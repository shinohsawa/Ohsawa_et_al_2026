#!/usr/bin/env python3
"""
3D Subcellular Object Relabeling Tool - Batch Processing

This tool relabels 3D subcellular objects based on which cell they belong to in 2D.
All subcellular objects within the same cell are merged into one discontinuous object
with the cell's ID, perfect for total volume analysis.

Workflow:
1. Scan folders for 3D and 2D segmentation files
2. Match files by name stem (e.g., "sample1_manual.npz" pairs with "sample1_2d.npz")
3. For each pair:
   - Load 3D subcellular segmentation and 2D cell segmentation
   - Project 3D objects to 2D (union of XY across all Z slices)
   - Find overlap between projected subcellular objects and cells
   - Create relabeling: all objects in same cell → same ID (the cell ID)
   - Apply relabeling to original 3D data
   - Save relabeled 3D, 2D projection, mapping, stats, and visualization

Usage:
1. Edit the CONFIGURATION section in main() to set:
   - subcell_dir: folder with 3D subcellular .npz files
   - cell_dir: folder with 2D cell .npz files (can be same)
   - subcell_3d_pattern: pattern to identify 3D files (e.g., "*manual*", "*3d*")
   - cell_2d_pattern: pattern to identify 2D files (e.g., "*2d*", "*cells*")
   - output_dir: where to save results
2. Run: python relabel_subcellular_objects.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import ndimage
from collections import defaultdict
import traceback

def find_file_pairs(subcell_dir, cell_dir, subcell_3d_pattern='*manual*', cell_2d_pattern='*2d*'):
    """
    Find matching pairs of 3D subcellular and 2D cell segmentation files
    
    Args:
        subcell_dir: Directory containing 3D subcellular segmentation files
        cell_dir: Directory containing 2D cell segmentation files
        subcell_3d_pattern: Pattern to identify 3D subcellular files (e.g., '*manual*', '*3d*')
        cell_2d_pattern: Pattern to identify 2D cell files (e.g., '*2d*', '*cells*')
    
    Returns:
        List of tuples: [(subcell_3d_file, cell_2d_file, common_stem), ...]
    """
    subcell_dir = Path(subcell_dir)
    cell_dir = Path(cell_dir)
    
    print(f"Searching for file pairs...")
    print(f"  3D subcellular directory: {subcell_dir}")
    print(f"  2D cell directory: {cell_dir}")
    print(f"  3D pattern: {subcell_3d_pattern}")
    print(f"  2D pattern: {cell_2d_pattern}")
    
    # Find all npz files matching patterns
    subcell_files = list(subcell_dir.glob(f"**/{subcell_3d_pattern}.npz"))
    cell_files = list(cell_dir.glob(f"**/{cell_2d_pattern}.npz"))
    
    print(f"\nFound {len(subcell_files)} 3D subcellular files")
    print(f"Found {len(cell_files)} 2D cell files")
    
    if len(subcell_files) == 0:
        print(f"WARNING: No 3D files found matching pattern '{subcell_3d_pattern}' in {subcell_dir}")
    if len(cell_files) == 0:
        print(f"WARNING: No 2D files found matching pattern '{cell_2d_pattern}' in {cell_dir}")
    
    # Match files by extracting common stem (before pattern markers)
    pairs = []
    
    for subcell_file in subcell_files:
        # Extract the base name and try to find matching cell file
        subcell_name = subcell_file.stem
        
        # Remove pattern-specific parts to get common stem
        # E.g., "3506_6h_EST__s1_segm_manual" -> "3506_6h_EST__s1"
        common_stem = subcell_name
        
        # Try to remove common suffixes
        for suffix in ['_segm_manual', '_manual', '_3d', '_subcell']:
            if suffix in common_stem:
                common_stem = common_stem.replace(suffix, '')
        
        # Find matching cell file
        best_match = None
        best_match_score = 0
        
        for cell_file in cell_files:
            cell_name = cell_file.stem
            
            # Calculate match score (how many characters match from the start)
            match_len = 0
            for i in range(min(len(common_stem), len(cell_name))):
                if common_stem[i] == cell_name[i]:
                    match_len += 1
                else:
                    break
            
            # Also check if common_stem is a substring
            if common_stem in cell_name or cell_name.startswith(common_stem):
                match_len = max(match_len, len(common_stem))
            
            if match_len > best_match_score:
                best_match_score = match_len
                best_match = cell_file
        
        if best_match and best_match_score > 5:  # Require at least 5 character match
            pairs.append((subcell_file, best_match, common_stem))
            print(f"  ✓ Matched: {subcell_file.name} <-> {best_match.name}")
        else:
            print(f"  ✗ No match found for: {subcell_file.name}")
    
    print(f"\n✓ Found {len(pairs)} file pairs")
    
    return pairs

def load_segmentation_files(subcell_3d_file, cell_2d_file):
    """Load and validate segmentation files"""
    print("Loading segmentation files...")
    
    # Load files
    subcell_data = np.load(subcell_3d_file)
    cell_data = np.load(cell_2d_file)
    
    # Extract main arrays (assume first key contains the segmentation)
    subcell_key = list(subcell_data.keys())[0]
    cell_key = list(cell_data.keys())[0]
    
    subcell_3d = subcell_data[subcell_key]
    cell_2d = cell_data[cell_key]
    
    print(f"Loaded 3D subcellular: {subcell_3d.shape} (key: '{subcell_key}')")
    print(f"Loaded 2D cells: {cell_2d.shape} (key: '{cell_key}')")
    
    # Validate dimensions
    if len(subcell_3d.shape) != 3:
        raise ValueError(f"Subcellular data must be 3D, got {len(subcell_3d.shape)}D")
    if len(cell_2d.shape) != 2:
        raise ValueError(f"Cell data must be 2D, got {len(cell_2d.shape)}D")
    
    # Check XY compatibility
    if subcell_3d.shape[1:] != cell_2d.shape:
        raise ValueError(f"XY dimensions don't match: 3D={subcell_3d.shape[1:]} vs 2D={cell_2d.shape}")
    
    print("✓ Dimensions validated")
    
    subcell_data.close()
    cell_data.close()
    
    return subcell_3d, cell_2d

def project_3d_to_2d(subcell_3d):
    """Project 3D objects to 2D by taking union across Z dimension"""
    print("Projecting 3D objects to 2D...")
    
    # Get unique object IDs (excluding background)
    object_ids = np.unique(subcell_3d)
    object_ids = object_ids[object_ids > 0]
    
    print(f"Found {len(object_ids)} unique objects")
    
    # Create 2D projection
    projection_2d = np.zeros(subcell_3d.shape[1:], dtype=subcell_3d.dtype)
    
    for obj_id in object_ids:
        # Find all pixels belonging to this object across all Z slices
        obj_mask_3d = subcell_3d == obj_id
        
        # Project to 2D by taking union (any Z slice where object exists)
        obj_mask_2d = np.any(obj_mask_3d, axis=0)
        
        # Set projected pixels to object ID
        projection_2d[obj_mask_2d] = obj_id
    
    print(f"✓ 2D projection created with {len(np.unique(projection_2d))-1} objects")
    
    return projection_2d

def find_object_cell_mapping(projected_subcells, cell_2d):
    """Find which cell each subcellular object belongs to"""
    print("Finding object-to-cell mapping...")
    
    object_ids = np.unique(projected_subcells)
    object_ids = object_ids[object_ids > 0]
    
    cell_ids = np.unique(cell_2d)
    cell_ids = cell_ids[cell_ids > 0]
    
    print(f"Mapping {len(object_ids)} objects to {len(cell_ids)} cells")
    
    mapping = {}
    overlap_stats = {}
    
    for obj_id in object_ids:
        # Get mask for this object
        obj_mask = projected_subcells == obj_id
        
        # Find overlapping cells
        overlapping_cells = cell_2d[obj_mask]
        overlapping_cells = overlapping_cells[overlapping_cells > 0]  # Remove background
        
        if len(overlapping_cells) == 0:
            # Object doesn't overlap with any cell - assign to background (0)
            mapping[obj_id] = 0
            overlap_stats[obj_id] = {'cell': 0, 'overlap_pixels': 0, 'total_pixels': np.sum(obj_mask)}
            print(f"  Object {obj_id}: no cell overlap -> background")
        else:
            # Find which cell has most overlap
            unique_cells, counts = np.unique(overlapping_cells, return_counts=True)
            best_cell = unique_cells[np.argmax(counts)]
            max_overlap = counts.max()
            total_pixels = np.sum(obj_mask)
            overlap_ratio = max_overlap / total_pixels
            
            mapping[obj_id] = best_cell
            overlap_stats[obj_id] = {
                'cell': best_cell, 
                'overlap_pixels': max_overlap, 
                'total_pixels': total_pixels,
                'overlap_ratio': overlap_ratio
            }
            
            print(f"  Object {obj_id}: -> Cell {best_cell} ({overlap_ratio:.1%} overlap)")
    
    return mapping, overlap_stats

def create_relabeling_scheme(mapping, cell_ids):
    """Create a systematic relabeling scheme"""
    print("Creating relabeling scheme...")
    
    # Group objects by their assigned cell
    cell_to_objects = defaultdict(list)
    for obj_id, cell_id in mapping.items():
        cell_to_objects[cell_id].append(obj_id)
    
    # Create new labeling scheme: All objects in same cell get the same ID (cell_id)
    # This creates discontinuous objects for total volume analysis
    relabeling_dict = {}
    
    for cell_id, objects in cell_to_objects.items():
        objects.sort()  # Ensure consistent ordering
        
        if cell_id == 0:
            # Background objects keep original IDs
            for i, obj_id in enumerate(objects):
                relabeling_dict[obj_id] = obj_id  # Keep original ID
        else:
            # All objects in the same cell get the cell ID
            for obj_id in objects:
                new_id = cell_id  # Simply use the cell ID
                relabeling_dict[obj_id] = new_id
                
        if cell_id == 0:
            print(f"  Background: {len(objects)} objects kept with original IDs")
        else:
            new_id = relabeling_dict[objects[0]]  # All get same ID
            print(f"  Cell {cell_id}: {len(objects)} objects -> All merged to ID {new_id}")
    
    return relabeling_dict

def apply_relabeling(subcell_3d, relabeling_dict):
    """Apply relabeling to 3D data"""
    print("Applying relabeling to 3D data...")
    
    relabeled_3d = np.zeros_like(subcell_3d)
    
    for old_id, new_id in relabeling_dict.items():
        mask = subcell_3d == old_id
        relabeled_3d[mask] = new_id
    
    original_objects = len(np.unique(subcell_3d)) - 1
    relabeled_objects = len(np.unique(relabeled_3d)) - 1
    
    print(f"✓ Relabeling complete: {original_objects} -> {relabeled_objects} objects")
    
    return relabeled_3d

def apply_relabeling_2d(projected_2d, relabeling_dict):
    """Apply relabeling to 2D projection"""
    print("Applying relabeling to 2D projection...")
    
    relabeled_2d = np.zeros_like(projected_2d)
    
    for old_id, new_id in relabeling_dict.items():
        mask = projected_2d == old_id
        relabeled_2d[mask] = new_id
    
    original_objects = len(np.unique(projected_2d)) - 1
    relabeled_objects = len(np.unique(relabeled_2d)) - 1
    
    print(f"✓ 2D relabeling complete: {original_objects} -> {relabeled_objects} objects")
    
    return relabeled_2d

def visualize_results(subcell_2d_proj, cell_2d, relabeled_2d_proj, output_dir, file_stem):
    """Create visualization of the relabeling process"""
    print("Creating visualization...")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Original projections
    axes[0, 0].imshow(subcell_2d_proj, cmap='tab20')
    axes[0, 0].set_title('Original Subcellular Objects (2D projection)')
    axes[0, 0].axis('off')
    
    axes[0, 1].imshow(cell_2d, cmap='tab10')
    axes[0, 1].set_title('Cell Segmentation (2D)')
    axes[0, 1].axis('off')
    
    # Overlay
    overlay = np.zeros((*cell_2d.shape, 3))
    overlay[:, :, 0] = (cell_2d > 0).astype(float) * 0.5  # Cells in red
    overlay[:, :, 1] = (subcell_2d_proj > 0).astype(float) * 0.7  # Objects in green
    axes[1, 0].imshow(overlay)
    axes[1, 0].set_title('Overlay (Red=Cells, Green=Objects)')
    axes[1, 0].axis('off')
    
    # Relabeled result
    axes[1, 1].imshow(relabeled_2d_proj, cmap='tab20')
    axes[1, 1].set_title('Relabeled Objects (Cell-based IDs)')
    axes[1, 1].axis('off')
    
    plt.tight_layout()
    
    viz_file = output_dir / f'{file_stem}_visualization.png'
    plt.savefig(viz_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Visualization saved to: {viz_file}")

def save_results(relabeled_3d, relabeled_2d_proj, relabeling_dict, overlap_stats, output_dir, file_stem):
    """Save results to files"""
    print("Saving results...")
    
    # Save relabeled 3D data
    relabeled_3d_file = output_dir / f'{file_stem}_relabeled_3d.npz'
    np.savez_compressed(relabeled_3d_file, segmentation=relabeled_3d)
    print(f"✓ Relabeled 3D data saved to: {relabeled_3d_file}")
    
    # Save relabeled 2D projection
    relabeled_2d_file = output_dir / f'{file_stem}_relabeled_2d.npz'
    np.savez_compressed(relabeled_2d_file, segmentation=relabeled_2d_proj)
    print(f"✓ Relabeled 2D projection saved to: {relabeled_2d_file}")
    
    # Save relabeling dictionary
    mapping_file = output_dir / f'{file_stem}_mapping.txt'
    with open(mapping_file, 'w') as f:
        f.write("# Original_ID -> New_ID (Cell_based)\n")
        for old_id, new_id in sorted(relabeling_dict.items()):
            f.write(f"{old_id} -> {new_id}\n")
    print(f"✓ Relabeling mapping saved to: {mapping_file}")
    
    # Save overlap statistics
    stats_file = output_dir / f'{file_stem}_overlap_stats.txt'
    with open(stats_file, 'w') as f:
        f.write("# Object_ID, Assigned_Cell, Overlap_Pixels, Total_Pixels, Overlap_Ratio\n")
        for obj_id, stats in sorted(overlap_stats.items()):
            f.write(f"{obj_id}, {stats['cell']}, {stats['overlap_pixels']}, "
                   f"{stats['total_pixels']}, {stats['overlap_ratio']:.3f}\n")
    print(f"✓ Overlap statistics saved to: {stats_file}")

def process_single_pair(subcell_3d_file, cell_2d_file, output_dir, file_stem):
    """Process a single pair of segmentation files"""
    print("\n" + "-"*60)
    print(f"Processing: {file_stem}")
    print("-"*60)
    
    try:
        # Step 1: Load data
        subcell_3d, cell_2d = load_segmentation_files(subcell_3d_file, cell_2d_file)
        
        # Step 2: Project 3D to 2D
        subcell_2d_proj = project_3d_to_2d(subcell_3d)
        
        # Step 3: Find object-cell mapping
        mapping, overlap_stats = find_object_cell_mapping(subcell_2d_proj, cell_2d)
        
        # Step 4: Create relabeling scheme
        cell_ids = np.unique(cell_2d)
        cell_ids = cell_ids[cell_ids > 0]
        relabeling_dict = create_relabeling_scheme(mapping, cell_ids)
        
        # Step 5: Apply relabeling to 3D data
        relabeled_3d = apply_relabeling(subcell_3d, relabeling_dict)
        
        # Step 6: Apply relabeling to 2D projection  
        relabeled_2d_proj = apply_relabeling_2d(subcell_2d_proj, relabeling_dict)
        
        # Visualize results
        visualize_results(subcell_2d_proj, cell_2d, relabeled_2d_proj, output_dir, file_stem)
        
        # Save results
        save_results(relabeled_3d, relabeled_2d_proj, relabeling_dict, overlap_stats, output_dir, file_stem)
        
        print(f"\n✓ Successfully processed: {file_stem}")
        print(f"  Original objects: {len(np.unique(subcell_3d))-1}")
        print(f"  Relabeled objects: {len(np.unique(relabeled_3d))-1}")
        print(f"  Cells found: {len(np.unique(cell_2d))-1}")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error processing {file_stem}: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    print("="*70)
    print("3D SUBCELLULAR OBJECT RELABELING TOOL - MULTI-POSITION MODE")
    print("="*70)

    # ==================== CONFIGURATION ====================
    base_dir = Path("/Users/sosawa/Documents/microscope_offline/251015/Images/output/corrected_tiff")

    # Position_1～Position_8 を処理
    position_start = 9
    position_end = 25

    # ファイル名パターン
    subcell_3d_pattern = "*segm_nuc_3D*"              # 3D核セグメンテーション
    cell_2d_pattern = "*segm_YeaZ_v2*"    # 2Dセグメンテーション

    # ========================================================

    print(f"Processing Positions {position_start} to {position_end}")
    print(f"Base directory: {base_dir}")
    print("="*70)

    total_success = 0
    total_failed = 0

    for pos_idx in range(position_start, position_end + 1):
        pos_dir = base_dir / f"Position_{pos_idx}" / "Images"
        output_dir = pos_dir  # ✅ 出力を Images 内に直接保存

        if not pos_dir.exists():
            print(f"✗ Skipping Position_{pos_idx} (folder not found: {pos_dir})")
            continue

        print(f"\n{'-'*60}")
        print(f"Scanning Position_{pos_idx}")
        print(f"Input folder: {pos_dir}")
        print(f"Output folder: {output_dir}")
        print(f"{'-'*60}")

        # ファイルペア探索
        file_pairs = find_file_pairs(
            subcell_dir=pos_dir,
            cell_dir=pos_dir,
            subcell_3d_pattern=subcell_3d_pattern,
            cell_2d_pattern=cell_2d_pattern
        )

        if len(file_pairs) == 0:
            print(f"✗ No matching file pairs found in Position_{pos_idx}")
            total_failed += 1
            continue

        success = 0
        failed = 0

        for subcell_file, cell_file, common_stem in file_pairs:
            try:
                ok = process_single_pair(subcell_file, cell_file, output_dir, common_stem)
                if ok:
                    success += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"✗ Error in Position_{pos_idx}: {e}")
                traceback.print_exc()
                failed += 1

        print(f"\nPosition_{pos_idx} Summary:")
        print(f"  Successful: {success}")
        print(f"  Failed: {failed}")
        total_success += success
        total_failed += failed

    print("\n" + "="*70)
    print("ALL POSITIONS COMPLETED")
    print("="*70)
    print(f"Total successful: {total_success}")
    print(f"Total failed: {total_failed}")
    print("="*70)


if __name__ == "__main__":
    main()
