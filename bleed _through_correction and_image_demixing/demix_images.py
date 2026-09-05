"""
Standalone bleedthrough correction script for ND2 images with cell segmentation masks.

This script:
1. Loads ND2 images and matching .npz segmentation masks
2. Pools voxels inside the mask from both channels
3. Fits a linear model (dependent ~ slope * independent + intercept)
   - Excludes low-intensity voxels (bottom 10% of independent channel max)
4. Generates scatter plot with fit overlay
5. Creates corrected dependent channel using the fitted model
6. Exports multi-channel OME-TIFF with original + corrected channels
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Tuple, Optional, List
import numpy as np


def load_mask(mask_path: Path) -> np.ndarray:
    """Load a segmentation mask from .npz and project to 2D if needed."""
    arr = np.load(mask_path)
    if hasattr(arr, 'files'):
        key = arr.files[0]
        data = arr[key]
    else:
        data = arr
    data = np.asarray(data)
    if data.ndim == 3:
        # Project across Z dimension (any positive voxel)
        projected = np.any(data > 0, axis=0)
    elif data.ndim == 2:
        projected = data > 0
    else:
        raise ValueError(f"Unsupported mask dimensionality: {data.shape}")
    return projected.astype(bool)


def find_masks(mask_folder: Path) -> dict[str, Path]:
    """Build a mapping of mask stem names to file paths."""
    mapping: dict[str, Path] = {}
    for path in sorted(mask_folder.glob("*.npz")):
        mapping[path.stem] = path
    return mapping


def gather_nd2_files(folder: Path, patterns: Optional[List[str]] = None) -> List[Path]:
    """Return ordered list of ND2 files matching optional patterns."""
    files: List[Path] = []
    if patterns:
        for pattern in patterns:
            files.extend(sorted(folder.glob(pattern)))
    else:
        files = sorted(folder.glob("*.nd2"))
    unique: List[Path] = []
    seen = set()
    for path in files:
        try:
            key = path.resolve()
        except FileNotFoundError:
            continue
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def match_mask_for_image(image_stem: str, mask_map: dict[str, Path]) -> Optional[Path]:
    """
    Find the best matching mask for an image by stem name.
    
    Tries exact match first, then prefix matching with preference for 
    masks that differ only by trailing underscores.
    """
    # Exact match
    exact = mask_map.get(image_stem)
    if exact:
        return exact
    
    # Try trimmed version (remove trailing underscores)
    trimmed = image_stem.rstrip("_")
    if trimmed != image_stem:
        exact_trim = mask_map.get(trimmed)
        if exact_trim:
            return exact_trim
    
    # Prefix matching
    candidates: List[Tuple[Path, str]] = []
    for mask_stem, path in mask_map.items():
        if mask_stem.startswith(image_stem):
            remainder = mask_stem[len(image_stem):]
            candidates.append((path, remainder))
        elif trimmed and mask_stem.startswith(trimmed):
            remainder = mask_stem[len(trimmed):]
            candidates.append((path, remainder))
    
    if not candidates:
        return None
    
    # Prefer masks with underscore suffix, then by shorter remainder
    candidates.sort(key=lambda item: (
        0 if item[1].startswith("_") or item[1] == "" else 1,
        len(item[1]),
    ))
    return candidates[0][0]


def read_nd2_image(nd2_path: Path) -> 'xarray.DataArray':
    """Read ND2 file using nd2 library and return xarray DataArray."""
    try:
        import nd2
        import xarray as xr
    except ImportError as e:
        raise ImportError(
            "Required libraries not found. Install with:\n"
            "  conda install -c conda-forge nd2 xarray\n"
            "or:\n"
            "  pip install nd2 xarray"
        ) from e
    
    with nd2.ND2File(nd2_path) as f:
        da = f.to_xarray(squeeze=False)
    return da


def select_channel(da: 'xarray.DataArray', channel_idx: int) -> 'xarray.DataArray':
    """Extract a single channel from multi-channel DataArray."""
    if 'C' not in da.dims:
        raise ValueError("No channel dimension 'C' found in data")
    channel_coord = da.coords['C'].values[channel_idx]
    return da.sel(C=channel_coord)


def apply_percentile_background_subtraction(
    img: np.ndarray, 
    percentile: float = 10.0
) -> np.ndarray:
    """Subtract per-plane background estimated by percentile."""
    img = img.astype(np.float32)
    
    # Squeeze out any singleton dimensions first
    img = np.squeeze(img)
    
    if img.ndim == 2:
        bg = float(np.percentile(img, percentile))
        return np.maximum(img - bg, 0.0)
    elif img.ndim == 3:
        # Z-stack: apply per Z-plane
        result = np.empty_like(img, dtype=np.float32)
        for zi in range(img.shape[0]):
            plane = img[zi]
            bg = float(np.percentile(plane, percentile))
            result[zi] = np.maximum(plane - bg, 0.0)
        return result
    else:
        raise ValueError(f"Expected 2D or 3D image after squeezing, got shape {img.shape}")


def pool_masked_voxels(
    dep_stack: np.ndarray,
    ind_stack: np.ndarray,
    mask_2d: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pool voxels from 3D stacks using 2D mask broadcast across Z.
    
    Returns:
        dep_values: 1D array of dependent channel values
        ind_values: 1D array of independent channel values
    """
    if dep_stack.ndim == 2:
        dep_stack = dep_stack[None, :, :]
    if ind_stack.ndim == 2:
        ind_stack = ind_stack[None, :, :]
    
    # Broadcast 2D mask to 3D
    mask_3d = np.broadcast_to(mask_2d, dep_stack.shape)
    
    dep_vals = dep_stack[mask_3d]
    ind_vals = ind_stack[mask_3d]
    
    # Keep only finite values
    valid = np.isfinite(dep_vals) & np.isfinite(ind_vals)
    return dep_vals[valid], ind_vals[valid]


def fit_linear_model_with_threshold(
    y: np.ndarray,
    x: np.ndarray,
    threshold_fraction: float = 0.1,
) -> Tuple[float, float, float, float]:
    """
    Fit linear model: y = slope * x + intercept
    
    Excludes points where x < threshold_fraction * max(x) to avoid 
    low-intensity noise affecting the fit.
    
    Returns:
        slope, intercept, r2, threshold_value
    """
    if x.size < 2:
        return float('nan'), float('nan'), float('nan'), float('nan')
    
    # Compute threshold
    x_max = float(x.max())
    threshold = threshold_fraction * x_max
    
    # Filter data
    mask = x >= threshold
    x_fit = x[mask]
    y_fit = y[mask]
    
    if x_fit.size < 2:
        return float('nan'), float('nan'), float('nan'), threshold
    
    # Fit linear model: y = slope * x + intercept
    A = np.vstack([x_fit, np.ones_like(x_fit)]).T
    slope, intercept = np.linalg.lstsq(A, y_fit, rcond=None)[0]
    
    # Compute R²
    y_pred = slope * x_fit + intercept
    ss_res = float(np.sum((y_fit - y_pred) ** 2))
    ss_tot = float(np.sum((y_fit - y_fit.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    
    return float(slope), float(intercept), float(r2), float(threshold)


def compute_offset_intercept(
    x: np.ndarray,
    y: np.ndarray,
    slope: float,
    intercept: float,
    offset_mode: str,
    offset_value: float,
    threshold: Optional[float] = None,
) -> Tuple[float, float, float]:
    """Return (offset_intercept, residual_std, offset_amount)."""
    mask = np.ones_like(x, dtype=bool)
    if threshold is not None:
        mask &= x >= threshold

    x_sel = x[mask]
    y_sel = y[mask]
    if x_sel.size < 2:
        x_sel = x
        y_sel = y

    residuals = y_sel - (slope * x_sel + intercept)
    if residuals.size:
        positive = residuals[residuals > 0]
        base_sample = positive if positive.size >= 2 else residuals
        residual_std = float(np.std(base_sample))
    else:
        residual_std = 0.0
    mode = offset_mode.lower()
    if mode == 'residual_std':
        offset = float(offset_value) * residual_std
    elif mode == 'fixed':
        offset = float(offset_value)
    else:
        offset = 0.0
    return intercept + offset, residual_std, offset


def create_scatter_plot(
    x: np.ndarray,
    y: np.ndarray,
    slope: float,
    intercept: float,
    r2: float,
    threshold: float,
    output_path: Path,
    max_points: int = 50000,
    image_name: str = "",
    offset_intercept: Optional[float] = None,
    offset_label: Optional[str] = None,
):
    """Generate scatter plot with linear fit and optional offset overlay."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not installed, skipping scatter plot generation")
        return
    
    # Subsample for plotting if needed
    if x.size > max_points:
        idx = np.random.choice(x.size, size=max_points, replace=False)
        x_plot = x[idx]
        y_plot = y[idx]
    else:
        x_plot = x
        y_plot = y
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Scatter plot
    ax.scatter(x_plot, y_plot, alpha=0.3, s=1, c='gray', rasterized=True)
    
    # Fit line
    x_line = np.array([x.min(), x.max()])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'r-', linewidth=2, 
            label=f'y = {slope:.4f}x + {intercept:.2f}\nR² = {r2:.4f}')

    if offset_intercept is not None and np.isfinite(offset_intercept):
        y_line_offset = slope * x_line + offset_intercept
        label = offset_label or f'Offset line (intercept {offset_intercept:.2f})'
        ax.plot(
            x_line,
            y_line_offset,
            color='royalblue',
            linestyle='--',
            linewidth=2.5,
            label=label,
        )
    
    # Threshold line
    ax.axvline(threshold, color='orange', linestyle='--', linewidth=1.5,
               label=f'Threshold = {threshold:.1f}')
    
    ax.set_xlabel('Independent channel (red)', fontsize=12)
    ax.set_ylabel('Dependent channel (green)', fontsize=12)
    ax.set_title(f'Bleedthrough correction fit\n{image_name}', fontsize=14)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Scatter plot: {output_path}")


def export_ome_tiff(
    original_data: 'xarray.DataArray',
    corrected_channel: np.ndarray,
    dep_channel_idx: int,
    output_path: Path,
    nd2_path: Path,
    compression: Optional[str] = None,
):
    """
    Export multi-channel OME-TIFF with corrected dependent channel.
    
    All original channels are preserved, with the dependent channel replaced
    by the corrected version.
    """
    try:
        import tifffile as tiff
    except ImportError:
        raise ImportError(
            "tifffile not installed. Install with:\n"
            "  conda install -c conda-forge tifffile\n"
            "or:\n"
            "  pip install tifffile"
        )
    
    # Ensure standard dimension order
    da = original_data
    
    # Handle position/series dimension (select first if present)
    if 'P' in da.dims:
        da = da.isel(P=0)
    
    if 'T' not in da.dims:
        da = da.expand_dims({'T': [0]})
    if 'Z' not in da.dims:
        da = da.expand_dims({'Z': [0]})
    if 'C' not in da.dims:
        da = da.expand_dims({'C': [0]})
    
    da = da.transpose('T', 'Z', 'C', 'Y', 'X')
    
    Tn, Zn, Cn, Yn, Xn = [da.sizes[d] for d in ('T', 'Z', 'C', 'Y', 'X')]
    
    # Get channel names
    if 'C' in original_data.coords:
        channel_names = [str(v) for v in original_data.coords['C'].values]
    else:
        channel_names = [f"Channel{i}" for i in range(Cn)]
    
    # Prepare output array
    dtype = np.asarray(da).dtype
    out_arr = np.empty((Tn, Zn, Cn, Yn, Xn), dtype=dtype)
    
    # Copy all channels, replacing dependent with corrected version
    for ti in range(Tn):
        for zi in range(Zn):
            for ci in range(Cn):
                if ci == dep_channel_idx:
                    # Use corrected channel
                    if corrected_channel.ndim == 2:
                        out_arr[ti, zi, ci] = corrected_channel
                    else:
                        out_arr[ti, zi, ci] = corrected_channel[zi]
                else:
                    # Copy original
                    out_arr[ti, zi, ci] = np.asarray(da.isel(T=ti, Z=zi, C=ci))
    
    # Build OME metadata
    ome_meta = {
        'axes': 'TZCYX',
        'Channel': {'Name': channel_names},
    }
    
    # Try to extract pixel sizes from ND2 metadata
    psx, psy, psz = None, None, None
    try:
        import nd2
        with nd2.ND2File(nd2_path) as f:
            vs = f.voxel_size()
            if isinstance(vs, dict):
                psx = vs.get('x') or vs.get('X')
                psy = vs.get('y') or vs.get('Y')
                psz = vs.get('z') or vs.get('Z')
            elif isinstance(vs, (list, tuple)) and len(vs) >= 2:
                if len(vs) == 2:
                    psy, psx = float(vs[0]), float(vs[1])
                    psz = None
                else:
                    psz, psy, psx = float(vs[-3]), float(vs[-2]), float(vs[-1])
    except Exception:
        pass
    
    if psx is not None:
        ome_meta['PhysicalSizeX'] = float(psx)
        ome_meta['PhysicalSizeXUnit'] = 'µm'
    if psy is not None:
        ome_meta['PhysicalSizeY'] = float(psy)
        ome_meta['PhysicalSizeYUnit'] = 'µm'
    if psz is not None:
        ome_meta['PhysicalSizeZ'] = float(psz)
        ome_meta['PhysicalSizeZUnit'] = 'µm'
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tiff.imwrite(
        str(output_path),
        out_arr,
        photometric='minisblack',
        metadata=ome_meta,
        compression=compression,
    )
    print(f"  OME-TIFF: {output_path}")


def process_single_image(
    nd2_path: Path,
    mask_path: Path,
    dep_channel_idx: int,
    ind_channel_idx: int,
    bg_percentile: float,
    threshold_fraction: float,
    output_folder: Path,
    scatter_output_folder: Path,
    tiff_compression: Optional[str],
) -> dict:
    """
    Process a single ND2 image with mask-based bleedthrough correction.
    
    Returns a dictionary with fit parameters and statistics.
    """
    print(f"\nProcessing: {nd2_path.name}")
    print(f"  Mask: {mask_path.name}")
    
    # Load data
    da = read_nd2_image(nd2_path)
    mask_2d = load_mask(mask_path)
    
    # Select first position/series if present
    if 'P' in da.dims:
        da = da.isel(P=0)
    
    # Select time point if present
    if 'T' in da.dims:
        da = da.isel(T=0)
    
    # Extract channels
    dep_da = select_channel(da, dep_channel_idx)
    ind_da = select_channel(da, ind_channel_idx)
    
    dep_arr = np.asarray(dep_da).astype(np.float32)
    ind_arr = np.asarray(ind_da).astype(np.float32)
    
    # Squeeze out singleton dimensions (e.g., position/series)
    dep_arr = np.squeeze(dep_arr)
    ind_arr = np.squeeze(ind_arr)
    
    # Apply background subtraction
    dep_processed = apply_percentile_background_subtraction(dep_arr, bg_percentile)
    ind_processed = apply_percentile_background_subtraction(ind_arr, bg_percentile)
    
    # Pool masked voxels
    dep_vals, ind_vals = pool_masked_voxels(dep_processed, ind_processed, mask_2d)
    
    print(f"  Pooled {len(dep_vals):,} voxels inside mask")
    
    # Fit linear model
    slope, intercept, r2, threshold = fit_linear_model_with_threshold(
        dep_vals, ind_vals, threshold_fraction
    )
    
    n_above_threshold = np.sum(ind_vals >= threshold)
    print(f"  Fit: slope={slope:.6f}, intercept={intercept:.2f}, R²={r2:.4f}")
    print(f"  Threshold: {threshold:.1f} ({n_above_threshold:,} voxels used for fit)")
    
    # Generate scatter plot
    scatter_path = scatter_output_folder / f"{nd2_path.stem}_scatter.pdf"
    create_scatter_plot(
        x=ind_vals,
        y=dep_vals,
        slope=slope,
        intercept=intercept,
        r2=r2,
        threshold=threshold,
        output_path=scatter_path,
        image_name=nd2_path.name,
        offset_intercept=None,
        offset_label=None,
    )
    
    # Apply correction to full dependent channel
    prediction = slope * ind_processed + intercept
    corrected = dep_arr - prediction
    
    # Clip to valid range if integer dtype
    original_dtype = np.asarray(dep_da).dtype
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        corrected = np.clip(np.rint(corrected), info.min, info.max).astype(original_dtype)
    else:
        corrected = corrected.astype(original_dtype)
    
    # Export OME-TIFF
    tiff_path = output_folder / f"{nd2_path.stem}_corrected.ome.tif"
    export_ome_tiff(
        da, corrected, dep_channel_idx, tiff_path, nd2_path, tiff_compression
    )
    
    return {
        'image': nd2_path.name,
        'mask': mask_path.name,
        'slope': slope,
        'intercept': intercept,
        'r2': r2,
        'threshold': threshold,
        'n_voxels_total': len(dep_vals),
        'n_voxels_fit': n_above_threshold,
    }


def main():
    """Main entry point - imports config and runs pipeline."""
    # Import configuration
    try:
        from config import CONFIG
    except ImportError:
        print("Error: config.py not found in the same directory!")
        print("Please ensure config.py exists and contains CONFIG dictionary.")
        sys.exit(1)
    
    # Extract parameters
    nd2_folder = Path(CONFIG['nd2_folder'])
    mask_folder = Path(CONFIG['mask_folder'])
    output_folder = Path(CONFIG['output_folder'])
    scatter_folder = Path(CONFIG['scatter_output_folder'])
    
    dep_ch = CONFIG['dependent_channel_index']
    ind_ch = CONFIG['independent_channel_index']
    bg_percentile = CONFIG['background_percentile']
    threshold_fraction = CONFIG['threshold_fraction']
    tiff_compression = CONFIG.get('tiff_compression')
    reference_config = CONFIG.get('reference_source')
    offset_mode = str(CONFIG.get('correction_offset_mode', 'none')).lower()
    offset_value = float(CONFIG.get('correction_offset_value', 0.0))
    valid_offset_modes = {'none', 'residual_std', 'fixed'}
    if offset_mode not in valid_offset_modes:
        print(f"Warning: unknown correction_offset_mode '{offset_mode}', defaulting to 'none'")
        offset_mode = 'none'
    
    # Validate inputs
    if not nd2_folder.exists():
        print(f"Error: ND2 folder not found: {nd2_folder}")
        sys.exit(1)
    if not mask_folder.exists():
        print(f"Error: Mask folder not found: {mask_folder}")
        sys.exit(1)
    
    # Find files
    nd2_files = gather_nd2_files(nd2_folder)
    if not nd2_files:
        print(f"No ND2 files found in {nd2_folder}")
        sys.exit(1)
    
    mask_map = find_masks(mask_folder)

    # Resolve reference selection
    if reference_config:
        ref_nd2_folder = Path(reference_config.get('nd2_folder', nd2_folder))
        ref_mask_folder = Path(reference_config.get('mask_folder', mask_folder))
        include_patterns = reference_config.get('include_patterns')
        if isinstance(include_patterns, str):
            include_patterns = [include_patterns]
        if not include_patterns:
            include_patterns = ['*.nd2']
    else:
        ref_nd2_folder = nd2_folder
        ref_mask_folder = mask_folder
        include_patterns = None

    if not ref_nd2_folder.exists():
        print(f"Error: Reference ND2 folder not found: {ref_nd2_folder}")
        sys.exit(1)
    if not ref_mask_folder.exists():
        print(f"Error: Reference mask folder not found: {ref_mask_folder}")
        sys.exit(1)

    if ref_mask_folder.resolve() == mask_folder.resolve():
        ref_mask_map = mask_map
        if not mask_map:
            print(f"No .npz mask files found in {mask_folder}")
            sys.exit(1)
    else:
        ref_mask_map = find_masks(ref_mask_folder)
        if not ref_mask_map:
            print(f"No .npz mask files found in {ref_mask_folder}")
            sys.exit(1)
        if not mask_map:
            print(f"Warning: no .npz mask files found in {mask_folder}")

    reference_nd2_files = gather_nd2_files(ref_nd2_folder, include_patterns)
    if not reference_nd2_files:
        print("No ND2 files matched the reference selection.")
        sys.exit(1)
    
    print(f"Found {len(nd2_files)} ND2 files for correction and {len(mask_map)} masks")
    print(f"\nConfiguration:")
    print(f"  Dependent channel: {dep_ch}")
    print(f"  Independent channel: {ind_ch}")
    print(f"  Background percentile: {bg_percentile}")
    print(f"  Threshold fraction: {threshold_fraction}")
    print(f"  Correction offset mode: {offset_mode}")
    print(f"  Correction offset value: {offset_value}")
    if reference_config:
        pattern_str = ', '.join(include_patterns) if include_patterns else '*.nd2'
        print("\nReference fit source:")
        print(f"  ND2 folder: {ref_nd2_folder}")
        print(f"  Mask folder: {ref_mask_folder}")
        print(f"  Patterns: {pattern_str}")
        print(f"  Selected images: {len(reference_nd2_files)}")
    else:
        print("\nReference fit source: using all ND2 files in nd2_folder")
    
    # Create output folders
    output_folder.mkdir(parents=True, exist_ok=True)
    scatter_folder.mkdir(parents=True, exist_ok=True)
    
    # ========================================================================
    # PHASE 1: Pool all voxels from all images to compute global fit
    # ========================================================================
    print("\n" + "="*60)
    print("PHASE 1: Pooling voxels from reference images")
    print("="*60)
    
    all_dep_vals = []
    all_ind_vals = []
    reference_image_names = []
    
    for nd2_path in reference_nd2_files:
        mask_path = match_mask_for_image(nd2_path.stem, ref_mask_map)
        if mask_path is None:
            print(f"Error: reference image {nd2_path.name} has no matching mask")
            sys.exit(1)
        
        try:
            print(f"\nLoading: {nd2_path.name}")
            print(f"  Mask: {mask_path.name}")
            
            # Load data
            da = read_nd2_image(nd2_path)
            mask_2d = load_mask(mask_path)
            
            # Select first position/series if present
            if 'P' in da.dims:
                da = da.isel(P=0)
            
            # Select time point if present
            if 'T' in da.dims:
                da = da.isel(T=0)
            
            # Extract channels
            dep_da = select_channel(da, dep_ch)
            ind_da = select_channel(da, ind_ch)
            
            dep_arr = np.asarray(dep_da).astype(np.float32)
            ind_arr = np.asarray(ind_da).astype(np.float32)
            
            # Squeeze out singleton dimensions
            dep_arr = np.squeeze(dep_arr)
            ind_arr = np.squeeze(ind_arr)
            
            # Apply background subtraction
            dep_processed = apply_percentile_background_subtraction(dep_arr, bg_percentile)
            ind_processed = apply_percentile_background_subtraction(ind_arr, bg_percentile)
            
            # Pool masked voxels
            dep_vals, ind_vals = pool_masked_voxels(dep_processed, ind_processed, mask_2d)
            
            print(f"  Pooled {len(dep_vals):,} voxels")
            
            all_dep_vals.append(dep_vals)
            all_ind_vals.append(ind_vals)
            reference_image_names.append(nd2_path.name)
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    if not all_dep_vals:
        print("\nNo data collected. Exiting.")
        sys.exit(1)
    
    # Concatenate all voxels
    print(f"\n" + "="*60)
    print("Computing global linear fit")
    print("="*60)
    
    dep_pooled = np.concatenate(all_dep_vals)
    ind_pooled = np.concatenate(all_ind_vals)
    
    print(f"Total pooled voxels: {len(dep_pooled):,}")
    
    # Fit global model
    slope, intercept, r2, threshold = fit_linear_model_with_threshold(
        dep_pooled, ind_pooled, threshold_fraction
    )

    # Compute residual statistics for optional offset line
    intercept_with_offset, residual_std, offset_amount = compute_offset_intercept(
        ind_pooled, dep_pooled, slope, intercept, offset_mode, offset_value, threshold
    )

    n_above_threshold = np.sum(ind_pooled >= threshold)
    print(f"\nGlobal fit:")
    print(f"  Slope: {slope:.6f}")
    print(f"  Intercept: {intercept:.2f}")
    print(f"  Residual std: {residual_std:.2f}")
    if offset_mode != 'none':
        print(f"  Offset mode: {offset_mode} ({offset_value})")
        print(f"  Offset applied: {offset_amount:.2f}")
        if np.isclose(offset_amount, 0.0):
            print("  Note: offset evaluated to ~0. Consider 'fixed' mode or larger value if you expect a visible shift.")
        print(f"  Intercept w/ offset: {intercept_with_offset:.2f}")
    print(f"  R²: {r2:.4f}")
    print(f"  Threshold: {threshold:.1f}")
    print(f"  Voxels used for fit: {n_above_threshold:,} / {len(ind_pooled):,}")
    
    # Generate global scatter plot
    print(f"\nGenerating global scatter plot...")
    global_scatter_path = scatter_folder / "global_fit_all_images.pdf"
    offset_label = None
    if offset_mode != 'none':
        offset_label = f"Offset line ({offset_mode}, +{offset_amount:.2f})"

    create_scatter_plot(
        x=ind_pooled,
        y=dep_pooled,
        slope=slope,
        intercept=intercept,
        r2=r2,
        threshold=threshold,
        output_path=global_scatter_path,
        image_name="All images combined",
        offset_intercept=intercept_with_offset if offset_mode != 'none' else None,
        offset_label=offset_label,
    )
    
    # ========================================================================
    # PHASE 2: Apply correction to each image
    # ========================================================================
    print("\n" + "="*60)
    print("PHASE 2: Applying correction to individual images")
    print("="*60)
    
    results = []
    
    for nd2_path in nd2_files:
        try:
            print(f"\nCorrecting: {nd2_path.name}")
            
            da = read_nd2_image(nd2_path)
            if 'P' in da.dims:
                da = da.isel(P=0)
            if 'T' in da.dims:
                da = da.isel(T=0)

            dep_da = select_channel(da, dep_ch)
            ind_da = select_channel(da, ind_ch)

            dep_arr_raw = np.squeeze(np.asarray(dep_da))
            ind_arr = np.squeeze(np.asarray(ind_da)).astype(np.float32)

            dep_arr = dep_arr_raw.astype(np.float32)
            ind_processed = apply_percentile_background_subtraction(ind_arr, bg_percentile)
            
            prediction = slope * ind_processed + intercept_with_offset
            corrected = dep_arr - prediction
            
            # Clip to valid range if integer dtype
            original_dtype = dep_arr_raw.dtype
            if np.issubdtype(original_dtype, np.integer):
                info = np.iinfo(original_dtype)
                corrected = np.clip(np.rint(corrected), info.min, info.max).astype(original_dtype)
            else:
                corrected = corrected.astype(original_dtype)
            
            # Export OME-TIFF
            tiff_path = output_folder / f"{nd2_path.stem}_corrected.ome.tif"
            export_ome_tiff(
                da, corrected, dep_ch, tiff_path, nd2_path, tiff_compression
            )

            mask_path = match_mask_for_image(nd2_path.stem, mask_map)
            if mask_path is None:
                print("  Warning: no matching mask found; proceeding without mask metadata")
            results.append({
                'image': nd2_path.name,
                'mask': mask_path.name if mask_path else None,
            })
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Successfully processed: {len(results)} images")
    print(f"\nGlobal fit applied to all images:")
    print(f"  Slope: {slope:.6f}")
    print(f"  Intercept: {intercept:.2f}")
    if offset_mode != 'none':
        print(f"  Offset mode: {offset_mode} ({offset_value})")
        print(f"  Offset applied: {offset_amount:.2f}")
        if np.isclose(offset_amount, 0.0):
            print("  Note: offset evaluated to ~0. Consider 'fixed' mode or larger value if you expect a visible shift.")
        print(f"  Intercept w/ offset: {intercept_with_offset:.2f}")
    print(f"  Residual std: {residual_std:.2f}")
    print(f"  R²: {r2:.4f}")

    if reference_image_names:
        print(f"\nReference images used for fit ({len(reference_image_names)}):")
        for name in reference_image_names:
            print(f"  - {name}")
    
    if results:
        print("\nCorrected images:")
        for res in results:
            mask_info = f" (mask: {res['mask']})" if res['mask'] else ""
            print(f"  - {res['image']}{mask_info}")
    
    print(f"\nOutputs saved to:")
    print(f"  - OME-TIFF files: {output_folder}")
    print(f"  - Scatter plots: {scatter_folder}")


if __name__ == '__main__':
    main()
