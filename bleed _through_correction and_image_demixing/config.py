"""
Configuration file for demixing pipeline.

Edit the values in the CONFIG dictionary below to customize the analysis.
"""

CONFIG = {
    # ============================================================================
    # INPUT/OUTPUT PATHS
    # ============================================================================
    
    # Folder containing ND2 image files
    'nd2_folder': '../microscopy/251002/Images/raw_microscopy_files/',
    
    # Folder containing .npz segmentation mask files
    # Masks should have filenames that match (or prefix-match) the ND2 files
    'mask_folder': '../microscopy/251002/Images/cell_masks/',
    
    # Output folder for corrected OME-TIFF files
    'output_folder': '../microscopy/251002/Images/output/corrected_tiff',
    
    # Output folder for scatter plots with fit overlays
    'scatter_output_folder': '../microscopy/251002/Images/output/scatter_plots',
    
    # Optional override for selecting reference images used to fit the
    # bleedthrough model. Leave as None to use all ND2 files from
    # nd2_folder. Example configuration:
    # 'reference_source': {
    #     'nd2_folder': '../reference_images',
    #     'mask_folder': '../reference_masks',
    #     'include_patterns': ['*_bleedthrough.nd2']
    # }
    'reference_source': {
        'include_patterns': ['251002_3584_ED_006.nd2']
    },

    
    # ============================================================================
    # CHANNEL SELECTION
    # ============================================================================
    
    # Index of the DEPENDENT channel (typically green, the channel being corrected)
    # This is the channel that suffers from bleedthrough contamination
    # Channel indices are 0-based (first channel = 0, second = 1, etc.)
    'dependent_channel_index': 1,
    
    # Index of the INDEPENDENT channel (typically red, the source of bleedthrough)
    # This is the channel that bleeds into the dependent channel
    'independent_channel_index': 2,
    
    
    # ============================================================================
    # PREPROCESSING PARAMETERS
    # ============================================================================
    
    # Background subtraction percentile (0-100)
    # Pixels below this percentile are considered background and subtracted
    # Lower values = more aggressive background removal
    # Typical range: 5-25
    # - Use ~10 for clean backgrounds
    # - Use ~25 for noisy backgrounds
    'background_percentile': 10.0,
    
    
    # ============================================================================
    # FIT PARAMETERS
    # ============================================================================
    
    # Threshold fraction for excluding low-intensity voxels from fit (0.0-1.0)
    # Voxels where independent channel < (threshold_fraction × max_intensity) are excluded
    # This prevents low-intensity noise from biasing the linear fit
    # 
    # How it works:
    # - 0.0 = include all voxels (not recommended, noise affects fit)
    # - 0.1 = exclude bottom 10% of intensity range (default, good starting point)
    # - 0.2 = exclude bottom 20% of intensity range (more aggressive)
    # 
    # Increase this if:
    # - Your scatter plot shows a curved/non-linear relationship at low intensities
    # - Background subtraction alone doesn't remove low-intensity noise
    # 
    # Decrease this if:
    # - You're losing too many valid signal voxels
    # - The fit line doesn't capture the bulk of your data
    'threshold_fraction': 0.1,

    # Optional offset applied to the correction model.
    # - 'none': no offset, uses the fitted line directly
    # - 'residual_std': shift by (correction_offset_value × residual std of fit)
    # - 'fixed': shift by a fixed intensity value (same units as the data)
    'correction_offset_mode': 'residual_std',
    'correction_offset_value': 1.0,
    
    
    # ============================================================================
    # EXPORT OPTIONS
    # ============================================================================
    
    # TIFF compression: None, 'zlib', or 'lzw'
    # - None: Fastest, largest files
    # - 'zlib': Good compression, moderate speed
    # - 'lzw': Less compression, faster than zlib
    # Set to None (no quotes) for no compression
    'tiff_compression': None,
}


# ============================================================================
# ADVANCED NOTES
# ============================================================================

# CHANNEL INDICES:
# To find the correct channel indices for your data:
# 1. Open an ND2 file in ImageJ/Fiji with Bio-Formats
# 2. Note which channel number corresponds to green (dependent) and red (independent)
# 3. Subtract 1 from the channel number (ImageJ uses 1-based indexing)
# 
# Example:
# - If green is "Channel 2" in ImageJ → use index 1
# - If red is "Channel 3" in ImageJ → use index 2

# MASK MATCHING:
# The script automatically matches masks to images by filename:
# - Exact match: "image.nd2" → "image.npz"
# - Prefix match: "image.nd2" → "image_segmentation.npz"
# - Ignores trailing underscores
# 
# If multiple masks match, it picks the one with the shortest suffix.

# BACKGROUND PERCENTILE:
# This parameter controls how much background signal is removed before fitting.
# The percentile is computed per Z-plane for Z-stacks.
# 
# If your scatter plot shows:
# - Upward curve at low intensities → try higher percentile (15-25)
# - Too few points remaining → try lower percentile (5-10)

# THRESHOLD FRACTION:
# After background subtraction, this parameter excludes low-intensity voxels
# from the linear regression fit. The threshold is computed as:
#   threshold = threshold_fraction × max(independent_channel)
# 
# Orange dashed line in scatter plot shows this threshold.
# Only points to the right of this line are used for fitting.

# REFERENCE SOURCE:
# Use reference_source to limit which images contribute to the global fit.
# Keys:
# - nd2_folder (optional): folder containing the reference ND2 files.
#   Defaults to nd2_folder when omitted.
# - mask_folder (optional): folder containing masks for the reference images.
#   Defaults to mask_folder.
# - include_patterns (optional): list of glob patterns evaluated inside nd2_folder.
#   Defaults to ['*.nd2']. Specify an explicit filename (e.g., ['my_image.nd2'])
#   to use a single image.
# Example:
# 'reference_source': {
#     'nd2_folder': '../reference_images',
#     'mask_folder': '../reference_masks',
#     'include_patterns': ['*_ref.nd2']
# }

# CORRECTION OFFSET:
# Use correction_offset_mode to shift the fitted line upward before subtraction.
# Modes:
# - 'none': no shift (default)
# - 'residual_std': shift by correction_offset_value × residual standard deviation
# - 'fixed': shift by a fixed intensity amount specified in correction_offset_value
# Example to subtract a line two residual standard deviations above the fit:
# 'correction_offset_mode': 'residual_std'
# 'correction_offset_value': 2.0
