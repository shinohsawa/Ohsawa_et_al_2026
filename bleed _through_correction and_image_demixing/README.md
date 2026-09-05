# Bleedthrough Correction Pipeline

Standalone tool for correcting spectral bleedthrough in multi-channel fluorescence microscopy images.

## Overview

This pipeline:
1. **Loads** ND2 images and matching segmentation masks (.npz format)
2. **Pools** voxels inside the mask from both channels
3. **Fits** a linear model excluding low-intensity noise (bottom 10% by default)
4. **Generates** scatter plots showing the fit quality
5. **Corrects** the dependent channel by subtracting predicted bleedthrough
6. **Exports** multi-channel OME-TIFF files with all original channels plus the corrected channel

## Quick Start

### 1. Set up conda environment

```bash
# Create environment from the included specification
conda env create -f environment.yml

# Activate the environment
conda activate demixing
```

Alternatively, create manually:
```bash
conda create -n demixing python=3.10 numpy scipy matplotlib nd2 xarray tifffile -c conda-forge
conda activate demixing
```

### 2. Configure the analysis

Edit `config.py` to set:
- **Input folders**: where to find ND2 images and masks
- **Output folders**: where to save corrected images and plots
- **Channel indices**: which channels to use (0-based indexing)
- **Processing parameters**: background subtraction and fit thresholds

**Important**: Check the channel indices match your data! See the comments in `config.py` for guidance.

### 3. Run the pipeline

```bash
# Activate the environment (if not already active)
conda activate demixing

# Run the analysis
python demix_images.py
```

The script will:
- Print progress for each image
- Show fit statistics (slope, intercept, R²)
- Save outputs to the configured folders

## Output Files

### Corrected Images
Location: `output/corrected_tiff/` (or as configured)

Multi-channel OME-TIFF files containing:
- All original channels (unchanged)
- Corrected dependent channel (replaces original dependent channel)

File naming: `<original_name>_corrected.ome.tif`

### Scatter Plots
Location: `output/scatter_plots/` (or as configured)

PDF files showing:
- Scatter plot of dependent vs independent channel intensities
- Fitted linear regression line (red)
- Threshold used for fitting (orange dashed line)
- Fit statistics (slope, intercept, R²)

File naming: `<original_name>_scatter.pdf`

## Configuration Details

### Channel Selection

```python
'dependent_channel_index': 1,    # Green channel (contaminated)
'independent_channel_index': 2,  # Red channel (source of bleedthrough)
```

**How to find channel indices:**
1. Open an ND2 file in ImageJ/Fiji
2. Note the channel numbers shown (1-based)
3. Subtract 1 for the config (0-based)

Example: If ImageJ shows "C2-green", use index `1` in config.

### Processing Parameters

#### Background Percentile
```python
'background_percentile': 10.0
```
- Controls how much background is removed before fitting
- Lower = more aggressive (use ~5-10 for clean data)
- Higher = more conservative (use ~15-25 for noisy data)

**Troubleshooting:**
- If scatter plot shows upward curve at low intensities → increase percentile
- If too few points remain → decrease percentile

#### Threshold Fraction
```python
'threshold_fraction': 0.1
```
- Excludes low-intensity voxels from the linear fit
- `0.1` = exclude bottom 10% of intensity range
- Prevents noise from biasing the fit

**Troubleshooting:**
- If fit line doesn't capture the data well → decrease threshold
- If low-intensity noise affects fit → increase threshold

### Reference Source (Advanced)

```python
'reference_source': None
```

Controls which images are used to calculate the global bleedthrough correction model.

**Default behavior (`None`):**
- Uses ALL ND2 files in `nd2_folder` to compute the fit
- All images contribute voxels to the pooled dataset for fitting

**Custom reference selection (dictionary):**
Restrict the fit to specific reference images while still correcting all images in the dataset.

```python
'reference_source': {
    'nd2_folder': '../reference_images',      # Optional: different folder
    'mask_folder': '../reference_masks',      # Optional: different folder
    'include_patterns': ['*_bleedthrough.nd2'] # Required: glob patterns
}
```

**Dictionary keys:**
- **`nd2_folder`** (optional, string): Path to folder containing reference ND2 files
  - If omitted, uses the main `nd2_folder` value
  - Use this when reference images are in a separate location
  
- **`mask_folder`** (optional, string): Path to folder containing reference masks
  - If omitted, uses the main `mask_folder` value
  - Must contain matching `.npz` masks for the reference ND2 files
  
- **`include_patterns`** (optional, list of strings): Glob patterns to select specific files
  - If omitted, defaults to `['*.nd2']` (all ND2 files in the reference folder)
  - Patterns are evaluated inside the `nd2_folder`
  - Examples:
    - `['sample1.nd2']` — single specific file
    - `['*_control.nd2']` — all files ending with `_control.nd2`
    - `['exp1_*.nd2', 'exp2_*.nd2']` — multiple patterns

**Common use cases:**

1. **Use one image from the main dataset as reference:**
   ```python
   'reference_source': {
       'include_patterns': ['my_reference_image.nd2']
   }
   ```
   - No need to specify folders (defaults to main paths)
   - Only this image is used to calculate slope/intercept
   - All images in `nd2_folder` are still corrected

2. **Use multiple images from the main dataset:**
   ```python
   'reference_source': {
       'include_patterns': ['*_control.nd2', 'baseline_*.nd2']
   }
   ```

3. **Use images from a separate reference folder:**
   ```python
   'reference_source': {
       'nd2_folder': '../bleedthrough_references',
       'mask_folder': '../bleedthrough_references/masks',
       'include_patterns': ['*.nd2']
   }
   ```

**Important notes:**
- Reference images MUST have matching masks in `mask_folder`
- The global fit is computed only from reference image voxels
- Correction is applied to ALL images in the main `nd2_folder`
- The pipeline will exit with error if a reference image has no matching mask

### Correction Offset (Advanced)

```python
'correction_offset_mode': 'none',
'correction_offset_value': 0.0,
```

Allows you to shift the bleedthrough correction line upward to be more conservative and avoid over-correction in the presence of noise.

**Why use an offset?**
- The linear fit represents the *mean* bleedthrough relationship
- Due to noise, some true signal may appear in the bleedthrough distribution
- An offset creates a "safety margin" by subtracting less than the fitted model predicts
- Helps prevent removing real signal that happens to correlate with bleedthrough

**`correction_offset_mode` options:**

1. **`'none'`** (default)
   - No offset applied
   - Uses the fitted intercept directly
   - Most aggressive correction
   - Best when noise is minimal

2. **`'residual_std'`**
   - Adds `correction_offset_value × σ` to the intercept
   - Where σ = standard deviation of fit residuals (reported in console output)
   - The offset is *data-driven* and scales with your noise level
   - **Recommended** for most cases with noisy data
   
   Example: `correction_offset_value: 1.0` shifts the line up by 1 standard deviation
   - ~68% of points fall below the correction line
   - Conservative correction that preserves most potential signal

3. **`'fixed'`**
   - Adds exactly `correction_offset_value` intensity units to the intercept
   - The offset is *absolute* and independent of data properties
   - Use when you know a specific intensity value to reserve
   
   Example: `correction_offset_value: 100.0` shifts intercept up by exactly 100 units

**`correction_offset_value` (float):**
- Scaling factor for `'residual_std'` mode (typically 0.5 to 3.0)
- Absolute intensity units for `'fixed'` mode
- Ignored when `correction_offset_mode: 'none'`

**Visual feedback:**
- The scatter plot shows TWO lines when offset is active:
  - **Red solid line** = best-fit line (not used for correction)
  - **Blue dashed line** = offset line (used for actual correction)
- Legend shows the offset amount applied
- Console output reports both intercepts and the offset amount

**Recommended values for `'residual_std'` mode:**
- `0.5` — minimal offset, slight safety margin
- `1.0` — moderate offset, balances correction vs. preservation
- `1.5-2.0` — conservative offset, prioritizes preserving signal
- `>2.0` — very conservative, may under-correct

**Examples:**

Conservative correction (1 standard deviation):
```python
'correction_offset_mode': 'residual_std',
'correction_offset_value': 1.0,
```

Very conservative (2 standard deviations):
```python
'correction_offset_mode': 'residual_std',
'correction_offset_value': 2.0,
```

Fixed offset of 200 intensity units:
```python
'correction_offset_mode': 'fixed',
'correction_offset_value': 200.0,
```

**Troubleshooting:**
- If offset is near zero, consider using `'fixed'` mode instead
- If you're over-correcting (losing real signal), increase the offset value
- If you're under-correcting (residual bleedthrough visible), decrease the offset value
- Check the scatter plot to visually verify the offset line position

## Mask Matching

The script automatically matches masks to images by filename:

| Image filename | Mask filename | Match type |
|---------------|---------------|------------|
| `image.nd2` | `image.npz` | Exact |
| `image.nd2` | `image_cells.npz` | Prefix |
| `image_.nd2` | `image.npz` | Trimmed |

If multiple masks match, it picks the one with the shortest suffix after the image name.

## Troubleshooting

### "No matching mask found"
- Check that mask filenames match or prefix-match the ND2 filenames
- Masks must be `.npz` format
- Check the `mask_folder` path in `config.py`

### "Channel index out of range"
- Open ND2 in ImageJ and count the channels
- Adjust `dependent_channel_index` and `independent_channel_index`
- Remember: indices are 0-based (first channel = 0)

### Poor fit quality (low R²)
- Check scatter plot to see the relationship between channels
- Adjust `background_percentile` if background isn't fully removed
- Adjust `threshold_fraction` to exclude more low-intensity noise
- Verify you selected the correct channel indices

### "Required libraries not found"
- Activate the conda environment: `conda activate demixing`
- Or install manually: `conda install numpy scipy matplotlib nd2 xarray tifffile -c conda-forge`

## Technical Details

### Algorithm

1. **Load & Select**: Read ND2, select first timepoint if present
2. **Background Subtraction**: Per-plane percentile-based subtraction
3. **Masking**: Broadcast 2D segmentation mask across Z-stack
4. **Pooling**: Extract all voxels inside mask from both channels
5. **Thresholding**: Compute `threshold = fraction × max(independent)`
6. **Fitting**: Linear regression on voxels above threshold
7. **Correction**: `corrected = dependent - (slope × independent + intercept)`
8. **Export**: Save as OME-TIFF preserving metadata

### Linear Model

```
dependent = slope × independent + intercept
corrected_dependent = dependent - (slope × independent + intercept)
```

The slope represents the bleedthrough coefficient (fraction of red signal appearing in green channel).

### File Formats

- **Input**: ND2 (Nikon microscopy format), NPZ (NumPy compressed arrays)
- **Output**: OME-TIFF (Open Microscopy Environment TIFF with metadata)

## Citation

If you use this pipeline, please cite:
- nd2 library: https://github.com/tlambert03/nd2
- OME-TIFF specification: https://docs.openmicroscopy.org/ome-model/

## Support

For issues or questions:
1. Check the scatter plots to diagnose fit quality
2. Review the configuration parameters in `config.py`
3. Verify input file formats and naming conventions
4. Check that the conda environment is activated

---
**Author**: Generated for spectral unmixing pipeline  
**Last updated**: 2025-10-03
