# Computational Spectral-Angular Reconstruction Demo

This repository contains MATLAB demos for 2D and 3D reconstruction using an experimentally calibrated transmission matrix and FISTA-based inverse reconstruction. The demos only use simulated objects and simulated measurements; no experimental measurement data are required.

Only the `mixed_l1l2` reconstruction method is included in this demo release. Other legacy regularizers from development versions were intentionally removed to keep the upload package compact and focused.

## Required Content

- Source code for 2D spectral-angular reconstruction and 3D hyperspectral reconstruction.
- A demo calibration dataset stored as web-uploadable binary parts:
  - `data/measured_response_data.mat.part01`
  - `data/measured_response_data.mat.part02`
  - Variables: `response_xyw`, `theta_x_data`, `theta_y_data`, `wavelength_data`
  - The full `data/measured_response_data.mat` file is rebuilt automatically on the first demo run.
- A small image pattern used by the 3D simulated object:
  - `char_A.jpeg`
- Demo scripts:
  - `demo_2D_reconstruction.m`
  - `demo_3D_reconstruction.m`
- MATLAB reconstruction source code:
  - `src/common`
  - `src/rec2d`
  - `src/rec3d`

## 1. System Requirements

### Software

Tested environment:

- MATLAB R2022b
- Image Processing Toolbox
- Parallel Computing Toolbox is optional, but recommended for GPU acceleration


### Hardware

No non-standard hardware is required for the demos.

Recommended hardware:

- 16 GB RAM or more
- NVIDIA GPU with CUDA support for faster reconstruction

The demos will fall back to CPU if no compatible GPU is found.

## 2. Installation Guide

1. Download or clone the repository.
2. Open MATLAB.
3. Set the MATLAB current folder to the repository root, where this README is located.
4. Confirm that the calibration parts exist:

   ```text
   data/measured_response_data.mat.part01
   data/measured_response_data.mat.part02
   ```

   The original full calibration file is larger than the GitHub web upload limit, so it is not stored directly. The helper function `src/common/ensure_demo_calibration_file.m` joins the parts and creates `data/measured_response_data.mat` automatically when either demo is run.

## 3. Demo

### 2D Reconstruction Demo

Run:

```matlab
demo_2D_reconstruction
```

This demo:

1. Loads the experimentally calibrated response matrix.
2. Interpolates the response to the original 2D reconstruction grid:
   - `theta_x = -60:1:60 deg`
   - `wavelength = 400:0.5:700 nm`
3. Creates a simulated spectral-angular object.
   - The angular profile has one compact nonzero center position.
   - All nonzero angular positions share the same smooth spectrum generated from three random control points; the final normalized spectrum is constrained to stay within `[0.35, 1]`.
4. Simulates one encoded measurement spectrum.
5. Runs FISTA reconstruction.
6. Shows a live MATLAB window with reconstruction and residual loss.
7. Saves PNG figures and a MAT result file.

Expected outputs:

```text
outputs/2D_demo/2D_ground_truth.png
outputs/2D_demo/2D_measurement.png
outputs/2D_demo/2D_reconstruction_summary.png
outputs/2D_demo/2D_demo_result.mat
```

Expected runtime on a normal desktop computer:

- GPU: several minutes
- CPU: tens of minutes, depending on CPU speed

### 3D Reconstruction Demo

Run:

```matlab
demo_3D_reconstruction
```

This demo:

1. Loads the experimentally calibrated 3D response matrix.
2. Interpolates the response to the original 3D reconstruction grid:
   - `theta_x = linspace(-25, 25, 50)`
   - `theta_y = linspace(-25, 25, 50)`
   - `wavelength = 400:0.5:700 nm`
   - 36 rotation measurements
3. Creates a simulated hyperspectral object.
   - The spatial pattern is read from `char_A.jpeg`.
   - Every angular position shares the same smooth spectrum generated from three random control points; the final normalized spectrum is constrained to stay within `[0.35, 1]`.
4. Simulates encoded spectra for each rotation angle.
5. Runs FISTA reconstruction.
6. Shows a live MATLAB window with slices, spectra, and residual loss.
7. Saves PNG figures and a MAT result file.


### Using Your Own 2D Object

In `demo_2D_reconstruction.m`, replace:

```matlab
xw_gt = create_demo_2d_object(theta_x_range, wavelength_range);
```

with your own nonnegative array:

```matlab
xw_gt(theta_x_index, wavelength_index)
```

Then run the script. The measurement is simulated by:

```matlab
measurement = Forward_model(xw_gt, xw_sampling_matrix);
```

### Using Your Own 3D Object

In `demo_3D_reconstruction.m`, replace:

```matlab
obj = create_demo_3d_object(theta_x_range, theta_y_range, wavelength_range);
```

with your own nonnegative hyperspectral cube:

```matlab
obj(theta_x_index, theta_y_index, wavelength_index)
```

Then run the script. The encoded measurements are simulated by:

```matlab
measurements = Forward_model(obj, xywm_weight);
```
