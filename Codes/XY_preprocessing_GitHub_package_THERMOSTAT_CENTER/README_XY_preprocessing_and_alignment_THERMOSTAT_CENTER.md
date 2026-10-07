# XY trajectory preprocessing and thermal-plate alignment

This code converts cleaned tracking output from multiple independent thermotaxis assay plates into the `X`, `Y`, and plate-ID matrices used by the downstream manuscript analyses.

## Why alignment is required

The agar/gel assay sits on a metal thermal slab with two fixed thermostats positioned at the left and right sides of the gel. The distance between these thermostats is fixed, but the **position of the metal slab within the camera image can vary between experiments**. As a result, the same physical position along the thermal gradient can occur at different raw image X coordinates across assay plates.

To place all experiments in a common coordinate system, the midpoint between the two thermostats is measured in the image for each assay plate. Each recording is then translated along X so that this midpoint is located at **X = 1350 pixels**.

The correction is:

```text
x_shift_px = 1350 - measured_thermostat_center_px
X_aligned  = X_raw + x_shift_px
```

For example:

```text
measured thermostat center = 1315 px
x_shift_px = 1350 - 1315 = +35 px
```

After correction, the thermostat midpoint is at X = 1350 px for every assay plate. Y coordinates are unchanged.

This is a translational registration only; trajectories are not rotated, rescaled, warped, or otherwise transformed.

## Files

- `preprocess_xy_align_thermal_gradient_THERMOSTAT_CENTER.m` — repository-facing preprocessing script
- `collect_xy_from_folders_correctX_plateID_GITHUB.m` — extraction/alignment function
- `plate_alignment_manifest_example_THERMOSTAT_CENTER.csv` — example plate manifest

## Alignment manifest

The manifest contains one row per independent biological assay plate:

```csv
plate_folder,thermostat_center_px
relative/path/to/plate_01,1315
relative/path/to/plate_02,1362
relative/path/to/plate_03,1350
```

The script computes the required X shift automatically from the measured thermostat midpoint.

Keeping the measured physical landmark in the manifest is preferable to storing only the final shift because it makes the alignment procedure directly auditable.

## Raw tracking format

Each assay folder is expected to contain exactly one file matching:

```text
*postclean_experiment_1*.mat
```

The MAT file must contain `experiment_1.track`. For each track, the code uses:

- `npts`
- `startFrame`
- `endFrame`
- `dq.sloc`

`dq.sloc(1,:)` is treated as X and `dq.sloc(2,:)` as Y.

The tracking output uses zero-based `startFrame`/`endFrame`; these are converted to MATLAB one-based indices when trajectories are placed into the common output timeline.

## Track processing

A track is retained when:

1. the required tracking fields are present;
2. `npts >= min_track_frames`;
3. at least two coordinate rows are available in `dq.sloc`;
4. at least part of the track lies within the requested assay timeline.

Each retained trajectory is stored as one column in the pooled matrices. Time points outside the observed lifetime of a track remain `NaN`. No interpolation is performed.

## Plate identity

Each manifest row represents one independent assay plate. Plate identity is preserved for every pooled trajectory column. The original within-plate track index and folder name are also retained.

## Main settings

Edit these values in the main preprocessing script:

```matlab
raw_data_root = "/path/to/raw_tracking_data";
manifest_file = "plate_alignment_manifest.csv";
output_dir = "./processed_data";
dataset_label = "N2";

reference_thermostat_center_px = 1350;
fps = 2;
total_minutes = 15;
min_track_frames = 300;
```

At 2 frames/s, a 300-frame minimum corresponds to 150 s of tracked data.

## Outputs

For `dataset_label = "N2"`, the script writes:

```text
X_N2.mat
Y_N2.mat
plate_N2.mat
N2_track_metadata.csv
N2_plate_alignment_summary.csv
N2_alignment_manifest_used.csv
N2_processed_XY_with_metadata.mat
N2_alignment_QC.pdf
```

The first three files are directly compatible with the downstream Figure 5 / Figure 6 analysis scripts.

The saved alignment manifest includes both the measured thermostat midpoint and the calculated X shift for every plate.

## Reproducibility summary

- Each assay folder is one independent biological plate.
- The midpoint between the two fixed thermostats is used as the physical alignment landmark.
- Every recording is translated so that this midpoint is at X = 1350 px.
- Alignment is performed only by adding a constant X offset.
- Y is unchanged.
- No rotation, scaling, nonlinear transformation, or interpolation is applied.
- Missing/unobserved tracking samples remain `NaN`.
- The exact plate-level alignment measurements and calculated offsets are saved with the processed dataset.
