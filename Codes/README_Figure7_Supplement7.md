# Figure 7 / Supplement 7 reproducibility analysis

This reduced MATLAB script contains only the tracking analyses used for the data-derived panels of Figure 7 and Supplement 7.

## Script

`Figure7_Supplement7_reproducibility.m`

## Panels generated

### Figure 7

- **B** — plate-mean physical position through the 30-min gradient-switch assay, shown as mean ± SEM across independent plates
- **C** — individual wild-type trajectories
- **D** — individual `unc-17` cKO trajectories

**Figure 7A is a schematic and is not generated from tracking data.**

### Supplement 7

- **A** — instantaneous animal-position histograms in consecutive 3-min windows
- **B** — individual `eat-4` cKO trajectories

### Supplement 7C

The supplied MATLAB source does not define the rule used to classify animals as **“Left”** versus **“Stayed”**, nor the statistical analysis used for that categorical endpoint. This reduced script therefore does **not** infer or reproduce Supplement 7C. That panel should be added only after its exact classification rule is explicitly documented.

## Expected input

The script expects:

```text
X_N2.mat
Y_N2.mat
plate_N2.mat

X_eat4.mat
Y_eat4.mat
plate_eat4.mat

X_unc17.mat
Y_unc17.mat
plate_unc17.mat
```

Rows are time points and columns are animal tracks. Each plate vector contains one independent assay-plate identifier per trajectory column.

## Assay geometry

The aligned physical arena is fixed at:

```text
left edge  = 900 px
right edge = 1800 px
```

Before the gradient switch:

```text
900 px  = 15 °C
1800 px = 25 °C
```

The switch begins at 15 min and is completed by approximately 18 min.

After the switch:

```text
900 px  = 23 °C
1800 px = 17 °C
```

All trajectory and mean-position analyses use physical X coordinate internally. Temperature values are used only to label the corresponding physical coordinate system. This prevents an artificial positional discontinuity when the gradient is inverted.

## Figure 7B

Each independent plate first contributes one mean-X trajectory by averaging all currently valid animal positions on that plate at each time point. Genotype means and SEM are then calculated across independent plates.

For statistical analysis, each plate contributes one mean physical X position in consecutive 3-min bins. A linear mixed-effects model is fitted:

```text
MeanX_px ~ Genotype * TimeBin + (1 | PlateUID)
```

The genotype-by-time-bin interaction tests whether the temporal position profile differs among genotypes while accounting for repeated observations from the same plate.

Planned pairwise genotype comparisons are also performed within each 3-min bin using Wilcoxon rank-sum tests on independent plate means. Holm correction is applied across the three genotype comparisons within each time bin.

## Supplement 7A

Every valid instantaneous animal position within each 3-min interval is included in the histogram. These distributions are descriptive; repeated video frames from the same animal are not treated as independent biological replicates.

Histogram bins are defined in fixed physical X coordinates. Before and after the switch, those same physical positions are labeled with the corresponding temperature calibration. During 15–18 min the X axis is shown in pixels because the gradient is actively changing and there is no single fixed temperature mapping for that interval.

## Individual trajectory panels

Each thin line corresponds to one tracked animal. Genotype-specific shade gradients are used only to help distinguish individual trajectories visually.

A valid spatial observation requires both X and Y coordinates to be finite. No interpolation is performed.

## Outputs

The script produces manuscript-panel PDFs plus source-data and statistics CSV files:

```text
Figure7B_plate_mean_position_SEM.pdf
Figure7C_wild_type_trajectories.pdf
Figure7D_unc17_trajectories.pdf
Supplement7A_instantaneous_position_histograms.pdf
Supplement7B_eat4_trajectories.pdf

Figure7B_plate_mean_position_source_data.csv
Figure7B_plate_mean_position_3min_source_data.csv
Figure7B_global_statistics.csv
Figure7B_3min_pairwise_statistics.csv
Supplement7A_instantaneous_position_source_data.csv
```

A compact MAT file containing the derived analysis objects is also saved.
