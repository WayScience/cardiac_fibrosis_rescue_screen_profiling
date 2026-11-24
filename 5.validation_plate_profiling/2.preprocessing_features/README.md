# Preprocessing features

In this module, we turn the CellProfiler SQLite outputs of the validation plates into QC-filtered, position-corrected **single-cell** and **bulk** (well-level) profiles.
Each plate is processed on its own; the bulk step also pools the wells of all plates.
For more detail on each step, see the [module 4 README](../../4.image_based_profiling/README.md), which uses the same approach for the screen.

| Notebook | What it does | Output (per plate) |
| --- | --- | --- |
| [`0.convert_cytotable.ipynb`](./0.convert_cytotable.ipynb) | Merge the objects in each SQLite file into single cells with CytoTable, and add the `Metadata_` prefix to the location, cell count, and image FileName/PathName columns | `data/converted_profiles/*_converted.parquet` |
| [`1.sc_quality_control.ipynb`](./1.sc_quality_control.ipynb) | Flag poor-quality cells with coSMicQC `label_outliers` (all cells are kept, with a `Metadata_cqc_<condition>_is_outlier` column per condition) | `data/qc_labeled_profiles/*_qc_labeled.parquet` and figures in `qc_figures/` |
| [`2.single_cell_processing.ipynb`](./2.single_cell_processing.ipynb) | Drop QC-failed cells, annotate with the plate map, standardize using all cells on the plate, and select features | `*_sc_annotated`, `*_sc_normalized`, `*_sc_feature_selected` |
| [`3.bulk_processing.ipynb`](./3.bulk_processing.ipynb) | Drop QC-failed cells, aggregate to wells (median), annotate, normalize to the failing DMSO controls (`mad_robustize`), and select features for each plate; then pool the normalized wells of all plates and select features once | `*_bulk_annotated`, `*_bulk_normalized`, `*_bulk_feature_selected`, and `bulk_feature_selected` (all plates pooled) |
| [`4a.position_correction/`](./4a.position_correction/) | Diagnose the plate-position effect, fit a correction, validate it on held-out plates, and draw the figures of the screen (analysis, not run by the pipeline script) | `position_correction_fit.npz`, tables, and figures |
| [`4b.apply_position_correction.ipynb`](./4b.apply_position_correction.ipynb) | Subtract the fitted position effect from every cell of the normalized profiles of step 2, and select features | `*_sc_position_corrected`, `*_sc_position_corrected_feature_selected` |
| [`5.bulk_position_corrected.ipynb`](./5.bulk_position_corrected.ipynb) | Aggregate the position-corrected cells to wells (median), pool the wells of all plates, and select features once | `bulk_position_corrected_feature_selected` (all plates pooled) |

Steps 4a, 4b, and 5 are the steps 3a, 3b, and 4 of the screen ([module 4 README](../../4.image_based_profiling/README.md)).
They follow the existing steps and do not change their outputs.

## Single-cell QC

coSMicQC flags a cell when its z-score (computed within the plate) passes the threshold for every feature in a condition:

| Condition | Features |
| --- | --- |
| Over-segmented nuclei | `Nuclei_AreaShape_Solidity` (low) and `Nuclei_Intensity_MassDisplacement_DNA` (high) |
| High intensity nuclei | `Nuclei_Intensity_MeanIntensity_DNA` (high) |
| Under-segmented (small) cells | `Cells_AreaShape_Area` (low) |
| Blurry cells | `Cells_Granularity_1_Actin` (low) |
| Mis-segmented cells | `Cells_Intensity_MeanIntensity_Actin` (low) and `Cells_AreaShape_FormFactor` (low) |

The thresholds for each plate are in [`sc_qc_thresholds.json`](./sc_qc_thresholds.json).
Every plate needs an entry there, because coSMicQC fails when a condition has no thresholds.
About 11–15% of cells fail QC on each plate.

## Plate-position correction

Where a well sits on the plate shifts its profile in a consistent direction, and the validation plates show this effect strongly (the first principal component of the wells is 53% explained by plate row, and 57% to 72% in the DMSO wells of one cell type alone).
We correct it as in the screen:

1. [`4a.position_correction/`](./4a.position_correction/) fits one tilt map (an additive row plus column effect for every feature) on the compound and DMSO wells of the five plates, with one effect for each of the 15 compounds and each kind of DMSO because each compound is on one plate only (the column effects carry a tiny ridge penalty, because the groups sit in blocks of columns), measures one amplitude for each plate and shrinks it toward the other plates as the screen does for its platemaps, and validates it on held-out plates.
   Its [README](./4a.position_correction/README.md) has the method, the figures, and the results.
2. [`4b.apply_position_correction.ipynb`](./4b.apply_position_correction.ipynb) subtracts the amplitude of its plate times the tilt of its well from every cell (all cells of a well shift by the same vector, controls included), then selects features with the settings of step 2.
3. [`5.bulk_position_corrected.ipynb`](./5.bulk_position_corrected.ipynb) takes the median of every well of the corrected cells, pools the wells of all plates, and selects features once.
   The single-cell profiles are already standardized on each plate, so it does not normalize again, as in the bulk step of the screen.
   Step 3 normalizes the bulk profiles of the uncorrected cells to the failing DMSO controls instead, and its outputs stay as they are.

The feature names are those of the validation plates (the d2 channel is `Golgi`).

## Where outputs are saved

All profiles are saved to `data/`, which git ignores.
The single-cell and bulk notebooks save each output locally with the `output_file` parameter of the pycytominer functions, then save a copy to bandicoot with pycytominer's `output()` function:

- local: `data/single_cell_profiles/` and `data/bulk_profiles/`
- bandicoot: the same folder names under `~/mnt/bandicoot/CFReT_screening_data/validation_profiles/`

Steps 4b and 5 save locally only, with pycytominer's `output()` function and the `output_file` parameter of `feature_select`, and do not write to bandicoot.
Their file names are new (`*_sc_position_corrected*` and `bulk_position_corrected_feature_selected`), so they do not overwrite the profiles of the earlier steps.
They read the profiles of step 2 from `data/single_cell_profiles/` when it exists and from bandicoot otherwise.

## Run the module

Run all steps in order with the bash script [preprocess_features.sh](./preprocess_features.sh):

```bash
source preprocess_features.sh
```

The script converts the notebooks to Python scripts and runs them.
QC runs once per plate, with the plate set through the `PLATE_ID` environment variable.
As scripts, the QC notebooks do not render the example cell images, so QC takes only seconds per plate.
To look at example outlier cells, open [`1.sc_quality_control.ipynb`](./1.sc_quality_control.ipynb) and set the default `plate_id` in the parameters cell.

The script runs steps 4b and 5 after step 3, and it does not run step 4a.
Step 4b reads the committed fit, [`4a.position_correction/position_correction_fit.npz`](./4a.position_correction/position_correction_fit.npz).
Run the notebooks of step 4a (see its [README](./4a.position_correction/README.md)) to refit the correction, for example when the normalized profiles change.
