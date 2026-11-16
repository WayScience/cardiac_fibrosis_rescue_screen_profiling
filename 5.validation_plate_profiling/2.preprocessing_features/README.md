# Preprocessing features

In this module, we turn the CellProfiler SQLite outputs of the validation plates into QC-filtered **single-cell** and **bulk** (well-level) profiles.
Each plate is processed on its own.
For more detail on each step, see the [module 4 README](../../4.image_based_profiling/README.md), which uses the same approach for the screen.

| Notebook | What it does | Output (per plate) |
| --- | --- | --- |
| [`0.convert_cytotable.ipynb`](./0.convert_cytotable.ipynb) | Merge the objects in each SQLite file into single cells with CytoTable, and add the `Metadata_` prefix to the location, cell count, and image FileName/PathName columns | `data/converted_profiles/*_converted.parquet` |
| [`1.sc_quality_control.ipynb`](./1.sc_quality_control.ipynb) | Flag poor-quality cells with coSMicQC (all cells are kept, with a `Metadata_cqc_failed_*` column per condition) | `data/qc_labeled_profiles/*_qc_labeled.parquet` and figures in `qc_figures/` |
| [`2.single_cell_processing.ipynb`](./2.single_cell_processing.ipynb) | Drop QC-failed cells, annotate with the plate map, standardize using all cells on the plate, and select features | `*_sc_annotated`, `*_sc_normalized`, `*_sc_feature_selected` |
| [`3.bulk_processing.ipynb`](./3.bulk_processing.ipynb) | Drop QC-failed cells, aggregate to wells (median), annotate, normalize to the failing DMSO controls (`mad_robustize`), select features, and spherize to the same controls | `*_bulk_annotated`, `*_bulk_normalized`, `*_bulk_feature_selected`, `*_bulk_spherized` |

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

## Where outputs are saved

Converted and QC-labeled profiles are saved to `data/`, which git ignores.
The single-cell and bulk notebooks save their profiles with `save_parquet()` from [`utils/io_utils.py`](../../utils/io_utils.py).
Set `save_to` at the top of each notebook to choose where the profiles go:

- `"local"`: `data/single_cell_profiles/` and `data/bulk_profiles/`
- `"bandicoot"`: the same folder names under `~/mnt/bandicoot/CFReT_screening_data/validation_profiles/`
- `"both"` (the current setting): both locations

## Run the module

Run all steps in order with the bash script [preprocess_features.sh](./preprocess_features.sh):

```bash
source preprocess_features.sh
```

The script converts the notebooks to Python scripts and runs them.
QC runs once per plate, with the plate set through the `PLATE_ID` environment variable.
As scripts, the QC notebooks do not render the example cell images, so QC takes only seconds per plate.
To look at example outlier cells, open [`1.sc_quality_control.ipynb`](./1.sc_quality_control.ipynb) and set the default `plate_id` in the parameters cell.
