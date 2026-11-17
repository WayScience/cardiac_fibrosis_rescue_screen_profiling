# Segmentation and feature extraction

In this module, we correct the illumination of the validation plate images, segment the cells, and extract morphology features with CellProfiler.
The pipeline is a modified version of the analysis pipeline from the [cellpainting_predicts_cardiac_fibrosis](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.

| Notebook | What it does | Output (per plate) |
| --- | --- | --- |
| [`0.create_loaddata_with_illum_csv.ipynb`](./0.create_loaddata_with_illum_csv.ipynb) | Add the paths to each channel's illumination function from [`0.illumination_correction`](../0.illumination_correction/) to that module's LoadData CSVs | `loaddata_csvs/loaddata_with_illum_<plate>.csv` |
| [`1.extract_features.ipynb`](./1.extract_features.ipynb) | Run the analysis pipeline on all plates in parallel | `cp_output/<plate>/cfret_features.sqlite`, outline images in the same folder, and a run log in `logs/` |

## CellProfiler pipeline

The pipeline in [`pipeline/analysis.cppipe`](./pipeline/analysis.cppipe) does the following for each image set:

1. Apply the illumination functions to the five channels (`Actin`, `DNA`, `ER`, `Golgi`, and `Mito`).
2. Segment three compartments:
   - **Nuclei** from the `DNA` channel
   - **Cells** by growing out from each nucleus in the `Actin` channel
   - **Cytoplasm** as the cell minus its nucleus
3. Measure features for each compartment and channel: colocalization, granularity, intensity, intensity distribution, neighbors, size and shape, and texture.
4. Save nuclei and cell outline images (`NucleiOutlines_*.tiff` and `CellsOutlines_*.tiff`), which the single-cell QC notebook uses to show example cells.
5. Export all measurements to one SQLite file per plate.

The `d2` channel is named `Golgi` here, while the screen pipeline calls it `PM`; see the [module README](../README.md#channels).

## Run the module

Run [`0.illumination_correction`](../0.illumination_correction/) first.
Then run the bash script [perform_extraction.sh](./perform_extraction.sh), which activates the `fibrosis_cp_env` conda environment, converts the notebooks to Python scripts, and runs them:

```bash
source perform_extraction.sh
```

The SQLite outputs in `cp_output/` are large, so git ignores them.
[`2.preprocessing_features`](../2.preprocessing_features/) converts them to parquet files.
The outputs for the deprecated preliminary plate are in `preliminary_validation_plate/`, which git ignores.
