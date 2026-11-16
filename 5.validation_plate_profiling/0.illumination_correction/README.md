# Illumination correction

In this module, we calculate an illumination correction (IC) function for each channel of each validation plate with CellProfiler.
The pipeline is a modified version of the IC pipeline from the [cellpainting_predicts_cardiac_fibrosis](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis) repository.
This module only calculates the functions; [`1.cellprofiler_processing`](../1.cellprofiler_processing/) applies them to the images before segmentation.

| Notebook | What it does | Output (per plate) |
| --- | --- | --- |
| [`0.create_loaddata_csv.ipynb`](./0.create_loaddata_csv.ipynb) | Find the raw images, and build a CellProfiler LoadData CSV with one row per image set (plate, well, and site) and one file name and path per channel | `loaddata_csvs/loaddata_<plate>.csv` |
| [`1.illum_correct.ipynb`](./1.illum_correct.ipynb) | Run the IC pipeline on all plates in parallel | `illum_directory/<plate>/<plate>_Illum<channel>.npy` and a run log in `logs/` |

## Channels

The LoadData CSVs name the five channels as follows (see the [module README](../README.md#channels) for how these differ from the screen):

| Image channel | CellProfiler image name |
| --- | --- |
| d0 | `OrigActin` |
| d1 | `OrigMito` |
| d2 | `OrigGolgi` |
| d3 | `OrigER` |
| d4 | `OrigDNA` |

## CellProfiler pipeline

The pipeline in [`pipeline/illum.cppipe`](./pipeline/illum.cppipe) does the following for each channel:

1. Downsize the images to speed up the calculation.
2. Calculate one illumination function from all images on the plate, using a median filter.
3. Resize the function back to the full image size and save it as a `.npy` file.

The pipeline also contains whole-image QC modules (`MeasureImageQuality` and `FlagImage`) that would flag over-saturated images, but they are turned off, so every image set is used.

## Run the module

The raw image directory is set in both notebooks (`IMAGE_INPUTS` and `base_dir`); change it if the images are mounted somewhere else on your machine.
Then run the bash script [perform_ic.sh](./perform_ic.sh), which activates the `fibrosis_cp_env` conda environment, converts the notebooks to Python scripts, and runs them:

```bash
source perform_ic.sh
```

The outputs for the deprecated preliminary plate are in `preliminary_validation_plate/`, which git ignores.
