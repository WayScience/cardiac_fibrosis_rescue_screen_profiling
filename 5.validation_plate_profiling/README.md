# Validation plate profiling

This module processes the validation plates with the same workflow as modules 2, 3, and 4, using one plate map per plate in [`metadata/platemaps/`](./metadata/platemaps/).
There are five validation plates. They share the same well layout (controls and compound wells in the same positions), but each plate tests three different compounds:

| Plate | Barcode | Compound A | Compound B | Compound C |
| --- | --- | --- | --- | --- |
| 1 | `CARD-CelIns-CX7_261030160001` | K310 | K543 | K278 |
| 2 | `CARD-CelIns-CX7_261031160001` | K756 | K020 | K845 |
| 3 | `CARD-CelIns-CX7_261101160003` | K463 | K208 | K779 |
| 4 | `CARD-CelIns-CX7_261102080001` | K221 | K318 | K645 |
| 5 | `CARD-CelIns-CX7_261102150001` | K504 | K687 | K191 |

The compounds on these plates are drug hits from the screen, selected with our optimized hit-calling methodology.

An earlier preliminary validation plate (`CARD-CelIns-CX7_260803130001`) is deprecated.
Its compounds were drug hits selected with an earlier hit-calling method that was not reliable, so it did not test the right set of hits.
Its local outputs are kept in each step's `preliminary_validation_plate/` folder, which git ignores.

## Channels

The validation plates use the same five-channel stain as the screen, but the d2 channel has a different name:

| Image channel | Stain | Name in the validation plates | Name in the screen (modules 2–4) |
| --- | --- | --- | --- |
| d0 | F-actin | `Actin` | `Actin` |
| d1 | Mitochondria | `Mito` | `Mito` |
| d2 | Golgi/plasma membrane | `Golgi` | `PM` |
| d3 | Endoplasmic reticulum | `ER` | `ER` |
| d4 | Nuclei | `DNA` | `DNA` |

So features from the d2 channel are named `*_Golgi` in the validation profiles (e.g., `Cells_Intensity_MeanIntensity_Golgi`) and `*_PM` in the screen profiles.
The validation profiles are analyzed on their own, so the names are left as they are; rename the channel first if you ever compare these profiles to the screen directly.
The blocklist in [`2.preprocessing_features/blocklist_features.txt`](./2.preprocessing_features/blocklist_features.txt) uses the `Golgi` name.

## Steps

| Step | Folder | What it does |
| --- | --- | --- |
| 0 | [`0.illumination_correction`](./0.illumination_correction/) | Create the LoadData files and correct illumination with CellProfiler |
| 1 | [`1.cellprofiler_processing`](./1.cellprofiler_processing/) | Segment cells and extract morphology features with CellProfiler |
| 2 | [`2.preprocessing_features`](./2.preprocessing_features/) | Convert to parquet, run single-cell QC, process single cells, and aggregate to bulk profiles |

[`metadata/`](./metadata/) holds the validation plate map and its visualization.

Run each step from its own folder with the bash script it contains:

```bash
# for example
cd 0.illumination_correction
source perform_ic.sh
```

Each step's README describes its outputs.
