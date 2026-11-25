# Validation plate profiling

This module processes the validation plates with the same workflow as modules 2, 3, and 4, using one plate map per plate in [`metadata/platemaps/`](./metadata/platemaps/).
There are five validation plates.
They share the same well layout: 36 DMSO control wells (18 non-failing and 18 failing) in blocks of columns that cover all six rows, and 24 compound wells in the other columns.
Each plate tests three different compounds in 8 wells each, and each compound is on one plate only:

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

## Plate maps and design

Each plate has one plate map, [`metadata/platemaps/<barcode>_platemap.csv`](./metadata/platemaps/), with one row per well (plate barcode and number, well position, heart number, cell type, heart failure type, and treatment).
All five plates have 60 wells (rows B to G, columns 2 to 11) and the same layout, in blocks of columns:

| Columns | Wells | Cells |
| --- | --- | --- |
| 2, 4, 10 | DMSO, 18 wells (3 per row) | failing: heart 25, dilated cardiomyopathy |
| 3, 9, 11 | DMSO, 18 wells (3 per row) | non-failing: heart 7 |
| 5 to 8 | compounds, 24 wells | failing: heart 25 |

The three compounds of a plate (A, B, and C in the table above) rotate through the 24 compound wells in the same pattern on every plate.
Each compound has 8 wells, which cover all four compound columns and all six rows (one or two wells per row).
Each compound is on one plate only, so its replicates are the 8 wells of that plate, and a plate has no replicate of its compounds on another plate.

What the design means for the analysis:

- **Controls for the hit calling.** The failing DMSO is the reference that compounds should move away from, and the non-failing DMSO is the state that they should move toward. The failing DMSO is also the negative control of the bulk normalization.
- **Rows are not confounded with cell type.** Both kinds of DMSO have three wells in every row, so a row effect can be estimated from the controls and does not look like a difference between non-failing and failing cells.
- **Identical layout on every plate.** The plates differ in their compounds only, so we can fit one position effect (row and column) to all plates and give each plate its own strength.
- **Columns are not covered by controls.** The compounds sit in columns 5 to 8, which have no DMSO, so the column effects of these columns come from the compound wells and are only partly identifiable. [`2.preprocessing_features/4a.position_correction`](./2.preprocessing_features/4a.position_correction/) explains how the correction handles this and how much of the position effect it removes.
- **Slots.** Each compound occupies the wells of one of three fixed sets of eight compound wells (slots A, B, and C) that are the same on every plate, and the compounds were dispensed into the slots at random. The compounds of slot A score lower than those of the other two on every plate, and the position correction cannot remove this: see [`2.preprocessing_features/4c.slot_a_adjustment`](./2.preprocessing_features/4c.slot_a_adjustment/README.md).

The plate maps enter the workflow in step 2 (the profiles get the treatment, cell type, and heart as metadata), in the position correction (the well row and column, and the group of each well), and in the [hit calling](https://github.com/WayScience/cardiac_fibrosis_rescue_screen_hit_calling) of the compounds (module `4.validation`).

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
| 2 | [`2.preprocessing_features`](./2.preprocessing_features/) | Convert to parquet, run single-cell QC, process single cells, correct the plate-position effect, and aggregate to bulk profiles; the optional step 4c writes separate profiles with the slot A adjustment |

[`metadata/`](./metadata/) holds the validation plate map and its visualization.

Run each step from its own folder with the bash script it contains:

```bash
# for example
cd 0.illumination_correction
source perform_ic.sh
```

Each step's README describes its outputs.
