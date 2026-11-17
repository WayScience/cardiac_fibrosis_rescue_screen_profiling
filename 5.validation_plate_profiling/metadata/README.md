# Validation plate metadata

This folder holds the plate map shared by all five validation plates.

| File | What it is |
| --- | --- |
| [`platemaps/platemap_validation.csv`](./platemaps/platemap_validation.csv) | Plate map with one row per well: well position, heart number, cell type, heart failure type, and treatment |
| [`edit_metadata_excel_file.ipynb`](./edit_metadata_excel_file.ipynb) | Converts the plate map Excel file into the plate map CSV |
| [`platemap_vis.ipynb`](./platemap_vis.ipynb) | Draws the example plate map figure in [`platemap_fig/`](./platemap_fig/) |
| [`metadata_env.yml`](./metadata_env.yml) | Conda environment for the two notebooks |

## Plate layout

The plate has 60 wells: 15 conditions with 4 wells each.

- 13 compounds on failing cells (heart 25, dilated cardiomyopathy)
- DMSO on failing cells (heart 25), the negative control used for bulk normalization
- DMSO on non-failing cells (heart 7), the healthy control

![example_platemap](./platemap_fig/example_platemap.png)
