# Plate-position correction of the validation plates (step 4a)

In the screen, the position of a well on the plate shifts its profile in a consistent, feature-specific direction (a "tilt") that has nothing to do with what is in the well.
The screen estimates the tilt, validates it on held-out platemaps, and subtracts it ([`4.image_based_profiling/3a.position_correction/`](../../../4.image_based_profiling/3a.position_correction/)).

This folder does the same for the five validation plates, before we call hits on them.
It fits the correction, validates it on held-out plates, and draws the same figures as the screen.
[`../4b.apply_position_correction.ipynb`](../4b.apply_position_correction.ipynb) applies the fit to the single cells (the step 3b of the screen), and [`../5.bulk_position_corrected.ipynb`](../5.bulk_position_corrected.ipynb) aggregates the corrected cells to bulk profiles (the bulk step of the screen).
Like in the screen, the pipeline script does not run this step: it reads the committed fit, `position_correction_fit.npz`.

The validation plates keep their own feature names (the d2 channel is `Golgi`), and the correction is fit on the validation plates only.

## The validation plates

All five replicate plates have 60 wells (rows B to G, columns 2 to 11):

| Wells | Cells | Where |
|---|---|---|
| 18 failing DMSO (negative control) | failing (heart 25) | columns 2, 4, and 10 (all six rows), identical positions on all plates |
| 18 healthy DMSO | non-failing (heart 7) | columns 3, 9, and 11 (all six rows), identical positions on all plates |
| 24 compound wells | failing (heart 25) | columns 5 to 8 (all six rows) |

There are 15 compounds (K020, K191, K208, K221, K278, K310, K318, K463, K504, K543, K645, K687, K756, K779, K845).
Each plate holds three of them in 8 wells each, and each compound is on one plate only, at positions that differ between plates.
Plate 1, for example, holds K278, K310, and K543.

This layout differs from the screen in three ways that matter for the correction:

- The screen places each compound once on each of 11 platemaps and the same controls everywhere, so compound effects average out of its tilt estimate.
  Here a compound is on one plate, so its effect has to be fit together with the tilt (one effect for each compound, `fit_compound_adjusted_tilt` in [`utils/validation_position_utils.py`](../../../utils/validation_position_utils.py)).
- The DMSO wells of both cell types are the only wells that repeat at the same position on five plates.
  Within a cell type, the difference between positions is position, and the difference between plates is noise.
  So the fit uses the DMSO wells too (`include_controls=True`), and the measurements of this folder that need replicates use only the DMSO wells.
- **The groups sit in blocks of columns.** Every row has the same number of wells of each kind (three failing DMSO and three healthy DMSO), so row is not confounded with the cell type, but the failing DMSO, the healthy DMSO, and the compounds each have their own columns.
  A group effect can therefore absorb the offset of its block of columns, and three column offsets cannot be told from the group effects (on one plate the design has 20 columns and rank 17).
  An unconstrained fit splits those offsets between group and column effects arbitrarily, which changes the difference between the healthy and the failing controls.
  The fit has no reference column and a tiny ridge penalty on the column effects, which sets what the data cannot determine to zero and leaves the effects that they do determine as they are.

## Notebooks

| Notebook | What it does | Screen counterpart |
|---|---|---|
| [`0.diagnose_position_effect.ipynb`](./0.diagnose_position_effect.ipynb) | Diagnoses the uncorrected plates: how much of each principal component plate row and column explain (in all wells, in the DMSO wells of one cell type, and in the DMSO wells without the cell-type effect), and whether the effect works through the number of cells | none (it motivates the correction) |
| [`1.fit_position_correction.ipynb`](./1.fit_position_correction.ipynb) | Estimates the tilt (with one effect for each compound and each kind of DMSO), writes `position_correction_fit.npz`, and measures the amplitude of each replicate plate and the effect on the controls | `0.fit_position_correction` |
| [`2.validate_position_correction.ipynb`](./2.validate_position_correction.ipynb) | Tests the correction on held-out plates with a row classifier, and checks on a simulation that the fit recovers a known tilt | `1.validate_position_correction` |
| [`3.effect_on_feature_space.ipynb`](./3.effect_on_feature_space.ipynb) | Computes the size of the tilt by feature family and the principal components of the wells before and after the correction | `2.effect_on_feature_space` |
| [`4.plot_effect_on_feature_space.ipynb`](./4.plot_effect_on_feature_space.ipynb) (R) | Draws the effect figures from those tables | `3.plot_effect_on_feature_space` |
| [`5.quantify_position_effect.ipynb`](./5.quantify_position_effect.ipynb) | Measures the position effect directly: the position variation of every feature within the DMSO groups, position against biology, whether the tilt carries over between cell types, and the DMSO wells on the layout of each plate, before and after the correction | none |
| [`6.compare_position_effects_in_validation_to_screen.ipynb`](./6.compare_position_effects_in_validation_to_screen.ipynb) | Puts the validation plates next to the screen: the variance of the principal components split into position and heart status, the survival of the healthy against failing difference, and the placement of compounds on the screen platemaps | none |

Run them in order; every notebook has a script version in [`nbconverted/`](./nbconverted/).
The first one that needs the single-cell profiles reads the five `<plate>_sc_normalized.parquet` files (about 1.4 GB, on bandicoot at `~/mnt/bandicoot/CFReT_screening_data/validation_profiles/single_cell_profiles/`, or in the folder that `VALIDATION_PROFILE_DIR` names) and caches the median profile of every well in `well_medians.parquet`, which git ignores.
Delete the cache whenever the single-cell profiles change.
The R notebook needs R with `ggplot2`, `dplyr`, `readr`, and `cowplot`.

Differences from the screen:

| | Screen | Validation plates |
|---|---|---|
| Tilt | one shared tilt, fit on compounds placed at random on 11 platemaps | one tilt fit on the compound wells and the DMSO wells of both cell types, with one effect for each compound and each kind of DMSO, and a ridge penalty on the column effects (blocks of columns cannot be told from the group effects) |
| Amplitude | one per platemap, shrunk toward the others | one per replicate plate, shrunk toward the others like the platemaps of the screen: 1.05, 0.92, 0.92, 0.95, and 0.82 for plates 1 to 5 (before shrinking: 1.10, 0.92, 0.92, 0.96, and 0.71) |
| Held-out validation | leave one platemap out | leave one plate out (the tilt is fit on four plates and applied to the fifth) |
| Figures of the principal components | the two components that plate row explains most | PC1, which plate row explains most and which is also the largest, and PC8, the next that plate row explains most |

## Figures

The figures are the five figures of the screen, for the validation plates.

| Figure | What it shows |
|---|---|
| [`figures/fit/tilt_size_and_amplitudes.png`](./figures/fit/tilt_size_and_amplitudes.png) | The size of the tilt at every well position, and the amplitude of the tilt in each replicate plate |
| [`figures/validation/row_prediction_validation.png`](./figures/validation/row_prediction_validation.png) | Held-out row prediction before and after the correction, per plate and pooled, against shuffled labels |
| [`figures/effect/tilt_size_by_feature_family.png`](./figures/effect/tilt_size_by_feature_family.png) | Which features the position moves most, by compartment, measurement type, and channel |
| [`figures/effect/well_pca_before_after_by_row.png`](./figures/effect/well_pca_before_after_by_row.png) | The wells in the principal components before and after the correction, with the displacement of the wells |
| [`figures/effect/plate_layout_top_row_pc_before_after.png`](./figures/effect/plate_layout_top_row_pc_before_after.png) | The component that plate row explains most on the plate layout, before and after, for each replicate plate |

[`figures/diagnostics/`](./figures/diagnostics/) has the figures of the diagnostics notebook. [`results/`](./results/) has the tables that later notebooks read (the fit, the effect on the feature space, and the row prediction); the other numbers are printed in the notebooks.

[`figures/quantify/`](./figures/quantify/) has the figures of the notebook that quantifies the position effect:

| Figure | What it shows |
|---|---|
| [`position_variance_by_feature.png`](./figures/quantify/position_variance_by_feature.png) | The share of the variation within the DMSO groups that is position, for every feature, uncorrected and after the in-sample and held-out corrections, and by channel and measurement type |
| [`position_vs_biology.png`](./figures/quantify/position_vs_biology.png) | The healthy against failing DMSO difference on each plate, with and without the correction, next to the difference that the tilt predicts between the healthy and the failing positions |
| [`tilt_cross_cell_type.png`](./figures/quantify/tilt_cross_cell_type.png) | Whether a tilt fit on the DMSO of one cell type predicts the position pattern of the other cell type |
| [`plate_layout_well_deviation.png`](./figures/quantify/plate_layout_well_deviation.png) | For each plate, the distance of every DMSO well from the mean profile of its cell type, before and after the held-out correction |
| [`tilt_terms_and_uncertainty.png`](./figures/quantify/tilt_terms_and_uncertainty.png) | The size of each of the 14 row and column terms of the tilt, with its uncertainty from leaving out each plate |

[`figures/compare_to_screen/`](./figures/compare_to_screen/) has the figures of `6.compare_position_effects_in_validation_to_screen.ipynb` (see Comparison with the screen below):

| Figure | What it shows |
|---|---|
| [`scree_position_and_heart_status.png`](./figures/compare_to_screen/scree_position_and_heart_status.png) | The variance of the first 20 principal components of the wells, with the part that plate position explains (red) and the part that heart status explains (blue), before and after the correction, for the screen (all wells, as in its own analysis) and the validation plates, on one y axis |
| [`healthy_failing_survival.png`](./figures/compare_to_screen/healthy_failing_survival.png) | The length of the healthy against failing DMSO difference on every plate before and after the correction, and the share that survives for every screen platemap and for the validation plates |
| [`screen_compound_placement.png`](./figures/compare_to_screen/screen_compound_placement.png) | The correlation between compound number and position on the plate for every screen platemap |

## What we found

**1. The validation plates have a strong position effect.**

- The first principal component of the wells (48% of the variance) is 53% explained by plate row and 20% by plate column (chance is about 2% and 3%).
  In the screen, plate row explained 25% of PC3 and 10% of PC2, the components that depend most on it, and 1% of PC1.
- It is not the effect of the cell type or of the compounds.
  Among the failing DMSO wells alone, which are one treatment on one cell type, PC1 (54% of their variance) is 69% explained by plate row and 7% by plate column.
  Among the healthy DMSO wells alone, PC1 (57%) is 76% explained by row and 4% by column.
  After we average the five plates and remove the mean of each cell type, what remains differs only between the 36 positions, and PC1 of that (79%) is still 84% explained by row and 6% by column.
- The effect is strongest in Actin features (median tilt size 0.34 standard deviations of a single cell, against 0.09 to 0.16 in the other channels) and in intensity and texture measurements, as in the screen.

**2. It goes partly through the number of cells.**
The median number of cells in a DMSO well falls from 257 in row B and 260 in row C to 226 in row F and 203 in row G.
Within a cell type, the cell number correlates with the row (-0.63 after the mean of the cell type is removed), and the cell type explains none of the variation in cell number.
In the failing DMSO wells, plate row explains 69% of PC1 and the log number of cells 5%; together 70%.
In the healthy DMSO wells row explains 76%, cell number 33%, and together 81%.
A correction by position removes the part that goes through position, and we do not model the cell number separately.

**3. The five plates do not show the effect equally.**
The amplitude of the tilt is 1.10 (plate 1), 0.92 (plate 2), 0.92 (plate 3), 0.96 (plate 4), and 0.71 (plate 5), with standard errors of 0.05 to 0.07.
Plate 5 is about five standard errors away from 1, and plate 1 about two, so one tilt for all plates is an approximation.
Like the screen, the correction uses each plate's amplitude shrunk toward the average of the others (1.05, 0.92, 0.92, 0.95, and 0.82).

**4. The correction works on held-out plates, but only partly.**
The table gives the leave-one-plate-out result: the AUC of a classifier that predicts the plate row of the wells of a plate that neither the tilt nor the classifier saw (0.5 is chance; shuffled labels, repeated several times, give the floor).

| | raw | corrected | shuffled labels |
|---|---|---|---|
| All wells: AUC over the rows | 0.85 | 0.78 | 0.50 (sd 0.02) |
| All wells: row B against the rest | 0.98 | 0.74 | 0.45 (sd 0.09) |
| DMSO wells: AUC over the rows | 0.87 | 0.89 | 0.51 (sd 0.02) |
| DMSO wells: row B against the rest | 0.99 | 0.85 | 0.44 (sd 0.09) |

- The correction removes 20% of the excess AUC over the rows in all wells and 46% of the excess for row B against the rest.
  On the DMSO wells alone the AUC over the rows does not fall.
  The DMSO positions are the same on every plate, so a classifier can learn the pattern of a well from the other plates, which a row plus column tilt does not contain; we have not tested whether this is the reason.
  In the screen, the same test took the pooled AUC from 0.76 to 0.52.
- On each plate separately, the AUC over the rows goes from 0.88 to 0.69 (plate 1), 0.88 to 0.79 (plate 2), 0.89 to 0.80 (plate 3), 0.84 to 0.85 (plate 4), and 0.81 to 0.81 (plate 5).
- On the plates themselves (in-sample), the correction leaves 12.6% of the total variance of the first 20 principal components as position, against 43.2% before.
  What remains is mostly column structure that is confounded with the groups (finding 5 and `3.effect_on_feature_space.ipynb`: after the correction, plate column explains 41% of PC1).
- In a simulation on the real layout, the fit with one effect for each group recovers a known tilt (median cosine 0.95), and the screen's way of fitting, without group effects, reaches 0.78.

**5. The correction does not change the healthy versus failing contrast.**
The length of the healthy minus failing DMSO difference is 12.8 to 14.5 on every plate before the correction and the same after it (the share that survives is 1.00), and the tilt itself predicts no difference between the positions of the healthy and the failing wells.
This is by design: the rows are balanced between the two cell types, and the fit does not split the column-block offsets between the group and the column effects.
In the screen the median over plates is 0.84, but it hides large differences between platemaps: the correction makes the difference 1.5 to 2.1 times longer on platemaps 1 to 4, and keeps 0.40 to 0.85 of it on platemaps 5 to 11 (see [`6.compare_position_effects_in_validation_to_screen.ipynb`](./6.compare_position_effects_in_validation_to_screen.ipynb)).

**6. The tilt of one cell type predicts the position pattern of the other well.**
Both cell types are at 18 positions, so we fitted the tilt on the failing DMSO only and compared its prediction for the 18 healthy positions with what the healthy wells show (and the reverse).
The cosine between the predicted and the observed pattern is 0.80 (fit on failing, predicts healthy) and 0.78 (the reverse), while two halves of the plates of the same cell type agree at 0.82 and 0.80, which is as high as any prediction can be at this noise.
The prediction for the healthy wells has about the right size (slope 1.09), and the prediction for the failing wells has 0.78 of it.
So position moves the two cell types in nearly the same way, and one shared tilt is a good approximation.

**7. Some terms of the tilt are poorly determined.**
Leaving out each plate in turn, the five terms of the rows C to G against row B are 5 to 7 times larger than their uncertainty.
None of the nine column terms is more than 2 times its uncertainty (1.0 to 2.0), so the column effects are not determined by these plates, in line with the fact that they are confounded with the groups and do not carry over between plates.

## Measuring the position effect directly

[`5.quantify_position_effect.ipynb`](./5.quantify_position_effect.ipynb) measures the position effect with four measurements that need no classifier.
The biology is removed by looking only at the DMSO wells of one cell type: they sit at the same 18 positions on each of the five plates, so the variation of a position mean is position, and the variation between plates at the same position is noise.

**1. Position variation of the features.**
Of the variation within the DMSO groups, the share that is position (with the chance level subtracted) has a median of 0.26 per feature before the correction, and 55% of the features are above 0.2.
After a tilt fit on all plates it is 0.10 (36% above 0.2), which a fit on its own data always favors.
After a tilt that never saw the plate (held-out) it is 0.05 (28% above 0.2).
For a median feature with a clear position effect (share above 0.2, 1,335 features), the held-out correction removes 85% of the position variation (quartiles 67% to 92%), and 6% of these features get worse.
Overall 14% of the position variation is left after the in-sample correction and 8% after the held-out one.
The held-out figure being lower than the in-sample one is not credible as a real gain: the measure subtracts a noise level estimated from the differences between plates, and the five held-out corrections differ a little between plates, which adds to that noise.
We take the two figures to be about equal, and the remaining position variation to be 8% to 14%.

**2. Position against biology.**
The length of the healthy against failing DMSO difference is 12.8 to 14.5 on the plates, before and after the correction (the share kept is 1.00), and the tilt itself differs by 0 between the healthy and the failing positions.

**3. Does the tilt carry over between cell types?**
See finding 6 above (`tilt_cross_cell_type.png`).

**4. Per plate.**
The mean distance of a DMSO well from the mean profile of its cell type is 19.4 in row B and 9.8 to 12.0 in the other rows before the correction, so the wells at the top edge of the plate sit furthest from the mean of their cell type.
After the held-out correction it is 12.6 in row B and 9.2 to 11.2 in the others (row B keeps 65% of its distance, rows C to G 81% to 96%).
The corrected distance is close to the noise between replicate plates, so only row B stands out, and it still does.
The same pattern holds on all five plates ([`plate_layout_well_deviation.png`](./figures/quantify/plate_layout_well_deviation.png)).

## Comparison with the screen

[`6.compare_position_effects_in_validation_to_screen.ipynb`](./6.compare_position_effects_in_validation_to_screen.ipynb) runs the same measurements on the screen and on the validation plates.

**What the tilt is.** The correction subtracts from every feature a position effect, the sum of one effect for the row and one for the column of a well.
The fit of a feature has one effect for each group (15 compounds and 2 kinds of DMSO), five row effects (rows C to G against row B), and ten column effects with a ridge penalty (reported against column 2): 32 parameters, of which three combinations are not determined by the layout.
It is fit on the 17 groups, with the replicate plates averaged (36 DMSO positions and 15 compounds at 8 positions each), once for each of the 2,453 features.
The group effects keep what a compound or a cell type does from being counted as position.

**Variance of the principal components.**
We fit the components as the screen analysis does (for the screen, on all 2,088 wells) and show the first 20, on one y axis.
In the screen, plate position (row plus column) is in PC2 and PC3: it explains 10% of PC2 and 27% of PC3 (row alone: 10% and 25%), and 4% of PC1.
In all, it explains 6.8% of the total variance in the first 20 principal components of the screen wells (85.7% of the total variance), and 1.4% after the correction.
On the validation plates it is in PC1, which holds 47.9% of the total variance, and position explains 34.4 points of it (21.1% of the total variance remain in PC1 after the correction, 8.1 points of it position).
In all, position explains 43.2% of the total variance in the first 20 components of the validation wells (94.6% of the total variance), and 12.6% after the correction.
Heart status (healthy against failing) explains 4.1% of the total variance before and after the correction on the validation plates, against 1.3% and 0.7% in the screen.
The 12.6% after the correction is in-sample, for the plates that the tilt was fit on, so it is an optimistic figure; the held-out figure is in `5.quantify_position_effect.ipynb`.

**What the tilt estimates well.**
The jackknife (leaving out each plate in turn) gives an uncertainty for each of the 14 position terms (finding 7).
The rows are well determined and the columns are not.

**What generalizes.**
For a tilt fit on all plates, 36% of the features keep a position share above 0.2, and for a tilt fit without the plate, 28% do (55% before the correction; see the caveat in the measurement above).
Of the 1,335 features with a clear position effect (a share above 0.2 before the correction), 113 (8%) keep at least 80% of it after the held-out correction: 48 DNA, 29 without a channel, 13 ER, 12 Golgi, 5 multi-channel, 5 Mito, and 1 Actin feature.

**Compound placement in the screen.**
The number of the compound identifier rises with the position on the plate on platemaps 5, 7, 8, and 10 (Spearman correlation 0.49 to 0.71, p < 0.001), and does not on platemaps 1, 2, 3, 6, 9, and 11.
The pathways do not differ by plate row (chi-square p = 0.62).
So compounds were not placed independently of position on some platemaps, which the tilt of the screen assumes.

## Assumptions and limits

- Position adds the same fixed vector to every cell of a well, the same for every compound and cell type (additive, no interaction of row and column, no effect that differs between compounds or cell types).
  Finding 6 shows that the last part is only approximately true between the two cell types.
- The five replicate plates share the DMSO layout, so they reduce noise, but the compounds are on one plate each.
  The compound effects are fit, not replicated, and the 8 positions of a compound (columns 5 to 8, in all six rows) are what separate its effect from position.
- The groups sit in blocks of columns (failing DMSO in columns 2, 4, and 10, healthy DMSO in 3, 9, and 11, compounds in 5 to 8), so a column effect is confounded with the cell type and the compounds, and only the differences between columns inside a block are determined.
  The column terms of the tilt are therefore not reliable, and the row terms carry the correction.
- The held-out test is held out by plate, so it tells how the correction works for a new plate on this layout, not for a new layout.
- We use wells with at least 50 cells, and we drop features with a missing value in any well.
- The principal-component figures fit the components on the uncorrected wells and project the corrected wells into them, as in the screen.

## Contents

```text
4a.position_correction/
├── 0.diagnose_position_effect.ipynb                          Diagnose the uncorrected plates
├── 1.fit_position_correction.ipynb                           Fit the tilt
├── 2.validate_position_correction.ipynb                      Test it on held-out plates
├── 3.effect_on_feature_space.ipynb                           Tables of the effect on the features and the wells
├── 4.plot_effect_on_feature_space.ipynb                      Figures of the effect (R)
├── 5.quantify_position_effect.ipynb                          Measure the position effect directly
├── 6.compare_position_effects_in_validation_to_screen.ipynb  Compare the validation plates with the screen
├── nbconverted/                                              Script version of each notebook
├── position_correction_fit.npz                               The tilt map and amplitude (committed; step 4b reads it)
├── results/                                                  Tables that later notebooks read
├── figures/                                                  Figures (diagnostics/, fit/, validation/, effect/, quantify/, compare_to_screen/)
└── well_medians.parquet                                      Median profile of every well (git ignores this cache)
```
