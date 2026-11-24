#!/usr/bin/env python
# coding: utf-8

# # Quantify the position effect
# 
# This notebook measures how large the position effect is on the validation plates and how much of it the correction removes, with four measurements.
# It changes no data: it reads the well medians that `1.fit_position_correction.ipynb` cached and the committed fit.
# 
# The DMSO wells are the backbone of the measurements, because they sit at the same positions on all five plates: 18 healthy and 18 failing positions, in blocks of columns that cover all six rows (failing in columns 2, 4, and 10, healthy in 3, 9, and 11).
# The variation of a DMSO position between the plates is noise, and the variation between the positions of one cell type is position.
# The compound wells sit on one plate each (three compounds per plate), so they cannot be replicated.
# 
# 1. **Position variation of every feature.** The share of the variation within each cell type of DMSO that is position, with the chance level subtracted (`vpu.position_variance_share`), on the uncorrected wells, after a tilt fit on all plates (in-sample, which a fit on its own data always favors), and after a tilt fit on the other four plates (held-out, `vpu.leave_one_plate_out_correction`).
# 2. **Position against biology.** The size of the healthy against failing DMSO difference on each plate, with and without the correction, next to the difference that the tilt predicts between the healthy and the failing positions.
# 3. **Does the tilt carry over between cell types?** The correction assumes position moves healthy and failing cells equally. Both cell types sit at many positions, so we fit the tilt on the failing DMSO only and check whether it predicts the position pattern of the healthy DMSO, and the other way round.
# 4. **Per plate.** The distance of every DMSO well from the mean of its cell type, on the layout of each plate, before and after the held-out correction.
# 
# Effect sizes are the length of a difference of profiles over all features, in standard deviations of a single cell.
# 
# **Output**
# 
# - `results/position_variance_by_feature.csv`, `position_variance_summary.csv`, `position_vs_biology.csv`, `position_vs_biology_plate_average.csv`, `tilt_cross_cell_type.csv`, `well_deviation_from_class_mean.csv`, and `well_deviation_by_row.csv`
# - `figures/quantify/`: one figure for each measurement

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append("../../../utils")
import position_correction_utils as pcu
import validation_position_utils as vpu


# ## Set paths and parameters

# In[2]:


figure_dir = pathlib.Path("figures") / "quantify"
figure_dir.mkdir(parents=True, exist_ok=True)
results_dir = pathlib.Path("results")
results_dir.mkdir(exist_ok=True)

fit_path = pathlib.Path("position_correction_fit.npz")
well_medians_cache = pathlib.Path("well_medians.parquet")
for path in (fit_path, well_medians_cache):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run 1.fit_position_correction.ipynb first.")

min_cells = 50  # minimum cells for a well to be used (as in the screen)
HEALTHY, FAILING = vpu.HEALTHY_CELL_TYPE, vpu.FAILING_CELL_TYPE
colors = {"raw": "#d62728", "in_sample": "#7f7f7f", "held_out": "#1f77b4"}


# ## Load the fit and the well medians

# In[3]:


fit = pcu.load_fit(fit_path)
well_medians, features = vpu.load_or_build_well_medians(well_medians_cache, ".")
features = [feature for feature in fit["features"] if feature in well_medians.columns]
feature_index = [list(fit["features"]).index(feature) for feature in features]
coef, center = fit["coef"][:, feature_index], fit["center"][feature_index]

wells = well_medians[well_medians["n"] >= min_cells].reset_index(drop=True)
rows, cols0 = pcu.well_index(wells["well"])
wells["row"], wells["column"] = rows, wells["well"].str[1:].astype(int)
raw = wells[features].to_numpy(float)
# the shift that the tilt predicts at every well (the correction subtracts it)
tilt_at_wells = pcu.tilt(coef, center, rows, cols0)
corrected = raw - tilt_at_wells
is_dmso = (wells["treatment"] == "DMSO").to_numpy()
is_healthy = (is_dmso & (wells["cell_type"] == HEALTHY)).to_numpy()
is_failing = (is_dmso & (wells["cell_type"] == FAILING)).to_numpy()
print("Wells:", len(wells), "| Features:", len(features), "| Plates:", wells["plate"].nunique())


# ## 1. Position variation of every feature
# 
# For each feature, the share of the variation within the two DMSO groups that is position.
# Zero is the chance level: no more variation between positions than the noise between plates predicts.

# In[4]:


shares = vpu.held_out_position_variance(well_medians, features, min_cells)
shares = pd.concat(
    [shares, pd.DataFrame([vpu.parse_feature_family(name) for name in shares["feature"]])],
    axis=1,
)
shares.to_csv(results_dir / "position_variance_by_feature.csv", index=False)
valid = shares.dropna(subset=["raw", "in_sample", "held_out"])

has_effect = valid["raw"] > 0.2
removed = 1 - valid.loc[has_effect, "held_out_excess"] / valid.loc[has_effect, "raw_excess"]
excess_left = {
    column: valid[f"{column}_excess"].sum() / valid["raw_excess"].sum()
    for column in ("raw", "in_sample", "held_out")
}
summary = pd.DataFrame(
    {
        "stage": ["raw", "in-sample correction", "held-out correction"],
        "median_share": [valid[c].median() for c in ("raw", "in_sample", "held_out")],
        "p90_share": [valid[c].quantile(0.9) for c in ("raw", "in_sample", "held_out")],
        "features_above_0.2": [(valid[c] > 0.2).mean() for c in ("raw", "in_sample", "held_out")],
        "position_variation_left": [excess_left[c] for c in ("raw", "in_sample", "held_out")],
    }
)
summary["features"] = len(valid)
summary.to_csv(results_dir / "position_variance_summary.csv", index=False)
print(summary.round(3).to_string(index=False))
print(
    f"\nFeatures with a clear position effect (share above 0.2): {int(has_effect.sum())}. "
    f"The held-out correction removes a median of {removed.median():.2f} of their position variation "
    f"(quartiles {removed.quantile(0.25):.2f} to {removed.quantile(0.75):.2f}); "
    f"{(removed < 0).mean():.0%} of these features get worse."
)


# In[5]:


fig, axes = plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw={"width_ratios": [1.1, 1, 1]})
labels = {"raw": "Uncorrected", "in_sample": "Corrected, tilt fit on all plates (in-sample)", "held_out": "Corrected, tilt fit on the other four plates (held-out)"}
ax = axes[0]
for column in ("raw", "in_sample", "held_out"):
    values = np.sort(valid[column].to_numpy())
    ax.plot(values, np.arange(1, len(values) + 1) / len(values), color=colors[column], label=labels[column], lw=2)
ax.axvline(0, color="black", lw=0.8, ls="--")
ax.set_xlabel("Share of the variation within the DMSO groups that is position")
ax.set_ylabel("Fraction of features")
ax.set_title("Position variation of the features")
ax.legend(fontsize=8, loc="lower right")
for ax, family, title in ((axes[1], "channel", "By channel"), (axes[2], "measurement", "By measurement type")):
    medians = valid.groupby(family)[["raw", "held_out"]].median()
    medians = medians[valid.groupby(family).size() >= 20].sort_values("raw")
    y = np.arange(len(medians))
    ax.barh(y + 0.2, medians["raw"], height=0.4, color=colors["raw"], label="Uncorrected")
    ax.barh(y - 0.2, medians["held_out"], height=0.4, color=colors["held_out"], label="Held-out correction")
    ax.set_yticks(y)
    ax.set_yticklabels(medians.index)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Median share that is position")
    ax.set_title(title)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(figure_dir / "position_variance_by_feature.png", dpi=150)
plt.show()


# ## 2. Position against biology
# 
# The healthy against failing DMSO difference on each plate (the length of the difference of the mean profiles of the two cell types), before and after the correction.
# The line is the difference that the tilt predicts between the positions of the healthy and the failing wells.
# If position were all of the difference, the corrected difference would be zero.

# In[6]:


records = []
for plate in sorted(wells["plate"].unique()):
    in_plate = (wells["plate"] == plate).to_numpy()
    h, f = in_plate & is_healthy, in_plate & is_failing
    records.append(
        {
            "plate": plate,
            "raw": np.linalg.norm(raw[h].mean(0) - raw[f].mean(0)),
            "corrected": np.linalg.norm(corrected[h].mean(0) - corrected[f].mean(0)),
        }
    )
position_vs_biology = pd.DataFrame(records)
position_vs_biology["share_kept"] = position_vs_biology["corrected"] / position_vs_biology["raw"]
tilt_healthy_vs_failing = float(np.linalg.norm(tilt_at_wells[is_healthy].mean(0) - tilt_at_wells[is_failing].mean(0)))
position_vs_biology["tilt_healthy_vs_failing"] = tilt_healthy_vs_failing
position_vs_biology.to_csv(results_dir / "position_vs_biology.csv", index=False)
print(position_vs_biology.round(2).to_string(index=False))

# the same difference for the profiles averaged over the five plates (less noise than a single plate)
pooled_raw = np.linalg.norm(raw[is_healthy].mean(0) - raw[is_failing].mean(0))
pooled_corrected = np.linalg.norm(corrected[is_healthy].mean(0) - corrected[is_failing].mean(0))
pd.DataFrame({"raw": [pooled_raw], "corrected": [pooled_corrected]}).to_csv(
    results_dir / "position_vs_biology_plate_average.csv", index=False
)
print(f"Plate-averaged healthy minus failing difference: {pooled_raw:.1f} before and {pooled_corrected:.1f} after the correction")


# In[7]:


fig, ax = plt.subplots(figsize=(8, 4.6))
x = np.arange(len(position_vs_biology))
ax.bar(x - 0.2, position_vs_biology["raw"], 0.4, color=colors["raw"], label="Uncorrected")
ax.bar(x + 0.2, position_vs_biology["corrected"], 0.4, color=colors["held_out"], label="Corrected")
ax.axhline(tilt_healthy_vs_failing, color="black", ls="--", label=f"Healthy against failing positions in the tilt ({tilt_healthy_vs_failing:.1f})")
ax.set_xticks(x)
ax.set_xticklabels([f"Plate {i + 1}" for i in x])
ax.set_ylabel("Length of the healthy minus failing DMSO difference")
ax.set_title("Healthy against failing DMSO, with and without the correction")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(figure_dir / "position_vs_biology.png", dpi=150)
plt.show()


# ## 3. Does the tilt carry over between cell types?
# 
# The correction subtracts one tilt from healthy and failing cells.
# That is right only if position moves both cell types in the same way.
# The two cell types sit at different positions (18 each), so we can test it: fit the tilt on the wells of one cell type only, predict the pattern across the positions of the other cell type, and compare it with the pattern that the wells show.
# 
# For each direction we report the cosine between the predicted and the observed pattern (1 is a perfect match), the slope of the observed on the predicted pattern (1 means the right size), and a ceiling for the cosine: the cosine between the patterns of two halves of the plates (plates 1, 3, 5 against plates 2, 4), which is as high as any prediction can be given the noise.

# In[8]:


dmso = wells[is_dmso]


def position_pattern(table):
    """Mean profile of every position (the plates averaged), centered over the positions."""
    means = table.groupby("well")[features].mean()
    return means.index, means.to_numpy(float) - means.to_numpy(float).mean(axis=0)


def cosine_and_slope(observed, predicted):
    """Cosine and slope between two patterns (positions by features)."""
    inner = float((observed * predicted).sum())
    return inner / np.sqrt((observed**2).sum() * (predicted**2).sum()), inner / float((predicted**2).sum())


plate_order = sorted(wells["plate"].unique())
records = []
for source, target in ((FAILING, HEALTHY), (HEALTHY, FAILING)):
    source_wells = dmso[dmso["cell_type"] == source]
    tilt_fit = vpu.fit_compound_adjusted_tilt(source_wells, features, min_cells, include_controls=True)
    target_wells = dmso[dmso["cell_type"] == target]
    positions, observed = position_pattern(target_wells)
    t_rows, t_cols = pcu.well_index(pd.Series(positions))
    predicted = pcu.tilt(tilt_fit["coef"], tilt_fit["center"], t_rows, t_cols)
    predicted = predicted - predicted.mean(axis=0)
    cosine, slope = cosine_and_slope(observed, predicted)
    half_a = target_wells[target_wells["plate"].isin(plate_order[0::2])]
    half_b = target_wells[target_wells["plate"].isin(plate_order[1::2])]
    _, pattern_a = position_pattern(half_a)
    _, pattern_b = position_pattern(half_b)
    ceiling, _ = cosine_and_slope(pattern_a, pattern_b)
    records.append({"fit_on": source, "predicts": target, "cosine": cosine, "slope": slope, "ceiling_cosine": ceiling, "n_positions": len(positions)})
cross = pd.DataFrame(records)
cross.to_csv(results_dir / "tilt_cross_cell_type.csv", index=False)
print(cross.round(2).to_string(index=False))


# In[9]:


fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
x = np.arange(len(cross))
ax.bar(x - 0.2, cross["cosine"], 0.4, color=colors["held_out"], label="Cosine of the prediction")
ax.bar(x + 0.2, cross["ceiling_cosine"], 0.4, color="#bdbdbd", label="Ceiling (two halves of the plates)")
ax.set_xticks(x)
ax.set_xticklabels([f"Fit on {a}\\npredicts {b}" for a, b in zip(cross["fit_on"], cross["predicts"])])
ax.set_ylim(0, 1.05)
ax.set_ylabel("Cosine between predicted and observed position pattern")
ax.set_title("Does the tilt carry over between cell types?")
ax.legend(fontsize=9, loc="lower right")
ax = axes[1]
source_wells = dmso[dmso["cell_type"] == FAILING]
tilt_fit = vpu.fit_compound_adjusted_tilt(source_wells, features, min_cells, include_controls=True)
positions, observed = position_pattern(dmso[dmso["cell_type"] == HEALTHY])
t_rows, t_cols = pcu.well_index(pd.Series(positions))
predicted = pcu.tilt(tilt_fit["coef"], tilt_fit["center"], t_rows, t_cols)
predicted = predicted - predicted.mean(axis=0)
direction = np.linalg.svd(predicted, full_matrices=False)[2][0]
row_colors = dict(zip(pcu.ROWS, plt.get_cmap("viridis")(np.linspace(1, 0, pcu.N_ROW_TERMS))))
for position, px, py in zip(positions, predicted @ direction, observed @ direction):
    ax.scatter(px, py, color=row_colors[position[0]], s=45)
limit = max(abs(predicted @ direction).max(), abs(observed @ direction).max()) * 1.1
ax.plot([-limit, limit], [-limit, limit], color="black", lw=0.8, ls="--")
ax.set_xlabel("Predicted by the tilt fit on the failing DMSO")
ax.set_ylabel("Observed in the healthy DMSO")
ax.set_title("Healthy DMSO positions on the main direction of the tilt (color: plate row B to G)")
ax.title.set_size(9)
fig.tight_layout()
fig.savefig(figure_dir / "tilt_cross_cell_type.png", dpi=150)
plt.show()


# ## 4. Per plate: distance of every DMSO well from its cell type
# 
# For every DMSO well, the distance between its profile and the mean profile of its cell type over all positions and plates.
# Position moves the wells of a cell type apart, so before the correction the wells at the edges of the plate sit far from the mean.
# After the held-out correction (the tilt is fit on the other four plates), what is left is noise and whatever the correction misses.
# The compound wells are not part of this measurement and are shown in grey.

# In[10]:


held_out_profiles = vpu.leave_one_plate_out_correction(wells, features, min_cells)


def distance_from_class_mean(values):
    """Distance of each DMSO well from the mean profile of its cell type."""
    result = np.full(len(wells), np.nan)
    for mask in (is_healthy, is_failing):
        result[mask] = np.linalg.norm(values[mask] - values[mask].mean(axis=0), axis=1)
    return result


deviation = wells[["plate", "well", "treatment", "cell_type"]].copy()
deviation["row"] = deviation["well"].str[0]
deviation["column"] = deviation["well"].str[1:].astype(int)
deviation["raw"] = distance_from_class_mean(raw)
deviation["held_out"] = distance_from_class_mean(held_out_profiles)
deviation = deviation[is_dmso]
deviation.to_csv(results_dir / "well_deviation_from_class_mean.csv", index=False)

by_row = deviation.groupby("row")[["raw", "held_out"]].mean()
by_row["share_kept"] = by_row["held_out"] / by_row["raw"]
by_row.to_csv(results_dir / "well_deviation_by_row.csv")
print("Mean distance from the mean of the cell type, by plate row")
print(by_row.round(2).to_string())


# In[11]:


plates = sorted(deviation["plate"].unique())
rows_label = list(pcu.ROWS)
column_numbers = list(range(2, 12))


def layout(plate, column_name):
    grid = np.full((len(rows_label), len(column_numbers)), np.nan)
    for _, well in deviation[deviation["plate"] == plate].iterrows():
        grid[rows_label.index(well["row"]), well["column"] - 2] = well[column_name]
    return grid


vmax = max(deviation["raw"].max(), deviation["held_out"].max())
fig, axes = plt.subplots(3, len(plates), figsize=(3.4 * len(plates), 8.6), squeeze=False)
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad("#dddddd")
diverging = plt.get_cmap("RdBu_r").copy()
diverging.set_bad("#dddddd")
for j, plate in enumerate(plates):
    raw_grid, held_grid = layout(plate, "raw"), layout(plate, "held_out")
    panels = [
        (raw_grid, cmap, 0, vmax, "Uncorrected"),
        (held_grid, cmap, 0, vmax, "Corrected (held-out)"),
        (raw_grid - held_grid, diverging, -vmax / 2, vmax / 2, "Uncorrected minus corrected"),
    ]
    for i, (grid, colormap, low, high, title) in enumerate(panels):
        ax = axes[i, j]
        image = ax.imshow(grid, cmap=colormap, vmin=low, vmax=high, aspect="equal")
        ax.set_xticks(range(len(column_numbers)))
        ax.set_xticklabels(column_numbers, fontsize=7)
        ax.set_yticks(range(len(rows_label)))
        ax.set_yticklabels(rows_label, fontsize=8)
        ax.set_title(f"Plate {j + 1}\n{title}" if i == 0 else title, fontsize=9)
        if j == len(plates) - 1:
            colorbar = fig.colorbar(image, ax=axes[i, :], shrink=0.8, pad=0.02)
            colorbar.set_label("Change" if i == 2 else "Distance", fontsize=8)
fig.suptitle("Distance of each DMSO well from the mean profile of its cell type, before and after the correction")
fig.savefig(figure_dir / "plate_layout_well_deviation.png", dpi=150, bbox_inches="tight")
plt.show()

