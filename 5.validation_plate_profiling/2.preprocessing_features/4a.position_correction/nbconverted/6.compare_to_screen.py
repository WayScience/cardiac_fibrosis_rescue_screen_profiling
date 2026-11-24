#!/usr/bin/env python
# coding: utf-8

# # Compare the position effect on the validation plates with the screen
# 
# This notebook puts the validation plates next to the screen, with the same measurements on both, to show why the plate-position correction removed most of the position effect in the screen and only part of it on the validation plates.
# It changes no data: it reads the well medians and the fit of this step (`1.fit_position_correction.ipynb`) and of the screen (`4.image_based_profiling/3a.position_correction`).
# 
# 1. **Variance by principal component, split into position and heart status.** The principal components of the well profiles, with the part of each component that plate position (row plus column) explains in red and the part that the heart status of the cells (healthy against failing) explains in blue, before and after the correction.
# 2. **Survival of the healthy against failing difference.** The length of the difference between the healthy and the failing DMSO profile on every plate, before and after the correction.
# 3. **Compound placement in the screen.** Whether the position of a compound on the screen platemaps is independent of the compound.
# 
# **Input**
# 
# - `well_medians.parquet` and `position_correction_fit.npz` of this folder (validation plates)
# - `../../../4.image_based_profiling/3a.position_correction/well_medians.parquet` and `position_correction_fit.npz` (screen; the cache of the well medians is git-ignored, so run the screen's fit notebook first when it is missing)
# - `../../../metadata/updated_platemaps/` (screen platemaps)
# 
# **Output**: the three figures in `figures/compare_to_screen/`, which is the only folder with screen data. The tables behind them are printed in the notebook.

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, spearmanr

sys.path.append("../../../utils")
import position_correction_effect_utils as pce
import position_correction_utils as pcu
import validation_position_utils as vpu


# ## Set paths and parameters

# In[2]:


figure_dir = pathlib.Path("figures") / "compare_to_screen"
figure_dir.mkdir(parents=True, exist_ok=True)

screen_dir = pathlib.Path("../../../4.image_based_profiling/3a.position_correction")
platemap_dir = pathlib.Path("../../../metadata/updated_platemaps")
for path in (
    pathlib.Path("well_medians.parquet"),
    pathlib.Path("position_correction_fit.npz"),
    screen_dir / "well_medians.parquet",
    screen_dir / "position_correction_fit.npz",
):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run the fit notebooks of the screen and of this step first.")

min_cells = 50  # minimum cells for a well to be used (as in the screen)
n_components = 20  # principal components to show
seed = 0

# the two data sets, with the label of the healthy cells in each
datasets = {
    "Screen": {
        "medians": pd.read_parquet(screen_dir / "well_medians.parquet"),
        "fit": pcu.load_fit(screen_dir / "position_correction_fit.npz"),
        "healthy": "healthy",
    },
    "Validation": {
        "medians": pd.read_parquet("well_medians.parquet"),
        "fit": pcu.load_fit("position_correction_fit.npz"),
        "healthy": vpu.HEALTHY_CELL_TYPE,
    },
}
colors = {"position": "#d62728", "status": "#1f77b4", "total": "#bdbdbd"}
for name, data in datasets.items():
    # the screen analysis of the principal components uses every well, whatever its number of cells
    data["all_wells"] = data["medians"].reset_index(drop=True)
    data["medians"] = data["medians"][data["medians"]["n"] >= min_cells].reset_index(drop=True)
    print(name, "| wells:", len(data["medians"]), "| plates:", data["medians"]["plate"].nunique(), "| features:", len(data["fit"]["features"]))


# ## 1. Variance by principal component: position and heart status
# 
# We fit the principal components on the uncorrected wells (controls included), as in `3.effect_on_feature_space.ipynb`, and project the corrected wells on the same components.
# For every component:
# 
# - the **grey bar** is the variance of the wells on the component, as a percentage of the total variance of the uncorrected wells;
# - the **red bar** is the part of it that plate position explains, from an additive row plus column regression (adjusted R squared times the variance);
# - the **blue bar** is the part that heart status explains (healthy against failing cells, R squared times the variance).
# 
# Heart status and position are partly confounded.
# On the screen, the healthy wells are all in row B and the failing control wells in row E, so the blue bar before the correction also contains position.
# On the validation plates the two cell types have the same number of wells in each of the six rows, so row does not confound heart status, but they sit in different blocks of columns (failing in columns 2, 4, and 10, healthy in 3, 9, and 11), so column does.
# The validation fit uses the DMSO wells of both cell types, and what remains in the blue bar afterwards is the difference between healthy and failing cells that the position terms do not explain.
# The red bar after the correction is the position variation that the fit leaves in the wells that it was fit on, which a fit on its own data always favors.
# The position variation that the correction leaves in wells of a plate that it never saw is in section 3.

# In[3]:


# the principal components of the screen are fit on all of its wells, as in its own analysis; the validation plates have 300 wells with at least 50 cells
partitions = {
    name: pce.partition_variance(
        data["all_wells"] if name == "Screen" else data["medians"], data["fit"], data["healthy"], n_components, seed
    )
    for name, data in datasets.items()
}
for name, table in partitions.items():
    print(f"\n{name}: share of the variance of the first {n_components} components, by stage")
    print(
        table.groupby("stage")[["variance_pct_of_total", "position_pct_of_total", "status_pct_of_total"]]
        .sum()
        .round(1)
        .to_string()
    )


# In[4]:


fig, axes = plt.subplots(len(datasets), 2, figsize=(17, 8.2), sharey=True)
for i, (name, table) in enumerate(partitions.items()):
    for j, stage in enumerate(("before", "after")):
        ax = axes[i, j]
        subset = table[table["stage"] == stage]
        x = np.arange(len(subset))
        ax.bar(x, subset["variance_pct_of_total"], width=0.8, color=colors["total"], label="Variance of the component")
        ax.bar(x - 0.19, subset["position_pct_of_total"], width=0.36, color=colors["position"], label="Explained by plate position (row + column)")
        ax.bar(x + 0.19, subset["status_pct_of_total"], width=0.36, color=colors["status"], label="Explained by heart status (healthy vs failing)")
        ax.set_xticks(x)
        ax.set_xticklabels(subset["component"], fontsize=9)
        ax.set_title(f"{name}, {'before' if stage == 'before' else 'after'} the correction", fontsize=11)
        if j == 0:
            ax.set_ylabel("% of the total variance")
axes[0, 0].legend(fontsize=8, loc="upper right")
fig.suptitle("Variance of the principal components of the wells, and how much of it position and heart status explain", fontsize=12)
fig.tight_layout()
fig.savefig(figure_dir / "scree_position_and_heart_status.png", dpi=150, bbox_inches="tight")
plt.show()


# ## 2. Survival of the healthy against failing difference
# 
# For every plate, the length of the difference between the mean healthy DMSO profile and the mean failing DMSO profile (over all features, in standard deviations of a single cell), before and after the correction.
# On the screen, the healthy controls are in row B and the failing controls in row E on every platemap, so the correction shifts the two controls differently.
# On the validation plates both cell types sit in every row, at the same positions on every plate, so the tilt shifts the two controls alike and leaves their difference as it is.
# The part of the difference that survives is what position does not explain.

# In[5]:


survival = {
    name: pce.healthy_failing_difference(data["medians"], data["fit"], data["healthy"])
    for name, data in datasets.items()
}
print(pd.DataFrame({name: table[["raw", "corrected", "share_that_survives"]].median() for name, table in survival.items()}).round(2).to_string())


# In[6]:


fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axes[0]
rng = np.random.default_rng(seed)
positions = {"Screen": (0, 1), "Validation": (2.6, 3.6)}
for name, table in survival.items():
    left, right = positions[name]
    for _, plate in table.iterrows():
        ax.plot([left, right], [plate["raw"], plate["corrected"]], color="#9e9e9e", lw=0.6 if name == "Screen" else 1.2, alpha=0.7)
    ax.scatter(np.full(len(table), left) + rng.uniform(-0.1, 0.1, len(table)), table["raw"], s=18, color=colors["position"], zorder=3)
    ax.scatter(np.full(len(table), right) + rng.uniform(-0.1, 0.1, len(table)), table["corrected"], s=18, color=colors["status"], zorder=3)
    ax.hlines([table["raw"].median(), table["corrected"].median()], [left - 0.3, right - 0.3], [left + 0.3, right + 0.3], color="black", lw=2.5, zorder=4)
ax.set_xticks([0, 1, 2.6, 3.6])
ax.set_xticklabels(["Screen\nbefore", "Screen\nafter", "Validation\nbefore", "Validation\nafter"])
ax.set_ylabel("Length of the healthy minus failing DMSO difference")
ax.set_title("Difference between healthy and failing DMSO, per plate (black: median)", fontsize=10)
ax = axes[1]
for k, (name, table) in enumerate(survival.items()):
    ax.scatter(np.full(len(table), k) + rng.uniform(-0.15, 0.15, len(table)), table["share_that_survives"], s=22, color=colors["status"], alpha=0.8)
    median = table["share_that_survives"].median()
    ax.hlines(median, k - 0.3, k + 0.3, color="black", lw=2.5)
    ax.text(k + 0.34, median, f"median {median:.2f}", va="center", fontsize=10)
ax.axhline(1, color="grey", ls="--", lw=1)
ax.set_xlim(-0.6, 1.9)
ax.set_ylim(0, 1.15)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Screen\n(44 plates)", "Validation\n(5 plates)"])
ax.set_ylabel("Share of the difference that survives the correction")
ax.set_title("Corrected over uncorrected length", fontsize=10)
fig.tight_layout()
fig.savefig(figure_dir / "healthy_failing_survival.png", dpi=150, bbox_inches="tight")
plt.show()


# ## 3. Compound placement in the screen
# 
# The tilt of the screen, and its use in the rescue score, assume that compounds were placed on the plates independently of the position.
# Each screen compound sits at one position of its platemap, so the assumption cannot be tested within a compound.
# We test two things: whether the compound classes (pathways) differ by plate row, and whether the compound identifier (a running number) follows the order of the wells on a platemap.

# In[7]:


frames = []
for path in sorted(platemap_dir.glob("*_with_pathways.csv")):
    frame = pd.read_csv(path)
    frame["platemap"] = int(path.name.split("Plate_")[1].split("_")[0])
    frames.append(frame)
platemaps = pd.concat(frames, ignore_index=True)
screen_compounds = platemaps[(platemaps["treatment"] != "DMSO") & platemaps["treatment"].notna()].copy()
screen_compounds["row_index"] = screen_compounds["well_position"].str[0].map({row: i for i, row in enumerate("BCDEFG")})
screen_compounds["column"] = screen_compounds["well_position"].str[1:].astype(int)
screen_compounds["reading_order"] = screen_compounds["row_index"] * 10 + screen_compounds["column"]
screen_compounds["compound_number"] = screen_compounds["treatment"].str.extract(r"(\d+)")[0].astype(float)

pathway_table = pd.crosstab(screen_compounds["Pathway"].fillna("none"), screen_compounds["well_position"].str[0])
chi2, p_pathway, dof, _ = chi2_contingency(pathway_table)
print(f"Pathway against plate row: chi-square p = {p_pathway:.2f} ({len(screen_compounds)} compounds, {pathway_table.shape[0]} pathways)")

placement = []
for platemap, group in screen_compounds.groupby("platemap"):
    result = spearmanr(group["compound_number"], group["reading_order"], nan_policy="omit")
    placement.append({"platemap": platemap, "n_compounds": len(group), "spearman_rho": result.statistic, "p_value": result.pvalue})
placement = pd.DataFrame(placement)
print(placement.round(3).to_string(index=False))


# In[8]:


fig, ax = plt.subplots(figsize=(9, 4.6))
order = placement.sort_values("platemap")
bar_colors = [colors["position"] if p < 0.01 else "#bdbdbd" for p in order["p_value"]]
ax.bar(order["platemap"].astype(str), order["spearman_rho"], color=bar_colors)
ax.axhline(0, color="black", lw=0.8)
ax.set_xlabel("Screen platemap")
ax.set_ylabel("Spearman correlation between compound number\nand position on the plate")
ax.set_title("Compound number follows the position on the plate in platemaps 5, 7, 8, and 10 (red: p < 0.01)", fontsize=10)
fig.tight_layout()
fig.savefig(figure_dir / "screen_compound_placement.png", dpi=150, bbox_inches="tight")
plt.show()

