"""Validate stage-specific source tables and render the complete Day-0 screen."""

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/task1_stage_specific"
F = ROOT / "figures/task1_stage_specific"
F.mkdir(parents=True, exist_ok=True)
region = pd.read_csv(R / "regional_tests.tsv", sep="\t")
peak = pd.read_csv(R / "peak_tests.tsv.gz", sep="\t")
lib = pd.read_csv(R / "peak_by_library.tsv.gz", sep="\t")
rna = pd.read_csv(R / "RNA_supplement.tsv", sep="\t")
members = pd.read_csv(R / "peak_membership.tsv", sep="\t")

# Independently reconcile stage totals with the existing per-library counts.
old = pd.read_csv(
    ROOT / "results/differentiation_fusion25/GSE208248_by_library.tsv.gz", sep="\t"
)
old = old[old.gate.eq("TSS_ge_3")]
c = lib.merge(
    old,
    on=["contrast", "gsm", "peak"],
    suffixes=("_new", "_old"),
    validate="one_to_one",
)
assert len(c) == len(lib)
for k in ["high_open", "low_open", "n_pairs"]:
    assert np.array_equal(c[k + "_new"], c[k + "_old"]), k
for _, d in peak.groupby(["stage", "contrast"]):
    assert len(d) == 565 and d.peak.is_unique
    p = [
        binomtest(int(h), int(h + l)).pvalue if h + l else 1
        for h, l in zip(d.high_only, d.low_only)
    ]
    assert np.allclose(p, d.p, atol=1e-14)
    assert np.allclose(multipletests(p, method="fdr_bh")[1], d.q565)
for _, d in region.groupby(["stage", "contrast", "mode"]):
    assert len(d) == 25 and d.gene.is_unique
    assert np.allclose((d.exceedances + 1) / (d.permutations + 1), d.p)
    assert np.allclose(multipletests(d.p, method="fdr_bh")[1], d.q25)
assert (
    peak[peak.stage.eq("stem") & peak.contrast.eq("2plus_vs_1")].n_pairs.eq(507).all()
)
(R / "validation.json").write_text(
    json.dumps(
        {
            "per_library_counts_match_previous": True,
            "all_peak_p_and_q_reproduced": True,
            "all_regional_p_and_q_reproduced": True,
            "primary_pairs": 507,
            "scope_genes": 25,
            "scope_peaks": 565,
        },
        indent=2,
    )
    + "\n"
)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
primary = region[
    region.stage.eq("stem")
    & region.contrast.eq("2plus_vs_1")
    & region["mode"].eq("opening")
].sort_values("p")
fig, axes = plt.subplots(
    1, 2, figsize=(12, 8), gridspec_kw={"width_ratios": [1, 1.15]}, layout="constrained"
)
ax = axes[0]
y = np.arange(len(primary))
ax.barh(
    y,
    -np.log10(primary.p),
    color=["#ce6a45" if g == "CSRP3" else "#8ca1ad" for g in primary.gene],
    height=0.7,
)
ax.set_yticks(y, primary.gene)
ax.invert_yaxis()
ax.set_xlim(0, 2.5)
for yy, row in zip(y, primary.itertuples()):
    ax.text(2.47, yy, f"{row.q25:.3f}", ha="right", va="center", fontsize=8)
ax.text(2.47, -1.1, "Regional q", ha="right", fontsize=9)
ax.set_xlabel("Regional opening -log10(p)")
ax.set_title("A  All 25 gene neighborhoods", loc="left", pad=22, fontweight="bold")
ax = axes[1]
d = peak[peak.stage.eq("stem") & peak.contrast.eq("2plus_vs_1")].copy()
d = d[np.isfinite(d.FC) & d.FC.gt(0)]
ax.scatter(np.log2(d.FC), -np.log10(d.p), s=19, c="#abb5bc", alpha=0.75, linewidths=0)
marks = {
    "chr22:36352375-36353254": ("MYH9 neighborhood", "#3a73a8", (-1.8, 4.6)),
    "chr11:19201752-19202603": ("CSRP3 proximal", "#ce6a45", (1.15, 2.6)),
    "chr11:19218592-19219518": ("CSRP3 distal", "#ce6a45", (-1.65, 2.6)),
}
for p, (label, col, xy) in marks.items():
    z = d[d.peak.eq(p)].iloc[0]
    x = np.log2(z.FC)
    yy = -np.log10(z.p)
    ax.scatter([x], [yy], s=65, c=col, edgecolors="white", zorder=3)
    ax.annotate(
        label,
        (x, yy),
        xytext=xy,
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": col, "lw": 0.8},
    )
ax.axvline(0, color="#777777", lw=0.7)
ax.set_xlim(-3.5, 3.5)
ax.set_ylim(-0.08, 5.1)
ax.set_xlabel("Peak log2(open-fraction ratio: high / low)")
ax.set_ylabel("Individual-peak -log10(p)")
ax.set_title("B  Candidate peak effects", loc="left", pad=22, fontweight="bold")
fig.suptitle(
    "Undifferentiated cultures | COQ8A >=2 vs 1 UMI | 507 matched pairs", fontsize=13
)
for ext in ["png", "pdf"]:
    fig.savefig(F / f"Figure_1_Day0_complete_screen.{ext}", dpi=180)
plt.close(fig)

fig, axes = plt.subplots(2, 1, figsize=(7.5, 8), layout="constrained")
colors = {
    "GSM6339597": "#2878a3",
    "GSM6339601": "#e08b42",
    "GSM6339599": "#2878a3",
    "GSM6339603": "#e08b42",
}
labels = {
    "GSM6339597": "Source 1, D0",
    "GSM6339601": "Source 2, D0",
    "GSM6339599": "Source 1, D7",
    "GSM6339603": "Source 2, D7",
}
for ax, stage, title in zip(
    axes, ["stem", "differentiated"], ["A  Undifferentiated", "B  Day 7"]
):
    d = lib[
        lib.stage.eq(stage)
        & lib.contrast.eq("2plus_vs_1")
        & lib.peak.eq("chr11:19201752-19202603")
    ]
    for row in d.itertuples():
        ax.plot(
            [0, 1],
            100 * np.array([row.low_open, row.high_open]) / row.n_pairs,
            "o-",
            color=colors[row.gsm],
            label=f"{labels[row.gsm]} (n={row.n_pairs} pairs)",
            lw=1.7,
            ms=7,
        )
    ax.set_xticks([0, 1], ["COQ8A low (1 UMI)", "COQ8A high (>=2 UMI)"])
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0, 10)
    ax.set_ylabel("Accessible nuclei (%)")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
fig.suptitle("CSRP3-proximal peak | chr11:19201752-19202603", fontsize=12)
for ext in ["png", "pdf"]:
    fig.savefig(F / f"Figure_2_CSRP3_stage_separated.{ext}", dpi=180)
plt.close(fig)
print("Validation passed; both figures written.")
