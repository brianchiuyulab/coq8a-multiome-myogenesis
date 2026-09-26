"""Plot external timing and the complete four-setting COQ8A module comparison."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

MODULES = [
    ("All_external_mapped", "All externally mapped"),
    ("All_opening", "All opening"),
    ("Early_by24h", "Early: half-rise by 24 h"),
    ("Middle_24to48h", "Middle: 24-48 h"),
    ("Late_after48h", "Late: after 48 h"),
    ("Low_change_profile", "Low-change profile"),
    ("Closing", "Closing"),
    ("Earliest_20pct", "Earliest 20%"),
    ("Earliest_25pct", "Earliest 25%"),
    ("Earliest_33pct", "Earliest 33%"),
]
COLORS = {
    "Early_by24h": "#137C8B",
    "Middle_24to48h": "#DC9B32",
    "Late_after48h": "#9E5A91",
}


def save(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=200, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--temporal", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(a.temporal / "temporal_effect_grid.tsv", sep="\t")
    m = d[
        (d.scope == "programme") & (d.promoter_kb == 2) & (d.reciprocal_overlap == 0.5)
    ]
    definitions = pd.read_csv(
        a.temporal / "temporal_set_definitions.tsv", sep="\t"
    ).set_index("set_id")
    mapped = pd.read_csv(a.temporal / "target_temporal_annotations.tsv.gz", sep="\t")
    opening = mapped[
        (mapped.promoter_kb == 2) & (mapped.reciprocal_overlap >= 0.5) & mapped.opening
    ].sort_values(["half_rise_hour", "peak"])
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(11, 6.5))
    gs = fig.add_gridspec(
        1,
        2,
        left=0.10,
        right=0.95,
        bottom=0.19,
        top=0.80,
        wspace=0.48,
        width_ratios=[1, 1.1],
    )
    fig.suptitle(
        "External ATAC timing defines the regions to test",
        x=0.06,
        y=0.96,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    count = lambda name: int(definitions.loc["P2_O50_" + name, "n_peaks"])
    fig.text(
        0.06,
        0.89,
        f'221 genes -> 5,097 candidate peaks -> {count("All_external_mapped"):,} external matches -> {count("All_opening")} opening peaks',
        fontsize=11,
    )
    ax = fig.add_subplot(gs[0, 0])
    v = opening[["p0", "p24", "p48", "p72"]].to_numpy()
    normalized = (v - v.min(axis=1, keepdims=True)) / np.maximum(
        np.ptp(v, axis=1)[:, None], 1e-12
    )
    im = ax.imshow(
        normalized,
        aspect="auto",
        cmap="viridis",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    ax.set_xticks(range(4), ["0", "24", "48", "72"])
    ax.set_xlabel("Hours after differentiation induction")
    ax.set_ylabel("Opening peaks, sorted by external half-rise time")
    ax.set_title("A  External time-course ATAC", loc="left", fontweight="bold", pad=13)
    n1 = count("Early_by24h")
    n2 = count("Middle_24to48h")
    n3 = count("Late_after48h")
    for y in [n1 - 0.5, n1 + n2 - 0.5]:
        ax.axhline(y, color="white", lw=1)
    ax.set_yticks(
        [max(0, n1 / 2), n1 + n2 / 2, n1 + n2 + n3 / 2],
        [f"Early\n{n1}", f"Middle\n{n2}", f"Late {n3}"],
        fontsize=9,
    )
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.18)
    cb.set_label("Within-peak scaled accessibility (0-1)", fontsize=9)
    ax = fig.add_subplot(gs[0, 1])
    for module, mask in [
        ("Early_by24h", opening.half_rise_hour <= 24),
        (
            "Middle_24to48h",
            (opening.half_rise_hour > 24) & (opening.half_rise_hour <= 48),
        ),
        ("Late_after48h", opening.half_rise_hour > 48),
    ]:
        subset = normalized[mask.to_numpy()]
        if len(subset):
            ax.plot(
                [0, 24, 48, 72],
                np.median(subset, axis=0),
                "o-",
                color=COLORS[module],
                lw=2,
                label=f'{module.split("_")[0]} ({len(subset)} peaks)',
            )
    ax.set_xticks([0, 24, 48, 72])
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("Hours after differentiation induction")
    ax.set_ylabel("Median within-peak scaled accessibility")
    ax.set_title(
        "B  Early, middle and late profiles", loc="left", fontweight="bold", pad=13
    )
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.text(
        0.06,
        0.05,
        "GSE109828: temporal definition; positive endpoint change in both experiments.\nCOQ8A expression and target-dataset effects were not used to select these regions.",
        fontsize=10,
    )
    save(fig, a.out / "Figure_external_timing")

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_axes([0.31, 0.17, 0.58, 0.66])
    fig.suptitle(
        "Does higher COQ8A accompany early-region accessibility?",
        x=0.05,
        y=0.97,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.05,
        0.91,
        "Same externally defined regions under all four COQ8A / TSS settings",
        fontsize=11,
    )
    settings = [
        ("TSS_ge_2", "2plus_vs_1"),
        ("TSS_ge_3", "2plus_vs_1"),
        ("TSS_ge_2", "3plus_vs_1"),
        ("TSS_ge_3", "3plus_vs_1"),
    ]
    matrix = np.full((len(MODULES), 4), np.nan)
    labels = {}
    for i, (module, _) in enumerate(MODULES):
        for j, (gate, contrast) in enumerate(settings):
            sub = m[
                (m.module == module)
                & (m.gate == gate)
                & (m.contrast == contrast)
                & (m.stage_analysis == "pooled")
            ]
            if sub.empty:
                continue
            row = sub.iloc[0]
            matrix[i, j] = row.fold_open
            labels[i, j] = f"{row.fold_open:.3f}x\np={row.p_pair:.3g}"
    extent = max(0.08, min(0.35, np.nanmax(np.abs(matrix - 1))))
    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=1 - extent, vcenter=1, vmax=1 + extent),
    )
    for (i, j), label in labels.items():
        ax.text(
            j,
            i,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color="white" if abs(matrix[i, j] - 1) > extent * 0.7 else "#202020",
        )
    ax.set_yticks(
        range(len(MODULES)),
        [f"{label}\n{count(module)} peaks" for module, label in MODULES],
        fontsize=10,
    )
    ax.set_xticks(
        range(4),
        [
            "COQ8A >=2 vs 1\nTSS >=2\n1,020 pairs",
            "COQ8A >=2 vs 1\nTSS >=3\n958 pairs",
            "COQ8A >=3 vs 1\nTSS >=2\n212 pairs",
            "COQ8A >=3 vs 1\nTSS >=3\n201 pairs",
        ],
        fontsize=9,
    )
    cb = fig.colorbar(image, ax=ax, orientation="vertical", fraction=0.05, pad=0.04)
    cb.set_label("High / low open-fraction ratio")
    fig.text(
        0.05,
        0.05,
        "Four libraries from two source lines; nuclei matched within library on RNA and ATAC depth.\nP: paired-nucleus tests. Complete within-setting q values and per-library effects are in source tables.",
        fontsize=10,
    )
    save(fig, a.out / "Figure_COQ8A_temporal_sensitivity")

    primary = m[(m.gate == "TSS_ge_3") & (m.contrast == "3plus_vs_1")]
    fig, axes = plt.subplots(
        1, 2, figsize=(12, 6.5), gridspec_kw={"width_ratios": [1, 1]}
    )
    fig.subplots_adjust(left=0.22, right=0.93, top=0.78, bottom=0.28, wspace=0.80)
    fig.suptitle(
        "Cell state and genes within the temporal accessibility signal",
        x=0.05,
        y=0.97,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.05,
        0.90,
        "COQ8A >=3 vs 1 UMI | TSS >=3 | fixed external timing and region mapping",
        fontsize=11,
    )
    names = ["Early_by24h", "Middle_24to48h", "All_opening", "Low_change_profile"]
    ax = axes[0]
    for j, (stage, color, label) in enumerate(
        [
            ("stem", "#137C8B", "Undifferentiated"),
            ("differentiated", "#CE7751", "Differentiated"),
            ("pooled", "#303642", "Pooled"),
        ]
    ):
        sub = (
            primary[(primary.stage_analysis == stage) & primary.module.isin(names)]
            .set_index("module")
            .reindex(names)
        )
        y = np.arange(len(names)) + (j - 1) * 0.2
        ax.errorbar(
            sub.delta_pp,
            y,
            xerr=[
                sub.delta_pp - sub.delta_ci_low_pp,
                sub.delta_ci_high_pp - sub.delta_pp,
            ],
            fmt="o",
            ms=5,
            color=color,
            label=label,
            capsize=2,
        )
    ax.axvline(0, color="#AAAAAA", lw=1)
    ax.invert_yaxis()
    ax.set_yticks(
        range(len(names)),
        ["Early by 24 h", "Middle: 24-48 h", "All opening", "Low-change"],
    )
    ax.set_xlabel("High minus low accessibility\n(percentage points; 95% CI)")
    ax.set_title("A  State-resolved effects", loc="left", fontweight="bold", pad=18)
    ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.45, -0.27),
        fontsize=9,
        ncol=3,
    )
    genes = d[
        (d.scope == "gene")
        & (d.module == "Middle_24to48h")
        & (d.gate == "TSS_ge_3")
        & (d.contrast == "3plus_vs_1")
        & (d.stage_analysis == "pooled")
    ].copy()
    ax = axes[1]
    genes["x"] = np.log2(genes.fold_open.replace(0, np.nan))
    genes["y"] = -np.log10(genes.p_pair.clip(lower=1e-300))
    ax.scatter(
        genes.x,
        genes.y,
        s=30 + 12 * genes.n_peaks,
        c=np.where(genes.q_within_family < 0.05, "#137C8B", "#ADB8C3"),
        edgecolor="white",
        lw=0.5,
    )
    ax.axvline(0, color="#AAAAAA", lw=1)
    ax.axhline(-np.log10(0.05), color="#AAAAAA", lw=1, ls="--")
    for j, r in enumerate(genes.sort_values("p_pair").head(5).itertuples()):
        ax.annotate(
            r.genes,
            (r.x, r.y),
            xytext=(4, 6 + 8 * (j % 2)),
            textcoords="offset points",
            fontsize=9,
        )
    ax.set_xlabel("log2(high / low open fraction)")
    ax.set_ylabel("-log10(paired-nucleus p)")
    ax.set_title(
        f"B  All {len(genes)} genes in the middle module",
        loc="left",
        fontweight="bold",
        pad=18,
    )
    fig.text(
        0.05,
        0.04,
        "Gene assignments use the original 221-gene candidate region map. Dot size: number of peaks.\nTeal gene dots: BH q<0.05 within this module; labels: five smallest p values. All genes remain in source tables.",
        fontsize=10,
    )
    save(fig, a.out / "Figure_temporal_states_and_genes")


if __name__ == "__main__":
    main()
