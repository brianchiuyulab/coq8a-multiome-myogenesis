"""Show phase-level effects and every middle-opening peak without target filtering."""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def save(fig, path):
    for ext in ["png", "pdf"]:
        fig.savefig(
            path.with_suffix("." + ext), dpi=180, facecolor="white", bbox_inches="tight"
        )
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--temporal", type=Path, required=True)
    p.add_argument("--figures", type=Path, required=True)
    a = p.parse_args()
    a.figures.mkdir(parents=True, exist_ok=True)
    allpeaks = pd.read_csv(a.temporal / "single_peak_effects.tsv", sep="\t")
    middle = allpeaks[allpeaks.phase.eq("Middle")].copy()
    for _, s in middle.groupby(["gate", "contrast"]):
        assert len(s) == 34
        middle.loc[s.index, "q_middle_34"] = multipletests(s.p_exact, method="fdr_bh")[
            1
        ]
    middle.to_csv(
        a.temporal / "middle_peak_followup.tsv", sep="\t", index=False, na_rep="NA"
    )
    effects = pd.read_csv(a.temporal / "temporal_effect_grid.tsv", sep="\t")
    libraries = pd.read_csv(a.temporal / "temporal_effects_by_library.tsv.gz", sep="\t")

    def primary(x):
        return x[(x.gate == "TSS_ge_3") & (x.contrast == "3plus_vs_1")]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.8), sharey=True, layout="constrained")
    phase_ids = ["P2_O50_Early_by24h", "P2_O50_Middle_24to48h", "P2_O50_Late_after48h"]
    phase_lib = primary(libraries).query("set_id in @phase_ids")
    ymax = phase_lib[["low_open_pct", "high_open_pct"]].to_numpy().max() * 1.15
    for ax, (label, suffix) in zip(
        axes,
        [
            ("Early", "Early_by24h"),
            ("Middle", "Middle_24to48h"),
            ("Late", "Late_after48h"),
        ],
    ):
        sid = "P2_O50_" + suffix
        row = (
            primary(effects).query('set_id==@sid and stage_analysis=="pooled"').iloc[0]
        )
        sub = primary(libraries).query("set_id==@sid").sort_values("gsm")
        for r, c in zip(sub.itertuples(), ["#176D9C", "#91BDD1", "#D87728", "#E6B987"]):
            ax.plot(
                [0, 1],
                [r.low_open_pct, r.high_open_pct],
                "o-",
                c=c,
                alpha=0.85,
                ms=5,
                lw=1,
            )
        ax.plot(
            [0, 1],
            [row.low_open_pct, row.high_open_pct],
            "D-",
            c="black",
            lw=2,
            ms=6,
            label="Pooled pairs",
        )
        ax.set_xticks([0, 1], ["COQ8A low\n1 UMI", "COQ8A high\n≥3 UMI"])
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(0, ymax)
        ax.set_title(
            f"{label}: {int(row.n_peaks)} peaks\nFC {row.fold_open:.3f} · p {row.p_pair:.4f}\n10-programme q {row.q_within_family:.4f}",
            fontsize=11,
        )
        ax.grid(axis="y", alpha=0.15)
    axes[0].set_ylabel("Mean fraction of accessible peaks per nucleus (%)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(
        "Externally timed regions: identical COQ8A groups and QC\nTSS ≥3; 201 matched pairs; colored lines = four libraries from two source lines",
        fontsize=12,
    )
    save(fig, a.figures / "Figure_three_temporal_phases")
    mid = primary(middle).sort_values(["p_exact", "peak"]).reset_index(drop=True)
    lib = primary(pd.read_csv(a.temporal / "single_peak_by_library.tsv.gz", sep="\t"))
    delta = (
        lib.pivot(index="peak", columns="gsm", values="delta_pp")
        .loc[mid.peak, sorted(lib.gsm.unique())]
        .to_numpy()
    )
    fig = plt.figure(figsize=(14, 13))
    grid = fig.add_gridspec(
        1,
        3,
        width_ratios=[2, 3, 3],
        left=0.31,
        right=0.98,
        top=0.89,
        bottom=0.1,
        wspace=0.12,
    )
    ax = fig.add_subplot(grid[0, 0])
    v = mid[["low_open_pct", "high_open_pct"]].to_numpy()
    im = ax.imshow(v, aspect="auto", cmap="YlOrRd", vmin=0, vmax=25)
    labels = [f"{r.nearby_candidate_genes}   {r.peak}" for r in mid.itertuples()]
    ax.set_yticks(range(34), labels, fontsize=8)
    ax.set_xticks([0, 1], ["Low", "High"])
    ax.set_title("Open nuclei (%)", fontsize=11)
    for i in range(34):
        for j in range(2):
            ax.text(
                j,
                i,
                f"{v[i,j]:.1f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if v[i, j] > 16 else "black",
            )
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, 0.045, pos.width, 0.012])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=8)
    ax = fig.add_subplot(grid[0, 1])
    lim = max(abs(delta.min()), abs(delta.max()))
    im = ax.imshow(
        delta,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
    )
    ax.set_yticks([])
    ax.set_xticks(
        range(4),
        ["Line1\nstem", "Line1\ndiff.", "Line2\nstem", "Line2\ndiff."],
        fontsize=8,
    )
    ax.set_title("High − low (percentage points)", fontsize=11)
    for i in range(34):
        for j in range(4):
            ax.text(
                j,
                i,
                f"{delta[i,j]:+.1f}",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if abs(delta[i, j]) > lim * 0.65 else "black",
            )
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, 0.045, pos.width, 0.012])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=8)
    ax = fig.add_subplot(grid[0, 2])
    ax.set_xlim(0, 3)
    ax.set_ylim(33.5, -0.5)
    ax.axis("off")
    for j, h in enumerate(["FC", "p", "q"]):
        ax.text(j + 0.48, -1, h, ha="center", fontsize=9, weight="bold")
    for i, r in enumerate(mid.itertuples()):
        if r.p_exact < 0.05:
            ax.axhspan(i - 0.5, i + 0.5, color="#F5E8C6", zorder=0)
        for j, t in enumerate(
            [f"{r.fold_open:.3f}", f"{r.p_exact:.3g}", f"{r.q_middle_34:.4f}"]
        ):
            ax.text(j + 0.48, i, t, ha="center", va="center", fontsize=9)
    fig.suptitle(
        "All 34 middle-opening peaks: COQ8A-associated accessibility",
        fontsize=16,
        y=0.97,
    )
    fig.text(0.31, 0.925, "COQ8A ≥3 vs 1 UMI; TSS ≥3; 201 matched pairs", fontsize=11)
    save(fig, a.figures / "Figure_middle_34_peak_heatmap")
    print(
        mid[
            [
                "nearby_candidate_genes",
                "fold_open",
                "p_exact",
                "q_middle_34",
                "q_229_peaks",
            ]
        ]
        .head(3)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
