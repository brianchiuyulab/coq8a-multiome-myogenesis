"""Display the complete temporal peak screen and the eight lowest raw p values."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


COLORS = {"Early": "#177D91", "Middle": "#D18B20", "Late": "#975695"}
SETTINGS = [
    ("TSS_ge_2", "2plus_vs_1"),
    ("TSS_ge_3", "2plus_vs_1"),
    ("TSS_ge_2", "3plus_vs_1"),
    ("TSS_ge_3", "3plus_vs_1"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temporal", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.temporal / "single_peak_effects.tsv", sep="\t")
    library = pd.read_csv(args.temporal / "single_peak_by_library.tsv.gz", sep="\t")
    primary = data[
        (data.gate == "TSS_ge_3") & (data.contrast == "3plus_vs_1")
    ].sort_values(["p_exact", "peak"])
    lead = primary.head(8)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(16, 10.5), layout="constrained")
    fig.suptitle(
        "Individual differentiation-opening peaks and COQ8A expression\n"
        "Primary display: COQ8A ≥3 vs 1 UMI · TSS ≥3 · 201 matched pairs · 4 libraries / 2 source lines",
        fontsize=16,
    )
    ax = axes[0, 0]
    for phase, color in COLORS.items():
        d = primary[primary.phase.eq(phase)]
        ax.scatter(
            d.delta_pp,
            -np.log10(d.p_exact),
            c=color,
            s=25,
            alpha=0.8,
            label=f"{phase} (n={len(d)})",
            linewidths=0,
        )
    ax.axvline(0, color=".6", lw=0.8)
    ax.axhline(-np.log10(0.05), color=".65", ls="--", lw=0.8)
    for row in primary.head(3).itertuples():
        ax.annotate(
            row.nearby_candidate_genes,
            (row.delta_pp, -np.log10(row.p_exact)),
            xytext=(-5, 9),
            textcoords="offset points",
            ha="right",
            fontsize=10,
        )
    ax.set(
        xlabel="Open-nucleus fraction: high − low (percentage points)",
        ylabel="−log10(raw paired exact p)",
        title="A   All 229 peaks, including 191 early peaks",
    )
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.text(
        0.02,
        0.68,
        "BH q < 0.05: 0 / 229\nBH q < 0.10: 3 / 229",
        transform=ax.transAxes,
        fontsize=10,
    )
    ax.margins(y=0.18)
    labels = [
        f"{r.nearby_candidate_genes} [{r.phase}]\n{r.peak}" for r in lead.itertuples()
    ]
    ax = axes[0, 1]
    for i, row in enumerate(lead.itertuples()):
        ax.plot([row.low_open_pct, row.high_open_pct], [i, i], c=".65", lw=1.5)
        ax.scatter(
            row.low_open_pct,
            i,
            c="#4B5D6B",
            s=38,
            label="COQ8A low" if i == 0 else None,
        )
        ax.scatter(
            row.high_open_pct,
            i,
            c="#C64B3C",
            s=38,
            label="COQ8A high" if i == 0 else None,
        )
        ax.text(
            26,
            i,
            f"{row.fold_open:.2f}×   p={row.p_exact:.2g}   q={row.q_229_peaks:.3f}",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(range(len(lead)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set(
        xlim=(0, 51),
        xlabel="Nuclei with an accessible peak (%)",
        title="B   Eight lowest p values, both directions retained",
    )
    ax.set_xticks([0, 5, 10, 15, 20, 25])
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax = axes[1, 0]
    matrix = np.array(
        [
            [
                data[
                    (data.gate == gate)
                    & (data.contrast == contrast)
                    & (data.peak == peak)
                ]
                .iloc[0]
                .fold_open
                for gate, contrast in SETTINGS
            ]
            for peak in lead.peak
        ]
    )
    ax.imshow(
        np.log2(matrix),
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5),
        aspect="auto",
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j,
                i,
                f"{matrix[i,j]:.2f}×",
                ha="center",
                va="center",
                color="white" if abs(np.log2(matrix[i, j])) > 1.6 else "black",
            )
    ax.set_yticks(range(len(lead)), lead.nearby_candidate_genes)
    ax.set_xticks(
        range(4),
        ["≥2 vs 1\nTSS ≥2", "≥2 vs 1\nTSS ≥3", "≥3 vs 1\nTSS ≥2", "≥3 vs 1\nTSS ≥3"],
    )
    ax.set_title("C   Same peaks across all four settings (high/low ratio)")
    ax = axes[1, 1]
    lib = library[(library.gate == "TSS_ge_3") & (library.contrast == "3plus_vs_1")]
    gsms = sorted(lib.gsm.unique())
    matrix = (
        lib.pivot(index="peak", columns="gsm", values="delta_pp")
        .loc[lead.peak, gsms]
        .to_numpy()
    )
    limit = max(10, np.max(np.abs(matrix)))
    im = ax.imshow(
        matrix,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        aspect="auto",
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j,
                i,
                f"{matrix[i,j]:+.1f}",
                ha="center",
                va="center",
                color="white" if abs(matrix[i, j]) > limit * 0.6 else "black",
            )
    ax.set_yticks(range(len(lead)), lead.nearby_candidate_genes)
    ax.set_xticks(
        range(4),
        [
            "Line 1\nstem",
            "Line 1\ndifferentiated",
            "Line 2\nstem",
            "Line 2\ndifferentiated",
        ],
    )
    ax.set_title("D   Library effects (percentage-point difference)")
    fig.colorbar(im, ax=ax, shrink=0.75, label="High − low (pp)")
    for suffix in ["png", "pdf"]:
        fig.savefig(
            args.out / ("Figure_temporal_single_peaks." + suffix),
            dpi=200,
            facecolor="white",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
