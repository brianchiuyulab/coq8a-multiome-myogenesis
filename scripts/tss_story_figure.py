"""Figure 1: fixed gene set, TSS-centred fragments, and transparent locus nomination.

The TSS display is descriptive. The exploratory ranking uses the separately
specified +/-100 kb candidate peaks and same-data RNA links. In particular,
the TSS heatmap does not select MYOD1 or constitute evidence for the six-peak
fold ratio.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


INK = "#202935"
MUTED = "#687582"
BLUE = "#2166ac"
RED = "#b2182b"
LIGHT = "#e4e9ee"


def _load_profile(tables: Path):
    data = pd.read_csv(tables / "tss_fragment_profile_221.tsv.gz", sep="\t")
    if len(data) != 22100 or data.gene.nunique() != 221:
        raise ValueError("Expected 221 genes x 100 TSS bins")
    if data.groupby("gene").size().ne(100).any():
        raise ValueError("Each gene must have exactly 100 bins")
    if not data.n_high.eq(958).all() or not data.n_low.eq(958).all():
        raise ValueError("TSS profile must use the 958 primary matched pairs")
    data = data.sort_values(["gene", "bin_start_bp"])
    ordered = (
        data.groupby("gene")[["high_cuts", "low_cuts"]]
        .sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index.tolist()
    )
    # Sorting uses only the pooled intensity, never high-minus-low direction.
    high = data.pivot(index="gene", columns="bin_start_bp", values="high_cuts").loc[
        ordered
    ]
    low = data.pivot(index="gene", columns="bin_start_bp", values="low_cuts").loc[
        ordered
    ]
    if high.columns.tolist() != low.columns.tolist():
        raise ValueError("High/low TSS bin positions differ")
    return ordered, high.to_numpy(dtype=float) / 958, low.to_numpy(dtype=float) / 958


def _panel(ax, letter, title):
    ax.set_title(
        f"{letter}  {title}", loc="left", fontsize=11, fontweight="bold", pad=10
    )


def make_story_figure(tables: Path, out: Path) -> None:
    order, high, low = _load_profile(tables)
    high = gaussian_filter1d(high, sigma=1, axis=1, mode="nearest")
    low = gaussian_filter1d(low, sigma=1, axis=1, mode="nearest")
    combined = np.concatenate([high.ravel(), low.ravel()])
    upper = float(np.quantile(combined, 0.995))
    if upper <= 0:
        raise ValueError("TSS profile has no signal")
    ranking = pd.read_csv(tables / "candidate_gene_ranking.tsv", sep="\t")
    if len(ranking) != 221 or ranking.iloc[0].gene != "MYOD1":
        raise ValueError("Unexpected exploratory ranking")
    sensitivity = pd.read_csv(tables / "myod1_locus_summary.tsv", sep="\t")
    extreme = sensitivity[
        (sensitivity.gate == "TSS_ge_3")
        & (sensitivity.contrast == "3plus_vs_1")
        & (sensitivity.region_set == "positive_link_nonMRF_tss2_6")
    ]
    if len(extreme) != 1 or not np.isclose(extreme.iloc[0].fold_open, 1.28, atol=0.001):
        raise ValueError("Expected the six-peak 1.280 sensitivity row")
    result = extreme.iloc[0]

    fig = plt.figure(figsize=(15.2, 10.4), facecolor="white")
    grid = GridSpec(
        3,
        3,
        figure=fig,
        height_ratios=[0.95, 1.1, 5.1],
        width_ratios=[1, 1, 1.15],
        hspace=0.34,
        wspace=0.34,
        left=0.075,
        right=0.96,
        top=0.89,
        bottom=0.105,
    )
    fig.suptitle(
        "From a fixed myogenesis gene set to an exploratory MYOD1 locus",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.075,
        0.935,
        "Same-nucleus RNA + ATAC  |  COQ8A >=2 versus 1 UMI  |  TSS QC >=3  |  958 depth-matched pairs",
        color=MUTED,
        fontsize=10,
    )

    ax = fig.add_subplot(grid[0, :])
    ax.axis("off")
    _panel(ax, "A", "Define the search before looking at COQ8A effects")
    steps = [
        ("221 genes", "Hallmark + Reactome"),
        ("7,132 peaks", "within 100 kb of a TSS"),
        ("6,871 tested", "prevalence-qualified"),
        ("0 FDR hits", "primary ATAC screen"),
    ]
    for i, (head, detail) in enumerate(steps):
        x = 0.025 + 0.25 * i
        ax.text(
            x,
            0.44,
            head,
            transform=ax.transAxes,
            fontsize=12,
            color=BLUE if i < 3 else INK,
            fontweight="bold",
        )
        ax.text(x, 0.13, detail, transform=ax.transAxes, fontsize=8.5, color=MUTED)
        if i < 3:
            ax.annotate(
                "",
                xy=(x + 0.225, 0.48),
                xytext=(x + 0.19, 0.48),
                xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "color": MUTED},
            )

    ax_meta = fig.add_subplot(grid[1, :2])
    x = np.arange(100) / 10 - 4.95
    ax_meta.plot(x, low.mean(axis=0), color=BLUE, lw=1.7, label="COQ8A low")
    ax_meta.plot(x, high.mean(axis=0), color=RED, lw=1.7, label="COQ8A high")
    ax_meta.axvline(0, color=MUTED, ls="--", lw=0.8)
    ax_meta.set_xlim(-5, 5)
    ax_meta.set_ylabel("Mean Tn5 cuts / nucleus", fontsize=8.5)
    ax_meta.legend(frameon=False, loc="upper right", fontsize=8, ncol=2)
    _panel(ax_meta, "B", "ATAC insertion profile around the 221 gene TSSs")
    ax_meta.text(
        0.01,
        0.91,
        f"High/low cuts across these windows: {high.sum() / low.sum():.3f}x",
        transform=ax_meta.transAxes,
        fontsize=8,
        color=MUTED,
        va="top",
    )

    ax_note = fig.add_subplot(grid[1, 2])
    ax_note.axis("off")
    ax_note.text(
        0,
        0.96,
        "TSS heatmap: descriptive view\n"
        "Locus ranking: all candidate peaks within 100 kb\n"
        "MYOD1 six-peak ratio: a separate, exploratory subset",
        transform=ax_note.transAxes,
        va="top",
        fontsize=9.3,
        color=INK,
        linespacing=1.65,
    )

    matrix_grid = grid[2, :2].subgridspec(1, 2, wspace=0.08)
    norm = Normalize(vmin=0, vmax=upper)
    axes = [fig.add_subplot(matrix_grid[0, i]) for i in range(2)]
    for ax_h, values, name in zip(axes, [low, high], ["COQ8A low", "COQ8A high"]):
        image = ax_h.imshow(
            values, aspect="auto", interpolation="nearest", cmap="Blues", norm=norm
        )
        ax_h.set_title(name, fontsize=10, pad=7)
        ax_h.set_xticks([0, 49.5, 99], ["-5", "TSS", "+5"])
        ax_h.set_xlabel("Distance from gene TSS (kb)", fontsize=8.5)
        ax_h.set_yticks([0, 54, 109, 164, 220], ["1", "55", "110", "165", "221"])
        ax_h.tick_params(axis="both", labelsize=8, length=0)
        ax_h.axvline(49.5, color="white", lw=0.6, alpha=0.8)
        ax_h.axhline(order.index("MYOD1"), color=RED, lw=0.65, alpha=0.9)
        for spine in ax_h.spines.values():
            spine.set_visible(False)
    axes[0].set_ylabel("221 genes, ordered by pooled TSS signal", fontsize=8.5)
    axes[1].set_yticklabels([])
    fig.colorbar(
        image,
        ax=axes,
        location="bottom",
        fraction=0.025,
        pad=0.10,
        label="Mean Tn5 cuts / nucleus / 100 bp (shared colour scale)",
    )
    axes[1].text(
        1.01,
        order.index("MYOD1") / 220,
        "MYOD1",
        transform=axes[1].transAxes,
        color=RED,
        fontsize=8,
        va="center",
    )

    right = grid[2, 2].subgridspec(2, 1, height_ratios=[1.65, 1], hspace=0.47)
    ax_rank = fig.add_subplot(right[0])
    top = ranking.head(10).iloc[::-1]
    colors = [RED if gene == "MYOD1" else "#aab8c7" for gene in top.gene]
    ax_rank.barh(top.gene, top.n_linked_coq_3of4, color=colors, height=0.65)
    ax_rank.set_xlim(0, 5.9)
    ax_rank.set_xlabel("Positive RNA link + ATAC rise in >=3/4 libraries", fontsize=8.5)
    ax_rank.tick_params(axis="both", labelsize=8)
    ax_rank.grid(axis="x", color=LIGHT, lw=0.7)
    ax_rank.set_axisbelow(True)
    for i, val in enumerate(top.n_linked_coq_3of4):
        ax_rank.text(val + 0.1, i, str(val), va="center", fontsize=8, color=INK)
    _panel(ax_rank, "C", "Exploratory ranking of all 221 genes")

    ax_zoom = fig.add_subplot(right[1])
    _panel(ax_zoom, "D", "MYOD1: zoom from 32 regions to 6")
    ax_zoom.text(
        0,
        0.94,
        "32 nearby -> 19 shared -> 11 linked (TSS QC >=2) -> 6 motif-low",
        transform=ax_zoom.transAxes,
        fontsize=8,
        color=INK,
    )
    ax_zoom.barh(
        [1, 0],
        [result.high_open_pct, result.low_open_pct],
        color=[RED, BLUE],
        height=0.5,
    )
    ax_zoom.set_yticks([1, 0], ["COQ8A >=3", "COQ8A =1"], fontsize=8)
    ax_zoom.set_xlim(0, 16.5)
    ax_zoom.set_ylim(-1.0, 1.75)
    ax_zoom.set_xticks([0, 5, 10, 15])
    ax_zoom.set_xlabel("Open nuclei at six peaks (%)", fontsize=8)
    ax_zoom.tick_params(axis="x", labelsize=8)
    ax_zoom.grid(axis="x", color=LIGHT, lw=0.7)
    ax_zoom.set_axisbelow(True)
    for y, value in [(1, result.high_open_pct), (0, result.low_open_pct)]:
        ax_zoom.text(value + 0.2, y, f"{value:.2f}%", va="center", fontsize=8)
    ax_zoom.text(
        0.98,
        0.10,
        f"Ratio {result.fold_open:.2f}x  |  {int(result.n_pairs)} pairs",
        transform=ax_zoom.transAxes,
        ha="right",
        fontsize=9,
        fontweight="bold",
        color=RED,
    )
    for side in ("top", "right", "left"):
        ax_zoom.spines[side].set_visible(False)

    fig.text(
        0.075,
        0.022,
        "TSS signal comes from indexed deduplicated ATAC fragments (GENCODE gene TSS, +/-5 kb, 100-bp bins; "
        "100-bp smoothing). Genes sorted by pooled signal. TSS display does not select MYOD1 or test the six peaks.",
        fontsize=8,
        color=MUTED,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(out.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
