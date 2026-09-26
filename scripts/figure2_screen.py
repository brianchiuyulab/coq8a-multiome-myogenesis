"""Figure 2: parallel complete screens for both COQ8A count thresholds.

The same fixed 221-gene search and QC rules are applied to both contrasts.
Locus analyses follow these full screens; the stronger-count subset is not
presented as independent replication of the broader contrast.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INK = "#202935"
MUTED = "#687582"
PALE = "#adb9c6"
BLUE = "#2166ac"
RED = "#b2182b"
GRID = "#e3e8ed"
GATE = "TSS_ge_3"
CONTRASTS = ("2plus_vs_1", "3plus_vs_1")


def read_screen(
    tables: Path, contrast: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    peaks = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peaks = peaks[
        (peaks.gate == GATE)
        & (peaks.contrast == contrast)
        & peaks.q_candidate_peaks.notna()
    ].copy()
    genes = pd.read_csv(tables / "gene_effects_pooled.tsv", sep="\t")
    genes = genes[
        (genes.gate == GATE) & (genes.contrast == contrast) & (genes.modality == "ATAC")
    ].copy()
    ranking = pd.read_csv(tables / f"candidate_gene_ranking_{contrast}.tsv", sep="\t")
    expected_peaks = 6_871 if contrast == "2plus_vs_1" else 3_731
    if len(peaks) != expected_peaks or len(genes) != 221 or len(ranking) != 221:
        raise ValueError(f"Incomplete fixed-gene screen for {contrast}")
    if set(genes.gene) != set(ranking.gene) or genes.q_221.lt(0.05).any():
        raise ValueError(f"Unexpected gene family or FDR result for {contrast}")
    return peaks, genes, ranking


def plot_peak_family(ax, peaks: pd.DataFrame, label: str, n_pairs: int) -> None:
    ax.scatter(
        peaks.delta_pp,
        -np.log10(peaks.p_pair_binomial.clip(lower=1e-300)),
        s=5,
        alpha=0.40,
        color=PALE,
        rasterized=True,
    )
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_title(
        f"{label}  {len(peaks):,} eligible peaks; {n_pairs} matched pairs",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.set_xlabel("High - low opening (percentage points)")
    ax.set_ylabel("-log10 nominal paired p")
    ax.text(
        0.02,
        0.97,
        f"Whole-peak-family minimum BH q = {peaks.q_candidate_peaks.min():.3f}",
        transform=ax.transAxes,
        va="top",
        color=MUTED,
        fontsize=8.1,
    )


def make_figure(tables: Path, out: Path) -> None:
    low_peaks, low_genes, low_rank = read_screen(tables, CONTRASTS[0])
    high_peaks, high_genes, high_rank = read_screen(tables, CONTRASTS[1])
    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}
    )
    fig = plt.figure(figsize=(13.6, 8.5), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1, 1],
        hspace=0.55,
        wspace=0.34,
        left=0.09,
        right=0.97,
        top=0.85,
        bottom=0.13,
    )
    fig.suptitle(
        "Complete 221-gene screens at both COQ8A thresholds",
        x=0.09,
        y=0.965,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.09,
        0.916,
        "Same four libraries, same TSS >=3 QC and search space | low group = 1 RNA UMI",
        fontsize=9.5,
        color=MUTED,
    )

    plot_peak_family(
        fig.add_subplot(grid[0, 0]), low_peaks, "A  COQ8A >=2 vs 1 UMI", 958
    )
    plot_peak_family(
        fig.add_subplot(grid[0, 1]), high_peaks, "B  COQ8A >=3 vs 1 UMI", 201
    )

    ax = fig.add_subplot(grid[1, 0])
    effects = (
        low_genes[["gene", "difference", "q_221"]]
        .rename(columns={"difference": "delta_2", "q_221": "q_2"})
        .merge(
            high_genes[["gene", "difference", "q_221"]].rename(
                columns={"difference": "delta_3", "q_221": "q_3"}
            ),
            on="gene",
            validate="one_to_one",
        )
    )
    ax.scatter(
        100 * effects.delta_2, 100 * effects.delta_3, s=23, color=PALE, alpha=0.70
    )
    for gene, color, offset in (("MYOD1", RED, (8, 6)), ("CKB", BLUE, (8, -12))):
        row = effects.set_index("gene").loc[gene]
        ax.scatter(100 * row.delta_2, 100 * row.delta_3, s=72, color=color, zorder=3)
        ax.annotate(
            gene,
            (100 * row.delta_2, 100 * row.delta_3),
            xytext=offset,
            textcoords="offset points",
            color=color,
            fontsize=8.5,
        )
    ax.axhline(0, color=MUTED, lw=0.8, ls="--")
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_title(
        "C  All 221 gene-region ATAC effects",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.set_xlabel(">=2 vs 1 difference (percentage points)")
    ax.set_ylabel(">=3 vs 1 difference (percentage points)")
    ax.text(
        0.02,
        0.97,
        f"Whole-gene-family minimum BH q: {effects.q_2.min():.3f} / {effects.q_3.min():.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8.1,
        color=MUTED,
    )

    ax = fig.add_subplot(grid[1, 1])
    rank = (
        low_rank[["gene", "n_linked_coq_3of4"]]
        .rename(columns={"n_linked_coq_3of4": "links_2"})
        .merge(
            high_rank[["gene", "n_linked_coq_3of4"]].rename(
                columns={"n_linked_coq_3of4": "links_3"}
            ),
            on="gene",
            validate="one_to_one",
        )
    )
    rank["max_links"] = rank[["links_2", "links_3"]].max(axis=1)
    shown = rank.sort_values(
        ["max_links", "links_3", "gene"], ascending=[False, False, True]
    ).head(10)
    shown = shown.iloc[::-1]
    y = np.arange(len(shown))
    ax.barh(y - 0.19, shown.links_2, height=0.34, color=BLUE, label=">=2 vs 1")
    ax.barh(y + 0.19, shown.links_3, height=0.34, color=RED, label=">=3 vs 1")
    ax.set_yticks(y, shown.gene)
    ax.set_xlim(0, 5.7)
    ax.set_xlabel("Linked peaks with ATAC rise in >=3/4 libraries (count)")
    ax.set_title(
        "D  Exploratory ranking across both contrasts",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=7.7, loc="lower right")
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for panel in fig.axes:
        panel.spines["top"].set_visible(False)
        panel.spines["right"].set_visible(False)
    fig.text(
        0.09,
        0.035,
        "A-C: complete test families at each count threshold. D: same-data ranking; RNA links learned in >=2 vs 1 nuclei. "
        "MYOD1 ranks 1st / 2nd; the >=3 MYOD1 locus effect is examined in Figure 3.",
        fontsize=7.8,
        color=MUTED,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(out.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
