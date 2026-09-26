"""Figure 2: complete primary ATAC test families and exploratory gene ranking.

Panels A and B display every eligible primary test. Panels C and D show the
subsequent same-data integration and ranking, not independent confirmation.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INK = "#202935"
MUTED = "#687582"
PALE = "#adb9c6"
RED = "#b2182b"
GRID = "#e3e8ed"
GATE = "TSS_ge_3"
CONTRAST = "2plus_vs_1"


def load_primary(tables: Path):
    peaks = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peaks = peaks[
        (peaks.gate == GATE)
        & (peaks.contrast == CONTRAST)
        & peaks.q_candidate_peaks.notna()
    ].copy()
    genes = pd.read_csv(tables / "gene_effects_pooled.tsv", sep="\t")
    genes = genes[
        (genes.gate == GATE) & (genes.contrast == CONTRAST) & (genes.modality == "ATAC")
    ].copy()
    rna = pd.read_csv(tables / "gene_effects_pooled.tsv", sep="\t")
    rna = rna[
        (rna.gate == GATE) & (rna.contrast == CONTRAST) & (rna.modality == "RNA")
    ].copy()
    ranking = pd.read_csv(tables / "candidate_gene_ranking.tsv", sep="\t")
    if (
        len(peaks) != 6_871
        or len(genes) != 221
        or len(rna) != 221
        or len(ranking) != 221
    ):
        raise ValueError(
            "Primary screening or ranking table has incomplete test families"
        )
    if set(genes.gene) != set(ranking.gene):
        raise ValueError("Ranked genes differ from the 221 primary gene regions")
    if peaks.q_candidate_peaks.lt(0.05).any() or genes.q_221.lt(0.05).any():
        raise ValueError("Figure caption assumes no primary ATAC FDR-positive result")
    return peaks, genes, rna, ranking


def make_figure(tables: Path, out: Path) -> None:
    peaks, genes, rna, ranking = load_primary(tables)
    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}
    )
    fig = plt.figure(figsize=(13.6, 8.5), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.05, 0.95],
        hspace=0.53,
        wspace=0.35,
        left=0.09,
        right=0.97,
        top=0.86,
        bottom=0.12,
    )
    fig.suptitle(
        "Myogenesis-space screen and candidate nomination",
        x=0.09,
        y=0.965,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.09,
        0.917,
        "COQ8A >=2 versus 1 UMI  |  TSS enrichment >=3  |  958 within-library matched pairs",
        fontsize=9.5,
        color=MUTED,
    )

    ax = fig.add_subplot(grid[0, 0])
    ax.scatter(
        peaks.delta_pp,
        -np.log10(peaks.p_pair_binomial.clip(lower=1e-300)),
        s=5,
        alpha=0.38,
        color=PALE,
        rasterized=True,
    )
    ax.axvline(0, color=MUTED, lw=0.85, ls="--")
    ax.set_xlabel("Peak opening difference (percentage points)")
    ax.set_ylabel("-log10 paired peak p")
    ax.set_title(
        "A  All 6,871 eligible peaks", loc="left", fontweight="bold", fontsize=11
    )
    ax.text(
        0.02,
        0.97,
        f"Minimum BH q = {peaks.q_candidate_peaks.min():.3f}; 0 peaks at q<0.05",
        transform=ax.transAxes,
        va="top",
        color=MUTED,
        fontsize=8.5,
    )

    ax = fig.add_subplot(grid[0, 1])
    ax.scatter(
        100 * genes.difference,
        -np.log10(genes.p_pair.clip(lower=1e-300)),
        s=22,
        alpha=0.65,
        color=PALE,
    )
    myod = genes.set_index("gene").loc["MYOD1"]
    ax.scatter(
        100 * myod.difference,
        -np.log10(myod.p_pair),
        s=70,
        color=RED,
        zorder=3,
    )
    ax.annotate(
        f"MYOD1\nq={myod.q_221:.3f}",
        (100 * myod.difference, -np.log10(myod.p_pair)),
        xytext=(10, 7),
        textcoords="offset points",
        fontsize=8.5,
        color=RED,
    )
    ax.axvline(0, color=MUTED, lw=0.85, ls="--")
    ax.set_xlabel("Gene-region ATAC difference (percentage points)")
    ax.set_ylabel("-log10 paired gene-region p")
    ax.set_title("B  All 221 gene regions", loc="left", fontweight="bold", fontsize=11)
    ax.text(
        0.02,
        0.97,
        f"Minimum BH q = {genes.q_221.min():.3f}; 0 regions at q<0.05",
        transform=ax.transAxes,
        va="top",
        color=MUTED,
        fontsize=8.5,
    )

    ax = fig.add_subplot(grid[1, 0])
    effects = (
        rna[["gene", "difference"]]
        .rename(columns={"difference": "rna_delta"})
        .merge(
            genes[["gene", "difference"]].rename(columns={"difference": "atac_delta"}),
            on="gene",
            validate="one_to_one",
        )
    )
    ax.scatter(
        effects.rna_delta, 100 * effects.atac_delta, s=19, alpha=0.65, color=PALE
    )
    myod_effect = effects.set_index("gene").loc["MYOD1"]
    ax.scatter(
        myod_effect.rna_delta, 100 * myod_effect.atac_delta, s=70, color=RED, zorder=3
    )
    ax.annotate(
        "MYOD1",
        (myod_effect.rna_delta, 100 * myod_effect.atac_delta),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=8.5,
        color=RED,
    )
    ax.axhline(0, color=MUTED, lw=0.8, ls="--")
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("Gene RNA difference (log1p CP10K)")
    ax.set_ylabel("Nearby ATAC difference (percentage points)")
    ax.set_title(
        "C  RNA and ATAC effects for all 221 genes",
        loc="left",
        fontweight="bold",
        fontsize=11,
    )

    ax = fig.add_subplot(grid[1, 1])
    shown = ranking.head(10).iloc[::-1].copy()
    colours = [RED if gene == "MYOD1" else PALE for gene in shown.gene]
    ax.barh(shown.gene, shown.n_linked_coq_3of4, height=0.7, color=colours)
    for index, row in enumerate(shown.itertuples()):
        ax.text(
            row.n_linked_coq_3of4 + 0.08,
            index,
            str(row.n_linked_coq_3of4),
            va="center",
            fontsize=8,
            color=INK,
        )
    ax.set_xlim(0, max(shown.n_linked_coq_3of4) + 0.65)
    ax.set_xlabel("Positive RNA links with ATAC rise in >=3/4 libraries (count)")
    ax.set_title(
        "D  Exploratory ranking - top 10 of 221",
        loc="left",
        fontweight="bold",
        fontsize=11,
    )
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8)
    for panel in fig.axes:
        panel.spines["top"].set_visible(False)
        panel.spines["right"].set_visible(False)
    fig.text(
        0.09,
        0.038,
        "A-B: complete test families; points show nominal paired p and labels show whole-family BH q.  "
        "C-D: same-data integration and nomination for Figure 3. Full 221-row ranking is released; source lines n=2.",
        color=MUTED,
        fontsize=8,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(out.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
