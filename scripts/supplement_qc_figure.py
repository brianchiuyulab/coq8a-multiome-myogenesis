"""Supplementary QC: TSS enrichment and paired RNA/ATAC depth balance."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INK = "#202935"
MUTED = "#687582"
BLUE = "#2166ac"
LIBRARIES = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]


def make_figure(tables: Path, out: Path) -> None:
    qc = pd.read_csv(tables / "qc_nuclei.tsv.gz", sep="\t")
    pairs = pd.read_csv(tables / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "2plus_vs_1")]
    if len(pairs) != 958:
        raise ValueError("QC supplement must use the 958 primary matched pairs")
    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}
    )
    fig, axes = plt.subplots(1, 3, figsize=(12.9, 4.2))
    for gsm in LIBRARIES:
        vals = qc.loc[qc.gsm == gsm, "tss_enrichment"].to_numpy()
        axes[0].hist(
            vals,
            bins=np.linspace(0, 12, 61),
            density=True,
            histtype="step",
            lw=1.4,
            label=gsm[-2:],
        )
    axes[0].axvline(3, color=INK, ls="--", lw=1)
    axes[0].set_xlabel("ATAC TSS enrichment")
    axes[0].set_ylabel("Nucleus density")
    axes[0].legend(frameon=False, fontsize=7, title="Library")
    axes[0].set_title("A  ATAC quality", loc="left", fontweight="bold", fontsize=10)

    indexed = qc.set_index(["gsm", "barcode"])
    for ax, column, title in [
        (axes[1], "total_rna_umi", "B  RNA depth after matching"),
        (axes[2], "total_open_peaks", "C  ATAC depth after matching"),
    ]:
        high = np.log1p(
            np.array(
                [
                    indexed.at[(g, b), column]
                    for g, b in zip(pairs.gsm, pairs.high_barcode)
                ]
            )
        )
        low = np.log1p(
            np.array(
                [
                    indexed.at[(g, b), column]
                    for g, b in zip(pairs.gsm, pairs.low_barcode)
                ]
            )
        )
        ax.scatter(low, high, s=4, alpha=0.22, color=BLUE)
        bounds = [min(low.min(), high.min()), max(low.max(), high.max())]
        ax.plot(bounds, bounds, "--", color=MUTED, lw=0.9)
        ax.set_xlabel("COQ8A low: log1p depth")
        ax.set_ylabel("COQ8A high: log1p depth")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        "Joint RNA/ATAC quality and matching", fontsize=13, fontweight="bold", y=1.05
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(out.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
