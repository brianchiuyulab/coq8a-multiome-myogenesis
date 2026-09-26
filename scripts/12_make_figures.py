"""Generate publication-style vector PDF and high-resolution PNG figures."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tss_story_figure import make_story_figure

BLUE = "#2166ac"
RED = "#b2182b"
GREEN = "#1b7837"
INK = "#202935"
GRAY = "#687582"
LIGHT = "#e5e9ed"
LIBRARIES = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(path.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)


def label(ax, letter, title):
    ax.text(
        -0.08,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
    )
    ax.set_title(title, loc="left", fontsize=10, pad=12, fontweight="bold")


def figure_one(t, out):
    """Show fixed 221-gene search, TSS fragments, and exploratory MYOD1 nomination."""
    make_story_figure(t, out / "Figure_1_global_discovery")


def supplementary_myod1(t, out):
    locus = pd.read_csv(t / "myod1_locus_summary.tsv", sep="\t")
    sample = pd.read_csv(t / "myod1_locus_by_library.tsv", sep="\t")
    peak = pd.read_csv(t / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peak = peak[(peak.gate == "TSS_ge_3") & (peak.contrast == "2plus_vs_1")]
    links = pd.read_csv(t / "peak_gene_links_tss2.tsv", sep="\t")
    region = pd.read_csv(t / "candidate_peak_gene.tsv.gz", sep="\t")
    region = region[(region.gene == "MYOD1") & (region.n_libraries == 4)]
    d = region.merge(peak[["peak", "delta_pp"]], on="peak")
    coord = d.peak.str.extract(r":(\d+)-(\d+)").astype(float)
    d["mid_kb"] = (coord[0] + coord[1]) / 2000
    linked = set(
        links[
            (links.gene == "MYOD1")
            & (links.q_all_links < 0.05)
            & (links.partial_r > 0)
            & (links.n_libraries == 4)
        ].peak
    )
    nonmrf = set(
        links[
            (links.gene == "MYOD1")
            & (links.q_all_links < 0.05)
            & (links.partial_r > 0)
            & (links.n_libraries == 4)
            & (links.mrf_max_score < 0.95)
        ].peak
    )
    fig = plt.figure(figsize=(12.6, 10.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.1], hspace=0.43, wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    label(ax, "A", "MYOD1 locus: all 19 four-library common peaks")
    colors = [RED if v > 0 else BLUE for v in d.delta_pp]
    ax.bar(d.mid_kb, d.delta_pp, width=0.85, color=colors, alpha=0.75)
    for r in d.itertuples():
        if r.peak in linked:
            ax.scatter(
                r.mid_kb,
                r.delta_pp,
                marker="o",
                facecolors="none",
                edgecolors=INK,
                s=93,
                lw=1.1,
                zorder=3,
            )
        if r.peak in nonmrf:
            ax.scatter(r.mid_kb, r.delta_pp, marker="*", color=INK, s=28, zorder=4)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.axvline(17719.565, color=GREEN, lw=1.3, linestyle="--")
    ax.text(
        17719.565,
        ax.get_ylim()[1] * 0.93,
        "MYOD1 TSS",
        rotation=90,
        ha="right",
        va="top",
        fontsize=8,
        color=GREEN,
    )
    ax.set_xlabel("chr11 coordinate (kb, hg38)")
    ax.set_ylabel("COQ8A ≥2 vs 1 UMI\nopen-fraction difference (pp)")
    ax.text(
        0.01,
        0.02,
        "Outline: TSS≥2 peak–RNA link  •  Star: maximum MRF PWM score <0.95",
        transform=ax.transAxes,
        fontsize=8,
        color=GRAY,
    )
    ax = fig.add_subplot(gs[1, 0])
    label(ax, "B", "Region sets and COQ8A count sensitivity")
    x = locus[locus.gate == "TSS_ge_3"]
    names = [
        "all_candidate_32",
        "common_all4_19",
        "positive_link_10",
        "positive_link_nonMRF_5",
        "positive_link_nonMRF_tss2_6",
    ]
    pretty = [
        "All 32 nearby",
        "19 shared",
        "10 linked (TSS≥3)",
        "5 MRF score<0.95 (TSS≥3 links)",
        "6 MRF score<0.95 (TSS≥2 links)",
    ]
    for j, contrast in enumerate(["2plus_vs_1", "3plus_vs_1"]):
        z = x[x.contrast == contrast].set_index("region_set").loc[names]
        y = np.arange(len(names)) + (j - 0.5) * 0.18
        ax.scatter(
            z.fold_open,
            y,
            s=46,
            color=[BLUE, RED][j],
            label=["COQ ≥2 vs 1", "COQ ≥3 vs 1"][j],
        )
        for k, r in enumerate(z.itertuples()):
            ax.text(
                r.fold_open + 0.008,
                y[k],
                f"p={r.p_pair:.3f}",
                fontsize=7.5,
                va="center",
                color=GRAY,
            )
    ax.axvline(1, color=GRAY, lw=0.8, linestyle="--")
    ax.set_yticks(range(len(names)), pretty)
    ax.invert_yaxis()
    ax.set_xlabel("High / low mean open fraction")
    ax.set_xlim(0.98, 1.37)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="x", color=LIGHT, lw=0.5)
    ax = fig.add_subplot(gs[1, 1])
    label(ax, "C", "Extreme contrast by library")
    s = (
        sample[
            (sample.gate == "TSS_ge_3")
            & (sample.contrast == "3plus_vs_1")
            & (sample.region_set == "positive_link_nonMRF_tss2_6")
        ]
        .set_index("gsm")
        .loc[LIBRARIES]
    )
    xx = np.arange(4)
    ax.errorbar(
        xx,
        s.delta_pp,
        yerr=[s.delta_pp - 100 * s.ci_low, 100 * s.ci_high - s.delta_pp],
        fmt="o",
        color=RED,
        capsize=3,
        markersize=5,
        lw=1,
    )
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xticks(
        xx,
        [
            f'{s.loc[gsm, "source"].replace("line", "L")} {s.loc[gsm, "stage"][:4]}\n'
            f'{s.loc[gsm, "n_pairs"]} pairs'
            for gsm in LIBRARIES
        ],
        fontsize=7,
    )
    ax.set_ylabel("6-region difference (percentage points)")
    fig.suptitle(
        "Exploratory MYOD1 locus and threshold sensitivity",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        -0.015,
        "Linked sets were learned in the same nuclei; all significance is exploratory at nucleus level.",
        ha="center",
        fontsize=8,
        color=GRAY,
    )
    save(fig, out / "Supplementary_Figure_MYOD1_Exploratory")


def supplementary_qc(t, out):
    qc = pd.read_csv(t / "qc_nuclei.tsv.gz", sep="\t")
    summary = pd.read_csv(t / "qc_by_library.tsv", sep="\t")
    pairs = pd.read_csv(t / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "2plus_vs_1")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    for i, gsm in enumerate(LIBRARIES):
        vals = qc[qc.gsm == gsm].tss_enrichment.to_numpy()
        axes[0].hist(
            vals,
            bins=np.linspace(0, 12, 61),
            density=True,
            histtype="step",
            lw=1.4,
            label=gsm[-2:],
        )
    axes[0].axvline(3, color=INK, linestyle="--", lw=1)
    axes[0].set_xlabel("TSS enrichment")
    axes[0].set_ylabel("Nucleus density")
    axes[0].legend(frameon=False, title="Library", fontsize=7)
    label(axes[0], "A", "ATAC quality by library")
    d = qc.set_index(["gsm", "barcode"])
    for ax, col, letter, title in [
        (axes[1], "total_rna_umi", "B", "RNA depth after matching"),
        (axes[2], "total_open_peaks", "C", "ATAC depth after matching"),
    ]:
        high = np.log1p(
            np.array([d.at[(g, b), col] for g, b in zip(pairs.gsm, pairs.high_barcode)])
        )
        low = np.log1p(
            np.array([d.at[(g, b), col] for g, b in zip(pairs.gsm, pairs.low_barcode)])
        )
        ax.scatter(low, high, s=3, alpha=0.20, color=BLUE)
        bounds = [min(low.min(), high.min()), max(low.max(), high.max())]
        ax.plot(bounds, bounds, "--", color=GRAY, lw=0.9)
        ax.set_xlabel("Low COQ8A: log1p depth")
        ax.set_ylabel("High COQ8A: log1p depth")
        label(ax, letter, title)
    fig.suptitle(
        "Joint RNA/ATAC quality and paired-depth balance",
        fontsize=13,
        fontweight="bold",
        y=1.05,
    )
    save(fig, out / "Supplementary_Figure_QC")


def figure_two(t, out):
    peaks = pd.read_csv(t / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peaks = peaks[
        (peaks.gate == "TSS_ge_3")
        & (peaks.contrast == "2plus_vs_1")
        & peaks.q_candidate_peaks.notna()
    ]
    genes = pd.read_csv(t / "gene_effects_pooled.tsv", sep="\t")
    genes = genes[
        (genes.gate == "TSS_ge_3")
        & (genes.contrast == "2plus_vs_1")
        & (genes.modality == "ATAC")
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    ax = axes[0]
    ax.scatter(
        peaks.delta_pp,
        -np.log10(np.maximum(peaks.p_pair_binomial, 1e-300)),
        color="#aeb9c5",
        s=7,
        alpha=0.45,
    )
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("COQ8A high − low open fraction (percentage points)")
    ax.set_ylabel("−log10 paired peak p")
    ax.text(
        0.02,
        0.97,
        f"{len(peaks):,} peaks tested; minimum BH q={peaks.q_candidate_peaks.min():.3f}",
        transform=ax.transAxes,
        fontsize=8,
        va="top",
        color=GRAY,
    )
    label(ax, "A", "All eligible candidate peaks")
    ax = axes[1]
    ax.scatter(
        100 * genes.difference,
        -np.log10(np.maximum(genes.p_pair, 1e-300)),
        color="#aeb9c5",
        s=18,
        alpha=0.65,
    )
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("COQ8A high − low nearby-peak mean (percentage points)")
    ax.set_ylabel("−log10 paired gene-region p")
    ax.text(
        0.02,
        0.97,
        f"{len(genes)} genes tested; minimum BH q={genes.q_221.min():.3f}",
        transform=ax.transAxes,
        fontsize=8,
        va="top",
        color=GRAY,
    )
    label(ax, "B", "All 221 gene regions")
    fig.suptitle(
        "Complete primary ATAC search space (TSS≥3, COQ8A≥2 vs 1 UMI)",
        fontsize=13,
        fontweight="bold",
        y=1.05,
    )
    save(fig, out / "Figure_2_Global_ATAC_Scan")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", type=Path, required=True)
    p.add_argument("--figures", type=Path, required=True)
    a = p.parse_args()
    style()
    figure_one(a.tables, a.figures / "main")
    figure_two(a.tables, a.figures / "main")
    supplementary_qc(a.tables, a.figures / "supplement")
    supplementary_myod1(a.tables, a.figures / "supplement")
    print("Wrote two primary figures and two supplements", flush=True)


if __name__ == "__main__":
    main()
