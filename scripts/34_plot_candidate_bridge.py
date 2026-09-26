"""Plot selected candidate loci alongside public and optionally private RNA.

This figure displays explicit illustrative loci, not a significance-filtered
discovery set. Private passage data and the resulting figure stay local.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LOCUS = {
    "CSRP3": "chr11:19201752-19202603",
    "CAV3": "chr3:8733438-8733968",
    "LDB3": "chr10:86674717-86675570",
    "TNNC1": "chr3:52446715-52447619",
    "SPARC": "chr5:151698875-151699764",
    "GNAO1": "chr16:56341742-56342652",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    all_rows = pd.read_csv(args.evidence, sep="\t")
    selected = (
        pd.concat(
            [
                all_rows.loc[all_rows.gene.eq(g) & all_rows.peak.eq(p)]
                for g, p in LOCUS.items()
            ]
        )
        .set_index("gene")
        .loc[list(LOCUS)]
    )
    assert len(selected) == 6 and selected.index.is_unique
    rna = pd.read_csv(args.tables / "gene_effects_pooled.tsv", sep="\t")
    rna = (
        rna.loc[
            rna.gate.eq("TSS_ge_3")
            & rna.contrast.eq("3plus_vs_1")
            & rna.modality.eq("RNA")
        ]
        .set_index("gene")
        .loc[list(LOCUS)]
    )
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42}
    )
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(13.6, 6.3),
        gridspec_kw={"width_ratios": [1.3, 1.3, 1]},
        layout="constrained",
    )
    y = np.arange(len(selected))
    palette = {"Early": "#287D8E", "Middle": "#D58B25", "Closing": "#876BA5"}
    colors = [palette[p] for p in selected.phase]
    for ax in axes:
        ax.set_ylim(len(selected) - 0.45, -0.65)
        ax.set_yticks(y)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.15)
        ax.set_axisbelow(True)
    ax = axes[0]
    for i, (_, row) in enumerate(selected.iterrows()):
        ax.plot([row.low_open_pct, row.high_open_pct], [i, i], color=colors[i], lw=2)
        ax.scatter(
            row.low_open_pct, i, color="white", edgecolor=colors[i], s=55, zorder=3
        )
        ax.scatter(row.high_open_pct, i, color=colors[i], s=55, zorder=3)
        ax.text(
            43,
            i,
            f"{row.fold_open:.2f}x\np={row.p_exact:.3g}\nq={row.q_410_peaks:.3g}",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(-1, 58)
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.set_yticklabels([g + "\n" + selected.loc[g, "phase"] for g in LOCUS])
    ax.set_xlabel("Nuclei with accessible peak (%)")
    ax.set_title(
        "A  Public multiome: ATAC\nOpen circle: low | Filled: high",
        loc="left",
        fontsize=11,
    )
    ax = axes[1]
    ax.axvline(0, color="#999999", lw=0.8)
    for i, (g, row) in enumerate(rna.iterrows()):
        if selected.loc[g, "n_gene_detected"] == 0:
            ax.text(
                0.02, i, "RNA not detected", va="center", fontsize=9, color="#666666"
            )
            continue
        ax.errorbar(
            row.difference,
            i,
            xerr=[[row.difference - row.ci_low], [row.ci_high - row.difference]],
            fmt="o",
            color=colors[i],
            capsize=3,
            markersize=5,
        )
        ax.text(
            0.40, i, f"p={row.p_pair:.3g}\nq={row.q_221:.3g}", va="center", fontsize=9
        )
    ax.set_xlim(-0.15, 0.61)
    ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
    ax.set_yticklabels([])
    ax.set_title(
        "B  Public multiome: RNA\nHigh minus low; 95% paired CI",
        loc="left",
        fontsize=11,
    )
    ax.set_xlabel("Mean log1p(CP10K) difference")
    ax = axes[2]
    ax.axvline(0, color="#999999", lw=0.8)
    for day, offset, marker in [(2, -0.12, "o"), (6, 0.12, "s")]:
        v = selected[f"P33_vs_P11_at_D{day}_log2FC"].to_numpy()
        ax.scatter(v, y + offset, color=colors, marker=marker, s=38, label=f"D{day}")
    ax.set_xticks([-4, -2, 0, 2], ["0.0625", "0.25", "1", "4"])
    ax.set_xlim(-4.7, 2.4)
    ax.set_yticklabels([])
    ax.set_xlabel("RNA fold change (log2 axis)")
    ax.set_title(
        "C  C2C12 RNA: P33 / P11\nCircle: D2 | Square: D6", loc="left", fontsize=11
    )
    fig.suptitle(
        "Candidate loci and RNA associations\nCOQ8A >=3 vs 1 UMI | TSS enrichment >=3 | 201 matched pairs",
        fontsize=14,
    )
    for suffix in ["png", "pdf"]:
        fig.savefig(
            args.out / f"Figure_candidate_bridge.{suffix}", dpi=180, facecolor="white"
        )
    plt.close(fig)
    selected.to_csv(args.out / "candidate_bridge_source.tsv", sep="\t", na_rep="NA")


if __name__ == "__main__":
    main()
