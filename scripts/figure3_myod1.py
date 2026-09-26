"""Figure 3: the nominated MYOD1 region and all key sensitivity challenges."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INK = "#202935"
MUTED = "#687582"
BLUE = "#2166ac"
RED = "#b2182b"
PALE = "#aeb9c5"
GRID = "#e4e9ee"
SIX = "positive_link_nonMRF_tss2_6"
LIBRARIES = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]
TSS_KB = 17719.565  # GENCODE v48 MYOD1 gene-feature TSS on hg38 chr11.


def load_locus(tables: Path):
    candidates = pd.read_csv(tables / "candidate_peak_gene.tsv.gz", sep="\t")
    candidates = candidates[
        (candidates.gene == "MYOD1") & (candidates.n_libraries == 4)
    ]
    effects = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    effects = effects[(effects.gate == "TSS_ge_3") & (effects.contrast == "3plus_vs_1")]
    links = pd.read_csv(tables / "peak_gene_links_tss2.tsv", sep="\t")
    links = links[
        (links.gene == "MYOD1")
        & (links.q_all_links < 0.05)
        & (links.partial_r > 0)
        & (links.n_libraries == 4)
    ]
    selected = set(links.loc[links.mrf_max_score < 0.95, "peak"])
    locus = candidates.merge(
        effects[["peak", "delta_pp"]], on="peak", validate="one_to_one"
    )
    coordinates = locus.peak.str.extract(r"^chr11:(\d+)-(\d+)$").astype(int)
    locus["mid_kb"] = (coordinates[0] + coordinates[1]) / 2000
    if len(locus) != 19 or len(selected) != 6 or not selected.issubset(set(locus.peak)):
        raise ValueError("MYOD1 common-peak or linked-six set changed")
    return locus.sort_values("mid_kb"), selected


def make_figure(tables: Path, out: Path) -> None:
    locus, selected = load_locus(tables)
    summary = pd.read_csv(tables / "myod1_locus_summary.tsv", sep="\t")
    by_library = pd.read_csv(tables / "myod1_locus_by_library.tsv", sep="\t")
    doublet = pd.read_csv(tables / "doublet_tail_sensitivity.tsv", sep="\t")
    primary = summary[
        (summary.gate == "TSS_ge_3")
        & (summary.contrast == "3plus_vs_1")
        & (summary.region_set == SIX)
    ].iloc[0]
    if (
        not np.isclose(primary.fold_open, 1.28, atol=0.001)
        or int(primary.n_pairs) != 201
    ):
        raise ValueError("Six-peak sensitivity result changed")

    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}
    )
    fig = plt.figure(figsize=(15.2, 9.2), facecolor="white")
    grid = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.0, 1.15],
        width_ratios=[1.45, 1, 1],
        left=0.09,
        right=0.96,
        top=0.86,
        bottom=0.12,
        hspace=0.44,
        wspace=0.38,
    )
    fig.suptitle(
        "Exploratory MYOD1 locus: which peaks, how large, how stable?",
        x=0.09,
        y=0.965,
        ha="left",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.09,
        0.913,
        f"Six selected peaks: {primary.high_open_pct:.2f}% vs {primary.low_open_pct:.2f}% open, "
        f"ratio {primary.fold_open:.2f}x in {int(primary.n_pairs)} pairs; paired p={primary.p_pair:.3f}. "
        "In B: blue >=2 UMI, red >=3 UMI.",
        color=MUTED,
        fontsize=9.5,
    )

    ax = fig.add_subplot(grid[0, :])
    colors = [
        RED if row.peak in selected else (BLUE if row.delta_pp < 0 else PALE)
        for row in locus.itertuples()
    ]
    ax.bar(locus.mid_kb, locus.delta_pp, width=0.75, color=colors, alpha=0.9)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.axvline(TSS_KB, color=INK, ls="--", lw=0.8)
    ax.text(TSS_KB + 0.6, ax.get_ylim()[1] * 0.85, "MYOD1 TSS", fontsize=8, color=INK)
    ax.set_ylabel("High - low open fraction\n(percentage points)")
    ax.set_xlabel("hg38 chr11 coordinate (kb)")
    ax.set_title(
        "A  All 19 four-library-common peaks near MYOD1",
        loc="left",
        fontweight="bold",
        fontsize=11,
    )
    ax.text(
        0.02,
        0.04,
        "Red: six positive peak-RNA links with MRF PWM score <0.95 at the TSS>=2 link gate; "
        "blue: lower opening; grey: other peaks.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
    )

    ax = fig.add_subplot(grid[1, 0])
    subset = summary[summary.gate == "TSS_ge_3"]
    names = [
        "all_candidate_32",
        "common_all4_19",
        "positive_link_10",
        "positive_link_nonMRF_5",
        SIX,
    ]
    labels = [
        "32 nearby",
        "19 shared",
        "10 links (TSS>=3)",
        "5 lower-MRF-score links",
        "6 lower-MRF-score links",
    ]
    for offset, contrast, color, legend in [
        (-0.09, "2plus_vs_1", BLUE, "COQ8A >=2 vs 1"),
        (+0.09, "3plus_vs_1", RED, "COQ8A >=3 vs 1"),
    ]:
        rows = subset[subset.contrast == contrast].set_index("region_set").loc[names]
        ax.scatter(
            rows.fold_open,
            np.arange(5) + offset,
            s=44,
            color=color,
            label=legend,
            zorder=3,
        )
    ax.axvline(1, color=MUTED, ls="--", lw=0.8)
    ax.set_yticks(range(5), labels)
    ax.invert_yaxis()
    ax.set_xlim(0.97, 1.35)
    ax.set_xlabel("High / low open-fraction ratio")
    ax.set_title(
        "B  Region and count thresholds", loc="left", fontweight="bold", fontsize=11
    )
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    ax = fig.add_subplot(grid[1, 1])
    rows = (
        by_library[
            (by_library.gate == "TSS_ge_3")
            & (by_library.contrast == "3plus_vs_1")
            & (by_library.region_set == SIX)
        ]
        .set_index("gsm")
        .loc[LIBRARIES]
    )
    y = np.arange(4)
    ax.errorbar(
        rows.delta_pp,
        y,
        xerr=[rows.delta_pp - 100 * rows.ci_low, 100 * rows.ci_high - rows.delta_pp],
        fmt="o",
        capsize=3,
        color=RED,
        markersize=5,
        lw=1.1,
    )
    ax.axvline(0, color=MUTED, ls="--", lw=0.8)
    ax.set_yticks(
        y,
        [
            f"L{i//2+1} {'stem' if i%2==0 else 'diff'}  ({int(row.n_pairs)})"
            for i, row in enumerate(rows.itertuples())
        ],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Six-peak high - low (percentage points)")
    ax.set_title(
        "C  Four library directions", loc="left", fontweight="bold", fontsize=11
    )
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)

    ax = fig.add_subplot(grid[1, 2])
    challenge = doublet[
        (doublet.contrast == "3plus_vs_1") & (doublet.region_set == SIX)
    ].sort_values("scrublet_tail_removed_pct")
    ax.plot(
        challenge.scrublet_tail_removed_pct,
        challenge.fold_open,
        marker="o",
        color=RED,
        lw=1.5,
        ms=5,
    )
    ax.axhline(1, color=MUTED, ls="--", lw=0.8)
    ax.set_ylim(1.0, 1.34)
    ax.set_xticks(challenge.scrublet_tail_removed_pct)
    ax.set_xlabel("Highest RNA doublet-score tail removed (%)")
    ax.set_ylabel("Six-peak high / low ratio")
    ax.set_title(
        "D  Doublet-score challenge", loc="left", fontweight="bold", fontsize=11
    )
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for row in challenge.itertuples():
        ax.annotate(
            f"n={row.n_pairs}",
            (row.scrublet_tail_removed_pct, row.fold_open),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            fontsize=7,
            color=MUTED,
        )

    for panel in fig.axes:
        panel.spines["top"].set_visible(False)
        panel.spines["right"].set_visible(False)
    fig.text(
        0.09,
        0.035,
        "MYOD1 is a same-data exploratory nomination. Linked subsets were learned in overlapping nuclei; "
        "pair-level p and intervals are descriptive, not independent-source replication (source lines n=2).",
        fontsize=8,
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
