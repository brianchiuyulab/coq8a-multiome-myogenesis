"""Figure 1: study design and 221-gene TSS-aligned ATAC fragment display."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


INK = "#202935"
MUTED = "#687582"
BLUE = "#2166ac"
RED = "#b2182b"
LIBRARIES = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]
LABELS = [
    "Line 1 - stem",
    "Line 1 - differentiated",
    "Line 2 - stem",
    "Line 2 - differentiated",
]


def read_profile(tables: Path):
    profile = pd.read_csv(tables / "tss_fragment_profile_221.tsv.gz", sep="\t")
    genes = pd.read_csv(tables / "gene_search_space.tsv", sep="\t")
    if len(profile) != 22_100 or set(profile.gene) != set(genes.gene):
        raise ValueError(
            "TSS profile must contain 100 bins for all 221 programme genes"
        )
    if profile.groupby("gene").size().ne(100).any():
        raise ValueError("Incomplete TSS bin grid")
    order = (
        profile.groupby("gene")[["high_cuts", "low_cuts"]]
        .sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index.tolist()
    )
    grid = profile.pivot(
        index="gene", columns="bin_start_bp", values=["low_cuts", "high_cuts"]
    )
    low = grid["low_cuts"].loc[order].to_numpy(dtype=float)
    high = grid["high_cuts"].loc[order].to_numpy(dtype=float)
    return order, low / 201, high / 201


def make_figure(tables: Path, out: Path) -> None:
    order, low, high = read_profile(tables)
    pairs = pd.read_csv(tables / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "3plus_vs_1")]
    sizes = pairs.groupby("gsm").size().reindex(LIBRARIES)
    if sizes.isna().any() or int(sizes.sum()) != 201:
        raise ValueError("Expected 201 main-contrast within-library matched pairs")
    ratio = high.sum() / low.sum()
    low = gaussian_filter1d(low, sigma=1, axis=1, mode="nearest")
    high = gaussian_filter1d(high, sigma=1, axis=1, mode="nearest")
    colour_max = float(np.quantile(np.r_[low.ravel(), high.ravel()], 0.995))

    plt.rcParams.update(
        {"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"}
    )
    fig = plt.figure(figsize=(12.8, 9.2), facecolor="white")
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.7, 1, 3.25],
        left=0.095,
        right=0.945,
        top=0.88,
        bottom=0.10,
        hspace=0.39,
        wspace=0.10,
    )
    fig.suptitle(
        "COQ8A groups and TSS accessibility in a fixed myogenesis programme",
        x=0.095,
        y=0.975,
        ha="left",
        fontsize=15.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.095,
        0.93,
        "Two source lines, four RNA/ATAC libraries | COQ8A >=3 versus 1 UMI | ATAC TSS QC >=3",
        fontsize=9.5,
        color=MUTED,
    )

    ax = fig.add_subplot(grid[0, :])
    ax.axis("off")
    ax.set_title(
        "A  Within-library matched comparison",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )
    for index, (gsm, label) in enumerate(zip(LIBRARIES, LABELS)):
        x = 0.01 + 0.25 * index
        ax.text(
            x,
            0.60,
            label,
            transform=ax.transAxes,
            fontsize=9,
            color=INK,
            fontweight="bold",
        )
        ax.text(
            x,
            0.27,
            f"{int(sizes[gsm])} high-low nucleus pairs",
            transform=ax.transAxes,
            fontsize=8.5,
            color=MUTED,
        )
    ax = fig.add_subplot(grid[1, :])
    position = np.arange(100) / 10 - 4.95
    ax.plot(position, low.mean(axis=0), color=BLUE, lw=1.6, label="COQ8A low")
    ax.plot(position, high.mean(axis=0), color=RED, lw=1.6, label="COQ8A high")
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(-5, 5)
    ax.set_ylabel("Mean Tn5 cuts / nucleus / 100 bp", fontsize=8.5)
    ax.set_title(
        "B  Mean ATAC insertion profile across all 221 gene TSSs",
        loc="left",
        fontsize=11,
        fontweight="bold",
    )
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
    ax.text(
        0.02,
        0.90,
        f"High/low total TSS cut ratio = {ratio:.3f}x",
        transform=ax.transAxes,
        fontsize=8.5,
        color=MUTED,
        va="top",
    )
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    norm = Normalize(vmin=0, vmax=colour_max)
    axes = [fig.add_subplot(grid[2, index]) for index in range(2)]
    for axis, values, label in zip(axes, [low, high], ["C  COQ8A low", "COQ8A high"]):
        image = axis.imshow(
            values, aspect="auto", cmap="Blues", norm=norm, interpolation="nearest"
        )
        axis.set_title(label, fontsize=10, pad=7)
        axis.set_xticks([0, 49.5, 99], ["-5", "TSS", "+5"])
        axis.set_xlabel("Distance from gene TSS (kb)", fontsize=8.5)
        axis.set_yticks([0, 54, 109, 164, 220], ["1", "55", "110", "165", "221"])
        axis.tick_params(length=0, labelsize=8)
        axis.axvline(49.5, color="white", lw=0.6)
        for side in axis.spines.values():
            side.set_visible(False)
    axes[0].set_ylabel("221 genes sorted by pooled TSS signal", fontsize=8.5)
    axes[1].set_yticklabels([])
    fig.colorbar(
        image,
        ax=axes,
        location="bottom",
        fraction=0.03,
        pad=0.10,
        label="Mean Tn5 cuts / nucleus / 100 bp (shared colour scale)",
    )
    fig.text(
        0.095,
        0.025,
        "GENCODE v48 gene TSS +/-5 kb; deduplicated ATAC insertions in 100-bp bins, 100-bp display smoothing. "
        "Rows are ordered by pooled signal, never by COQ8A effect. This descriptive profile does not select candidate loci.",
        fontsize=7.7,
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
