"""Display the 221-gene TSS context using shared signal scales and row order."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.tables / "tss_fragment_profile_221.tsv.gz", sep="\t")
    assert len(data) == 22100 and data.n_high.eq(201).all() and data.n_low.eq(201).all()
    order = (
        data.groupby("gene")[["high_cuts", "low_cuts"]]
        .sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
    )
    matrix = data.pivot(
        index="gene", columns="bin_start_bp", values=["low_cuts", "high_cuts"]
    )
    low = gaussian_filter1d(
        matrix["low_cuts"].loc[order].to_numpy() / 201 * 100, 1, axis=1, mode="nearest"
    )
    high = gaussian_filter1d(
        matrix["high_cuts"].loc[order].to_numpy() / 201 * 100, 1, axis=1, mode="nearest"
    )
    limit = np.quantile(np.r_[low.ravel(), high.ravel()], 0.995)
    pos = (matrix["low_cuts"].columns.to_numpy() + 50) / 1000
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42}
    )
    fig = plt.figure(figsize=(8.3, 8.6), layout="constrained")
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 3.3, 0.13])
    profile = fig.add_subplot(grid[0, :])
    profile.plot(pos, low.mean(axis=0), color="#A0B5C9", lw=2, label="COQ8A low: 1 UMI")
    profile.plot(
        pos, high.mean(axis=0), color="#174E83", lw=2, label="COQ8A high: >=3 UMI"
    )
    profile.set_ylabel("Mean insertions\nper 100 nuclei / 100 bp")
    profile.spines[["top", "right"]].set_visible(False)
    profile.set_xlim(-5, 5)
    profile.set_xticks([-5, 0, 5], ["-5 kb", "TSS", "+5 kb"])
    profile.legend(frameon=False, fontsize=9)
    for j, (values, title) in enumerate([(low, "COQ8A low"), (high, "COQ8A high")]):
        ax = fig.add_subplot(grid[1, j])
        im = ax.imshow(
            values,
            cmap="Blues",
            vmin=0,
            vmax=limit,
            aspect="auto",
            interpolation="nearest",
            extent=[-5, 5, 221, 0],
        )
        ax.axvline(0, color="#666666", ls=":", lw=0.7)
        ax.set_xticks([-5, 0, 5], ["-5 kb", "TSS", "+5 kb"])
        ax.set_xlabel("Distance from TSS")
        ax.set_title(title)
        if j == 0:
            ax.set_ylabel("221 myogenesis genes")
            ax.set_yticks(
                [0.5, 49.5, 99.5, 149.5, 220.5], ["1", "50", "100", "150", "221"]
            )
        else:
            ax.set_yticks([])
    cax = fig.add_subplot(grid[2, :])
    fig.colorbar(
        im,
        cax=cax,
        orientation="horizontal",
        extend="max",
        label="ATAC insertions per 100 nuclei / 100 bp",
    )
    fig.suptitle(
        "Myogenesis: TSS-aligned chromatin accessibility\n201 depth-matched pairs | TSS enrichment >=3",
        fontsize=14,
    )
    for suffix in ["png", "pdf"]:
        fig.savefig(
            args.out / f"Figure_221_gene_TSS_heatmap.{suffix}",
            dpi=180,
            facecolor="white",
        )
    pd.DataFrame({"row": np.arange(1, 222), "gene": order}).to_csv(
        args.out / "TSS_heatmap_row_order.tsv", sep="\t", index=False
    )
    print(
        "221 rows; common pooled-signal ordering; shared color scale; 100 bp display smoothing."
    )
    print("Whole-window high/low ratio:", data.high_cuts.sum() / data.low_cuts.sum())


if __name__ == "__main__":
    main()
