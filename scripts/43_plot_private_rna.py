"""Plot existing C2C12 passage/time RNA estimates to an explicitly private output."""

import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-results", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    if a.out.resolve().is_relative_to(root):
        p.error("Private RNA output must be outside the public repository")
    a.out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(a.gene_results, low_memory=False)
    d["gene"] = d.symbol.str.upper()
    d = d.set_index("gene")
    genes = ["COQ8A", "CSRP3", "CAV3", "MYOD1", "CACNA1H"]
    cols = [f"P{p}_D{day}" for p in [11, 22, 33] for day in [0, 1, 2, 6]]
    data = d.loc[genes, cols].astype(float)
    relative = data.subtract(data.P11_D0, axis=0)
    relative.to_csv(a.out / "candidate_rna_heatmap_source.tsv", sep="\t")
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42}
    )
    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.29, top=0.78)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#DEDEDE")
    im = ax.imshow(
        np.ma.masked_invalid(relative.to_numpy()),
        aspect="auto",
        cmap=cmap,
        vmin=-8.5,
        vmax=8.5,
    )
    ax.set_yticks(range(5), genes)
    ax.set_xticks(range(12), [f"D{day}" for _ in range(3) for day in [0, 1, 2, 6]])
    for x, label in [(1.5, "P11"), (5.5, "P22"), (9.5, "P33")]:
        ax.text(x, -0.9, label, ha="center", fontsize=12, fontweight="bold")
    for x in [3.5, 7.5]:
        ax.axvline(x, color="white", lw=3)
    for y in range(5):
        for x in range(12):
            v = relative.iloc[y, x]
            if np.isfinite(v):
                ax.text(
                    x,
                    y,
                    f"{v:+.1f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white" if abs(v) > 5 else "#222222",
                )
    ax.set_xlabel("Differentiation day")
    cax = fig.add_axes([0.26, 0.11, 0.55, 0.035])
    fig.colorbar(
        im,
        cax=cax,
        orientation="horizontal",
        label="Mean log₂(normalized count + 1), relative to P11 D0",
    )
    fig.text(
        0.15,
        0.94,
        "C2C12 passage and differentiation context",
        fontsize=15,
        fontweight="bold",
    )
    for ext in ["png", "pdf"]:
        fig.savefig(
            a.out / ("Figure_C2C12_candidate_context." + ext),
            dpi=220,
            facecolor="white",
        )
    plt.close(fig)
    print("Saved private figure and exact plotted values.")


if __name__ == "__main__":
    main()
