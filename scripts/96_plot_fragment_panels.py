"""TSS/peak-centered ATAC profiles and locus tracks from fixed matched nuclei."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
src = root / "results/fixed_day0_D206"
out = root / "figures/fixed_day0_D206"
d = pd.read_csv(src / "fragment_profiles.tsv.gz", sep="\t")
windows = pd.read_csv(src / "profile_windows.tsv", sep="\t")
scope = set(
    pd.read_csv(root / "results/differentiation_fusion25/genes.tsv", sep="\t").gene
)
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)
colors = {"low": "#777777", "high": "#397db8"}


def save(fig, name):
    fig.savefig(out / (name + ".png"), dpi=400, facecolor="white")
    fig.savefig(out / (name + ".pdf"))
    fig.savefig(out / (name + ".svg"))
    plt.close(fig)


def matrix(kind, subset=None):
    z = d[d.kind.eq(kind)]
    if subset is not None:
        z = z[z.id.isin(subset)]
    # Normalize each source per nucleus first, then average the two sources equally.
    z = z.groupby(["id", "group", "bin"]).per100_nuclei_per100bp.mean().reset_index()
    mats = {
        g: z[z.group.eq(g)].pivot(
            index="id", columns="bin", values="per100_nuclei_per100bp"
        )
        for g in ["low", "high"]
    }
    order = (
        ((mats["low"] + mats["high"]) / 2).mean(1).sort_values(ascending=False).index
    )
    return {g: m.loc[order].to_numpy() for g, m in mats.items()}, order


fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(
    2,
    5,
    width_ratios=[1, 1, 0.28, 1, 1],
    height_ratios=[1, 3],
    hspace=0.08,
    wspace=0.18,
)
for kind, subset, cols, half, title in [
    ("TSS", scope, [0, 1], 5000, "25 differentiation / fusion genes"),
    ("peak", None, [3, 4], 2000, "All 565 candidate peaks"),
]:
    mats, order = matrix(kind, subset)
    limit = np.quantile(np.concatenate(list(mats.values())), 0.99)
    ymax = max(v.mean(0).max() for v in mats.values()) * 1.12
    for col, g in zip(cols, ["low", "high"]):
        top = fig.add_subplot(gs[0, col])
        bot = fig.add_subplot(gs[1, col])
        v = mats[g]
        xs = -half + (np.arange(v.shape[1]) + 0.5) * (2 * half / v.shape[1])
        top.plot(xs, v.mean(0), c=colors[g], lw=1)
        top.set_ylim(0, ymax)
        top.set_xlim(-half, half)
        top.set_xticks([])
        top.set_title(
            f"COQ8A {g}\n" + ("0 UMI" if g == "low" else "≥2 UMI"), fontsize=7
        )
        if g == "low":
            top.set_ylabel("Mean insertions", fontsize=6.5)
        else:
            top.tick_params(axis="y", labelleft=False)
        im = bot.imshow(
            v,
            aspect="auto",
            extent=[-half, half, len(v), 0],
            cmap="Blues",
            vmin=0,
            vmax=limit,
            interpolation="nearest",
        )
        bot.set_xticks(
            [-half, 0, half],
            [
                f"-{half//1000} kb",
                "TSS" if kind == "TSS" else "Peak center",
                f"+{half//1000} kb",
            ],
        )
        bot.set_xlabel("Genomic position")
        bot.get_xticklabels()[0].set_ha("left")
        bot.get_xticklabels()[-1].set_ha("right")
        if kind == "TSS" and g == "low":
            bot.set_yticks(np.arange(len(order)) + 0.5, order, fontsize=6)
        elif g == "low":
            bot.set_yticks([0, len(v)], [1, len(v)])
            bot.set_ylabel("Peaks, common row order")
        else:
            bot.set_yticks([])
        if g == "high":
            cb = fig.colorbar(
                im, cax=bot.inset_axes([0, -0.18, 1, 0.025]), orientation="horizontal"
            )
            cb.set_label("Insertions / 100 nuclei / 100 bp", fontsize=5.5)
    fig.text(0.29 if kind == "TSS" else 0.76, 0.97, title, ha="center", fontsize=7)
fig.subplots_adjust(top=0.86, bottom=0.16, left=0.12, right=0.98)
save(fig, "Figure_4_position_aligned_ATAC")

selected = {
    "MYOD1": "chr11:17653349-17654252",
    "CSRP3": "chr11:19218592-19219518",
    "CAV3": "chr3:8757835-8758737",
}
fig, axs = plt.subplots(3, 1, figsize=(7.2, 5.4), layout="constrained")
for ax, (gene, peak) in zip(axs, selected.items()):
    win = windows[(windows.kind == "locus") & (windows.id == gene)].iloc[0]
    z = (
        d[(d.kind == "locus") & (d.id == gene)]
        .groupby(["group", "bin"])
        .per100_nuclei_per100bp.mean()
        .reset_index()
    )
    top = 0
    for g in ["low", "high"]:
        zz = z[z.group == g].sort_values("bin")
        left = win.start + zz.bin.to_numpy() * 100
        x = (left + np.minimum(100, win.end - left) / 2) / 1e6
        y = gaussian_filter1d(zz.per100_nuclei_per100bp.to_numpy(), 1)
        ax.plot(x, y, c=colors[g], lw=0.8, label=f"COQ8A {g}")
        ax.fill_between(x, 0, y, color=colors[g], alpha=0.15)
        top = max(top, y.max())
    lo, hi = map(int, peak.split(":")[1].split("-"))
    ax.axvspan(lo / 1e6, hi / 1e6, color="#edb34c", alpha=0.55)
    ax.axvline(win.center / 1e6, c="#444444", ls=":", lw=0.7)
    ax.text(win.center / 1e6, top * 1.1, f"{gene} TSS", ha="center", fontsize=6.5)
    ax.set_ylim(0, top * 1.3)
    ax.set_xlim(win.start / 1e6, win.end / 1e6)
    ax.ticklabel_format(useOffset=False, axis="x")
    ax.set_xlabel(f"{win.chrom} position (Mb)")
    ax.set_ylabel("ATAC insertions\n/ 100 nuclei / 100 bp", fontsize=6.5)
    ax.set_title(gene, loc="left", fontsize=7)
    ax.legend(loc="upper center", frameon=False, fontsize=6)
save(fig, "Figure_5_local_ATAC_tracks")
print("Saved position-aligned heatmaps and locus tracks")
