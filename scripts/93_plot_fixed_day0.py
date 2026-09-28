"""Render candidate screening, source contrasts, and local RNA target panels."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/fixed_day0_D206"
OUTPUT = ROOT / "figures/fixed_day0_D206"
BLUE, RED, GRAY = "#326b9b", "#bf5145", "#b8bdc3"


def read(name):
    return pd.read_csv(SOURCE / name, sep="\t")


def save(figure, name):
    OUTPUT.mkdir(exist_ok=True)
    figure.savefig(OUTPUT / f"{name}.png", dpi=400, facecolor="white")
    figure.savefig(OUTPUT / f"{name}.pdf", facecolor="white")
    figure.savefig(OUTPUT / f"{name}.svg", facecolor="white")
    plt.close(figure)


def panel(ax, letter, title):
    ax.set_title(title, loc="left", fontsize=7, pad=9)
    ax.text(
        -0.06 if ax.get_position().width > 0.7 else -0.14,
        1.05,
        letter,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=9,
    )


def screen(peaks, hits, labels):
    fig = plt.figure(figsize=(7.2, 4.9))
    grid = fig.add_gridspec(
        2, 2, height_ratios=[0.50, 2.8], width_ratios=[1, 1.4], hspace=0.48, wspace=0.42
    )
    ax = fig.add_subplot(grid[0, :])
    ax.axis("off")
    steps = [
        "221 myogenesis\ngenes",
        "25 differentiation /\nfusion genes",
        "565 shared\nATAC peaks",
        "Day 0\n527 matched pairs",
        "6 regions\nq < 0.05",
    ]
    for i, text in enumerate(steps):
        ax.text(
            (i + 0.5) / 5,
            0.45,
            text,
            ha="center",
            va="center",
            fontsize=6.5,
            bbox=dict(
                boxstyle="round,pad=.45",
                facecolor="#f0f4f7",
                edgecolor="#b8c7d4",
                linewidth=0.5,
            ),
            transform=ax.transAxes,
        )
        if i < 4:
            ax.annotate(
                "",
                xy=((i + 1.07) / 5, 0.45),
                xytext=((i + 0.93) / 5, 0.45),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#788692", linewidth=0.6),
            )
    panel(ax, "a", "External candidate definition and COQ8A comparison")
    ax = fig.add_subplot(grid[1, 0])
    values = (
        np.column_stack(
            [
                peaks.source1_low / peaks.source1_pairs,
                peaks.source1_high / peaks.source1_pairs,
                peaks.source2_low / peaks.source2_pairs,
                peaks.source2_high / peaks.source2_pairs,
            ]
        )
        * 100
    )
    order = np.argsort(values.mean(1))[::-1]
    im = ax.imshow(
        values[order],
        aspect="auto",
        cmap="Blues",
        vmin=0,
        vmax=100,
        interpolation="nearest",
    )
    ax.set_xticks(
        range(4),
        ["Source 1\nLow", "Source 1\nHigh", "Source 2\nLow", "Source 2\nHigh"],
        fontsize=6,
    )
    ax.set_yticks([0, 564], [1, 565])
    ax.set_ylabel("Peaks ordered by mean accessibility")
    ax.axvline(1.5, color="white", linewidth=1.5)
    panel(ax, "b", "All 565 candidate peaks")
    colorbar = fig.colorbar(
        im, ax=ax, orientation="horizontal", pad=0.22, fraction=0.045, aspect=25
    )
    colorbar.set_label("Accessible nuclei (%)", fontsize=6.5)
    ax = fig.add_subplot(grid[1, 1])
    x = np.log2((peaks.high + 0.5) / (peaks.low + 0.5))
    y = -np.log10(peaks.p)
    ax.scatter(x, y, s=7, color=GRAY, alpha=0.65, linewidths=0)
    for name, color, mask in [
        ("Opening", RED, peaks.q565.lt(0.05) & peaks.FC.gt(1)),
        ("Closing", BLUE, peaks.q565.lt(0.05) & peaks.FC.lt(1)),
    ]:
        ax.scatter(x[mask], y[mask], s=15, color=color, label=name, linewidths=0)
    offsets = {
        "MYOD1": (3, 9),
        "CSRP3": (8, -4),
        "CACNA1H": (8, 9),
        "BOC": (0, 10),
        "NOTCH1": (-33, -10),
        "CAV3": (-14, -23),
    }
    for row in hits[hits.q25.lt(0.05)].itertuples():
        j = peaks.index[peaks.peak.eq(row.top_peak)][0]
        ax.annotate(
            row.gene,
            (x[j], y[j]),
            xytext=offsets[row.gene],
            textcoords="offset points",
            fontsize=6.5,
            arrowprops=dict(arrowstyle="-", color="#666666", linewidth=0.5),
        )
    ax.axvline(0, color="#777777", linewidth=0.6)
    ax.set_xlabel("log₂ accessibility ratio")
    ax.set_ylabel("−log₁₀ paired peak p")
    ax.set_ylim(-0.08, 4.7)
    ax.set_xlim(-2.5, 3.3)
    ax.legend(
        title="Peak q < 0.05",
        loc="upper right",
        bbox_to_anchor=(1, 0.57),
        frameon=False,
        fontsize=6,
        title_fontsize=6,
    )
    panel(ax, "c", "Complete peak screen")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.15)
    save(fig, "Figure_1_screen")


def locus_panels(peaks, hits, labels, primary=True):
    subset = hits[hits.q25.lt(0.05) if primary else hits.q25.ge(0.05)].sort_values(
        "q25"
    )
    order = subset.top_peak.drop_duplicates().tolist()
    rows, cols = (2, 3) if primary else (1, 3)
    fig, axes = plt.subplots(rows, cols, figsize=(7.2, 4.4 if primary else 3.2))
    counts = read("selected_ATAC_by_source.tsv")
    for ax, peak in zip(np.asarray(axes).flat, order):
        data = counts[counts.peak.eq(peak)].sort_values("gsm")
        value = peaks.set_index("peak").loc[peak]
        region = subset[subset.top_peak.eq(peak)]
        for i, row in enumerate(data.itertuples()):
            y = np.array([row.low, row.high]) / row.pairs * 100
            ax.plot([0, 1], y, "o-", color=[BLUE, RED][i], linewidth=1, markersize=3)
        ax.set_xticks([0, 1], ["Low", "High"])
        ax.set_xlim(-0.18, 1.18)
        ax.set_ylim(0, max(5, ax.get_ylim()[1] * 1.17))
        ax.set_ylabel("Accessible nuclei (%)", fontsize=6.5)
        ax.set_title(
            f"{labels[peak]}\nFC {value.FC:.2f}   q {region.q25.min():.3g}",
            fontsize=7,
            pad=6,
        )
        ax.text(0.5, -0.28, peak, transform=ax.transAxes, ha="center", fontsize=5.5)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=color,
            linewidth=1,
            markersize=3,
            label=f"Source {i+1} ({n} pairs)",
        )
        for i, (color, n) in enumerate([(BLUE, 226), (RED, 301)])
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.52, 0.005),
        fontsize=6.5,
    )
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        top=0.91 if primary else 0.83,
        bottom=0.30 if not primary else 0.15,
        hspace=0.9,
        wspace=0.58,
    )
    save(fig, "Figure_2_selected_loci" if primary else "Figure_S2_secondary_loci")


def rna_targets(hits, labels):
    evidence = read("target_evidence.tsv")
    shown = []
    order = hits.sort_values("q25").top_peak.drop_duplicates()
    for peak in order:
        data = evidence[evidence.peak.eq(peak)]
        extra = (
            data[data.r.notna()]
            .assign(abs_r=lambda x: x.r.abs())
            .sort_values("abs_r", ascending=False)
            .head(2)
            .gene
        )
        wanted = set(labels[peak].split("/")) | set(extra)
        shown.append(data[data.gene.isin(wanted)].sort_values("gene"))
    data = pd.concat(shown, ignore_index=True)
    matched = read("links_combined.tsv").query("population == 'matched'")
    for model, prefix in [
        ("depth_COQ", "matched"),
        ("depth_COQ_state", "matched_state"),
    ]:
        block = matched[matched.model.eq(model)][
            ["peak", "gene", "r", "q_all_local_links"]
        ].rename(columns={"r": prefix + "_r", "q_all_local_links": prefix + "_q"})
        data = data.merge(block, on=["peak", "gene"], how="left", validate="one_to_one")
    columns = [
        "peak",
        "gene",
        "r",
        "q_all_local_links",
        "state_r",
        "state_q",
        "matched_r",
        "matched_q",
        "matched_state_r",
        "matched_state_q",
        "FC",
        "p_RNA",
        "q_local_RNAs",
    ]
    data[columns].to_csv(
        SOURCE / "figure3_display_rows.tsv", sep="\t", index=False, na_rep="NA"
    )
    fig = plt.figure(figsize=(7.2, 6.1))
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.08], wspace=0.24)
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    rowlabels = [f"{labels[row.peak]} → {row.gene}" for row in data.itertuples()]
    fig.text(0.035, 0.95, "Neighborhood → candidate RNA", fontsize=6.5)
    breaks = np.flatnonzero(data.peak.to_numpy()[1:] != data.peak.to_numpy()[:-1]) + 0.5
    for k, (rs, qs, title) in enumerate(
        [
            (["r", "state_r"], ["q_all_local_links", "state_q"], "All eligible nuclei"),
            (
                ["matched_r", "matched_state_r"],
                ["matched_q", "matched_state_q"],
                "Matched nuclei",
            ),
        ]
    ):
        ax = axes[k]
        values = data[rs].to_numpy(float)
        im = ax.imshow(
            np.ma.masked_invalid(values),
            cmap="RdBu_r",
            vmin=-0.18,
            vmax=0.18,
            aspect="auto",
            interpolation="nearest",
        )
        ax.set_facecolor("#eeeeee")
        for i in range(len(data)):
            for j in range(2):
                value = values[i, j]
                star = "*" if data.iloc[i][qs[j]] < 0.05 else ""
                text = f"{value:.3f}{star}" if np.isfinite(value) else "NA"
                ax.text(
                    j,
                    i,
                    text,
                    ha="center",
                    va="center",
                    fontsize=5.8,
                    color="white" if abs(value) > 0.11 else "black",
                )
        for boundary in breaks:
            ax.axhline(boundary, color="white", linewidth=1)
        ax.set_xticks(
            [0, 1],
            ["Depth + COQ8A", "+ cell state"],
            rotation=30,
            ha="right",
            fontsize=6,
        )
        ax.set_yticks(range(len(data)), rowlabels if k == 0 else [], fontsize=6.2)
        panel(ax, "ab"[k], title)
    ax = axes[2]
    for i, row in enumerate(data.itertuples()):
        if np.isfinite(row.FC) and row.FC > 0:
            color = RED if row.FC > 1 else BLUE
            ax.scatter(
                np.log2(row.FC),
                i,
                s=17,
                facecolor=color if row.q_local_RNAs < 0.05 else "white",
                edgecolor=color,
                linewidth=0.8,
            )
            ax.annotate(
                f"{row.FC:.2f}",
                (np.log2(row.FC), i),
                xytext=(5, 0),
                textcoords="offset points",
                va="center",
                fontsize=6,
            )
    for boundary in breaks:
        ax.axhline(boundary, color="#eeeeee", linewidth=0.7, zorder=0)
    ax.axvline(0, color=GRAY, linewidth=0.6, zorder=0)
    ax.set_ylim(len(data) - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlim(-2.3, 2.5)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xlabel("RNA log₂ FC (high / low)", fontsize=6.5)
    panel(ax, "c", "Same-pair RNA")
    cb = fig.colorbar(
        im,
        cax=fig.add_axes([0.34, 0.06, 0.24, 0.018]),
        orientation="horizontal",
        ticks=[-0.15, 0, 0.15],
    )
    cb.set_label("Partial correlation r", fontsize=6.5)
    fig.text(0.33, 0.13, "* Link q < 0.05", fontsize=6.5)
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markeredgecolor="#555555",
            markerfacecolor=face,
            markersize=3.8,
            label=text,
        )
        for face, text in [("#555555", "RNA q < 0.05"), ("white", "RNA q ≥ 0.05")]
    ]
    fig.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(0.985, 0.055),
        frameon=False,
        fontsize=6.5,
    )
    fig.subplots_adjust(left=0.31, right=0.98, bottom=0.23, top=0.93)
    save(fig, "Figure_3_RNA_targets")


def all_regions(regions):
    values = regions.pivot(index="gene", columns="mode", values="q25")[
        ["opening", "closing"]
    ]
    values = values.loc[values.min(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(3.5, 5.8))
    im = ax.imshow(-np.log10(values), cmap="Blues", vmin=0, vmax=3, aspect="auto")
    for i in range(len(values)):
        for j in range(2):
            value = values.iloc[i, j]
            ax.text(
                j,
                i,
                f"{value:.3f}",
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if value < 0.03 else "black",
            )
    ax.set_yticks(range(len(values)), values.index)
    ax.set_xticks([0, 1], ["Opening", "Closing"])
    ax.set_title("All 25 candidate neighborhoods", fontsize=7)
    fig.colorbar(
        im, ax=ax, fraction=0.04, pad=0.08, shrink=0.5, label="−log₁₀ regional q"
    )
    fig.subplots_adjust(left=0.24, right=0.82, bottom=0.08, top=0.95)
    save(fig, "Figure_S1_all_regions")


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.linewidth": 0.6,
            "lines.linewidth": 0.8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    peaks, regions, hits = (
        read("all_peaks.tsv"),
        read("all_regions.tsv"),
        read("selected_regions.tsv"),
    )
    labels = (
        hits.groupby("top_peak", sort=False).gene.apply(lambda x: "/".join(x)).to_dict()
    )
    screen(peaks, hits, labels)
    locus_panels(peaks, hits, labels)
    locus_panels(peaks, hits, labels, primary=False)
    rna_targets(hits, labels)
    all_regions(regions)


if __name__ == "__main__":
    main()
