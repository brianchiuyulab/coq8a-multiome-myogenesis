"""Render manuscript figures from frozen analysis tables and fragment profiles."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

LOW, HIGH = "#9BA8B4", "#145A86"
PHASES = ["Early", "Middle", "Late", "Closing"]
PC = {"Early": "#268B9B", "Middle": "#A675B5", "Late": "#D2B35A", "Closing": "#C46D49"}
GENES = ["CSRP3", "CAV3", "MYOD1", "CACNA1H"]
LIBS = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]
LL = ["Line 1 / undiff.", "Line 1 / D7", "Line 2 / undiff.", "Line 2 / D7"]


def read(root, path):
    return pd.read_csv(root / path, sep="\t")


def main_only(d):
    return d[(d.gate == "TSS_ge_3") & (d.contrast == "3plus_vs_1")].copy()


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=220, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def letter(ax, text):
    ax.text(
        -0.10,
        1.04,
        text,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=14,
        va="bottom",
    )


def phase_blocks(ax, frame, label=False):
    offset = 0
    for phase in PHASES:
        n = int(frame.phase.eq(phase).sum())
        if not n:
            continue
        ax.axhline(offset - 0.5, color="white", lw=1.4)
        ax.add_patch(
            Rectangle(
                (-0.035, offset - 0.5),
                0.018,
                n,
                transform=ax.get_yaxis_transform(),
                color=PC[phase],
                clip_on=False,
                lw=0,
            )
        )
        if label:
            ax.text(
                -0.065,
                offset + (n - 1) / 2,
                f"{phase} ({n})",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=9,
                color=PC[phase],
            )
        offset += n


def external_figure(root, out):
    d = read(root, "results/figure_source/external_410_regions.tsv")
    fig = plt.figure(figsize=(15.4, 8.3))
    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[0.55, 3],
        width_ratios=[1.2, 1, 1],
        left=0.12,
        right=0.97,
        bottom=0.11,
        top=0.94,
        wspace=0.58,
        hspace=0.36,
    )
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    letter(ax, "A")
    labels = [
        "Myogenesis\n221 genes",
        "Common regions\n5,097 peaks",
        "External dynamics\n410 peaks",
        "Functional subset\n54 peaks",
    ]
    for i, label in enumerate(labels):
        x = 0.02 + i * 0.25
        ax.text(
            x,
            0.52,
            label,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=11,
            bbox=dict(boxstyle="round,pad=.5", fc="#F1F5F8", ec="#D2DBE2"),
        )
        if i < 3:
            ax.annotate(
                "",
                xy=(x + 0.235, 0.52),
                xytext=(x + 0.187, 0.52),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#71818D"),
            )
    ax = fig.add_subplot(gs[1, 0])
    letter(ax, "B")
    v = d[["p0", "p24", "p48", "p72"]].to_numpy()
    scaled = (v - v.min(1)[:, None]) / np.maximum(np.ptp(v, axis=1)[:, None], 1e-12)
    im = ax.imshow(
        scaled, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="nearest"
    )
    ax.set_xticks(range(4), ["0", "24", "48", "72"])
    ax.set_xlabel("Differentiation time (h)")
    ax.set_yticks([])
    ax.set_title("External differentiation dynamics", pad=13)
    phase_blocks(ax, d, label=True)
    cb = fig.colorbar(
        im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.13, shrink=0.9
    )
    cb.set_label("Accessibility (row-scaled)")
    sub = gs[1, 1].subgridspec(4, 1, hspace=0.95)
    for i, phase in enumerate(PHASES):
        ax = fig.add_subplot(sub[i])
        s = d[d.phase == phase]
        mean = s[["p0", "p24", "p48", "p72"]].mean().to_numpy() * 100
        ax.plot([0, 24, 48, 72], mean, "o-", color=PC[phase], ms=4, lw=1.8)
        ax.set_xlim(-3, 75)
        ax.set_xticks([0, 24, 48, 72])
        ax.set_ylim(bottom=0)
        ax.set_title(
            f"{phase} ({len(s)} regions)",
            loc="left",
            fontsize=10,
            color=PC[phase],
            pad=3,
        )
        ax.set_ylabel("Open (%)", fontsize=9)
        if i == 0:
            letter(ax, "C")
        if i == 3:
            ax.set_xlabel("Differentiation time (h)")
    effects = main_only(read(root, "results/temporal/temporal_effect_grid.tsv"))
    bylib = main_only(read(root, "results/temporal/temporal_effects_by_library.tsv.gz"))
    ids = [
        "P2_O50_Early_by24h",
        "P2_O50_Middle_24to48h",
        "P2_O50_Late_after48h",
        "P2_O50_Closing",
    ]
    target = gs[1, 2].subgridspec(4, 1, hspace=0.95)
    for i, (phase, sid) in enumerate(zip(PHASES, ids)):
        ax = fig.add_subplot(target[i])
        a = bylib[bylib.set_id == sid].set_index("gsm").loc[LIBS]
        r = effects[
            (effects.set_id == sid) & (effects.stage_analysis == "pooled")
        ].iloc[0]
        for j, row in enumerate(a.itertuples()):
            ax.plot(
                [row.low_open_pct, row.high_open_pct], [j, j], color="#BCC5CC", lw=1
            )
            ax.scatter(row.low_open_pct, j, color=LOW, s=20, zorder=3)
            ax.scatter(row.high_open_pct, j, color=HIGH, s=20, zorder=3)
        ax.set_yticks(range(4), ["L1 U", "L1 D7", "L2 U", "L2 D7"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(left=0)
        ax.set_title(
            f"FC {r.fold_open:.3f} · p {r.p_pair:.3g}\nq {r.q_within_family:.3g}",
            loc="left",
            fontsize=9,
            pad=3,
        )
        if i == 0:
            letter(ax, "D")
            ax.plot([], [], "o", color=LOW, label="Low", ms=4)
            ax.plot([], [], "o", color=HIGH, label="High", ms=4)
            fig.legend(
                *ax.get_legend_handles_labels(),
                loc="upper right",
                bbox_to_anchor=(0.965, 0.815),
                frameon=False,
                fontsize=8,
                ncol=2,
            )
        if i == 3:
            ax.set_xlabel("Mean accessible peaks (%)")
    fig.text(0.42, 0.77, "GSE109828: differentiation", fontsize=10, color="#52606A")
    fig.text(0.73, 0.77, "GSE208248: COQ8A high / low", fontsize=10, color="#52606A")
    save(fig, out / "main/Figure_1_External_definition")


def pooled_profiles(root):
    p = read(root, "results/figure_source/fragment_profiles.tsv.gz")
    return p.groupby(
        ["window_id", "kind", "group", "bin_start", "bin_end", "relative_bp"],
        as_index=False,
    ).agg(insertions=("insertions", "sum"), n_nuclei=("n_nuclei", "sum"))


def functional_figure(root, out):
    order = read(root, "results/figure_source/functional_54_order.tsv")
    profiles = pooled_profiles(root)
    profiles = profiles[profiles.kind == "peak"]
    matrices = []
    for group in ["low", "high"]:
        p = profiles[profiles.group == group].copy()
        p["value"] = p.insertions / p.n_nuclei * 100
        mat = (
            p.pivot(index="window_id", columns="relative_bp", values="value")
            .loc[order.region_id]
            .to_numpy()
        )
        matrices.append(gaussian_filter1d(mat, 1, axis=1, mode="nearest"))
    limit = np.quantile(np.concatenate([x.ravel() for x in matrices]), 0.995)
    effects = main_only(
        read(root, "results/tables/candidate_peak_effects_by_library.tsv.gz")
    )
    delta = effects.pivot(index="peak", columns="gsm", values="delta_pp").loc[
        order.peak, LIBS
    ]
    fig = plt.figure(figsize=(12.5, 11.5))
    gs = fig.add_gridspec(
        2,
        4,
        width_ratios=[1, 1, 1.1, 0.4],
        height_ratios=[1, 0.035],
        left=0.25,
        right=0.97,
        bottom=0.11,
        top=0.91,
        wspace=0.18,
        hspace=0.17,
    )
    for i in range(2):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(
            matrices[i],
            aspect="auto",
            cmap="Blues",
            vmin=0,
            vmax=limit,
            interpolation="nearest",
        )
        ax.set_xticks([0, 19.5, 39], ["−2", "0", "+2"])
        ax.set_xlabel("From peak centre (kb)")
        ax.set_title("COQ8A low" if i == 0 else "COQ8A high", pad=13)
        if i == 0:
            labels = [
                f"{r.nearby_candidate_genes}  {r.region_id}" for r in order.itertuples()
            ]
            ax.set_yticks(range(54), labels, fontsize=7)
            ax.text(
                -0.10,
                1.015,
                "A",
                transform=ax.transAxes,
                fontweight="bold",
                fontsize=14,
            )
            phase_blocks(ax, order)
        else:
            ax.set_yticks([])
            phase_blocks(ax, order)
    cax = fig.add_subplot(gs[1, :2])
    fig.colorbar(
        im,
        cax=cax,
        orientation="horizontal",
        extend="max",
        label="ATAC insertions per 100 nuclei / 100 bp",
    )
    ax = fig.add_subplot(gs[0, 2])
    letter(ax, "B")
    lim = 20
    im = ax.imshow(
        delta.to_numpy(),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-lim,
        vmax=lim,
        interpolation="nearest",
    )
    ax.set_title("Accessibility: high − low", pad=13)
    ax.set_yticks([])
    ax.set_xticks(
        range(4),
        ["L1 U", "L1 D7", "L2 U", "L2 D7"],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    phase_blocks(ax, order)
    cax = fig.add_subplot(gs[1, 2])
    fig.colorbar(
        im,
        cax=cax,
        orientation="horizontal",
        extend="both",
        label="High − low (percentage points)",
    )
    ax = fig.add_subplot(gs[0, 3])
    letter(ax, "C")
    vals = np.log2(order.fold_open.to_numpy())
    vals = np.ma.masked_invalid(np.clip(vals, -3, 3))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#D4D4D4")
    im = ax.imshow(
        vals[:, None],
        aspect="auto",
        cmap=cmap,
        vmin=-3,
        vmax=3,
        interpolation="nearest",
    )
    ax.set_title("High / low", pad=13)
    ax.set_xticks([0], ["FC"])
    ax.set_yticks([])
    phase_blocks(ax, order)
    cax = fig.add_subplot(gs[1, 3])
    fig.colorbar(
        im,
        cax=cax,
        orientation="horizontal",
        extend="both",
        label="log₂ FC",
        ticks=[-3, 0, 3],
    )
    fig.text(0.25, 0.975, "Functional myogenic regions", fontsize=16, fontweight="bold")
    fig.text(
        0.25,
        0.946,
        "54 regions · 201 matched pairs · 2 source lines",
        fontsize=10,
        color="#52606A",
    )
    save(fig, out / "main/Figure_2_Functional_accessibility")


def locus_figure(root, out):
    prof = pooled_profiles(root)
    windows = read(root, "results/figure_source/fragment_windows.tsv").set_index(
        "window_id"
    )
    tx = read(root, "results/figure_source/focal_transcripts.tsv").set_index("gene")
    exons = read(root, "results/figure_source/focal_exons.tsv")
    stats = read(root, "results/figure_source/focal_statistics.tsv").set_index("gene")
    atac = main_only(
        read(root, "results/tables/candidate_peak_effects_by_library.tsv.gz")
    )
    rna = main_only(read(root, "results/tables/gene_effects_by_library.tsv.gz"))
    rna = rna[rna.modality == "RNA"]
    fig = plt.figure(figsize=(13, 11.6))
    gs = fig.add_gridspec(
        4,
        3,
        width_ratios=[1.65, 1, 1],
        left=0.085,
        right=0.97,
        bottom=0.08,
        top=0.885,
        hspace=0.95,
        wspace=0.53,
    )
    for i, gene in enumerate(GENES):
        w = windows.loc[gene]
        s = stats.loc[gene]
        ax = fig.add_subplot(gs[i, 0])
        for group, color in [("low", LOW), ("high", HIGH)]:
            p = prof[(prof.window_id == gene) & (prof.group == group)].sort_values(
                "bin_start"
            )
            v = p.insertions / p.n_nuclei * 100 * 100 / (p.bin_end - p.bin_start)
            y = gaussian_filter1d(v.to_numpy(), 1, mode="nearest")
            ax.plot(
                (p.bin_start + p.bin_end) / 2 / 1000,
                y,
                color=color,
                lw=1.5,
                label=group.capitalize(),
            )
            ax.fill_between(
                (p.bin_start + p.bin_end) / 2 / 1000, 0, y, color=color, alpha=0.10
            )
        lo, hi = map(int, w.peak.split(":")[1].split("-"))
        ax.axvspan(lo / 1000, hi / 1000, color="#D9B36A", alpha=0.32, lw=0)
        ax.set_xlim(w.start / 1000, w.end / 1000)
        ax.set_ylim(bottom=0)
        ax.set_title(gene, loc="left", fontweight="bold", fontsize=12)
        ax.set_ylabel("ATAC insertions\nper 100 nuclei / 100 bp", fontsize=8)
        ax.set_xlabel(f"{w.chrom} position (kb)", fontsize=8, labelpad=28)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        tr = tx.loc[gene]
        transform = ax.get_xaxis_transform()
        start = max(w.start, tr.start) / 1000
        end = min(w.end, tr.end) / 1000
        ax.plot(
            [start, end],
            [-0.23, -0.23],
            transform=transform,
            color="#3B4146",
            lw=0.7,
            clip_on=False,
        )
        for r in exons[exons.transcript == tr.transcript].itertuples():
            left = max(w.start, r.start) / 1000
            right = min(w.end, r.end) / 1000
            if left < right:
                ax.add_patch(
                    Rectangle(
                        (left, -0.26),
                        right - left,
                        0.06,
                        transform=transform,
                        color="#3B4146",
                        clip_on=False,
                    )
                )
        direction = 1 if tr.strand == "+" else -1
        midpoint = (start + end) / 2
        ax.annotate(
            "",
            xy=(midpoint + direction * (end - start) * 0.12, -0.23),
            xytext=(midpoint - direction * (end - start) * 0.12, -0.23),
            xycoords=transform,
            arrowprops=dict(arrowstyle="->", color="#3B4146", lw=0.7),
            annotation_clip=False,
        )
        if i == 0:
            letter(ax, "A")
            ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax = fig.add_subplot(gs[i, 1])
        a = atac[atac.peak == w.peak].set_index("gsm").loc[LIBS]
        for j, r in enumerate(a.itertuples()):
            lo2, hi2 = r.low_open / r.n_pairs * 100, r.high_open / r.n_pairs * 100
            ax.plot([lo2, hi2], [j, j], color="#BCC5CC", lw=1)
            ax.scatter(lo2, j, color=LOW, s=28, zorder=3)
            ax.scatter(hi2, j, color=HIGH, s=28, zorder=3)
        ax.set_yticks(range(4), LL, fontsize=7)
        ax.invert_yaxis()
        xmax = max(
            (a.high_open / a.n_pairs * 100).max(),
            (a.low_open / a.n_pairs * 100).max(),
            1,
        )
        ax.set_xlim(-0.04 * xmax, 1.10 * xmax)
        ax.set_xlabel("Open nuclei (%)", fontsize=9)
        ax.set_title(
            f"FC {s.fold_open:.2f} · p {s.p_exact:.3g}\nq {s.q_functional_union:.3g}",
            fontsize=9,
            loc="left",
        )
        if i == 0:
            letter(ax, "B")
        ax = fig.add_subplot(gs[i, 2])
        a = rna[rna.gene == gene].set_index("gsm").loc[LIBS]
        for j, r in enumerate(a.itertuples()):
            ax.plot([r.low_mean, r.high_mean], [j, j], color="#BCC5CC", lw=1)
            ax.scatter(r.low_mean, j, color=LOW, s=28, zorder=3)
            ax.scatter(r.high_mean, j, color=HIGH, s=28, zorder=3)
        ax.set_yticks(range(4), LL, fontsize=7)
        ax.invert_yaxis()
        xmax = max(a.high_mean.max(), a.low_mean.max(), 0.001)
        ax.set_xlim(-0.04 * xmax, 1.10 * xmax)
        ax.set_xlabel("Mean RNA logₑ(1 + CP10k)", fontsize=9)
        ax.set_title(
            f"Δ {s.rna_difference:.3f} · p {s.rna_p_pair:.3g}\nq {s.rna_q_221:.3g}",
            fontsize=9,
            loc="left",
        )
        if i == 0:
            letter(ax, "C")
    fig.text(
        0.085,
        0.975,
        "Candidate loci and paired RNA measurements",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(0.085, 0.944, "ATAC fragment signal", fontsize=11, color="#52606A")
    fig.text(0.493, 0.944, "ATAC accessibility", fontsize=11, color="#52606A")
    fig.text(0.787, 0.944, "RNA expression", fontsize=11, color="#52606A")
    save(fig, out / "main/Figure_3_Candidate_loci")


def sensitivity_figure(root, out):
    temporal = read(root, "results/temporal/temporal_effect_grid.tsv")
    temporal = temporal.query(
        "scope=='programme' and promoter_kb==2 and reciprocal_overlap==0.5 and stage_analysis=='pooled'"
    )
    ids = [
        "P2_O50_Early_by24h",
        "P2_O50_Middle_24to48h",
        "P2_O50_Late_after48h",
        "P2_O50_Closing",
    ]
    functional = read(root, "results/functional/functional_effects.tsv")
    settings = [
        ("TSS_ge_2", "2plus_vs_1"),
        ("TSS_ge_3", "2plus_vs_1"),
        ("TSS_ge_2", "3plus_vs_1"),
        ("TSS_ge_3", "3plus_vs_1"),
    ]
    fig, axes = plt.subplots(
        1, 2, figsize=(12, 6.5), gridspec_kw=dict(width_ratios=[1, 1.5])
    )
    fig.subplots_adjust(left=0.09, right=0.97, top=0.86, bottom=0.22, wspace=0.56)
    for panel, (frame, rowids, labels, title) in enumerate(
        [
            (
                temporal,
                ids,
                ["Early (191)", "Middle (34)", "Late (4)", "Closing (181)"],
                "Temporal regions",
            ),
            (
                functional,
                sorted(functional.set_id.unique()),
                None,
                "Functional regions",
            ),
        ]
    ):
        values = []
        pv = []
        for sid in rowids:
            sub = frame[frame.set_id == sid].set_index(["gate", "contrast"])
            values.append([sub.loc[z, "fold_open"] for z in settings])
            pv.append([sub.loc[z, "p_pair"] for z in settings])
        values = np.array(values)
        pv = np.array(pv)
        ax = axes[panel]
        im = ax.imshow(
            np.log2(values), cmap="RdBu_r", vmin=-0.35, vmax=0.35, aspect="auto"
        )
        for y in range(len(rowids)):
            for x in range(4):
                ax.text(
                    x,
                    y,
                    f"{values[y,x]:.3f}\np={pv[y,x]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=8,
                )
        if labels is None:
            labels = [s.replace("__", " / ").replace("_", " ") for s in rowids]
        ax.set_yticks(range(len(rowids)), labels, fontsize=8)
        ax.set_xticks(
            range(4),
            ["TSS≥2\n≥2 vs 1", "TSS≥3\n≥2 vs 1", "TSS≥2\n≥3 vs 1", "TSS≥3\n≥3 vs 1"],
            fontsize=8,
        )
        ax.set_title(title)
        letter(ax, chr(65 + panel))
    cax = fig.add_axes([0.37, 0.08, 0.32, 0.025])
    fig.colorbar(im, cax=cax, orientation="horizontal", label="log₂ accessibility FC")
    save(fig, out / "supplement/Figure_S2_Module_sensitivity")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    a = p.parse_args()
    out = a.root / "figures"
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
        }
    )
    external_figure(a.root, out)
    functional_figure(a.root, out)
    locus_figure(a.root, out)
    sensitivity_figure(a.root, out)
    print("Rendered three main figures and module sensitivity supplement.")


if __name__ == "__main__":
    main()
