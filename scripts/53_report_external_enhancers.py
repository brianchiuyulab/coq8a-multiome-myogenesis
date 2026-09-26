"""Report the independent enhancer annotation and distinguish time from COQ8A."""

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(a):
    out = a.root / "results/external_enhancers"
    out.mkdir(parents=True, exist_ok=True)
    peak = "chr11:17649919-17650798"
    if a.work:
        for f in a.work.glob("*.tsv"):
            shutil.copy2(f, out / f.name)
        shutil.copy2(
            a.work / "reference_manifest.json", out / "reference_manifest.json"
        )
        temporal = pd.read_csv(
            a.root / "results/temporal/target_temporal_annotations.tsv.gz", sep="\t"
        )
        temporal = temporal[(temporal.peak == peak) & (temporal.promoter_kb == 2)]
        temporal.to_csv(out / "MYOD1_external_time_profile.tsv", sep="\t", index=False)
        states = pd.read_csv(a.hmm, sep="\t", header=None)
        local = states[
            (states[0] == "chr11") & (states[1] < 17678000) & (states[2] > 17669000)
        ].iloc[:, :4]
        local.columns = ["chrom", "start", "end", "state"]
        local.to_csv(out / "MYOD1_HSMM_state_segments.tsv", sep="\t", index=False)
        rows = []
        for stage, name in [
            ("Myoblast", "GSM1218849_MB135GMMD.peak.txt.gz"),
            ("Myotube", "GSM1218850_MB135DMMD.peak.txt.gz"),
        ]:
            chip = pd.read_csv(a.chip / name, sep="\t", header=None)
            chip = chip[
                (chip[0] == "chr11") & (chip[1] < 17678000) & (chip[2] > 17669000)
            ]
            rows.extend(
                dict(stage=stage, chrom=r[0], start=r[1], end=r[2])
                for _, r in chip.iterrows()
            )
        pd.DataFrame(rows).to_csv(
            out / "MYOD1_local_ChIP_intervals.tsv", sep="\t", index=False
        )
        c = pd.read_csv(a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
        m = pd.read_csv(
            a.root / "results/temporal/temporal_set_memberships.tsv.gz", sep="\t"
        )
        early = m[m.set_id.eq("P2_O50_Early_by24h")].merge(
            c[c.gene.eq("MYOD1")], on="peak"
        )
        effects = pd.read_csv(
            a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t"
        )
        early.merge(effects, on="peak").to_csv(
            out / "MYOD1_early_peaks.tsv", sep="\t", index=False
        )
    states = pd.read_csv(out / "MYOD1_HSMM_state_segments.tsv", sep="\t")
    chip = pd.read_csv(out / "MYOD1_local_ChIP_intervals.tsv", sep="\t")
    time = pd.read_csv(out / "MYOD1_external_time_profile.tsv", sep="\t").iloc[0]
    ann = pd.read_csv(out / "external_peak_annotations.tsv", sep="\t")
    target = ann[ann.peak.eq(peak)].iloc[0]
    assert target.strong_enhancer_fraction == 1 and target.MYOD_chip_calls > 0
    effects = pd.read_csv(a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t")
    effect = effects[effects.peak.eq(peak)].iloc[0]
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42}
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(13, 3.8), gridspec_kw={"width_ratios": [1.4, 1, 1]}
    )
    ax = axes[0]
    ax.axvspan(target.start19 / 1e6, target.end19 / 1e6, color="#EAE4D4", zorder=0)
    for r in states.itertuples():
        if r.state in ["4_Strong_Enhancer", "5_Strong_Enhancer"]:
            ax.broken_barh(
                [(r.start / 1e6, (r.end - r.start) / 1e6)],
                (2.75, 0.45),
                facecolors="#DAA34A",
            )
    for stage, y, color in [
        ("Myoblast", 1.75, "#3E7586"),
        ("Myotube", 0.75, "#65A3A4"),
    ]:
        for r in chip[chip.stage.eq(stage)].itertuples():
            ax.broken_barh(
                [(r.start / 1e6, (r.end - r.start) / 1e6)], (y, 0.45), facecolors=color
            )
    ax.broken_barh(
        [(target.start19 / 1e6, (target.end19 - target.start19) / 1e6)],
        (-0.25, 0.45),
        facecolors="#333333",
    )
    ax.set(
        xlim=(17.669, 17.678),
        ylim=(-0.65, 3.6),
        yticks=[0, 1, 2, 3],
        yticklabels=[
            "Target ATAC peak",
            "MYOD ChIP: myotube",
            "MYOD ChIP: myoblast",
            "HSMM strong enhancer",
        ],
        xlabel="chr11 position (Mb; hg19)",
        title="Independent regulatory annotation",
    )
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax = axes[1]
    hours = np.array([0, 24, 48, 72])
    values = np.array([time[f"p{h}"] for h in hours]) * 100
    ax.plot(hours, values, "o-", color="#3E7586", label="Experiment 1")
    ax.plot(
        [0, 72],
        [time.rep_p0 * 100, time.rep_p72 * 100],
        "s--",
        color="#969696",
        label="Experiment 2",
    )
    ax.set(
        xticks=hours,
        ylim=(0, 40),
        xlabel="Differentiation time (h)",
        ylabel="Depth-standardized detection (%)",
        title="External time course: GSE109828",
    )
    ax.legend(frameon=False, fontsize=8)
    ax = axes[2]
    values = np.array([effect.low_open, effect.high_open]) / effect.n_pairs * 100
    ax.bar([0, 1], values, color=["#8AA8BC", "#C47439"], width=0.6)
    ax.set(
        xticks=[0, 1],
        xticklabels=["COQ8A low", "COQ8A high"],
        ylim=(0, 33),
        ylabel="Nuclei with peak detected (%)",
        title="Matched comparison: GSE208248",
    )
    ax.text(
        0.5,
        27,
        f"FC = {effect.fold_open:.2f}\np = {effect.p_pair_binomial:.3f}",
        ha="center",
    )
    for ax, letter in zip(axes, "ABC"):
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.17, 1.09, letter, transform=ax.transAxes, weight="bold", fontsize=14)
    fig.subplots_adjust(left=0.16, right=0.985, bottom=0.19, top=0.83, wspace=0.55)
    folder = a.root / "figures/external_enhancers"
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        fig.savefig(folder / ("MYOD1_external_annotation_and_time." + ext), dpi=220)
    print(
        "Focal peak fully covered by independent strong-enhancer annotation; figure rendered."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path)
    p.add_argument("--hmm", type=Path)
    p.add_argument("--chip", type=Path)
    main(p.parse_args())
