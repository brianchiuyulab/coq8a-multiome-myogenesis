"""Validate, package and visualize the temporally unstratified exploration."""

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

FOCAL = {
    "MYOD1": "chr11:17649919-17650798",
    "CSRP3": "chr11:19201752-19202603",
    "CAV3": "chr3:8733438-8733968",
}


def validate_and_rna(a):
    genes = (a.work / "genes.txt").read_text().splitlines()
    peaks = (a.work / "peaks.txt").read_text().splitlines()
    pairs = pd.read_csv(a.root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "3plus_vs_1")]
    totals = {
        k: np.zeros(len(peaks), dtype=int)
        for k in ["high_open", "low_open", "high_only", "low_only"]
    }
    high = []
    low = []
    librows = []
    for gsm, ps in pairs.groupby("gsm"):
        meta = pd.read_csv(a.work / (gsm + "_meta.tsv"), sep="\t").set_index("barcode")
        hi = meta.index.get_indexer(ps.high_barcode)
        lo = meta.index.get_indexer(ps.low_barcode)
        assert min(hi.min(), lo.min()) >= 0
        atac = sparse.load_npz(a.work / (gsm + "_atac.npz"))
        h = atac[:, hi].toarray()
        l = atac[:, lo].toarray()
        for key, values in [
            ("high_open", h),
            ("low_open", l),
            ("high_only", h * (1 - l)),
            ("low_only", l * (1 - h)),
        ]:
            totals[key] += values.sum(1).astype(int)
        raw = sparse.load_npz(a.work / (gsm + "_rna.npz"))
        yh = np.log1p(
            raw[:, hi].toarray().astype(float)
            / meta.iloc[hi].total_rna_umi.to_numpy()[None, :]
            * 10000
        )
        yl = np.log1p(
            raw[:, lo].toarray().astype(float)
            / meta.iloc[lo].total_rna_umi.to_numpy()[None, :]
            * 10000
        )
        high.append(yh)
        low.append(yl)
        librows.append(
            pd.DataFrame({"gene": genes, "gsm": gsm, "difference": (yh - yl).mean(1)})
        )
    effects = (
        pd.read_csv(a.work / "ATAC_all_5097.tsv", sep="\t").set_index("peak").loc[peaks]
    )
    for key, values in totals.items():
        assert np.array_equal(values, effects[key])
    discord = totals["high_only"] + totals["low_only"]
    p = np.minimum(
        1,
        2
        * stats.binom.cdf(
            np.minimum(totals["high_only"], totals["low_only"]), discord, 0.5
        ),
    )
    assert np.allclose(p, effects.p_pair_binomial)
    h = np.hstack(high)
    l = np.hstack(low)
    delta = h - l
    with np.errstate(invalid="ignore", divide="ignore"):
        t = delta.mean(1) / (delta.std(1, ddof=1) / np.sqrt(delta.shape[1]))
        pv = 2 * stats.t.sf(np.abs(t), delta.shape[1] - 1)
    pv = np.where(np.isfinite(pv), pv, 1)
    per = pd.concat(librows)
    signs = per.groupby("gene").difference.apply(lambda x: int(sum(x > 0)))
    rna = pd.DataFrame(
        dict(
            gene=genes,
            n_pairs=delta.shape[1],
            high_mean=h.mean(1),
            low_mean=l.mean(1),
            difference=delta.mean(1),
            p_pair=pv,
            positive_libraries=[signs[g] for g in genes],
        )
    )
    rna["q_measured_targets"] = multipletests(pv, method="fdr_bh")[1]
    rna.to_csv(a.work / "RNA_COQ8A_high_low.tsv", sep="\t", index=False)
    old = pd.read_csv(a.root / "results/tables/gene_effects_pooled.tsv", sep="\t")
    old = old[
        (old.gate == "TSS_ge_3")
        & (old.contrast == "3plus_vs_1")
        & (old.modality == "RNA")
    ]
    matched = rna.merge(old, on="gene", suffixes=("_new", "_old"))
    assert np.allclose(matched.difference_new, matched.difference_old, atol=1e-7)
    valid = matched.p_pair_old.notna()
    assert np.allclose(
        matched.loc[valid, "p_pair_new"], matched.loc[valid, "p_pair_old"], atol=1e-7
    )
    validation = json.loads((a.work / "validation.json").read_text())
    validation.update(
        ATAC_peaks_reproduced=len(peaks),
        RNA_gene_effects_reproduced=len(matched),
        RNA_high_low_genes=len(rna),
    )
    (a.work / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")


def package(a, out):
    validate_and_rna(a)
    files = [
        "ATAC_all_5097.tsv",
        "count_followup_selection.tsv",
        "count_candidate_pairs.tsv",
        "count_models_by_library.tsv",
        "count_models_combined.tsv",
        "scent_style_bootstrap_followup.tsv",
        "bootstrap_targets.tsv",
        "correlations_combined.tsv.gz",
        "correlations_by_library.tsv.gz",
        "programme_genes.txt",
        "validation.json",
        "RNA_COQ8A_high_low.tsv",
    ]
    for name in files:
        shutil.copy2(a.work / name, out / name)
    pd.read_csv(a.work / "candidate_pairs.tsv", sep="\t").to_csv(
        out / "candidate_pairs.tsv.gz", sep="\t", index=False
    )
    sensitive = out / "detection_1pct"
    sensitive.mkdir(exist_ok=True)
    for name in [
        "count_candidate_pairs.tsv",
        "count_models_by_library.tsv",
        "count_models_combined.tsv",
        "scent_style_bootstrap_followup.tsv",
    ]:
        file = a.work / "detection_1pct" / name
        if file.exists():
            shutil.copy2(file, sensitive / name)
    # Independent annotation uses author peak calls, without conditioning on
    # the target accessibility or RNA result.
    mapping = pd.read_csv(
        a.root / "results/temporal/candidate_hg19_mapping.tsv", sep="\t"
    ).set_index("peak")
    external = []
    for stage, name in [
        ("myoblast", "GSM1218849_MB135GMMD.peak.txt.gz"),
        ("myotube_72h", "GSM1218850_MB135DMMD.peak.txt.gz"),
    ]:
        chip = pd.read_csv(a.external / name, sep="\t", header=None)
        for gene, peak in FOCAL.items():
            r = mapping.loc[peak]
            s = chip[
                (chip[0] == r.chrom19) & (chip[1] < r.end19) & (chip[2] > r.start19)
            ]
            external.append(
                dict(gene=gene, peak=peak, stage=stage, overlap_peak_calls=len(s))
            )
    pd.DataFrame(external).to_csv(
        out / "focal_MYOD1_ChIP_overlap.tsv", sep="\t", index=False
    )


def report(a, out):
    effects = pd.read_csv(out / "ATAC_all_5097.tsv", sep="\t")
    counts = pd.read_csv(out / "count_models_combined.tsv", sep="\t")
    cor = pd.read_csv(out / "correlations_combined.tsv.gz", sep="\t")
    rna = pd.read_csv(out / "RNA_COQ8A_high_low.tsv", sep="\t")
    main = counts[counts.model.eq("technical_state_coq")].rename(
        columns={"p": "p_peak_RNA", "positive_libraries": "positive_link_libraries"}
    )
    chain = main.merge(
        effects[
            ["peak", "fold_open", "p_pair_binomial", "positive_libraries", "delta_pp"]
        ],
        on="peak",
    ).rename(columns={"positive_libraries": "positive_ATAC_libraries"})
    chain = chain.merge(
        rna[["gene", "difference", "p_pair", "positive_libraries"]], on="gene"
    ).rename(columns={"positive_libraries": "positive_RNA_libraries"})
    chain.to_csv(out / "complete_followup_evidence.tsv", sep="\t", index=False)
    coherent = chain[
        (chain.fold_open > 1)
        & (chain.p_pair_binomial < 0.05)
        & (chain.RNA_ratio > 1)
        & (chain.p_peak_RNA < 0.05)
        & (chain.difference > 0)
        & (chain.p_pair < 0.05)
    ]
    coherent.sort_values("p_pair_binomial").to_csv(
        out / "nominal_positive_chain.tsv", sep="\t", index=False
    )
    focalrows = []
    for gene, peak in FOCAL.items():
        e = effects[effects.peak.eq(peak)].iloc[0]
        rr = rna[rna.gene.eq(gene)].iloc[0]
        det = 0.01 if gene == "CSRP3" else 0.05
        cc = cor[
            (cor.gene == gene)
            & (cor.peak == peak)
            & (cor.model == "technical_state_coq")
            & (cor.detection_fraction == det)
        ].iloc[0]
        focalrows.append(
            dict(
                gene=gene,
                peak=peak,
                ATAC_FC=e.fold_open,
                p_ATAC=e.p_pair_binomial,
                high_open=int(e.high_open),
                low_open=int(e.low_open),
                n_pairs=int(e.n_pairs),
                RNA_delta=rr.difference,
                p_RNA=rr.p_pair,
                r=cc.r,
                p_link=cc.p,
                detection_fraction=det,
                n_link_libraries=cc.n_libraries,
            )
        )
    focal = pd.DataFrame(focalrows)
    focal.to_csv(out / "focal_figure_source.tsv", sep="\t", index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 4.2), gridspec_kw={"width_ratios": [1.05, 1, 1.15]}
    )
    ax = axes[0]
    color = np.where(
        effects.p_pair_binomial < 0.05,
        np.where(effects.delta_pp > 0, "#BA622E", "#437DA2"),
        "#C3C8CA",
    )
    ax.scatter(
        effects.delta_pp,
        -np.log10(effects.p_pair_binomial.clip(lower=1e-300)),
        c=color,
        s=8,
        alpha=0.65,
        linewidths=0,
    )
    ax.axhline(-np.log10(0.05), color="#999999", ls=":", lw=1)
    offsets = {"MYOD1": (3, 8), "CSRP3": (3, -15), "CAV3": (3, -15)}
    for r in focal.itertuples():
        row = effects[effects.peak.eq(r.peak)].iloc[0]
        ax.scatter(
            row.delta_pp, -np.log10(row.p_pair_binomial), s=24, c="#202020", zorder=4
        )
        ax.annotate(
            r.gene,
            (row.delta_pp, -np.log10(row.p_pair_binomial)),
            xytext=offsets[r.gene],
            textcoords="offset points",
            fontsize=8,
        )
    ax.set(
        xlabel="Accessibility difference (percentage points)",
        ylabel="−log₁₀ nominal p",
        title="All 5,097 candidate peaks",
    )
    ax = axes[1]
    x = np.arange(3)
    width = 0.32
    ax.bar(
        x - width / 2,
        focal.low_open / focal.n_pairs * 100,
        width,
        color="#8AA8BC",
        label="COQ8A low",
    )
    ax.bar(
        x + width / 2,
        focal.high_open / focal.n_pairs * 100,
        width,
        color="#C47439",
        label="COQ8A high",
    )
    for i, r in enumerate(focal.itertuples()):
        height = max(r.high_open, r.low_open) / r.n_pairs * 100
        ax.text(
            i,
            height + 1,
            f"{r.ATAC_FC:.2f}×\np={r.p_ATAC:.3f}",
            ha="center",
            fontsize=8,
        )
    ax.set(
        xticks=x,
        xticklabels=focal.gene,
        ylim=(0, 32),
        ylabel="Nuclei with peak detected (%)",
        title="COQ8A high / low accessibility",
    )
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax = axes[2]
    plotrows = []
    for i, gene in enumerate(["MYOD1", "CAV3"]):
        for j, (model, color) in enumerate(
            [("technical", "#7096AE"), ("technical_state_coq", "#BA622E")]
        ):
            r = counts[
                (counts.gene == gene)
                & (counts.peak == FOCAL[gene])
                & (counts.model == model)
            ].iloc[0]
            y = i * 2.2 + j * 0.65
            ax.errorbar(
                r.RNA_ratio,
                y,
                xerr=[[r.RNA_ratio - r.CI_low], [r.CI_high - r.RNA_ratio]],
                fmt="o",
                color=color,
                capsize=3,
                label=["Technical", "+ COQ8A and muscle score"][j] if i == 0 else None,
            )
            ax.text(2.02, y, f"{r.RNA_ratio:.2f}×", va="center", fontsize=8)
            plotrows.append(
                dict(
                    gene=gene,
                    model=model,
                    **{k: r[k] for k in ["RNA_ratio", "CI_low", "CI_high", "p"]},
                )
            )
    ax.axvline(1, color="#999999", ls=":", lw=1)
    ax.set(
        yticks=[0.325, 2.525],
        yticklabels=["MYOD1", "CAV3"],
        ylim=(3.2, -0.7),
        xlim=(0.96, 2.25),
        xticks=[1, 1.25, 1.5, 1.75, 2],
        xlabel="RNA count ratio: peak detected / undetected",
        title="Peak–RNA coupling",
    )
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    for ax, letter in zip(axes, "ABC"):
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.15, 1.08, letter, transform=ax.transAxes, weight="bold", fontsize=14)
    fig.subplots_adjust(left=0.065, right=0.98, bottom=0.20, top=0.83, wspace=0.42)
    folder = a.root / "figures/unstratified"
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        fig.savefig(folder / ("Unstratified_peak_RNA_exploration." + ext), dpi=220)
    pd.DataFrame(plotrows).to_csv(
        out / "count_figure_source.tsv", sep="\t", index=False
    )
    print(focal.to_string(index=False))
    print("Nominal positive chains:", len(coherent), "genes:", coherent.gene.unique())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path)
    p.add_argument("--external", type=Path)
    a = p.parse_args()
    out = a.root / "results/unstratified"
    out.mkdir(parents=True, exist_ok=True)
    if a.work:
        package(a, out)
    report(a, out)
