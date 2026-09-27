"""Validate the union analysis against matrices and render cohort summaries."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats, io


def main(root, oldwork, newwork):
    out = root / "results/differentiation_fusion25"
    figdir = root / "figures/differentiation_fusion25"
    figdir.mkdir(parents=True, exist_ok=True)
    read = lambda p: pd.read_csv(p, sep="\t")
    old = read(out / "GSE208248_primary.tsv")
    genes = read(out / "genes.tsv").gene.tolist()
    peaks = (oldwork / "peaks.txt").read_text().splitlines()
    rnagenes = (oldwork / "genes.txt").read_text().splitlines()
    pi = pd.Index(peaks).get_indexer(old.peak)
    gi = pd.Index(rnagenes).get_indexer(genes)
    assert min(pi) >= 0 and min(gi) >= 0
    pairs = read(root / "results/tables/matched_pairs.tsv.gz")
    pairs = pairs[pairs.gate.eq("TSS_ge_3") & pairs.contrast.eq("3plus_vs_1")]
    totals = {
        k: np.zeros(len(pi), dtype=int)
        for k in ["high_open", "low_open", "high_only", "low_only"]
    }
    modules, rhi, rlo = [], [], []
    for gsm, d in pairs.groupby("gsm"):
        meta = read(oldwork / f"{gsm}_meta.tsv").set_index("barcode")
        hi, lo = meta.index.get_indexer(d.high_barcode), meta.index.get_indexer(
            d.low_barcode
        )
        assert min(hi) >= 0 and min(lo) >= 0
        x = sparse.load_npz(oldwork / f"{gsm}_atac.npz")[pi]
        h, l = x[:, hi].toarray(), x[:, lo].toarray()
        for key, val in [
            ("high_open", h),
            ("low_open", l),
            ("high_only", h * (1 - l)),
            ("low_only", l * (1 - h)),
        ]:
            totals[key] += val.sum(axis=1).astype(int)
        modules.extend(
            dict(gsm=gsm, pair=k, high=h[:, k].mean(), low=l[:, k].mean())
            for k in range(len(d))
        )
        raw = sparse.load_npz(oldwork / f"{gsm}_rna.npz")[gi]
        rhi.append(
            raw[:, hi].toarray()
            / meta.iloc[hi].total_rna_umi.to_numpy()[None, :]
            * 10000
        )
        rlo.append(
            raw[:, lo].toarray()
            / meta.iloc[lo].total_rna_umi.to_numpy()[None, :]
            * 10000
        )
    for k, v in totals.items():
        assert np.array_equal(v, old[k].to_numpy()), k
    module = pd.DataFrame(modules)
    module.to_csv(out / "GSE208248_module_pairs.tsv", sep="\t", index=False)
    mh, ml = float(module.high.mean()), float(module.low.mean())
    module_summary = dict(
        n_peaks=len(pi),
        n_pairs=len(module),
        high_fraction=mh,
        low_fraction=ml,
        FC=mh / ml,
        difference_percentage_points=100 * (mh - ml),
        paired_nucleus_t_p=stats.ttest_rel(module.high, module.low).pvalue,
    )
    (out / "GSE208248_module_summary.json").write_text(
        json.dumps(module_summary, indent=2) + "\n"
    )
    mlibrary = module.groupby("gsm")[["high", "low"]].mean().reset_index()
    mlibrary["FC"] = mlibrary.high / mlibrary.low
    mlibrary.to_csv(out / "GSE208248_module_by_library.tsv", sep="\t", index=False)
    rh, rl = np.hstack(rhi).mean(1), np.hstack(rlo).mean(1)
    rna = pd.DataFrame(
        dict(
            gene=genes,
            high_mean_CP10k=rh,
            low_mean_CP10k=rl,
            RNA_CP10k_FC=np.divide(
                rh, rl, out=np.full(len(genes), np.nan), where=rl > 0
            ),
        )
    )
    rna = rna.merge(read(out / "GSE208248_RNA.tsv"), on="gene")
    rna.to_csv(out / "GSE208248_RNA_with_linear_FC.tsv", sep="\t", index=False)

    new = read(out / "GSE240061_all_settings.tsv.gz")
    targets = new[new.caliper.ne("none") & new.q_union.lt(0.1)]
    columns = read(root / "results/replication_sensitivity/pseudobulk_columns.tsv")
    features = (newwork / "all_ATAC_features.txt").read_text().splitlines()
    matrix = np.memmap(
        newwork / "sensitivity_ATAC_pseudobulk.bin",
        dtype="<i4",
        mode="r",
        shape=(len(columns), len(features)),
    )
    nf = read(root / "results/replication_sensitivity/ATAC_normalization_factors.tsv")
    columns = columns.merge(
        nf[["column", "TMM_factor"]], on="column", how="left", validate="one_to_one"
    )
    support = []
    for d in targets.itertuples():
        ix = features.index(d.peak)
        c = columns[columns.config.eq(d.config)].copy()
        c["peak"] = d.peak
        c["count"] = matrix[c.column.to_numpy(), ix]
        c["CPM"] = c["count"] / (c.raw_ATAC_size * c.TMM_factor) * 1e6
        support.append(c)
    raw = pd.concat(support, ignore_index=True)
    raw.to_csv(out / "GSE240061_q10_matched_donor_counts.tsv", sep="\t", index=False)
    rfeatures = (newwork / "candidate_RNA_genes.txt").read_text().splitlines()
    rcounts = io.mmread(newwork / "sensitivity_RNA_pseudobulk.mtx").tocsr()
    c058 = columns[columns.config.eq("C058")][
        ["donor", "group", "n_nuclei", "column"]
    ].copy()
    c058["MYOD1_RNA_count"] = (
        rcounts[rfeatures.index("MYOD1"), c058.column.to_numpy()].toarray().ravel()
    )
    c058.to_csv(out / "GSE240061_C058_MYOD1_RNA_counts.tsv", sep="\t", index=False)

    ann = read(root / "results/candidate_identity/all_nearby_transcript_TSS.tsv")
    close = ann[ann.peak.eq("chr16:1311478-1312392")].copy()
    close["new_peak"] = "chr16:1308794-1309860"
    close["new_midpoint_distance_bp"] = (
        close.tss_0based - ((1308794 + 1309860) / 2)
    ).abs()
    close["new_peak_contains_TSS"] = close.tss_0based.between(1308794, 1309859)
    close.sort_values("new_midpoint_distance_bp").to_csv(
        out / "GSE240061_primary_closed_peak_TSS.tsv", sep="\t", index=False
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    subsets = [
        (old, "GSE208248 · 565 peaks", "p_recomputed"),
        (
            new[new.config.eq("C285")],
            "GSE240061 · Slow, Post · 261 tested peaks",
            "PValue",
        ),
    ]
    for ax, (d, title, p) in zip(axes, subsets):
        finite = d.FC.gt(0) & np.isfinite(d.FC)
        d = d[finite]
        colors = np.where(
            d.q_union.lt(0.1), "#AA3377", np.where(d[p].lt(0.05), "#0077AA", "#A5ADB4")
        )
        ax.scatter(
            np.log2(d.FC),
            -np.log10(d[p].clip(lower=1e-300)),
            c=colors,
            s=18,
            alpha=0.8,
            rasterized=True,
        )
        ax.axvline(0, c="#777777", lw=0.7)
        ax.axhline(-np.log10(0.05), c="#777777", lw=0.7, ls=":")
        ax.set(
            xlabel="Accessibility log2 fold change (COQ8A high / low)",
            ylabel="−log10(raw p)",
            title=title,
        )
    axes[0].set_title(
        "GSE208248 · 565 peaks\n201 matched pairs · 2 source lines", fontsize=11
    )
    axes[1].set_title(
        "GSE240061 · Slow, Post · 261 tested peaks\n307 matched pairs · 6 donors",
        fontsize=11,
    )
    for ax, letter in zip(axes, ["A", "B"]):
        ax.text(
            -0.12, 1.05, letter, transform=ax.transAxes, fontweight="bold", fontsize=15
        )
    from matplotlib.lines import Line2D

    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor=c, label=l)
            for c, l in [
                ("#A5ADB4", "p ≥ 0.05"),
                ("#0077AA", "p < 0.05"),
                ("#AA3377", "q < 0.10"),
            ]
        ],
        loc="outside lower center",
        ncol=3,
        frameon=False,
    )
    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"mapped_union_peak_screen.{ext}", dpi=230)
    plt.close(fig)
    (out / "validation.json").write_text(
        json.dumps(
            dict(
                old_primary_matrix_counts_reproduced=len(pi),
                old_rna_genes=len(rna),
                old_pairs=len(module),
                new_matched_q10_count_rows=len(raw),
            ),
            indent=2,
        )
        + "\n"
    )
    print(module_summary)
    print(mlibrary.to_string(index=False))
    print(close.sort_values("new_midpoint_distance_bp").head(3).to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--old-work", type=Path, required=True)
    p.add_argument("--new-work", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.old_work, a.new_work)
