"""Integrate native-peak donor tests with RNA and render complete screens."""

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main(root, work):
    out = root / "results/differentiation_fusion25/native_GSE240061"
    figdir = root / "figures/differentiation_fusion25"
    read = lambda p: pd.read_csv(p, sep="\t")
    spec = importlib.util.spec_from_file_location(
        "union", root / "scripts/76_differentiation_fusion_analysis.py"
    )
    union = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(union)
    configs = read(
        root / "results/differentiation_fusion25/GSE240061_configurations.tsv"
    )
    membership = read(out / "native_peak_membership.tsv")
    names = (
        membership.groupby("peak")
        .gene.agg(lambda x: ";".join(sorted(set(x))))
        .rename("scope_genes")
    )
    tests = read(out / "ATAC_native_models.tsv.gz")
    assert not tests.duplicated(["config", "peak"]).any()
    for _, d in tests.groupby("config"):
        assert np.allclose(
            d.q_native_union, multipletests(d.PValue, method="fdr_bh")[1]
        )
        assert len(d) == d.n_native_tested.iloc[0]
    tests = tests.merge(configs.drop(columns="n_donors"), on="config").merge(
        names, on="peak"
    )
    tests["q_union"] = tests.q_native_union
    union.save(tests, out / "ATAC_all_settings_annotated.tsv.gz")
    union.save(tests[tests.primary], out / "ATAC_primary.tsv")
    summ = union.summary(tests, ["config"], "PValue").merge(
        configs, on="config", how="right"
    )
    summ = summ.merge(read(out / "model_status.tsv"), on="config", how="left")
    union.save(summ, out / "summary.tsv")
    rna = read(
        root / "results/differentiation_fusion25/GSE240061_RNA_all_settings.tsv.gz"
    )
    rna = rna[["config", "gene", "FC", "PValue", "q_union"]].rename(
        columns={"FC": "RNA_FC", "PValue": "RNA_p", "q_union": "RNA_q"}
    )
    links = read(out / "peak_RNA_links.tsv.gz")
    links = links[links.model.eq("technical_COQ")].rename(
        columns={
            "q_family": "link_q",
            "p_donor_t": "link_p",
            "n_donors": "link_n_donors",
        }
    )
    linked = tests.merge(membership[["peak", "gene"]].drop_duplicates(), on="peak")
    linked = linked.merge(
        links[
            [
                "peak",
                "gene",
                "celltype",
                "time",
                "r_equal_donor",
                "link_p",
                "link_q",
                "link_n_donors",
            ]
        ],
        on=["peak", "gene", "celltype", "time"],
        how="left",
    )
    linked = linked.merge(rna, on=["config", "gene"], how="left")
    union.save(linked, out / "ATAC_RNA_all_settings.tsv.gz")
    union.save(linked[linked.primary], out / "ATAC_RNA_primary.tsv")
    union.save(
        linked[linked.caliper.ne("none") & linked.q_union.lt(0.1)].sort_values(
            "PValue"
        ),
        out / "matched_q10_candidates.tsv",
    )

    columns = read(root / "results/replication_sensitivity/pseudobulk_columns.tsv")
    norm = read(root / "results/replication_sensitivity/ATAC_normalization_factors.tsv")
    columns = columns.merge(
        norm[["column", "TMM_factor"]], on="column", how="left", validate="one_to_one"
    )
    features = (work / "all_ATAC_features.txt").read_text().splitlines()
    mat = np.memmap(
        work / "sensitivity_ATAC_pseudobulk.bin",
        dtype="<i4",
        mode="r",
        shape=(len(columns), len(features)),
    )
    support = []
    for t in tests[tests.caliper.ne("none") & tests.q_union.lt(0.1)].itertuples():
        c = columns[columns.config.eq(t.config)].copy()
        c["peak"] = t.peak
        c["count"] = mat[c.column.to_numpy(), features.index(t.peak)]
        c["CPM"] = c["count"] / (c.raw_ATAC_size * c.TMM_factor) * 1e6
        support.append(c)
    union.save(
        pd.concat(support, ignore_index=True), out / "matched_q10_donor_counts.tsv"
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    def panel(ax, d, p, q, title):
        d = d[d.FC.gt(0) & np.isfinite(d.FC)]
        color = np.where(
            d[q].lt(0.1), "#AA3377", np.where(d[p].lt(0.05), "#0077AA", "#A5ADB4")
        )
        ax.scatter(
            np.log2(d.FC),
            -np.log10(d[p].clip(lower=1e-300)),
            c=color,
            s=16,
            alpha=0.8,
            rasterized=True,
        )
        ax.axvline(0, c="#777777", lw=0.7)
        ax.axhline(-np.log10(0.05), c="#777777", lw=0.7, ls=":")
        ax.set(
            title=title,
            xlabel="Accessibility log2 fold change (high / low)",
            ylabel="−log10(raw p)",
        )

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c, label=l)
        for c, l in [
            ("#A5ADB4", "p ≥ 0.05"),
            ("#0077AA", "p < 0.05"),
            ("#AA3377", "q < 0.10"),
        ]
    ]
    old = read(root / "results/differentiation_fusion25/GSE208248_primary.tsv")
    crossmap = read(
        root / "results/differentiation_fusion25/GSE240061_peak_membership.tsv"
    )[["author_peak", "target_peak"]].drop_duplicates()
    cross = tests[tests.caliper.ne("none") & tests.q_union.lt(0.1)].merge(
        crossmap, left_on="peak", right_on="author_peak", how="left"
    )
    oldcmp = old[["peak", "FC", "p_recomputed", "q_union"]].rename(
        columns={
            "peak": "target_peak",
            "FC": "old_primary_FC",
            "p_recomputed": "old_primary_p",
            "q_union": "old_primary_q",
        }
    )
    union.save(
        cross.merge(oldcmp, on="target_peak", how="left"),
        out / "matched_candidates_old_comparison.tsv",
    )
    slow = tests[tests.config.eq("C285")]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    panel(
        ax[0],
        old,
        "p_recomputed",
        "q_union",
        "GSE208248 · 565 peaks\n201 matched pairs · 2 source lines",
    )
    panel(
        ax[1],
        slow,
        "PValue",
        "q_union",
        f"GSE240061 · Slow, Post · {len(slow)} peaks\n307 matched pairs · 6 donors",
    )
    for a, label in zip(ax, ["A", "B"]):
        a.text(
            -0.12, 1.05, label, transform=a.transAxes, fontweight="bold", fontsize=15
        )
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"union_peak_screen.{ext}", dpi=230)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.3), layout="constrained")
    for ax, (cfg, d) in zip(axes.flat, tests[tests.primary].groupby("config")):
        panel(
            ax,
            d,
            "PValue",
            "q_union",
            f"{d.celltype.iloc[0]} · {d.time.iloc[0]} · {len(d)} peaks\n{d.n_donors.iloc[0]} donors · {d.n_high.iloc[0]} pairs",
        )
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"native_primary_contexts.{ext}", dpi=220)
    plt.close(fig)
    print(summ[summ.primary].to_string(index=False))
    show = [
        "config",
        "peak",
        "scope_genes",
        "FC",
        "PValue",
        "q_union",
        "n_donors",
        "n_high",
    ]
    print(
        "PRIMARY",
        tests[tests.primary & tests.q_union.lt(0.1)][show].to_string(index=False),
    )
    print(
        "MATCHED OPEN",
        tests[tests.caliper.ne("none") & tests.FC.gt(1) & tests.q_union.lt(0.1)]
        .sort_values("PValue")[show]
        .head(25)
        .to_string(index=False),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work)
