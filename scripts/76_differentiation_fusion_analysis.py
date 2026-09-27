"""Analyze the frozen differentiation/fusion union in both multiome cohorts.

Existing per-peak models retain their original inputs and normalization. This
script recomputes the old paired tests from discordant counts, establishes the
union test families, and combines accessibility with separately measured RNA.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def read(path):
    return pd.read_csv(path, sep="\t")


def save(frame, path):
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        compression=(
            {"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None
        ),
    )


def adjust(frame, keys, p):
    frame = frame.copy()
    frame["q_union"] = np.nan
    frame["n_tests_union"] = 0
    for _, d in frame.groupby(keys, dropna=False):
        valid = d[p].notna()
        idx = d.index[valid]
        frame.loc[idx, "q_union"] = multipletests(d.loc[idx, p], method="fdr_bh")[1]
        frame.loc[d.index, "n_tests_union"] = len(idx)
    return frame


def summary(frame, keys, p):
    rows = []
    for k, d in frame.groupby(keys):
        k = k if isinstance(k, tuple) else (k,)
        rows.append(
            dict(zip(keys, k))
            | dict(
                n_peaks=len(d),
                nominal_open=int(((d.FC > 1) & (d[p] < 0.05)).sum()),
                nominal_closed=int(((d.FC < 1) & (d[p] < 0.05)).sum()),
                q10_open=int(((d.FC > 1) & (d.q_union < 0.1)).sum()),
                q10_closed=int(((d.FC < 1) & (d.q_union < 0.1)).sum()),
                q05_open=int(((d.FC > 1) & (d.q_union < 0.05)).sum()),
                q05_closed=int(((d.FC < 1) & (d.q_union < 0.05)).sum()),
                min_p=d[p].min(),
                min_q=d.q_union.min(),
            )
        )
    return pd.DataFrame(rows)


def main(root, work):
    out = root / "results/differentiation_fusion25"
    out.mkdir(parents=True, exist_ok=True)
    members = read(
        root / "results/chr16_RNA_followup/differentiation_fusion_membership.tsv"
    )
    genes = set(members.gene)
    assert len(genes) == 25
    save(members, out / "genes.tsv")
    near = read(root / "results/tables/candidate_peak_gene.tsv.gz")
    oldmap = near[
        near.n_libraries.eq(4) & near.gene.isin(genes) & near.nearest_tss_bp.le(100000)
    ].copy()
    save(oldmap, out / "GSE208248_peak_membership.tsv")
    names = (
        oldmap.groupby("peak")
        .gene.agg(lambda x: ";".join(sorted(set(x))))
        .rename("scope_genes")
    )
    old = read(root / "results/tables/candidate_peak_effects_pooled.tsv.gz")
    old = old[old.peak.isin(names.index)].copy()
    old["p_recomputed"] = [
        binomtest(int(h), int(h + l), 0.5).pvalue if h + l else 1.0
        for h, l in zip(old.high_only, old.low_only)
    ]
    assert np.allclose(old.p_recomputed, old.p_pair_binomial, atol=1e-12)
    old["FC"] = old.high_open / old.low_open.replace(0, np.nan)
    old["high_open_fraction"] = old.high_open / old.n_pairs
    old["low_open_fraction"] = old.low_open / old.n_pairs
    old = adjust(old, ["gate", "contrast"], "p_recomputed").merge(names, on="peak")
    save(old, out / "GSE208248_all_settings.tsv.gz")
    anchor = old[old.gate.eq("TSS_ge_3") & old.contrast.eq("3plus_vs_1")].sort_values(
        "p_recomputed"
    )
    save(anchor, out / "GSE208248_primary.tsv")
    save(
        summary(old, ["gate", "contrast"], "p_recomputed"),
        out / "GSE208248_summary.tsv",
    )
    bylib = read(root / "results/tables/candidate_peak_effects_by_library.tsv.gz")
    save(bylib[bylib.peak.isin(names.index)], out / "GSE208248_by_library.tsv.gz")
    rna = read(root / "results/unstratified/RNA_COQ8A_high_low.tsv")
    rna = rna[rna.gene.isin(genes)].copy()
    rna["q_union"] = multipletests(rna.p_pair, method="fdr_bh")[1]
    save(rna, out / "GSE208248_RNA.tsv")
    links = read(root / "results/unstratified/correlations_combined.tsv.gz")
    links = links[links.peak.isin(names.index)].copy()
    links = adjust(links, ["model", "detection_fraction"], "p")
    save(links, out / "GSE208248_all_local_RNA_links.tsv.gz")
    chain = oldmap.merge(anchor, on="peak", suffixes=("_membership", ""))
    link = links[links.model.eq("technical_coq") & links.detection_fraction.eq(0.05)][
        ["peak", "gene", "r", "p", "q_union", "n_libraries"]
    ]
    link = link.rename(
        columns={"p": "link_p", "q_union": "link_q", "n_libraries": "link_n_libraries"}
    )
    rr = rna.rename(
        columns={
            "p_pair": "RNA_p",
            "q_union": "RNA_q",
            "difference": "RNA_log1p_difference",
        }
    )
    chain = chain.merge(link, on=["peak", "gene"], how="left").merge(
        rr[["gene", "RNA_p", "RNA_q", "RNA_log1p_difference"]], on="gene", how="left"
    )
    sensitivity_link = links[
        links.model.eq("technical_coq") & links.detection_fraction.eq(0.01)
    ][["peak", "gene", "r", "p", "q_union"]]
    sensitivity_link = sensitivity_link.rename(
        columns={"r": "link_1pct_r", "p": "link_1pct_p", "q_union": "link_1pct_q"}
    )
    chain = chain.merge(
        sensitivity_link, on=["peak", "gene"], how="left", validate="one_to_one"
    )
    save(chain.sort_values("p_recomputed"), out / "GSE208248_primary_ATAC_RNA.tsv")

    overlap = read(work / "candidate_peak_overlap.tsv")
    overlap = overlap[
        overlap.target_overlap_fraction.ge(0.5)
        & overlap.author_overlap_fraction.ge(0.5)
    ]
    newmap = overlap.merge(oldmap, left_on="target_peak", right_on="peak").drop(
        columns="peak"
    )
    save(newmap, out / "GSE240061_peak_membership.tsv")
    newnames = (
        newmap.groupby("author_peak")
        .gene.agg(lambda x: ";".join(sorted(set(x))))
        .rename("scope_genes")
    )
    configs = read(root / "results/replication_sensitivity/configurations.tsv")
    configs["primary"] = (
        configs.qc.eq("author")
        & configs.rule.eq("3plus_vs_1")
        & configs.caliper.eq("0.1")
    )
    save(configs, out / "GSE240061_configurations.tsv")
    new = read(root / "results/replication_sensitivity/ATAC_donor_models.tsv.gz")
    new = new[new.peak.isin(newnames.index)].copy()
    new = adjust(new, ["config"], "PValue").merge(
        newnames, left_on="peak", right_index=True
    )
    new = new.merge(configs.drop(columns="n_donors"), on="config")
    save(new, out / "GSE240061_all_settings.tsv.gz")
    primary = new[new.primary].sort_values(["config", "PValue"])
    save(primary, out / "GSE240061_primary.tsv")
    summ = summary(new, ["config"], "PValue").merge(configs, on="config", how="right")
    status = read(root / "results/replication_sensitivity/model_status.tsv")
    summ = summ.merge(status[["config", "status"]], on="config", how="left")
    save(summ, out / "GSE240061_summary.tsv")
    nrna = read(root / "results/replication_sensitivity/RNA_donor_models.tsv.gz")
    nrna = adjust(nrna[nrna.gene.isin(genes)], ["config"], "PValue")
    save(nrna, out / "GSE240061_RNA_all_settings.tsv.gz")
    nl = read(root / "results/replication_links/peak_RNA_links.tsv.gz")
    nl = nl[nl.peak.isin(newnames.index)].copy()
    nl = adjust(nl, ["celltype", "time", "model"], "p_donor_t")
    save(nl, out / "GSE240061_RNA_links.tsv.gz")
    ng = (
        newmap[["author_peak", "gene"]]
        .drop_duplicates()
        .rename(columns={"author_peak": "peak"})
    )
    nc = new.merge(ng, on="peak")
    nlr = nl[nl.model.eq("technical_COQ")][
        [
            "peak",
            "gene",
            "celltype",
            "time",
            "r_equal_donor",
            "p_donor_t",
            "q_union",
            "n_donors",
        ]
    ]
    nlr = nlr.rename(
        columns={
            "q_union": "link_q",
            "p_donor_t": "link_p",
            "n_donors": "link_n_donors",
        }
    )
    nrr = nrna[["config", "gene", "FC", "PValue", "q_union"]].rename(
        columns={"FC": "RNA_FC", "PValue": "RNA_p", "q_union": "RNA_q"}
    )
    nc = nc.merge(nlr, on=["peak", "gene", "celltype", "time"], how="left").merge(
        nrr, on=["config", "gene"], how="left"
    )
    save(nc, out / "GSE240061_ATAC_RNA_all_settings.tsv.gz")
    save(nc[nc.primary], out / "GSE240061_primary_ATAC_RNA.tsv")
    save(
        nc[(nc.FC > 1) & (nc.q_union < 0.1)].sort_values("PValue"),
        out / "GSE240061_open_q10_sensitivity.tsv",
    )
    manifest = dict(
        n_genes=len(genes),
        old_peaks=len(names),
        new_mapped_peaks=len(newnames),
        old_primary_pairs=int(anchor.n_pairs.iloc[0]),
        old_primary_sources=2,
        new_primary_configs=configs.loc[configs.primary, "config"].tolist(),
        new_estimable_configs=int(new.config.nunique()),
        old_q="BH over all distinct union peaks including non-discordant peaks; per setting",
        new_q="BH over count-supported mapped union peaks; per cell type/time/configuration",
        RNA_q="BH over measured union genes; per configuration",
        links="Old: all available local targets +/-500kb; new: stored 221-gene-neighborhood pairs only; not all possible targets",
        selection="Union agreed after earlier exploration; prospective for this report, not preregistered or independent confirmation",
        computation="Old p recomputed from paired discordant counts. New fitted coefficients/p retained; union BH recomputed. Dispersion and normalization background unchanged.",
    )
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(
        "OLD",
        anchor[["peak", "scope_genes", "FC", "p_recomputed", "q_union"]]
        .head(12)
        .to_string(index=False),
    )
    print("NEW PRIMARY", summ[summ.primary].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.work)
