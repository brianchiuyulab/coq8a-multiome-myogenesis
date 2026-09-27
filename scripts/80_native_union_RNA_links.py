"""Test native union peak–RNA pairs within donor and cell-type/time strata."""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse, stats
from statsmodels.stats.multitest import multipletests


def main(root, work):
    out = root / "results/differentiation_fusion25/native_GSE240061"
    spec = importlib.util.spec_from_file_location(
        "links", root / "scripts/69_replication_continuous_links.py"
    )
    links = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(links)
    pairs = pd.read_csv(out / "native_peak_membership.tsv", sep="\t")[
        ["peak", "gene", "index"]
    ].drop_duplicates()
    native = pairs[["peak", "index"]].drop_duplicates().sort_values("index")
    pnames = native.peak.tolist()
    shape = tuple(map(int, (work / "all_ATAC_shape.txt").read_text().split()))
    x = np.memmap(work / "all_ATAC_x.bin", dtype="<i4", mode="r")
    i = np.memmap(work / "all_ATAC_i.bin", dtype="<i4", mode="r")
    p = np.memmap(work / "all_ATAC_p.bin", dtype="<i4", mode="r")
    full = sparse.csc_matrix((x, i, p), shape=shape, copy=False)
    atac = full[native["index"].to_numpy(), :].tocsc()
    atac.data[:] = 1
    pair_table = pd.read_csv(work / "sensitivity_pairs.tsv.gz", sep="\t")
    nucleus_score = np.asarray(atac.mean(axis=0)).ravel()
    pair_table["high_fraction"] = nucleus_score[pair_table.high_index.to_numpy()]
    pair_table["low_fraction"] = nucleus_score[pair_table.low_index.to_numpy()]
    module_donors = (
        pair_table.groupby(["config", "donor"])[["high_fraction", "low_fraction"]]
        .mean()
        .reset_index()
    )
    module_donors.to_csv(out / "module_by_donor.tsv", sep="\t", index=False)
    module_rows = []
    for cfg, d in module_donors.groupby("config"):
        if len(d) < 2:
            continue
        module_rows.append(
            dict(
                config=cfg,
                n_donors=len(d),
                high_fraction=d.high_fraction.mean(),
                low_fraction=d.low_fraction.mean(),
                FC=d.high_fraction.mean() / d.low_fraction.mean(),
                p_donor=stats.ttest_rel(d.high_fraction, d.low_fraction).pvalue,
            )
        )
    modules = pd.DataFrame(module_rows).merge(
        pd.read_csv(
            root / "results/differentiation_fusion25/GSE240061_configurations.tsv",
            sep="\t",
        ).drop(columns="n_donors"),
        on="config",
    )
    for _, d in modules.groupby(["qc", "rule", "caliper"]):
        valid = d.p_donor.notna()
        modules.loc[d.index[valid], "q_contexts"] = multipletests(
            d.loc[valid, "p_donor"], method="fdr_bh"
        )[1]
    modules.to_csv(out / "module_effects.tsv", sep="\t", index=False)
    bcs = (work / "all_ATAC_barcodes.txt").read_text().splitlines()
    rbcs = (work / "barcodes.txt").read_text().splitlines()
    ridx = pd.Index(rbcs).get_indexer(bcs)
    assert min(ridx) >= 0
    genes = (work / "candidate_RNA_genes.txt").read_text().splitlines()
    rna = io.mmread(work / "candidate_RNA_counts.mtx").tocsc()[:, ridx]
    meta = (
        pd.read_csv(work / "nucleus_inventory.tsv", sep="\t")
        .set_index("barcode")
        .loc[bcs]
    )
    meta["donor"] = meta.Sample.str.split("_").str[0]
    meta["celltype"] = (
        meta.refined_annotations_wknn_0_8
        if "refined_annotations_wknn_0_8" in meta
        else meta["refined_annotations_wknn_0.8"]
    )
    prog = set(
        (root / "results/unstratified/programme_genes.txt").read_text().splitlines()
    ) & set(genes)
    progi = [genes.index(g) for g in prog]
    groups = {
        g: pd.Index(pnames).get_indexer(d.peak)
        for g, d in pairs.groupby("gene")
        if g in genes
    }
    rows = []
    for (ct, time, donor), d in meta[
        meta.celltype.isin(["Satellite Cells", "Fast", "Slow", "Intermediate"])
    ].groupby(["celltype", "Time", "donor"]):
        if len(d) < 8:
            continue
        idx = meta.index.get_indexer(d.index)
        counts = rna[:, idx].toarray().T.astype(float)
        y = np.log1p(counts / d.raw_RNA_library_size.to_numpy()[:, None] * 10000)
        a = atac[:, idx].toarray().T.astype(float)
        state = y[:, progi].mean(1)
        coq = y[:, genes.index("COQ8A")]
        tech = np.column_stack(
            [
                np.ones(len(d)),
                np.log1p(d.raw_RNA_library_size),
                np.log1p(d.nFeature_ATAC),
                d["percent.mt"],
                coq,
            ]
        )
        for gene, js in groups.items():
            j = genes.index(gene)
            if np.count_nonzero(counts[:, j]) < max(3, np.ceil(0.01 * len(d))):
                continue
            yg = y[:, j]
            adjusted_state = (
                (state * len(prog) - yg) / (len(prog) - 1) if gene in prog else state
            )
            for model, z in [
                ("technical_COQ", tech),
                ("technical_COQ_state", np.column_stack([tech, adjusted_state])),
            ]:
                rr = links.correlation(
                    links.residual(a[:, js], z), links.residual(yg, z)[:, None]
                )
                good = (
                    (a[:, js].sum(0) >= 3)
                    & ((len(d) - a[:, js].sum(0)) >= 3)
                    & np.isfinite(rr)
                )
                rows.extend(
                    dict(
                        celltype=ct,
                        time=time,
                        donor=donor,
                        model=model,
                        peak=pnames[js[k]],
                        gene=gene,
                        n_nuclei=len(d),
                        r=rr[k],
                    )
                    for k in np.flatnonzero(good)
                )
        print(ct, time, donor, len(d), flush=True)
    per = pd.DataFrame(rows)
    per.to_csv(
        out / "peak_RNA_by_donor.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    combined = links.combine(per, ["celltype", "time", "model", "peak", "gene"])
    combined.to_csv(
        out / "peak_RNA_links.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work)
