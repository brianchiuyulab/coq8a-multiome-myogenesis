"""Donor-level continuous COQ8A accessibility and candidate RNA-link analyses."""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import io, stats
from statsmodels.stats.multitest import multipletests


def residual(x, z):
    return x - z @ np.linalg.lstsq(z, x, rcond=None)[0]


def correlation(x, y):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (x * y).sum(axis=0) / np.sqrt((x * x).sum(axis=0) * (y * y).sum(axis=0))


def combine(frame, keys):
    rows = []
    for k, d in frame.groupby(keys):
        rs = d.r.to_numpy()
        rs = rs[np.isfinite(rs)]
        n = len(rs)
        if n < 2:
            continue
        z = np.arctanh(np.clip(rs, -0.999999, 0.999999))
        se = z.std(ddof=1) / np.sqrt(n)
        p = (
            2 * stats.t.sf(abs(z.mean() / se), n - 1)
            if se > 0
            else (1.0 if z.mean() == 0 else 0.0)
        )
        loo = [np.tanh(np.delete(z, i).mean()) for i in range(n)]
        rows.append(
            dict(zip(keys, k if isinstance(k, tuple) else (k,)))
            | dict(
                n_donors=n,
                r_equal_donor=np.tanh(z.mean()),
                p_donor_t=p,
                positive_donors=int((rs > 0).sum()),
                negative_donors=int((rs < 0).sum()),
                loo_r_min=min(loo),
                loo_r_max=max(loo),
            )
        )
    result = pd.DataFrame(rows)
    if len(result):
        for _, d in result.groupby(["celltype", "time", "model"]):
            result.loc[d.index, "q_family"] = multipletests(
                d.p_donor_t, method="fdr_bh"
            )[1]
    return result


def main(root, work, out, minimum_events=3):
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(work / "nucleus_inventory.tsv", sep="\t").set_index("barcode")
    bcs = (work / "barcodes.txt").read_text().splitlines()
    meta = meta.loc[bcs]
    genes = (work / "candidate_RNA_genes.txt").read_text().splitlines()
    peaks = (work / "candidate_ATAC_peaks.txt").read_text().splitlines()
    rna = io.mmread(work / "candidate_RNA_counts.mtx").tocsc()
    atac = io.mmread(work / "candidate_ATAC_counts.mtx").tocsc()
    atac.data[:] = 1
    overlap = pd.read_csv(work / "candidate_peak_overlap.tsv", sep="\t")
    overlap = overlap[
        overlap.target_overlap_fraction.ge(0.5)
        & overlap.author_overlap_fraction.ge(0.5)
    ]
    near = pd.read_csv(root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
    pairs = overlap.merge(
        near[["peak", "gene"]], left_on="target_peak", right_on="peak"
    )[["author_peak", "gene"]].drop_duplicates()
    pairs = pairs[pairs.gene.isin(genes)]
    pnames = sorted(pairs.author_peak.unique())
    pidx = pd.Index(peaks).get_indexer(pnames)
    pair_groups = {
        g: pd.Index(pnames).get_indexer(d.author_peak) for g, d in pairs.groupby("gene")
    }
    prog = set(
        (root / "results/unstratified/programme_genes.txt").read_text().splitlines()
    ) & set(genes)
    progidx = [genes.index(g) for g in prog]
    meta["donor"] = meta.Sample.str.split("_").str[0]
    meta["celltype"] = meta["refined_annotations_wknn_0.8"]
    coqrows = []
    linkrows = []
    for (ct, time, donor), d in meta[
        meta.celltype.isin(["Satellite Cells", "Fast", "Slow", "Intermediate"])
    ].groupby(["celltype", "Time", "donor"]):
        idx = meta.index.get_indexer(d.index)
        n = len(d)
        if n < 8:
            continue
        counts = rna[:, idx].toarray().T.astype(float)
        y = np.log1p(counts / d.raw_RNA_library_size.to_numpy()[:, None] * 10000)
        state = y[:, progidx].mean(1)
        coq = y[:, genes.index("COQ8A")]
        a = atac[pidx][:, idx].toarray().T.astype(float)
        tech = np.column_stack(
            [
                np.ones(n),
                np.log1p(d.raw_RNA_library_size),
                np.log1p(d.nFeature_ATAC),
                d["percent.mt"],
            ]
        )
        for model, z in [
            ("technical", tech),
            ("technical_state", np.column_stack([tech, state])),
        ]:
            coqr = residual(coq, z)[:, None]
            ar = residual(a, z)
            rs = correlation(ar, coqr)
            eligible = (
                (a.sum(0) >= minimum_events)
                & ((n - a.sum(0)) >= minimum_events)
                & (np.count_nonzero(coq) >= minimum_events)
            )
            coqrows.extend(
                dict(
                    celltype=ct,
                    time=time,
                    donor=donor,
                    model=model,
                    peak=pnames[j],
                    n_nuclei=n,
                    n_covariates=z.shape[1],
                    n_COQ_detected=int(np.count_nonzero(coq)),
                    r=rs[j],
                )
                for j in np.flatnonzero(eligible)
            )
        for gene, js in pair_groups.items():
            j = genes.index(gene)
            if (counts[:, j] > 0).sum() < max(minimum_events, np.ceil(0.01 * n)):
                continue
            yg = y[:, j]
            leave_state = (
                (state * len(prog) - yg) / (len(prog) - 1) if gene in prog else state
            )
            for model, z in [
                ("technical_COQ", np.column_stack([tech, coq])),
                ("technical_COQ_state", np.column_stack([tech, coq, leave_state])),
            ]:
                yr = residual(yg, z)[:, None]
                ar = residual(a[:, js], z)
                rs = correlation(ar, yr)
                valid = (a[:, js].sum(0) >= minimum_events) & (
                    (n - a[:, js].sum(0)) >= minimum_events
                )
                linkrows.extend(
                    dict(
                        celltype=ct,
                        time=time,
                        donor=donor,
                        model=model,
                        peak=pnames[pj],
                        gene=gene,
                        n_nuclei=n,
                        n_covariates=z.shape[1],
                        r=rs[k],
                    )
                    for k, pj in enumerate(js)
                    if valid[k]
                )
        print(ct, time, donor, n, flush=True)
    coq = pd.DataFrame(coqrows)
    links = pd.DataFrame(linkrows)
    coq.to_csv(
        out / "COQ_ATAC_by_donor.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    links.to_csv(
        out / "peak_RNA_by_donor.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    combine(coq, ["celltype", "time", "model", "peak"]).to_csv(
        out / "COQ_ATAC_continuous.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    combine(links, ["celltype", "time", "model", "peak", "gene"]).to_csv(
        out / "peak_RNA_links.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--minimum-events", type=int, default=3)
    a = p.parse_args()
    main(a.root, a.work, a.out, a.minimum_events)
