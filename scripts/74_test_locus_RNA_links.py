"""Test all available local RNA targets of the chr16 peak in fixed contexts."""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import io, sparse, stats
from statsmodels.stats.multitest import multipletests


def main(root, work, locus, out):
    out.mkdir(parents=True, exist_ok=True)
    genes = (locus / "locus_RNA_genes.txt").read_text().splitlines()
    bcs = (work / "all_ATAC_barcodes.txt").read_text().splitlines()
    rb = (locus / "locus_RNA_barcodes.txt").read_text().splitlines()
    rna = io.mmread(locus / "locus_RNA.mtx").tocsr()[:, pd.Index(rb).get_indexer(bcs)]
    membership = sparse.load_npz(work / "sensitivity_membership.npz")
    io.mmwrite(locus / "locus_RNA_pseudobulk.mtx", rna @ membership)
    meta = (
        pd.read_csv(work / "nucleus_inventory.tsv", sep="\t")
        .set_index("barcode")
        .loc[bcs]
    )
    meta["donor"] = meta.Sample.str.split("_").str[0]
    cg = (work / "candidate_RNA_genes.txt").read_text().splitlines()
    cb = (work / "barcodes.txt").read_text().splitlines()
    order = pd.Index(cb).get_indexer(bcs)
    cr = io.mmread(work / "candidate_RNA_counts.mtx").tocsr()[:, order]
    peaks = (work / "candidate_ATAC_peaks.txt").read_text().splitlines()
    peak = "chr16-1311478-1312392"
    ac = (
        io.mmread(work / "candidate_ATAC_counts.mtx")
        .tocsr()[peaks.index(peak), order]
        .toarray()
        .ravel()
    )
    opened = (ac > 0).astype(float)
    pairs = pd.read_csv(work / "sensitivity_pairs.tsv.gz", sep="\t")
    contexts = {
        "Fast_Post_all_author": np.flatnonzero(
            meta["refined_annotations_wknn_0.8"].eq("Fast") & meta.Time.eq("Post")
        )
    }
    for cfg in ["C197", "C198", "C199"]:
        p = pairs[pairs.config.eq(cfg)]
        contexts[cfg] = np.unique(np.r_[p.high_index, p.low_index])
    pg = set(
        (root / "results/unstratified/programme_genes.txt").read_text().splitlines()
    ) & set(cg)
    programme = np.log1p(
        cr[[cg.index(g) for g in pg]].toarray().T
        / meta.raw_RNA_library_size.to_numpy()[:, None]
        * 10000
    ).mean(1)
    coq = np.log1p(
        cr[cg.index("COQ8A")].toarray().ravel()
        / meta.raw_RNA_library_size.to_numpy()
        * 10000
    )
    y = np.log1p(
        rna.toarray().T / meta.raw_RNA_library_size.to_numpy()[:, None] * 10000
    )
    counts = rna.toarray().T
    rows, availability = [], []
    for context, ix in contexts.items():
        for donor, d in meta.iloc[ix].groupby("donor"):
            ii = meta.index.get_indexer(d.index)
            n = len(ii)
            z = np.column_stack(
                [
                    np.ones(n),
                    np.log1p(d.raw_RNA_library_size),
                    np.log1p(d.nFeature_ATAC),
                    d["percent.mt"],
                    coq[ii],
                ]
            )
            for j, gene in enumerate(genes):
                detected = int((counts[ii, j] > 0).sum())
                availability.append(
                    dict(
                        context=context,
                        donor=donor,
                        gene=gene,
                        n_nuclei=n,
                        RNA_detected=detected,
                        RNA_counts=int(counts[ii, j].sum()),
                        peak_detected=int(opened[ii].sum()),
                    )
                )
                if (
                    n < 8
                    or detected < max(3, np.ceil(0.01 * n))
                    or min(opened[ii].sum(), n - opened[ii].sum()) < 3
                ):
                    continue
                for model, cov in [
                    ("technical_COQ", z),
                    ("technical_COQ_programme", np.column_stack([z, programme[ii]])),
                ]:
                    if n <= np.linalg.matrix_rank(cov) + 2:
                        continue
                    ar = (
                        opened[ii]
                        - cov @ np.linalg.lstsq(cov, opened[ii], rcond=None)[0]
                    )
                    yr = y[ii, j] - cov @ np.linalg.lstsq(cov, y[ii, j], rcond=None)[0]
                    r = float(np.corrcoef(ar, yr)[0, 1])
                    rows.append(
                        dict(
                            context=context,
                            donor=donor,
                            model=model,
                            gene=gene,
                            n_nuclei=n,
                            r=r,
                        )
                    )
    donor = pd.DataFrame(rows)
    donor.to_csv(out / "peak_RNA_by_donor.tsv", sep="\t", index=False)
    pd.DataFrame(availability).to_csv(
        out / "RNA_detection_by_donor.tsv", sep="\t", index=False
    )
    results = []
    for keys, d in donor.groupby(["context", "model", "gene"]):
        d = d[np.isfinite(d.r)]
        if len(d) < 2:
            continue
        z = np.arctanh(np.clip(d.r, -0.999999, 0.999999))
        se = z.std(ddof=1) / np.sqrt(len(z))
        p = 2 * stats.t.sf(abs(z.mean() / se), len(z) - 1) if se else 1.0
        results.append(
            dict(zip(["context", "model", "gene"], keys))
            | dict(
                n_donors=len(z),
                r_equal_donor=np.tanh(z.mean()),
                p=p,
                positive_donors=int((d.r > 0).sum()),
            )
        )
    result = pd.DataFrame(results)
    for _, d in result.groupby(["context", "model"]):
        result.loc[d.index, "q_local_targets"] = multipletests(d.p, method="fdr_bh")[1]
        result.loc[d.index, "n_local_targets"] = len(d)
    result.to_csv(out / "peak_RNA_links.tsv", sep="\t", index=False)
    inv = pd.read_csv(locus / "locus_RNA_inventory.tsv", sep="\t")
    inv.to_csv(out / "locus_RNA_inventory.tsv", sep="\t", index=False)
    print(
        result[result.context.eq("C198")]
        .sort_values(["model", "p"])
        .to_string(index=False)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--locus", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work, a.locus, a.out)
