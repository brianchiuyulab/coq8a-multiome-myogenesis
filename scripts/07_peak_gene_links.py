"""Explore all candidate peak–gene RNA links within the 221-gene search space.

Primary gate/contrast only. Links are conditioned on COQ8A group and both
depths, then combined across libraries. Same-data links are exploratory.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, read_barcodes


def residualize(x, design):
    return x - design @ np.linalg.lstsq(design, x, rcond=None)[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--tables", type=Path, required=True)
    p.add_argument("--motif-scores", type=Path)
    p.add_argument("--gate", choices=["TSS_ge_2", "TSS_ge_3"], default="TSS_ge_3")
    a = p.parse_args()
    suffix = "_tss2" if a.gate == "TSS_ge_2" else ""
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == a.gate) & (pairs.contrast == "2plus_vs_1")]
    cells = pd.read_csv(a.tables / "qc_nuclei.tsv.gz", sep="\t")
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    consensus = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t").set_index("peak")
    candidate_rows = []
    for gsm, ps in pairs.groupby("gsm"):
        n = len(ps)
        barcodes = ps.high_barcode.tolist() + ps.low_barcode.tolist()
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            R, P, rna_names, peak_names = read_barcodes(h5, barcodes)
        rna_idx = {g: i for i, g in enumerate(rna_names)}
        peak_idx = {str(q): i for i, q in enumerate(peak_names)}
        q = candidates.copy()
        target = q.peak if gsm == "GSM6339597" else q.peak.map(consensus[gsm + "_peak"])
        q["local_peak"] = target.map(peak_idx)
        q["local_gene"] = q.gene.map(rna_idx)
        q = q.dropna(subset=["local_peak", "local_gene"]).copy()
        q[["local_peak", "local_gene"]] = q[["local_peak", "local_gene"]].astype(int)
        unique_peaks = np.sort(q.local_peak.unique())
        unique_genes = np.sort(q.local_gene.unique())
        peak_lookup = {v: j for j, v in enumerate(unique_peaks)}
        gene_lookup = {v: j for j, v in enumerate(unique_genes)}
        X = P[:, unique_peaks].toarray().astype(np.float32)
        raw = R[:, unique_genes].toarray().astype(np.float32)
        depth = (
            cells[cells.gsm == gsm].set_index("barcode").loc[barcodes, "total_rna_umi"].to_numpy()
        )
        atac_depth = (
            cells[cells.gsm == gsm]
            .set_index("barcode")
            .loc[barcodes, "total_open_peaks"]
            .to_numpy()
        )
        Y = np.log1p(10000 * raw / depth[:, None])
        design = np.column_stack(
            [np.ones(2 * n), np.r_[np.ones(n), np.zeros(n)], np.log1p(depth), np.log1p(atac_depth)]
        )
        design[:, 2:] = (design[:, 2:] - design[:, 2:].mean(axis=0)) / design[:, 2:].std(axis=0)
        XR, YR = residualize(X, design), residualize(Y, design)
        xs = np.sqrt((XR * XR).sum(axis=0))
        ys = np.sqrt((YR * YR).sum(axis=0))
        pi = q.local_peak.map(peak_lookup).to_numpy()
        gi = q.local_gene.map(gene_lookup).to_numpy()
        corr = np.full(len(q), np.nan)
        for left in range(0, len(q), 500):
            sl = slice(left, min(left + 500, len(q)))
            num = np.einsum("ij,ij->j", XR[:, pi[sl]], YR[:, gi[sl]])
            den = xs[pi[sl]] * ys[gi[sl]]
            corr[sl] = np.divide(num, den, out=np.full(len(num), np.nan), where=den > 0)
        q["gsm"] = gsm
        q["n_pairs"] = n
        q["n_peak_open"] = X[:, pi].sum(axis=0).astype(int)
        q["n_gene_detected"] = (raw[:, gi] > 0).sum(axis=0).astype(int)
        q["partial_r"] = corr
        q.loc[(q.n_peak_open < 10) | (q.n_gene_detected < 20), "partial_r"] = np.nan
        candidate_rows.append(
            q[
                [
                    "gsm",
                    "peak",
                    "gene",
                    "nearest_tss_bp",
                    "n_pairs",
                    "n_peak_open",
                    "n_gene_detected",
                    "partial_r",
                ]
            ]
        )
        print(gsm, "tested", int(q.partial_r.notna().sum()), "of", len(q), flush=True)
    by_library = pd.concat(candidate_rows, ignore_index=True)
    by_library.to_csv(
        a.tables / f"peak_gene_links_by_library{suffix}.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    rows = []
    for (peak, gene), sub in by_library.groupby(["peak", "gene"]):
        s = sub.dropna(subset=["partial_r"])
        if len(s) < 3:
            continue
        z = np.arctanh(np.clip(s.partial_r.to_numpy(), -0.999, 0.999))
        w = (2 * s.n_pairs.to_numpy() - 5).astype(float)
        mean_z = np.average(z, weights=w)
        p_cell = 2 * stats.norm.sf(abs(mean_z * np.sqrt(w.sum())))
        rows.append(
            {
                "peak": peak,
                "gene": gene,
                "nearest_tss_bp": int(s.nearest_tss_bp.iloc[0]),
                "n_libraries": len(s),
                "positive_link_libraries": int((s.partial_r > 0).sum()),
                "n_peak_open": int(s.n_peak_open.sum()),
                "n_gene_detected": int(s.n_gene_detected.sum()),
                "partial_r": float(np.tanh(mean_z)),
                "p_link_cell": float(p_cell),
            }
        )
    links = pd.DataFrame(rows)
    links["q_all_links"] = multipletests(links.p_link_cell, method="fdr_bh")[1]
    if a.motif_scores:
        scores = pd.read_csv(a.motif_scores, sep="\t")
        cols = ["MYOD1_score", "MYOG_score", "MYF5_score", "MYF6_score"]
        scores["mrf_max_score"] = scores[cols].max(axis=1)
        links = links.merge(
            scores[["peak", "mrf_max_score"]], on="peak", how="left", validate="many_to_one"
        )
    links.sort_values("q_all_links").to_csv(
        a.tables / f"peak_gene_links{suffix}.tsv", sep="\t", index=False
    )
    if a.gate == "TSS_ge_2":
        print(
            "TSS>=2 links tested",
            len(links),
            "positive q<.05",
            int(((links.q_all_links < 0.05) & (links.partial_r > 0)).sum()),
            flush=True,
        )
        return
    print(
        "links tested",
        len(links),
        "positive q<.05",
        int(((links.q_all_links < 0.05) & (links.partial_r > 0)).sum()),
    )


if __name__ == "__main__":
    main()
