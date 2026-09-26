"""Extract all-QC nuclei and reassess middle-region cis associations by library."""

import argparse
import gzip
import json
import re
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmwrite
from statsmodels.stats.multitest import multipletests


def extract(path, genes, peaks, cells):
    with h5py.File(path) as h:
        m = h["matrix"]
        names = np.char.decode(m["features"]["name"][:])
        bars = np.char.decode(m["barcodes"][:])
        index = {s: i for i, s in enumerate(names)}
        wanted = [index[g] for g in genes] + [index[p] for p in peaks]
        lookup = np.full(len(names), -1)
        lookup[wanted] = np.arange(len(wanted))
        mito = np.char.startswith(names, "MT-")
        ptr = m["indptr"][:]
        blocks = []
        mitos = []
        for start in range(0, len(bars), 256):
            end = min(start + 256, len(bars))
            lo, hi = int(ptr[start]), int(ptr[end])
            ids = m["indices"][lo:hi]
            values = m["data"][lo:hi]
            cols = np.repeat(np.arange(end - start), np.diff(ptr[start : end + 1]))
            mitos.extend(
                np.bincount(
                    cols[mito[ids]], weights=values[mito[ids]], minlength=end - start
                )
            )
            keep = lookup[ids] >= 0
            blocks.append(
                sparse.coo_matrix(
                    (values[keep], (lookup[ids[keep]], cols[keep])),
                    shape=(len(wanted), end - start),
                ).tocsr()
            )
        matrix = sparse.hstack(blocks, format="csr")
        ix = pd.Index(bars).get_indexer(cells)
        assert (ix >= 0).all()
        return matrix[:, ix], np.asarray(mitos)[ix]


def prepare(a):
    tables = a.root / "results/tables"
    temporal = a.root / "results/temporal"
    members = pd.read_csv(temporal / "temporal_set_memberships.tsv.gz", sep="\t")
    peaks = sorted(
        members.loc[members.set_id.eq("P2_O50_Middle_24to48h"), "peak"].unique()
    )
    assert len(peaks) == 34
    positions = {
        p: (p.split(":")[0], np.mean(list(map(int, p.split(":")[1].split("-")))))
        for p in peaks
    }
    rows = []
    with gzip.open(a.gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            z = line.rstrip().split("\t")
            if z[2] != "transcript":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', z[8]))
            tss = int(z[3] if z[6] == "+" else z[4]) - 1
            for p, (chrom, center) in positions.items():
                if z[0] == chrom and abs(tss - center) <= 500000:
                    rows.append(
                        (
                            p,
                            attrs["gene_name"],
                            attrs["gene_id"],
                            attrs["gene_type"],
                            abs(tss - center),
                        )
                    )
    candidates = (
        pd.DataFrame(
            rows, columns=["peak", "gene", "gene_id", "gene_type", "distance_bp"]
        )
        .groupby(["peak", "gene", "gene_id", "gene_type"], as_index=False)
        .distance_bp.min()
    )
    first = next(a.h5.glob("GSM6339597_*filtered_feature_bc_matrix.h5"))
    with h5py.File(first) as h:
        f = h["matrix/features"]
        names = np.char.decode(f["name"][:])
        types = f["feature_type"][:]
        rna_names = names[types == b"Gene Expression"]
        counts = pd.Series(rna_names).value_counts()
        ambiguous = set(counts[counts > 1].index)
        pd.Series(sorted(ambiguous)).to_csv(
            a.out / "ambiguous_rna_symbols_excluded.txt", index=False, header=False
        )
        rna_names = np.array([g for g in rna_names if g not in ambiguous])
    candidates = candidates[candidates.gene.isin(rna_names)].drop_duplicates(
        ["peak", "gene"]
    )
    myo = pd.read_csv(
        a.root / "reference/myogenesis_221_gene_sources.tsv", sep="\t"
    ).gene.tolist()
    scoregenes = sorted(set(myo) - set(candidates.gene) - {"COQ8A"})
    genes = sorted(
        set(candidates.gene)
        | set(scoregenes)
        | {"COQ8A", "PAX7", "MYOD1", "MYOG", "MYH3", "PTPRC", "PECAM1", "PDGFRA"}
    )
    genes = [g for g in genes if g in rna_names]
    candidates.to_csv(a.out / "candidate_pairs.tsv", sep="\t", index=False)
    pd.Series(scoregenes).to_csv(
        a.out / "programme_genes.txt", index=False, header=False
    )
    pd.Series(genes).to_csv(a.out / "genes.txt", index=False, header=False)
    pd.Series(peaks).to_csv(a.out / "peaks.txt", index=False, header=False)
    mapping = pd.read_csv(tables / "consensus_peak_map.tsv.gz", sep="\t").set_index(
        "peak"
    )
    qc = pd.read_csv(tables / "qc_nuclei.tsv.gz", sep="\t")
    qc = qc[qc.pass_tss3].copy()
    allmeta = []
    allr = []
    allp = []
    for gsm, meta in qc.groupby("gsm", sort=True):
        meta = meta.sort_values("barcode").copy()
        local = [
            p if gsm == "GSM6339597" else mapping.loc[p, gsm + "_peak"] for p in peaks
        ]
        matrix, mito = extract(
            next(a.h5.glob(gsm + "_*filtered_feature_bc_matrix.h5")),
            genes,
            local,
            meta.barcode.tolist(),
        )
        r = matrix[: len(genes)].astype(float)
        p = matrix[len(genes) :].astype(float)
        p.data[:] = 1
        gi = {g: i for i, g in enumerate(genes)}
        assert np.array_equal(
            r[gi["COQ8A"]].toarray().ravel(), meta.COQ8A_umi.to_numpy()
        )
        score = r[[gi[g] for g in scoregenes]].toarray()
        score = np.log1p(score / meta.total_rna_umi.to_numpy()[None, :] * 10000).mean(
            axis=0
        )
        meta["myogenesis_score"] = score
        meta["percent_mito"] = mito / meta.total_rna_umi.to_numpy() * 100
        meta["cell"] = gsm + "|" + meta.barcode
        allmeta.append(meta)
        allr.append(r)
        allp.append(p)
        sparse.save_npz(a.out / (gsm + "_rna.npz"), r)
        sparse.save_npz(a.out / (gsm + "_atac.npz"), p)
        meta.to_csv(a.out / (gsm + "_meta.tsv"), sep="\t", index=False)
        print(gsm, len(meta), "nuclei extracted", flush=True)
    metadata = pd.concat(allmeta, ignore_index=True)
    metadata.to_csv(a.out / "metadata.tsv", sep="\t", index=False)
    mmwrite(a.out / "rna.mtx", sparse.hstack(allr).tocoo())
    mmwrite(a.out / "atac.mtx", sparse.hstack(allp).tocoo())
    print(
        "Candidate pairs",
        len(candidates),
        "genes",
        len(genes),
        "independent score genes",
        len(scoregenes),
        flush=True,
    )


def analyze(a):
    cand = pd.read_csv(a.out / "candidate_pairs.tsv", sep="\t")
    genes = (a.out / "genes.txt").read_text().splitlines()
    peaks = (a.out / "peaks.txt").read_text().splitlines()
    gi = {g: i for i, g in enumerate(genes)}
    pi = {p: i for i, p in enumerate(peaks)}
    ix = cand.peak.map(pi).to_numpy()
    iy = cand.gene.map(gi).to_numpy()
    outputs = []
    pairs = pd.read_csv(a.root / "results/tables/matched_pairs.tsv.gz", sep="\t").query(
        "gate=='TSS_ge_3' and contrast=='3plus_vs_1'"
    )
    for file in sorted(a.out.glob("GSM*_meta.tsv")):
        gsm = file.name.split("_")[0]
        meta = pd.read_csv(file, sep="\t")
        raw = sparse.load_npz(a.out / (gsm + "_rna.npz")).toarray().T
        x = sparse.load_npz(a.out / (gsm + "_atac.npz")).toarray().T
        y = np.log1p(raw / meta.total_rna_umi.to_numpy()[:, None] * 10000)
        for population in ["all_qc", "matched"]:
            ids = (
                np.arange(len(meta))
                if population == "all_qc"
                else np.flatnonzero(
                    meta.barcode.isin(
                        set(pairs.loc[pairs.gsm.eq(gsm), "high_barcode"])
                        | set(pairs.loc[pairs.gsm.eq(gsm), "low_barcode"])
                    )
                )
            )
            m = meta.iloc[ids]
            xx = x[ids]
            yy = y[ids]
            rr = raw[ids]
            for model in [
                "technical",
                "technical_coq",
                "technical_state",
                "technical_state_coq",
            ]:
                cov = [np.log1p(m.total_rna_umi), np.log1p(m.total_open_peaks)]
                if "state" in model:
                    cov.append(m.myogenesis_score.to_numpy())
                if "coq" in model:
                    cov.append(
                        np.log1p(m.COQ8A_umi).to_numpy()
                        if population == "all_qc"
                        else (m.COQ8A_umi >= 3).to_numpy()
                    )
                cov = np.column_stack(cov)
                cov = (cov - cov.mean(0)) / np.maximum(cov.std(0), 1e-10)
                cov = np.column_stack([np.ones(len(ids)), cov])
                xr = xx - cov @ np.linalg.lstsq(cov, xx, rcond=None)[0]
                yr = yy - cov @ np.linalg.lstsq(cov, yy, rcond=None)[0]
                numerator = np.einsum("ij,ij->j", xr[:, ix], yr[:, iy])
                den = np.linalg.norm(xr[:, ix], axis=0) * np.linalg.norm(
                    yr[:, iy], axis=0
                )
                r = np.divide(
                    numerator, den, out=np.full(len(cand), np.nan), where=den > 1e-12
                )
                d = cand.copy()
                d["gsm"] = gsm
                d["population"] = population
                d["model"] = model
                d["n"] = len(ids)
                d["covariate_rank"] = np.linalg.matrix_rank(cov)
                d["r"] = r
                d["peak_detected"] = xx[:, ix].sum(0).astype(int)
                d["gene_detected"] = (rr[:, iy] > 0).sum(0)
                outputs.append(d)
        print(gsm, "correlations completed", flush=True)
    per = pd.concat(outputs, ignore_index=True)
    per.to_csv(a.out / "correlations_by_library.tsv.gz", sep="\t", index=False)
    rows = []
    for (pop, model, peak, gene), group in per.groupby(
        ["population", "model", "peak", "gene"], sort=False
    ):
        for detection in ([0.01, 0.05] if pop == "all_qc" else [0]):
            good = group[
                (group.peak_detected >= 10)
                & (group.gene_detected >= 20)
                & (group.peak_detected / group.n > detection)
                & (group.gene_detected / group.n > detection)
                & group.r.notna()
            ]
            rec = dict(
                population=pop,
                model=model,
                peak=peak,
                gene=gene,
                detection_fraction=detection,
                n_libraries=len(good),
                positive_libraries=int((good.r > 0).sum()),
                n_nuclei=int(good.n.sum()),
                gene_type=group.gene_type.iloc[0],
            )
            if len(good) >= 3:
                z = np.arctanh(np.clip(good.r, -0.999999, 0.999999))
                w = (
                    good.n - 5
                    if pop == "matched" and model == "technical_coq"
                    else good.n - good.covariate_rank - 2
                )
                mz = np.average(z, weights=w)
                se = 1 / np.sqrt(w.sum())
                rec.update(r=np.tanh(mz), p=2 * stats.norm.sf(abs(mz / se)))
                loo = [
                    np.tanh(
                        np.average(
                            np.delete(z.to_numpy(), i),
                            weights=np.delete(w.to_numpy(), i),
                        )
                    )
                    for i in range(len(good))
                ]
                rec["loo_r_min"] = min(loo)
                rec["loo_r_max"] = max(loo)
            rows.append(rec)
    result = pd.DataFrame(rows)
    result["q"] = np.nan
    for _, g in result.groupby(["population", "model", "detection_fraction"]):
        ok = g.p.notna()
        result.loc[g.index[ok], "q"] = multipletests(g.loc[ok, "p"], method="fdr_bh")[1]
    result.to_csv(a.out / "correlations_combined.tsv", sep="\t", index=False, na_rep="NA")
    old = pd.read_csv(
        a.root / "results/temporal/dynamic_cis_links.tsv", sep="\t"
    ).query("gate=='TSS_ge_3' and contrast=='3plus_vs_1'")
    check = result.query("population=='matched' and model=='technical_coq'").merge(
        old, on=["peak", "gene"]
    )
    mask = check.partial_r.notna()
    assert np.allclose(check.loc[mask, "r"], check.loc[mask, "partial_r"], atol=1e-7)
    (a.out / "validation.json").write_text(
        json.dumps(
            {
                "matched_correlations_reproduced": int(mask.sum()),
                "candidate_pairs": len(cand),
                "all_qc_nuclei": int(
                    pd.read_csv(a.out / "metadata.tsv", sep="\t").shape[0]
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--h5", type=Path, required=True)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reuse", action="store_true")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if not a.reuse:
        prepare(a)
    analyze(a)
