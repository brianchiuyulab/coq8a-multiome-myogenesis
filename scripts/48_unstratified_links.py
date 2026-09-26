"""Screen all four-library myogenesis peaks without temporal selection."""

import argparse
import gzip
import importlib.util
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmwrite
from statsmodels.stats.multitest import multipletests


def prepare(a):
    spec = importlib.util.spec_from_file_location(
        "extract_middle", a.root / "scripts/44_expand_middle_links.py"
    )
    extraction = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(extraction)
    original = pd.read_csv(
        a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    peaks = sorted(original.loc[original.n_libraries.eq(4), "peak"].unique())
    assert len(peaks) == 5097
    genomic = {}
    for p in peaks:
        chrom, pos = p.split(":")
        center = np.mean(list(map(int, pos.split("-"))))
        genomic.setdefault(chrom, []).append((center, p))
    genomic = {c: sorted(v) for c, v in genomic.items()}
    centers = {c: np.array([x[0] for x in v]) for c, v in genomic.items()}
    first = next(a.h5.glob("GSM6339597_*filtered_feature_bc_matrix.h5"))
    with h5py.File(first) as h:
        features = h["matrix/features"]
        names = np.char.decode(features["name"][:])
        names = names[features["feature_type"][:] == b"Gene Expression"]
        frequency = pd.Series(names).value_counts()
        valid = set(frequency[frequency.eq(1)].index)
    pairs = {}
    with gzip.open(a.gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            z = line.rstrip().split("\t")
            if z[2] != "transcript" or z[0] not in genomic:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', z[8]))
            gene = attrs["gene_name"]
            if gene not in valid:
                continue
            tss = int(z[3] if z[6] == "+" else z[4]) - 1
            c = centers[z[0]]
            lo, hi = np.searchsorted(c, [tss - 500000, tss + 500000], side="left")
            hi = np.searchsorted(c, tss + 500000, side="right")
            for center, peak in genomic[z[0]][lo:hi]:
                key = (peak, gene)
                distance = abs(center - tss)
                if key not in pairs or distance < pairs[key][2]:
                    pairs[key] = (attrs["gene_id"], attrs["gene_type"], distance)
    candidates = pd.DataFrame(
        [(p, g, *v) for (p, g), v in pairs.items()],
        columns=["peak", "gene", "gene_id", "gene_type", "distance_bp"],
    )
    candidates = candidates.sort_values(["peak", "gene"])
    scoregenes = (a.middle / "programme_genes.txt").read_text().splitlines()
    genes = sorted(set(candidates.gene) | set(scoregenes) | {"COQ8A"})
    candidates.to_csv(a.out / "candidate_pairs.tsv", sep="\t", index=False)
    (a.out / "peaks.txt").write_text("\n".join(peaks) + "\n")
    (a.out / "genes.txt").write_text("\n".join(genes) + "\n")
    (a.out / "programme_genes.txt").write_text("\n".join(scoregenes) + "\n")
    mapping = pd.read_csv(
        a.root / "results/tables/consensus_peak_map.tsv.gz", sep="\t"
    ).set_index("peak")
    print(
        "Candidates",
        len(candidates),
        "peaks",
        len(peaks),
        "genes",
        len(genes),
        flush=True,
    )
    allr, allp, allmeta = [], [], []
    for file in sorted(a.middle.glob("GSM*_meta.tsv")):
        gsm = file.name.split("_")[0]
        meta = pd.read_csv(file, sep="\t")
        local = [
            p if gsm == "GSM6339597" else mapping.loc[p, gsm + "_peak"] for p in peaks
        ]
        matrix, mito = extraction.extract(
            next(a.h5.glob(gsm + "_*filtered_feature_bc_matrix.h5")),
            genes,
            local,
            meta.barcode.tolist(),
        )
        r = matrix[: len(genes)].astype(np.float32)
        p = matrix[len(genes) :].astype(np.float32)
        p.data[:] = 1
        assert np.array_equal(r[genes.index("COQ8A")].toarray().ravel(), meta.COQ8A_umi)
        assert np.allclose(mito / meta.total_rna_umi * 100, meta.percent_mito)
        sparse.save_npz(a.out / (gsm + "_rna.npz"), r)
        sparse.save_npz(a.out / (gsm + "_atac.npz"), p)
        meta.to_csv(a.out / (gsm + "_meta.tsv"), sep="\t", index=False)
        allr.append(r)
        allp.append(p)
        allmeta.append(meta)
        print(gsm, len(meta), "extracted", flush=True)
    mmwrite(a.out / "rna.mtx", sparse.hstack(allr).tocoo())
    mmwrite(a.out / "atac.mtx", sparse.hstack(allp).tocoo())
    pd.concat(allmeta).to_csv(a.out / "metadata.tsv", sep="\t", index=False)


def analyze(a):
    candidates = pd.read_csv(a.out / "candidate_pairs.tsv", sep="\t")
    genes = (a.out / "genes.txt").read_text().splitlines()
    peaks = (a.out / "peaks.txt").read_text().splitlines()
    scoregenes = set((a.out / "programme_genes.txt").read_text().splitlines())
    gi = {g: i for i, g in enumerate(genes)}
    pi = {p: i for i, p in enumerate(peaks)}
    outputs = []
    for file in sorted(a.out.glob("GSM*_meta.tsv")):
        gsm = file.name.split("_")[0]
        meta = pd.read_csv(file, sep="\t")
        raw = sparse.load_npz(a.out / (gsm + "_rna.npz"))
        x = sparse.load_npz(a.out / (gsm + "_atac.npz")).toarray().T.astype(np.float64)
        gene_detected = np.asarray((raw > 0).sum(axis=1)).ravel()
        peak_detected = x.sum(axis=0).astype(int)
        y = raw.toarray().T.astype(np.float64)
        y /= meta.total_rna_umi.to_numpy()[:, None]
        y *= 10000
        np.log1p(y, out=y)
        del raw
        ip = candidates.peak.map(pi).to_numpy()
        ig = candidates.gene.map(gi).to_numpy()
        eligible = (
            (peak_detected[ip] / len(meta) > 0.01)
            & (gene_detected[ig] / len(meta) > 0.01)
            & (gene_detected[ig] >= 20)
            & (peak_detected[ip] >= 10)
        )
        subset = candidates.loc[eligible].copy()
        ip = ip[eligible]
        ig = ig[eligible]
        in_score = subset.gene.isin(scoregenes).to_numpy(dtype=float)
        ss = meta.myogenesis_score.to_numpy() * len(scoregenes)
        for coq in [False, True]:
            terms = [np.log1p(meta.total_rna_umi), np.log1p(meta.total_open_peaks)]
            if coq:
                terms.append(np.log1p(meta.COQ8A_umi))
            cov = np.column_stack(terms)
            cov = (cov - cov.mean(0)) / cov.std(0)
            cov = np.column_stack([np.ones(len(meta)), cov])
            # Projection sufficient statistics avoid retaining residual copies.
            basis = np.linalg.qr(cov, mode="reduced")[0]
            qx = basis.T @ x
            qy = basis.T @ y
            qs = basis.T @ ss
            vx = np.einsum("ij,ij->j", x, x) - np.einsum("ij,ij->j", qx, qx)
            vy = np.einsum("ij,ij->j", y, y) - np.einsum("ij,ij->j", qy, qy)
            xsum = x.T @ ss - qx.T @ qs
            ysum = y.T @ ss - qy.T @ qs
            vss = ss @ ss - qs @ qs
            base = np.full(len(subset), np.nan)
            adjusted = base.copy()
            for left in range(0, len(subset), 128):
                sl = slice(left, left + 128)
                px = ip[sl]
                gy = ig[sl]
                e = in_score[sl]
                xy = np.einsum("ij,ij->j", x[:, px], y[:, gy]) - np.einsum(
                    "ij,ij->j", qx[:, px], qy[:, gy]
                )
                base[sl] = xy / np.sqrt(vx[px] * vy[gy])
                xs = xsum[px] - e * xy
                ys = ysum[gy] - e * vy[gy]
                vs = vss - 2 * e * ysum[gy] + e * vy[gy]
                adjusted[sl] = (xy - xs * ys / vs) / np.sqrt(
                    np.maximum(vx[px] - xs**2 / vs, 0)
                    * np.maximum(vy[gy] - ys**2 / vs, 0)
                )
            for state, values in [(False, base), (True, adjusted)]:
                d = subset.copy()
                d["gsm"] = gsm
                d["n"] = len(meta)
                d["r"] = values
                d["model"] = (
                    "technical" + ("_state" if state else "") + ("_coq" if coq else "")
                )
                d["rank"] = cov.shape[1] + int(state)
                d["peak_detected"] = peak_detected[ip]
                d["gene_detected"] = gene_detected[ig]
                outputs.append(d)
        del x, y
        print(gsm, "evaluated", len(subset), "pairs", flush=True)
    per = pd.concat(outputs, ignore_index=True)
    per.to_csv(
        a.out / "correlations_by_library.tsv.gz", sep="\t", index=False, na_rep="NA"
    )
    rows = []
    for detection in [0.05, 0.01]:
        eligible = per[
            (per.peak_detected / per.n > detection)
            & (per.gene_detected / per.n > detection)
            & per.r.notna()
        ]
        for (model, peak, gene), g in eligible.groupby(
            ["model", "peak", "gene"], sort=False
        ):
            if len(g) < 3:
                continue
            z = np.arctanh(np.clip(g.r.to_numpy(), -0.999999, 0.999999))
            w = (g.n - g["rank"] - 2).to_numpy()
            mz = np.average(z, weights=w)
            stat = mz * np.sqrt(w.sum())
            loo = [
                np.tanh(np.average(np.delete(z, i), weights=np.delete(w, i)))
                for i in range(len(g))
            ]
            rows.append(
                dict(
                    model=model,
                    peak=peak,
                    gene=gene,
                    detection_fraction=detection,
                    n_libraries=len(g),
                    positive_libraries=int(sum(g.r > 0)),
                    r=np.tanh(mz),
                    p=2 * stats.norm.sf(abs(stat)),
                    loo_r_min=min(loo),
                    loo_r_max=max(loo),
                )
            )
    combined = pd.DataFrame(rows)
    combined["q"] = np.nan
    for _, g in combined.groupby(["model", "detection_fraction"]):
        combined.loc[g.index, "q"] = multipletests(g.p, method="fdr_bh")[1]
    combined.to_csv(
        a.out / "correlations_combined.tsv.gz", sep="\t", index=False, na_rep="NA"
    )
    # Every eligible original middle link must agree because its target was
    # already excluded from the fixed programme score.
    old = pd.read_csv(a.middle / "correlations_combined.tsv", sep="\t")
    old = old[(old.population == "all_qc") & old.r.notna()]
    check = combined.merge(
        old,
        on=["peak", "gene", "model", "detection_fraction"],
        suffixes=("_new", "_old"),
    )
    assert len(check) == len(old)
    assert np.allclose(check.r_new, check.r_old, atol=1e-7)
    effects = pd.read_csv(
        a.root / "results/tables/candidate_peak_effects_pooled.tsv.gz", sep="\t"
    )
    effects = effects[
        (effects.gate == "TSS_ge_3")
        & (effects.contrast == "3plus_vs_1")
        & effects.peak.isin(peaks)
    ].copy()
    assert (
        len(effects) == 5097
        and effects.n_libraries.eq(4).all()
        and effects.n_pairs.eq(201).all()
    )
    effects["fold_open"] = effects.high_open / effects.low_open.replace(0, np.nan)
    effects["q_all_5097"] = multipletests(effects.p_pair_binomial, method="fdr_bh")[1]
    effects.to_csv(a.out / "ATAC_all_5097.tsv", sep="\t", index=False, na_rep="NA")
    validation = {
        "peaks": len(peaks),
        "candidate_pairs": len(candidates),
        "middle_correlations_reproduced": len(check),
        "all_qc_nuclei": 32977,
        "detection_families": combined.groupby(["model", "detection_fraction"])
        .size()
        .to_dict(),
    }
    validation["detection_families"] = {
        str(k): v for k, v in validation["detection_families"].items()
    }
    (a.out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--h5", type=Path, required=True)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--middle", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reuse", action="store_true")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if not a.reuse:
        prepare(a)
    analyze(a)
