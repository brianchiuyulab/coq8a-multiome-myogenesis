"""Label-blind peak matching and knowledge-based myogenesis target space."""

import argparse
import gzip
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SAMPLES = ["GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603"]
CHR = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def coordinates(names):
    x = pd.Series(names).str.extract(r"^([^:]+):(\d+)-(\d+)$")
    assert x.notna().all().all()
    return x[0].to_numpy(), x[1].astype(int).to_numpy(), x[2].astype(int).to_numpy()


def read_peaks(root, gsm):
    path = next(root.glob(f"{gsm}_*filtered_feature_bc_matrix.h5"))
    with h5py.File(path) as h:
        f = h["matrix/features"]
        is_peak = f["feature_type"][:] == b"Peaks"
        return np.char.decode(f["name"][:][is_peak], "utf-8")


def map_one(ref, target):
    rc, rs, re_ = coordinates(ref)
    tc, ts, te = coordinates(target)
    match = np.full(len(ref), -1, dtype=np.int32)
    match_score = np.zeros(len(ref), dtype=np.float32)
    for chrom in sorted(set(rc) & set(tc)):
        r = np.flatnonzero(rc == chrom)
        t = np.flatnonzero(tc == chrom)
        mid = (ts[t] + te[t]) / 2
        order = np.argsort(mid)
        t, mid = t[order], mid[order]
        pos = np.searchsorted(mid, (rs[r] + re_[r]) / 2)
        for offset in [-3, -2, -1, 0, 1, 2]:
            ix = np.clip(pos + offset, 0, len(t) - 1)
            cand = t[ix]
            overlap = np.maximum(0, np.minimum(re_[r], te[cand]) - np.maximum(rs[r], ts[cand]))
            reciprocal = np.minimum(overlap / (re_[r] - rs[r]), overlap / (te[cand] - ts[cand]))
            better = reciprocal > match_score[r]
            match[r[better]] = cand[better]
            match_score[r[better]] = reciprocal[better]
    match[match_score < 0.5] = -1
    # A target peak cannot serve two anchor peaks in one library.
    mapped = np.flatnonzero(match >= 0)
    candidate = pd.DataFrame({"ref": mapped, "target": match[mapped], "score": match_score[mapped]})
    candidate = candidate.sort_values(["score", "ref"], ascending=[False, True]).drop_duplicates(
        "target"
    )
    keep = np.full(len(ref), -1, dtype=np.int32)
    keep[candidate.ref.to_numpy()] = candidate.target.to_numpy()
    return keep


def blacklist_flags(chrom, start, end, bed):
    blocked = np.zeros(len(chrom), dtype=bool)
    with gzip.open(bed, "rt") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            c, lo, hi = fields[0], int(fields[1]), int(fields[2])
            idx = np.flatnonzero(chrom == c)
            blocked[idx] |= (start[idx] < hi) & (end[idx] > lo)
    return blocked


def tss_positions(gtf, genes):
    positions = {}
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript" or fields[0] not in CHR:
                continue
            attr = fields[8]
            if 'gene_type "protein_coding"' not in attr:
                continue
            m = re.search(r'gene_name "([^"]+)"', attr)
            if not m or m.group(1) not in genes:
                continue
            tss = int(fields[3] if fields[6] == "+" else fields[4])
            positions.setdefault(fields[0], set()).add((tss, m.group(1)))
    return {c: sorted(v) for c, v in positions.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5-root", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    peaks = {gsm: read_peaks(a.h5_root, gsm) for gsm in SAMPLES}
    anchor = peaks[SAMPLES[0]]
    mapping = pd.DataFrame({"peak": anchor})
    for gsm in SAMPLES[1:]:
        match = map_one(anchor, peaks[gsm])
        mapping[gsm + "_peak"] = np.where(match >= 0, peaks[gsm][np.maximum(match, 0)], "")
        print(gsm, "reciprocal one-to-one matches", int((match >= 0).sum()), flush=True)
    mapping["n_libraries"] = 1 + (mapping[[g + "_peak" for g in SAMPLES[1:]]] != "").sum(axis=1)
    chrom, start, end = coordinates(anchor)
    mapping["chrom"] = chrom
    mapping["start"] = start
    mapping["end"] = end
    mapping["width"] = end - start
    mapping["blacklisted"] = blacklist_flags(chrom, start, end, a.blacklist)
    mapping["coordinate_qc"] = (
        mapping.chrom.isin(CHR) & mapping.width.between(200, 2000) & ~mapping.blacklisted
    )
    mapping.to_csv(
        a.out / "consensus_peak_map.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    genes = pd.read_csv(a.gene_list, sep="\t")
    assert len(genes) == 221 and genes.gene.nunique() == 221
    tss = tss_positions(a.gtf, set(genes.gene))
    rows = []
    for p in mapping.itertuples():
        if not p.coordinate_qc or p.n_libraries < 3 or p.chrom not in tss:
            continue
        midpoint = (p.start + p.end) // 2
        positions = tss[p.chrom]
        for site, gene in positions:
            distance = abs(site - midpoint)
            if distance <= 100_000:
                rows.append((p.peak, gene, distance, p.n_libraries))
    candidates = pd.DataFrame(rows, columns=["peak", "gene", "nearest_tss_bp", "n_libraries"])
    candidates = candidates.sort_values("nearest_tss_bp").drop_duplicates(["peak", "gene"])
    candidates.to_csv(
        a.out / "candidate_peak_gene.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    counts = (
        candidates.groupby("gene")
        .agg(
            n_candidate_peaks=("peak", "nunique"),
            n_common4=("n_libraries", lambda z: int((z == 4).sum())),
        )
        .reset_index()
    )
    genes.merge(counts, on="gene", how="left").fillna(
        {"n_candidate_peaks": 0, "n_common4": 0}
    ).to_csv(a.out / "gene_search_space.tsv", sep="\t", index=False)
    print(
        "candidate peak-gene pairs",
        len(candidates),
        "genes mapped",
        candidates.gene.nunique(),
        flush=True,
    )


if __name__ == "__main__":
    main()
