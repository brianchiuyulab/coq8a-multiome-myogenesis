"""Recompute the complete fixed Day-0 comparison without scanning parameters."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DESIGN = dict(
    config="D206",
    TSS=3,
    population="all_culture",
    contrast="2plus_vs_zero",
    caliper=0.30,
    matching="depth",
)


def load_grid():
    spec = importlib.util.spec_from_file_location(
        "day0", ROOT / "scripts/88_day0_sensitivity_grid.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(args):
    out = ROOT / "results/fixed_day0_D206"
    out.mkdir(exist_ok=True)
    universe = set(
        pd.read_csv(ROOT / "reference/myogenesis_221_gene_sources.tsv", sep="\t").gene
    )
    external = json.loads((ROOT / "reference/muscle_subprogrammes.json").read_text())
    genes = sorted(
        universe
        & (
            set(external["GOBP_MYOBLAST_DIFFERENTIATION"]["geneSymbols"])
            | set(external["GOBP_MYOBLAST_FUSION"]["geneSymbols"])
        )
    )
    candidates = pd.read_csv(
        ROOT / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    member = candidates[
        candidates.n_libraries.eq(4)
        & candidates.nearest_tss_bp.le(100000)
        & candidates.gene.isin(genes)
    ].copy()
    names = sorted(member.peak.unique())
    assert len(genes) == 25 and len(names) == 565
    member.to_csv(out / "peak_membership.tsv", sep="\t", index=False)
    index = pd.Index(names)
    groups = {gene: index.get_indexer(d.peak) for gene, d in member.groupby("gene")}
    grid = load_grid()
    if not args.reuse_cache:
        grid.prepare(ROOT, args.h5, args.work, names, genes)
    hs, ls, counts, pair_rows = [], [], [], []
    for gsm in grid.SAMPLES:
        meta = pd.read_csv(args.work / f"{gsm}_meta.tsv", sep="\t")
        atac = sparse.load_npz(args.work / f"{gsm}_atac.npz").toarray().T
        pairs, _ = grid.match(meta, DESIGN)
        hi, lo = np.asarray(pairs).T
        h, l = atac[hi], atac[lo]
        hs.append(h)
        ls.append(l)
        counts.append((h.sum(0), l.sum(0), len(pairs)))
        pair_rows.extend(
            dict(
                config="D206",
                gsm=gsm,
                high_barcode=meta.iloc[i].barcode,
                low_barcode=meta.iloc[j].barcode,
            )
            for i, j in pairs
        )
    h, l = np.vstack(hs), np.vstack(ls)
    assert len(h) == 527
    high, low = h.sum(0), l.sum(0)
    high_only, low_only = ((h > 0) & (l == 0)).sum(0), ((l > 0) & (h == 0)).sum(0)
    p = np.array(
        [
            stats.binomtest(int(x), int(x + y)).pvalue if x + y else 1
            for x, y in zip(high_only, low_only)
        ]
    )
    q = multipletests(p, method="fdr_bh")[1]
    fc = np.divide(high, low, out=np.full(len(names), np.nan), where=low > 0)
    up = (counts[0][0] > counts[0][1]) & (counts[1][0] > counts[1][1])
    down = (counts[0][0] < counts[0][1]) & (counts[1][0] < counts[1][1])
    peaks = pd.DataFrame(
        dict(
            config="D206",
            peak=names,
            n_pairs=len(h),
            high=high.astype(int),
            low=low.astype(int),
            FC=fc,
            p=p,
            q565=q,
            high_only=high_only,
            low_only=low_only,
            both_up=up,
            both_down=down,
        )
    )
    for i, (hc, lc, number) in enumerate(counts, 1):
        peaks[f"source{i}_high"] = hc.astype(int)
        peaks[f"source{i}_low"] = lc.astype(int)
        peaks[f"source{i}_pairs"] = number
    print(
        "Recomputing 25 regions in both directions: 2,000,000 permutations", flush=True
    )
    regional = grid.permutation(h - l, groups, genes, 2000000, 20261134)
    for row in regional:
        j = row.pop("top_index")
        row.update(DESIGN)
        row.update(
            n_pairs=len(h),
            top_peak=names[j],
            top_peak_FC=fc[j],
            top_peak_p=p[j],
            top_peak_q565=q[j],
            both_sources_same_direction=bool(
                up[j] if row["mode"] == "opening" else down[j]
            ),
            n_peaks=len(groups[row["gene"]]),
        )
    grid.save(peaks, out / "all_peaks.tsv")
    grid.save(pd.DataFrame(regional), out / "all_regions.tsv")
    grid.save(pd.DataFrame(pair_rows), out / "matched_pairs.tsv")
    print("Fixed-scope statistics complete", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--reuse-cache", action="store_true")
    main(parser.parse_args())
