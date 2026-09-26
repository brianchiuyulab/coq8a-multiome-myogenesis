"""Annotate HSMM enhancer selectivity against eight external cell references."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import sparse
from statsmodels.stats.multitest import multipletests


def main(a):
    a.out.mkdir(parents=True, exist_ok=True)
    a.references.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location(
        "region_tests", a.root / "scripts/54_enhancer_region_tests.py"
    )
    region = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(region)
    ann = pd.read_csv(
        a.root / "results/external_enhancers/external_peak_annotations.tsv", sep="\t"
    )
    ann = ann[ann.strong_enhancer].copy()
    names = ["Gm12878", "H1hesc", "Hepg2", "Hmec", "Huvec", "K562", "Nhek", "Nhlf"]

    def reference(name):
        filename = f"wgEncodeBroadHmm{name}HMM.bed.gz"
        url = (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/"
            + filename
        )
        path = a.references / filename
        if not path.exists():
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            path.write_bytes(response.content)
        states = pd.read_csv(path, sep="\t", header=None).iloc[:, :4]
        states.columns = ["chrom", "start", "end", "state"]
        states = states[states.state.isin(["4_Strong_Enhancer", "5_Strong_Enhancer"])]
        refs = dict(tuple(states.groupby("chrom")))
        overlap = []
        for r in ann.itertuples():
            g = refs.get(r.chrom19)
            overlap.append(
                False
                if g is None
                else bool(((g.start < r.end19) & (g.end > r.start19)).any())
            )
        return (
            name,
            overlap,
            dict(
                name=name, url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest()
            ),
        )

    manifest = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for name, overlap, source in pool.map(reference, names):
            ann[name + "_strong"] = overlap
            manifest.append(source)
            print(name, sum(overlap), flush=True)
    ann["n_other_strong"] = ann[[n + "_strong" for n in names]].sum(1)
    ann.to_csv(a.out / "peak_selectivity.tsv", sep="\t", index=False)
    sets = {
        "HSMM_strong_all": ann.peak.tolist(),
        "HSMM_strong_absent_in_8_others": ann.loc[
            ann.n_other_strong.eq(0), "peak"
        ].tolist(),
        "HSMM_strong_at_most_1_other": ann.loc[
            ann.n_other_strong.le(1), "peak"
        ].tolist(),
    }
    pairs = pd.read_csv(a.root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[
        (pairs.gate == "TSS_ge_3") & (pairs.contrast == "3plus_vs_1")
    ].sort_values(["gsm", "pair"])
    peaks = (a.cache / "peaks.txt").read_text().splitlines()
    pi = {p: i for i, p in enumerate(peaks)}
    highs, lows = [], []
    for gsm, g in pairs.groupby("gsm"):
        meta = pd.read_csv(a.cache / f"{gsm}_meta.tsv", sep="\t").set_index("barcode")
        x = sparse.load_npz(a.cache / f"{gsm}_atac.npz")
        hi, lo = meta.index.get_indexer(g.high_barcode), meta.index.get_indexer(
            g.low_barcode
        )
        assert min(hi.min(), lo.min()) >= 0
        highs.append(x[:, hi].toarray().T)
        lows.append(x[:, lo].toarray().T)
    h, l = np.vstack(highs), np.vstack(lows)
    effects = pd.read_csv(a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t")
    original = pd.read_csv(
        a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    genes = (
        original.groupby("peak")
        .gene.agg(lambda g: ";".join(sorted(set(g))))
        .rename("nearby_genes")
    )
    rows, perpeak, bylib = [], [], []
    for name, members in sets.items():
        idx = [pi[p] for p in members]
        hc, lc = h[:, idx].sum(1).astype(int), l[:, idx].sum(1).astype(int)
        delta = hc - lc
        rows.append(
            dict(
                set_id=name,
                n_peaks=len(members),
                FC=hc.sum() / lc.sum(),
                p=region.paired_exact_score(delta),
                contains_MYOD1_focal="chr11:17649919-17650798" in members,
            )
        )
        for gsm in pairs.gsm.unique():
            mask = pairs.gsm.to_numpy() == gsm
            bylib.append(
                dict(
                    set_id=name,
                    gsm=gsm,
                    n_pairs=int(mask.sum()),
                    high=hc[mask].mean() / len(members),
                    low=lc[mask].mean() / len(members),
                    FC=hc[mask].sum() / lc[mask].sum(),
                )
            )
        e = (
            effects[effects.peak.isin(members)]
            .copy()
            .merge(genes, on="peak", how="left")
        )
        e["set_id"] = name
        e["q_within_set"] = multipletests(e.p_pair_binomial, method="fdr_bh")[1]
        perpeak.append(e)
    summary = pd.DataFrame(rows)
    summary["q_three_sets"] = multipletests(summary.p, method="fdr_bh")[1]
    summary.to_csv(a.out / "selectivity_set_results.tsv", sep="\t", index=False)
    pd.concat(perpeak).to_csv(
        a.out / "selectivity_peak_results.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame(bylib).to_csv(
        a.out / "selectivity_by_library.tsv", sep="\t", index=False
    )
    (a.out / "references.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(ann[ann.peak.eq("chr11:17649919-17650798")].to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--references", type=Path, required=True)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())
