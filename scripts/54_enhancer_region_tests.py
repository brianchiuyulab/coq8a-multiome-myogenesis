"""Test external enhancer blocks and gene neighborhoods without a TF filter."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


def merge_states(path):
    states = pd.read_csv(path, sep="\t", header=None).iloc[:, :4]
    states.columns = ["chrom", "start", "end", "state"]
    states = states[states.state.isin(["4_Strong_Enhancer", "5_Strong_Enhancer"])]
    blocks = []
    for chrom, group in states.groupby("chrom"):
        local = []
        for row in group.sort_values("start").itertuples():
            if local and row.start <= local[-1][1]:
                local[-1][1] = max(row.end, local[-1][1])
            else:
                local.append([row.start, row.end])
        blocks.extend((chrom, s, e, f"{chrom}:{s}-{e}") for s, e in local)
    return pd.DataFrame(blocks, columns=["chrom", "start", "end", "block"])


def paired_exact_score(delta):
    """Exact two-sided sign-flip distribution conditional on pair magnitudes."""
    delta = np.asarray(delta, dtype=int)
    weights = np.abs(delta[delta != 0])
    distribution = np.ones(1)
    for weight in weights:
        new = np.zeros(len(distribution) + weight)
        new[: len(distribution)] += 0.5 * distribution
        new[weight:] += 0.5 * distribution
        distribution = new
    signed_sum = 2 * np.arange(len(distribution)) - weights.sum()
    assert np.isclose(distribution.sum(), 1)
    return min(1.0, float(distribution[np.abs(signed_sum) >= abs(delta.sum())].sum()))


def main(a):
    a.out.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(
        a.root / "results/external_enhancers/external_peak_annotations.tsv", sep="\t"
    )
    ann = ann[ann.strong_enhancer].copy()
    assert len(ann) == 1777
    if a.selectivity is not None:
        selective = pd.read_csv(a.selectivity, sep="\t")
        ann = ann[ann.peak.isin(selective.loc[selective.n_other_strong.eq(0), "peak"])]
        assert len(ann) == 401
    blocks = merge_states(a.hmm)
    chroms = dict(tuple(blocks.groupby("chrom")))
    membership = []
    for r in ann.itertuples():
        g = chroms[r.chrom19]
        g = g[(g.start < r.end19) & (g.end > r.start19)].copy()
        g["overlap_bp"] = np.minimum(g.end, r.end19) - np.maximum(g.start, r.start19)
        winner = g.sort_values(["overlap_bp", "start"], ascending=[False, True]).iloc[0]
        membership.append(
            dict(
                peak=r.peak,
                block=winner.block,
                overlap_bp=winner.overlap_bp,
                overlapping_blocks=len(g),
            )
        )
    membership = pd.DataFrame(membership)
    membership.to_csv(a.out / "block_membership.tsv", sep="\t", index=False)
    blocks[blocks.block.isin(membership.block)].to_csv(
        a.out / "external_blocks_hg19.tsv", sep="\t", index=False
    )
    candidates = pd.read_csv(
        a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    gene_members = candidates[candidates.peak.isin(ann.peak)][
        ["gene", "peak"]
    ].drop_duplicates()
    gene_members.to_csv(
        a.out / "gene_neighborhood_membership.tsv", sep="\t", index=False
    )
    sets = [
        ("enhancer_block", name, sorted(set(g.peak)))
        for name, g in membership.groupby("block")
    ]
    sets += [
        ("gene_neighborhood", name, sorted(set(g.peak)))
        for name, g in gene_members.groupby("gene")
    ]
    peaks = (a.cache / "peaks.txt").read_text().splitlines()
    idx = {p: i for i, p in enumerate(peaks)}
    pairs = pd.read_csv(a.root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[
        (pairs.gate == "TSS_ge_3") & (pairs.contrast == "3plus_vs_1")
    ].sort_values(["gsm", "pair"])
    assert len(pairs) == 201
    for gsm, g in pairs.groupby("gsm"):
        assert not set(g.high_barcode) & set(g.low_barcode)
        assert g.high_barcode.is_unique and g.low_barcode.is_unique
    highs, lows = [], []
    for gsm, g in pairs.groupby("gsm", sort=True):
        meta = pd.read_csv(a.cache / f"{gsm}_meta.tsv", sep="\t").set_index("barcode")
        x = sparse.load_npz(a.cache / f"{gsm}_atac.npz")
        hi, lo = meta.index.get_indexer(g.high_barcode), meta.index.get_indexer(
            g.low_barcode
        )
        assert min(hi.min(), lo.min()) >= 0
        highs.append(x[:, hi].toarray().T)
        lows.append(x[:, lo].toarray().T)
    high, low = np.vstack(highs), np.vstack(lows)
    assert set(np.unique(high)) <= {0, 1}
    effects = pd.read_csv(
        a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t"
    ).set_index("peak")
    assert np.array_equal(high.sum(0), effects.loc[peaks, "high_open"])
    assert np.array_equal(low.sum(0), effects.loc[peaks, "low_open"])
    focal = "chr11:17649919-17650798"
    librows = []
    for gsm in pairs.gsm.unique():
        mask = pairs.gsm.to_numpy() == gsm
        h, l = high[mask, idx[focal]], low[mask, idx[focal]]
        plus, minus = int(((h == 1) & (l == 0)).sum()), int(((h == 0) & (l == 1)).sum())
        librows.append(
            dict(
                gsm=gsm,
                source=pairs.loc[pairs.gsm.eq(gsm), "source"].iloc[0],
                stage=pairs.loc[pairs.gsm.eq(gsm), "stage"].iloc[0],
                n_pairs=int(mask.sum()),
                high_open=int(h.sum()),
                low_open=int(l.sum()),
                high_only=plus,
                low_only=minus,
                FC=h.sum() / l.sum() if l.sum() else np.nan,
                p=stats.binomtest(plus, plus + minus).pvalue if plus + minus else 1,
            )
        )
    pd.DataFrame(librows).to_csv(a.out / "MYOD1_by_library.tsv", sep="\t", index=False)
    summary, libeffects, sourceeffects = [], [], []
    for begin in range(0, len(sets), 128):
        subset = sets[begin : begin + 128]
        h = np.column_stack(
            [high[:, [idx[p] for p in members]].mean(1) for _, _, members in subset]
        ).astype(float)
        l = np.column_stack(
            [low[:, [idx[p] for p in members]].mean(1) for _, _, members in subset]
        ).astype(float)
        delta = h - l
        for j, (kind, name, members) in enumerate(subset):
            pvals = np.sort(effects.loc[members, "p_pair_binomial"].to_numpy())
            simes = min(
                1.0, float(np.min(pvals * len(pvals) / np.arange(1, len(pvals) + 1)))
            )
            indices = [idx[p] for p in members]
            count_delta = high[:, indices].sum(1).astype(int) - low[:, indices].sum(
                1
            ).astype(int)
            pp = paired_exact_score(count_delta)
            if len(members) == 1:
                assert np.isclose(pp, pvals[0])
            for gsm in pairs.gsm.unique():
                mask = pairs.gsm.to_numpy() == gsm
                libeffects.append(
                    dict(
                        kind=kind,
                        region=name,
                        gsm=gsm,
                        n_pairs=int(mask.sum()),
                        high=float(h[mask, j].mean()),
                        low=float(l[mask, j].mean()),
                        delta_pp=float(100 * delta[mask, j].mean()),
                    )
                )
            source_deltas = []
            for source in pairs.source.unique():
                gsms = pairs.loc[pairs.source.eq(source), "gsm"].unique()
                sh = np.mean([h[pairs.gsm.to_numpy() == gsm, j].mean() for gsm in gsms])
                sl = np.mean([l[pairs.gsm.to_numpy() == gsm, j].mean() for gsm in gsms])
                source_deltas.append(sh - sl)
                sourceeffects.append(
                    dict(
                        kind=kind,
                        region=name,
                        source=source,
                        high=sh,
                        low=sl,
                        FC=sh / sl if sl else np.nan,
                        delta_pp=100 * (sh - sl),
                    )
                )
            summary.append(
                dict(
                    kind=kind,
                    region=name,
                    n_peaks=len(members),
                    n_pairs=len(pairs),
                    high=float(h[:, j].mean()),
                    low=float(l[:, j].mean()),
                    FC=(
                        float(h[:, j].mean() / l[:, j].mean())
                        if l[:, j].mean()
                        else np.nan
                    ),
                    delta_pp=float(100 * delta[:, j].mean()),
                    p_mean_score=pp,
                    p_simes=simes,
                    positive_sources=int(np.sum(np.array(source_deltas) > 0)),
                    contains_MYOD1_focal=focal in members,
                    nearby_genes=";".join(
                        sorted(
                            set(
                                gene_members.loc[
                                    gene_members.peak.isin(members), "gene"
                                ]
                            )
                        )
                    ),
                )
            )
    summary = pd.DataFrame(summary)
    for kind, g in summary.groupby("kind"):
        for p in ["p_mean_score", "p_simes"]:
            summary.loc[g.index, "q_" + p[2:]] = multipletests(g[p], method="fdr_bh")[1]
    summary.sort_values(["kind", "p_mean_score"]).to_csv(
        a.out / "region_results.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame(libeffects).to_csv(
        a.out / "region_effects_by_library.tsv", sep="\t", index=False
    )
    pd.DataFrame(sourceeffects).to_csv(
        a.out / "region_effects_by_source.tsv", sep="\t", index=False, na_rep="NA"
    )
    (a.out / "analysis_parameters.json").write_text(
        json.dumps(
            dict(
                gate="TSS_ge_3",
                contrast="3plus_vs_1",
                n_pairs=201,
                n_sources=2,
                n_libraries=4,
                start_peaks=len(ann),
                selectivity=(
                    "No strong enhancer overlap in eight other references"
                    if a.selectivity
                    else "None"
                ),
                merge_rule="Overlapping or directly touching external states 4/5; zero gap",
                assignment="Largest overlap with external block; genomic-start tie break",
                test="Exact two-sided paired sign-flip of summed binary accessibility; dynamic programming",
                inference="Conditional matched-nucleus association; not donor-population inference",
                correction="BH separately for each region definition and endpoint; no best-method selection",
                gene_assignment="Original candidate neighborhoods, not established enhancer targets",
            ),
            indent=2,
        )
        + "\n"
    )
    print(pd.DataFrame(librows).to_string(index=False))
    print(summary.groupby("kind").size())
    print(summary[summary.contains_MYOD1_focal].to_string(index=False))
    print(summary.groupby("kind")[["q_mean_score", "q_simes"]].min())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--hmm", type=Path, required=True)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--selectivity", type=Path)
    main(p.parse_args())
