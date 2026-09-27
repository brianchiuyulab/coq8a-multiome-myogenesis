"""Test sparse local accessibility changes with paired max-statistic permutations.

The unit of hypothesis is one frozen gene neighborhood. Each permutation swaps
COQ8A-high/low labels within matched nucleus pairs and preserves correlations
between peaks. This is a nucleus-association test, not a donor-level test.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


def main(root, work, permutations, seed):
    out = root / "results/local_signal_decision"
    out.mkdir(parents=True, exist_ok=True)
    membership = pd.read_csv(
        root / "results/differentiation_fusion25/GSE208248_peak_membership.tsv",
        sep="\t",
    )
    names = sorted(membership.peak.unique())
    genes = sorted(membership.gene.unique())
    group = {
        g: pd.Index(names).get_indexer(d.peak) for g, d in membership.groupby("gene")
    }
    allpeaks = (work / "peaks.txt").read_text().splitlines()
    pi = pd.Index(allpeaks).get_indexer(names)
    pairs = pd.read_csv(root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[pairs.gate.eq("TSS_ge_3")]
    tables = {}
    for gsm in sorted(pairs.gsm.unique()):
        meta = pd.read_csv(work / f"{gsm}_meta.tsv", sep="\t").set_index("barcode")
        matrix = sparse.load_npz(work / f"{gsm}_atac.npz")[pi]
        tables[gsm] = (meta, matrix)
    rng = np.random.default_rng(seed)
    rows = []
    with threadpool_limits(limits=2):
        for contrast, ps in pairs.groupby("contrast"):
            diffs = []
            for gsm, d in ps.groupby("gsm"):
                meta, matrix = tables[gsm]
                hi, lo = meta.index.get_indexer(d.high_barcode), meta.index.get_indexer(
                    d.low_barcode
                )
                assert min(hi) >= 0 and min(lo) >= 0
                diffs.append(
                    (matrix[:, hi] - matrix[:, lo]).toarray().T.astype(np.float32)
                )
            delta = np.vstack(diffs)
            discordance = np.square(delta).sum(0)
            scale = np.sqrt(np.maximum(discordance, 1))
            standardized = delta / scale
            observed = standardized.sum(0)
            modes = {"opening": lambda x: x, "two_sided": np.abs}
            counters = {mode: np.zeros(len(genes), dtype=int) for mode in modes}
            obs = {
                mode: np.array(
                    [max(0, transform(observed[group[g]]).max()) for g in genes]
                )
                for mode, transform in modes.items()
            }
            for start in range(0, permutations, 500):
                batch = min(500, permutations - start)
                signs = (
                    rng.integers(0, 2, size=(batch, len(delta)), dtype=np.int8) * 2 - 1
                ).astype(np.float32)
                simulated = signs @ standardized
                for mode, transform in modes.items():
                    values = transform(simulated)
                    for j, g in enumerate(genes):
                        best = np.maximum(0, values[:, group[g]].max(1))
                        counters[mode][j] += np.count_nonzero(
                            best >= obs[mode][j] - 1e-5
                        )
            for mode, transform in modes.items():
                p = (counters[mode] + 1) / (permutations + 1)
                q = multipletests(p, method="fdr_bh")[1]
                for j, g in enumerate(genes):
                    ids = group[g]
                    focal = ids[np.argmax(transform(observed[ids]))]
                    rows.append(
                        dict(
                            gate="TSS_ge_3",
                            contrast=contrast,
                            mode=mode,
                            gene=g,
                            n_peaks=len(ids),
                            n_pairs=len(delta),
                            max_statistic=obs[mode][j],
                            permutation_p=p[j],
                            q_25_genes=q[j],
                            top_peak=names[focal],
                            top_peak_signed_statistic=observed[focal],
                            permutation_exceedances=int(counters[mode][j]),
                            MC_p_lower=stats.beta.ppf(
                                0.025,
                                counters[mode][j] + 1,
                                permutations - counters[mode][j] + 1,
                            ),
                            MC_p_upper=stats.beta.ppf(
                                0.975,
                                counters[mode][j] + 1,
                                permutations - counters[mode][j] + 1,
                            ),
                            permutations=permutations,
                            seed=seed,
                        )
                    )
            print(contrast, "completed", permutations, "permutations", flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(out / "GSE208248_gene_region_tests.tsv", sep="\t", index=False)
    print(
        result.groupby(["contrast", "mode"])
        .agg(min_p=("permutation_p", "min"), min_q=("q_25_genes", "min"))
        .to_string()
    )
    export_rna(root, work)


def export_rna(root, work):
    out = root / "results/local_signal_decision"
    genes = pd.read_csv(
        root / "results/differentiation_fusion25/genes.tsv", sep="\t"
    ).gene.tolist()
    gi = pd.Index((work / "genes.txt").read_text().splitlines()).get_indexer(genes)
    pairs = pd.read_csv(root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[pairs.gate.eq("TSS_ge_3")]
    for contrast, ps in pairs.groupby("contrast"):
        highs, lows = [], []
        for gsm, d in ps.groupby("gsm"):
            meta = pd.read_csv(work / f"{gsm}_meta.tsv", sep="\t").set_index("barcode")
            hi, lo = meta.index.get_indexer(d.high_barcode), meta.index.get_indexer(
                d.low_barcode
            )
            assert min(hi) >= 0 and min(lo) >= 0
            x = sparse.load_npz(work / f"{gsm}_rna.npz")[gi]
            highs.append(
                x[:, hi].toarray()
                / meta.iloc[hi].total_rna_umi.to_numpy()[None, :]
                * 10000
            )
            lows.append(
                x[:, lo].toarray()
                / meta.iloc[lo].total_rna_umi.to_numpy()[None, :]
                * 10000
            )
        high, low = np.hstack(highs), np.hstack(lows)
        p = np.nan_to_num(
            stats.ttest_rel(np.log1p(high), np.log1p(low), axis=1).pvalue, nan=1.0
        )
        fc = np.divide(
            high.mean(1),
            low.mean(1),
            out=np.full(len(genes), np.nan),
            where=low.mean(1) > 0,
        )
        result = pd.DataFrame(
            dict(
                gene=genes,
                contrast=contrast,
                gate="TSS_ge_3",
                n_pairs=high.shape[1],
                RNA_mean_CP10k_high=high.mean(1),
                RNA_mean_CP10k_low=low.mean(1),
                RNA_FC=fc,
                RNA_p=p,
                RNA_q25=multipletests(p, method="fdr_bh")[1],
            )
        )
        result.to_csv(
            out / f"GSE208248_{contrast}_RNA.tsv", sep="\t", index=False, na_rep="NA"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=2000000)
    parser.add_argument("--seed", type=int, default=20260928)
    a = parser.parse_args()
    main(a.root, a.work, a.permutations, a.seed)
