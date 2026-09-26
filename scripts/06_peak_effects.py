"""Outcome screen across every knowledge-based candidate peak, then gene rank."""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, read_barcodes


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5-root', type=Path, required=True)
    p.add_argument('--tables', type=Path, required=True)
    a = p.parse_args()
    pairs = pd.read_csv(a.tables / 'matched_pairs.tsv.gz', sep='\t')
    candidate = pd.read_csv(a.tables / 'candidate_peak_gene.tsv.gz', sep='\t')
    mapping = pd.read_csv(a.tables / 'consensus_peak_map.tsv.gz', sep='\t').set_index('peak')
    anchors = sorted(candidate.peak.unique())
    by_library = []
    for gsm, sample_pairs in pairs.groupby('gsm'):
        barcodes = sorted(set(sample_pairs.high_barcode) | set(sample_pairs.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            _, P, _, peak_names = read_barcodes(h5, barcodes)
        bidx = {b: i for i, b in enumerate(barcodes)}
        pidx = {str(p): i for i, p in enumerate(peak_names)}
        local = [pidx.get(anchor if gsm == 'GSM6339597' else mapping.at[anchor, gsm + '_peak'], -1)
                 for anchor in anchors]
        valid = np.asarray(local) >= 0
        available = np.asarray(anchors)[valid]
        X = P[:, np.asarray(local)[valid]].tocsr()
        for (gate, contrast), sub in sample_pairs.groupby(['gate', 'contrast']):
            hi = np.asarray([bidx[b] for b in sub.high_barcode])
            lo = np.asarray([bidx[b] for b in sub.low_barcode])
            high, low = X[hi], X[lo]
            h = np.asarray(high.sum(axis=0)).ravel().astype(int)
            l = np.asarray(low.sum(axis=0)).ravel().astype(int)
            both = np.asarray(high.multiply(low).sum(axis=0)).ravel().astype(int)
            plus, minus = h - both, l - both
            disc = plus + minus
            row = pd.DataFrame({'gsm': gsm, 'gate': gate, 'contrast': contrast,
                                'peak': available, 'n_pairs': len(sub),
                                'high_open': h, 'low_open': l,
                                'high_only': plus, 'low_only': minus,
                                'delta_pp': 100 * (h - l) / len(sub)})
            by_library.append(row)
        print(gsm, 'peaks', len(available), flush=True)
    lib = pd.concat(by_library, ignore_index=True)
    lib.to_csv(a.tables / 'candidate_peak_effects_by_library.tsv.gz', sep='\t', index=False, compression='gzip')
    rows = []
    for (gate, contrast, peak), sub in lib.groupby(['gate', 'contrast', 'peak']):
        if len(sub) < 3:
            continue
        h, l = int(sub.high_open.sum()), int(sub.low_open.sum())
        plus, minus = int(sub.high_only.sum()), int(sub.low_only.sum())
        disc = plus + minus
        pval = min(1., 2 * stats.binom.cdf(min(plus, minus), disc, .5)) if disc else 1.
        rows.append({'gate': gate, 'contrast': contrast, 'peak': peak,
                     'n_libraries': len(sub), 'n_pairs': int(sub.n_pairs.sum()),
                     'high_open': h, 'low_open': l, 'high_only': plus, 'low_only': minus,
                     'delta_pp': 100 * (h - l) / sub.n_pairs.sum(),
                     'positive_libraries': int((sub.delta_pp > 0).sum()),
                     'p_pair_binomial': pval})
    out = pd.DataFrame(rows)
    out['q_candidate_peaks'] = np.nan
    for _, sub in out.groupby(['gate', 'contrast']):
        test = (sub.high_open + sub.low_open) >= 20
        if test.any():
            out.loc[sub.index[test], 'q_candidate_peaks'] = multipletests(
                sub.loc[test, 'p_pair_binomial'], method='fdr_bh')[1]
    out.to_csv(a.tables / 'candidate_peak_effects_pooled.tsv.gz', sep='\t', index=False, compression='gzip')
    links = pd.read_csv(a.tables / 'peak_gene_links.tsv', sep='\t')
    primary = out[(out.gate == 'TSS_ge_3') & (out.contrast == '2plus_vs_1')]
    linked = links.merge(primary[['peak', 'delta_pp', 'positive_libraries', 'q_candidate_peaks']],
                         on='peak', how='left', validate='many_to_one')
    linked['linked_positive'] = (linked.q_all_links < .05) & (linked.partial_r > 0) & (linked.n_libraries == 4)
    linked['linked_coq_positive'] = linked.linked_positive & (linked.delta_pp > 0) & (linked.positive_libraries >= 3)
    rank = linked.groupby('gene').agg(n_linked=('linked_positive', 'sum'),
                                      n_linked_coq_3of4=('linked_coq_positive', 'sum')).reset_index()
    full = pd.read_csv(a.tables / 'candidate_gene_ranking.tsv', sep='\t')
    full = full.drop(columns=['n_positive_links']).merge(rank, on='gene', how='left')
    full[['n_linked', 'n_linked_coq_3of4']] = full[['n_linked', 'n_linked_coq_3of4']].fillna(0).astype(int)
    full = full.sort_values(['n_linked_coq_3of4', 'n_linked', 'difference'], ascending=False)
    full.to_csv(a.tables / 'candidate_gene_ranking.tsv', sep='\t', index=False)
    print('primary ATAC peak q<.05', int((primary.q_candidate_peaks < .05).sum()))
    print(full.head(12)[['gene', 'n_tested_links', 'n_linked', 'n_linked_coq_3of4',
                         'difference', 'q_221']].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
