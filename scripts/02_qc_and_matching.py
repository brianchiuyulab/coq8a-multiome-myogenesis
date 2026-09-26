"""Joint RNA/ATAC QC and within-library depth matching, independent of outcomes."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

CONTRASTS = [('2plus_vs_1', 2, 1), ('3plus_vs_1', 3, 1)]
GATES = [('TSS_ge_2', 2), ('TSS_ge_3', 3)]
CALIPER = 0.10


def depth_match(cells, high_cut, low_count):
    hi = cells[cells.COQ8A_umi >= high_cut].reset_index(drop=True)
    lo = cells[cells.COQ8A_umi == low_count].reset_index(drop=True)
    if len(hi) == 0 or len(lo) == 0:
        return []
    hx = np.log1p(hi[['total_rna_umi', 'total_open_peaks']].to_numpy())
    lx = np.log1p(lo[['total_rna_umi', 'total_open_peaks']].to_numpy())
    nn = NearestNeighbors(n_neighbors=min(100, len(lo)), algorithm='kd_tree').fit(lx)
    dist, candidates = nn.kneighbors(hx)
    used, out = set(), []
    for i in np.argsort(-dist[:, 0]):
        for j in candidates[i]:
            if int(j) not in used and np.max(np.abs(hx[i] - lx[j])) <= CALIPER:
                used.add(int(j))
                out.append((hi.iloc[i], lo.iloc[j]))
                break
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nuclei', type=Path, required=True)
    parser.add_argument('--fragment-qc', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(args.nuclei, sep='\t')
    a['qc_rna_depth'] = a.total_rna_umi >= 500
    a['qc_atac_depth'] = a.total_open_peaks >= 500
    a = a[a.qc_rna_depth & a.qc_atac_depth].copy()
    # Upper-tail cap is computed after both lower-depth filters, within library.
    cap = a.groupby('gsm')[['total_rna_umi', 'total_open_peaks']].transform(lambda z: z.quantile(.95))
    a['qc_depth_upper'] = (a.total_rna_umi <= cap.total_rna_umi) & (a.total_open_peaks <= cap.total_open_peaks)
    a = a[a.qc_depth_upper].copy()
    for kind in ['tss', 'nucleosome', 'blacklist']:
        files = sorted(args.fragment_qc.glob(f'GSM6339*_fragment_{kind}_qc.tsv'))
        assert len(files) == 4, (kind, files)
        q = pd.concat((pd.read_csv(p, sep='\t') for p in files), ignore_index=True)
        a = a.merge(q, on=['gsm', 'barcode'], validate='one_to_one')
    assert a.tss_enrichment.notna().all()
    a['frip'] = a.atac_peak_region_fragments / a.atac_fragments
    a['blacklist_ratio'] = a.blacklist_fragments / a.atac_fragments
    a['pass_atac_other'] = ((a.nucleosome_signal < 4) & (a.frip >= .25) &
                            (a.blacklist_ratio < .05))
    a['pass_tss2'] = a.pass_atac_other & (a.tss_enrichment >= 2)
    a['pass_tss3'] = a.pass_atac_other & (a.tss_enrichment >= 3)
    a.to_csv(out / 'qc_nuclei.tsv.gz', sep='\t', index=False, compression='gzip')
    rows, pair_rows = [], []
    for gsm, lib in a.groupby('gsm'):
        rows.append({'gsm': gsm, 'source': lib.source.iloc[0], 'stage': lib.stage.iloc[0],
                     'n_pre_fragment_qc': len(lib), 'n_tss2': int(lib.pass_tss2.sum()),
                     'n_tss3': int(lib.pass_tss3.sum()),
                     'median_tss': lib.tss_enrichment.median(),
                     'median_frip': lib.frip.median(),
                     'median_nucleosome': lib.nucleosome_signal.median()})
        for gate, threshold in GATES:
            eligible = lib[lib[f'pass_tss{threshold}']]
            for contrast, high_cut, low_count in CONTRASTS:
                for pair_id, (h, l) in enumerate(depth_match(eligible, high_cut, low_count)):
                    pair_rows.append({'gsm': gsm, 'source': h.source, 'stage': h.stage,
                                      'gate': gate, 'contrast': contrast, 'pair': pair_id,
                                      'high_barcode': h.barcode, 'low_barcode': l.barcode,
                                      'high_coq_umi': int(h.COQ8A_umi),
                                      'low_coq_umi': int(l.COQ8A_umi),
                                      'abs_log_rna_depth_difference': abs(np.log1p(h.total_rna_umi) - np.log1p(l.total_rna_umi)),
                                      'abs_log_atac_depth_difference': abs(np.log1p(h.total_open_peaks) - np.log1p(l.total_open_peaks))})
    pd.DataFrame(rows).to_csv(out / 'qc_by_library.tsv', sep='\t', index=False)
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(out / 'matched_pairs.tsv.gz', sep='\t', index=False, compression='gzip')
    print(pd.DataFrame(rows).to_string(index=False))
    print(pairs.groupby(['gate', 'contrast', 'gsm']).size().to_string())


if __name__ == '__main__':
    main()
