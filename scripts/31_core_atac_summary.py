"""Summarize existing main-contrast ATAC effects without refitting or reselection."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(args.tables / 'gene_search_space.tsv', sep='\t')
    effects = pd.read_csv(args.tables / 'gene_effects_pooled.tsv', sep='\t')
    programmes = pd.read_csv(args.tables / 'programme_effects_pooled.tsv', sep='\t')

    def main_atac(frame):
        return frame.loc[
            frame.gate.eq('TSS_ge_3')
            & frame.contrast.eq('3plus_vs_1')
            & frame.modality.eq('ATAC')
        ].copy()

    effects = main_atac(effects)
    programmes = main_atac(programmes)
    members = {
        'Hallmark_myogenesis': genes.loc[genes.Hallmark_Myogenesis, 'gene'],
        'Reactome_myogenesis': genes.loc[genes.Reactome_Myogenesis, 'gene'],
        'MRF_loci': ['MYOD1', 'MYOG', 'MYF5', 'MYF6'],
        'MEF2_loci': ['MEF2A', 'MEF2C', 'MEF2D'],
    }
    for index, row in programmes.iterrows():
        subset = effects.loc[effects.gene.isin(members[row.programme])]
        assert len(subset) == row.n_genes_min == row.n_genes_max
        programmes.loc[index, 'fold_open'] = subset.high_mean.mean() / subset.low_mean.mean()
    programmes.to_csv(args.out / 'core_programme_ATAC_summary.tsv', sep='\t', index=False)
    core = effects.loc[effects.gene.isin(members['MRF_loci'] + members['MEF2_loci'] + ['MEF2B'])].copy()
    core['fold_open'] = core.high_mean / core.low_mean
    core.to_csv(args.out / 'core_TF_ATAC_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
