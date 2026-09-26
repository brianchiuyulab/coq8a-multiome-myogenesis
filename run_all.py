"""Run the fixed GSE208248 COQ8A same-nucleus RNA/ATAC workflow.

Raw 10x H5 files and the GENCODE GTF are supplied by the caller. The smaller
fragment QC, motif and Scrublet reference tables are versioned in reference/.
"""
import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5-root', type=Path, required=True,
                   help='Directory with four GEO filtered_feature_bc_matrix.h5 files and barcode metrics')
    p.add_argument('--gtf', type=Path, required=True, help='GENCODE v48 hg38 annotation.gtf.gz')
    p.add_argument('--out-root', type=Path, default=ROOT,
                   help='Output repository root; default is this directory')
    args = p.parse_args()
    t = args.out_root / 'results' / 'tables'
    f = args.out_root / 'figures'
    t.mkdir(parents=True, exist_ok=True)
    f.mkdir(parents=True, exist_ok=True)
    r = ROOT / 'reference'
    h = args.h5_root.resolve()
    g = args.gtf.resolve()
    steps = [
        ('01_extract_nuclei.py', '--h5-root', h, '--out', t / 'nuclei.tsv.gz'),
        ('02_qc_and_matching.py', '--nuclei', t / 'nuclei.tsv.gz', '--fragment-qc', r / 'fragment_qc', '--out', t),
        ('03_define_regions.py', '--h5-root', h, '--gtf', g, '--gene-list', r / 'myogenesis_221_gene_sources.tsv',
         '--blacklist', r / 'hg38-blacklist.v2.bed.gz', '--out', t),
        ('04_gene_programme_effects.py', '--h5-root', h, '--tables', t, '--genes', r / 'myogenesis_221_gene_sources.tsv'),
        ('05_peak_gene_links.py', '--h5-root', h, '--tables', t, '--motif-scores', r / 'jaspar2024_peak_scores.tsv.gz'),
        ('05_peak_gene_links.py', '--h5-root', h, '--tables', t, '--motif-scores', r / 'jaspar2024_peak_scores.tsv.gz',
         '--gate', 'TSS_ge_2'),
        ('06_peak_effects.py', '--h5-root', h, '--tables', t),
        ('10_link_gate_reconciliation.py', '--tables', t),
        ('11_primary_decision.py', '--tables', t),
        ('07_locus_sensitivity.py', '--h5-root', h, '--tables', t),
        ('08_doublet_sensitivity.py', '--tables', t, '--scrublet', r / 'rna_scrublet_qc.tsv.gz'),
        ('09_make_figures.py', '--tables', t, '--figures', f),
    ]
    for name, *arguments in steps:
        cmd = [sys.executable, str(ROOT / 'scripts' / name), *map(str, arguments)]
        print('RUN', name, flush=True)
        subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
