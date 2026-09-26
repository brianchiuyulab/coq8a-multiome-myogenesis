"""COQ8A-associated RNA and common-peak accessibility across all 221 genes.

The primary comparison is TSS>=3, COQ8A>=2 versus 1 UMI. Three other gate/
count combinations are reported as sensitivity. All p values use paired nuclei
and are descriptive, because only two source cell lines are available.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, paired_stats, read_barcodes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--tables", type=Path, required=True)
    p.add_argument("--genes", type=Path, required=True)
    a = p.parse_args()
    genes = pd.read_csv(a.genes, sep="\t")
    names = genes.gene.tolist()
    gene_idx = {gene: j for j, gene in enumerate(names)}
    groups = {
        "Hallmark_myogenesis": genes.loc[genes.Hallmark_Myogenesis, "gene"].tolist(),
        "Reactome_myogenesis": genes.loc[genes.Reactome_Myogenesis, "gene"].tolist(),
        "MRF_loci": ["MYOD1", "MYOG", "MYF5", "MYF6"],
        "MEF2_loci": ["MEF2A", "MEF2C", "MEF2D"],
    }
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    nuclei = pd.read_csv(a.tables / "qc_nuclei.tsv.gz", sep="\t")
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    candidates = candidates[candidates.n_libraries == 4]
    consensus = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t")
    cmap = consensus.set_index("peak")
    assigned = candidates.groupby("gene").peak.apply(list).to_dict()
    gene_rows, programme_rows, pair_programme = [], [], []
    pooled = {
        (gate, contrast, modality, gene): []
        for gate in ["TSS_ge_2", "TSS_ge_3"]
        for contrast in ["2plus_vs_1", "3plus_vs_1"]
        for modality in ["RNA", "ATAC"]
        for gene in names
    }
    pooled_programme = {}
    for gsm, gsm_pairs in pairs.groupby("gsm"):
        barcodes = sorted(set(gsm_pairs.high_barcode) | set(gsm_pairs.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            R, P, rna_names, peak_names = read_barcodes(h5, barcodes)
        barcode_idx = {b: i for i, b in enumerate(barcodes)}
        rna_idx = {g: i for i, g in enumerate(rna_names)}
        local_peak_idx = {str(x): i for i, x in enumerate(peak_names)}
        depth = (
            nuclei[nuclei.gsm == gsm].set_index("barcode").loc[barcodes, "total_rna_umi"].to_numpy()
        )
        rna_cols = [rna_idx.get(g, -1) for g in names]
        assert min(rna_cols) >= 0, f"{gsm}: one or more genes absent from RNA matrix"
        Y = R[:, rna_cols].toarray()
        Y = np.log1p(10000 * Y / depth[:, None])
        peak_row, gene_col, weights = [], [], []
        n_regions = np.zeros(len(names), dtype=int)
        for gene, anchor_peaks in assigned.items():
            j = gene_idx[gene]
            local = []
            for anchor in anchor_peaks:
                target = anchor if gsm == "GSM6339597" else cmap.at[anchor, gsm + "_peak"]
                if isinstance(target, str) and target in local_peak_idx:
                    local.append(local_peak_idx[target])
            local = sorted(set(local))
            n_regions[j] = len(local)
            for i in local:
                peak_row.append(i)
                gene_col.append(j)
                weights.append(1 / len(local))
        region_map = sparse.csr_matrix(
            (weights, (peak_row, gene_col)), shape=(len(peak_names), len(names))
        )
        A = (P @ region_map).toarray()
        for (gate, contrast), sub in gsm_pairs.groupby(["gate", "contrast"]):
            hi = np.asarray([barcode_idx[b] for b in sub.high_barcode])
            lo = np.asarray([barcode_idx[b] for b in sub.low_barcode])
            n = len(hi)
            assert n == len(lo)
            for modality, values in [("RNA", Y), ("ATAC", A)]:
                diff = values[hi] - values[lo]
                for j, gene in enumerate(names):
                    if modality == "ATAC" and not n_regions[j]:
                        continue
                    s = paired_stats(diff[:, j])
                    gene_rows.append(
                        {
                            "gsm": gsm,
                            "source": sub.source.iloc[0],
                            "stage": sub.stage.iloc[0],
                            "gate": gate,
                            "contrast": contrast,
                            "modality": modality,
                            "gene": gene,
                            "n_regions": int(n_regions[j]),
                            "high_mean": float(values[hi, j].mean()),
                            "low_mean": float(values[lo, j].mean()),
                            **s,
                        }
                    )
                    pooled[(gate, contrast, modality, gene)].append(diff[:, j])
                for programme, gene_names in groups.items():
                    use = [
                        gene_idx[g]
                        for g in gene_names
                        if modality == "RNA" or n_regions[gene_idx[g]] > 0
                    ]
                    d = diff[:, use].mean(axis=1)
                    s = paired_stats(d)
                    programme_rows.append(
                        {
                            "gsm": gsm,
                            "source": sub.source.iloc[0],
                            "stage": sub.stage.iloc[0],
                            "gate": gate,
                            "contrast": contrast,
                            "modality": modality,
                            "programme": programme,
                            "n_genes": len(use),
                            **s,
                        }
                    )
                    pooled_programme.setdefault((gate, contrast, modality, programme), []).append(d)
                    pair_programme.extend(
                        {
                            "gsm": gsm,
                            "gate": gate,
                            "contrast": contrast,
                            "modality": modality,
                            "programme": programme,
                            "pair": int(pid),
                            "difference": float(delta),
                        }
                        for pid, delta in zip(sub.pair, d)
                    )
        print(gsm, "effects complete", flush=True)
    gd = pd.DataFrame(gene_rows)
    gd.to_csv(
        a.tables / "gene_effects_by_library.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(programme_rows).to_csv(
        a.tables / "programme_effects_by_library.tsv", sep="\t", index=False
    )
    pd.DataFrame(pair_programme).to_csv(
        a.tables / "pair_programme_effects.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    gene_summary = []
    for (gate, contrast, modality, gene), arrays in pooled.items():
        if not arrays:
            continue
        d = np.concatenate(arrays)
        sample = gd[
            (gd.gate == gate)
            & (gd.contrast == contrast)
            & (gd.modality == modality)
            & (gd.gene == gene)
        ]
        high = np.average(sample.high_mean, weights=sample.n_pairs)
        low = np.average(sample.low_mean, weights=sample.n_pairs)
        gene_summary.append(
            {
                "gate": gate,
                "contrast": contrast,
                "modality": modality,
                "gene": gene,
                "n_regions": int(sample.n_regions.max()),
                "high_mean": high,
                "low_mean": low,
                "positive_libraries": int((sample.difference > 0).sum()),
                **paired_stats(d),
            }
        )
    gene_summary = pd.DataFrame(gene_summary)
    gene_summary["q_221"] = np.nan
    for _, sub in gene_summary.groupby(["gate", "contrast", "modality"]):
        gene_summary.loc[sub.index, "q_221"] = multipletests(sub.p_pair.fillna(1), method="fdr_bh")[
            1
        ]
    gene_summary.to_csv(a.tables / "gene_effects_pooled.tsv", sep="\t", index=False)
    programme_summary = []
    for (gate, contrast, modality, programme), arrays in pooled_programme.items():
        d = np.concatenate(arrays)
        library = pd.DataFrame(programme_rows)
        sample = library[
            (library.gate == gate)
            & (library.contrast == contrast)
            & (library.modality == modality)
            & (library.programme == programme)
        ]
        programme_summary.append(
            {
                "gate": gate,
                "contrast": contrast,
                "modality": modality,
                "programme": programme,
                "n_genes_min": int(sample.n_genes.min()),
                "n_genes_max": int(sample.n_genes.max()),
                "positive_libraries": int((sample.difference > 0).sum()),
                **paired_stats(d),
            }
        )
    programme_summary = pd.DataFrame(programme_summary)
    programme_summary["q_4_programmes"] = np.nan
    for _, sub in programme_summary.groupby(["gate", "contrast", "modality"]):
        programme_summary.loc[sub.index, "q_4_programmes"] = multipletests(
            sub.p_pair.fillna(1), method="fdr_bh"
        )[1]
    programme_summary.to_csv(a.tables / "programme_effects_pooled.tsv", sep="\t", index=False)
    print(
        programme_summary[
            (programme_summary.gate == "TSS_ge_3") & (programme_summary.contrast == "2plus_vs_1")
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
