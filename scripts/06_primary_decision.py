"""Summarize the fixed 221-gene, TSS>=3, COQ8A>=3-vs-1 scan.

This runs before interpreting any MYOD1 subset. It reads the complete tested
families and records whether any ATAC result survives the stated FDR rule.
"""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", type=Path, required=True)
    a = p.parse_args()
    genes = pd.read_csv(a.tables / "gene_effects_pooled.tsv", sep="\t")
    genes = genes[(genes.gate == "TSS_ge_3") & (genes.contrast == "3plus_vs_1")]
    peaks = pd.read_csv(a.tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peaks = peaks[(peaks.gate == "TSS_ge_3") & (peaks.contrast == "3plus_vs_1")]
    tested_peaks = peaks[peaks.q_candidate_peaks.notna()]
    programme = pd.read_csv(a.tables / "programme_effects_pooled.tsv", sep="\t")
    programme = programme[(programme.gate == "TSS_ge_3") & (programme.contrast == "3plus_vs_1")]
    atac = genes[genes.modality == "ATAC"]
    rna = genes[genes.modality == "RNA"]
    rows = [
        ("primary_gate", "TSS_ge_3"),
        ("primary_contrast", "3plus_vs_1"),
        ("n_candidate_genes", len(atac)),
        ("n_atac_gene_regions_q_lt_0_05", int((atac.q_221 < 0.05).sum())),
        ("min_atac_gene_region_q", atac.q_221.min()),
        ("n_candidate_peaks", len(peaks)),
        ("n_tested_candidate_peaks", len(tested_peaks)),
        ("n_candidate_peaks_q_lt_0_05", int((tested_peaks.q_candidate_peaks < 0.05).sum())),
        ("min_candidate_peak_q", tested_peaks.q_candidate_peaks.min()),
        ("n_rna_genes_q_lt_0_05", int((rna.q_221 < 0.05).sum())),
    ]
    for modality in ["RNA", "ATAC"]:
        for name in ["Hallmark_myogenesis", "Reactome_myogenesis", "MRF_loci", "MEF2_loci"]:
            row = programme[(programme.modality == modality) & (programme.programme == name)].iloc[
                0
            ]
            rows.append((f"{modality}_{name}_difference", row.difference))
            rows.append((f"{modality}_{name}_q_4", row.q_4_programmes))
            rows.append((f"{modality}_{name}_positive_libraries", row.positive_libraries))
    out = pd.DataFrame(rows, columns=["metric", "value"])
    out.to_csv(a.tables / "primary_decision_summary.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
