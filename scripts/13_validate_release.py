"""Check the analysis tables and figure inventory for internal consistency."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PRIMARY_GATE = "TSS_ge_3"
PRIMARY_CONTRAST = "2plus_vs_1"
SIX_REGION_SET = "positive_link_nonMRF_tss2_6"


def check_pairs(tables: Path) -> int:
    nuclei = pd.read_csv(tables / "qc_nuclei.tsv.gz", sep="\t")
    pairs = pd.read_csv(tables / "matched_pairs.tsv.gz", sep="\t")
    if nuclei.duplicated(["gsm", "barcode"]).any():
        raise ValueError("QC nucleus identifiers are not unique")
    nucleus_index = nuclei.set_index(["gsm", "barcode"])
    for (gate, contrast, gsm), group in pairs.groupby(["gate", "contrast", "gsm"]):
        if (
            group.high_barcode.duplicated().any()
            or group.low_barcode.duplicated().any()
        ):
            raise ValueError(f"Nucleus reused within {gsm} {gate} {contrast}")
        if set(group.high_barcode) & set(group.low_barcode):
            raise ValueError(f"High/low nucleus overlap in {gsm} {gate} {contrast}")
        high = nucleus_index.loc[[(gsm, barcode) for barcode in group.high_barcode]]
        low = nucleus_index.loc[[(gsm, barcode) for barcode in group.low_barcode]]
        threshold = 3 if contrast == "3plus_vs_1" else 2
        tss_gate = 3 if gate == PRIMARY_GATE else 2
        if not (high.COQ8A_umi.ge(threshold).all() and low.COQ8A_umi.eq(1).all()):
            raise ValueError(
                f"COQ8A group definition failed for {gsm} {gate} {contrast}"
            )
        if not (
            high.tss_enrichment.ge(tss_gate).all()
            and low.tss_enrichment.ge(tss_gate).all()
        ):
            raise ValueError(f"TSS gate failed for {gsm} {gate} {contrast}")
        for field in ("abs_log_rna_depth_difference", "abs_log_atac_depth_difference"):
            if group[field].gt(0.1000001).any():
                raise ValueError(f"Matching caliper failed for {gsm} {gate} {contrast}")
    primary = pairs[(pairs.gate == PRIMARY_GATE) & (pairs.contrast == PRIMARY_CONTRAST)]
    return len(primary)


def check_primary(tables: Path) -> tuple[int, int, int]:
    genes = pd.read_csv(tables / "gene_effects_pooled.tsv", sep="\t")
    genes = genes[(genes.gate == PRIMARY_GATE) & (genes.contrast == PRIMARY_CONTRAST)]
    atac_genes = genes[genes.modality == "ATAC"]
    rna_genes = genes[genes.modality == "RNA"]
    if not (len(atac_genes) == len(rna_genes) == 221):
        raise ValueError(
            "Primary RNA and ATAC gene families must each have 221 members"
        )
    peaks = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    peaks = peaks[(peaks.gate == PRIMARY_GATE) & (peaks.contrast == PRIMARY_CONTRAST)]
    tested = peaks[peaks.q_candidate_peaks.notna()]
    summary = pd.read_csv(tables / "primary_decision_summary.tsv", sep="\t").set_index(
        "metric"
    )
    expected = {
        "n_candidate_genes": len(atac_genes),
        "n_tested_candidate_peaks": len(tested),
        "n_atac_gene_regions_q_lt_0_05": int(atac_genes.q_221.lt(0.05).sum()),
        "n_candidate_peaks_q_lt_0_05": int(tested.q_candidate_peaks.lt(0.05).sum()),
    }
    for metric, observed in expected.items():
        if int(summary.at[metric, "value"]) != observed:
            raise ValueError(
                f"Primary summary disagrees with complete test family: {metric}"
            )
    return len(atac_genes), len(tested), expected["n_candidate_peaks_q_lt_0_05"]


def check_six_region_effect(tables: Path) -> tuple[int, float]:
    detail = pd.read_csv(tables / "myod1_locus_pair_scores.tsv.gz", sep="\t")
    detail = detail[
        (detail.gate == PRIMARY_GATE)
        & (detail.contrast == "3plus_vs_1")
        & (detail.region_set == SIX_REGION_SET)
    ]
    summary = pd.read_csv(tables / "myod1_locus_summary.tsv", sep="\t")
    row = summary[
        (summary.gate == PRIMARY_GATE)
        & (summary.contrast == "3plus_vs_1")
        & (summary.region_set == SIX_REGION_SET)
    ].squeeze()
    if len(detail) != int(row.n_pairs) or int(row.n_regions) != 6:
        raise ValueError("Six-region effect has inconsistent pair or region counts")
    fold = detail.high_open.mean() / detail.low_open.mean()
    delta = 100 * (detail.high_open.mean() - detail.low_open.mean())
    if not np.isclose(fold, row.fold_open, atol=1e-10):
        raise ValueError("Six-region fold differs from pair-level data")
    if not np.isclose(delta, row.delta_pp, atol=1e-10):
        raise ValueError(
            "Six-region percentage-point effect differs from pair-level data"
        )
    return len(detail), fold


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    n_pairs = check_pairs(args.tables)
    n_genes, n_peaks, n_significant = check_primary(args.tables)
    n_extreme, fold = check_six_region_effect(args.tables)
    figures = [
        args.figures / "main/Figure_1_global_discovery.pdf",
        args.figures / "main/Figure_2_Global_ATAC_Scan.pdf",
        args.figures / "supplement/Supplementary_Figure_MYOD1_Exploratory.pdf",
        args.figures / "supplement/Supplementary_Figure_QC.pdf",
    ]
    missing = [str(path) for path in figures if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing figures: {missing}")
    print(
        f"Validated {n_pairs} primary pairs, {n_genes} genes, {n_peaks} tested peaks "
        f"({n_significant} FDR-positive); six-region fold {fold:.3f} in {n_extreme} pairs."
    )


if __name__ == "__main__":
    main()
