"""Derive fragment-QC barcode and TSS inputs from the public primary files.

The barcode table contains nuclei that pass the two lower-depth thresholds and
the within-library 95th-percentile upper-depth caps. Fragment-level QC is
computed only for this fixed set to avoid processing excluded barcodes.
"""

import argparse
import gzip
from pathlib import Path

import pandas as pd


CANONICAL_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def prepare_barcodes(nuclei_path: Path) -> pd.DataFrame:
    nuclei = pd.read_csv(nuclei_path, sep="\t")
    eligible = nuclei[(nuclei.total_rna_umi >= 500) & (nuclei.total_open_peaks >= 500)].copy()
    caps = eligible.groupby("gsm")[["total_rna_umi", "total_open_peaks"]].transform(
        lambda values: values.quantile(0.95)
    )
    eligible = eligible[
        (eligible.total_rna_umi <= caps.total_rna_umi)
        & (eligible.total_open_peaks <= caps.total_open_peaks)
    ]
    return eligible[["gsm", "barcode"]].sort_values(["gsm", "barcode"])


def prepare_tss(gtf_path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf_path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if (
                len(fields) != 9
                or fields[2] != "gene"
                or fields[0] not in CANONICAL_CHROMOSOMES
                or 'gene_type "protein_coding"' not in fields[8]
            ):
                continue
            position = int(fields[3] if fields[6] == "+" else fields[4])
            rows.append((fields[0], position))
    return pd.DataFrame(rows, columns=["chr", "tss_1based"]).drop_duplicates()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nuclei", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    barcodes = prepare_barcodes(args.nuclei)
    tss = prepare_tss(args.gtf)
    barcodes.to_csv(args.out / "gse208248_qc_barcodes.tsv", sep="\t", index=False)
    tss.to_csv(args.out / "gencode_v48_protein_coding_gene_tss.tsv", sep="\t", index=False)
    print(f"Prepared {len(barcodes):,} barcodes and {len(tss):,} gene TSS entries")


if __name__ == "__main__":
    main()
