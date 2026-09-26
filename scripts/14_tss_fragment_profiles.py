"""Tabulate TSS-aligned ATAC insertion profiles for the fixed 221-gene set.

This visual analysis uses the *same* primary matched nuclei as the peak tests,
but fetches deduplicated fragments directly rather than plotting 10x peak calls.
It is descriptive and does not participate in gene or peak nomination.

Requires pysam and the four indexed GSE208248 ATAC fragment files. A gene's
GENCODE v48 gene-feature TSS supplies one fixed alignment point; the primary
search itself considers all transcript TSSs within +/-100 kb.
"""

import argparse
import csv
import gzip
import io
import re
from collections import defaultdict
from pathlib import Path

import pysam


FRAGMENT_GSM = {
    "GSM6339597": "GSM6339598",
    "GSM6339599": "GSM6339600",
    "GSM6339601": "GSM6339602",
    "GSM6339603": "GSM6339604",
}
HALF_WINDOW = 5000
BIN_WIDTH = 100
N_BINS = 2 * HALF_WINDOW // BIN_WIDTH
CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def gene_tss(gtf: Path, selected: set[str]) -> dict[str, tuple[str, int, str]]:
    result = {}
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in CHROMOSOMES:
                continue
            if 'gene_type "protein_coding"' not in fields[8]:
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match or match.group(1) not in selected:
                continue
            tss_0based = int(fields[3] if fields[6] == "+" else fields[4]) - 1
            name = match.group(1)
            if name in result:
                raise ValueError(f"Multiple gene-feature TSS records for {name}")
            result[name] = (fields[0], tss_0based, fields[6])
    if set(result) != selected:
        raise ValueError(f"Missing GENCODE gene TSS: {sorted(selected - set(result))}")
    return result


def matched_barcodes(path: Path) -> dict[str, dict[str, int]]:
    result = defaultdict(dict)
    counts = defaultdict(int)
    with gzip.open(path, "rt") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["gate"] != "TSS_ge_3" or row["contrast"] != "2plus_vs_1":
                continue
            gsm = row["gsm"]
            for group, field in ((1, "high_barcode"), (0, "low_barcode")):
                barcode = row[field]
                if barcode in result[gsm]:
                    raise ValueError(f"Reused matched barcode: {gsm} {barcode}")
                result[gsm][barcode] = group
                counts[gsm, group] += 1
    if set(result) != set(FRAGMENT_GSM) or sum(counts[gsm, 1] for gsm in result) != 958:
        raise ValueError("Expected the four libraries and 958 primary pairs")
    for gsm in result:
        if counts[gsm, 1] != counts[gsm, 0]:
            raise ValueError(f"Unequal matched groups in {gsm}")
    return dict(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with args.genes.open(newline="") as stream:
        selected = {row["gene"] for row in csv.DictReader(stream, delimiter="\t")}
    if len(selected) != 221:
        raise ValueError("Expected fixed 221-gene Hallmark/Reactome union")
    tss = gene_tss(args.gtf, selected)
    barcodes = matched_barcodes(args.pairs)
    counts = {(gene, group): [0] * N_BINS for gene in selected for group in (0, 1)}

    for gsm, fragment_gsm in FRAGMENT_GSM.items():
        files = list(args.fragments_dir.glob(f"{fragment_gsm}_*fragments.tsv.gz"))
        if len(files) != 1 or not Path(str(files[0]) + ".tbi").exists():
            raise ValueError(f"Expected one indexed fragment file for {fragment_gsm}")
        with pysam.TabixFile(str(files[0])) as fragments:
            for gene in sorted(selected):
                chrom, site, strand = tss[gene]
                start = max(0, site - HALF_WINDOW - 10)
                end = site + HALF_WINDOW + 10
                for line in fragments.fetch(chrom, start, end):
                    fields = line.split("\t")
                    group = barcodes[gsm].get(fields[3])
                    if group is None:
                        continue
                    # One row is one PCR-deduplicated fragment. Use its two
                    # Tn5 insertion sites rather than its reported read support.
                    for cut in (int(fields[1]) + 4, int(fields[2]) - 5):
                        relative = cut - site if strand == "+" else site - cut
                        if -HALF_WINDOW <= relative < HALF_WINDOW:
                            index = (relative + HALF_WINDOW) // BIN_WIDTH
                            counts[gene, group][index] += 1
        print(f"Counted {gsm} fragments across 221 TSS windows", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(args.out), mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(
                [
                    "gene",
                    "bin_start_bp",
                    "bin_end_bp",
                    "high_cuts",
                    "low_cuts",
                    "n_high",
                    "n_low",
                ]
            )
            for gene in sorted(selected):
                for index in range(N_BINS):
                    writer.writerow(
                        [
                            gene,
                            -HALF_WINDOW + index * BIN_WIDTH,
                            -HALF_WINDOW + (index + 1) * BIN_WIDTH,
                            counts[gene, 1][index],
                            counts[gene, 0][index],
                            958,
                            958,
                        ]
                    )
    print(f"Wrote {len(selected)} genes x {N_BINS} bins: {args.out}", flush=True)


if __name__ == "__main__":
    main()
