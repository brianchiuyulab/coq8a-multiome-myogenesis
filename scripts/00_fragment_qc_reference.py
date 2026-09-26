"""Per-nucleus fragment QC from GSE208248 ATAC fragments (run with Linux pysam).

TSS score follows Signac fast TSSEnrichment geometry: per-base insertions in
the central 1001 bp divided by mean per-base insertion signal in the two
outermost 100-bp flanks of a +/-1000-bp TSS window. TSS annotations are
GENCODE v48 protein-coding gene TSSs (one gene-level TSS per gene).
Nucleosome signal is exact mononucleosomal [147,294) / nucleosome-free <147
fragment ratio for selected barcodes, rather than Signac's optional sampling.
"""

from collections import defaultdict
import csv
import gzip
from pathlib import Path
import sys
import time
import argparse

import pysam

ATAC_TO_GEX = {
    "GSM6339598": "GSM6339597",
    "GSM6339600": "GSM6339599",
    "GSM6339602": "GSM6339601",
    "GSM6339604": "GSM6339603",
}


def main(atac_gsm, mode, R, fragment_root, out_root):
    gex_gsm = ATAC_TO_GEX[atac_gsm]
    with (R / "gse208248_qc_barcodes.tsv").open() as f:
        cells = {
            row["barcode"] for row in csv.DictReader(f, delimiter="\t") if row["gsm"] == gex_gsm
        }
    path = next(fragment_root.glob(atac_gsm + "*fragments.tsv.gz"))
    tbx = pysam.TabixFile(str(path), index=str(path) + ".tbi")
    print(atac_gsm, mode, "selected cells", len(cells), "contigs", len(tbx.contigs), flush=True)
    if mode == "tss":
        tss = []
        with (R / "gencode_v48_protein_coding_gene_tss.tsv").open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["chr"] in tbx.contigs:
                    tss.append((row["chr"], int(row["tss_1based"]) - 1))
        center = defaultdict(int)
        flank = defaultdict(int)
        for i, (chrom, t) in enumerate(tss):
            for row in tbx.fetch(chrom, max(0, t - 1000), t + 1001):
                fields = row.split("\t")
                bc = fields[3]
                if bc not in cells:
                    continue
                for cut in (int(fields[1]), int(fields[2]) - 1):
                    dist = abs(cut - t)
                    if dist <= 500:
                        center[bc] += 1
                    if 901 <= dist <= 1000:
                        flank[bc] += 1
            if (i + 1) % 2000 == 0:
                print(atac_gsm, "TSS", i + 1, "/", len(tss), flush=True)
        flank_mean = sum(flank.values()) / (200 * len(cells))
        out = out_root / (gex_gsm + "_fragment_tss_qc.tsv")
        with out.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                [
                    "gsm",
                    "barcode",
                    "n_tss",
                    "center_insertions",
                    "flank_insertions",
                    "tss_enrichment",
                ]
            )
            for bc in sorted(cells):
                expected = flank[bc] / 200 if flank[bc] else flank_mean
                score = (center[bc] / 1001) / expected if expected > 0 else 0
                w.writerow([gex_gsm, bc, len(tss), center[bc], flank[bc], score])
        print(atac_gsm, "TSS done", out, "flank mean", flank_mean, flush=True)
    elif mode == "nuc":
        free = defaultdict(int)
        mono = defaultdict(int)
        n = 0
        for chrom in tbx.contigs:
            for row in tbx.fetch(chrom):
                fields = row.split("\t")
                bc = fields[3]
                if bc not in cells:
                    continue
                length = int(fields[2]) - int(fields[1])
                if length < 147:
                    free[bc] += 1
                elif length < 294:
                    mono[bc] += 1
                n += 1
                if n % 10000000 == 0:
                    print(atac_gsm, "selected fragments", n, flush=True)
        out = out_root / (gex_gsm + "_fragment_nucleosome_qc.tsv")
        with out.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                ["gsm", "barcode", "nucleosome_free", "mononucleosomal", "nucleosome_signal"]
            )
            for bc in sorted(cells):
                score = mono[bc] / free[bc] if free[bc] else float("nan")
                w.writerow([gex_gsm, bc, free[bc], mono[bc], score])
        print(atac_gsm, "nucleosome done", out, "selected fragment rows", n, flush=True)
    elif mode == "blacklist":
        counts = defaultdict(int)
        with gzip.open(R / "hg38-blacklist.v2.bed.gz", "rt") as f:
            intervals = [
                (x[0], int(x[1]), int(x[2]))
                for line in f
                if (x := line.rstrip().split("\t")) and len(x) >= 3
            ]
        for chrom, start, end in intervals:
            if chrom not in tbx.contigs:
                continue
            for row in tbx.fetch(chrom, start, end):
                fields = row.split("\t")
                bc = fields[3]
                if bc in cells:
                    counts[bc] += 1
        out = out_root / (gex_gsm + "_fragment_blacklist_qc.tsv")
        with out.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["gsm", "barcode", "blacklist_fragments"])
            for bc in sorted(cells):
                w.writerow([gex_gsm, bc, counts[bc]])
        print(atac_gsm, "blacklist done", out, flush=True)
    else:
        raise ValueError(mode)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("atac_gsm", choices=ATAC_TO_GEX)
    p.add_argument("mode", choices=["tss", "nuc", "blacklist"])
    p.add_argument(
        "--reference-root",
        type=Path,
        required=True,
        help="Directory with gse208248_qc_barcodes.tsv, GENCODE TSS TSV and blacklist",
    )
    p.add_argument("--fragment-root", type=Path, required=True)
    p.add_argument("--out-root", type=Path, required=True)
    a = p.parse_args()
    a.out_root.mkdir(parents=True, exist_ok=True)
    main(a.atac_gsm, a.mode, a.reference_root, a.fragment_root, a.out_root)
