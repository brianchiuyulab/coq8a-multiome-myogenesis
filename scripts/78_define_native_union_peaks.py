"""Define the 25-gene neighborhood directly in the GSE240061 author peak set."""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def main(root, work, gtf):
    out = root / "results/differentiation_fusion25/native_GSE240061"
    out.mkdir(parents=True, exist_ok=True)
    genes = set(
        pd.read_csv(root / "results/differentiation_fusion25/genes.tsv", sep="\t").gene
    )
    names = (work / "all_ATAC_features.txt").read_text().splitlines()
    peaks = pd.Series(names).str.extract(r"^(.+)-(\d+)-(\d+)$")
    peaks.columns = ["chrom", "start", "end"]
    peaks[["start", "end"]] = peaks[["start", "end"]].astype(int)
    peaks["peak"] = names
    peaks["index"] = np.arange(len(peaks))
    valid = peaks.chrom.isin([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
    valid &= (peaks.end - peaks.start).between(200, 2000)
    with gzip.open(root / "reference/hg38-blacklist.v2.bed.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            z = line.split()
            if len(z) >= 3:
                valid &= ~(
                    peaks.chrom.eq(z[0])
                    & peaks.start.lt(int(z[2]))
                    & peaks.end.gt(int(z[1]))
                )
    peaks = peaks[valid].copy()
    peaks["midpoint"] = (peaks.start + peaks.end) // 2
    tss = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            z = line.rstrip().split("\t")
            if z[2] != "transcript":
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', z[8]))
            if (
                a.get("gene_name") not in genes
                or a.get("gene_type") != "protein_coding"
            ):
                continue
            tss.append(
                dict(
                    chrom=z[0],
                    gene=a["gene_name"],
                    tss_1based=int(z[3] if z[6] == "+" else z[4]),
                )
            )
    tss = pd.DataFrame(tss).drop_duplicates()
    rows = []
    for (chrom, gene), d in tss.groupby(["chrom", "gene"]):
        p = peaks[peaks.chrom.eq(chrom)].copy()
        p["nearest_tss_bp"] = np.min(
            np.abs(p.midpoint.to_numpy()[:, None] - d.tss_1based.to_numpy()[None, :]),
            axis=1,
        )
        p = p[p.nearest_tss_bp.le(100000)]
        p["gene"] = gene
        rows.append(p)
    membership = pd.concat(rows, ignore_index=True)
    mapped = set(
        pd.read_csv(
            root / "results/differentiation_fusion25/GSE240061_peak_membership.tsv",
            sep="\t",
        ).author_peak
    )
    membership["previously_mapped"] = membership.peak.isin(mapped)
    membership.to_csv(out / "native_peak_membership.tsv", sep="\t", index=False)
    tss.to_csv(out / "gene_TSS.tsv", sep="\t", index=False)
    print(
        "Native peaks",
        membership.peak.nunique(),
        "genes",
        membership.gene.nunique(),
        "shared mapped",
        membership[membership.previously_mapped].peak.nunique(),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--gtf", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work, a.gtf)
