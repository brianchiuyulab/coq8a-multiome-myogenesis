"""Annotate native candidate intervals against unrestricted GENCODE transcripts."""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def main(root, gtf):
    out = root / "results/differentiation_fusion25/native_GSE240061"
    m = pd.read_csv(out / "native_peak_membership.tsv", sep="\t")
    chromosomes = set(m.chrom)
    rows = []
    bodies = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            z = line.rstrip().split("\t")
            if z[0] not in chromosomes or z[2] not in ["gene", "transcript"]:
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', z[8]))
            lo, hi = int(z[3]) - 1, int(z[4])
            if z[2] == "gene":
                bodies.append((z[0], lo, hi, a.get("gene_name", a["gene_id"])))
            else:
                rows.append(
                    (
                        z[0],
                        lo if z[6] == "+" else hi - 1,
                        a.get("gene_name", a["gene_id"]),
                        a.get("transcript_type"),
                        a.get("transcript_id"),
                    )
                )
    tss = pd.DataFrame(
        rows, columns=["chrom", "tss", "gene", "transcript_type", "transcript"]
    )
    body = pd.DataFrame(bodies, columns=["chrom", "start", "end", "gene"])
    result = []
    for chrom, d in m.drop_duplicates("peak").groupby("chrom"):
        ts = tss[tss.chrom.eq(chrom)].sort_values(["tss", "gene"])
        g = body[body.chrom.eq(chrom)]
        for p in d.itertuples():
            distances = (ts.tss - (p.start + p.end) / 2).abs()
            closest = ts[distances.eq(distances.min())]
            contained = ts[ts.tss.ge(p.start) & ts.tss.lt(p.end)]
            overlaps = g[g.start.lt(p.end) & g.end.gt(p.start)]
            result.append(
                dict(
                    peak=p.peak,
                    nearest_TSS_genes=";".join(sorted(set(closest.gene))),
                    nearest_TSS_distance_bp=distances.min(),
                    contained_TSS_genes=";".join(sorted(set(contained.gene))),
                    contained_transcript_types=";".join(
                        sorted(set(contained.transcript_type))
                    ),
                    overlapping_gene_bodies=";".join(sorted(set(overlaps.gene))),
                )
            )
    pd.DataFrame(result).replace("", pd.NA).to_csv(
        out / "unrestricted_locus_annotation.tsv", sep="\t", index=False, na_rep="NA"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--gtf", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.gtf)
