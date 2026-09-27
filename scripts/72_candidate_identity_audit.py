"""Separate biological-set membership from unrestricted locus annotation."""

import argparse
import gzip
import re
from pathlib import Path
import pandas as pd


def main(root, gtf):
    out = root / "results/candidate_identity"
    out.mkdir(parents=True, exist_ok=True)
    universe = pd.read_csv(root / "results/tables/gene_search_space.tsv", sep="\t")
    functions = pd.read_csv(root / "results/functional/functional_genes.tsv", sep="\t")
    near = pd.read_csv(root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
    near = near[near.n_libraries.eq(4)]
    members = pd.read_csv(
        root / "results/cav3_scope_grid/scope_memberships.tsv.gz", sep="\t"
    )
    rows = []
    for gene in ["MYOD1", "MYOG", "MYF5", "CAV3", "CACNA1H"]:
        p = near[near.gene.eq(gene)]
        rows.append(
            dict(
                gene=gene,
                in_221=gene in set(universe.gene),
                functions=";".join(functions.loc[functions.gene.eq(gene), "function"]),
                n_common_peaks=p.peak.nunique(),
                minimum_candidate_TSS_distance_bp=p.nearest_tss_bp.min(),
                fusion_500_peaks=members.loc[
                    members.family.eq("Fusion__TSS500__All") & members.gene.eq(gene),
                    "peak",
                ].nunique(),
                fusion_500_selective_peaks=members.loc[
                    members.family.eq("Fusion__TSS500__Muscle_selective")
                    & members.gene.eq(gene),
                    "peak",
                ].nunique(),
            )
        )
    pd.DataFrame(rows).to_csv(out / "gene_scope_membership.tsv", sep="\t", index=False)
    loci = [("chr16", 1311478, 1312392), ("chr3", 8733222, 8734203)]
    tss_rows, body_rows = [], []
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            v = line.rstrip().split("\t")
            if v[2] not in ["gene", "transcript"] or v[0] not in ["chr16", "chr3"]:
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]+)"', v[8]))
            start, end = int(v[3]) - 1, int(v[4])
            tss = start if v[6] == "+" else end - 1
            for chrom, lo, hi in loci:
                if chrom != v[0]:
                    continue
                peak = f"{chrom}:{lo}-{hi}"
                if v[2] == "gene" and min(hi, end) > max(lo, start):
                    body_rows.append(
                        dict(
                            peak=peak,
                            gene=attributes.get("gene_name"),
                            gene_start=start,
                            gene_end=end,
                            gene_type=attributes.get("gene_type"),
                        )
                    )
                distance = abs((lo + hi) / 2 - tss)
                if v[2] == "transcript" and distance <= 200000:
                    tss_rows.append(
                        dict(
                            peak=peak,
                            gene=attributes.get("gene_name"),
                            transcript=attributes.get("transcript_id"),
                            transcript_type=attributes.get("transcript_type"),
                            tss_0based=tss,
                            midpoint_distance_bp=distance,
                            peak_contains_TSS=lo <= tss < hi,
                        )
                    )
    tss = pd.DataFrame(tss_rows).sort_values(["peak", "midpoint_distance_bp", "gene"])
    tss.to_csv(out / "all_nearby_transcript_TSS.tsv", sep="\t", index=False)
    tss.drop_duplicates(["peak", "gene"]).to_csv(
        out / "nearest_TSS_per_gene.tsv", sep="\t", index=False
    )
    pd.DataFrame(body_rows).to_csv(
        out / "overlapping_gene_bodies.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(rows).to_string(index=False))
    print(tss.drop_duplicates("peak").to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--gtf", type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.gtf)
