"""Assemble fixed figure sources from complete analysis outputs and GENCODE."""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

FOCAL = {
    "CSRP3": "chr11:19201752-19202603",
    "CAV3": "chr3:8733438-8733968",
    "MYOD1": "chr11:17696914-17697783",
    "CACNA1H": "chr16:1152563-1153401",
}
PHASES = ["Early", "Middle", "Late", "Closing"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--gtf", type=Path, required=True)
    a = parser.parse_args()
    out = a.root / "results/figure_source"
    out.mkdir(parents=True, exist_ok=True)
    temporal = a.root / "results/temporal"
    functional = a.root / "results/functional"
    effects = pd.read_csv(temporal / "bidirectional_peak_effects.tsv", sep="\t")
    main_effects = effects.query("gate=='TSS_ge_3' and contrast=='3plus_vs_1'").copy()
    annotation = pd.read_csv(temporal / "target_temporal_annotations.tsv.gz", sep="\t")
    annotation = annotation.query("promoter_kb==2 and reciprocal_overlap>=0.5")
    full = main_effects.merge(annotation, on="peak", validate="one_to_one")
    full["phase_order"] = full.phase.map({p: i for i, p in enumerate(PHASES)})
    full = full.sort_values(
        ["phase_order", "half_rise_hour", "peak"], na_position="last"
    )
    assert len(full) == 410
    full.to_csv(out / "external_410_regions.tsv", sep="\t", index=False)
    members = pd.read_csv(functional / "functional_peak_memberships.tsv", sep="\t")
    selected = full[full.peak.isin(members.peak)].copy()
    assert len(selected) == 54
    selected["region_id"] = [f"R{i:02d}" for i in range(1, 55)]
    selected.to_csv(out / "functional_54_order.tsv", sep="\t", index=False)
    windows = []
    for r in selected.itertuples():
        chrom, lo, hi = re.split("[:-]", r.peak)
        center = (int(lo) + int(hi)) // 2
        windows.append(
            dict(
                window_id=r.region_id,
                kind="peak",
                peak=r.peak,
                gene=r.nearby_candidate_genes,
                chrom=chrom,
                start=center - 2000,
                end=center + 2000,
                center=center,
                bin_bp=100,
            )
        )
    transcripts, exonrows = [], []
    with gzip.open(a.gtf, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] not in ["transcript", "exon"]:
                continue
            gene = re.search(r'gene_name "([^"]+)"', fields[8])
            if gene is None or gene.group(1) not in FOCAL:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', fields[8]).group(1)
            record = dict(
                gene=gene.group(1),
                transcript=tid,
                chrom=fields[0],
                start=int(fields[3]) - 1,
                end=int(fields[4]),
                strand=fields[6],
            )
            (transcripts if fields[2] == "transcript" else exonrows).append(record)
    tx = pd.DataFrame(transcripts)
    tx["tss"] = np.where(tx.strand.eq("+"), tx.start, tx.end - 1)
    chosen = []
    for gene, peak in FOCAL.items():
        chrom, lo, hi = re.split("[:-]", peak)
        center = (int(lo) + int(hi)) // 2
        sub = tx[(tx.gene == gene) & (tx.chrom == chrom)].copy()
        sub["distance"] = abs(sub.tss - center)
        row = sub.sort_values(["distance", "transcript"]).iloc[0].to_dict()
        row["peak"] = peak
        chosen.append(row)
        lo2 = min(int(lo), int(row["tss"])) - 5000
        hi2 = max(int(hi), int(row["tss"])) + 5000
        windows.append(
            dict(
                window_id=gene,
                kind="locus",
                peak=peak,
                gene=gene,
                chrom=chrom,
                start=lo2,
                end=hi2,
                center=center,
                bin_bp=100,
            )
        )
    chosen = pd.DataFrame(chosen)
    chosen.to_csv(out / "focal_transcripts.tsv", sep="\t", index=False)
    exons = pd.DataFrame(exonrows)
    exons[exons.transcript.isin(chosen.transcript)].to_csv(
        out / "focal_exons.tsv", sep="\t", index=False
    )
    pd.DataFrame(windows).to_csv(out / "fragment_windows.tsv", sep="\t", index=False)
    focal = pd.read_csv(temporal / "dynamic_evidence.tsv.gz", sep="\t")
    focal = focal.query("gate=='TSS_ge_3' and contrast=='3plus_vs_1'")
    focal = focal[[FOCAL.get(r.gene) == r.peak for r in focal.itertuples()]].copy()
    q = pd.read_csv(functional / "functional_peak_effects.tsv", sep="\t")
    q = q.query("gate=='TSS_ge_3' and contrast=='3plus_vs_1'")
    focal.merge(
        q[["peak", "q_functional_union"]], on="peak", validate="one_to_one"
    ).to_csv(out / "focal_statistics.tsv", sep="\t", index=False)
    print("Prepared 410 external regions, 54 functional regions and four focal loci.")


if __name__ == "__main__":
    main()
