"""Report all matched candidates through the same ATAC/RNA evidence steps."""

import argparse
from pathlib import Path

import pandas as pd


def main(root):
    out = root / "results/local_signal_decision"
    out.mkdir(parents=True, exist_ok=True)
    base = root / "results/differentiation_fusion25/native_GSE240061"
    x = pd.read_csv(base / "ATAC_RNA_all_settings.tsv.gz", sep="\t")
    matched = x[x.caliper.ne("none")].copy()
    matched["ATAC_open_nominal"] = matched.FC.gt(1) & matched.PValue.lt(0.05)
    matched["RNA_up_nominal"] = matched.RNA_FC.gt(1) & matched.RNA_p.lt(0.05)
    matched["positive_link_nominal"] = matched.r_equal_donor.gt(0) & matched.link_p.lt(
        0.05
    )
    criteria = {
        "matched_evaluable_peak_gene_settings": pd.Series(True, index=matched.index),
        "ATAC_open_p005": matched.ATAC_open_nominal,
        "ATAC_open_p005_and_RNA_up_p005": matched.ATAC_open_nominal
        & matched.RNA_up_nominal,
        "ATAC_RNA_p005_and_positive_link_direction": matched.ATAC_open_nominal
        & matched.RNA_up_nominal
        & matched.r_equal_donor.gt(0),
        "ATAC_RNA_link_all_p005": matched.ATAC_open_nominal
        & matched.RNA_up_nominal
        & matched.positive_link_nominal,
    }
    counts = []
    for label, mask in criteria.items():
        d = matched[mask]
        counts.append(
            dict(
                step=label,
                setting_peak_gene_rows=len(d),
                distinct_peaks=d.peak.nunique(),
                distinct_genes=d.gene.nunique(),
            )
        )
    pd.DataFrame(counts).to_csv(
        out / "GSE240061_evidence_steps.tsv", sep="\t", index=False
    )
    d = matched[matched.ATAC_open_nominal & matched.RNA_up_nominal].sort_values(
        "PValue"
    )
    annotation = pd.read_csv(base / "unrestricted_locus_annotation.tsv", sep="\t")
    d = d.merge(annotation, on="peak", how="left")
    d.to_csv(
        out / "GSE240061_ATAC_and_RNA_nominal_candidates.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    rows = []
    for (ct, time, peak, gene), d in matched.groupby(
        ["celltype", "time", "peak", "gene"]
    ):
        b = d.loc[d.PValue.idxmin()]
        rows.append(
            dict(
                celltype=ct,
                time=time,
                peak=peak,
                gene=gene,
                n_evaluable_settings=len(d),
                n_open_p005=int(d.ATAC_open_nominal.sum()),
                n_open_q01=int((d.FC.gt(1) & d.q_union.lt(0.1)).sum()),
                n_ATAC_RNA_p005=int((d.ATAC_open_nominal & d.RNA_up_nominal).sum()),
                n_all_three_p005=int(
                    (
                        d.ATAC_open_nominal & d.RNA_up_nominal & d.positive_link_nominal
                    ).sum()
                ),
                best_raw_p_config=b.config,
                best_raw_p_FC=b.FC,
                best_raw_p=b.PValue,
                best_raw_p_q=b.q_union,
                peak_RNA_r=b.r_equal_donor,
                peak_RNA_p=b.link_p,
            )
        )
    pd.DataFrame(rows).to_csv(
        out / "GSE240061_all_candidate_stability.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    myod = matched[
        matched.peak.eq("chr11-17719024-17719942")
        & matched.celltype.eq("Satellite Cells")
        & matched.time.eq("Post")
    ]
    myod.to_csv(
        out / "MYOD1_promoter_all_matched_settings.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    csrp = matched[
        matched.peak.eq("chr11-19201660-19202638")
        & matched.gene.eq("CSRP3")
        & matched.celltype.eq("Fast")
        & matched.time.eq("Pre")
    ]
    csrp.to_csv(
        out / "CSRP3_new_promoter_all_matched_settings.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    print(pd.DataFrame(counts).to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    main(p.parse_args().root)
