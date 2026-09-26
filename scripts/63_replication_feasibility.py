"""Inventory replication cell types and unchanged within-sample depth matching."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main(root, work, out):
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(work / "nucleus_inventory.tsv", sep="\t")
    geo = pd.read_csv(out / "GEO_sample_manifest.tsv", sep="\t")
    parsed = geo.sample_id.str.extract(r"^\d+-([A-Z])-M([12])_")
    geo["Sample"] = parsed[0] + "_M" + parsed[1]
    design = geo[["Sample", "subject_id", "time", "group", "Sex"]].drop_duplicates()
    assert design.Sample.is_unique and len(design) == 12
    meta = meta.merge(
        design, on="Sample", validate="many_to_one", suffixes=("", "_GEO")
    )
    assert meta.Time.eq(meta.time).all() and meta.Sex.eq(meta.Sex_GEO).all()
    assert meta.barcode.is_unique and len(meta) == 37154
    assert np.array_equal(meta.raw_RNA_library_size, meta.nCount_RNA)
    assert np.array_equal(meta.raw_ATAC_library_size, meta.nCount_ATAC)
    meta["celltype"] = meta["refined_annotations_wknn_0.8"]
    meta["COQ8A_umi"] = meta.COQ8A_raw_UMI
    meta["total_rna_umi"] = meta.raw_RNA_library_size
    meta["total_open_peaks"] = meta.nFeature_ATAC
    rows = []
    for keys, d in meta.groupby(["Sample", "subject_id", "Time", "Group", "celltype"]):
        rows.append(
            dict(zip(["Sample", "donor", "time", "group", "celltype"], keys))
            | dict(
                n_nuclei=len(d),
                n_COQ8A_zero=int(d.COQ8A_umi.eq(0).sum()),
                n_COQ8A_eq1=int(d.COQ8A_umi.eq(1).sum()),
                n_COQ8A_ge2=int(d.COQ8A_umi.ge(2).sum()),
                n_COQ8A_ge3=int(d.COQ8A_umi.ge(3).sum()),
                n_CAV3_detected=int(d.CAV3_raw_UMI.gt(0).sum()),
                median_RNA_depth=float(d.total_rna_umi.median()),
                median_open_peaks=float(d.total_open_peaks.median()),
            )
        )
    inventory = pd.DataFrame(rows)
    inventory.to_csv(out / "celltype_donor_group_inventory.tsv", sep="\t", index=False)
    summary = (
        inventory.groupby(["time", "celltype"])[
            [
                "n_nuclei",
                "n_COQ8A_zero",
                "n_COQ8A_eq1",
                "n_COQ8A_ge2",
                "n_COQ8A_ge3",
                "n_CAV3_detected",
            ]
        ]
        .sum()
        .reset_index()
    )
    spec = importlib.util.spec_from_file_location(
        "original_matching", root / "scripts/02_qc_and_matching.py"
    )
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    pairs = []
    match_rows = []
    for (sample, ctype), d in meta.groupby(["Sample", "celltype"]):
        for cut in [2, 3]:
            matched = original.depth_match(d, cut, 1)
            match_rows.append(
                dict(
                    Sample=sample,
                    celltype=ctype,
                    contrast=f"{cut}plus_vs_1",
                    donor=d.subject_id.iloc[0],
                    time=d.Time.iloc[0],
                    n_pairs=len(matched),
                )
            )
            for high, low in matched:
                pairs.append(
                    dict(
                        Sample=sample,
                        celltype=ctype,
                        contrast=f"{cut}plus_vs_1",
                        donor=high.subject_id,
                        time=high.Time,
                        high_barcode=high.barcode,
                        low_barcode=low.barcode,
                        abs_log_RNA=np.abs(
                            np.log1p(high.total_rna_umi) - np.log1p(low.total_rna_umi)
                        ),
                        abs_log_open_peaks=np.abs(
                            np.log1p(high.total_open_peaks)
                            - np.log1p(low.total_open_peaks)
                        ),
                    )
                )
    pd.DataFrame(pairs).to_csv(
        work / "author_QC_depth_matched_pairs.tsv.gz", sep="\t", index=False
    )
    matches = pd.DataFrame(match_rows)
    matches.to_csv(out / "depth_matching_feasibility.tsv", sep="\t", index=False)
    total = (
        matches.groupby(["time", "celltype", "contrast"])
        .agg(
            n_pairs=("n_pairs", "sum"),
            n_donors_with_pairs=("n_pairs", lambda x: int((x > 0).sum())),
        )
        .reset_index()
    )
    total.to_csv(out / "depth_matching_summary.tsv", sep="\t", index=False)
    sensitivity = []
    satellite = meta[meta.celltype.eq("Satellite Cells")]
    for caliper in [0.10, 0.20, 0.30]:
        original.CALIPER = caliper
        for sample, d in satellite.groupby("Sample"):
            for cut in [2, 3]:
                matched = original.depth_match(d, cut, 1)
                errors = [
                    max(
                        abs(np.log1p(h.total_rna_umi) - np.log1p(l.total_rna_umi)),
                        abs(
                            np.log1p(h.total_open_peaks) - np.log1p(l.total_open_peaks)
                        ),
                    )
                    for h, l in matched
                ]
                sensitivity.append(
                    dict(
                        Sample=sample,
                        donor=d.subject_id.iloc[0],
                        time=d.Time.iloc[0],
                        caliper=caliper,
                        contrast=f"{cut}plus_vs_1",
                        n_pairs=len(matched),
                        max_log_depth_mismatch=max(errors) if errors else np.nan,
                    )
                )
    pd.DataFrame(sensitivity).to_csv(
        out / "satellite_matching_sensitivity.tsv", sep="\t", index=False, na_rep="NA"
    )
    summary.to_csv(out / "celltype_group_summary.tsv", sep="\t", index=False)
    overlap = pd.read_csv(work / "candidate_peak_overlap.tsv", sep="\t")
    a = pd.read_csv(root / "results/enhancer_regions/peak_selectivity.tsv", sep="\t")
    f = pd.read_csv(
        root / "results/candidate_provenance/external_fusion_promoter_candidates.tsv",
        sep="\t",
    )
    scope_rows = []
    for name, peaks in [
        ("all5097", None),
        ("HSMM1777", set(a.peak)),
        ("selective401", set(a.loc[a.n_other_strong.eq(0), "peak"])),
        ("fusion_promoter31", set(f.peak)),
    ]:
        sub = overlap if peaks is None else overlap[overlap.target_peak.isin(peaks)]
        strict = sub[
            sub.target_overlap_fraction.ge(0.5) & sub.author_overlap_fraction.ge(0.5)
        ]
        scope_rows.append(
            dict(
                scope=name,
                any_overlap_old_peaks=sub.target_peak.nunique(),
                reciprocal50_old_peaks=strict.target_peak.nunique(),
                reciprocal50_author_peaks=strict.author_peak.nunique(),
            )
        )
    pd.DataFrame(scope_rows).to_csv(
        out / "candidate_coordinate_coverage.tsv", sep="\t", index=False
    )
    overlap[overlap.target_peak.eq("chr3:8733438-8733968")].to_csv(
        out / "CAV3_interval_correspondence.tsv", sep="\t", index=False
    )
    report = dict(
        n_nuclei=len(meta),
        n_donors=meta.subject_id.nunique(),
        n_samples=meta.Sample.nunique(),
        RNA_features=36601,
        ATAC_features=144663,
        author_annotation_column="refined_annotations_wknn_0.8",
        n_celltypes=meta.celltype.nunique(),
        same_nucleus_barcodes_verified=True,
        QC="Author-retained nuclei; deposited metadata lack TSS enrichment, FRiP, blacklist and nucleosome metrics; no new TSS>=3 claim",
        matching="Original log1p RNA UMI and detected-peak caliper 0.10; without replacement, within sample and cell type",
        status="Feasibility and coordinate audit; differential-accessibility p/q not yet calculated",
        matching_scope="No additional old-library depth caps or reconstructed fragment QC imposed at this feasibility stage",
    )
    (out / "feasibility.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        summary[summary.celltype.isin(["Satellite Cells", "Fast", "Slow"])].to_string(
            index=False
        )
    )
    print(
        total[total.celltype.isin(["Satellite Cells", "Fast", "Slow"])].to_string(
            index=False
        )
    )
    print(pd.DataFrame(scope_rows).to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work, a.out)
