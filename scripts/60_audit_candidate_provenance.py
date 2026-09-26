"""Compare all candidate neighborhoods with identical evidence columns."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main(root):
    out = root / "results/candidate_provenance"
    out.mkdir(parents=True, exist_ok=True)
    peaks = pd.read_csv(root / "results/unstratified/ATAC_all_5097.tsv", sep="\t")
    assert len(peaks) == 5097 and peaks.peak.is_unique
    annotations = pd.read_csv(
        root / "results/enhancer_regions/peak_selectivity.tsv", sep="\t"
    )
    peaks = peaks.merge(
        annotations[["peak", "n_other_strong"]],
        on="peak",
        how="left",
        validate="one_to_one",
    )
    peaks["HSMM_strong"] = peaks.peak.isin(annotations.peak)
    peaks["muscle_selective"] = peaks.HSMM_strong & peaks.n_other_strong.eq(0)
    for name, mask in [
        ("all5097", np.ones(len(peaks), dtype=bool)),
        ("HSMM1777", peaks.HSMM_strong),
        ("selective401", peaks.muscle_selective),
    ]:
        peaks.loc[mask, "p_rank_" + name] = peaks.loc[mask, "p_pair_binomial"].rank(
            method="min"
        )
        peaks.loc[mask, "q_" + name] = multipletests(
            peaks.loc[mask, "p_pair_binomial"], method="fdr_bh"
        )[1]
    source = pd.read_csv(
        root / "results/candidate_comparison/peak_effects_by_source.tsv.gz", sep="\t"
    )
    wide = source.pivot(
        index="peak", columns="source", values="FC_equal_stage"
    ).add_prefix("source_FC_")
    peaks = peaks.merge(wide, on="peak", validate="one_to_one")
    neighbors = pd.read_csv(
        root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    neighbors = neighbors[neighbors.n_libraries.eq(4)].drop(columns="n_libraries")
    evidence = neighbors.merge(peaks, on="peak", validate="many_to_one")
    links = pd.read_csv(
        root / "results/unstratified/correlations_combined.tsv.gz", sep="\t"
    )
    links = links[
        links.model.eq("technical_state_coq") & links.detection_fraction.eq(0.05)
    ]
    links = links[
        [
            "peak",
            "gene",
            "n_libraries",
            "r",
            "p",
            "q",
            "positive_libraries",
            "loo_r_min",
            "loo_r_max",
        ]
    ]
    links = links.rename(
        columns={c: "link_" + c for c in links.columns if c not in ["peak", "gene"]}
    )
    evidence = evidence.merge(
        links, on=["peak", "gene"], how="left", validate="one_to_one"
    )
    rna = pd.read_csv(root / "results/unstratified/RNA_COQ8A_high_low.tsv", sep="\t")
    rna = rna[["gene", "difference", "p_pair", "q_measured_targets"]].rename(
        columns={
            "difference": "RNA_high_low_delta",
            "p_pair": "RNA_high_low_p",
            "q_measured_targets": "RNA_high_low_q",
        }
    )
    evidence = evidence.merge(rna, on="gene", how="left", validate="many_to_one")
    time = pd.read_csv(
        root / "results/external_c2c12/all_gene_time_contrasts.tsv.gz", sep="\t"
    )
    time["gene_human_symbol"] = time.gene.str.upper()
    unique = time.groupby("gene_human_symbol").ID.nunique()
    time = time[time.gene_human_symbol.isin(unique[unique.eq(1)].index)]
    for hour in [24, 48, 96]:
        sub = time[time.time_h.eq(hour)][["gene_human_symbol", "FC", "padj"]].rename(
            columns={
                "gene_human_symbol": "gene",
                "FC": f"external_RNA_FC_{hour}h",
                "padj": f"external_RNA_q_{hour}h",
            }
        )
        evidence = evidence.merge(sub, on="gene", how="left", validate="many_to_one")
    functions = pd.read_csv(root / "results/functional/functional_genes.tsv", sep="\t")
    labels = functions.groupby("gene")["function"].agg(
        lambda x: ";".join(sorted(set(x)))
    )
    evidence["external_GO_functions"] = evidence.gene.map(labels)
    fusion_genes = set(functions.loc[functions["function"].eq("Fusion"), "gene"])
    fusion_promoters = evidence[
        evidence.gene.isin(fusion_genes) & evidence.nearest_tss_bp.le(2000)
    ]
    fusion_promoters[["peak", "gene", "nearest_tss_bp"]].to_csv(
        out / "external_fusion_promoter_candidates.tsv", sep="\t", index=False
    )
    evidence["prior_count_followup"] = False
    follow = pd.read_csv(
        root / "results/unstratified/count_followup_selection.tsv", sep="\t"
    )
    evidence = evidence.merge(
        follow[["peak", "gene", "followup_reason"]],
        on=["peak", "gene"],
        how="left",
        validate="one_to_one",
    )
    evidence["prior_count_followup"] = evidence.followup_reason.notna()
    evidence["nominal_positive_ATAC"] = evidence.fold_open.gt(
        1
    ) & evidence.p_pair_binomial.lt(0.05)
    evidence["positive_full_family_RNA_link"] = evidence.link_r.gt(
        0
    ) & evidence.link_q.lt(0.1)
    evidence["positive_COQ_RNA_nominal"] = evidence.RNA_high_low_delta.gt(
        0
    ) & evidence.RNA_high_low_p.lt(0.05)
    evidence["external_RNA_48h_up"] = evidence.external_RNA_FC_48h.gt(
        1
    ) & evidence.external_RNA_q_48h.lt(0.1)
    evidence.sort_values(["p_pair_binomial", "peak", "gene"]).to_csv(
        out / "all_candidate_evidence.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    rows = []
    for scope, mask in [
        ("all5097", evidence.peak.notna()),
        ("HSMM1777", evidence.HSMM_strong),
        ("selective401", evidence.muscle_selective),
    ]:
        sub = evidence[mask].copy()
        for step in [
            "all",
            "nominal_positive_ATAC",
            "positive_full_family_RNA_link",
            "positive_COQ_RNA_nominal",
            "external_RNA_48h_up",
        ]:
            if step != "all":
                sub = sub[sub[step]]
            rows.append(
                dict(
                    scope=scope,
                    cumulative_step=step,
                    n_peaks=sub.peak.nunique(),
                    n_pairs=len(sub),
                    n_genes=sub.gene.nunique(),
                    genes=";".join(sorted(sub.gene.unique())) or "NA",
                )
            )
        sub.to_csv(
            out / (scope + "_nominal_chain.tsv"), sep="\t", index=False, na_rep="NA"
        )
    pd.DataFrame(rows).to_csv(
        out / "descriptive_gate_counts.tsv", sep="\t", index=False
    )
    focal = evidence[evidence.gene.eq("CAV3")].sort_values("p_pair_binomial")
    focal.to_csv(out / "CAV3_all_neighbors.tsv", sep="\t", index=False, na_rep="NA")
    manifest = dict(
        status="retrospective provenance audit, not a new confirmatory screen",
        fixed_contrast="TSS>=3; COQ8A>=3 versus 1 UMI; 201 pairs",
        source_of_focal_choice="Prior focal gene list in script 50; lowest ATAC p among eligible original CAV3 pairs; displayed explicitly in scripts 51 and 58",
        gates="Descriptive cumulative intersections; no combined p value and no new FDR guarantee",
        external_RNA_time="48h retained from prior reference comparison; 24h and 96h also supplied",
        independent_replication="GSE240061 outcomes not inspected before this audit",
    )
    manifest["independent_replication_secondary_scope"] = dict(
        definition="GO myoblast fusion genes within the existing 221-gene universe; peak midpoint within 2 kb of an annotated transcript TSS; no target-effect filter",
        genes=sorted(fusion_genes),
        n_genes=len(fusion_genes),
        n_peaks=int(fusion_promoters.peak.nunique()),
        use="Secondary function-based replication family; neither a redefinition of the original discovery scope nor evidence that every promoter responds to COQ8A",
    )
    (out / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(rows).drop(columns="genes").to_string(index=False))
    print(
        focal[
            [
                "peak",
                "p_pair_binomial",
                "p_rank_HSMM1777",
                "p_rank_selective401",
                "followup_reason",
            ]
        ]
        .head(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    main(parser.parse_args().root)
