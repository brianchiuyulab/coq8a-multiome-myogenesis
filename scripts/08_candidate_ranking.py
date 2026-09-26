"""Rank exploratory loci after the complete primary ATAC test families.

This ranking combines same-data peak–RNA links and the direction of COQ8A
accessibility effects. It nominates candidates for follow-up; it does not
provide an independent statistical test of the nominated locus.
"""

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    args = parser.parse_args()
    tables = args.tables

    links = pd.read_csv(tables / "peak_gene_links.tsv", sep="\t")
    effects = pd.read_csv(tables / "candidate_peak_effects_pooled.tsv.gz", sep="\t")
    search = pd.read_csv(tables / "gene_search_space.tsv", sep="\t")
    gene_effects = pd.read_csv(tables / "gene_effects_pooled.tsv", sep="\t")
    for contrast in ("2plus_vs_1", "3plus_vs_1"):
        peak_effects = effects[
            (effects.gate == "TSS_ge_3") & (effects.contrast == contrast)
        ]
        linked = links.merge(
            peak_effects[["peak", "delta_pp", "positive_libraries"]],
            on="peak",
            how="left",
            validate="many_to_one",
        )
        linked["positive_link"] = (
            (linked.q_all_links < 0.05)
            & (linked.partial_r > 0)
            & (linked.n_libraries == 4)
        )
        linked["positive_link_and_effect"] = (
            linked.positive_link
            & (linked.delta_pp > 0)
            & (linked.positive_libraries >= 3)
        )
        link_summary = (
            linked.groupby("gene")
            .agg(
                n_tested_links=("peak", "size"),
                n_positive_links=("positive_link", "sum"),
                n_linked_coq_3of4=("positive_link_and_effect", "sum"),
                strongest_link_r=("partial_r", "max"),
            )
            .reset_index()
            .rename(columns={"n_positive_links": "n_linked"})
        )
        atac = gene_effects[
            (gene_effects.gate == "TSS_ge_3")
            & (gene_effects.contrast == contrast)
            & (gene_effects.modality == "ATAC")
        ]
        ranking = search.merge(link_summary, on="gene", how="left").merge(
            atac[["gene", "difference", "p_pair", "q_221", "positive_libraries"]],
            on="gene",
            how="left",
            validate="one_to_one",
        )
        count_columns = ["n_tested_links", "n_linked", "n_linked_coq_3of4"]
        ranking[count_columns] = ranking[count_columns].fillna(0).astype(int)
        ranking = ranking.sort_values(
            ["n_linked_coq_3of4", "n_linked", "difference", "gene"],
            ascending=[False, False, False, True],
        )
        ranking = ranking[
            [
                "gene",
                "Hallmark_Myogenesis",
                "Reactome_Myogenesis",
                "named_core_TF",
                "n_candidate_peaks",
                "n_common4",
                "n_tested_links",
                "strongest_link_r",
                "difference",
                "p_pair",
                "q_221",
                "positive_libraries",
                "n_linked",
                "n_linked_coq_3of4",
            ]
        ]
        ranking.to_csv(
            tables / f"candidate_gene_ranking_{contrast}.tsv", sep="\t", index=False
        )
        print(
            f"{contrast}: {ranking.iloc[0].gene} leads the exploratory ranking",
            flush=True,
        )


if __name__ == "__main__":
    main()
