"""Audit why TSS>=2 and TSS>=3 MYOD1 link discovery nominate 6 vs 5 peaks."""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables", type=Path, required=True)
    a = p.parse_args()
    rows = []
    for gate, suffix in [("TSS_ge_2", "_tss2"), ("TSS_ge_3", "")]:
        links = pd.read_csv(a.tables / f"peak_gene_links{suffix}.tsv", sep="\t")
        by_lib = pd.read_csv(a.tables / f"peak_gene_links_by_library{suffix}.tsv.gz", sep="\t")
        links = links[links.gene == "MYOD1"].copy()
        links["linked_four_libraries"] = (
            (links.q_all_links < 0.05) & (links.partial_r > 0) & (links.n_libraries == 4)
        )
        fourth = by_lib[(by_lib.gene == "MYOD1") & (by_lib.gsm == "GSM6339603")]
        merged = links.merge(
            fourth[["peak", "n_peak_open", "n_gene_detected", "partial_r"]],
            on="peak",
            how="left",
            suffixes=("", "_GSM6339603"),
        )
        merged.insert(0, "link_discovery_gate", gate)
        rows.append(merged)
    out = pd.concat(rows, ignore_index=True)
    out.sort_values(["peak", "link_discovery_gate"]).to_csv(
        a.tables / "myod1_link_gate_reconciliation.tsv", sep="\t", index=False
    )
    check = out[out.peak == "chr11:17652327-17652866"]
    print(
        check[
            [
                "link_discovery_gate",
                "peak",
                "n_libraries",
                "linked_four_libraries",
                "n_peak_open_GSM6339603",
                "partial_r_GSM6339603",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
