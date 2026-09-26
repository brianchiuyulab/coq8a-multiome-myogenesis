"""Join public candidate neighborhoods to an existing private RNA result table.

No RNA model is refit. Contrasts and original genome-wide FDRs are retained.
Use a private output directory: these combined results must not be published
unless the owner has separately authorized release of the RNA experiment.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--rna-contrasts", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(args.evidence, sep="\t")
    rna = pd.read_csv(args.rna_contrasts)
    rna["gene"] = rna.symbol.str.upper()
    rna = rna.loc[rna.gene.isin(evidence.gene)].copy()
    assert not rna.duplicated(["gene", "contrast"]).any()
    rna["mapping_method"] = "case-insensitive symbol match; not peak orthology"
    rna.to_csv(
        args.out / "dynamic_existing_RNA_contrasts.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    selected = [
        "P11_D1_vs_D0",
        "P11_D2_vs_D0",
        "P11_D6_vs_D0",
        "P22_vs_P11_at_D2",
        "P33_vs_P11_at_D2",
        "P33_vs_P11_at_D6",
    ]
    assert set(selected).issubset(set(rna.contrast))
    for contrast in selected:
        sub = rna.loc[rna.contrast.eq(contrast), ["gene", "log2FC", "fdr"]].copy()
        sub["FC"] = np.exp2(sub.log2FC)
        sub = sub.rename(
            columns={c: contrast + "_" + c for c in ["log2FC", "fdr", "FC"]}
        )
        evidence = evidence.merge(sub, on="gene", how="left", validate="many_to_one")
    evidence.to_csv(
        args.out / "dynamic_candidate_evidence.tsv", sep="\t", index=False, na_rep="NA"
    )
    print(
        "Candidate peak-gene rows:",
        len(evidence),
        "; RNA genes with contrasts:",
        rna.gene.nunique(),
    )


if __name__ == "__main__":
    main()
