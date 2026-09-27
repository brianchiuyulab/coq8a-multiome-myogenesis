"""Cross-reference the frozen functional genes against owner-provided RNA results.

Output must be outside the public repository. This compares gene symbols, not
cross-species enhancer coordinates, and retains the original contrast labels.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main(root, source, out):
    root = root.resolve()
    out = out.resolve()
    if out == root or root in out.parents:
        raise ValueError("Private output must be outside the public repository")
    out.mkdir(parents=True, exist_ok=True)
    genes = set(
        pd.read_csv(root / "results/differentiation_fusion25/genes.tsv", sep="\t").gene
    )
    data = pd.read_csv(source)
    data["human_scope_symbol"] = data.symbol.str.upper()
    data = data[data.human_scope_symbol.isin(genes)].copy()
    data["FC"] = np.exp2(data.log2FC)
    data["p_rounded_to_zero_in_source"] = data.p_value.eq(0)
    data["fdr_rounded_to_zero_in_source"] = data.fdr.eq(0)
    data.to_csv(out / "private_union_RNA_contrasts.tsv", sep="\t", index=False)
    summary = data[data.contrast.str.match(r"P(11|22|33)_D[12]_vs_D0$")]
    summary.to_csv(out / "private_D1_D2_vs_D0.tsv", sep="\t", index=False)
    print(
        summary[
            summary.human_scope_symbol.isin(["CAV3", "CSRP3", "MYOD1", "MYF5", "MYOG"])
        ][["symbol", "contrast", "FC", "p_value", "fdr"]].to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.root, args.source, args.out)
