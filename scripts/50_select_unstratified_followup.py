"""Select and document count-model follow-ups after the complete cis screen."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
p.add_argument("--work", type=Path, required=True)
a = p.parse_args()
links = pd.read_csv(a.work / "correlations_combined.tsv.gz", sep="\t")
main = links[(links.model == "technical_coq") & (links.detection_fraction == 0.05)]
effects = pd.read_csv(a.work / "ATAC_all_5097.tsv", sep="\t")
original = pd.read_csv(a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
focal_genes = ["MYOD1", "CSRP3", "CAV3", "CACNA1H"]
focal_pairs = original[(original.n_libraries == 4) & original.gene.isin(focal_genes)][
    ["peak", "gene"]
]
focal_pairs["prior_pair"] = True
focal_peaks = set(focal_pairs.peak)
joined = main.merge(
    effects[["peak", "fold_open", "p_pair_binomial", "positive_libraries"]],
    on="peak",
    suffixes=("_link", "_ATAC"),
)
joined = joined.merge(focal_pairs, on=["peak", "gene"], how="left")
joined["prior_pair"] = joined.prior_pair.eq(True)
selected = joined[
    (
        (joined.p < 0.05)
        & ((joined.p_pair_binomial < 0.05) | joined.peak.isin(focal_peaks))
    )
    | joined.prior_pair
].copy()
selected["followup_reason"] = np.select(
    [selected.prior_pair, selected.p_pair_binomial < 0.05],
    ["prior_focal_gene_pair", "nominal_ATAC_and_link"],
    default="prior_focal_locus_nominal_link",
)
selected.to_csv(
    a.work / "count_followup_selection.tsv", sep="\t", index=False, na_rep="NA"
)
candidates = pd.read_csv(a.work / "candidate_pairs.tsv", sep="\t")
candidates.merge(selected[["peak", "gene"]], on=["peak", "gene"]).to_csv(
    a.work / "count_candidate_pairs.tsv", sep="\t", index=False
)
targets = (
    selected[selected.prior_pair & selected.gene.isin(["MYOD1", "CSRP3", "CAV3"])]
    .sort_values(["p_pair_binomial", "peak"])
    .drop_duplicates("gene")
)
targets[["gene", "peak"]].to_csv(
    a.work / "bootstrap_targets.tsv", sep="\t", index=False
)
validation = json.loads((a.work / "validation.json").read_text())
validation["count_followup_pairs"] = len(selected)
validation["bootstrap_targets"] = targets[["gene", "peak"]].to_dict(orient="records")
(a.work / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
print("Count follow-up pairs:", len(selected))
print(
    targets[["gene", "peak", "fold_open", "p_pair_binomial", "r", "p"]].to_string(
        index=False
    )
)

# The previously discussed CSRP3 interval has sparse RNA and is evaluated
# separately at the already specified 1% detection sensitivity.
sensitivity = a.work / "detection_1pct"
sensitivity.mkdir(exist_ok=True)
csrp3 = candidates[
    (candidates.gene == "CSRP3") & (candidates.peak == "chr11:19201752-19202603")
]
assert len(csrp3) == 1
csrp3.to_csv(sensitivity / "count_candidate_pairs.tsv", sep="\t", index=False)
csrp3[["gene", "peak"]].to_csv(
    sensitivity / "bootstrap_targets.tsv", sep="\t", index=False
)
