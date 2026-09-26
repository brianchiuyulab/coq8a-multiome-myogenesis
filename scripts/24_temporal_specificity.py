"""Direct module contrasts and cell-state interactions for temporal effects."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from multiome_core import paired_stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--temporal", type=Path, required=True)
    a = p.parse_args()
    scores = pd.read_csv(a.temporal / "matched_pair_temporal_scores.tsv.gz", sep="\t")
    definitions = pd.read_csv(a.temporal / "tested_set_definitions.tsv", sep="\t")
    programmes = set(definitions.loc[definitions.scope == "programme", "set_id"])
    scores = scores[scores.set_id.isin(programmes)]
    interactions = []
    for key, frame in scores.groupby(["set_id", "gate", "contrast"]):
        y = (frame.high_open - frame.low_open).to_numpy()
        design = pd.DataFrame(
            {
                "intercept": np.ones(len(frame)),
                "undifferentiated": (frame.stage == "stem").to_numpy(dtype=float),
                "source_line2": (frame.source == "line2").to_numpy(dtype=float),
            }
        )
        fit = sm.OLS(y, design).fit(cov_type="HC3")
        ci = fit.conf_int().loc["undifferentiated"]
        interactions.append(
            dict(
                set_id=key[0],
                gate=key[1],
                contrast=key[2],
                stem_minus_diff_delta_pp=100 * fit.params.undifferentiated,
                ci_low_pp=100 * ci.iloc[0],
                ci_high_pp=100 * ci.iloc[1],
                p_interaction=fit.pvalues.undifferentiated,
                n_pairs=len(frame),
            )
        )
    interactions = pd.DataFrame(interactions).merge(definitions, on="set_id")
    interactions["q_within_family"] = np.nan
    for _, ix in interactions.groupby(
        ["gate", "contrast", "promoter_kb", "reciprocal_overlap"]
    ).groups.items():
        interactions.loc[ix, "q_within_family"] = multipletests(
            interactions.loc[ix, "p_interaction"], method="fdr_bh"
        )[1]
    interactions.to_csv(a.temporal / "state_interactions.tsv", sep="\t", index=False)
    comparisons = []
    for stage in ["pooled", "stem", "differentiated"]:
        data = scores if stage == "pooled" else scores[scores.stage == stage]
        for (gate, contrast), frame in data.groupby(["gate", "contrast"]):
            for comparator in ["Early_by24h", "Low_change_profile"]:
                first = frame[frame.set_id == "P2_O50_Middle_24to48h"]
                second = frame[frame.set_id == "P2_O50_" + comparator]
                joint = first.merge(
                    second, on=["gsm", "pair", "source"], suffixes=("_first", "_second")
                )
                values = joint[
                    [
                        "high_open_first",
                        "low_open_first",
                        "high_open_second",
                        "low_open_second",
                    ]
                ].to_numpy()
                h1, l1, h2, l2 = values.mean(axis=0)
                result = paired_stats(
                    (values[:, 0] - values[:, 1]) / l1
                    - (values[:, 2] - values[:, 3]) / l2
                )
                rng = np.random.default_rng(208248)
                boot = []
                groups = [
                    np.flatnonzero(joint.gsm.to_numpy() == gsm)
                    for gsm in joint.gsm.unique()
                ]
                for _ in range(2000):
                    idx = np.concatenate(
                        [rng.choice(g, len(g), replace=True) for g in groups]
                    )
                    bh1, bl1, bh2, bl2 = values[idx].mean(axis=0)
                    if min(bl1, bl2) > 0:
                        boot.append(bh1 / bl1 - bh2 / bl2)
                comparisons.append(
                    dict(
                        gate=gate,
                        contrast=contrast,
                        stage_analysis=stage,
                        comparator=comparator,
                        n_pairs=len(joint),
                        middle_minus_comparator_fold=h1 / l1 - h2 / l2,
                        bootstrap_ci_low=np.quantile(boot, 0.025),
                        bootstrap_ci_high=np.quantile(boot, 0.975),
                        p_scaled_difference=result["p_pair"],
                    )
                )
    pd.DataFrame(comparisons).to_csv(
        a.temporal / "module_specificity.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
