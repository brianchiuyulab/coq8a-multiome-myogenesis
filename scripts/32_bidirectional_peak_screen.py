"""Evaluate all external opening and closing peaks in fixed COQ8A comparisons.

The new family contains 410 regions: 191 early, 34 middle, four late-opening,
and 181 closing peaks. Original 229-peak results are retained and reconciled.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--temporal", type=Path, required=True)
    args = parser.parse_args()
    members = pd.read_csv(args.temporal / "temporal_set_memberships.tsv.gz", sep="\t")
    phases = {
        "Early": "Early_by24h",
        "Middle": "Middle_24to48h",
        "Late": "Late_after48h",
        "Closing": "Closing",
    }
    family = pd.concat(
        [
            members.loc[members.set_id.eq("P2_O50_" + suffix), ["peak"]].assign(
                phase=phase
            )
            for phase, suffix in phases.items()
        ],
        ignore_index=True,
    )
    assert family.peak.is_unique and len(family) == 410
    counts = pd.read_csv(
        args.tables / "candidate_peak_effects_by_library.tsv.gz", sep="\t"
    )
    counts = counts.merge(family, on="peak", validate="many_to_one")
    assert len(counts) == 410 * 4 * 4
    genes = pd.read_csv(args.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    genes = genes.groupby("peak").gene.agg(lambda x: ";".join(sorted(set(x))))
    rows = []
    for (gate, contrast, peak, phase), subset in counts.groupby(
        ["gate", "contrast", "peak", "phase"]
    ):
        n, high, low, plus, minus = [
            int(subset[c].sum())
            for c in ["n_pairs", "high_open", "low_open", "high_only", "low_only"]
        ]
        assert high - low == plus - minus
        rows.append(
            dict(
                gate=gate,
                contrast=contrast,
                peak=peak,
                phase=phase,
                nearby_candidate_genes=genes.loc[peak],
                n_pairs=n,
                high_open=high,
                low_open=low,
                high_only=plus,
                low_only=minus,
                high_open_pct=100 * high / n,
                low_open_pct=100 * low / n,
                delta_pp=100 * (high - low) / n,
                fold_open=high / low if low else (np.inf if high else np.nan),
                p_exact=binomtest(plus, plus + minus).pvalue if plus + minus else 1.0,
                positive_libraries=int((subset.delta_pp > 0).sum()),
                negative_libraries=int((subset.delta_pp < 0).sum()),
                tied_libraries=int((subset.delta_pp == 0).sum()),
            )
        )
    effects = pd.DataFrame(rows)
    for _, subset in effects.groupby(["gate", "contrast"]):
        assert len(subset) == 410
        effects.loc[subset.index, "q_410_peaks"] = multipletests(
            subset.p_exact, method="fdr_bh"
        )[1]
    for _, subset in effects.groupby(["gate", "contrast", "phase"]):
        effects.loc[subset.index, "q_within_phase_peaks"] = multipletests(
            subset.p_exact, method="fdr_bh"
        )[1]
    previous = pd.read_csv(args.temporal / "single_peak_effects.tsv", sep="\t")
    check = effects.merge(
        previous, on=["gate", "contrast", "peak"], suffixes=("", "_old")
    )
    assert len(check) == 229 * 4
    for column in ["high_open", "low_open", "p_exact", "fold_open"]:
        assert np.allclose(check[column], check[column + "_old"], equal_nan=True)
    effects = effects.merge(
        previous[["gate", "contrast", "peak", "q_229_peaks"]],
        on=["gate", "contrast", "peak"],
        how="left",
        validate="one_to_one",
    )
    modules = pd.read_csv(args.temporal / "temporal_effect_grid.tsv", sep="\t")
    for (gate, contrast, phase), subset in effects.groupby(
        ["gate", "contrast", "phase"]
    ):
        module = modules.loc[
            (modules.set_id == "P2_O50_" + phases[phase])
            & (modules.gate == gate)
            & (modules.contrast == contrast)
            & (modules.stage_analysis == "pooled")
        ].iloc[0]
        assert np.isclose(
            subset.high_open.sum() / subset.low_open.sum(), module.fold_open
        )
    effects.sort_values(["gate", "contrast", "p_exact", "peak"]).to_csv(
        args.temporal / "bidirectional_peak_effects.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    audit = dict(
        peaks=410,
        comparisons=4,
        opening_tests_reconciled=len(check),
        module_reconciliations=16,
        test="two-sided exact paired binomial",
        correction="BH within each comparison; 410 total and phase-specific families reported",
    )
    (args.temporal / "bidirectional_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    primary = effects.loc[
        (effects.gate == "TSS_ge_3") & (effects.contrast == "3plus_vs_1")
    ]
    print(
        primary.loc[primary.phase.eq("Closing")]
        .sort_values("p_exact")
        .head(15)
        .to_string(index=False)
    )
    print(json.dumps(audit))


if __name__ == "__main__":
    main()
