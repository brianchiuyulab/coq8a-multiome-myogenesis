"""Test every externally defined opening peak against COQ8A within matched nuclei.

The 229-peak family is fixed by external timing, with no target-count filter.
Early, middle and late peaks share one BH family per COQ8A/TSS setting.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


PHASES = {
    "Early": "Early_by24h",
    "Middle": "Middle_24to48h",
    "Late": "Late_after48h",
}


def write(frame, path):
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", compression=compression)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--temporal", type=Path, required=True)
    args = parser.parse_args()
    members = pd.read_csv(args.temporal / "temporal_set_memberships.tsv.gz", sep="\t")
    parts = []
    for phase, suffix in PHASES.items():
        subset = members.loc[members.set_id.eq("P2_O50_" + suffix), ["peak"]].copy()
        subset["phase"] = phase
        parts.append(subset)
    family = pd.concat(parts, ignore_index=True)
    assert not family.peak.duplicated().any()
    assert family.phase.value_counts().to_dict() == {
        "Early": 191,
        "Middle": 34,
        "Late": 4,
    }
    assert set(family.peak) == set(
        members.loc[members.set_id.eq("P2_O50_All_opening"), "peak"]
    )
    counts = pd.read_csv(
        args.tables / "candidate_peak_effects_by_library.tsv.gz", sep="\t"
    )
    counts = counts.merge(family, on="peak", validate="many_to_one")
    assert not counts.duplicated(["gate", "contrast", "gsm", "peak"]).any()
    assert len(counts) == 229 * 4 * 4
    pairs = pd.read_csv(args.tables / "matched_pairs.tsv.gz", sep="\t")
    libraries = pairs[["gsm", "source", "stage"]].drop_duplicates()
    counts = counts.merge(libraries, on="gsm", validate="many_to_one")
    expected_n = pairs.groupby(["gate", "contrast", "gsm"]).size()
    for key, subset in counts.groupby(["gate", "contrast", "gsm"]):
        assert (subset.n_pairs == expected_n.loc[key]).all()
    assert (
        counts.high_open - counts.low_open == counts.high_only - counts.low_only
    ).all()
    assert (counts.high_only + counts.low_only <= counts.n_pairs).all()
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
        p = binomtest(plus, plus + minus, p=0.5).pvalue if plus + minus else 1.0
        rows.append(
            dict(
                gate=gate,
                contrast=contrast,
                peak=peak,
                phase=phase,
                nearby_candidate_genes=genes.loc[peak],
                n_pairs=n,
                n_libraries=len(subset),
                high_open=high,
                low_open=low,
                high_only=plus,
                low_only=minus,
                high_open_pct=100 * high / n,
                low_open_pct=100 * low / n,
                delta_pp=100 * (high - low) / n,
                fold_open=high / low if low else (np.inf if high else np.nan),
                p_exact=p,
                positive_libraries=int((subset.delta_pp > 0).sum()),
                negative_libraries=int((subset.delta_pp < 0).sum()),
                tied_libraries=int((subset.delta_pp == 0).sum()),
            )
        )
    effects = pd.DataFrame(rows)
    summary = []
    for (gate, contrast), subset in effects.groupby(["gate", "contrast"]):
        assert len(subset) == 229
        q = multipletests(subset.p_exact, method="fdr_bh")[1]
        effects.loc[subset.index, "q_229_peaks"] = q
        for phase in PHASES:
            ix = subset.phase.eq(phase)
            summary.append(
                dict(
                    gate=gate,
                    contrast=contrast,
                    phase=phase,
                    n_peaks=int(ix.sum()),
                    nominal_p_lt_005=int((subset.loc[ix, "p_exact"] < 0.05).sum()),
                    nominal_positive=int(
                        (
                            (subset.loc[ix, "p_exact"] < 0.05)
                            & (subset.loc[ix, "delta_pp"] > 0)
                        ).sum()
                    ),
                    q_lt_005=int((q[ix] < 0.05).sum()),
                    q_lt_010=int((q[ix] < 0.10).sum()),
                )
            )
    # Sum single-peak counts back to the independently stored module effects.
    modules = pd.read_csv(args.temporal / "temporal_effect_grid.tsv", sep="\t")
    checked = 0
    for (gate, contrast), subset in effects.groupby(["gate", "contrast"]):
        for phase, suffix in {**PHASES, "All": "All_opening"}.items():
            selected = subset if phase == "All" else subset[subset.phase.eq(phase)]
            module = modules[
                (modules.scope == "programme")
                & (modules.stage_analysis == "pooled")
                & (modules.set_id == "P2_O50_" + suffix)
                & (modules.gate == gate)
                & (modules.contrast == contrast)
            ].iloc[0]
            assert np.isclose(
                selected.high_open.sum() / selected.low_open.sum(),
                module.fold_open,
                rtol=1e-10,
            )
            assert np.isclose(
                100 * selected.high_open.sum() / (len(selected) * module.n_pairs),
                module.high_open_pct,
            )
            checked += 1
    effects = effects.sort_values(["gate", "contrast", "p_exact", "peak"])
    write(effects, args.temporal / "single_peak_effects.tsv")
    write(
        counts.sort_values(["gate", "contrast", "peak", "gsm"]),
        args.temporal / "single_peak_by_library.tsv.gz",
    )
    write(pd.DataFrame(summary), args.temporal / "single_peak_screen_summary.tsv")
    audit = dict(
        status="passed",
        family_size=229,
        settings=4,
        peak_tests=916,
        by_library_counts=len(counts),
        module_reconciliations=checked,
        target_count_filter="none",
        test="two-sided exact paired binomial",
        adjustment="BH across all 229 peaks separately in each setting",
    )
    (args.temporal / "single_peak_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit))
    primary = effects[(effects.gate == "TSS_ge_3") & (effects.contrast == "3plus_vs_1")]
    print(primary.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
