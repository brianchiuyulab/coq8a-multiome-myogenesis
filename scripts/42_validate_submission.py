"""Validate figure sources, matched comparisons and reported statistical families."""

from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / "results/figure_source"
    order = pd.read_csv(source / "functional_54_order.tsv", sep="\t")
    assert len(order) == 54 and order.peak.is_unique and order.region_id.is_unique
    assert order.phase.value_counts().to_dict() == {
        "Closing": 25,
        "Early": 24,
        "Middle": 5,
    }
    fx = pd.read_csv(ROOT / "results/functional/functional_peak_effects.tsv", sep="\t")
    for _, d in fx.groupby(["gate", "contrast"]):
        assert set(d.peak) == set(order.peak)
        p = [
            (
                binomtest(int(r.high_only), int(r.high_only + r.low_only)).pvalue
                if r.high_only + r.low_only
                else 1
            )
            for r in d.itertuples()
        ]
        assert np.allclose(p, d.p_exact)
        assert np.allclose(multipletests(p, method="fdr_bh")[1], d.q_functional_union)
    profiles = pd.read_csv(source / "fragment_profiles.tsv.gz", sep="\t")
    manifest = json.loads((source / "fragment_profile_manifest.json").read_text())
    pairs = pd.read_csv(ROOT / "results/tables/matched_pairs.tsv.gz", sep="\t").query(
        "gate=='TSS_ge_3' and contrast=='3plus_vs_1'"
    )
    identities = "\n".join(
        sorted(
            "|".join(getattr(r, k) for k in ["gsm", "high_barcode", "low_barcode"])
            for r in pairs.itertuples()
        )
    )
    assert (
        hashlib.sha256(identities.encode()).hexdigest()
        == manifest["paired_barcode_sha256"]
    )
    assert (
        hashlib.sha256(
            (source / "fragment_windows.tsv").read_text().encode()
        ).hexdigest()
        == manifest["windows_sha256"]
    )
    counts = profiles.groupby(["window_id", "group", "bin_start"]).n_nuclei.sum()
    assert counts.eq(201).all()
    assert (profiles.insertions >= 0).all()
    assert np.allclose(
        profiles.insertions_per_100_nuclei_100bp,
        profiles.insertions
        / profiles.n_nuclei
        * 100
        * 100
        / (profiles.bin_end - profiles.bin_start),
    )
    stats = pd.read_csv(source / "focal_statistics.tsv", sep="\t")
    assert set(stats.gene) == {"CSRP3", "CAV3", "MYOD1", "CACNA1H"}
    mainfx = fx.query("gate=='TSS_ge_3' and contrast=='3plus_vs_1'")
    joined = stats.merge(
        mainfx, on="peak", suffixes=("", "_check"), validate="one_to_one"
    )
    for col in ["fold_open", "p_exact", "q_functional_union"]:
        assert np.allclose(joined[col], joined[col + "_check"])
    expected = [
        "main/Figure_1_External_definition",
        "main/Figure_2_Functional_accessibility",
        "main/Figure_3_Candidate_loci",
        "supplement/Figure_221_gene_TSS_heatmap",
        "supplement/Figure_S2_Module_sensitivity",
    ]
    for name in expected:
        for ext in ["png", "pdf"]:
            path = ROOT / "figures" / (name + "." + ext)
            assert path.is_file() and path.stat().st_size > 1000
    result = dict(
        status="passed",
        functional_peaks=54,
        single_peak_tests_checked=len(fx),
        main_matched_pairs=201,
        fragment_windows=profiles.window_id.nunique(),
        main_figures=3,
        supplement_figures=2,
    )
    (source / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
