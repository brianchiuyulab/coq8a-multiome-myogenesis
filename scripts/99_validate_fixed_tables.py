"""Validate published statistics, fragment normalization, and table integrity."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/fixed_day0_D206"


def read(name):
    return pd.read_csv(SOURCE / name, sep="\t")


def main():
    peaks, regions, pairs = (
        read("all_peaks.tsv"),
        read("all_regions.tsv"),
        read("matched_pairs.tsv"),
    )
    assert len(peaks) == 565 and peaks.peak.is_unique
    assert len(regions) == 50 and not regions.duplicated(["gene", "mode"]).any()
    assert pairs.groupby("gsm").size().to_dict() == {
        "GSM6339597": 226,
        "GSM6339601": 301,
    }
    for _, part in pairs.groupby("gsm"):
        assert len(set(part.high_barcode) | set(part.low_barcode)) == 2 * len(part)
    for _, part in regions.groupby("mode"):
        assert np.allclose((part.exceedances + 1) / (part.permutations + 1), part.p)
        assert np.allclose(multipletests(part.p, method="fdr_bh")[1], part.q25)
    p = np.array(
        [
            stats.binomtest(int(a), int(a + b)).pvalue if a + b else 1
            for a, b in zip(peaks.high_only, peaks.low_only)
        ]
    )
    assert np.allclose(p, peaks.p)
    assert np.allclose(multipletests(p, method="fdr_bh")[1], peaks.q565)
    assert np.allclose(
        peaks.loc[peaks.low.gt(0), "FC"], (peaks.high / peaks.low)[peaks.low.gt(0)]
    )
    assert (peaks.high == peaks.source1_high + peaks.source2_high).all()
    assert (peaks.low == peaks.source1_low + peaks.source2_low).all()
    selected = read("selected_regions.tsv")
    assert (
        (selected.q25 < 0.05).sum() == 6
        and len(selected) == 10
        and selected.top_peak.nunique() == 9
    )
    by_source = (
        read("selected_ATAC_by_source.tsv").groupby("peak")[["high", "low"]].sum()
    )
    assert np.array_equal(
        by_source, peaks.set_index("peak").loc[by_source.index, ["high", "low"]]
    )
    candidates = read("local_RNA_candidates.tsv")
    assert (
        len(candidates) == 170
        and candidates.gene.nunique() == 170
        and candidates.distance_bp.le(500000).all()
    )
    links = read("links_combined.tsv")
    for _, part in links.groupby(["population", "model"]):
        assert len(part) == 170
        expected = multipletests(part.p.fillna(1), method="fdr_bh")[1]
        assert np.allclose(
            expected[part.p.notna()], part.loc[part.p.notna(), "q_all_local_links"]
        )
    profiles, windows = read("fragment_profiles.tsv.gz"), read("profile_windows.tsv")
    assert (
        len(windows) == 789
        and not profiles.duplicated(["id", "kind", "gsm", "group", "bin"]).any()
    )
    merged = profiles.merge(windows, on=["id", "kind"], validate="many_to_one")
    n_bins = np.ceil((merged.end - merged.start) / merged.bin_bp).astype(int)
    original_bin = np.where(merged.strand.eq("-"), n_bins - 1 - merged.bin, merged.bin)
    width = np.minimum(
        merged.bin_bp, merged.end - merged.start - original_bin * merged.bin_bp
    )
    denominator = merged.gsm.map(pairs.groupby("gsm").size())
    expected = merged.insertions / denominator * 100 * 100 / width
    assert np.allclose(expected, merged.per100_nuclei_per100bp)
    for gsm in ["GSM6339597", "GSM6339601"]:
        folder = ROOT / "reference/fragment_qc"
        tss = pd.read_csv(folder / f"{gsm}_fragment_tss_qc.tsv", sep="\t")
        flank = tss.flank_insertions / 200
        expected = tss.center_insertions / 1001 / flank.where(flank.gt(0), flank.mean())
        assert np.allclose(expected, tss.tss_enrichment)
        nuc = pd.read_csv(folder / f"{gsm}_fragment_nucleosome_qc.tsv", sep="\t")
        expected = nuc.mononucleosomal / nuc.nucleosome_free.replace(0, np.nan)
        assert np.allclose(expected, nuc.nucleosome_signal, equal_nan=True)
    manifest = [
        {
            "file": p.relative_to(ROOT).as_posix(),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for p in sorted(SOURCE.glob("*.tsv*"))
    ]
    audit = SOURCE / "audit"
    audit.mkdir(exist_ok=True)
    pd.DataFrame(manifest).to_csv(
        audit / "source_table_sha256.tsv", sep="\t", index=False
    )
    report = dict(
        passed=True,
        peaks=565,
        regions_per_direction=25,
        pairs=527,
        sources=2,
        local_RNA_pairs=170,
        fragment_profile_rows=len(profiles),
        TSS_and_nucleosome_formulas=True,
        fragment_normalization=True,
    )
    (audit / "published_table_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
