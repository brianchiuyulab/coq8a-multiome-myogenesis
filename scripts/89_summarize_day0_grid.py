"""Audit the sensitivity grid, refine subgroup maxima, and report all settings."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(a):
    root = a.root
    out = root / "results/day0_sensitivity_grid"
    spec = importlib.util.spec_from_file_location(
        "grid", root / "scripts/88_day0_sensitivity_grid.py"
    )
    grid = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grid)
    reg = pd.read_csv(out / "regional_screen.tsv.gz", sep="\t")
    refined = pd.read_csv(out / "regional_refined.tsv", sep="\t")
    peaks = pd.read_csv(out / "peak_screen.tsv.gz", sep="\t")
    pairs = pd.read_csv(out / "matched_pairs.tsv.gz", sep="\t")
    configs = pd.read_csv(out / "configurations.tsv", sep="\t")
    status = pd.read_csv(out / "status.tsv", sep="\t")
    member = pd.read_csv(
        root / "results/task1_stage_specific/peak_membership.tsv", sep="\t"
    )
    names = sorted(member.peak.unique())
    genes = sorted(member.gene.unique())
    groups = {g: pd.Index(names).get_indexer(d.peak) for g, d in member.groupby("gene")}
    # Additional precision for the best setting in each population at TSS >=3.
    selected = reg[reg.TSS.eq(3) & reg.both_sources_same_direction].sort_values(
        ["q25", "p", "config", "gene"]
    )
    extra = sorted(
        set(selected.drop_duplicates(["population", "mode"]).config)
        - set(refined.config)
    )
    (out / "population_precision_refinement.json").write_text(
        json.dumps(extra, indent=2) + "\n"
    )
    if extra:
        data = {
            gsm: (
                pd.read_csv(a.work / f"{gsm}_meta.tsv", sep="\t").set_index("barcode"),
                sparse.load_npz(a.work / f"{gsm}_atac.npz").toarray().T,
            )
            for gsm in grid.SAMPLES
        }
        rows = []
        for config in extra:
            ds = []
            for gsm, (meta, atac) in data.items():
                d = pairs[pairs.config.eq(config) & pairs.gsm.eq(gsm)]
                hi = meta.index.get_indexer(d.high_barcode)
                lo = meta.index.get_indexer(d.low_barcode)
                assert min(hi) >= 0 and min(lo) >= 0
                ds.append(atac[hi] - atac[lo])
            rr = grid.permutation(
                np.vstack(ds), groups, genes, 2000000, 20260928 + int(config[1:])
            )
            c = configs[configs.config.eq(config)].iloc[0].to_dict()
            for row in rr:
                j = row.pop("top_index")
                row.update(c)
                row["top_peak"] = names[j]
                rows.append(row)
            print("subgroup precision", config, flush=True)
        refined = pd.concat([refined, pd.DataFrame(rows)], ignore_index=True)
        grid.save(refined, out / "regional_refined.tsv")
    effective = reg.copy()
    for _, row in refined.iterrows():
        mask = (
            effective.config.eq(row.config)
            & effective["mode"].eq(row["mode"])
            & effective.gene.eq(row.gene)
        )
        assert mask.sum() == 1
        for col in [
            "p",
            "q25",
            "exceedances",
            "permutations",
            "seed",
            "MC_lower",
            "MC_upper",
        ]:
            effective.loc[mask, col] = row[col]
    grid.save(effective, out / "regional_results.tsv.gz")
    best = effective[effective.both_sources_same_direction].sort_values(
        ["q25", "p", "config", "gene"]
    )
    grid.save(
        best.drop_duplicates(["population", "mode", "gene"]),
        out / "best_per_population_gene.tsv",
    )

    # Check all exact peak tests and both BH families, not only the winners.
    for cfg, d in peaks.groupby("config"):
        assert len(d) == 565
        pv = [
            stats.binomtest(int(h), int(h + l)).pvalue if h + l else 1.0
            for h, l in zip(d.high_only, d.low_only)
        ]
        assert np.allclose(pv, d.p)
        assert np.allclose(multipletests(pv, method="fdr_bh")[1], d.q565)
    for _, d in effective.groupby(["config", "mode"]):
        assert len(d) == 25
        assert np.allclose((d.exceedances + 1) / (d.permutations + 1), d.p)
        assert np.allclose(multipletests(d.p, method="fdr_bh")[1], d.q25)
    old = pd.read_csv(
        root / "results/task1_stage_specific/matched_pairs.tsv.gz", sep="\t"
    )
    old = old[old.stage.eq("stem") & old.contrast.eq("2plus_vs_1")]
    base = pairs[pairs.config.eq("D170")]
    key = ["gsm", "high_barcode", "low_barcode"]
    assert set(map(tuple, old[key].to_numpy())) == set(map(tuple, base[key].to_numpy()))
    for cfg, d in pairs.groupby("config"):
        for gsm, g in d.groupby("gsm"):
            assert g.high_barcode.is_unique and g.low_barcode.is_unique
            assert not (set(g.high_barcode) & set(g.low_barcode))
    audit = dict(
        configurations=len(configs),
        evaluated=int(status.status.eq("evaluated").sum()),
        unavailable=int(status.status.ne("evaluated").sum()),
        baseline_pairs_reproduced=True,
        all_exact_p_and_BH_reproduced=True,
        regional_BH_reproduced=True,
        refined_configs=int(refined.config.nunique()),
    )
    (out / "validation.json").write_text(json.dumps(audit, indent=2) + "\n")
    counts = []
    for (pop, mode), d in effective.groupby(["population", "mode"]):
        z = d[d.q25.lt(0.1) & d.both_sources_same_direction]
        counts.append(
            dict(
                population=pop,
                mode=mode,
                n_evaluable_configurations=d.config.nunique(),
                n_configs_q01=z.config.nunique(),
                n_genes_q01=z.gene.nunique(),
            )
        )
    grid.save(pd.DataFrame(counts), out / "signal_inventory.tsv")

    figdir = root / "figures/day0_sensitivity_grid"
    figdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42}
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), layout="constrained")
    panels = [
        (
            "ADAM12",
            "opening",
            "myogenic_detected",
            "3plus_vs_1",
            "Myogenic-marker positive | >=3 vs 1",
        ),
        (
            "MYOD1",
            "opening",
            "all_culture",
            "2plus_vs_zero",
            "All D0 cultures | >=2 vs 0",
        ),
        (
            "NOTCH1",
            "closing",
            "all_culture",
            "2plus_vs_zero",
            "All D0 cultures | >=2 vs 0",
        ),
        (
            "MYOG",
            "closing",
            "PAX7_detected",
            "2plus_vs_le1",
            "PAX7 detected | >=2 vs <=1",
        ),
    ]
    for ax, (gene, mode, pop, contrast, subtitle) in zip(axes.flat, panels):
        d = effective[
            effective.TSS.eq(3)
            & effective.gene.eq(gene)
            & effective["mode"].eq(mode)
            & effective.population.eq(pop)
            & effective.contrast.eq(contrast)
        ]
        values = np.zeros((2, 4))
        qs = np.ones((2, 4))
        for i, mt in enumerate(["depth", "depth_state"]):
            for j, cal in enumerate([0.05, 0.1, 0.2, 0.3]):
                row = d[d.matching.eq(mt) & d.caliper.eq(cal)].iloc[0]
                values[i, j] = np.log2(row.top_peak_FC)
                qs[i, j] = row.q25
        im = ax.imshow(values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
        for i in range(2):
            for j in range(4):
                ax.text(
                    j,
                    i,
                    f"FC {2**values[i,j]:.2f}\nq {qs[i,j]:.3f}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    color="white" if abs(values[i, j]) > 1.8 else "black",
                )
        ax.set_xticks(range(4), ["0.05", "0.10", "0.20", "0.30"])
        ax.set_yticks([0, 1], ["Depth", "Depth + state"])
        ax.set_xlabel("Log-depth matching caliper")
        ax.set_title(f"{gene} neighborhood: {mode}\n{subtitle}", fontsize=11, pad=10)
    fig.colorbar(
        im, ax=axes, shrink=0.7, label="Top-peak log2 accessibility ratio (high / low)"
    )
    fig.suptitle(
        "Day-0 sensitivity | TSS >=3 | regional q across 25 neighborhoods", fontsize=13
    )
    for ext in ["png", "pdf"]:
        fig.savefig(figdir / f"Parameter_sensitivity.{ext}", dpi=180)
    plt.close(fig)
    print(json.dumps(audit), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    main(p.parse_args())
