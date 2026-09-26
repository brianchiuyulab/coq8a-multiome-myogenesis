"""Enumerate externally annotated GSE208248 scopes without changing stored tests."""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def main(root):
    out = root / "results/cav3_scope_grid"
    out.mkdir(parents=True, exist_ok=True)
    near = pd.read_csv(root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
    near = near[near.n_libraries.eq(4)]
    eff = pd.read_csv(
        root / "results/tables/candidate_peak_effects_pooled.tsv.gz", sep="\t"
    )
    eff = eff[eff.n_libraries.eq(4)].copy()
    eff["FC"] = eff.high_open / eff.low_open.replace(0, np.nan)
    ann = pd.read_csv(
        root / "results/external_enhancers/external_peak_annotations.tsv", sep="\t"
    )
    selective = pd.read_csv(
        root / "results/enhancer_regions/peak_selectivity.tsv", sep="\t"
    )
    funcs = pd.read_csv(root / "results/functional/functional_genes.tsv", sep="\t")
    scopes = {"Myogenesis221": set(near.gene)}
    scopes.update({f: set(d.gene) for f, d in funcs.groupby("function")})
    scopes["Differentiation_or_fusion"] = scopes["Differentiation"] | scopes["Fusion"]
    regulatory = {
        "All": set(near.peak),
        "HSMM": set(ann.loc[ann.strong_enhancer, "peak"]),
        "MYOD_bound_HSMM": set(ann.loc[ann.MYOD_bound_enhancer, "peak"]),
        "Muscle_selective": set(selective.loc[selective.n_other_strong.eq(0), "peak"]),
    }
    rows = []
    members = []
    full = []
    for scope, genes in scopes.items():
        for distance in [500, 2000, 100000]:
            pairs = near[near.gene.isin(genes) & near.nearest_tss_bp.le(distance)]
            for annotation, peaks in regulatory.items():
                selected = pairs[pairs.peak.isin(peaks)]
                chosen = set(selected.peak)
                if not chosen:
                    continue
                family = f"{scope}__TSS{distance}__{annotation}"
                members.extend(
                    dict(family=family, peak=p, gene=g)
                    for p, g in selected[["peak", "gene"]].itertuples(
                        index=False, name=None
                    )
                )
                for (gate, contrast), d in eff[eff.peak.isin(chosen)].groupby(
                    ["gate", "contrast"]
                ):
                    d = d.copy()
                    d["q_scope"] = multipletests(d.p_pair_binomial, method="fdr_bh")[1]
                    d["p_rank"] = d.p_pair_binomial.rank(method="min")
                    d["family"] = family
                    full.append(d)
                    cav = set(selected.loc[selected.gene.eq("CAV3"), "peak"])
                    c = d[d.peak.isin(cav)]
                    focal = d[d.peak.eq("chr3:8733438-8733968")]
                    rows.append(
                        dict(
                            family=family,
                            gene_scope=scope,
                            tss_distance=distance,
                            annotation=annotation,
                            gate=gate,
                            contrast=contrast,
                            n_peaks=len(d),
                            n_genes=selected.gene.nunique(),
                            n_open_q01=int(((d.FC > 1) & (d.q_scope < 0.1)).sum()),
                            n_closed_q01=int(((d.FC < 1) & (d.q_scope < 0.1)).sum()),
                            n_CAV3_peaks=len(c),
                            CAV3_min_q=float(c.q_scope.min()) if len(c) else np.nan,
                            CAV3_focal_FC=(
                                float(focal.FC.iloc[0]) if len(focal) else np.nan
                            ),
                            CAV3_focal_p=(
                                float(focal.p_pair_binomial.iloc[0])
                                if len(focal)
                                else np.nan
                            ),
                            CAV3_focal_q=(
                                float(focal.q_scope.iloc[0]) if len(focal) else np.nan
                            ),
                            CAV3_focal_rank=(
                                float(focal.p_rank.iloc[0]) if len(focal) else np.nan
                            ),
                        )
                    )
    grid = pd.DataFrame(rows)
    grid.to_csv(out / "scope_grid.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(members).to_csv(
        out / "scope_memberships.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.concat(full).to_csv(
        out / "all_peak_tests.tsv.gz",
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0},
    )
    manifest = dict(
        gene_sets=list(scopes),
        distances_bp=[500, 2000, 100000],
        external_annotations=list(regulatory),
        settings=grid[["gate", "contrast"]].drop_duplicates().to_dict("records"),
        n_scope_comparisons=len(grid),
        selection="Retrospective sensitivity exploration; all external subsets tested and retained; no one-gene CAV3-only family",
        q_definition="BH across all distinct member peaks within each setting; not across sensitivity settings",
        biological_sources=2,
    )
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        grid[grid.CAV3_focal_q.notna()]
        .sort_values(["CAV3_focal_q", "n_peaks"])
        .head(12)
        .to_string(index=False)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    main(p.parse_args().root)
