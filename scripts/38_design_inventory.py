"""Export input, membership and sample-size inventories without refitting effects."""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/design"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / "results/tables/gene_search_space.tsv", sep="\t")
    candidates = pd.read_csv(
        ROOT / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    common = candidates[candidates.n_libraries.eq(4)]
    temporal = pd.read_csv(
        ROOT / "results/temporal/bidirectional_peak_effects.tsv", sep="\t"
    )
    phases = temporal[["peak", "phase"]].drop_duplicates()
    rows = []
    for name, selected in [
        ("union", genes),
        ("intersection", genes[genes.Hallmark_Myogenesis & genes.Reactome_Myogenesis]),
    ]:
        members = common[common.gene.isin(selected.gene)]
        dynamic = phases[phases.peak.isin(members.peak)]
        selected.to_csv(OUT / f"{name}_genes.tsv", sep="\t", index=False)
        members.to_csv(OUT / f"{name}_common_peak_gene.tsv", sep="\t", index=False)
        rows.append(
            dict(
                scope=name,
                n_genes=len(selected),
                n_common_peaks=members.peak.nunique(),
                n_early=dynamic.phase.eq("Early").sum(),
                n_middle=dynamic.phase.eq("Middle").sum(),
                n_late=dynamic.phase.eq("Late").sum(),
                n_closing=dynamic.phase.eq("Closing").sum(),
            )
        )
    scope = pd.DataFrame(rows)
    assert scope.n_genes.tolist() == [221, 6]
    assert scope.n_common_peaks.tolist() == [5097, 146]
    scope.to_csv(OUT / "union_intersection_inventory.tsv", sep="\t", index=False)
    nuclei = pd.read_csv(ROOT / "results/tables/nuclei.tsv.gz", sep="\t")
    qc = pd.read_csv(ROOT / "results/tables/qc_by_library.tsv", sep="\t")
    qc = qc.merge(
        nuclei.groupby("gsm").size().rename("n_input_filtered_barcodes"), on="gsm"
    )
    pairs = pd.read_csv(ROOT / "results/tables/matched_pairs.tsv.gz", sep="\t")
    main_pairs = pairs[pairs.gate.eq("TSS_ge_3") & pairs.contrast.eq("3plus_vs_1")]
    qc = qc.merge(
        main_pairs.groupby("gsm").size().rename("n_main_matched_pairs"), on="gsm"
    )
    qc.to_csv(OUT / "target_sample_inventory.tsv", sep="\t", index=False)
    ext = []
    for experiment in [1, 2]:
        cells = pd.read_csv(
            ROOT / f"results/temporal/HSMM{experiment}_cells.tsv.gz", sep="\t"
        )
        for hour, frame in cells.groupby("hour"):
            ext.append(
                dict(
                    experiment=experiment,
                    gsm=f"GSM{2970929+experiment}",
                    hour=hour,
                    n_input_cells=len(frame),
                    n_primary_retained=int(frame.pass_2kb.sum()),
                    n_sensitivity_retained=int(frame.pass_1kb.sum()),
                )
            )
    pd.DataFrame(ext).to_csv(
        OUT / "external_sample_inventory.tsv", sep="\t", index=False
    )
    print(scope.to_string(index=False))
    print(
        qc[
            ["gsm", "n_input_filtered_barcodes", "n_tss3", "n_main_matched_pairs"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
