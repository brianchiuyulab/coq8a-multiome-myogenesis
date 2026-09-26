"""Export reproducible candidate summaries, diagnostics and scientific figures."""

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/replication_sensitivity"
OUT = ROOT / "results/replication_summary"
FIG = ROOT / "figures/replication"
CAV = "chr3-8733222-8734203"
OLD = "chr3:8733438-8733968"
FOCAL = [CAV, "chr1-211133720-211134678", "chr16-1311478-1312392"]


def read(path):
    return pd.read_csv(path, sep="\t")


def save(df, name):
    df.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 220,
        }
    )


def export(fig, name):
    fig.savefig(FIG / (name + ".png"), bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    configs = read(R / "configurations.tsv")
    a = read(R / "ATAC_donor_models.tsv.gz")
    rna = read(R / "RNA_donor_models.tsv.gz")
    cav = a[a.peak.eq(CAV)].merge(configs.drop(columns="n_donors"), on="config")
    save(cav, "CAV3_all_count_models.tsv")
    focal = a[a.peak.isin(FOCAL)].merge(configs.drop(columns="n_donors"), on="config")
    save(focal, "focal_all_count_models.tsv")
    save(
        rna[rna.gene.isin(["CAV3", "CACNA1H", "KCNH1"])].merge(
            configs.drop(columns="n_donors"), on="config"
        ),
        "focal_RNA_models.tsv",
    )
    status = read(R / "model_status.tsv")
    an = a.groupby("config").size()
    rn = rna.groupby("config").size()
    status["n_ATAC_written"] = status.config.map(an).fillna(0).astype(int)
    status["n_RNA_written"] = status.config.map(rn).fillna(0).astype(int)
    status["ATAC_status"] = np.where(status.n_ATAC_written.gt(0), "ok", status.status)
    status["RNA_status"] = np.where(status.n_RNA_written.gt(0), "ok", status.status)
    save(status, "model_status_reconciled.tsv")
    save(
        cav[
            (cav.qc == "author")
            & (cav.rule == "3plus_vs_1")
            & (cav.caliper == "0.1")
            & (cav.time == "Pre")
        ],
        "CAV3_anchor_results.tsv",
    )
    best = cav[cav.FC.gt(1) & cav.n_donors.ge(3)].sort_values("PValue")
    save(best, "CAV3_favorable_ranked.tsv")
    links = []
    donors = []
    for folder, label in [
        ("replication_links", "minimum_3_events"),
        ("replication_links_relaxed", "minimum_1_event"),
    ]:
        d = read(ROOT / "results" / folder / "COQ_ATAC_continuous.tsv.gz")
        d = d[d.peak.isin(FOCAL)].copy()
        d["eligibility"] = label
        links.append(d)
        x = read(ROOT / "results" / folder / "COQ_ATAC_by_donor.tsv.gz")
        x = x[x.peak.isin(FOCAL)].copy()
        x["eligibility"] = label
        donors.append(x)
    save(pd.concat(links), "focal_continuous_sensitivity.tsv")
    ds = pd.concat(donors)
    save(ds, "focal_continuous_by_donor.tsv")
    pl = read(ROOT / "results/replication_links_relaxed/peak_RNA_links.tsv.gz")
    save(pl[pl.peak.isin(FOCAL)], "focal_peak_RNA_links.tsv")
    # Precision-aware Fisher-z meta-analysis, including low-information donors.
    meta = []
    for keys, d in ds[ds.eligibility.eq("minimum_1_event")].groupby(
        ["celltype", "time", "model", "peak"]
    ):
        d = d[np.isfinite(d.r) & (d.n_nuclei > d.n_covariates + 2)]
        if len(d) < 2:
            continue
        z = np.arctanh(np.clip(d.r.to_numpy(), -0.999999, 0.999999))
        v = 1 / (d.n_nuclei.to_numpy() - d.n_covariates.to_numpy() - 2)
        w = 1 / v
        mu = (w * z).sum() / w.sum()
        q = (w * (z - mu) ** 2).sum()
        tau = max(0, (q - (len(d) - 1)) / (w.sum() - (w * w).sum() / w.sum()))
        w = 1 / (v + tau)
        mu = (w * z).sum() / w.sum()
        se = np.sqrt(max(1, (w * (z - mu) ** 2).sum() / (len(d) - 1)) / w.sum())
        p = 2 * stats.t.sf(abs(mu / se), len(d) - 1)
        ci = mu + np.array([-1, 1]) * stats.t.ppf(0.975, len(d) - 1) * se
        meta.append(
            dict(zip(["celltype", "time", "model", "peak"], keys))
            | dict(
                n_donors=len(d),
                r=np.tanh(mu),
                lower=np.tanh(ci[0]),
                upper=np.tanh(ci[1]),
                p=p,
                tau2=tau,
            )
        )
    save(pd.DataFrame(meta), "focal_precision_meta.tsv")
    # Validate emitted family adjustments independently of the fitting script.
    sc = read(R / "tested_peak_scopes.tsv").set_index("peak")
    for cfg in ["C296", "C297", "C198"]:
        d = a[a.config.eq(cfg)]
        for family in ["all_candidates", "HSMM", "muscle_selective", "fusion_promoter"]:
            z = d[sc.loc[d.peak, family].to_numpy()]
            np.testing.assert_allclose(
                z["q_" + family],
                multipletests(z.PValue, method="fdr_bh")[1],
                rtol=1e-10,
                atol=1e-12,
            )
    old = read(ROOT / "results/cav3_scope_grid/all_peak_tests.tsv.gz")
    old = old[
        (old.family == "Fusion__TSS500__Muscle_selective")
        & (old.gate == "TSS_ge_2")
        & (old.contrast == "3plus_vs_1")
    ]
    assert len(old) == 2
    np.testing.assert_allclose(
        old.q_scope, multipletests(old.p_pair_binomial, method="fdr_bh")[1]
    )
    np.testing.assert_allclose(old.FC, old.high_open / old.low_open)
    save(old, "old_selected_two_promoters.tsv")
    plot_old(old)
    plot_grid(cav)
    plot_new(cav)
    plot_continuous(ds)
    plot_distal(focal)
    metrics = {
        "configurations": len(configs),
        "ATAC_estimable": int(status.n_ATAC_written.gt(0).sum()),
        "RNA_estimable": int(status.n_RNA_written.gt(0).sum()),
        "CAV3_count_estimable": len(cav),
        "CAV3_positive_nominal": int((cav.FC.gt(1) & cav.PValue.lt(0.05)).sum()),
        "CAV3_positive_nominal_matched": int(
            (cav.FC.gt(1) & cav.PValue.lt(0.05) & cav.caliper.ne("none")).sum()
        ),
        "CAV3_q01_fusion_positive": int(
            (cav.FC.gt(1) & cav.q_fusion_promoter.lt(0.1)).sum()
        ),
        "checks": "Selected BH families and old exact FC/q reconciled; full pseudobulk sums verified in scripts 66/67.",
    }
    (OUT / "summary.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


def plot_old(old):
    members = read(ROOT / "results/cav3_scope_grid/scope_memberships.tsv.gz")
    d = members[members.family.eq("Fusion__TSS500__All")][
        ["peak", "gene"]
    ].drop_duplicates()
    d = (
        d.groupby("peak")
        .gene.agg(", ".join)
        .reset_index()
        .sort_values(["gene", "peak"])
    )
    ann = read(
        ROOT / "results/external_enhancers/external_peak_annotations.tsv"
    ).set_index("peak")
    sel = read(ROOT / "results/enhancer_regions/peak_selectivity.tsv").set_index("peak")
    mat = np.column_stack(
        [
            ann.reindex(d.peak).strong_enhancer.fillna(False).astype(int),
            ann.reindex(d.peak).MYOD_bound_enhancer.fillna(False).astype(int),
            sel.reindex(d.peak).n_other_strong.eq(0).astype(int),
        ]
    )
    save(
        d.assign(HSMM=mat[:, 0], MYOD_bound=mat[:, 1], muscle_selective=mat[:, 2]),
        "fusion_promoter_annotation_matrix.tsv",
    )
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.4, 1], height_ratios=[1, 1], wspace=0.65, hspace=0.5
    )
    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(
        mat, cmap=ListedColormap(["#edf1f4", "#236f93"]), vmin=0, vmax=1, aspect="auto"
    )
    ax.set_yticks(
        range(len(d)),
        [
            f'{g}  |  {p.replace("chr","")}'
            for p, g in d.itertuples(index=False, name=None)
        ],
        fontsize=8,
    )
    ax.set_xticks(
        range(3),
        ["HSMM\nregulatory state", "MYOD1\noccupancy", "Muscle\nselectivity"],
        fontsize=8,
    )
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    ax.set_title(
        "A   Fusion-associated promoters", loc="left", pad=43, fontweight="bold"
    )
    ax.set_xlabel("Blue: external annotation present", labelpad=12)
    for i, g in enumerate(d.gene):
        if g in ["CAV3", "KCNH1"]:
            ax.get_yticklabels()[i].set_color("#a4482d")
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    steps = [
        "221 myogenesis genes",
        "12 fusion-annotated genes",
        "21 shared promoter peaks\n(midpoint within 500 bp of TSS)",
        "2 muscle-selective promoters\nCAV3 and KCNH1",
    ]
    for i, t in enumerate(steps):
        y = 0.9 - i * 0.26
        ax.text(
            0.5,
            y,
            t,
            ha="center",
            va="center",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=.5", fc="#edf4f7", ec="#90afbd"),
            transform=ax.transAxes,
        )
        if i < 3:
            ax.annotate(
                "",
                xy=(0.5, y - 0.17),
                xytext=(0.5, y - 0.065),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#526d7a"),
            )
    ax.set_title("B   External functional scope", loc="left", fontweight="bold", pad=16)
    ax = fig.add_subplot(gs[1, 1])
    low = old.low_open / old.n_pairs * 100
    high = old.high_open / old.n_pairs * 100
    for j, row in enumerate(old.itertuples()):
        ax.plot([j - 0.16, j + 0.16], [low.iloc[j], high.iloc[j]], color="#c3cbd0")
        ax.scatter(
            j - 0.16,
            low.iloc[j],
            color="#7b8e9c",
            s=48,
            label="COQ8A low" if j == 0 else None,
        )
        ax.scatter(
            j + 0.16,
            high.iloc[j],
            color="#b95538",
            s=48,
            label="COQ8A high" if j == 0 else None,
        )
        ax.text(
            j,
            max(low.iloc[j], high.iloc[j]) + 1.6,
            f"FC {row.FC:.2f}\np = {row.p_pair_binomial:.3f}\nq = {row.q_scope:.3f}",
            ha="center",
            fontsize=9,
        )
    ax.set_xticks(range(len(old)), ["CAV3" if p == OLD else "KCNH1" for p in old.peak])
    ax.set_ylabel("Nuclei with detected peak (%)")
    ax.set_ylim(0, max(high.max(), low.max()) + 12)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.set_title(
        "C   GSE208248 accessibility\nTSS ≥2; COQ8A ≥3 vs 1 UMI",
        loc="left",
        fontweight="bold",
        pad=14,
    )
    export(fig, "Figure_1_CAV3_functional_scope")


def plot_grid(cav):
    rules = [
        "2plus_vs_1",
        "3plus_vs_1",
        "4plus_vs_1",
        "3plus_vs_1to2",
        "positive_q25",
        "positive_q33",
        "detected_vs_zero",
    ]
    labels = ["≥2 / 1", "≥3 / 1", "≥4 / 1", "≥3 / 1–2", "Q25", "Q33", ">0 / 0"]
    groups = [
        (ct, t)
        for ct in ["Satellite Cells", "Fast", "Slow", "Intermediate"]
        for t in ["Pre", "Post"]
    ]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True, layout="constrained")
    for ax, cal, title in zip(
        axs,
        ["none", "0.1"],
        [
            "A   All eligible nuclei; RNA-depth covariate",
            "B   Depth-matched nuclei; caliper 0.10",
        ],
    ):
        mat = np.full((8, 7), np.nan)
        pv = mat.copy()
        ns = mat.copy()
        for i, (ct, t) in enumerate(groups):
            for j, rule in enumerate(rules):
                d = cav[
                    (cav.celltype == ct)
                    & (cav.time == t)
                    & (cav.rule == rule)
                    & (cav.qc == "author")
                    & (cav.caliper == cal)
                ]
                if len(d):
                    mat[i, j] = d.logFC.iloc[0]
                    pv[i, j] = d.PValue.iloc[0]
                    ns[i, j] = d.n_donors.iloc[0]
        cm = plt.colormaps["RdBu_r"].copy()
        cm.set_bad("#e5e5e5")
        im = ax.imshow(mat, cmap=cm, vmin=-2, vmax=2, aspect="auto")
        for i in range(8):
            for j in range(7):
                if np.isfinite(mat[i, j]):
                    ax.text(
                        j,
                        i,
                        f"{2**mat[i,j]:.2f}" + ("*" if pv[i, j] < 0.05 else ""),
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if abs(mat[i, j]) > 1.4 else "black",
                    )
        ax.set_xticks(range(7), labels, rotation=35, ha="right")
        ax.set_yticks(
            range(8),
            [f'{"MuSC" if ct=="Satellite Cells" else ct} | {t}' for ct, t in groups],
        )
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
        ax.set_xlabel("COQ8A high / low definition")
    fig.colorbar(
        im,
        ax=axs,
        label="CAV3 normalized-count log₂ fold change",
        shrink=0.8,
        extend="both",
    )
    fig.suptitle(
        "CAV3 accessibility across cell states and expression definitions\nCell values: fold change; * nominal donor-model p < 0.05; grey: not estimable",
        fontsize=11,
    )
    export(fig, "Figure_2_CAV3_sensitivity")


def plot_new(cav):
    d = read(R / "CAV3_detection_by_donor.tsv")
    d = d[d.config.eq("C297")]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6), layout="constrained")
    ax = axs[0]
    for j, row in enumerate(d.itertuples()):
        ax.plot(
            [j - 0.15, j + 0.15],
            [row.low_fraction * 100, row.high_fraction * 100],
            color="#b5bfc6",
        )
        ax.scatter(
            j - 0.15,
            row.low_fraction * 100,
            color="#7b8e9c",
            s=40,
            label="COQ8A low" if j == 0 else None,
        )
        ax.scatter(
            j + 0.15,
            row.high_fraction * 100,
            color="#b95538",
            s=40,
            label="COQ8A high" if j == 0 else None,
        )
    ax.set_xticks(
        range(len(d)), [f"{r.donor}\n{r.n_high} pairs" for r in d.itertuples()]
    )
    ax.set_ylim(0, 100)
    ax.set_ylabel("Nuclei with detected CAV3 peak (%)")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.set_title(
        "A   Six paired donors\nPooled detection FC 1.32; sign-flip p = 0.031",
        loc="left",
        fontweight="bold",
        fontsize=10,
    )
    ax = axs[1]
    for i, cfg in enumerate(["C296", "C297", "C298", "C299"]):
        r = cav[cav.config.eq(cfg)].iloc[0]
        ax.bar(
            i, r.FC, color=["#83929b", "#b95538", "#c88b75", "#dab7a9"][i], width=0.6
        )
        ax.text(
            i,
            r.FC + 0.12,
            f"{r.FC:.2f}×\np = {r.PValue:.4f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.axhline(1, color="#4a565e", lw=0.8, ls="--")
    ax.set_ylim(0, 4.3)
    ax.set_xticks(
        range(4),
        [
            "All nuclei\n+ RNA covariate",
            "Matched\n0.10",
            "Matched\n0.20",
            "Matched\n0.30",
        ],
    )
    ax.tick_params(axis="x", labelsize=8)
    ax.set_ylabel("Normalized peak-count fold change")
    ax.set_title(
        "B   Dependence on depth handling", loc="left", fontweight="bold", fontsize=10
    )
    fig.suptitle(
        "GSE240061 | Slow myonuclei, post sampling\nCOQ8A-positive nuclei: upper vs lower expression quartile",
        fontsize=11,
    )
    export(fig, "Figure_3_CAV3_donor_followup")


def plot_continuous(ds):
    d = ds[
        (ds.eligibility == "minimum_1_event")
        & (ds.celltype == "Satellite Cells")
        & (ds.time == "Pre")
        & (ds.model == "technical")
        & (ds.peak == CAV)
    ].sort_values("donor")
    fig, ax = plt.subplots(figsize=(6.5, 4), layout="constrained")
    ax.scatter(
        d.r,
        range(len(d)),
        s=55,
        color=["#b95538" if r >= 0 else "#346f99" for r in d.r],
    )
    ax.axvline(0, color="#777", ls="--", lw=0.8)
    ax.set_yticks(
        range(len(d)),
        [
            f"{r.donor}: {r.n_nuclei} nuclei, {int(r.n_COQ_detected)} COQ8A+"
            for r in d.itertuples()
        ],
    )
    ax.set_xlabel("Within-donor partial correlation")
    ax.set_title(
        "Baseline MuSC: COQ8A and CAV3 accessibility",
        loc="left",
        fontweight="bold",
        fontsize=11,
    )
    ax.invert_yaxis()
    export(fig, "Supplement_MuSC_donor_sensitivity")


def plot_distal(focal):
    peak = "chr16-1311478-1312392"
    d = read(ROOT / "results/replication_stress/candidate_counts_by_donor.tsv")
    d = d[(d.config == "C198") & (d.peak == peak)]
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
    ax = axs[0]
    for i, (donor, g) in enumerate(d.groupby("donor")):
        g = g.set_index("group")
        ys = [g.loc["low", "CPM"], g.loc["high", "CPM"]]
        ax.plot([i - 0.15, i + 0.15], ys, color="#bdc7cd")
        ax.scatter(
            i - 0.15, ys[0], color="#7b8e9c", label="COQ8A low" if i == 0 else None
        )
        ax.scatter(
            i + 0.15, ys[1], color="#b95538", label="COQ8A high" if i == 0 else None
        )
    ax.set_xticks(range(6), sorted(d.donor.unique()))
    ax.set_ylabel("TMM-normalized peak counts per million")
    ax.set_ylim(0, 25)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(
        "A   Six donors, 490 matched pairs\nFC 3.01; p = 8.19 × 10⁻⁶; q = 0.022",
        loc="left",
        fontsize=10,
        fontweight="bold",
    )
    ax = axs[1]
    for j, qc in enumerate(["author", "trim95"]):
        dd = focal[
            (focal.peak == peak)
            & (focal.celltype == "Fast")
            & (focal.time == "Post")
            & (focal.rule == "2plus_vs_1")
            & (focal.qc == qc)
            & (focal.caliper.isin(["0.1", "0.2", "0.3"]))
        ].sort_values("caliper")
        xs = np.arange(3) + (0.06 if j else -0.06)
        ax.plot(
            xs,
            dd.FC,
            "o-",
            label="Author QC" if j == 0 else "Depth cap95",
            color=["#607f95", "#b95538"][j],
        )
        for x, r in zip(xs, dd.itertuples()):
            ax.text(
                x,
                r.FC + (0.15 if j else -0.25),
                f"q={r.q_all_candidates:.3f}",
                ha="center",
                fontsize=8,
            )
    ax.axhline(1, color="#777", ls="--", lw=0.8)
    ax.set_ylim(0.8, 3.8)
    ax.set_xlim(-0.30, 2.30)
    ax.set_xticks(range(3), ["0.10", "0.20", "0.30"])
    ax.set_xlabel("Depth-matching caliper")
    ax.set_ylabel("Normalized peak-count fold change")
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 0.08), fontsize=8)
    ax.set_title(
        "B   QC and matching sensitivity", loc="left", fontsize=10, fontweight="bold"
    )
    fig.suptitle(
        "Fast myonuclei, post sampling | COQ8A ≥2 vs 1 UMI\nhg38 chr16:1311478–1312392",
        fontsize=11,
    )
    export(fig, "Supplement_distal_candidate")


if __name__ == "__main__":
    main()
