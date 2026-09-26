"""Summarize matched detection effects and donor-level sensitivity statistics."""

from pathlib import Path
import argparse
import itertools
import gzip
import numpy as np
import pandas as pd
from scipy import io, stats
from statsmodels.stats.multitest import multipletests


def bh(p):
    return multipletests(np.where(np.isfinite(p), p, 1), method="fdr_bh")[1]


def donor_test(delta):
    n = delta.shape[1]
    mean = delta.mean(1)
    if n < 2:
        return np.full(len(mean), np.nan), np.full(len(mean), np.nan)
    sd = delta.std(1, ddof=1) / np.sqrt(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 2 * stats.t.sf(abs(mean / sd), n - 1)
    p = np.where((sd == 0) & (mean == 0), 1, p)
    signs = np.asarray(list(itertools.product([-1, 1], repeat=n)), dtype=float)
    perm = delta @ signs.T / n
    exact = (np.abs(perm) >= np.abs(mean[:, None]) - 1e-14).mean(1)
    return p, exact


def main(root, work, out):
    configs = pd.read_csv(out / "configurations.tsv", sep="\t")
    columns = pd.read_csv(out / "pseudobulk_columns.tsv", sep="\t")
    scope = pd.read_csv(out / "tested_peak_scopes.tsv", sep="\t")
    peaks = (work / "candidate_ATAC_peaks.txt").read_text().splitlines()
    idx = pd.Index(peaks).get_indexer(scope.peak)
    open_counts = np.load(work / "sensitivity_open_counts.npy")[idx]
    binary = io.mmread(work / "candidate_ATAC_counts.mtx").tocsc()[idx]
    binary.data[:] = 1
    barcodes = (work / "barcodes.txt").read_text().splitlines()
    full_barcodes = (work / "all_ATAC_barcodes.txt").read_text().splitlines()
    binary = binary[:, pd.Index(barcodes).get_indexer(full_barcodes)]
    pairs = pd.read_csv(work / "sensitivity_pairs.tsv.gz", sep="\t")
    family_cols = list(scope.columns[1:])
    modules = []
    donor_rows = []
    focal = []
    with gzip.open(out / "detection_sensitivity.tsv.gz", "wt", newline="") as dest:
        first = True
        for cfg in configs.itertuples():
            c = columns[columns.config.eq(cfg.config)]
            if c.empty:
                continue
            hi = c[c.group.eq("high")].sort_values("donor")
            lo = c[c.group.eq("low")].sort_values("donor")
            assert hi.donor.tolist() == lo.donor.tolist()
            nh = hi.n_nuclei.to_numpy()
            nl = lo.n_nuclei.to_numpy()
            hc = open_counts[:, hi.column]
            lc = open_counts[:, lo.column]
            h = hc / nh
            l = lc / nl
            delta = h - l
            pt, pp = donor_test(delta)
            with np.errstate(divide="ignore", invalid="ignore"):
                fc = h.mean(1) / l.mean(1)
                fcpooled = (hc.sum(1) / nh.sum()) / (lc.sum(1) / nl.sum())
            df = pd.DataFrame(
                dict(
                    config=cfg.config,
                    peak=scope.peak,
                    n_donors=len(hi),
                    n_high=nh.sum(),
                    n_low=nl.sum(),
                    FC_equal_donor=fc,
                    FC_pooled=fcpooled,
                    delta_pp=100 * delta.mean(1),
                    positive_donors=(delta > 1e-14).sum(1),
                    negative_donors=(delta < -1e-14).sum(1),
                    high_open_total=hc.sum(1),
                    low_open_total=lc.sum(1),
                    p_donor_t=pt,
                    p_donor_signflip=pp,
                )
            )
            pair = pairs[pairs.config.eq(cfg.config)]
            df["p_paired_nuclei"] = np.nan
            if len(pair):
                high = binary[:, pair.high_index]
                low = binary[:, pair.low_index]
                common = np.asarray(high.multiply(low).sum(1)).ravel()
                a = np.asarray(high.sum(1)).ravel() - common
                b = np.asarray(low.sum(1)).ravel() - common
                pv = np.minimum(1, 2 * stats.binom.cdf(np.minimum(a, b), a + b, 0.5))
                df["p_paired_nuclei"] = pv
            for family in family_cols:
                mask = scope[family].to_numpy(bool)
                for test in ["donor_t", "donor_signflip", "paired_nuclei"]:
                    df["q_" + test + "_" + family] = np.nan
                    if df.loc[mask, "p_" + test].notna().any():
                        df.loc[mask, "q_" + test + "_" + family] = bh(
                            df.loc[mask, "p_" + test]
                        )
                mh = h[mask].mean(0)[None, :]
                ml = l[mask].mean(0)[None, :]
                mt, mp = donor_test(mh - ml)
                modules.append(
                    dict(
                        config=cfg.config,
                        module=family,
                        n_peaks=int(mask.sum()),
                        n_donors=len(hi),
                        FC_equal_donor=(
                            float(mh.mean() / ml.mean()) if ml.mean() else np.nan
                        ),
                        delta_pp=float(100 * (mh - ml).mean()),
                        p_donor_t=mt[0],
                        p_donor_signflip=mp[0],
                        positive_donors=int((mh > ml).sum()),
                    )
                )
            cav = np.flatnonzero(scope.peak.eq("chr3-8733222-8734203"))
            if len(cav):
                j = cav[0]
                focal.append(df.iloc[[j]])
                donor_rows.extend(
                    dict(
                        config=cfg.config,
                        donor=donor,
                        n_high=int(nh[k]),
                        n_low=int(nl[k]),
                        high_open=int(hc[j, k]),
                        low_open=int(lc[j, k]),
                        high_fraction=h[j, k],
                        low_fraction=l[j, k],
                        delta_pp=100 * delta[j, k],
                    )
                    for k, donor in enumerate(hi.donor)
                )
            df.to_csv(dest, sep="\t", index=False, header=first, na_rep="NA")
            first = False
    pd.concat(focal).to_csv(
        out / "CAV3_detection_grid.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame(donor_rows).to_csv(
        out / "CAV3_detection_by_donor.tsv", sep="\t", index=False
    )
    module = pd.DataFrame(modules)
    for _, d in module.groupby("config"):
        for test in ["donor_t", "donor_signflip"]:
            module.loc[d.index, "q_" + test] = (
                bh(d["p_" + test]) if d["p_" + test].notna().any() else np.nan
            )
    module.to_csv(out / "module_sensitivity.tsv", sep="\t", index=False, na_rep="NA")
    print(
        "Detection and module effects complete:", len(configs), "settings", flush=True
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work, a.out)
