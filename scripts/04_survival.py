"""Step 4: time-to-event analysis by molecular subtype (Kaplan-Meier, log-rank, Cox).

Endpoint: progression-free interval (PFI) from TCGA-CDR. Liu et al. (2018, Cell) recommend PFI for
UCEC: follow-up is short and few patients die, so overall survival (OS) is underpowered. OS is
run too, as a secondary check. Times are in days, converted to years.

Cox model: PFI ~ subtype (reference CN-low) + age (per 10 years) + FIGO stage (I vs II vs III/IV).
The proportional-hazards assumption is checked with scaled Schoenfeld residuals.

Outputs (results/): logrank_pfi.tsv, cox_pfi.tsv, cox_os.tsv
Figures:  figures/km_pfi_subtype.png, figures/km_os_subtype.png, figures/cox_pfi_forest.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ORDER = ["POLE", "MSI", "CN_LOW", "CN_HIGH"]
COLORS = {"POLE": "#2a9d8f", "MSI": "#e9c46a", "CN_LOW": "#2c6fb0", "CN_HIGH": "#c0392b"}

cohort = pd.read_csv(RES / "cohort.tsv", sep="\t", index_col=0)


def km_plot(endpoint, fname):
    d = cohort.dropna(subset=[endpoint, f"{endpoint}.time"]).copy()
    d["years"] = d[f"{endpoint}.time"] / 365.25
    fig, ax = plt.subplots(figsize=(6, 5))
    fitters = []
    for s in ORDER:
        g = d[d.subtype == s]
        kmf = KaplanMeierFitter(label=f"{s} (n={len(g)})")
        kmf.fit(g.years, g[endpoint])
        kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS[s], show_censors=True,
                                   censor_styles={"ms": 3})
        fitters.append(kmf)
    lr = multivariate_logrank_test(d.years, d.subtype, d[endpoint])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_xlabel("years from diagnosis")
    ax.set_ylabel(f"{endpoint} probability")
    ax.set_title(f"TCGA-UCEC {endpoint} by molecular subtype (log-rank p = {lr.p_value:.1e})", fontsize=10)
    add_at_risk_counts(*fitters, ax=ax, rows_to_show=["At risk"], fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=200)
    print(f"{endpoint}: n={len(d)}, events={int(d[endpoint].sum())}, 4-group log-rank p={lr.p_value:.2e}")
    for kmf in fitters:
        print(f"  {kmf._label:<16} 5-yr {endpoint} = {kmf.predict(5):.2f}")
    return d


d = km_plot("PFI", "km_pfi_subtype.png")
pw = pairwise_logrank_test(d.years, d.subtype, d.PFI).summary
pw["p_bonferroni"] = (pw.p * len(pw)).clip(upper=1)
pw.to_csv(RES / "logrank_pfi.tsv", sep="\t")
print(pw[["test_statistic", "p", "p_bonferroni"]].round(4).to_string())
km_plot("OS", "km_os_subtype.png")


def cox(endpoint):
    d = cohort.dropna(subset=[endpoint, f"{endpoint}.time", "age", "stage"]).copy()
    d["years"] = d[f"{endpoint}.time"] / 365.25
    d = d[d.years > 0]
    X = pd.DataFrame({"years": d.years, "event": d[endpoint].astype(int),
                      "age_per10y": d.age / 10})
    for s in ["POLE", "MSI", "CN_HIGH"]:
        X[f"subtype_{s}"] = (d.subtype == s).astype(int)
    X["stage_II"] = (d.stage == "II").astype(int)
    X["stage_III_IV"] = d.stage.isin(["III", "IV"]).astype(int)
    cph = CoxPHFitter()
    cph.fit(X, duration_col="years", event_col="event")
    print(f"\nCox {endpoint} (n={len(X)}, events={X.event.sum()}), reference: CN-low, stage I")
    out = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    print(out.round(3).to_string())
    print(f"concordance = {cph.concordance_index_:.3f}")
    return cph, X, out


cph, X, out = cox("PFI")
out.to_csv(RES / "cox_pfi.tsv", sep="\t")
print("\nProportional-hazards check (scaled Schoenfeld residuals):")
ph = __import__("lifelines.statistics", fromlist=["proportional_hazard_test"]).proportional_hazard_test(
    cph, X, time_transform="rank")
print(ph.summary[["test_statistic", "p"]].round(3).to_string())
ph.summary.to_csv(RES / "cox_pfi_ph_test.tsv", sep="\t")
cox("OS")[2].to_csv(RES / "cox_os.tsv", sep="\t")

fig, ax = plt.subplots(figsize=(6, 3.2))
cph.plot(hazard_ratios=True, ax=ax)
ax.set_title("Cox model, PFI (reference: CN-low, FIGO stage I)", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "cox_pfi_forest.png", dpi=200)
