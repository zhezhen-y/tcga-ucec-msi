"""Step 5: batch effects and confounding. Does the contributing institution leave a mark on expression,
and is it tangled with MSI status?

TCGA tumors come from ~30 tissue source sites (TSS), each with its own collection and handling. Three checks:
  1. PCA of expression: how much of each leading component does site explain, versus MSI status?
  2. Is MSI prevalence different across sites (chi-square)? If so, a classifier could learn site instead of biology.
  3. A mixed model, expression ~ MSI status + (1 | site), for MLH1 and the top DE genes: does the site random
     intercept change the MSI effect or its standard error compared with ordinary least squares?

Outputs (results/): pca_variance_by_site.tsv, mixed_model_site.tsv
Figures:  figures/pca_site_msi.png
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")  # MixedLM warns when the site variance sits at the zero boundary, which is the finding
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"

cohort = pd.read_csv(RES / "cohort.tsv", sep="\t", index_col=0)
counts = pd.read_csv(RES / "counts.tsv.gz", sep="\t", index_col=0)[cohort.index]
logcpm = np.log2(counts / counts.sum() * 1e6 + 1).T  # samples x genes
top = logcpm.var().sort_values(ascending=False).index[:2000]
X = StandardScaler().fit_transform(logcpm[top])

# ---------------- 1. PCA: site vs MSI ----------------
pca = PCA(n_components=10, random_state=0)
pcs = pd.DataFrame(pca.fit_transform(X), index=cohort.index, columns=[f"PC{i+1}" for i in range(10)])
d = cohort.join(pcs)
d["msi"] = d.msi_assay.map({"MSI-H": 1, "MSS": 0})
big_sites = d.tss.value_counts()[lambda s: s >= 5].index
rows = []
for i in range(10):
    pc = f"PC{i+1}"
    r2_site = smf.ols(f"{pc} ~ C(tss)", d[d.tss.isin(big_sites)]).fit().rsquared
    r2_msi = smf.ols(f"{pc} ~ msi", d.dropna(subset=["msi"])).fit().rsquared
    r2_sub = smf.ols(f"{pc} ~ C(subtype)", d).fit().rsquared
    rows.append({"PC": pc, "var_explained": pca.explained_variance_ratio_[i], "R2_site": r2_site, "R2_msi_assay": r2_msi, "R2_subtype": r2_sub})
tab = pd.DataFrame(rows)
tab.to_csv(RES / "pca_variance_by_site.tsv", sep="\t", index=False)
print("Share of each PC explained by site (R2, sites with >= 5 tumors), by MSI assay status, and by molecular subtype:")
print(tab.round(3).to_string(index=False))

# ---------------- 2. MSI prevalence by site ----------------
lab = d.dropna(subset=["msi"])
ct = pd.crosstab(lab.tss, lab.msi)
ct = ct[ct.sum(axis=1) >= 10]
chi2, p, _, _ = stats.chi2_contingency(ct)
print(f"\nMSI-H prevalence across {len(ct)} sites with >= 10 labeled tumors: chi-square p = {p:.3f}; "
      f"range {ct[1].div(ct.sum(axis=1)).min():.2f} to {ct[1].div(ct.sum(axis=1)).max():.2f}")

# ---------------- 3. mixed model for MLH1 and the top DE genes ----------------
de = pd.read_csv(RES / "de_msi_vs_cnlow.tsv", sep="\t", index_col=0)
genes = ["MLH1"] + [g for g in de.index[:6] if g != "MLH1"][:5]
rows = []
for g in genes:
    lab = lab.assign(expr=logcpm.loc[lab.index, g].values)
    ols = smf.ols("expr ~ msi", lab).fit()
    mm = smf.mixedlm("expr ~ msi", lab, groups=lab["tss"]).fit(reml=True)
    var_site, var_res = float(mm.cov_re.iloc[0, 0]), float(mm.scale)
    rows.append({"gene": g, "ols_coef": ols.params["msi"], "ols_se": ols.bse["msi"],
                 "mixed_coef": mm.params["msi"], "mixed_se": mm.bse["msi"],
                 "site_variance": var_site, "residual_variance": var_res, "ICC_site": var_site / (var_site + var_res)})
mmtab = pd.DataFrame(rows)
mmtab.to_csv(RES / "mixed_model_site.tsv", sep="\t", index=False)
print("\nMSI-H effect on log2 CPM, OLS vs mixed model with a site random intercept (ICC = share of residual variance between sites):")
print(mmtab.round(3).to_string(index=False))

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
ax = axes[0]
for lbl, c in [("MSS", "#7f8c8d"), ("MSI-H", "#c0392b")]:
    g = d[d.msi_assay == lbl]
    ax.scatter(g.PC1, g.PC2, s=10, c=c, alpha=0.7, label=f"{lbl} (n={len(g)})", rasterized=True)
ax.scatter(d[d.msi_assay.isna()].PC1, d[d.msi_assay.isna()].PC2, s=6, c="lightgray", alpha=0.5, label="no assay call", rasterized=True)
ax.set_title(f"By MSI assay status (R2 on PC1 = {tab.R2_msi_assay[0]:.2f})", fontsize=9)
ax.legend(fontsize=7)
ax = axes[1]
top_sites = d.tss.value_counts().index[:6]
cmap = plt.get_cmap("tab10")
for i, s in enumerate(top_sites):
    g = d[d.tss == s]
    ax.scatter(g.PC1, g.PC2, s=10, color=cmap(i), alpha=0.75, label=f"site {s} (n={len(g)})", rasterized=True)
ax.scatter(d[~d.tss.isin(top_sites)].PC1, d[~d.tss.isin(top_sites)].PC2, s=6, c="lightgray", alpha=0.5, label="other sites", rasterized=True)
ax.set_title(f"By tissue source site (R2 on PC1 = {tab.R2_site[0]:.2f})", fontsize=9)
ax.legend(fontsize=7)
for ax in axes:
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.0%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.0%})")
fig.suptitle("TCGA-UCEC expression, top 2,000 variable genes", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "pca_site_msi.png", dpi=200)
