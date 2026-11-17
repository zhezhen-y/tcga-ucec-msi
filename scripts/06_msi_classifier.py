"""Step 6: predict the clinical MSI call from RNA-seq, and ask whether machine learning adds anything
over a single gene.

Label: the clinical MSI assay call (MSI-H vs MSS; MSI-L and indeterminate excluded), which is
independent of the expression data. Four models, from simplest to most flexible:
  MLH1 alone        logistic regression on one gene (the known silenced gene in sporadic MSI)
  10-gene panel     genes picked by F-test *inside each training fold*, then logistic regression
  elastic net       penalized logistic regression on all 3,000 most variable genes, C tuned by inner CV
  random forest     500 trees on the same genes

Evaluation: repeated stratified 5-fold CV (3 repeats), nested for the elastic-net penalty, reporting
AUROC, average precision and balanced accuracy. Then leave-one-site-out: each tissue source site with
>= 10 labeled tumors is held out in turn, so a model that had learned site signatures would fall apart.

Outputs (results/): classifier_cv.tsv, classifier_site_holdout.tsv, elastic_net_coefficients.tsv
Figures:  figures/roc_msi_classifier.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
SEED = 7

cohort = pd.read_csv(RES / "cohort.tsv", sep="\t", index_col=0)
lab = cohort[cohort.msi_assay.isin(["MSI-H", "MSS"])]
counts = pd.read_csv(RES / "counts.tsv.gz", sep="\t", index_col=0)[lab.index]
logcpm = np.log2(counts / counts.sum() * 1e6 + 1).T
genes = logcpm.var().sort_values(ascending=False).index[:3000]  # unsupervised filter, no label used
X = logcpm[genes].values
y = (lab.msi_assay == "MSI-H").astype(int).values
sites = lab.tss.values
print(f"labeled tumors: {len(y)} ({y.sum()} MSI-H, {len(y) - y.sum()} MSS); features: {X.shape[1]} genes")

mlh1 = list(genes).index("MLH1")


class ColumnPicker:
    """Keeps one column; lets the single-gene model run through the same pipeline code."""
    def __init__(self, col): self.col = col
    def fit(self, X, y=None): return self
    def transform(self, X): return X[:, [self.col]]
    def get_params(self, deep=True): return {"col": self.col}
    def set_params(self, **p): self.col = p.get("col", self.col); return self


inner = StratifiedKFold(3, shuffle=True, random_state=SEED)
models = {
    "MLH1 alone": Pipeline([("pick", ColumnPicker(mlh1)), ("sc", StandardScaler()), ("lr", LogisticRegression())]),
    "10-gene panel (in-fold selection)": Pipeline([("sel", SelectKBest(f_classif, k=10)), ("sc", StandardScaler()), ("lr", LogisticRegression())]),
    "elastic net, 3,000 genes": GridSearchCV(
        Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, max_iter=3000, tol=1e-3))]),
        {"lr__C": [0.01, 0.03, 0.1, 0.3]}, cv=inner, scoring="roc_auc", n_jobs=4),
    "random forest, 3,000 genes": RandomForestClassifier(n_estimators=500, min_samples_leaf=2, random_state=SEED, n_jobs=4),
}

# ---------------- repeated nested CV ----------------
outer = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=SEED)
rows, pooled = [], {m: (np.zeros(len(y)), np.zeros(len(y))) for m in models}  # summed prob, count
for fold, (tr, te) in enumerate(outer.split(X, y)):
    for name, model in models.items():
        m = clone(model).fit(X[tr], y[tr])
        p = m.predict_proba(X[te])[:, 1]
        rows.append({"model": name, "fold": fold, "auroc": roc_auc_score(y[te], p), "auprc": average_precision_score(y[te], p),
                     "balanced_acc": balanced_accuracy_score(y[te], p >= 0.5),
                     "C": m.best_params_["lr__C"] if isinstance(m, GridSearchCV) else np.nan})
        pooled[name][0][te] += p
        pooled[name][1][te] += 1
cv = pd.DataFrame(rows)
cv.to_csv(RES / "classifier_cv.tsv", sep="\t", index=False)
summ = cv.groupby("model")[["auroc", "auprc", "balanced_acc"]].agg(["mean", "std"]).round(3)
print("\nRepeated 5-fold CV (3 repeats), mean and sd over 15 folds:")
print(summ.to_string())
chosen_C = cv.dropna(subset=["C"]).C.median()
print(f"elastic-net C chosen by inner CV: median {chosen_C}, values {sorted(cv.C.dropna().unique())}")

# ---------------- leave-one-site-out ----------------
big = [s for s, n in pd.Series(sites).value_counts().items() if n >= 10]
rows = []
hold = {m: np.full(len(y), np.nan) for m in models}
en_fixed = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=chosen_C, max_iter=3000, tol=1e-3))])
for s in big:
    te, tr = np.where(sites == s)[0], np.where(sites != s)[0]
    for name, model in models.items():
        m = en_fixed if name.startswith("elastic") else model
        p = clone(m).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        hold[name][te] = p
        rows.append({"model": name, "site": s, "n": len(te), "n_msi_h": int(y[te].sum()),
                     "auroc": roc_auc_score(y[te], p) if 0 < y[te].sum() < len(te) else np.nan,
                     "balanced_acc": balanced_accuracy_score(y[te], p >= 0.5)})
loso = pd.DataFrame(rows)
loso.to_csv(RES / "classifier_site_holdout.tsv", sep="\t", index=False)
mask = ~np.isnan(hold["MLH1 alone"])
print(f"\nLeave-one-site-out over {len(big)} sites ({mask.sum()} tumors), pooled held-out predictions:")
for name in models:
    print(f"  {name:<36} AUROC {roc_auc_score(y[mask], hold[name][mask]):.3f}   balanced acc {balanced_accuracy_score(y[mask], hold[name][mask] >= 0.5):.3f}")

# ---------------- what the elastic net uses ----------------
en = clone(en_fixed).fit(X, y)
coef = pd.Series(en.named_steps["lr"].coef_[0], index=genes)
nz = coef[coef != 0].sort_values(key=np.abs, ascending=False)
nz.to_frame("coef").to_csv(RES / "elastic_net_coefficients.tsv", sep="\t")
print(f"\nelastic net (C={chosen_C}) keeps {len(nz)} of {len(genes)} genes; largest coefficients (negative = lower in MSI-H):")
print(nz.head(12).round(3).to_string())

# ---------------- ROC figure ----------------
fig, ax = plt.subplots(figsize=(5.2, 5))
for name, (sp, n) in pooled.items():
    p = sp / n
    fpr, tpr, _ = roc_curve(y, p)
    ax.plot(fpr, tpr, lw=1.5, label=f"{name} (AUROC {roc_auc_score(y, p):.3f})")
ax.plot([0, 1], [0, 1], ls="--", lw=0.7, c="k")
ax.set_xlabel("false positive rate")
ax.set_ylabel("true positive rate")
ax.set_title("MSI-H vs MSS from RNA-seq, pooled cross-validated predictions", fontsize=9)
ax.legend(fontsize=7, loc="lower right")
fig.tight_layout()
fig.savefig(FIG / "roc_msi_classifier.png", dpi=200)
