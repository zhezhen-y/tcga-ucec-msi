"""Step 2: differential expression, MSI vs CN-low endometrioid tumors (PyDESeq2).

Why this contrast: MSI and CN-low tumors are both mostly endometrioid and both copy-number quiet;
the main thing separating them is mismatch-repair deficiency. Restricting to endometrioid
histology keeps histology from confounding the comparison (CN-high tumors are mostly serous).

Outputs (results/): de_msi_vs_cnlow.tsv  (baseMean, log2FC, shrunk log2FC, Wald stat, p, padj)
Figures:  figures/volcano_msi_vs_cnlow.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
FIG.mkdir(exist_ok=True)

cohort = pd.read_csv(RES / "cohort.tsv", sep="\t", index_col=0)
counts = pd.read_csv(RES / "counts.tsv.gz", sep="\t", index_col=0)

meta = cohort[(cohort.subtype.isin(["MSI", "CN_LOW"]))
              & (cohort.histology == "Endometrioid endometrial adenocarcinoma")].copy()
meta["condition"] = meta.subtype.map({"MSI": "MSI", "CN_LOW": "CNlow"})
print(meta.condition.value_counts().to_string())

X = counts[meta.index].T  # PyDESeq2 wants samples x genes
dds = DeseqDataSet(counts=X, metadata=meta[["condition"]], design="~condition",
                   refit_cooks=True, inference=DefaultInference(n_cpus=8), quiet=True)
dds.deseq2()

stats = DeseqStats(dds, contrast=["condition", "MSI", "CNlow"], alpha=0.05, quiet=True)
stats.summary()
res = stats.results_df.copy()
stats.lfc_shrink(coeff="condition[T.MSI]")  # apeGLM-style shrinkage for ranking and plotting
res["log2FC_shrunk"] = stats.results_df["log2FoldChange"]
res = res.sort_values("pvalue")
res.to_csv(RES / "de_msi_vs_cnlow.tsv", sep="\t")

sig = res[(res.padj < 0.05) & (res.log2FC_shrunk.abs() >= 1)]
print(f"genes tested: {res.padj.notna().sum()}  |  padj<0.05 & |shrunk log2FC|>=1: "
      f"{len(sig)} ({(sig.log2FC_shrunk > 0).sum()} up in MSI, {(sig.log2FC_shrunk < 0).sum()} down)")
print("MMR genes:\n" + res.loc[[g for g in ["MLH1", "PMS2", "MSH2", "MSH6", "EPM2AIP1"] if g in res.index],
                               ["baseMean", "log2FC_shrunk", "padj"]].to_string())
print("top 15 by p:\n" + res.head(15)[["baseMean", "log2FC_shrunk", "padj"]].to_string())

# --- Volcano ---
fig, ax = plt.subplots(figsize=(6, 5))
y = -np.log10(res.padj.clip(lower=1e-300))
ax.scatter(res.log2FC_shrunk, y, s=3, c="lightgray", rasterized=True)
ax.scatter(sig.log2FC_shrunk, -np.log10(sig.padj.clip(lower=1e-300)), s=4,
           c=np.where(sig.log2FC_shrunk > 0, "#c0392b", "#2c6fb0"), rasterized=True)
for g in ["MLH1", "EPM2AIP1", "CXCL10", "GBP1", "IDO1", "CXCL9"]:
    if g in res.index:
        ax.annotate(g, (res.at[g, "log2FC_shrunk"], y[g]), fontsize=8,
                    xytext=(4, 2), textcoords="offset points")
ax.axhline(-np.log10(0.05), lw=0.5, ls="--", c="k")
ax.set_xlabel("log2 fold change, MSI vs CN-low (shrunk)")
ax.set_ylabel("-log10 adjusted p")
ax.set_title(f"TCGA-UCEC endometrioid: MSI (n={sum(meta.condition=='MSI')}) "
             f"vs CN-low (n={sum(meta.condition=='CNlow')})", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "volcano_msi_vs_cnlow.png", dpi=200)
