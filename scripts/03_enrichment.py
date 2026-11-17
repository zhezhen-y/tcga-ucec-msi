"""Step 3: gene set enrichment (GSEA) and pathway over-representation (ORA), MSI vs CN-low.

GSEA (preranked): every tested gene, ranked by the DESeq2 Wald statistic, against MSigDB Hallmark
  (50 sets) and Reactome. No significance cutoff needed, so it catches coordinated small shifts.
ORA: hypergeometric test of the significant up- and down-regulated genes separately against
  Reactome and KEGG (MEDICUS), with the tested genes as the background universe. The background
  matters: using the whole genome instead inflates every pathway that is simply expressed in uterus.

Outputs (results/): gsea_hallmark.tsv, gsea_reactome.tsv, ora_up.tsv, ora_down.tsv
Figures:  figures/gsea_hallmark.png, figures/gsea_ifng_running.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA, RES, FIG = ROOT / "data", ROOT / "results", ROOT / "figures"
GMT = {k: str(DATA / f"{k}.v2023.2.Hs.symbols.gmt") for k in ["h.all", "c2.cp.reactome", "c2.cp.kegg_medicus"]}

de = pd.read_csv(RES / "de_msi_vs_cnlow.tsv", sep="\t", index_col=0).dropna(subset=["stat", "padj"])

# ---------------- GSEA preranked on the Wald statistic ----------------
rnk = de["stat"].sort_values(ascending=False)
gsea_out = {}
for name, gmt in [("hallmark", GMT["h.all"]), ("reactome", GMT["c2.cp.reactome"])]:
    pre = gp.prerank(rnk=rnk, gene_sets=gmt, min_size=15, max_size=500,
                     permutation_num=1000, seed=7, threads=8, outdir=None, verbose=False)
    r = pre.res2d.copy()
    r["NES"] = r["NES"].astype(float)
    r["FDR q-val"] = r["FDR q-val"].astype(float)
    r = r.sort_values("NES", ascending=False)
    r.to_csv(RES / f"gsea_{name}.tsv", sep="\t", index=False)
    gsea_out[name] = (pre, r)
    s = r[r["FDR q-val"] < 0.05]
    print(f"GSEA {name}: {len(s)} sets FDR<0.05 ({(s.NES > 0).sum()} up in MSI, {(s.NES < 0).sum()} down)")
    print(s.head(8)[["Term", "NES", "FDR q-val"]].to_string(index=False))
    print(s.tail(5)[["Term", "NES", "FDR q-val"]].to_string(index=False))

# Hallmark bar chart
r = gsea_out["hallmark"][1]
r = r[r["FDR q-val"] < 0.05].copy()
r["Term"] = r.Term.str.replace("HALLMARK_", "").str.replace("_", " ").str.lower()
fig, ax = plt.subplots(figsize=(6, 0.22 * len(r) + 1))
ax.barh(r.Term[::-1], r.NES[::-1], color=np.where(r.NES[::-1] > 0, "#c0392b", "#2c6fb0"))
ax.axvline(0, c="k", lw=0.5)
ax.set_xlabel("normalized enrichment score (MSI vs CN-low)")
ax.set_title("Hallmark gene sets, GSEA FDR < 0.05", fontsize=10)
ax.tick_params(axis="y", labelsize=7)
fig.tight_layout()
fig.savefig(FIG / "gsea_hallmark.png", dpi=200)

# Running-score plot for the interferon-gamma response set
pre = gsea_out["hallmark"][0]
term = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
if term in pre.results:
    gp.gseaplot(rank_metric=pre.ranking, term=term, ofname=str(FIG / "gsea_ifng_running.png"),
                **pre.results[term])

# ---------------- ORA: up and down lists separately ----------------
universe = list(de.index)
sig = de[(de.padj < 0.05) & (de.log2FC_shrunk.abs() >= 1)]
for direction, genes in [("up", sig.index[sig.log2FC_shrunk > 0]), ("down", sig.index[sig.log2FC_shrunk < 0])]:
    frames = []
    for name in ["c2.cp.reactome", "c2.cp.kegg_medicus"]:
        enr = gp.enrich(gene_list=list(genes), gene_sets=GMT[name], background=universe, outdir=None)
        frames.append(enr.results.assign(library=name))
    ora = pd.concat(frames).sort_values("Adjusted P-value")
    ora.to_csv(RES / f"ora_{direction}.tsv", sep="\t", index=False)
    print(f"ORA {direction} ({len(genes)} genes): {(ora['Adjusted P-value'] < 0.05).sum()} pathways adj p<0.05")
    print(ora.head(8)[["Term", "Overlap", "Adjusted P-value"]].to_string(index=False))
