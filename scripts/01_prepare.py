"""Step 1: assemble a tidy cohort table and a raw count matrix for TCGA-UCEC.

Inputs (data/, from the UCSC Xena GDC and PanCanAtlas hubs):
  TCGA-UCEC.star_counts.tsv.gz      STAR gene counts, stored as log2(count + 1)
  gencode.v36.annotation.gtf.gene.probemap   Ensembl ID -> gene symbol
  TCGASubtype.20170308.tsv.gz       PanCanAtlas molecular subtypes
  Survival_SupplementalTable_S1_20171025_xena_sp   TCGA-CDR clinical endpoints (Liu et al. 2018, Cell)

Outputs (results/):
  cohort.tsv     one row per primary tumor: subtype, age, stage, histology, OS and PFI
  counts.tsv.gz  integer counts, expressed genes x cohort samples
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(exist_ok=True)

# --- Molecular subtypes (TCGA 2013 endometrial classes, PanCanAtlas calls) ---
sub = pd.read_csv(DATA / "TCGASubtype.20170308.tsv.gz", sep="\t")
sub = sub[sub.Subtype_Selected.astype(str).str.startswith("UCEC.")]
sub["subtype"] = sub.Subtype_Selected.str.replace("UCEC.", "", regex=False)
sub = sub[sub.subtype != "NA"]
sub["patient"] = sub.sampleID.str[:12]
sub = sub[sub.sampleID.str[13:15] == "01"]  # primary tumor only
sub = sub[["patient", "subtype"]].drop_duplicates("patient")

# --- Clinical endpoints: TCGA-CDR, the curated resource recommended for survival work ---
cdr = pd.read_csv(DATA / "Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t")
cdr = cdr[cdr["cancer type abbreviation"] == "UCEC"].copy()
cdr = cdr.rename(columns={"_PATIENT": "patient",
                          "age_at_initial_pathologic_diagnosis": "age",
                          "clinical_stage": "stage_raw",
                          "histological_type": "histology"})
cdr = cdr.drop_duplicates("patient")


def figo_major(s):
    """'Stage IIIC1' -> 'III'. Endometrial staging in TCGA is FIGO clinical stage."""
    if not isinstance(s, str) or not s.startswith("Stage "):
        return np.nan
    roman = s.split()[1].rstrip("ABC12")
    return roman if roman in {"I", "II", "III", "IV"} else np.nan


cdr["stage"] = cdr.stage_raw.map(figo_major)
cols = ["patient", "age", "stage", "histology", "OS", "OS.time", "PFI", "PFI.time"]
cohort = sub.merge(cdr[cols], on="patient", how="inner")

# --- Counts: undo the hub's log2(x + 1), map Ensembl IDs to gene symbols ---
log_counts = pd.read_csv(DATA / "TCGA-UCEC.star_counts.tsv.gz", sep="\t", index_col=0)
log_counts = log_counts.loc[:, log_counts.columns.str[13:15] == "01"]
log_counts.columns = log_counts.columns.str[:12]
log_counts = log_counts.loc[:, ~log_counts.columns.duplicated()]
counts = np.rint(np.power(2.0, log_counts) - 1).astype(int)

probemap = pd.read_csv(DATA / "gencode.v36.annotation.gtf.gene.probemap", sep="\t", index_col=0)
counts = counts[counts.index.isin(probemap.index)]
counts.index = probemap.loc[counts.index, "gene"]
counts = counts.groupby(level=0).sum()  # collapse PAR_Y duplicates and repeated symbols

# Keep genes expressed (>= 10 counts) in at least 10% of samples; DESeq2 does its own filtering later.
keep_samples = sorted(set(counts.columns) & set(cohort.patient))
counts = counts[keep_samples]
counts = counts[(counts >= 10).sum(axis=1) >= 0.1 * len(keep_samples)]
cohort = cohort[cohort.patient.isin(keep_samples)].set_index("patient").loc[keep_samples]

cohort.to_csv(OUT / "cohort.tsv", sep="\t")
counts.to_csv(OUT / "counts.tsv.gz", sep="\t", compression="gzip")

print(f"cohort: {len(cohort)} tumors with RNA-seq and a molecular subtype")
print(cohort.subtype.value_counts().to_string())
print(f"genes kept: {counts.shape[0]}")
print(cohort[["OS", "PFI"]].sum().rename("events").to_string())
