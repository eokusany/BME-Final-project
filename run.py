#!/usr/bin/env python3
"""
Mendelian Randomization Analysis Script
Implements two-sample MR using IVW and MR-Egger methods

Methods follow Burgess et al. (2023) and Sanderson et al. (2022)
- IVW: Inverse-variance weighted meta-analysis of Wald ratios
- MR-Egger: Weighted regression with intercept to test for pleiotropy
- Allele harmonization: Following TwoSampleMR package approach

Author: Emmanuel Okusanya
Date: 2024
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm

# Open output file for writing results
output_file = "mr_results.txt"
f = open(output_file, 'w')

def print_and_write(*args, **kwargs):
    """Print to console and write to file"""
    message = ' '.join(str(arg) for arg in args)
    print(*args, **kwargs)
    f.write(message + '\n')

print_and_write("=" * 70)
print_and_write("Mendelian Randomization Analysis Results")
print_and_write("=" * 70)
print_and_write()

# ========= 1. LOAD + STANDARDISE =========

exposure_path = "exposure.txt"
outcome_path  = "outcome.txt"

exposure_raw = pd.read_csv(exposure_path, sep="\t")
outcome_raw  = pd.read_csv(outcome_path, sep="\t")

# Rename columns into a standard format
exposure = exposure_raw.rename(columns={
    "rsid": "snp",
    "effect_allele": "ea_exposure",
    "other_allele": "oa_exposure",
    "beta": "beta_exposure",
    "se": "se_exposure",
    "pval": "pval_exposure",
    "n": "n_exposure"
})

# Preserve gene_name if present
has_gene_name = "gene_name" in exposure_raw.columns
if has_gene_name:
    exposure["gene_name"] = exposure_raw["gene_name"]

outcome = outcome_raw.rename(columns={
    "rsid": "snp",
    "effect_allele": "ea_outcome",
    "other_allele": "oa_outcome",
    "beta": "beta_outcome",
    "se": "se_outcome",
    "pval": "pval_outcome",
    "n": "n_outcome"
})

# Keep only relevant columns and drop missing (but allow missing n and gene_name)
exposure_cols = ["snp", "ea_exposure", "oa_exposure",
                 "beta_exposure", "se_exposure", "pval_exposure", "n_exposure"]
if has_gene_name:
    exposure_cols.append("gene_name")
exposure = exposure[exposure_cols].dropna(subset=["snp", "ea_exposure", "oa_exposure", "beta_exposure", "se_exposure", "pval_exposure"])

outcome = outcome[["snp", "ea_outcome", "oa_outcome",
                   "beta_outcome", "se_outcome", "pval_outcome", "n_outcome"]].dropna(subset=["snp", "ea_outcome", "oa_outcome", "beta_outcome", "se_outcome", "pval_outcome"])

print_and_write(f"Exposure SNPs: {exposure.shape[0]}")
print_and_write(f"Outcome SNPs:  {outcome.shape[0]}")

# Report gene information if available
if has_gene_name:
    unique_genes = exposure["gene_name"].nunique()
    print_and_write(f"Unique genes in exposure: {unique_genes}")
    gene_counts = exposure["gene_name"].value_counts()
    print_and_write(f"Genes represented: {', '.join(sorted(gene_counts.index))}")


# ========= 2. MERGE + HARMONISE ALLELES =========
# Allele harmonization is critical for MR - exposure and outcome must use same reference allele
# Following approach from TwoSampleMR R package (Hemani et al. 2018)

merged = exposure.merge(outcome, on="snp", how="inner")
print_and_write(f"Merged SNPs:   {merged.shape[0]}")

# Define helper for palindromic SNPs
# Palindromic SNPs (A/T, C/G) are ambiguous and should be excluded conservatively
# This was causing issues early on - learned from MR tutorials to exclude these
palindromic_pairs = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def is_palindromic(a1, a2):
    return (a1, a2) in palindromic_pairs

def harmonise_row(row):
    ea_exp, oa_exp = row["ea_exposure"], row["oa_exposure"]
    ea_out, oa_out = row["ea_outcome"], row["oa_outcome"]
    beta_out = row["beta_outcome"]

    # 1) Drop palindromic SNPs (conservative)
    if is_palindromic(ea_exp, oa_exp):
        return pd.Series({"keep": False})

    # 2) Alleles already aligned
    if (ea_exp == ea_out) and (oa_exp == oa_out):
        return pd.Series({
            "beta_out_harmonised": beta_out,
            "ea": ea_exp,
            "oa": oa_exp,
            "keep": True
        })

    # 3) Alleles reversed -> flip outcome beta
    if (ea_exp == oa_out) and (oa_exp == ea_out):
        return pd.Series({
            "beta_out_harmonised": -beta_out,
            "ea": ea_exp,
            "oa": oa_exp,
            "keep": True
        })

    # 4) Anything else -> drop
    return pd.Series({"keep": False})


harm = merged.apply(harmonise_row, axis=1)
merged_h = pd.concat([merged, harm], axis=1)
merged_h = merged_h[merged_h["keep"]].copy()

# Use harmonised outcome beta and unified allele labels
merged_h["beta_outcome"] = merged_h["beta_out_harmonised"]

print_and_write(f"SNPs after harmonisation: {merged_h.shape[0]}")


# ========= 3. WALD RATIOS =========
# Each SNP provides an estimate of the causal effect: beta_outcome / beta_exposure
# This is the Wald ratio estimator (standard MR approach)

# Avoid division by almost-zero beta_exposure
# Had to add this after getting NaN values - learned from debugging
min_beta_exp = 1e-6
merged_h = merged_h[np.abs(merged_h["beta_exposure"]) > min_beta_exp].copy()
print_and_write(f"SNPs after beta_exposure filter: {merged_h.shape[0]}")

beta_x = merged_h["beta_exposure"].values
beta_y = merged_h["beta_outcome"].values
se_x   = merged_h["se_exposure"].values
se_y   = merged_h["se_outcome"].values

# Wald ratio for each SNP
beta_mr = beta_y / beta_x

# Delta-method SE for the ratio
# Standard error propagation formula - learned this from stats courses
se_mr = np.sqrt(
    (se_y ** 2) / (beta_x ** 2) +
    (beta_y ** 2) * (se_x ** 2) / (beta_x ** 4)
)

merged_h["beta_mr"] = beta_mr
merged_h["se_mr"] = se_mr


# ========= 4. IVW MR =========
# Inverse-variance weighted meta-analysis combines Wald ratios
# Weights are inverse of variance (1/SE^2) - gives more weight to precise estimates
# This is the standard MR approach (Burgess et al. 2023)

w = 1.0 / (merged_h["se_mr"] ** 2)

beta_ivw = np.sum(w * merged_h["beta_mr"]) / np.sum(w)
se_ivw = np.sqrt(1.0 / np.sum(w))
z_ivw = beta_ivw / se_ivw

# P-value calculation
# Initially was getting p-values showing as 0.0 for very significant results
# Found solution online: use log survival function for extreme z-scores to avoid underflow
# This was a learning moment - didn't realize numerical precision could be an issue!
if np.abs(z_ivw) > 6:
    log_p = stats.norm.logsf(np.abs(z_ivw)) + np.log(2)
    p_ivw = np.exp(log_p)
else:
    p_ivw = 2 * (1 - stats.norm.cdf(np.abs(z_ivw)))

print_and_write()
print_and_write("=== IVW MR result (expression -> outcome) ===")
print_and_write(f"beta_IVW = {beta_ivw:.4f}")
print_and_write(f"SE_IVW   = {se_ivw:.4f}")
print_and_write(f"z        = {z_ivw:.3f}")
print_and_write(f"p-value  = {p_ivw:.3e}")


# ========= 5. MR-EGGER (SENSITIVITY) =========
# MR-Egger regression tests for directional pleiotropy
# Intercept != 0 suggests pleiotropy (Burgess & Thompson 2017) [21]
# Slope provides pleiotropy-adjusted causal estimate
# Using statsmodels WLS (weighted least squares) - learned this from MR tutorials

X = merged_h["beta_exposure"].values
Y = merged_h["beta_outcome"].values
W = 1.0 / (merged_h["se_outcome"].values ** 2)  # Weight by inverse variance of outcome

X_design = sm.add_constant(X)  # adds intercept term
egger_model = sm.WLS(Y, X_design, weights=W).fit()

intercept, slope = egger_model.params
se_intercept, se_slope = egger_model.bse

z_slope = slope / se_slope
# Use log survival function for very small p-values to avoid numerical precision issues
if np.abs(z_slope) > 6:
    log_p = stats.norm.logsf(np.abs(z_slope)) + np.log(2)
    p_slope = np.exp(log_p)
else:
    p_slope = 2 * (1 - stats.norm.cdf(np.abs(z_slope)))

z_intercept = intercept / se_intercept
if np.abs(z_intercept) > 6:
    log_p = stats.norm.logsf(np.abs(z_intercept)) + np.log(2)
    p_intercept = np.exp(log_p)
else:
    p_intercept = 2 * (1 - stats.norm.cdf(np.abs(z_intercept)))

print_and_write()
print_and_write("=== MR-Egger result ===")
print_and_write(f"slope (causal)      = {slope:.4f} (SE = {se_slope:.4f}, p = {p_slope:.3e})")
print_and_write(f"intercept (pleio)   = {intercept:.4f} (SE = {se_intercept:.4f}, p = {p_intercept:.3e})")
print_and_write(f"Number of SNPs used = {merged_h.shape[0]}")

# Report final gene breakdown if available
if has_gene_name and "gene_name" in merged_h.columns:
    final_genes = merged_h["gene_name"].value_counts()
    print_and_write()
    print_and_write("=== Final SNPs by Gene ===")
    for gene, count in sorted(final_genes.items()):
        print_and_write(f"  {gene}: {count} SNPs")

# Close the output file
f.close()
