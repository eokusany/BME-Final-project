# Causal Confidence in Genes Linked with Anthracycline Induced Cardiotoxicity (AIC)

Mendelian Randomization analysis pipeline for assessing causal relationships between gene expression and anthracycline-induced cardiotoxicity.

## Overview

This project implements a two-sample Mendelian Randomization (MR) framework to evaluate causal confidence for genes associated with Anthracycline Induced Cardiotoxicity (AIC). The analysis uses:

- **Exposure Data**: GTEx Heart Left Ventricle eQTL data (gene expression)
- **Outcome Data**: GWAS summary statistics for cardiac traits (LVEF, heart failure, etc.)

## Project Structure

- `run.py` - Main MR analysis script (IVW and MR-Egger)
- `process_gtex_data.py` - Processes GTEx eQTL data for MR
- `get_real_data.py` - Helper script to format downloaded GWAS/eQTL data
- `download_data.py` - Guide for downloading data from public databases
- `causal_learn_example.py` - Causal discovery using causal-learn library

## Data Sources

### Exposure (eQTL)
- **GTEx Portal**: Heart Left Ventricle cis-eQTLs
- Format: rsid, effect_allele, other_allele, beta, se, pval, n

### Outcome (GWAS)
- **IEU OpenGWAS**: LVEF, heart failure, cardiac troponin
- **GWAS Catalog**: Anthracycline cardiotoxicity studies
- Format: rsid, effect_allele, other_allele, beta, se, pval, n

## Usage

### 1. Prepare Data

Process GTEx exposure data:
```bash
python3 process_gtex_data.py
```

Format outcome GWAS data:
```bash
python3 get_real_data.py
```

### 2. Run MR Analysis

```bash
python3 run.py
```

Results are saved to `mr_results.txt` and include:
- IVW MR estimates
- MR-Egger regression (pleiotropy testing)
- SNP harmonization statistics

## Requirements

```bash
pip install pandas numpy scipy statsmodels causal-learn
```

## Methods

- **Inverse-Variance Weighted (IVW)**: Main causal estimate
- **MR-Egger**: Sensitivity analysis for pleiotropy
- **Allele Harmonization**: Automatic alignment of exposure/outcome alleles
- **Palindromic SNP Filtering**: Conservative removal of ambiguous SNPs

## Output

The analysis produces:
- Causal effect estimates (beta, SE, p-value)
- Pleiotropy assessment (MR-Egger intercept)
- Harmonization statistics
- Number of SNPs used in analysis

## References

- Burgess S, Scott RA, Timpson NJ, et al. Using published data in Mendelian randomization: a blueprint for efficient identification of causal risk factors. Eur J Epidemiol. 2015; 30(7): 543–552.
- Sanderson, E., Glymour, M.M., Holmes, M.V. et al. Mendelian randomization. Nat Rev Methods Primers 2, 6 (2022).

## Authors

- Pulkit Digani (digani@ualberta.ca)
- Emmanuel Okusanya (eokusany@ualberta.ca)

University of Alberta, Canada
