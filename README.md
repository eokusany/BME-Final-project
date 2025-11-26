# Causal Confidence in Genes Linked with Anthracycline Induced Cardiotoxicity (AIC)

Mendelian Randomization analysis pipeline for assessing causal relationships between gene expression and anthracycline-induced cardiotoxicity.

## Overview

This project implements a two-sample Mendelian Randomization (MR) framework to evaluate causal confidence for genes associated with Anthracycline Induced Cardiotoxicity (AIC). The analysis uses:

- **Exposure Data**: GTEx Heart Left Ventricle eQTL data (gene expression)
- **Outcome Data**: GWAS summary statistics for cardiac traits (LVEF, heart failure, etc.)

## Project Structure

- `run.py` - Main MR analysis script (IVW and MR-Egger)
- `causal_graph.py` - Causal graph discovery using PC algorithm and MR results
- `exposure.txt` - Processed exposure data (eQTL summary statistics)
- `outcome.txt` - Processed outcome data (GWAS summary statistics)

## Data Sources

### Exposure (eQTL)
- **GTEx Portal**: Heart Left Ventricle cis-eQTLs
- Format: rsid, effect_allele, other_allele, beta, se, pval, n

### Outcome (GWAS)
- **IEU OpenGWAS**: LVEF, heart failure, cardiac troponin
- **GWAS Catalog**: Anthracycline cardiotoxicity studies
- Format: rsid, effect_allele, other_allele, beta, se, pval, n

## Usage

### Prerequisites

Ensure you have the required data files:
- `exposure.txt` - Exposure data (eQTL summary statistics)
- `outcome.txt` - Outcome data (GWAS summary statistics)

### 1. Run MR Analysis

```bash
python3 run.py
```

Results are saved to `mr_results.txt` and include:
- IVW MR estimates
- MR-Egger regression (pleiotropy testing)
- SNP harmonization statistics
- Gene-level SNP counts

### 2. Generate Causal Graphs

```bash
python3 causal_graph.py
```

This generates:
- `causal_graph_mr_improved.png` - MR-based causal graph (gene → LVEF edges based on MR significance)
- `causal_graph_pc_improved.png` - PC algorithm causal graph (constraint-based causal discovery)
- `causal_graph_mr_edges_improved.csv` - MR graph edge data with statistics
- `causal_graph_pc_edges_improved.csv` - PC algorithm edge data
- `gene_mr_summary.csv` - Gene-level MR summary statistics

## Requirements

```bash
pip install pandas numpy scipy statsmodels causal-learn
```

## Methods

### Mendelian Randomization
- **Inverse-Variance Weighted (IVW)**: Main causal estimate using meta-analysis of Wald ratios
- **MR-Egger**: Sensitivity analysis for pleiotropy detection
- **Allele Harmonization**: Automatic alignment of exposure/outcome alleles
- **Palindromic SNP Filtering**: Conservative removal of ambiguous SNPs (A/T, C/G)

### Causal Graph Discovery
- **PC Algorithm**: Constraint-based causal discovery using causal-learn library
- **MR-based Graph**: Directed graph using MR significance and effect sizes
- **Bootstrap Sampling**: Creates multiple observations for PC algorithm from SNP-level data

## Output

### MR Analysis (`run.py`)
- Causal effect estimates (beta, SE, p-value)
- Pleiotropy assessment (MR-Egger intercept)
- Harmonization statistics
- Number of SNPs used in analysis
- Gene-level SNP breakdown

### Causal Graphs (`causal_graph.py`)
- MR-based causal graph visualization
- PC algorithm causal structure visualization
- Edge data files with MR statistics
- Gene-level MR summary statistics

## References

- Sanderson E, Glymour MM, Holmes MV, et al. Mendelian Randomization. Nat Rev Methods Primers 2, 6 (2022). [10]
- Burgess S, Davey Smith G, Davies NM, et al. Guidelines for performing Mendelian randomization investigations. Wellcome Open Res. 2023; 4: 186. [16]
- Burgess S, Scott RA, Timpson NJ, Davey Smith G, Thompson SG. Using published data in Mendelian randomization: a blueprint for efficient identification of causal risk. Eur J Epidemiol. 2015;30(7):543-552. [17]
- rondolab. MR-PRESSO (2023). https://github.com/rondolab/MR-PRESSO. [20]
- Burgess S, Thompson SG. Interpreting findings from Mendelian randomization using the MR-Egger method. Eur J Epidemiol. 2017;32(5):377-389. [21]

## Authors

- Pulkit Digani (digani@ualberta.ca)
- Emmanuel Okusanya (eokusany@ualberta.ca)

University of Alberta, Canada
