# Data Download Instructions for AIC Mendelian Randomization

This guide explains how to replace the synthetic data files with real eQTL and GWAS data.

## Required Data Format

Both exposure and outcome files must have these columns (tab-separated):
- `rsid`: SNP identifier (e.g., rs12345)
- `effect_allele`: Effect allele (A, T, G, or C)
- `other_allele`: Other allele (A, T, G, or C)
- `beta`: Effect size
- `se`: Standard error
- `pval`: P-value
- `n`: Sample size

## Step 1: Download Exposure Data (eQTL - Gene Expression)

### Option A: GTEx Portal (Recommended for Heart Tissue)

1. Visit: https://www.gtexportal.org/home/eqtlDashboardPage
2. Select tissue: **Heart - Left Ventricle** or **Heart - Atrial Appendage**
3. For specific genes: Enter gene symbol (e.g., from your 80 AIC-linked genes)
4. For all cis-eQTLs: Download bulk cis-eQTL data
5. Save as: `exposure_eqtl_gtex.txt`

**GTEx Data Format:**
- Usually includes: variant_id, gene_id, tss_distance, ma_samples, ma_count, maf, pval_nominal, slope, slope_se
- You may need to convert variant_id to rsid using dbSNP
- Use `slope` as `beta` and `slope_se` as `se`

### Option B: eQTLGen (Blood eQTLs)

1. Visit: https://www.eqtlgen.org/cis-eqtls.html
2. Download: "cis-eQTLs (all significant)"
3. Filter for your genes of interest
4. Save as: `exposure_eqtl_eqtlgen.txt`

**eQTLGen Format:**
- Usually includes: Pvalue, FDR, NrSamples, Gene, GeneSymbol, etc.
- May need column renaming

## Step 2: Download Outcome Data (GWAS - AIC/Cardiac Traits)

### Option A: IEU OpenGWAS (Recommended)

1. Visit: https://gwas.mrcieu.ac.uk/
2. Search for one of these traits:
   - "left ventricular ejection fraction" (LVEF)
   - "heart failure"
   - "cardiac troponin"
   - "QT interval"
   - "anthracycline" (if available)
3. Note the GWAS ID (format: `ieu-a-XXXX` or `ebi-a-XXXX`)
4. Click "Download" to get summary statistics
5. Save as: `outcome_gwas_ieu.txt`

**Common IEU OpenGWAS IDs:**
- LVEF: Search for "ejection fraction"
- Heart Failure: Search for "heart failure"
- Note: You may need to search and find the most relevant dataset

### Option B: GWAS Catalog

1. Visit: https://www.ebi.ac.uk/gwas/
2. Search: "anthracycline cardiotoxicity" or related cardiac traits
3. Download summary statistics if available
4. Save as: `outcome_gwas_catalog.txt`

## Step 3: Format Your Data

Use the provided script to format your downloaded data:

```bash
python3 get_real_data.py
```

Or format manually using Python:

```python
from get_real_data import format_data_for_mr

# Format exposure data
format_data_for_mr(
    input_file="exposure_eqtl_gtex.txt",
    output_file="synthetic_exposure_eqtl.txt",  # Replaces synthetic
    data_type="exposure"
)

# Format outcome data
format_data_for_mr(
    input_file="outcome_gwas_ieu.txt",
    output_file="synthetic_outcome_lvef.txt",  # Replaces synthetic
    data_type="outcome"
)
```

## Step 4: Verify Data Format

Check that your files have the correct format:

```bash
head -5 synthetic_exposure_eqtl.txt
head -5 synthetic_outcome_lvef.txt
```

Should show:
```
rsid	effect_allele	other_allele	beta	se	pval	n
rs12345	A	G	0.15	0.02	1e-10	500
...
```

## Step 5: Run Analysis

Once data is formatted and in place:

```bash
python3 run.py
```

## Troubleshooting

### Missing rsid column
- If you have variant IDs (chr:pos), convert to rsid using:
  - dbSNP VCF files
  - `bcftools annotate`
  - Online tools like Ensembl VEP

### Different column names
- The `format_data_for_mr()` function handles common column name variations
- If it doesn't work, manually rename columns to match required format

### Missing sample size (n)
- If sample size is not in the data, you may need to:
  - Add it manually if you know it
  - Leave as NaN (some analyses can proceed without it)

### Allele orientation
- Ensure alleles are on the same strand (forward/reverse)
- The harmonization step in `run.py` will handle some mismatches

## Notes

- **Cis-eQTLs only**: For exposure data, focus on cis-eQTLs (variants near the gene)
- **LD Clumping**: You may want to perform LD clumping on exposure SNPs before MR
- **Multiple genes**: If analyzing multiple genes, you may need separate analyses or combine data appropriately
- **Sample overlap**: Check if exposure and outcome datasets have sample overlap (affects MR assumptions)

## Quick Start (If You Have Data Ready)

If you already have formatted data files:

```bash
# Copy your files to replace synthetic data
cp your_exposure_data.txt synthetic_exposure_eqtl.txt
cp your_outcome_data.txt synthetic_outcome_lvef.txt

# Run analysis
python3 run.py
```

