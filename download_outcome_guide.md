# Downloading GWAS Outcome Data - Step by Step Guide

## Quick Method: Using IEU OpenGWAS Website

### Step 1: Find a Dataset

1. **Visit**: https://gwas.mrcieu.ac.uk/
2. **Search for**: "left ventricular ejection fraction" or "heart failure"
3. **Look for datasets** with IDs like `ieu-a-XXXX` or `ebi-a-XXXX`
4. **Click on a dataset** that looks relevant

### Step 2: Download Summary Statistics

1. On the dataset page, look for a **"Download"** button or link
2. Download the summary statistics file (usually TSV or CSV format)
3. Save it to your BME REPORT folder

### Step 3: Format the Data

Run the formatting script:

```bash
python3 get_real_data.py
```

Or use Python directly:

```python
from get_real_data import format_data_for_mr

format_data_for_mr(
    input_file="your_downloaded_file.tsv",
    output_file="synthetic_outcome_lvef.txt",
    data_type="outcome"
)
```

## Alternative: Using R TwoSampleMR (If R is installed)

If you have R installed, this is often easier:

```r
# Install if needed
install.packages("TwoSampleMR")

# Load library
library(TwoSampleMR)

# Search for available outcomes
ao <- available_outcomes()

# Search for LVEF or heart failure
lvef <- subset(ao, grepl("ejection|heart failure|LVEF|cardiac", trait, ignore.case=TRUE))

# View results
print(lvef[, c("id", "trait", "sample_size", "nsnp")])

# Once you find an ID (e.g., "ieu-a-XXXX"), you can extract data
# First, get your exposure SNPs
exposure_snps <- read.table("synthetic_exposure_eqtl.txt", header=TRUE, sep="\t")$rsid

# Extract outcome data for those SNPs
outcome_dat <- extract_outcome_data(
    snps = exposure_snps,
    outcomes = "ieu-a-XXXX"  # Replace with actual ID
)

# Format and save
write.table(outcome_dat, "synthetic_outcome_lvef.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

## Common IEU OpenGWAS Dataset IDs (Examples - Verify on website)

These are examples - you should verify on the IEU OpenGWAS website:

- **Heart Failure**: Search for "heart failure" 
- **LVEF**: Search for "ejection fraction"
- **Cardiac Traits**: Search for "cardiac" or "cardiovascular"

**Note**: Dataset IDs change, so always verify on the website.

## What to Look For

Your downloaded file should ideally have columns like:
- SNP identifier (rsid, variant_id, etc.)
- Effect allele (EA, A1, alt, etc.)
- Other allele (OA, A2, ref, etc.)
- Beta/effect size
- Standard error (SE)
- P-value
- Sample size (N)

The formatting script will handle common column name variations.

## Troubleshooting

If the formatting script doesn't work:
1. Check that your file has the required columns
2. Open the file and note the exact column names
3. Manually rename columns to match: rsid, effect_allele, other_allele, beta, se, pval, n
4. Save as tab-separated (.txt) file

