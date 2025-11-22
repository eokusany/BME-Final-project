#!/usr/bin/env python3
"""
Manual instructions and helper for downloading outcome GWAS data
Since automatic download may be limited, this provides clear steps
"""

import pandas as pd
import os

def create_manual_instructions():
    """
    Create a guide for manually downloading outcome data
    """
    print("=" * 70)
    print("Downloading GWAS Outcome Data - Manual Guide")
    print("=" * 70)
    
    print("\nSince automatic API access may be limited, here are the best options:")
    
    print("\n" + "=" * 70)
    print("OPTION 1: IEU OpenGWAS (Recommended)")
    print("=" * 70)
    print("\n1. Visit: https://gwas.mrcieu.ac.uk/")
    print("2. Search for: 'left ventricular ejection fraction'")
    print("   OR search for: 'heart failure'")
    print("3. Look for datasets with format: ieu-a-XXXX or ebi-a-XXXX")
    print("4. Click on a dataset")
    print("5. Click 'Download' to get summary statistics")
    print("6. Save the file")
    
    print("\n" + "=" * 70)
    print("OPTION 2: Use R TwoSampleMR Package")
    print("=" * 70)
    print("\nIf you have R installed, you can use:")
    print("""
# Install if needed
install.packages("TwoSampleMR")

# Load library
library(TwoSampleMR)

# Search for available outcomes
ao <- available_outcomes()

# Search for LVEF or heart failure
lvef <- subset(ao, grepl("ejection|heart failure|LVEF", trait, ignore.case=TRUE))

# View results
print(lvef[, c("id", "trait", "sample_size")])

# Download data for a specific ID (replace with actual ID)
# dat <- extract_outcome_data(snps=exposure_snps, outcomes="ieu-a-XXXX")
""")
    
    print("\n" + "=" * 70)
    print("OPTION 3: GWAS Catalog")
    print("=" * 70)
    print("\n1. Visit: https://www.ebi.ac.uk/gwas/")
    print("2. Search for: 'left ventricular ejection fraction'")
    print("   OR: 'heart failure'")
    print("   OR: 'anthracycline cardiotoxicity'")
    print("3. Find studies and check if summary statistics are available")
    print("4. Download if available")
    
    print("\n" + "=" * 70)
    print("OPTION 4: Use Pre-formatted Example")
    print("=" * 70)
    print("\nFor testing purposes, I can create a formatted template")
    print("that you can replace with real data later.")
    
    response = input("\nWould you like me to create a template file for now? (y/n): ").strip().lower()
    
    if response == 'y':
        create_template_outcome()
    else:
        print("\nOnce you download your data, use get_real_data.py to format it.")


def create_template_outcome():
    """
    Create a template outcome file with instructions
    """
    print("\nCreating template outcome file...")
    
    # Create a small example dataset that matches the format
    # This is just for structure - user should replace with real data
    template_data = {
        'rsid': ['rs12345', 'rs67890', 'rs11111'],
        'effect_allele': ['A', 'G', 'T'],
        'other_allele': ['G', 'C', 'C'],
        'beta': [0.1, -0.15, 0.08],
        'se': [0.02, 0.03, 0.015],
        'pval': [1e-8, 1e-6, 1e-7],
        'n': [10000, 10000, 10000]
    }
    
    df = pd.DataFrame(template_data)
    output_file = "synthetic_outcome_lvef.txt"
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✓ Created template: {output_file}")
    print("\n⚠️  WARNING: This is a TEMPLATE with example data.")
    print("   You MUST replace this with real GWAS data before running MR analysis!")
    print("\nThe file has the correct format - just replace the data with your real GWAS results.")


def try_download_known_datasets():
    """
    Try to download from known IEU OpenGWAS dataset IDs
    """
    # Some potential IDs (these may need to be verified)
    # These are examples - actual IDs need to be found on the website
    known_ids = [
        # These are placeholders - need to find actual IDs
    ]
    
    print("\nTo find actual dataset IDs:")
    print("1. Go to https://gwas.mrcieu.ac.uk/")
    print("2. Search and note the ID from the URL or dataset page")
    print("3. Use that ID with download scripts")


if __name__ == "__main__":
    create_manual_instructions()

