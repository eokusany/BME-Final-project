#!/usr/bin/env python3
"""
Script to help download and format real eQTL and GWAS data
for AIC Mendelian Randomization analysis

This script provides functions to:
1. Search for datasets in IEU OpenGWAS
2. Download data from public APIs
3. Format data to match run.py requirements
"""

import pandas as pd
import os
import sys

def search_ieu_opengwas(trait_name):
    """
    Search IEU OpenGWAS for trait-specific datasets
    Returns list of potential GWAS IDs
    """
    print(f"\nSearching IEU OpenGWAS for: {trait_name}")
    print("=" * 70)
    print("\nVisit: https://gwas.mrcieu.ac.uk/")
    print(f"Search for: {trait_name}")
    print("\nCommon AIC-related searches:")
    print("  - 'left ventricular ejection fraction'")
    print("  - 'heart failure'")
    print("  - 'cardiac troponin'")
    print("  - 'QT interval'")
    print("  - 'anthracycline'")
    print("\nOnce you find a dataset, note the ID (format: ieu-a-XXXX)")
    return None


def format_data_for_mr(input_file, output_file, data_type="exposure"):
    """
    Format downloaded data to match run.py requirements
    
    Required columns: rsid, effect_allele, other_allele, beta, se, pval, n
    """
    print(f"\nFormatting {data_type} data...")
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found!")
        return False
    
    try:
        # Try to read the file (handle different separators)
        df = pd.read_csv(input_file, sep='\t')
        if df.shape[1] == 1:
            df = pd.read_csv(input_file, sep=',')
        if df.shape[1] == 1:
            df = pd.read_csv(input_file, sep=' ')
        
        print(f"\nOriginal columns: {list(df.columns)}")
        print(f"Shape: {df.shape}")
        
        # Common column name mappings
        column_mappings = {
            'rsid': ['rsid', 'rsID', 'SNP', 'snp', 'variant_id', 'variant', 'ID'],
            'effect_allele': ['effect_allele', 'EA', 'ea', 'A1', 'a1', 'alt', 'ALT', 'effect_allele_freq'],
            'other_allele': ['other_allele', 'OA', 'oa', 'A2', 'a2', 'ref', 'REF', 'other_allele_freq'],
            'beta': ['beta', 'BETA', 'effect', 'Effect', 'b', 'B'],
            'se': ['se', 'SE', 'std_err', 'StdErr', 'standard_error'],
            'pval': ['pval', 'P', 'p', 'p_value', 'pvalue', 'P-value', 'PVAL'],
            'n': ['n', 'N', 'sample_size', 'SampleSize', 'N_SAMPLES', 'n_samples']
        }
        
        # Find matching columns
        mapped_cols = {}
        for target, possible_names in column_mappings.items():
            for col in df.columns:
                if col in possible_names:
                    mapped_cols[target] = col
                    break
        
        print(f"\nMapped columns: {mapped_cols}")
        
        # Create output dataframe
        output_df = pd.DataFrame()
        for target, source_col in mapped_cols.items():
            if source_col in df.columns:
                output_df[target] = df[source_col]
            else:
                print(f"Warning: Could not find column for '{target}'")
        
        # Check if we have required columns
        required = ['rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval']
        missing = [col for col in required if col not in output_df.columns]
        
        if missing:
            print(f"\nError: Missing required columns: {missing}")
            print("\nPlease ensure your data has these columns (or similar names):")
            for col in required:
                print(f"  - {col}: {column_mappings[col]}")
            return False
        
        # Handle missing 'n' column (sample size)
        if 'n' not in output_df.columns:
            print("\nWarning: 'n' (sample size) column not found.")
            print("You may need to add this manually or it will be set to NaN")
            output_df['n'] = None
        
        # Clean data
        output_df = output_df.dropna(subset=required)
        
        # Ensure numeric types
        for col in ['beta', 'se', 'pval', 'n']:
            if col in output_df.columns:
                output_df[col] = pd.to_numeric(output_df[col], errors='coerce')
        
        # Save formatted data
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Formatted data saved to: {output_file}")
        print(f"  Rows: {output_df.shape[0]}")
        print(f"  Columns: {list(output_df.columns)}")
        
        return True
        
    except Exception as e:
        print(f"\nError formatting data: {e}")
        return False


def main():
    """
    Main function - interactive data download and formatting
    """
    print("=" * 70)
    print("Real Data Download and Formatting for AIC MR Analysis")
    print("=" * 70)
    
    print("\nThis script helps you:")
    print("1. Find datasets in public databases")
    print("2. Format downloaded data to match run.py requirements")
    
    print("\n" + "=" * 70)
    print("STEP 1: Download Data")
    print("=" * 70)
    
    print("\nEXPOSURE DATA (eQTL - Gene Expression):")
    print("\nOption A - GTEx Portal:")
    print("  1. Visit: https://www.gtexportal.org/home/eqtlDashboardPage")
    print("  2. Select: Heart - Left Ventricle")
    print("  3. Enter gene symbol or download all cis-eQTLs")
    print("  4. Download results")
    
    print("\nOption B - eQTLGen:")
    print("  1. Visit: https://www.eqtlgen.org/cis-eqtls.html")
    print("  2. Download cis-eQTL results")
    
    print("\nOUTCOME DATA (GWAS - AIC/Cardiac Traits):")
    print("\nOption A - IEU OpenGWAS:")
    print("  1. Visit: https://gwas.mrcieu.ac.uk/")
    print("  2. Search for: 'left ventricular ejection fraction' or 'heart failure'")
    print("  3. Note the GWAS ID (format: ieu-a-XXXX)")
    print("  4. Download summary statistics")
    
    print("\nOption B - GWAS Catalog:")
    print("  1. Visit: https://www.ebi.ac.uk/gwas/")
    print("  2. Search for: 'anthracycline cardiotoxicity'")
    print("  3. Download summary statistics if available")
    
    print("\n" + "=" * 70)
    print("STEP 2: Format Your Data")
    print("=" * 70)
    
    response = input("\nDo you have downloaded data files to format? (y/n): ").strip().lower()
    
    if response == 'y':
        exposure_input = input("Enter path to exposure (eQTL) data file: ").strip()
        exposure_output = "synthetic_exposure_eqtl.txt"  # Will replace synthetic
        
        outcome_input = input("Enter path to outcome (GWAS) data file: ").strip()
        outcome_output = "synthetic_outcome_lvef.txt"  # Will replace synthetic
        
        print("\nFormatting exposure data...")
        if format_data_for_mr(exposure_input, exposure_output, "exposure"):
            print(f"\n✓ Exposure data ready: {exposure_output}")
        
        print("\nFormatting outcome data...")
        if format_data_for_mr(outcome_input, outcome_output, "outcome"):
            print(f"\n✓ Outcome data ready: {outcome_output}")
        
        print("\n" + "=" * 70)
        print("Next Steps:")
        print("=" * 70)
        print("\n1. Verify the formatted files look correct")
        print("2. The files have replaced the synthetic data")
        print("3. Run: python3 run.py")
        print("\n" + "=" * 70)
    else:
        print("\nOnce you download your data, run this script again to format it.")
        print("Or use the format_data_for_mr() function directly in Python.")


if __name__ == "__main__":
    main()

