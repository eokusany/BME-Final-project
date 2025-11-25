#!/usr/bin/env python3
"""
Process GTEx Heart Left Ventricle eQTL data for Mendelian Randomization
Converts GTEx format to MR-required format
"""

import pandas as pd
import numpy as np

# Curated gene lists for focused analysis (30-40 SNPs)
POSITIVE_CONTROLS = ['SLC28A1', 'SLC28A3', 'RARG', 'CBR3', 'NOS3']

GENE_LIST = [
    'PRDM2', 'WDR4', 'ZNF521', 'SP4', 'RIN3', 'ABCC1', 'ABCC9', 'ABCC5', 
    'HNMT', 'SLC22A17', 'ERCC2', 'MYH7', 'CYP2J2', 'COL1A2', 'SPG7', 
    'GPX3', 'ABCB4', 'GSTP1', 'PLCE1', 'GSTM1', 'CELF4', 'CYBA', 'HAS3', 
    'MLH1', 'POR', 'RAC2', 'CAT', 'ABCC2', 'ATP2B1', 'CBR1', 'ERBB2'
]

NEGATIVE_CONTROLS = ['KRT1', 'KRT10', 'INS', 'AMY1']

HOUSEKEEPING_GENES = ['GAPDH', 'ACTB', 'RPLP0']

# Combine all genes into a single set for filtering
CURATED_GENES = set(POSITIVE_CONTROLS + GENE_LIST + NEGATIVE_CONTROLS + HOUSEKEEPING_GENES)

def process_gtex_egenes(input_file, output_file, pval_threshold=5e-8, min_maf=0.01, filter_genes=None, top_snp_per_gene=True):
    """
    Process GTEx egenes file to MR format
    
    Parameters:
    -----------
    input_file : str
        Path to GTEx egenes file
    output_file : str
        Path to output file (will replace synthetic_exposure_eqtl.txt)
    pval_threshold : float
        P-value threshold for filtering (default: 5e-8 for genome-wide significance)
    min_maf : float
        Minimum minor allele frequency (default: 0.01)
    filter_genes : set or None
        Set of gene names to filter by. If None, uses CURATED_GENES. If empty set, no filtering.
    top_snp_per_gene : bool
        If True, select only the top (lowest p-value) SNP per gene. Default: True
    """
    print("=" * 70)
    print("Processing GTEx Heart Left Ventricle eQTL Data")
    print("=" * 70)
    
    # Read the GTEx file
    print(f"\nReading: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total rows: {df.shape[0]}")
    print(f"Columns: {list(df.columns)}")
    
    # Check for required columns
    required_cols = ['rs_id_dbSNP151_GRCh38p7', 'ref', 'alt', 'slope', 'slope_se', 'pval_nominal']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"\nError: Missing required columns: {missing_cols}")
        return False
    
    # Filter out rows with missing rsid
    print(f"\nFiltering data...")
    initial_count = df.shape[0]
    df = df[df['rs_id_dbSNP151_GRCh38p7'].notna()].copy()
    print(f"  - Removed {initial_count - df.shape[0]} rows with missing rsid")
    
    # Filter by p-value threshold
    df = df[df['pval_nominal'] <= pval_threshold].copy()
    print(f"  - After p-value filter (p <= {pval_threshold}): {df.shape[0]} rows")
    
    # Filter by MAF if available
    if 'maf' in df.columns:
        df = df[df['maf'] >= min_maf].copy()
        print(f"  - After MAF filter (MAF >= {min_maf}): {df.shape[0]} rows")
    
    # Filter by curated gene list if specified
    if filter_genes is None:
        filter_genes = CURATED_GENES
    
    if filter_genes and 'gene_name' in df.columns:
        initial_gene_count = df.shape[0]
        df = df[df['gene_name'].isin(filter_genes)].copy()
        print(f"  - After gene filter ({len(filter_genes)} genes): {df.shape[0]} rows")
        if initial_gene_count > df.shape[0]:
            print(f"  - Removed {initial_gene_count - df.shape[0]} rows not matching curated genes")
        
        # Report which genes were found
        found_genes = df['gene_name'].unique()
        missing_genes = filter_genes - set(found_genes)
        print(f"  - Found SNPs for {len(found_genes)} genes: {sorted(found_genes)}")
        if missing_genes:
            print(f"  - Warning: No SNPs found for {len(missing_genes)} genes: {sorted(missing_genes)}")
    
    # Create output dataframe with required columns
    output_df = pd.DataFrame()
    
    # Map columns
    output_df['rsid'] = df['rs_id_dbSNP151_GRCh38p7'].astype(str)
    
    # Add gene_name if available (useful for tracking)
    if 'gene_name' in df.columns:
        output_df['gene_name'] = df['gene_name'].astype(str)
    
    # Determine effect allele and other allele
    # In GTEx, 'alt' is typically the effect allele, 'ref' is the other allele
    # But we need to check the slope direction
    # For now, we'll use alt as effect_allele and ref as other_allele
    # The harmonization step in run.py will handle any mismatches
    output_df['effect_allele'] = df['alt'].astype(str)
    output_df['other_allele'] = df['ref'].astype(str)
    
    # Beta and SE
    output_df['beta'] = df['slope'].astype(float)
    output_df['se'] = df['slope_se'].astype(float)
    output_df['pval'] = df['pval_nominal'].astype(float)
    
    # Sample size - GTEx v8 typically has ~700 samples for heart tissue
    # We can estimate from minor_allele_samples if available, or use a constant
    if 'minor_allele_samples' in df.columns and 'maf' in df.columns:
        # Estimate sample size: n = minor_allele_samples / (2 * maf) for heterozygotes
        # More accurate: use the total number of samples
        # GTEx v8 Heart Left Ventricle has ~387 samples
        # But we can try to estimate from the data
        estimated_n = df['minor_allele_samples'].max() / df['maf'].max() if df['maf'].max() > 0 else 387
        output_df['n'] = int(estimated_n)
        print(f"  - Estimated sample size: {int(estimated_n)}")
    else:
        # GTEx v8 Heart Left Ventricle typical sample size
        output_df['n'] = 387
        print(f"  - Using default sample size: 387 (GTEx v8 Heart Left Ventricle)")
    
    # Remove any rows with invalid data
    output_df = output_df[
        (output_df['beta'].notna()) & 
        (output_df['se'].notna()) & 
        (output_df['se'] > 0) &
        (output_df['pval'].notna()) &
        (output_df['effect_allele'].isin(['A', 'T', 'G', 'C'])) &
        (output_df['other_allele'].isin(['A', 'T', 'G', 'C']))
    ].copy()
    
    print(f"  - After data quality filter: {output_df.shape[0]} rows")
    
    # Remove duplicates (keep first occurrence if same SNP appears multiple times)
    initial_dup = output_df.shape[0]
    output_df = output_df.drop_duplicates(subset=['rsid'], keep='first')
    if initial_dup > output_df.shape[0]:
        print(f"  - Removed {initial_dup - output_df.shape[0]} duplicate SNPs")
    
    # If gene_name is available and top_snp_per_gene is True, select top SNP per gene
    if top_snp_per_gene and 'gene_name' in output_df.columns:
        initial_count = output_df.shape[0]
        # Sort by p-value, then take first (lowest p-value) per gene
        output_df = output_df.sort_values('pval').groupby('gene_name').first().reset_index()
        if initial_count > output_df.shape[0]:
            print(f"  - Selected top SNP per gene: {output_df.shape[0]} SNPs from {output_df['gene_name'].nunique()} genes")
    
    # Sort by p-value
    output_df = output_df.sort_values('pval')
    
    # Select and reorder columns (include gene_name if available)
    base_cols = ['rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'n']
    if 'gene_name' in output_df.columns:
        base_cols.insert(1, 'gene_name')  # Insert after rsid
    output_df = output_df[base_cols]
    
    # Save to file
    output_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✓ Processed data saved to: {output_file}")
    print(f"  Final SNPs: {output_df.shape[0]}")
    
    # Report gene breakdown if gene_name is available
    if 'gene_name' in output_df.columns:
        gene_counts = output_df['gene_name'].value_counts()
        print(f"\nSNPs per gene:")
        for gene, count in gene_counts.items():
            print(f"  - {gene}: {count} SNPs")
    
    print(f"\nSummary statistics:")
    print(f"  - P-value range: {output_df['pval'].min():.2e} to {output_df['pval'].max():.2e}")
    print(f"  - Beta range: {output_df['beta'].min():.4f} to {output_df['beta'].max():.4f}")
    print(f"  - Mean SE: {output_df['se'].mean():.4f}")
    
    print(f"\nFirst few rows:")
    print(output_df.head(10).to_string())
    
    return True


if __name__ == "__main__":
    input_file = "Heart_Left_Ventricle.v8.egenes.txt"
    output_file = "synthetic_exposure_eqtl.txt"  # Replace synthetic file
    
    print("\nNote: This will replace the synthetic exposure data file.")
    print("Make sure you have a backup if needed.\n")
    print(f"Filtering by curated gene list ({len(CURATED_GENES)} genes):")
    print(f"  - Positive controls: {len(POSITIVE_CONTROLS)} genes")
    print(f"  - Main gene list: {len(GENE_LIST)} genes")
    print(f"  - Negative controls: {len(NEGATIVE_CONTROLS)} genes")
    print(f"  - Housekeeping genes: {len(HOUSEKEEPING_GENES)} genes")
    print()
    
    success = process_gtex_egenes(
        input_file=input_file,
        output_file=output_file,
        pval_threshold=5e-4,  # More lenient threshold to get 30-40 SNPs
        min_maf=0.01,  # Minimum MAF of 1%
        filter_genes=CURATED_GENES,  # Filter by curated genes
        top_snp_per_gene=True  # Select top SNP per gene for better coverage
    )
    
    if success:
        print("\n" + "=" * 70)
        print("Processing complete!")
        print("=" * 70)
        print(f"\nNext steps:")
        print(f"1. Verify the data looks correct")
        print(f"2. Make sure you have outcome (GWAS) data ready")
        print(f"3. Run: python3 run.py")
    else:
        print("\nProcessing failed. Please check the error messages above.")

