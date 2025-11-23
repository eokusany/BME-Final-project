#!/usr/bin/env python3
"""
Process VCF file from IEU OpenGWAS to MR format
Converts VCF format to required MR format: rsid, effect_allele, other_allele, beta, se, pval, n
"""

import gzip
import pandas as pd
import sys
import os

def read_vcf_header(vcf_file):
    """
    Read VCF header to understand the format
    """
    header_lines = []
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line.strip())
            elif line.startswith('#CHROM'):
                # This is the column header line
                columns = line.strip().split('\t')
                return header_lines, columns
    return header_lines, None


def process_vcf_to_mr(vcf_file, output_file="synthetic_outcome_lvef.txt"):
    """
    Convert VCF file to MR format
    VCF files from IEU OpenGWAS may have summary statistics in INFO field
    """
    print("=" * 70)
    print("Processing VCF file to MR format")
    print("=" * 70)
    print(f"\nInput: {vcf_file}")
    print(f"Output: {output_file}")
    
    if not os.path.exists(vcf_file):
        print(f"\nError: File '{vcf_file}' not found!")
        return False
    
    # Read header to understand format
    print("\nReading VCF header...")
    header_lines, columns = read_vcf_header(vcf_file)
    
    if columns:
        print(f"Columns: {columns[:10]}...")  # Show first 10 columns
    
    # Check if this is a standard VCF or has summary stats
    # IEU OpenGWAS VCF files might have summary stats in INFO field
    print("\nReading VCF data...")
    
    data_rows = []
    with gzip.open(vcf_file, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            
            if i % 100000 == 0 and i > 0:
                print(f"  Processed {i:,} variants...")
            
            fields = line.strip().split('\t')
            
            if len(fields) < 8:  # VCF should have at least 8 standard fields
                continue
            
            # Standard VCF fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, DATA
            if len(fields) < 10:
                continue
                
            chrom = fields[0]
            pos = fields[1]
            variant_id = fields[2]  # Usually rsid or chr:pos format
            ref = fields[3]
            alt = fields[4]
            info = fields[7] if len(fields) > 7 else ""
            format_field = fields[8] if len(fields) > 8 else ""
            data_field = fields[9] if len(fields) > 9 else ""
            
            # Parse FORMAT and DATA fields (not INFO)
            # FORMAT: "ES:SE:LP:AF:SI:ID"
            # DATA: "-0.088471:0.0838036:0.537602:0.599372:0.466089:rs1264289758"
            # ES = Effect size (beta)
            # SE = Standard error
            # LP = -log10(p-value), so pval = 10^(-LP)
            # SS = Sample size (may be in separate field or INFO)
            
            beta = None
            se = None
            pval = None
            n = None
            
            if format_field and data_field:
                format_keys = format_field.split(':')
                data_values = data_field.split(':')
                
                # Create dictionary mapping format keys to values
                format_dict = dict(zip(format_keys, data_values))
                
                # Extract values
                if 'ES' in format_dict and format_dict['ES'] != '.':
                    try:
                        beta = float(format_dict['ES'])
                    except:
                        pass
                
                if 'SE' in format_dict and format_dict['SE'] != '.':
                    try:
                        se = float(format_dict['SE'])
                    except:
                        pass
                
                if 'LP' in format_dict and format_dict['LP'] != '.':
                    try:
                        lp = float(format_dict['LP'])  # -log10(p-value)
                        pval = 10 ** (-lp)  # Convert back to p-value
                    except:
                        pass
                
                if 'SS' in format_dict and format_dict['SS'] != '.':
                    try:
                        n = float(format_dict['SS'])
                    except:
                        pass
            
            # Also check INFO field for sample size or other metadata
            info_dict = {}
            if info:
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                        if key in ['SS', 'N', 'sample_size'] and not n:
                            try:
                                n = float(value)
                            except:
                                pass
            
            # Only include rows with essential data
            if variant_id and ref and alt:
                # If we don't have beta/se/pval in INFO, we might need to check FORMAT/GT fields
                # But for summary stats VCF, they should be in INFO
                
                data_rows.append({
                    'rsid': variant_id if variant_id.startswith('rs') else f"{chrom}:{pos}",
                    'effect_allele': alt.split(',')[0],  # Take first ALT allele
                    'other_allele': ref,
                    'beta': beta,
                    'se': se,
                    'pval': pval,
                    'n': n,
                    'info': info  # Keep for debugging
                })
    
    print(f"\nTotal variants read: {len(data_rows):,}")
    
    # Convert to DataFrame
    df = pd.DataFrame(data_rows)
    
    # Filter rows with required data
    print("\nFiltering data...")
    initial_count = len(df)
    
    # We need at least rsid, alleles, and one of beta/se/pval
    df = df[df['rsid'].notna()].copy()
    df = df[df['effect_allele'].notna()].copy()
    df = df[df['other_allele'].notna()].copy()
    
    # Try to keep rows with at least p-value (most important)
    df_with_pval = df[df['pval'].notna()].copy()
    
    if len(df_with_pval) > 0:
        df = df_with_pval
        print(f"  - Kept {len(df):,} variants with p-values")
    else:
        print("  - Warning: No p-values found in INFO field")
        print("  - Checking if data is in different format...")
    
    # If we still don't have beta/se, we might need to estimate or use placeholder
    # But let's see what we have first
    print(f"\nData summary:")
    print(f"  - Variants with beta: {df['beta'].notna().sum():,}")
    print(f"  - Variants with SE: {df['se'].notna().sum():,}")
    print(f"  - Variants with pval: {df['pval'].notna().sum():,}")
    print(f"  - Variants with sample size: {df['n'].notna().sum():,}")
    
    # Show sample of INFO field to understand format
    if len(df) > 0:
        print(f"\nSample INFO fields (first 3):")
        for i, info in enumerate(df['info'].head(3)):
            print(f"  {i+1}. {info[:200]}...")  # First 200 chars
    
    # Prepare output
    output_df = pd.DataFrame()
    output_df['rsid'] = df['rsid']
    output_df['effect_allele'] = df['effect_allele']
    output_df['other_allele'] = df['other_allele']
    output_df['beta'] = df['beta']
    output_df['se'] = df['se']
    output_df['pval'] = df['pval']
    output_df['n'] = df['n']
    
    # Remove rows missing critical data
    # At minimum, we need rsid, alleles, and pval
    output_df = output_df[
        (output_df['rsid'].notna()) &
        (output_df['effect_allele'].notna()) &
        (output_df['other_allele'].notna()) &
        (output_df['pval'].notna())
    ].copy()
    
    print(f"\nFinal variants for MR: {len(output_df):,}")
    
    if len(output_df) == 0:
        print("\nError: No valid data extracted from VCF file!")
        print("\nThe VCF file might be in a different format than expected.")
        print("Please check the file format or try downloading summary statistics")
        print("in TSV/CSV format instead of VCF.")
        return False
    
    # Save
    output_df = output_df[['rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'n']]
    output_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✓ Saved to: {output_file}")
    print(f"\nFirst few rows:")
    print(output_df.head(10).to_string())
    
    return True


if __name__ == "__main__":
    vcf_file = "ieu-b-5120.vcf.gz"
    output_file = "synthetic_outcome_lvef.txt"
    
    if len(sys.argv) > 1:
        vcf_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    success = process_vcf_to_mr(vcf_file, output_file)
    
    if success:
        print("\n" + "=" * 70)
        print("SUCCESS!")
        print("=" * 70)
        print(f"\nOutcome data ready: {output_file}")
        print("You can now run: python3 run.py")
    else:
        print("\nProcessing failed. Please check the error messages above.")

