#!/usr/bin/env python3
"""
Download GWAS outcome data from IEU OpenGWAS
Searches for and downloads LVEF, heart failure, or related cardiac traits
"""

import pandas as pd
import requests
import json
import os
import time

def search_ieu_opengwas(query):
    """
    Search IEU OpenGWAS for datasets
    """
    print(f"Searching IEU OpenGWAS for: {query}")
    
    # IEU OpenGWAS API endpoint
    api_url = "https://api.gwas-api.mrcieu.ac.uk/gwasinfo/search"
    
    try:
        response = requests.get(api_url, params={"query": query}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            return data
        else:
            print(f"API returned status code: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error searching API: {e}")
        return None


def get_gwas_data(gwas_id, output_file):
    """
    Download GWAS summary statistics from IEU OpenGWAS
    """
    print(f"\nAttempting to download GWAS ID: {gwas_id}")
    
    # Try different API endpoints
    endpoints = [
        f"https://api.gwas-api.mrcieu.ac.uk/associations/{gwas_id}",
        f"https://gwas-api.mrcieu.ac.uk/associations/{gwas_id}",
    ]
    
    for endpoint in endpoints:
        try:
            print(f"Trying: {endpoint}")
            response = requests.get(endpoint, timeout=60)
            if response.status_code == 200:
                # Try to parse as JSON first
                try:
                    data = response.json()
                    if isinstance(data, list) and len(data) > 0:
                        df = pd.DataFrame(data)
                        return df
                except:
                    # If not JSON, try as TSV/CSV
                    from io import StringIO
                    df = pd.read_csv(StringIO(response.text), sep='\t')
                    return df
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    return None


def format_gwas_data(df, output_file):
    """
    Format downloaded GWAS data to MR format
    """
    print(f"\nFormatting GWAS data...")
    print(f"Columns: {list(df.columns)}")
    print(f"Shape: {df.shape}")
    
    # Common column mappings for IEU OpenGWAS
    column_mappings = {
        'rsid': ['rsid', 'rsID', 'SNP', 'snp', 'variant', 'variant_id', 'id', 'ID'],
        'effect_allele': ['ea', 'EA', 'effect_allele', 'A1', 'a1', 'alt', 'ALT', 'effect_allele_freq'],
        'other_allele': ['oa', 'OA', 'other_allele', 'A2', 'a2', 'ref', 'REF', 'other_allele_freq'],
        'beta': ['beta', 'BETA', 'effect', 'Effect', 'b', 'B'],
        'se': ['se', 'SE', 'std_err', 'StdErr', 'standard_error'],
        'pval': ['p', 'P', 'pval', 'p_value', 'pvalue', 'P-value', 'PVAL'],
        'n': ['n', 'N', 'sample_size', 'SampleSize', 'N_SAMPLES', 'n_samples']
    }
    
    # Find matching columns
    mapped_cols = {}
    for target, possible_names in column_mappings.items():
        for col in df.columns:
            if col.lower() in [n.lower() for n in possible_names]:
                mapped_cols[target] = col
                break
    
    print(f"Mapped columns: {mapped_cols}")
    
    # Create output dataframe
    output_df = pd.DataFrame()
    for target, source_col in mapped_cols.items():
        if source_col in df.columns:
            output_df[target] = df[source_col]
    
    # Check required columns
    required = ['rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval']
    missing = [col for col in required if col not in output_df.columns]
    
    if missing:
        print(f"\nError: Missing required columns: {missing}")
        print("Available columns:", list(df.columns))
        return False
    
    # Handle sample size
    if 'n' not in output_df.columns:
        # Try to estimate or use a default
        if 'n_total' in df.columns:
            output_df['n'] = df['n_total']
        elif 'n_cases' in df.columns and 'n_controls' in df.columns:
            output_df['n'] = df['n_cases'] + df['n_controls']
        else:
            # Use a reasonable default or median if available
            output_df['n'] = 10000  # Placeholder
            print("Warning: Sample size not found, using placeholder value")
    
    # Clean data
    output_df = output_df.dropna(subset=required)
    
    # Ensure numeric types
    for col in ['beta', 'se', 'pval', 'n']:
        if col in output_df.columns:
            output_df[col] = pd.to_numeric(output_df[col], errors='coerce')
    
    # Remove invalid rows
    output_df = output_df[
        (output_df['se'] > 0) &
        (output_df['pval'].notna()) &
        (output_df['beta'].notna())
    ].copy()
    
    # Save
    output_df = output_df[['rsid', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'n']]
    output_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n✓ Formatted data saved to: {output_file}")
    print(f"  Rows: {output_df.shape[0]}")
    
    return True


def main():
    """
    Main function to search and download outcome data
    """
    print("=" * 70)
    print("Downloading GWAS Outcome Data from IEU OpenGWAS")
    print("=" * 70)
    
    # Common cardiovascular trait searches
    searches = [
        "left ventricular ejection fraction",
        "heart failure",
        "ejection fraction",
        "cardiac function"
    ]
    
    print("\nSearching for available datasets...")
    
    # Try to find a suitable dataset
    found_datasets = []
    
    for query in searches:
        print(f"\nSearching: '{query}'")
        results = search_ieu_opengwas(query)
        
        if results:
            if isinstance(results, list):
                found_datasets.extend(results[:5])  # Limit to first 5
            elif isinstance(results, dict) and 'data' in results:
                found_datasets.extend(results['data'][:5])
        
        time.sleep(1)  # Be polite to the API
    
    if not found_datasets:
        print("\n" + "=" * 70)
        print("Could not automatically find datasets via API.")
        print("=" * 70)
        print("\nManual download instructions:")
        print("\n1. Visit: https://gwas.mrcieu.ac.uk/")
        print("2. Search for one of these traits:")
        print("   - 'left ventricular ejection fraction'")
        print("   - 'heart failure'")
        print("   - 'ejection fraction'")
        print("3. Find a dataset and note its ID (format: ieu-a-XXXX)")
        print("4. Download the summary statistics")
        print("5. Use get_real_data.py to format it")
        print("\nAlternatively, you can use the TwoSampleMR R package:")
        print("  library(TwoSampleMR)")
        print("  ao <- available_outcomes()")
        print("  # Search for LVEF or heart failure")
        print("  dat <- extract_instruments(outcomes='ieu-a-XXXX')")
        print("  # Or download full summary stats")
        return
    
    print(f"\nFound {len(found_datasets)} potential datasets")
    
    # Try to download the first suitable one
    for i, dataset in enumerate(found_datasets[:3]):  # Try first 3
        if isinstance(dataset, dict):
            gwas_id = dataset.get('id', '')
            trait = dataset.get('trait', 'Unknown')
            
            print(f"\n{i+1}. Trying: {trait} (ID: {gwas_id})")
            
            if gwas_id:
                df = get_gwas_data(gwas_id, None)
                if df is not None and not df.empty:
                    output_file = "synthetic_outcome_lvef.txt"
                    if format_gwas_data(df, output_file):
                        print(f"\n✓ Successfully downloaded and formatted: {output_file}")
                        print(f"  Dataset: {trait}")
                        print(f"  ID: {gwas_id}")
                        return
    
    print("\n" + "=" * 70)
    print("Automatic download was not successful.")
    print("=" * 70)
    print("\nPlease download manually from:")
    print("  https://gwas.mrcieu.ac.uk/")
    print("\nThen use get_real_data.py to format it.")


if __name__ == "__main__":
    try:
        import requests
    except ImportError:
        print("Installing requests package...")
        import subprocess
        subprocess.check_call(["pip3", "install", "requests"])
        import requests
    
    main()

