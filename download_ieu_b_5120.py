#!/usr/bin/env python3
"""
Download and format ieu-b-5120 outcome dataset
"""

import pandas as pd
import requests
import os

DATASET_ID = "ieu-b-5120"
DATASET_URL = f"https://gwas.mrcieu.ac.uk/datasets/{DATASET_ID}/"

def get_dataset_info():
    """
    Get information about the dataset
    """
    print("=" * 70)
    print(f"Dataset: {DATASET_ID}")
    print("=" * 70)
    print(f"\nDataset page: {DATASET_URL}")
    
    try:
        # Try to get info from API
        url = f"https://gwas-api.mrcieu.ac.uk/gwasinfo/{DATASET_ID}"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            info = response.json()
            print(f"\nTrait: {info.get('trait', 'Unknown')}")
            print(f"Sample size: {info.get('sample_size', 'Unknown')}")
            print(f"Number of SNPs: {info.get('nsnp', 'Unknown')}")
            return info
    except Exception as e:
        print(f"\nCould not retrieve API info: {e}")
    
    print("\n" + "=" * 70)
    print("DOWNLOAD INSTRUCTIONS")
    print("=" * 70)
    print("\n1. Visit:", DATASET_URL)
    print("2. Click the 'Download' button to get summary statistics")
    print("3. Save the file (e.g., 'ieu_b_5120_raw.tsv' or similar)")
    print("4. Run this script again with the file path, OR")
    print("5. Use: python3 get_real_data.py to format it")
    
    return None


def format_ieu_b_5120(input_file, output_file="synthetic_outcome_lvef.txt"):
    """
    Format the downloaded ieu-b-5120 file for MR analysis
    """
    print(f"\n{'='*70}")
    print(f"Formatting {input_file} for MR analysis")
    print(f"{'='*70}")
    
    if not os.path.exists(input_file):
        print(f"\nError: File '{input_file}' not found!")
        print("\nPlease download the file first from:")
        print(f"  {DATASET_URL}")
        return False
    
    # Use the formatting function from get_real_data.py
    from get_real_data import format_data_for_mr
    
    success = format_data_for_mr(input_file, output_file, "outcome")
    
    if success:
        print(f"\n{'='*70}")
        print("SUCCESS!")
        print(f"{'='*70}")
        print(f"\nFormatted data saved to: {output_file}")
        print("This file will be used by run.py for MR analysis")
        print("\nYou can now run: python3 run.py")
    
    return success


if __name__ == "__main__":
    import sys
    
    # Get dataset info
    get_dataset_info()
    
    # If file path provided as argument, format it
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        format_ieu_b_5120(input_file)
    else:
        print("\n" + "=" * 70)
        print("USAGE")
        print("=" * 70)
        print("\nAfter downloading the file, run:")
        print(f"  python3 download_ieu_b_5120.py <downloaded_file_path>")
        print("\nOr use the interactive formatter:")
        print("  python3 get_real_data.py")

