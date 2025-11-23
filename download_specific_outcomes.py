#!/usr/bin/env python3
"""
Download specific GWAS outcome datasets from IEU OpenGWAS
Uses provided dataset IDs to download and format outcome data
"""

import pandas as pd
import requests
import json
import time
import os

# Dataset IDs provided by user
DATASET_IDS = [
    "ieu-b-5120",
    "ieu-b-5121", 
    "finn-b-I9_HEARTFAIL_ALLCAUSE",
    "finn-b-I9_HEARTFAIL_AND_HYPERTCARDIOM",
    "prot-a-3066"
]

def get_dataset_info(dataset_id):
    """
    Get information about a dataset from IEU OpenGWAS
    """
    try:
        # Try to get dataset info
        url = f"https://gwas-api.mrcieu.ac.uk/gwasinfo/{dataset_id}"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
    except Exception as e:
        print(f"  Error getting info: {e}")
    return None


def download_gwas_data(dataset_id, output_file):
    """
    Download GWAS summary statistics for specific SNPs
    Note: Full downloads may require API access or manual download
    """
    print(f"\n{'='*70}")
    print(f"Dataset: {dataset_id}")
    print(f"{'='*70}")
    
    # Get dataset info first
    info = get_dataset_info(dataset_id)
    if info:
        print(f"Trait: {info.get('trait', 'Unknown')}")
        print(f"Sample size: {info.get('sample_size', 'Unknown')}")
        print(f"Number of SNPs: {info.get('nsnp', 'Unknown')}")
    
    # Try to download via API
    # Method 1: Try associations endpoint (requires SNP list)
    # Method 2: Direct download link
    print(f"\nDownload options:")
    print(f"1. Manual download: https://gwas.mrcieu.ac.uk/datasets/{dataset_id}/")
    print(f"2. API endpoint: https://gwas-api.mrcieu.ac.uk/associations/{dataset_id}")
    
    # For now, provide instructions for manual download
    # Full API access may require authentication or specific endpoints
    return None


def format_downloaded_data(input_file, dataset_id, output_file):
    """
    Format downloaded GWAS data to MR format
    """
    print(f"\nFormatting data from {input_file}...")
    
    try:
        # Try different separators
        for sep in ['\t', ',', ' ']:
            try:
                df = pd.read_csv(input_file, sep=sep, low_memory=False)
                if df.shape[1] > 1:
                    break
            except:
                continue
        
        print(f"  Columns: {list(df.columns)}")
        print(f"  Shape: {df.shape}")
        
        # Use the formatting function from get_real_data.py
        from get_real_data import format_data_for_mr
        return format_data_for_mr(input_file, output_file, "outcome")
        
    except Exception as e:
        print(f"Error formatting: {e}")
        return False


def main():
    """
    Main function to download and process outcome datasets
    """
    print("=" * 70)
    print("Downloading Specific IEU OpenGWAS Outcome Datasets")
    print("=" * 70)
    
    print(f"\nDatasets to process: {len(DATASET_IDS)}")
    for i, dataset_id in enumerate(DATASET_IDS, 1):
        print(f"  {i}. {dataset_id}")
    
    print("\n" + "=" * 70)
    print("INSTRUCTIONS FOR MANUAL DOWNLOAD")
    print("=" * 70)
    
    print("\nFor each dataset, you can:")
    print("\n1. Visit the dataset page directly:")
    for dataset_id in DATASET_IDS:
        print(f"   https://gwas.mrcieu.ac.uk/datasets/{dataset_id}/")
    
    print("\n2. Click 'Download' to get summary statistics")
    print("3. Save the file with a descriptive name")
    print("4. Use format_downloaded_data() or get_real_data.py to format it")
    
    print("\n" + "=" * 70)
    print("USING R TwoSampleMR (Recommended)")
    print("=" * 70)
    
    print("\nIf you have R installed, this is the easiest method:")
    print("""
library(TwoSampleMR)

# Load your exposure SNPs
exposure_snps <- read.table("synthetic_exposure_eqtl.txt", header=TRUE, sep="\\t")$rsid

# Download outcome data for each dataset
datasets <- c("ieu-b-5120", "ieu-b-5121", "finn-b-I9_HEARTFAIL_ALLCAUSE", 
              "finn-b-I9_HEARTFAIL_AND_HYPERTCARDIOM", "prot-a-3066")

for (id in datasets) {
  cat("Downloading:", id, "\\n")
  outcome_dat <- extract_outcome_data(snps=exposure_snps, outcomes=id)
  
  if (!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
    output_file <- paste0("outcome_", gsub("-", "_", id), ".txt")
    write.table(outcome_dat, output_file, sep="\\t", row.names=FALSE, quote=FALSE)
    cat("Saved:", output_file, "\\n")
  }
}
""")
    
    print("\n" + "=" * 70)
    print("AUTOMATIC DOWNLOAD ATTEMPT")
    print("=" * 70)
    
    # Try to get info for each dataset
    for dataset_id in DATASET_IDS:
        info = get_dataset_info(dataset_id)
        if info:
            print(f"\n✓ {dataset_id}: {info.get('trait', 'Unknown')}")
            print(f"  Sample size: {info.get('sample_size', 'N/A')}")
        else:
            print(f"\n✗ {dataset_id}: Could not retrieve info")
        time.sleep(0.5)  # Be polite to API
    
    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print("\n1. Download data manually from the links above, OR")
    print("2. Use R TwoSampleMR script above, OR")
    print("3. Contact me with the downloaded file paths to format them")
    print("\nOnce you have the data files, run:")
    print("  python3 get_real_data.py")
    print("  OR")
    print("  from get_real_data import format_data_for_mr")
    print("  format_data_for_mr('your_file.tsv', 'synthetic_outcome_lvef.txt', 'outcome')")


if __name__ == "__main__":
    try:
        import requests
    except ImportError:
        print("Installing requests...")
        import subprocess
        subprocess.check_call(["pip3", "install", "requests", "--quiet"])
        import requests
    
    main()

