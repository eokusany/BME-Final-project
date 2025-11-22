#!/usr/bin/env python3
"""
Data Download Script for AIC Mendelian Randomization Analysis
Downloads eQTL (exposure) and GWAS (outcome) data from public databases
"""

import pandas as pd
import os
from pathlib import Path

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("Note: 'requests' module not installed. Some download features may be limited.")

def download_gtex_eqtl(gene_name=None, tissue="Heart_Left_Ventricle", output_file="exposure_eqtl.txt"):
    """
    Download eQTL data from GTEx Portal
    Note: GTEx requires manual download or API access
    For now, this provides instructions and template
    """
    print("=" * 70)
    print("GTEx eQTL Data Download")
    print("=" * 70)
    print("\nGTEx Portal: https://www.gtexportal.org/home/eqtlDashboardPage")
    print("\nInstructions:")
    print("1. Go to GTEx Portal eQTL Dashboard")
    print("2. Select tissue: Heart - Left Ventricle or Heart - Atrial Appendage")
    print("3. Enter gene symbol or browse by gene")
    print("4. Download cis-eQTL results (all significant associations)")
    print("5. Format should be: rsid, effect_allele, other_allele, beta, se, pval, n")
    print("\nAlternatively, use GTEx API:")
    print("  API endpoint: https://gtexportal.org/api/v2/association/singleTissueEqtl")
    print("\n" + "=" * 70)
    
    # Create a template file with instructions
    template = """# GTEx eQTL Data Template
# Replace this file with actual GTEx data
# Required columns: rsid, effect_allele, other_allele, beta, se, pval, n
# 
# Example format:
# rsid	effect_allele	other_allele	beta	se	pval	n
# rs12345	A	G	0.15	0.02	1e-10	500
"""
    
    if not os.path.exists(output_file):
        with open(output_file, 'w') as f:
            f.write(template)
        print(f"\nTemplate file created: {output_file}")
        print("Please replace with actual GTEx data")
    
    return output_file


def download_ieu_opengwas(trait_id=None, output_file="outcome_lvef.txt"):
    """
    Download GWAS data from IEU OpenGWAS
    Uses the ieugwasr R package via rpy2 or direct API calls
    """
    print("\n" + "=" * 70)
    print("IEU OpenGWAS Data Download")
    print("=" * 70)
    
    # Common AIC-related GWAS IDs from IEU OpenGWAS
    # You may need to search for specific IDs
    potential_ids = {
        "LVEF": "ieu-a-XXXX",  # Replace with actual ID
        "Heart Failure": "ieu-a-XXXX",
        "Cardiac Troponin": "ieu-a-XXXX",
        "QT Interval": "ieu-a-XXXX"
    }
    
    print("\nIEU OpenGWAS: https://gwas.mrcieu.ac.uk/")
    print("\nTo download data:")
    print("1. Search for your trait (e.g., 'left ventricular ejection fraction', 'anthracycline')")
    print("2. Note the GWAS ID (format: ieu-a-XXXX or ebi-a-XXXX)")
    print("3. Use the API or R package 'ieugwasr' to download")
    print("\nPython method (requires requests):")
    print("  API: https://api.gwas-api.mrcieu.ac.uk/")
    print("\nR method (recommended):")
    print("  install.packages('ieugwasr')")
    print("  library(ieugwasr)")
    print("  dat <- associations(variants=variants, id='ieu-a-XXXX')")
    print("\n" + "=" * 70)
    
    # Try to use requests if we have a specific ID
    if trait_id:
        try:
            api_url = f"https://api.gwas-api.mrcieu.ac.uk/associations/{trait_id}"
            print(f"\nAttempting to download from: {api_url}")
            # Note: This may require authentication or specific format
            print("Note: Direct API access may require authentication")
        except Exception as e:
            print(f"Error: {e}")
    
    # Create template
    template = """# IEU OpenGWAS Outcome Data Template
# Replace this file with actual GWAS summary statistics
# Required columns: rsid, effect_allele, other_allele, beta, se, pval, n
#
# Example format:
# rsid	effect_allele	other_allele	beta	se	pval	n
# rs12345	A	G	0.10	0.015	1e-8	1000
"""
    
    if not os.path.exists(output_file):
        with open(output_file, 'w') as f:
            f.write(template)
        print(f"\nTemplate file created: {output_file}")
        print("Please replace with actual IEU OpenGWAS data")
    
    return output_file


def download_eqtlgen(output_file="exposure_eqtlgen.txt"):
    """
    Download eQTL data from eQTLGen
    """
    print("\n" + "=" * 70)
    print("eQTLGen Data Download")
    print("=" * 70)
    print("\neQTLGen: https://www.eqtlgen.org/")
    print("\nInstructions:")
    print("1. Visit https://www.eqtlgen.org/cis-eqtls.html")
    print("2. Download cis-eQTL results (all significant)")
    print("3. Format the data to match required columns")
    print("4. Required columns: rsid, effect_allele, other_allele, beta, se, pval, n")
    print("\n" + "=" * 70)
    
    template = """# eQTLGen Data Template
# Replace this file with actual eQTLGen data
# Required columns: rsid, effect_allele, other_allele, beta, se, pval, n
"""
    
    if not os.path.exists(output_file):
        with open(output_file, 'w') as f:
            f.write(template)
        print(f"\nTemplate file created: {output_file}")
        print("Please replace with actual eQTLGen data")
    
    return output_file


def download_gwas_catalog(output_file="outcome_gwas_catalog.txt"):
    """
    Download GWAS data from GWAS Catalog
    """
    print("\n" + "=" * 70)
    print("GWAS Catalog Data Download")
    print("=" * 70)
    print("\nGWAS Catalog: https://www.ebi.ac.uk/gwas/")
    print("\nInstructions:")
    print("1. Search for 'anthracycline cardiotoxicity' or related terms")
    print("2. Download summary statistics if available")
    print("3. Format to match required columns")
    print("4. Required columns: rsid, effect_allele, other_allele, beta, se, pval, n")
    print("\n" + "=" * 70)
    
    template = """# GWAS Catalog Data Template
# Replace this file with actual GWAS Catalog data
# Required columns: rsid, effect_allele, other_allele, beta, se, pval, n
"""
    
    if not os.path.exists(output_file):
        with open(output_file, 'w') as f:
            f.write(template)
        print(f"\nTemplate file created: {output_file}")
        print("Please replace with actual GWAS Catalog data")
    
    return output_file


def main():
    """
    Main function to guide data download
    """
    print("=" * 70)
    print("AIC Mendelian Randomization - Data Download Guide")
    print("=" * 70)
    print("\nThis script provides instructions for downloading real data.")
    print("You need to download data from the following sources:\n")
    
    print("\n1. EXPOSURE DATA (eQTL - Gene Expression):")
    print("   - GTEx: Heart Left Ventricle or Heart Atrial Appendage")
    print("   - eQTLGen: Blood eQTLs")
    print("   - Format: rsid, effect_allele, other_allele, beta, se, pval, n")
    
    print("\n2. OUTCOME DATA (GWAS - AIC/Cardiac Traits):")
    print("   - IEU OpenGWAS: LVEF, Heart Failure, Cardiac Troponin, QT Interval")
    print("   - GWAS Catalog: Anthracycline cardiotoxicity")
    print("   - Format: rsid, effect_allele, other_allele, beta, se, pval, n")
    
    print("\n" + "=" * 70)
    
    # Create template files
    exposure_file = download_gtex_eqtl()
    outcome_file = download_ieu_opengwas()
    
    print("\n" + "=" * 70)
    print("Next Steps:")
    print("=" * 70)
    print("\n1. Download actual data from the sources above")
    print(f"2. Replace '{exposure_file}' with your eQTL data")
    print(f"3. Replace '{outcome_file}' with your GWAS data")
    print("4. Ensure data format matches: rsid, effect_allele, other_allele, beta, se, pval, n")
    print("5. Update run.py to use the new file names")
    print("6. Run: python3 run.py")
    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()

