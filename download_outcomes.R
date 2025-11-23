# R script to download outcome GWAS data from IEU OpenGWAS
# Using TwoSampleMR package

# Install if needed
# install.packages("TwoSampleMR")

library(TwoSampleMR)

# Load exposure SNPs (from GTEx processed data)
exposure_snps <- read.table("synthetic_exposure_eqtl.txt", header=TRUE, sep="\t")$rsid

cat("Loaded", length(exposure_snps), "exposure SNPs\n\n")

# Dataset IDs to download
datasets <- c(
  "ieu-b-5120",
  "ieu-b-5121", 
  "finn-b-I9_HEARTFAIL_ALLCAUSE",
  "finn-b-I9_HEARTFAIL_AND_HYPERTCARDIOM",
  "prot-a-3066"
)

# Download outcome data for each dataset
for (id in datasets) {
  cat("=", rep("=", 50), "\n", sep="")
  cat("Downloading:", id, "\n")
  cat("=", rep("=", 50), "\n", sep="")
  
  tryCatch({
    # Extract outcome data for our exposure SNPs
    outcome_dat <- extract_outcome_data(snps=exposure_snps, outcomes=id)
    
    if (!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
      # Create output filename
      output_file <- paste0("outcome_", gsub("-", "_", id), ".txt")
      
      # Format for MR analysis (rsid, effect_allele, other_allele, beta, se, pval, n)
      formatted <- data.frame(
        rsid = outcome_dat$SNP,
        effect_allele = outcome_dat$effect_allele.outcome,
        other_allele = outcome_dat$other_allele.outcome,
        beta = outcome_dat$beta.outcome,
        se = outcome_dat$se.outcome,
        pval = outcome_dat$pval.outcome,
        n = outcome_dat$samplesize.outcome
      )
      
      # Write to file
      write.table(formatted, output_file, sep="\t", row.names=FALSE, quote=FALSE)
      cat("✓ Saved:", output_file, "\n")
      cat("  SNPs:", nrow(formatted), "\n")
      cat("  Sample size:", unique(formatted$n)[1], "\n\n")
    } else {
      cat("✗ No data retrieved for", id, "\n\n")
    }
  }, error = function(e) {
    cat("✗ Error downloading", id, ":", conditionMessage(e), "\n\n")
  })
  
  # Be polite to the API
  Sys.sleep(1)
}

cat("\n", rep("=", 52), "\n", sep="")
cat("Download complete!\n")
cat("Files saved with prefix: outcome_\n")
cat("You can now use these files with run.py\n")
cat(rep("=", 52), "\n", sep="")

