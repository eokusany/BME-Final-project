#==== Installation of necessary packages - It is necessary to have an updated version of R and RTools =====

install.packages(c("devtools", "knitr", "rmarkdown", "ggplot2", "dplyr", "data.table", "ieugwasr", "vcfR", "coloc"))

# ====== Download needed packages, click 1 when asked ==========
library(devtools)
install_github(c("MRCIEU/TwoSampleMR","MRCIEU/MRInstruments"))
devtools::install_github("adletaw/captioner")
install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
devtools::install_github("rondolab/MR-PRESSO")
install_github("linjf15/MR_tricks")

# ============= Loading the necessary libraries ==============

library(MRInstruments)
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)
library(ggplot2)
library(data.table)
library(MRPRESSO)
library(vcfR)
library(coloc)
library(gridExtra)

# ============= Gene Lists ==============
positive_controls = c(RARG, CBR3, SPG7, MLH1)
gene_list = c(NOS3, PRDM2, WDR4, ZNF521, SP4, RIN3, ABCC1, ABCC9, ABCC5, HNMT, SLC22A17, ERCC2, MYH7, CYP2J2, COL1A2, SPG7, GPX3, ABCB4, GSTP1, PLCE1, GSTM1, CELF4, CYBA, HAS3, MLH1, POR, RAC2, CAT, ABCC2, ATP2B1, CBR1, ERBB2)
negative_controls = c(KRT1, KRT10, INS, AMY1)
housekeeping_genes = c(GAPDH, ACTB, RPLP0)

positive_variants = c(rs1854063, rs1332536, rs11140622)


# ============= Exposure Data ==============

gtex_heart <- data.table::fread("Heart_Left_Ventricle.v8.egenes.txt.gz")
gtex_data <- data.table::fread("Heart_Left_Ventricle.v8.signif_variant_gene_pairs.txt.gz")
ieu_data <- read.vcfR("ieu-b-5121.vcf.gz")
eqtl_sig <- data.table::fread("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
maf_eqtlgen <- data.table::fread("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")


filtered_exposure <- eqtl_sig %>% 
  filter(GeneSymbol %in% c("RARG", "CBR3", "SPG7", "MLH1", "CYBA", "POR", "SLC28A1", "SLC22A17", "SLC28A3", "ABCC2", "ABCC5", "HNMT", "GSTM1", "GAPDH", "ACTB", "RPLP0", "KRT1", "KRT10"))

filtered_exposure <- filtered_exposure %>%
  rename(
    effect_allele.exposure = AssessedAllele,
    other_allele.exposure = OtherAllele
  )


exposure <- merge(filtered_exposure, maf_eqtlgen, by="SNP")

filtered_exp <- data.frame(exposure)
filtered_exp$se <- (1/sqrt(2* filtered_exp$NrSamples * filtered_exp$AlleleB_all * (1 - filtered_exp$AlleleB_all)))
filtered_exp$beta <- (filtered_exp$Zscore * filtered_exp$se)


exp <- format_data(
  filtered_exp,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "AlleleB_all",
  z_col = "Zscore",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "PValue",
  gene_col = "GeneSymbol",
  samplesize_col = "NrSamples",
  chr_col = "SNPChr",
  pos_col = "SNPPos",
)


# ========== Instrument strength ================
exp$fstats <- (exp$beta.exposure^2)/(exp$se.exposure^2)

strong_exposure <- exp %>%
  filter(fstats > 50)


# =========== Clumping ============ -> Useless
# strong_exposure$rsid <- strong_exposure$SNP
#strong_exposure$pval <- strong_exposure$pval.exposure

#exposure_clumped <- ld_clump(strong_exposure, clump_kb = 10000, clump_r2 = 0.001, opengwas_jwt = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJkaWdhbmlAdWFsYmVydGEuY2EiLCJpYXQiOjE3NjMyNDgxNjAsImV4cCI6MTc2NDQ1Nzc2MH0.fVHt8-xwxQCqty4pPwMgxIragr-VC8JA7i4LKJNuaTf3CrtpEb7XMUkc7K3nucUyJkOq_WauhpI0Zj5MCPWCXVxGvMC6_iXIq7BRErycFGq-Dx0g3c_992M_rBFOL_UngCQeIERF-qf49L_tDnwH8Z1YhObnwYGu-gPCI1qtMnuhmNWqsXmtSEwBBOfr8ehw57R_tGEOWBeu5LbDxEzO4zX2AOMXgXIefuI5FK-y-0PnGR3XhxqF8qWqSpAjU4nv2fOi7H2u4oOruX4Q83occmLlKNgQso7iPVhhFLjnI-VssGJx7q8QxUveixBubHDwO61Da3ZQxAXyoPfCpLpJKg")

#strong_exposure <- strong_exposure %>% select(-fstats)


#============= Outcome data - Males ===============
LVEF_gwas_dat <- extract_outcome_data(snps = strong_exposure$SNP, outcomes = "ieu-b-5121", opengwas_jwt = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJkaWdhbmlAdWFsYmVydGEuY2EiLCJpYXQiOjE3NjMyNDgxNjAsImV4cCI6MTc2NDQ1Nzc2MH0.fVHt8-xwxQCqty4pPwMgxIragr-VC8JA7i4LKJNuaTf3CrtpEb7XMUkc7K3nucUyJkOq_WauhpI0Zj5MCPWCXVxGvMC6_iXIq7BRErycFGq-Dx0g3c_992M_rBFOL_UngCQeIERF-qf49L_tDnwH8Z1YhObnwYGu-gPCI1qtMnuhmNWqsXmtSEwBBOfr8ehw57R_tGEOWBeu5LbDxEzO4zX2AOMXgXIefuI5FK-y-0PnGR3XhxqF8qWqSpAjU4nv2fOi7H2u4oOruX4Q83occmLlKNgQso7iPVhhFLjnI-VssGJx7q8QxUveixBubHDwO61Da3ZQxAXyoPfCpLpJKg")


#============= Harmonize data ==============================
harmonised_data <- harmonise_data(strong_exposure, LVEF_gwas_dat, action = 2)

harmonised_data$samplesize.outcome = 17899  # Yang Sample Size


#============= Single SNP Filter ==============================
single_s <- mr_singlesnp(harmonised_data, all_method = "mr_egger_regression")
meo <- merge(single_s, strong_exposure, by="SNP")

mr_forest_plot(single_s)

meo2 <- meo %>%
  filter(p < 0.05)

meo3 <- merge(harmonised_data, meo2, by="SNP")
meo3 <- meo3 %>%
  rename(
    effect_allele.exposure = effect_allele.exposure.x,
    other_allele.exposure = other_allele.exposure.x,
    beta.exposure = beta.exposure.x,
    eaf.exposure = eaf.exposure.x,
    id.outcome = id.outcome.x,
    outcome = outcome.x,
    chr.exposure = chr.exposure.x,
    pos.exposure = pos.exposure.x,
    z.exposure = z.exposure.x,
    gene.exposure = gene.exposure.x,
    pval.exposure.x = pval.exposure.x,
    samplesize.exposure = samplesize.exposure.x,
    pval_origin.exposure = pval_origin.exposure.x,
    se.exposure = se.exposure.x,
    fstats = fstats.x,
  )

single <- mr_singlesnp(meo3, all_method = "mr_egger_regression")
mr_forest_plot(single)

meo4 <- merge(single, strong_exposure, by="SNP")
steiger_data <- directionality_test(harmonised_data)

head(steiger_data)

# ============== Variants of interest ===========
CBR3 = c(rs2835287, rs62229301, rs8129035)
SPG7 = c(rs67689854, rs9930567, rs55637757)
MLH1 = c(rs925649, rs7624304, rs6550471)

CYBA = c(rs12933512, rs12933505, rs4468599)
POR = c(rs6971988, rs12535057, rs60791363)

SLC22A17 = c(rs72686239)
SLC28A3 = c(rs4877275, rs11140622, rs1332536)

ABCC2 = c(rs11190259, rs12218468, rs11190170)
ABCC5 = c(rs113927964, rs112975088, rs111476402)
HNMT = c(rs72854949, rs10496794, rs72854953)
GSTM1 = c(rs10127988, rs1799875, rs11587687)

GAPDH = c(rs116979198, rs77843719, rs117684743)

KRT1 = c(rs11608629, rs11170178, rs676255)
KRT10 = c(rs7214830, rs76843554)


meo3_1 <- meo3 %>%
  filter(
    SNP %in% c("rs2835287", "rs62229301", "rs8129035", "rs67689854", "rs9930567", "rs55637757", "rs925649", "rs7624304", "rs6550471", 
               "rs12933512", "rs12933505", "rs4468599", "rs6971988", "rs12535057", "rs60791363",
               "rs72686239", "rs4877275", "rs11140622", "rs1332536",
               "rs11190259", "rs12218468", "rs11190170", "rs113927964", "rs112975088", "rs111476402", "rs72854949", "rs10496794", "rs72854953", "rs10127988", "rs1799875", "rs11587687",
               "rs116979198", "rs77843719", "rs117684743",
               "rs11608629", "rs11170178", "rs676255", "rs7214830", "rs76843554")
  )

mr_report(harmonised_data,
          output_path = ".",
          output_type = "html",
          author = "Pulkit",
          study = "study",
          path = system.file("reports", package = "TwoSampleMR"))


#============= Next Steps ==============================


# Steiger Filtering
## Test if variance explained in exposure > outcome
steiger_data <- directionality_test(meo3)



####################################

# Complete analysis function
run_complete_mr_analysis <- function(harmonised_data, gene_name) {
  
  tmp_data <- harmonised_data %>%
    filter(
      gene.exposure == gene_name
    )
  
  # 1. Main MR methods
  mr_results <- mr(tmp_data, method_list = c("mr_ivw", "mr_egger_regression", 
                                                    "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
  
  # 2. Sensitivity analyses
  heterogeneity <- mr_heterogeneity(tmp_data)
  pleiotropy <- mr_pleiotropy_test(tmp_data)
  loo <- mr_leaveoneout(tmp_data)
  single_snp <- mr_singlesnp(tmp_data)
  
  # 3. Generate plots
  scatter_plot <- mr_scatter_plot(mr_results, tmp_data)
  forest_plot <- mr_forest_plot(single_snp)
  funnel_plot <- mr_funnel_plot(single_snp)
  
  # 4. Compile results
  results <- list(
    mr_results = mr_results,
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy,
    leave_one_out = loo,
    single_snp = single_snp,
    plots = list(scatter = scatter_plot, forest = forest_plot, funnel = funnel_plot)
  )
  
  return(results)
}

########## Run for a specific gene ###########
CBR3_results <- run_complete_mr_analysis(harmonised_data, "CBR3")
SPG7_results <- run_complete_mr_analysis(harmonised_data, "SPG7")
MLH1_results <- run_complete_mr_analysis(harmonised_data, "MLH1")

CYBA_results <- run_complete_mr_analysis(harmonised_data, "CYBA")
POR_results <- run_complete_mr_analysis(harmonised_data, "POR")

SLC22A17_results <- run_complete_mr_analysis(harmonised_data, "SLC22A17")
SLC28A3_results <- run_complete_mr_analysis(harmonised_data, "SLC28A3")

ABCC2_results <- run_complete_mr_analysis(harmonised_data, "ABCC2")
ABCC5_results <- run_complete_mr_analysis(harmonised_data, "ABCC5")
HNMT_results <- run_complete_mr_analysis(harmonised_data, "HNMT")
GSTM1_results <- run_complete_mr_analysis(harmonised_data, "GSTM1")

GAPDH_results <- run_complete_mr_analysis(harmonised_data, "GAPDH")

KRT1_results <- run_complete_mr_analysis(harmonised_data, "KRT1")
KRT10_results <- run_complete_mr_analysis(harmonised_data, "KRT10")

RARG_results <- run_complete_mr_analysis(harmonised_data, "RARG")
RPLP0_results <- run_complete_mr_analysis(harmonised_data, "RPLP0")


########## Run for a specific gene ###########
sp1 <- CBR3_results$plots$scatter[[1]]
  
sp1 <- sp1 + labs(caption = "CBR3") 
sp1 <- sp1 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp2 <- SPG7_results$plots$scatter[[1]]

sp2 <- sp2 + labs(caption = "SPG7") 
sp2 <- sp2 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp3 <- MLH1_results$plots$scatter[[1]]

sp3 <- sp3 + labs(caption = "MLH1") 
sp3 <- sp3 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp4 <- CYBA_results$plots$scatter[[1]]

sp4 <- sp4 + labs(caption = "CYBA") 
sp4 <- sp4 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp5 <- POR_results$plots$scatter[[1]]

sp5 <- sp5 + labs(caption = "POR") 
sp5 <- sp5 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp6 <- SLC22A17_results$plots$scatter[[1]]

sp6 <- sp6 + labs(caption = "SLC22A17") 
sp6 <- sp6 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp7 <- SLC28A3_results$plots$scatter[[1]]

sp7 <- sp7 + labs(caption = "SLC28A3") 
sp7 <- sp7 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp8 <- ABCC2_results$plots$scatter[[1]]

sp8 <- sp8 + labs(caption = "ABCC2") 
sp8 <- sp8 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp9 <- ABCC5_results$plots$scatter[[1]]

sp9 <- sp9 + labs(caption = "ABCC5") 
sp9 <- sp9 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp10 <- HNMT_results$plots$scatter[[1]]

sp10 <- sp10 + labs(caption = "HNMT") 
sp10 <- sp10 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp11 <- GSTM1_results$plots$scatter[[1]]

sp11 <- sp11 + labs(caption = "GSTM1") 
sp11 <- sp11 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp12 <- GAPDH_results$plots$scatter[[1]]

sp12 <- sp12 + labs(caption = "GAPDH") 
sp12 <- sp12 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp13 <- KRT1_results$plots$scatter[[1]]

sp13 <- sp13 + labs(caption = "KRT1") 
sp13 <- sp13 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp14 <- KRT10_results$plots$scatter[[1]]

sp14 <- sp14 + labs(caption = "KRT10") 
sp14 <- sp14 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp15 <- RARG_results$plots$scatter[[1]]

sp15 <- sp15 + labs(caption = "RARG") 
sp15 <- sp15 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp16 <- RPLP0_results$plots$scatter[[1]]

sp16 <- sp16 + labs(caption = "RPLP0") 
sp16 <- sp16 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


grid.arrange(sp1, sp2, sp3, sp15)

grid.arrange(sp4, sp5)

grid.arrange(sp6, sp7)

grid.arrange(sp8, sp9, sp10, sp11)

grid.arrange(sp12, sp13, sp14, sp16)

######## Funnel Plots ###########
sp1 <- CBR3_results$plots$forest[[1]]

sp1 <- sp1 + labs(caption = "CBR3") 
sp1 <- sp1 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp2 <- SPG7_results$plots$forest[[1]]

sp2 <- sp2 + labs(caption = "SPG7") 
sp2 <- sp2 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp3 <- MLH1_results$plots$forest[[1]]

sp3 <- sp3 + labs(caption = "MLH1") 
sp3 <- sp3 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp4 <- CYBA_results$plots$forest[[1]]

sp4 <- sp4 + labs(caption = "CYBA") 
sp4 <- sp4 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp5 <- POR_results$plots$forest[[1]]

sp5 <- sp5 + labs(caption = "POR") 
sp5 <- sp5 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp6 <- SLC22A17_results$plots$forest[[1]]

sp6 <- sp6 + labs(caption = "SLC22A17") 
sp6 <- sp6 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp7 <- SLC28A3_results$plots$forest[[1]]

sp7 <- sp7 + labs(caption = "SLC28A3") 
sp7 <- sp7 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp8 <- ABCC2_results$plots$forest[[1]]

sp8 <- sp8 + labs(caption = "ABCC2") 
sp8 <- sp8 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp9 <- ABCC5_results$plots$forest[[1]]

sp9 <- sp9 + labs(caption = "ABCC5") 
sp9 <- sp9 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp10 <- HNMT_results$plots$forest[[1]]

sp10 <- sp10 + labs(caption = "HNMT") 
sp10 <- sp10 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp11 <- GSTM1_results$plots$forest[[1]]

sp11 <- sp11 + labs(caption = "GSTM1") 
sp11 <- sp11 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp12 <- GAPDH_results$plots$forest[[1]]

sp12 <- sp12 + labs(caption = "GAPDH") 
sp12 <- sp12 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp13 <- KRT1_results$plots$forest[[1]]

sp13 <- sp13 + labs(caption = "KRT1") 
sp13 <- sp13 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp14 <- KRT10_results$plots$forest[[1]]

sp14 <- sp14 + labs(caption = "KRT10") 
sp14 <- sp14 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp15 <- RARG_results$plots$forest[[1]]

sp15 <- sp15 + labs(caption = "RARG") 
sp15 <- sp15 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


sp16 <- RPLP0_results$plots$forest[[1]]

sp16 <- sp16 + labs(caption = "RPLP0") 
sp16 <- sp16 + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2), face = "bold"))


grid.arrange(sp1, sp2, sp3, sp15)

grid.arrange(sp4, sp5)

grid.arrange(sp6, sp7)

grid.arrange(sp8, sp9, sp10, sp11)

grid.arrange(sp12, sp13, sp14, sp16)

########### Tables ###########
CBR3_results$mr_results
MLH1_results$mr_results
SPG7_results$mr_results
RARG_results$mr_results

CYBA_results$mr_results
POR_results$mr_results

SLC22A17_results$mr_results
SLC28A3_results$mr_results

ABCC2_results$mr_results
ABCC5_results$mr_results
HNMT_results$mr_results
GSTM1_results$mr_results

GAPDH_results$mr_results
KRT1_results$mr_results
KRT10_results$mr_results
RPLP0_results$mr_results



