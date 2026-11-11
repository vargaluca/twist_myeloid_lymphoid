library(openxlsx)
library(dbplyr)
library(tidyverse)
library(maftools)

args <- commandArgs(trailingOnly = T)
file <- args[1]
outfile <- args[2]

maffile <- read.maf(file)
gene_csv <- read_csv(
  "/disk/work/diagnostics/TWLM/lv1/twist_myeloid_lymphoid/bed_files/updated_2024/Twist_myeloid+lymphoid_genlista_update_2407.csv",
  col_names = T)
genelist <- gene_csv$gene

# covered_regions <- read_tsv("/disk/work/diagnostics/DLBC/lv1/DLBCL_genclass_K21/bed_files/S3405393/DLBCL_final_comprehensive_panel_1_Targets.txt", col_names = c("target"))
# genelist_table <- covered_regions %>% 
#   filter(grepl("#", target) == FALSE) %>% 
#   filter(grepl("chr", target) == FALSE)
# genelist <- genelist_table$target

filtered_maf <- subsetMaf(maffile, genes = genelist)
maf_data <- filtered_maf@data
maf_silent <- filtered_maf@maf.silent
maf_all <- rbind(maf_data, maf_silent, fill = T)
  
maf_all <- maf_all %>%
  mutate(VAF = t_alt_count / t_depth * 100)

write.xlsx(maf_all, outfile)
