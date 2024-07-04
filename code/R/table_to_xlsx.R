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

filtered_maf <- subsetMaf(maffile, genes = genelist)
maf_data <- filtered_maf@data
maf_silent <- filtered_maf@maf.silent
maf_all <- rbind(maf_data, maf_silent, fill = T)
  
maf_all <- maf_all %>%
  mutate(VAF = t_alt_count / t_depth * 100)

write.xlsx(maf_all, outfile)
