library("tidyverse")
library("writexl")
library("vroom")

files <- fs::dir_ls("/disk/work/diagnostics/QAML/gb1/annotated_vcf/run_29/vep", recurse = TRUE, glob = "*.tab")

tables <- vroom(files, id = "path", col_names = FALSE, skip = 1) 


dataset1 <- tables %>%
  mutate(Allele = str_split_fixed(X102, "\\|", 54)[, 1]) %>%
  mutate(Consequence = str_split_fixed(X102, "\\|", 54)[, 2]) %>%
  mutate(IMPACT = str_split_fixed(X102, "\\|", 54)[, 3]) %>%
  mutate(SYMBOL = str_split_fixed(X102, "\\|", 54)[, 4]) %>%
  mutate(Gene = str_split_fixed(X102, "\\|", 54)[, 5]) %>%
  mutate(Feature_type = str_split_fixed(X102, "\\|", 54)[, 6]) %>%
  mutate(Feature = str_split_fixed(X102, "\\|", 54)[, 7]) %>%
  mutate(BIOTYPE = str_split_fixed(X102, "\\|", 54)[, 8]) %>%
  mutate(EXON = str_split_fixed(X102, "\\|", 54)[, 9]) %>%
  mutate(INTRON = str_split_fixed(X102, "\\|", 54)[, 10]) %>%
  mutate(HGVSc = str_split_fixed(X102, "\\|", 54)[, 11]) %>%
  mutate(HGVSp = str_split_fixed(X102, "\\|", 54)[, 12]) %>%
  mutate(cDNA_position = str_split_fixed(X102, "\\|", 54)[, 13]) %>%
  mutate(CDS_position = str_split_fixed(X102, "\\|", 54)[, 14]) %>%
  mutate(Protein_position = str_split_fixed(X102, "\\|", 54)[, 15]) %>%
  mutate(Amino_acids= str_split_fixed(X102, "\\|", 54)[, 16]) %>%
  mutate(Codons = str_split_fixed(X102, "\\|", 54)[, 17]) %>%
  mutate(Existing_variation = str_split_fixed(X102, "\\|", 54)[, 18]) %>%
  mutate(REF_ALLELE = str_split_fixed(X102, "\\|", 54)[, 19]) %>%
  mutate(DISTANCE = str_split_fixed(X102, "\\|", 54)[, 20]) %>%
  mutate(STRAND = str_split_fixed(X102, "\\|", 54)[, 21]) %>%
  mutate(FLAGS = str_split_fixed(X102, "\\|", 54)[, 22]) %>%
  mutate(SYMBOL_SOURCE = str_split_fixed(X102, "\\|", 54)[, 23]) %>%
  mutate(HGNC_ID = str_split_fixed(X102, "\\|", 54)[, 24]) %>%
  mutate(CANONICAL = str_split_fixed(X102, "\\|", 54)[, 25]) %>%
  mutate(CCDS = str_split_fixed(X102, "\\|", 54)[, 26]) %>%
  mutate(ENSP = str_split_fixed(X102, "\\|", 54)[, 27]) %>%
  mutate(REFSEQ_MATCH = str_split_fixed(X102, "\\|", 54)[, 28]) %>%
  mutate(SOURCE = str_split_fixed(X102, "\\|", 54)[, 29]) %>%
  mutate(REFSEQ_OFFSET = str_split_fixed(X102, "\\|", 54)[, 30]) %>%
  mutate(GIVEN_REF = str_split_fixed(X102, "\\|", 54)[, 31]) %>%
  mutate(USED_REF = str_split_fixed(X102, "\\|", 54)[, 32]) %>%
  mutate(BAM_EDIT = str_split_fixed(X102, "\\|", 54)[, 33]) %>%
  mutate(HGVS_OFFSET = str_split_fixed(X102, "\\|", 54)[, 34]) %>%
  mutate(AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 35])) %>%
  mutate(AFR_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 36])) %>%
  mutate(AMR_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 37])) %>%
  mutate(EAS_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 38])) %>%
  mutate(EUR_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 39])) %>%
  mutate(SAS_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 40])) %>%
  mutate(AA_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 41])) %>%
  mutate(EA_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 42])) %>%
  mutate(gnomAD_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 43])) %>%
  mutate(gnomAD_AFR_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 44])) %>%
  mutate(gnomAD_AMR_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 45])) %>%
  mutate(gnomAD_ASJ_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 46])) %>%
  mutate(gnomAD_EAS_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 47])) %>%
  mutate(gnomAD_FIN_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 48])) %>%
  mutate(gnomAD_NFE_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 49])) %>%
  mutate(gnomAD_OTH_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 50])) %>%
  mutate(gnomAD_SAS_AF = as.numeric(str_split_fixed(X102, "\\|", 54)[, 51])) %>%
  mutate(CLIN_SIG = str_split_fixed(X102, "\\|", 54)[, 52]) %>%
  mutate(SOMATIC = str_split_fixed(X102, "\\|", 54)[, 53]) %>%
  mutate(PHENO = str_split_fixed(X102, "\\|", 55)[, 54]) %>%
  mutate(PUBMED = str_split_fixed(X102, "\\|", 55)[, 55])

VEP_data <- dataset1 

Final_dataset <- VEP_data %>%
  mutate(Sample_ID_pre_3 = str_split_fixed(path, "\\/", 10)[, 10]) %>%
  mutate(Sample_ID_pre_2 = str_split_fixed(Sample_ID_pre_3, "\\.", 3)[, 2]) %>%
  mutate(Sample_ID = str_split_fixed(Sample_ID_pre_2, "\\_", 2)[, 1]) %>%
  select(-Sample_ID_pre_3, -Sample_ID_pre_2, -path, -X16) 

write_xlsx(Final_dataset, "/disk/work/diagnostics/QALL/gb1/annotation_final_results/QAML_run_29_full_annotation-2022-11-28.xlsx")
write_xlsx(dataset2, "/disk/work/diagnostics/QAML/gb1/annotated_vcf/run_27/QAML_run_27_HEADER-2022-10-07.xlsx")

dataset2 <- tables[1,]

