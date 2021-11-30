library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggforce)

GCT <- read.delim(file="~/Desktop/GCGR_pub/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip=2)

GCT_GCGR <- filter(GCT, Description == "GCGR") %>% 
  as_tibble() %>%  
  t(.) %>% 
  as.data.frame() 

rown <- gsub(".", '-', rownames(GCT_GCGR), fixed = T)

rownames(GCT_GCGR) <- rown
GCT_GCGR <- t(GCT_GCGR) %>% as_tibble()

Tissue_table <- read_tsv("~/Desktop/GCGR_pub/study_design.txt")
Tissue_pheno <- read_tsv("~/Desktop/GCGR_pub/study_design-pheno.txt")

sudden_death <- filter(Tissue_pheno, DTHHRDY == 1 | DTHHRDY == 2 | DTHHRDY == 3) %>% dplyr::select(.,SUBJID)
annotation_usefull <- dplyr::select(Tissue_table, SAMPID, SMTS, SMTSD)

GCGR_sudden_death <- dplyr::select(GCT_GCGR, 1) 

for(i in sudden_death$SUBJID){
  GCGR_join <- dplyr::select(GCT_GCGR, starts_with(i)) 
  GCGR_sudden_death <- bind_cols(x= GCGR_sudden_death, y= GCGR_join)
}

GCGR_sudden_death <- dplyr::select(GCGR_sudden_death, -1) %>% as_tibble %>% t(.) %>% 
  as.data.frame() %>% rownames_to_column("SAMPID") %>% as_tibble() %>% dplyr::rename(., TPM = V1)

GCGR_sudden_death$TPM <- as.numeric(GCGR_sudden_death$TPM)

GCGR_sudden_death_metajoin <- inner_join(x = GCGR_sudden_death, y = annotation_usefull, by = "SAMPID")
distribution_tissues <- dplyr::count(GCGR_sudden_death_metajoin, SMTS)
distribution_tissues_fine <- dplyr::count(GCGR_sudden_death_metajoin, SMTSD)

GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Cervix Uteri")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTSD != "Bladder")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Breast")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Ovary")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Uterus")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Vagina")
GCGR_sudden_death_metajoin <- filter(GCGR_sudden_death_metajoin, SMTS != "Blood")

distribution_tissues

ggplot(GCGR_sudden_death_metajoin, aes(y=TPM, x=reorder(SMTS, -TPM))) + 
  geom_boxplot(aes(), show.legend = FALSE) + 
  geom_jitter(aes(color = SMTS), alpha = 0.3, width= 0.1, show.legend = FALSE) +
  xlab("") +
  theme_pubr() +
  theme(axis.text.x=element_text(angle=50,hjust=1))
