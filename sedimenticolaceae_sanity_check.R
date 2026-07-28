library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(microViz)
library(phyloseq)
library(tibble)



phylo_elevation <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07_16S_raw_phyloseq.rds") %>%
  tax_filter(min_prevalence = 10,
             prev_detection_threshold = 5,
             min_total_abundance = 5,
             min_sample_abundance = 5) %>%
  subset_taxa(!(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))

phylo_elevation_prop <- transform_sample_counts(phylo_elevation, function(otu) otu/sum(otu))
sulfox <- subset_taxa(phylo_elevation_prop, Genus %in% c("Sedimenticola", "Candidatus Thiodiazotropha"))
sulfox@otu_table


sulfox@otu_table %>% data.frame() %>% writexl::write_xlsx("sulfox_relative_abundance.xlsx")


