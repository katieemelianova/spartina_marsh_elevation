
library(patchwork)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)


samplesheet <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07.xlsx") %>% dplyr::select(`JMF sample ID`, `User sample ID`)
sample_annotation <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/Sequencing/March2025_LeFaouAmplicon/LeFaou_samplesheet.xlsx") %>% dplyr::select(`User sample ID`, `Sample description`)
sample_mapping <- left_join(samplesheet, sample_annotation, by="User sample ID") %>% rename("JMF_sample_ID" = "JMF sample ID")


phylo_elevation <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07_16S_raw_phyloseq.rds") %>%
                   tax_filter(min_prevalence = 10,
                   prev_detection_threshold = 5,
                   min_total_abundance = 5,
                   min_sample_abundance = 5) %>%
  ps_join(sample_mapping, by = "JMF_sample_ID") %>%
  subset_taxa(!(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))



phylo_elevation %<>% subset_samples(!(is.na(Sample.description)) & !(Sample.description %in% c("Sedment unknown", "Root unknown")))


phylo_elevation_prop <- transform_sample_counts(phylo_elevation, function(otu) otu/sum(otu))
ord.nmds.bray_elevation <- ordinate(phylo_elevation_prop, method="NMDS", distance="bray")

plot_ordination(phylo_elevation_prop, ord.nmds.bray_elevation, color="Sample.description", title="Bray NMDS") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  ggtitle("")











# geometric means function
# doing it this way because see here:
# https://github.com/joey711/phyloseq/issues/445
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

run_deseq <- function(phylo_object, design_term){
  formula_parsed<-paste("~", design_term)
  phylo_deseq <- phylo_object %>% phyloseq_to_deseq2(as.formula(formula_parsed))
  phylo_deseq_geomeans <- apply(counts(phylo_deseq), 1, gm_mean)
  phylo_deseq = estimateSizeFactors(phylo_deseq, geoMeans = phylo_deseq_geomeans)
  phylo_deseq = DESeq(phylo_deseq, fitType="local")
  return(phylo_deseq)
}

annotate_deseq_results <- function(deseq_result, phylo_object){
  annotated <- deseq_result %>% 
    data.frame() %>% 
    arrange(padj) %>% 
    rownames_to_column(var="amplicon") %>% 
    left_join(tax_table(phylo_object) %>% data.frame() %>% rownames_to_column(var="amplicon"))
}



get_pos_neg_only_abundances <- function(da_annot){
  all_negative <- da_annot %>% 
    dplyr::select(log2FoldChange, Family) %>%
    group_by(Family) %>%
    summarise(all(log2FoldChange < 0)) %>%
    filter(`all(log2FoldChange < 0)` == TRUE) %>% 
    pull(Family)
  
  all_positive <- da_annot %>% 
    dplyr::select(log2FoldChange, Family) %>%
    group_by(Family) %>%
    summarise(all(log2FoldChange > 0)) %>%
    filter(`all(log2FoldChange > 0)` == TRUE) %>% 
    pull(Family)
  
  return(c(all_negative, all_positive))
}


################################################################
#             sediment differential abundance.                 #
################################################################

sediment_da <- subset_samples(phylo_elevation, Sample.description %in% c("Sediment marsh", "Sediment dry")) %>% run_deseq("0 +Sample.description")

#results(dds, contrast = c("condition", "treated", "untreated"))
#In this case, treated (numerator) is compared to untreated (denominator/baseline)
# so a negative fold change means that something is lower in dry compared to marsh
# so negative is seaward and positive is landweard
sediment_da_annot <- results(sediment_da, contrast = list("Sample.descriptionSediment.dry", "Sample.descriptionSediment.marsh")) %>% annotate_deseq_results(phylo_elevation)

sediment_neg_pos <- get_pos_neg_only_abundances(sediment_da_annot)

sediment_plot <- sediment_da_annot %>% 
  dplyr::select(log2FoldChange, Family) %>% 
  drop_na() %>%
  filter(Family %in% sediment_neg_pos & abs(log2FoldChange) > 2) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10) %>%
  mutate(whatever = case_when(log2FoldChange < 0 ~ "seaward",
                   log2FoldChange > 0 ~ "landward")) %>%
  ggplot(aes(y=reorder(Family, log2FoldChange), x=log2FoldChange, fill=whatever)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=17),
        legend.title = element_blank(),
        legend.position = "none") + 
  ylab("") +
  xlab("Log2 Fold Change") +
  scale_fill_manual(values = c("brown3", "dodgerblue"), labels=c("High Marsh", "Low Marsh"))


################################################################
#             root differential abundance.                 #
################################################################

root_da <- subset_samples(phylo_elevation, Sample.description %in% c("Root marsh", "Root dry")) %>% run_deseq("0 + Sample.description")

#results(dds, contrast = c("condition", "treated", "untreated"))
#In this case, treated (numerator) is compared to untreated (denominator/baseline)
# so a negative fold change means that something is lower in dry compared to marsh
# so negative is seaward and positive is landweard
root_da_annot <- results(root_da, contrast = list("Sample.descriptionRoot.dry", "Sample.descriptionRoot.marsh")) %>% annotate_deseq_results(phylo_elevation)

root_neg_pos <- get_pos_neg_only_abundances(root_da_annot)

root_plot <- root_da_annot %>% 
  dplyr::select(log2FoldChange, Family) %>% 
  drop_na() %>%
  filter(Family %in% root_neg_pos & abs(log2FoldChange) > 2 & Family != "Gammaproteobacteria Incertae Sedis Unknown Family") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10) %>%
  mutate(whatever = case_when(log2FoldChange < 0 ~ "seaward",
                              log2FoldChange > 0 ~ "landward")) %>%
  ggplot(aes(y=reorder(Family, log2FoldChange), x=log2FoldChange, fill=whatever)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=17),
        legend.title = element_blank()) + 
  ylab("") +
  xlab("Log2 Fold Change") +
  scale_fill_manual(values = c("brown3", "dodgerblue"), labels=c("High Marsh", "Low Marsh")) +
  xlim()


png("differential_abundance.png", height = 500, width=1500)
(sediment_plot | (root_plot))
dev.off()
  



