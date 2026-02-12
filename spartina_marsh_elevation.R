



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
    filter(padj < 0.05) %>% 
    rownames_to_column(var="amplicon") %>% 
    left_join(tax_table(phylo_object) %>% data.frame() %>% rownames_to_column(var="amplicon"))
}


tst <- subset_samples(phylo_elevation, Sample.description %in% c("Sediment marsh", "Sediment dry")) %>% run_deseq("0 +Sample.description")

resultsNames(tst)
results(tst)

#results(dds, contrast = c("condition", "treated", "untreated"))
#In this case, treated (numerator) is compared to untreated (denominator/baseline)
# so a negative fold change means that something is lower in dry compared to marsh
# so negative is seaward and positive is landweard
tst_annot <- results(tst, contrast = list("Sample.descriptionSediment.dry", "Sample.descriptionSediment.marsh")) %>% annotate_deseq_results(phylo_elevation)




all_negative <- tst_annot %>% 
  dplyr::select(log2FoldChange, Family) %>%
  group_by(Family) %>%
  summarise(all(log2FoldChange < 0)) %>%
  filter(`all(log2FoldChange < 0)` == TRUE) %>% 
  pull(Family)

all_positive <- tst_annot %>% 
  dplyr::select(log2FoldChange, Family) %>%
  group_by(Family) %>%
  summarise(all(log2FoldChange > 0)) %>%
  filter(`all(log2FoldChange > 0)` == TRUE) %>% 
  pull(Family)


# this should come up with no overlap (sanity check)
intersect(all_negative, all_positive)

tst_annot %>% 
  dplyr::select(log2FoldChange, Family) %>% 
  drop_na() %>%
  filter(Family %in% c(all_negative, all_positive)) %>%
  mutate(whatever = case_when(log2FoldChange < 0 ~ "landward",
                   log2FoldChange > 0 ~ "seaward")) %>%
  ggplot(aes(x=reorder(Family, log2FoldChange), y=log2FoldChange, fill=whatever)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
