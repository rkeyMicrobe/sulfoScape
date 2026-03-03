lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD SCRIPT PACKAGES
library("tidyverse")
library("cowplot")
library("ggalluvial")
library("feather")

dat_dir = "data/metaT/"
res_dir = "results/metaT/metaT_rentrez/"

# Load in Datasets

## Taxonomy
attach_tax <- function(top_df, tax_df) {
     top_df %>%
          select(qseqid, taxid) %>%
          left_join(tax_df, by = c("taxid" = "taxID")) %>%
          select(-kingdom, -cruise) %>%
          rename(contig = qseqid)
}
pa_tax <- attach_tax(read.csv(paste0(dat_dir, "G2PA/G2PA_topHit_wTax.csv")) %>% tibble, 
                     read.csv(paste0(res_dir, "G2PA_taxTable_Clean.csv")) %>% tibble)
ns_tax <- attach_tax(read.csv(paste0(dat_dir, "G2NS/G2NS_topHit_wTax.csv")) %>% tibble, 
                     read.csv(paste0(res_dir, "G2NS_taxTable_Clean.csv")) %>% tibble)

## Hits w/ cutoffs considered
getHits <- function(data = NULL){
     cuts <- read.csv(paste0(dat_dir, "/query_eValues.csv")) %>% tibble %>%
          rename(QueryName = queryName, e_cut = eValue_cutoffs) %>%
          select(QueryName, e_cut)
     
     df1 <- data %>% tibble %>% 
          inner_join(cuts, by = "QueryName") %>%              
          mutate(pass = Evalue <= e_cut) %>% 
          select(TargetName, QueryName, Sequence, Evalue, e_cut, pass) %>% 
          rename(contig = TargetName, gene = QueryName) %>% 
          filter(pass == "TRUE")
     return(df1)
}
pa_hit <- getHits(data = read.csv(paste0(dat_dir, "G2PA/all_G2PA_hmmHits.csv"))) 
ns_hit <- getHits(data = read.csv(paste0(dat_dir, "G2NS/all_G2NS_hmmHits.csv")))

### Apply taxonomy onto hit info
pa_hit <- pa_hit %>% 
     left_join(., pa_tax, by = "contig") %>% 
     select(-Sequence) %>% 
     mutate(contig = sub("_[0-9]$", "", contig))
ns_hit <- ns_hit %>% 
     left_join(., ns_tax, by = "contig") %>% 
     select(-Sequence) %>% 
     mutate(contig = sub("_[0-9]$", "", contig))

## Counts
list <- unique(pa_hit$contig)
pa_cnt <- read.csv(paste0(dat_dir, "G2PA/G2PA_hitCounts.csv")) %>% tibble %>% 
     filter(target_id %in% list) %>% 
     select(-length) %>% rename(contig = target_id)
list <- unique(ns_hit$contig)
ns_cnt <- read.csv(paste0(dat_dir, "G2NS/G2NS_hitCounts.csv")) %>% tibble %>% 
     filter(nt_id %in% list) %>% 
     select(-`Unnamed..0`, -nt_id_sample) %>% 
     rename(sampleID = sample_name, contig = nt_id)

### Clean NS to match PA structure
ns_cnt <- ns_cnt %>% 
     select(contig, sampleID, transcripts_L_A, transcripts_L_B, transcripts_L_C) %>% 
     mutate(sampleID = str_replace(sampleID, "\\.2um$", ".0_2um")) %>% 
     pivot_longer(., cols = -c(contig, sampleID), names_to = "rep", values_to = "counts") %>% 
     mutate(rep = str_extract(rep, "[A-Z]$"),
            sampleID = paste0(sampleID, ".", rep)) %>% 
     select(-rep) %>% 
     pivot_wider(., names_from = "sampleID", values_from = "counts") %>% 
     mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
### Apply norm factors to PA
nf <- read.csv("data/metaT/normFactors/NPac.G2PA.norm.factors.csv") %>% tibble %>% 
     rename(sampleID = sample_name) %>% 
     select(sampleID, NORM_FACTOR) %>%
     deframe()     
pa_cnt <- pa_cnt %>%
     mutate(across(
          .cols = -contig,                       
          .fns  = ~ .x * nf[cur_column()],         
          .names = "{.col}"                       
     ))

#------------------------------------------------------------------------------------

pa_hit 
ns_hit

# What genes contain redudant contigs? 
make_gene_tables <- function(hit_tbl, dataset_label) {
     
     contig_gene_summary <- hit_tbl %>%
          distinct(contig, gene) %>%              
          group_by(contig) %>%
          summarise(
               n_genes = n(),                          
               genes   = paste(sort(gene), collapse = ", "),
               .groups = "drop"
          ) %>%
          filter(n_genes >= 2)
     
     gene_table <- contig_gene_summary %>%
          separate_rows(genes, sep = ", ") %>%
          count(genes, name = "n_multi_contigs") %>%
          arrange(desc(n_multi_contigs)) %>%
          mutate(dataset = dataset_label, .before = 1)
     
     combo_table <- contig_gene_summary %>%
          count(genes, name = "n_contigs_with_combo") %>%
          arrange(desc(n_contigs_with_combo)) %>%
          mutate(dataset = dataset_label, .before = 1)
     
     list(
          contig_gene_summary = contig_gene_summary,
          gene_table          = gene_table,
          combo_table         = combo_table
     )
}
pa_res <- make_gene_tables(pa_hit, "PA")
pa_res$gene_table
pa_res$combo_table

ns_res <- make_gene_tables(ns_hit, "NS")
ns_res$gene_table
data.frame(ns_res$combo_table)

# Pick "best" gene per contig (ie. Lowest Eval)
temp_best <- pa_hit %>%
     arrange(contig, Evalue) %>%   
     group_by(contig) %>%
     slice(1) %>%
     ungroup()
list <- unique(temp_best$contig)
pa_master <- pa_cnt %>% 
     pivot_longer(., cols = -contig, names_to = "sampleID", values_to = "transcripts_L") %>% 
     left_join(., temp_best, by = "contig") %>% 
     filter(!is.na(gene)) 
  
temp_best <- ns_hit %>%
     arrange(contig, Evalue) %>%   
     group_by(contig) %>%
     slice(1) %>%
     ungroup()
ns_master <- ns_cnt %>% 
     pivot_longer(., cols = -contig, names_to = "sampleID", values_to = "transcripts_L") %>% 
     left_join(., temp_best, by = "contig") %>% 
     filter(!is.na(gene)) 

write_feather(pa_master, "results/metaT/pa_master.feather")
write_feather(ns_master, "results/metaT/ns_master.feather")

################################################################################
################################################################################

# Dominance across samples
df1 <- pa_master %>%
     group_by(gene) %>%
     summarize(total = sum(transcripts_L)) %>%
     arrange(desc(total)) %>%
     mutate(rank = row_number())
df2 <- ns_master %>%
     group_by(gene) %>%
     summarize(total = sum(transcripts_L)) %>%
     arrange(desc(total)) %>%
     mutate(rank = row_number())

p1 <- ggplot(df1, aes(rank, total, label = gene)) +
     geom_point() +
     geom_text(angle = 45, hjust = -.4, nudge_y = 0.1, size = 2, check_overlap = TRUE) +
     #scale_y_log10() +
     labs(title = "Dominance Ranks Across Smps", subtitle = "G2PA") +
     theme_classic()
p2 <- ggplot(df2, aes(rank, total, label = gene)) +
     geom_point() +
     geom_text(angle = 45, hjust = -.4, nudge_y = 0.1, size = 2, check_overlap = TRUE) +
     #scale_y_log10() +
     labs(title = " ", subtitle = "G2NS") +
     theme_classic()
plot_grid(p1, p2)

# Who's more dominant in which fraction? 
dominance_compare <- df1 %>%
     select(gene, total_pa = total) %>%
     inner_join(df2 %>% select(gene, total_ns = total), by = "gene") %>%
     mutate(pa_ns_ratio = log10(total_pa / total_ns)) %>%
     arrange(pa_ns_ratio)
ggplot(dominance_compare, aes(x = reorder(gene, pa_ns_ratio),
                              y = pa_ns_ratio,
                              fill = pa_ns_ratio > 0)) +
     geom_col() +
     coord_flip() +
     scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue3"),
                       labels = c("NS-dominant", "PA-dominant")) +
     labs(x = "gene",
          y = "log10(PA / NS)",
          fill = "Dominance",
          title = "PA vs NS Comparison, Bias") +
     theme_classic()

################################################################################
################################################################################
library(RColorBrewer)

plotTaxRank <- function(data = NULL,
                        type_cruise = NULL,
                        type_domain = NULL,
                        rank_col = "phylum",
                        N = 12) {
     
     rank_var <- sym(rank_col)   # turn string into symbol for tidy eval
     
     # Top N taxa at the chosen rank
     topN <- data %>%
          filter(!is.na(!!rank_var)) %>% 
          filter(domain == type_domain) %>%
          group_by(!!rank_var) %>%
          summarize(total = sum(transcripts_L), .groups = "drop") %>%
          arrange(desc(total)) %>%
          slice_head(n = N) %>%
          pull(!!rank_var) %>%
          as.character()
     
     # Palette (10 taxa + Other = black)
     pal <- c(
          "#BC80BD",
          "#BEBADA",
          "#FB8072",
          "#8DD3C7",
          "#80B1D3",
          "#B3DE69",
          "#FCCDE5",
          "#FDB462",
          "#CCEBC5",
          "#FFED6F",
          "brown",
          "gray",
          "black"
     )
     names(pal) <- c(topN, "Other")
     
     # Collapse to Top N + Other at chosen rank
     tax_topN <- data %>%
          filter(!is.na(!!rank_var)) %>% 
          filter(domain == type_domain) %>%
          mutate(
               rank2 = ifelse(
                    as.character(!!rank_var) %in% topN,
                    as.character(!!rank_var),
                    "Other"
               )
          ) %>%
          group_by(gene, rank2) %>%
          summarize(total = sum(transcripts_L), .groups = "drop")
     
     # Order stack within each gene by abundance
     tax_topN2 <- tax_topN %>%
          group_by(gene) %>%
          mutate(rank2 = fct_reorder(rank2, total, .fun = identity, .desc = TRUE)) %>%
          ungroup()
     
     # Nicely formatted label for legend
     rank_label <- tools::toTitleCase(rank_col)
     
     ggplot(tax_topN2, aes(x = gene, y = total, fill = rank2)) +
          geom_col() +
          scale_fill_manual(values = pal) +
          theme_cowplot() +
          labs(title = paste0("Top ", rank_label, " across Genes, ", type_cruise),
               y = "Total transcripts",
               x = "Gene",
               fill = rank_label) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
p1 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Bacteria", rank_col = "phylum")
p2 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Eukaryota", rank_col = "phylum")
plot_grid(p1, p2, ncol = 1)

p1 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Bacteria", rank_col = "class")
p2 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Eukaryota", rank_col = "class")
plot_grid(p1, p2, ncol = 1)

p1 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Bacteria", rank_col = "order")
p2 <- plotTaxRank(data = pa_master, type_cruise = "G2PA", type_domain = "Eukaryota", rank_col = "order")
plot_grid(p1, p2, ncol = 1)


p1 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Bacteria", rank_col = "phylum")
p2 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Eukaryota", rank_col = "phylum")
plot_grid(p1, p2, ncol = 1)

p1 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Bacteria", rank_col = "class")
p2 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Eukaryota", rank_col = "class")
plot_grid(p1, p2, ncol = 1)

p1 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Bacteria", rank_col = "order")
p2 <- plotTaxRank(data = ns_master, type_cruise = "G2NS", type_domain = "Eukaryota", rank_col = "order")
plot_grid(p1, p2, ncol = 1)




















