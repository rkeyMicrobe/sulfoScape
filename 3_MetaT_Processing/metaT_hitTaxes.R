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

##########################

# Read in taxonomy
pa_tax <- read.csv(paste0(res_dir, "G2PA_taxTable_Clean.csv")) %>% tibble
ns_tax <- read.csv(paste0(res_dir, "G2NS_taxTable_Clean.csv")) %>% tibble

# Read in hits, no e-value filter
pa_hit <- read.csv(paste0(dat_dir, "G2PA/all_G2PA_hmmHits.csv")) %>% tibble 
ns_hit <- read.csv(paste0(dat_dir, "G2NS/all_G2NS_hmmHits.csv")) %>% tibble

# Read in top-Hits, one contig per target
pa_top <- read.csv(paste0(dat_dir, "G2PA/G2PA_topHit_wTax.csv")) %>% tibble
ns_top <- read.csv(paste0(dat_dir, "G2NS/G2NS_topHit_wTax.csv")) %>% tibble

############################

# Make main dfs

fetchMaster = function(topDF = NULL, hitDF = NULL, taxDF = NULL){
    df =  topDF %>% 
          rename(TargetName = qseqid, taxID = taxid) %>% 
          select(TargetName, taxID) %>% 
          left_join(., hitDF %>% 
                         select(TargetName, QueryName, Evalue, BitScore, Sequence),
                    by = "TargetName") %>% 
          left_join(., taxDF %>% 
                         select(-cruise, -kingdom),
                    by = "taxID") %>% 
          relocate(., Sequence, .after = strain)
    return(df)
}

pa_df <- fetchMaster(topDF = pa_top, hitDF = pa_hit, taxDF = pa_tax)
ns_df <- fetchMaster(topDF = ns_top, hitDF = ns_hit, taxDF = ns_tax)

length(unique(pa_df$QueryName))
length(unique(ns_df$QueryName))

pa_df
ns_df

############################

# Split genes into own files for evalue eval

## G2PA
out_dir <- "results/metaT/hits_wTax/G2PA"
pa_splits = pa_df %>% select(TargetName, taxID, QueryName, Evalue,
                             BitScore, domain, class, Sequence) %>% 
     rename(contigID = TargetName, Gene = QueryName) %>% 
     group_split(Gene)
pa_splits %>%
     walk(function(d) {
          d <- d %>% arrange(Evalue)
          gene <- unique(d$Gene)
          gene_sane <- gsub("[^A-Za-z0-9_.-]", "_", gene)
          outfile <- file.path(out_dir, paste0(gene_sane, "_hits_wTax_G2PA.csv"))
          write_csv(d, outfile)
     })

## G2NS
out_dir <- "results/metaT/hits_wTax/G2NS"
ns_splits = ns_df %>% select(TargetName, taxID, QueryName, Evalue,
                             BitScore, domain, class, Sequence) %>% 
     rename(contigID = TargetName, Gene = QueryName) %>% 
     group_split(Gene)
ns_splits %>%
     walk(function(d) {
          d <- d %>% arrange(Evalue)
          gene <- unique(d$Gene)
          gene_sane <- gsub("[^A-Za-z0-9_.-]", "_", gene)
          outfile <- file.path(out_dir, paste0(gene_sane, "_hits_wTax_G2NS.csv"))
          write_csv(d, outfile)
     })

# Apply cut-offs
cutoffs <- 
     read.csv(paste0(dat_dir, "query_eValues_v2.csv")) %>% tibble %>% 
     select(queryName, eValue, eValue_original) %>%
     rename(Gene = queryName) %>% 
     mutate(cutoffs = case_when(is.na(eValue) ~ eValue_original, T ~ eValue) 
     ) %>% 
     select(Gene, cutoffs)


## G2NS Cut Offs
out_dir <- "results/metaT/hits_wTax/G2NS_cutoffs"

ns_splits %>%
     walk(function(d) {
          
          gene <- unique(d$Gene)
          
          d_filt <- d %>%
               left_join(cutoffs, by = "Gene") %>%
               filter(Evalue <= cutoffs) %>%
               arrange(Evalue)
          
          gene_sane <- gsub("[^A-Za-z0-9_.-]", "_", gene)
          outfile <- file.path(out_dir, paste0(gene_sane, "_hits_wTax_G2NS_cutoff.csv")
          )
          
          write_csv(d_filt, outfile)
     })

## G2PA Cut offs
out_dir <- "results/metaT/hits_wTax/G2PA_cutoffs"

pa_splits %>%
     walk(function(d) {
          
          gene <- unique(d$Gene)
          
          d_filt <- d %>%
               left_join(cutoffs, by = "Gene") %>%
               filter(Evalue <= cutoffs) %>%
               arrange(Evalue)
          
          gene_sane <- gsub("[^A-Za-z0-9_.-]", "_", gene)
          outfile <- file.path(out_dir, paste0(gene_sane, "_hits_wTax_G2PA_cutoff.csv")
          )
          
          write_csv(d_filt, outfile)
     })

############################

# Apply cut-offs
cutoffs <- 
     read.csv(paste0(dat_dir, "query_eValues_v2.csv")) %>% tibble %>% 
     select(queryName, eValue, eValue_original) %>%
     rename(QueryName = queryName) %>% 
     mutate(cutoffs = case_when(is.na(eValue) ~ eValue_original, T ~ eValue) 
     ) %>% 
     select(QueryName, cutoffs)

pa_cuts <- pa_df %>%
     inner_join(cutoffs, by = "QueryName") %>%
     filter(Evalue <= cutoffs) %>% 
     arrange(QueryName, Evalue) %>% 
     relocate(., cutoffs, .before = Evalue) %>% 
     mutate(TargetName = str_sub(TargetName, end = -3)) %>% 
# Best Gene per Contig
     group_by(TargetName) %>%
     slice_min(Evalue, n = 1, with_ties = FALSE) %>%
     ungroup()

ns_cuts <- ns_df %>%
     inner_join(cutoffs, by = "QueryName") %>%
     filter(Evalue <= cutoffs) %>% 
     arrange(QueryName, Evalue) %>% 
     relocate(., cutoffs, .before = Evalue) %>% 
# Best Gene per Contig
     group_by(TargetName) %>%
     slice_min(Evalue, n = 1, with_ties = FALSE) %>%
     ungroup()

# Compare number of contigs before and after cuts

before_counts <- ns_df %>%
     group_by(QueryName) %>%
     summarise(before = n_distinct(TargetName), .groups = "drop")
after_counts <- ns_cuts %>%
     group_by(QueryName) %>%
     summarise(after = n_distinct(TargetName), .groups = "drop")

contig_change <- before_counts %>%
     left_join(after_counts, by = "QueryName") %>%
     left_join(cutoffs,     by = "QueryName") %>%   # add evalue cutoff
     mutate(
          before = replace_na(before, 0),
          after  = replace_na(after, 0)
     ) %>%
     arrange(desc(after)) %>%
     mutate(QueryName = factor(QueryName, levels = unique(QueryName))) %>%
     pivot_longer(
          cols = c(before, after, cutoffs),
          names_to = "stage",
          values_to = "value"
     ) %>%
     mutate(stage = factor(stage, levels = c("before", "after", "cutoffs")))

ggplot(contig_change,
       aes(x = stage, y = QueryName)) +
     geom_tile(fill = "white", color = "black") +
     geom_text(aes(label = value), color = "black") +
     labs(x = NULL, y = "Gene", title = "G2NS. Number of Contigs") +
     theme_cowplot() +
     theme(legend.position = "none")

#-------------------------------------------------------
# Append Counts
## PA
list <- unique(pa_cuts$TargetName); head(list)
pa_cnt <- read.csv(paste0(dat_dir, "G2PA/G2PA_hitCounts.csv")) %>% tibble
head(pa_cnt$target_id)

pa_cnt <- pa_cnt %>% 
     filter(target_id %in% list) %>% 
     select(-length) %>% rename(TargetName = target_id)

nf <- read.csv("data/metaT/normFactors/NPac.G2PA.norm.factors.csv") %>% tibble %>% 
     rename(sampleID = sample_name) %>% 
     select(sampleID, NORM_FACTOR) %>%
     deframe()     
pa_cnt <- pa_cnt %>%
     mutate(across(
          .cols = -TargetName,                       
          .fns  = ~ .x * nf[cur_column()],         
          .names = "{.col}"                       
     ))

## NS
list <- unique(ns_cuts$TargetName); head(list)
ns_cnt <- read.csv(paste0(dat_dir, "G2NS/G2NS_hitCounts.csv")) %>% tibble
head(ns_cnt$nt_id)

ns_cnt <- ns_cnt %>% 
     filter(nt_id %in% list) %>% 
     select(-`Unnamed..0`, -nt_id_sample) %>% 
     rename(sampleID = sample_name, TargetName = nt_id) %>% 
     select(TargetName, sampleID, transcripts_L_A, transcripts_L_B, transcripts_L_C) %>% 
     mutate(sampleID = str_replace(sampleID, "\\.2um$", ".0_2um")) %>% 
     pivot_longer(., cols = -c(TargetName, sampleID), names_to = "rep", values_to = "counts") %>% 
     mutate(rep = str_extract(rep, "[A-Z]$"),
            sampleID = paste0(sampleID, ".", rep)) %>% 
     select(-rep) %>% 
     pivot_wider(., names_from = "sampleID", values_from = "counts") %>% 
     mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Check
pa_cnt
ns_cnt

pa_cuts 
ns_cuts

pa_master <- pa_cnt %>% 
     pivot_longer(., cols = -TargetName, names_to = "sampleID", values_to = "transcripts_L") %>% 
     left_join(., pa_cuts, by = "TargetName") %>% 
     mutate(taxID = as.character(taxID)) %>% 
     rename(contig = TargetName, gene = QueryName, eVal_cutOff = cutoffs)

ns_master <- ns_cnt %>% 
     pivot_longer(., cols = -TargetName, names_to = "sampleID", values_to = "transcripts_L") %>% 
     left_join(., ns_cuts, by = "TargetName") %>% 
     mutate(taxID = as.character(taxID)) %>% 
     rename(contig = TargetName, gene = QueryName, eVal_cutOff = cutoffs)

#########################################################################
#########################################################################
#########################################################################

tally <- pa_master %>%
     group_by(gene, class) %>%
     summarise(n_contigs = n_distinct(contig), .groups = "drop")

conc <- pa_master %>%
     filter(!is.na(class)) %>%
     group_by(gene, class) %>%
     summarise(transcripts_L = mean(transcripts_L, na.rm = TRUE),
               .groups = "drop")

make_gene_plot <- function(g) {
     
     g_tally <- tally %>%
          filter(gene == g) %>%
          arrange(desc(n_contigs)) %>%
          mutate(class = factor(class, levels = unique(class)))
     
     g_conc <- conc %>%
          filter(gene == g) %>%
          mutate(class = factor(class, levels = levels(g_tally$class)))
     
     plot_df <- g_tally %>%
          select(class, n_contigs) %>%
          left_join(
               g_conc %>% select(class, transcripts_L),
               by = "class"
          ) %>%
          pivot_longer(
               cols = c(n_contigs, transcripts_L),
               names_to = "metric",
               values_to = "value"
          )
     
     ggplot(plot_df, aes(x = class, y = value)) +
          geom_col() +
          coord_flip() +
          facet_wrap(~ metric, nrow = 1, scales = "free_x") +
          labs(x = "class", y = NULL, title = g) +
          theme_cowplot()
}

genes <- tally %>%
     distinct(gene) %>%
     pull(gene)

pdf("g2pa_plots.pdf", width = 10, height = 5)
for (g in genes) {
     print(make_gene_plot(g))
}
dev.off()


tally <- ns_master %>%
     group_by(gene, class) %>%
     summarise(n_contigs = n_distinct(contig), .groups = "drop")

conc <- ns_master %>%
     filter(!is.na(class)) %>%
     group_by(gene, class) %>%
     summarise(transcripts_L = sum(transcripts_L, na.rm = TRUE),
               .groups = "drop")

genes <- tally %>%
     distinct(gene) %>%
     pull(gene)

pdf("g2ns_plots.pdf", width = 10, height = 5)
for (g in genes) {
     print(make_gene_plot(g))
}
dev.off()
