lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD SCRIPT PACKAGES
library("tidyverse")
library("cowplot")
library("rentrez")
library("taxize")
library("ggalluvial")
dat_dir = "data/metaT/"
res_dir = "results/metaT/data"

cruise = "G2PA"
metaT_dir = "results/metaT/metaT_rentrez/"
d = paste0(dat_dir, cruise, "/")

hits <- read.csv(paste0(d, cruise, "_topHit_wTax_MM.csv")) %>% tibble %>% 
     rename(taxID = taxid, TargetName = qseqid)
hmmerHits <- read.csv(paste0(d, "all_", cruise, "_hmmHits.csv")) %>% tibble

cutoffs <- read.csv(paste0(dat_dir, "query_eValues_v2.csv")) %>% tibble %>% 
     select(queryName, eValue, eValue_original) %>% 
     distinct %>% 
     mutate(eValue_cutoffs = case_when(is.na(eValue) ~ eValue_original,
                                       T ~ eValue) 
     ) %>% 
     select(queryName, eValue_cutoffs)

# Quick View
hits_plot <- hits %>% mutate(mlog10_e = -log10(evalue))

ggplot(hits_plot, aes(mlog10_e)) +
     geom_histogram(bins = 60) +
     labs(x = "-log10(E-value)", y = "Count", title = "BLASTp E-values") +
     theme_classic()
hits %>%
     summarise(
          min_evalue = min(evalue, na.rm = TRUE),
          max_evalue = max(evalue, na.rm = TRUE),
          n_hits_allGenes     = n(),
          n_taxIDs   = n_distinct(taxID)
     )

hits %>%
     summarise(
          n_hits_allGenes = n(),
          n_taxIDs        = n_distinct(taxID),
          min_positive    = min(evalue[evalue > 0], na.rm = TRUE),
          max_evalue      = max(evalue, na.rm = TRUE)
     )

#----------------------------------------------------------------------------

# Apply filters to hmmerHits to reduce workload on hit tax ID table

## Apply eValue Cut-Offs
getHits <- function(data = NULL, cuts = NULL){
     filterVals <- cuts %>% tibble %>%
          rename(QueryName = queryName, e_cut = eValue_cutoffs) %>%
          select(QueryName, e_cut)
     
     df1 <- data %>% tibble %>% 
          inner_join(filterVals, by = "QueryName") %>%              
          mutate(pass = Evalue <= e_cut) %>% 
          select(TargetName, QueryName, Sequence, Evalue, e_cut, pass) %>% 
          rename(contig = TargetName, gene = QueryName) %>% 
          filter(pass == "TRUE")
     return(df1)
}
hmmerHits_Clean <- getHits(data = hmmerHits, cuts = cutoffs)

## Remove Redundancy: Parse out which contigs have > 2 genes associated
## Then choose the gene with the best evalue score
tempContigs <- hmmerHits_Clean %>% 
     distinct(contig, gene) %>%              
     group_by(contig) %>%
     summarise(
          n_genes = n(),                          
          genes   = paste(sort(gene), collapse = ", "),
          .groups = "drop"
     ) %>%
     filter(n_genes >= 2)
tempContigs

hmmerHits_Clean <- hmmerHits_Clean %>% 
     arrange(contig, Evalue) %>%   
     group_by(contig) %>%
     slice(1) %>%
     ungroup()

list <- unique(hmmerHits_Clean$contig)

## Apply list to the hit dataframe
hits_clean <- hits %>% filter(TargetName %in% list)

#----------------------------------------------------------------------------

# Taxonomy Fetch
taxIDs <- hits_clean %>% select(taxID) %>% distinct

Sys.setenv(ENTREZ_KEY = "ffff45ab08d02129f62745870b2c36145607")
options(ENTREZ_KEY = "ffff45ab08d02129f62745870b2c36145607") 
rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))

# taxize
ids <- taxIDs$taxID %>%
     unique() %>%
     na.omit()

safe_one <- function(id) {
     out <- try(classification(id, db = "ncbi"), silent = TRUE)
     
     # if the query fails or returns nothing, return an NA row with this taxID
     if (inherits(out, "try-error") || length(out) == 0 || !is.data.frame(out[[1]])) {
          return(tibble(
               name  = NA_character_,
               rank  = NA_character_,
               id    = NA_character_,
               taxID = as.integer(id)
          ))
     }
     
     df <- out[[1]]
     df$taxID <- as.integer(id)
     as_tibble(df)
}

# one taxonomy df per taxID, no flattening / nesting issues
cl_df <- purrr::map_dfr(ids, safe_one)

pick1 <- function(x) if (length(x)) x[1] else NA_character_

ncbi_taxa <- cl_df %>%
     mutate(rank = tolower(rank)) %>%
     group_by(taxID) %>%
     summarise(
          domain  = pick1(name[rank %in% c("superkingdom", "domain")]),
          kingdom = pick1(name[rank == "kingdom"]),
          phylum  = pick1(name[rank %in% c("phylum", "division")]),
          class   = pick1(name[rank == "class"]),
          order   = pick1(name[rank == "order"]),
          family  = pick1(name[rank == "family"]),
          genus   = pick1(name[rank == "genus"]),
          species = pick1(name[rank == "species"]),
          strain  = pick1(name[rank == "strain"]),
          .groups = "drop"
     )
ncbi_taxa # Go on NCBI Ta Browser, check taxID to make sure they align

#----------------------------------------------------------------------------


seqtemp <- hmmerHits_Clean %>% select(contig, Sequence) %>% distinct
hittemp <- hits_clean %>% select(TargetName, taxID) %>% distinct %>% 
     left_join(., seqtemp, by = c("TargetName" = "contig"))

#####



write.csv(ncbi_taxa %>% distinct %>% mutate(cruise = cruise), 
          paste0(metaT_dir, cruise, "_unclean_taxTable.csv"), row.names = F)

# Finish up
ncbi_taxa <- read.csv(paste0(metaT_dir, cruise, "_unclean_taxTable.csv"))

ranks <- c("domain","kingdom","phylum","class","order","family","genus","species","strain")
N <- nrow(ncbi_taxa)

ncbi_taxa2 <- ncbi_taxa %>% mutate(across(all_of(ranks), ~na_if(., "")))

ncbi_taxa2

write.csv(ncbi_taxa2 %>% distinct %>% mutate(cruise = cruise), 
          paste0(metaT_dir, cruise, "_taxTable.csv"), row.names = F)

#----------------------------------------------------------------------------
# What is missing at each level?

na_summary <- ncbi_taxa2 %>%
     summarise(across(all_of(ranks), ~sum(is.na(.)))) %>%
     pivot_longer(everything(), names_to = "rank", values_to = "n_missing") %>%
     mutate(
          n_present = N - n_missing,
          pct_missing = round(100 * n_missing / N, 1)
     )

plot_dat <- na_summary %>%
     mutate(rank = factor(rank, levels = ranks))

ggplot(plot_dat, aes(rank, pct_missing)) +
     geom_col() +
     geom_text(aes(label = paste0(pct_missing, "%")),
               vjust = -0.3, size = 3) +
     labs(x = "Rank", y = "% missing", title = "Annotation Coverage Gaps") +
     theme_classic() +
     coord_cartesian(ylim = c(0, max(plot_dat$pct_missing) * 1.12))

plot_dat2 <- na_summary %>%
     select(rank, n_missing, n_present) %>%
     pivot_longer(-rank, names_to = "status", values_to = "n") %>%
     mutate(rank = factor(rank, levels = ranks),
            status = factor(status, levels = c("n_present","n_missing"),
                            labels = c("Present","Missing")))

ggplot(plot_dat2, aes(rank, n, fill = status)) +
     geom_col(position = "dodge") +
     labs(x = "Rank", y = "Count", title = "Taxonomy annotations by rank") +
     theme_classic()

#----------------------------------------------------------------------------
# Apply hmmer thresholds
blastHits <- hits

ncbi_taxa2
unique(hmmerHits$QueryName)
cut_tab <- cutoffs %>%
     rename(QueryName = queryName, e_cut = eValue_cutoffs) %>%
     select(QueryName, e_cut)

hmmerHits_scored <- hmmerHits %>%
     inner_join(cut_tab, by = "QueryName") %>%              # only genes with a cutoff
     mutate(pass = Evalue <= e_cut)

# keep only hits that meet the cutoff
hmmerHits_filt <- hmmerHits_scored %>%
     filter(pass) %>%
     select(-pass)

cutoff_summary <- hmmerHits_scored %>%
     count(QueryName, pass) %>%
     tidyr::pivot_wider(names_from = pass, values_from = n,
                        values_fill = 0, names_prefix = "n_") %>%
     rename(n_fail = n_FALSE, n_pass = n_TRUE) %>% 
     mutate(cruise = cruise)
cutoff_summary

write.csv(cutoff_summary, paste0(metaT_dir, cruise, "_cutOff_summary_v2.csv"), row.names = F)


#----------------------------------------------------------------------------
# Revisit old plots... 
hmmerHits <- hmmerHits_filt %>% 
     select(TargetName, QueryName) %>% distinct
hits <- blastHits %>% left_join(., ncbi_taxa, by = "taxID") %>% 
     select(TargetName, taxID, domain:strain)

hits <- hits %>% left_join(., hmmerHits, by = "TargetName") %>% 
     drop_na(QueryName)

hits %>%
     count(QueryName, name = "n_rows") %>%
     arrange(desc(n_rows))

hits %>%
     filter(domain == "Eukaryota") %>% 
     mutate(phylum = if_else(is.na(phylum) | phylum == "NA", "Unassigned", as.character(phylum))) %>%
     group_by(QueryName, phylum) %>%
     summarise(n = n(), .groups = "drop") %>%
     ggplot(aes(x = QueryName, y = n, fill = phylum)) +
     geom_col(position = "stack") +
     scale_fill_viridis_d() +
     labs(x = "Gene", y = "Number of hits", fill = "Eukaryote Phylum") +
     theme_classic(12) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

flows <- hits %>%
     mutate(
          domain = if_else(is.na(phylum) | domain=="NA","Unassigned", domain),
          phylum = if_else(is.na(phylum) | phylum=="NA","Unassigned", phylum),
          class  = if_else(is.na(class)  | class =="NA","Unassigned", class),
          order  = if_else(is.na(order)  | order =="NA","Unassigned", order),
          family  = if_else(is.na(family)  | order =="NA","Unassigned", family),
          genus  = if_else(is.na(genus)  | order =="NA","Unassigned", genus),
          species  = if_else(is.na(species)  | order =="NA","Unassigned", species)
     ) %>% 
     filter(phylum != "Unassigned") %>% 
     count(QueryName, taxID, domain, phylum, class, order, family, genus, species, name = "n") 


p <- flows %>% filter(domain == "Bacteria") %>% 
     filter(QueryName == "ClaA") %>% 
     ggplot(.,
            aes(axis1 = QueryName, axis2 = phylum, axis3 = class, 
                axis4 = order, axis5 = family, axis6 = genus,
                y = n)) +
     geom_alluvium(aes(fill = QueryName), alpha = 0.7, width = .2) +
     geom_stratum(width = .2, color = "grey40", aes(fill = after_stat(stratum))) +
     geom_text(stat = "stratum", size = 3, color = "black",
               aes(label = after_stat(stratum)), min.y = 3) +
     scale_x_discrete(limits = c("Gene", "Phylum","Class"), expand = c(.05,.05)) +
     guides(fill = "none") +
     labs(y = "Hits") +
     theme_classic(12)
p


x <- flows %>% filter(domain == "Bacteria") %>% 
     filter(QueryName == "ClaA") %>% 
     select(-domain)




