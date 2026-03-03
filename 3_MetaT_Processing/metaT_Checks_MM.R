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
library("feather")
dat_dir = "data/metaT/"
res_dir = "results/metaT/data/"
plo_dir = "results/metaT/figs/"
tab_dir = "results/metaT/tabs/"

cruise = "G2PA"
metaT_dir = "results/metaT/metaT_rentrez/"
d = paste0(dat_dir, cruise, "/")

# Load in Dataframes
hmmerHits <- read.csv(paste0(d, "all_", cruise, "_hmmHits.csv")) %>% tibble

cutoffs <- read.csv(paste0(dat_dir, "query_eValues_v2.csv")) %>% tibble %>% 
     select(queryName, eValue, eValue_original) %>% 
     distinct %>% 
     mutate(eValue_cutoffs = case_when(is.na(eValue) ~ eValue_original,
                                       T ~ eValue) 
     ) %>% 
     select(queryName, eValue_cutoffs)

tax <- read.csv(paste0(d, cruise, "_topHit_wTax_MM.csv")) %>% tibble

# Apply CutOffs
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
### Then choose the gene with the best evalue score
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

# Subset Taxonomy further... 
tax_clean <- tax %>% filter(qseqid %in% list) %>% 
     rename(contig = qseqid, taxID = taxid) %>% 
     select(contig, taxID) %>% distinct
length(list)
nrow(tax_clean)
tax_clean

################################################################################
################################################################################

# Efetch 
## Taxonomy Fetch
taxIDs <- tax_clean %>% select(taxID) %>% distinct

Sys.setenv(ENTREZ_KEY = "ffff45ab08d02129f62745870b2c36145607")
options(ENTREZ_KEY = "ffff45ab08d02129f62745870b2c36145607") 
rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY"))

# taxize
ids <- taxIDs$taxID %>%
     unique() %>%
     na.omit()
length(ids)

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

# Make Master, uncleaned
temp <- hmmerHits_Clean %>% select(contig, gene, Sequence) %>% distinct
temp <- tax_clean %>% select(contig, taxID) %>% distinct %>% 
     left_join(., temp, by =  "contig") %>% 
     left_join(., ncbi_taxa, by = "taxID")
ranks <- c("domain","kingdom","phylum","class","order","family","genus","species","strain")
N <- nrow(ncbi_taxa)
temp <- temp %>% mutate(across(all_of(ranks), ~na_if(., "")))


unclean <- temp

################################################################################
################################################################################

# Cleaning
# Amplicon 2019. Pr2v5 and Silva 138.1
asv_euk <- read_feather(paste0(dat_dir, "g123_18S_taxonomy.feather")) %>% tibble %>% 
     select(-featureID, -specific) %>% distinct
asv_pro <- read_feather(paste0(dat_dir, "g123_16S_taxonomy.feather")) %>% tibble %>% 
     select(-featureID, -specific) %>% distinct

# Genus helps family
cleanTax <- function(data = NULL, helperDF = NULL){
     
     data2 <- data %>%
          mutate(
               genus = if_else(
                    is.na(genus) &
                         !is.na(species) &
                         grepl("^[A-Z][a-z]+\\s", species),
                    sub("^([A-Z][a-z]+).*", "\\1", species),
                    genus
               )
          ) %>% 
          mutate(phylum = ifelse(phylum == "Candidatus", NA, phylum),
                 phylum = gsub("^Candidatus ", "", phylum))
     
     lookup <- helperDF %>% select(genus, species) %>% 
          distinct(species, genus) %>% 
          filter(!is.na(species)) %>% 
          rename(help = genus)
     temp <- data2 %>%
          left_join(., lookup, by = "species") %>%
          mutate(genus = if_else(is.na(genus) & !is.na(help), help, genus)
          )
     lookup <- helperDF %>% select(family, genus) %>% 
          distinct(genus, family) %>% 
          filter(!is.na(genus)) %>% 
          rename(help = family)
     temp <- temp %>% select(-help) %>% 
          left_join(., lookup, by = "genus") %>%
          mutate(family = if_else(is.na(family) & !is.na(help), help, family)
          )
     # Family helps Order
     lookup <- helperDF %>% select(order, family) %>% 
          distinct(family, order) %>% 
          filter(!is.na(family)) %>% 
          rename(help = order)
     temp <- temp %>% select(-help) %>% 
          left_join(., lookup, by = "family") %>%
          mutate(order = if_else(is.na(order) & !is.na(help), help, order)
          ) 
     # order helps class
     lookup <- helperDF %>% select(class, order) %>% 
          distinct(order, class) %>% 
          filter(!is.na(order)) %>% 
          rename(help = class)
     temp <- temp %>% select(-help) %>% 
          left_join(., lookup, by = "order") %>%
          mutate(class = if_else(is.na(class) & !is.na(help), help, class)
          ) 
     # class helps phylum
     lookup <- helperDF %>% select(phylum, class) %>% 
          distinct(class, phylum) %>% 
          filter(!is.na(class)) %>% 
          rename(help = phylum) %>% 
          arrange(class)
     temp <- temp %>% select(-help) %>% 
          left_join(., lookup, by = "class") %>%
          mutate(phylum = if_else(is.na(phylum) & !is.na(help), help, phylum)
          ) %>% 
          select(-help)  
     return(temp)
}
clean <- cleanTax(data = unclean, helperDF = asv_euk) 
clean <- cleanTax(data = clean, helperDF = asv_pro) 

################################################################################
################################################################################

# PLOT: MISSINGNESS ACROSS RANKS
N = nrow(clean)
plot_dat <- clean %>%
     select(contig, domain:strain) %>% 
     summarise(across(all_of(ranks), ~sum(is.na(.)))) %>%
     select(-kingdom) %>% 
     pivot_longer(everything(), names_to = "rank", values_to = "n_missing") %>%
     mutate(
          n_present = N - n_missing,
          pct_missing = round(100 * n_missing / N, 1)
     ) %>%
     mutate(rank = factor(rank, levels = ranks)) %>% 
     mutate(total = n_missing + n_present) 

total_n = unique(plot_dat$total)
type = cruise

p1 <- ggplot(plot_dat, aes(rank, pct_missing)) +
     geom_col() +
     geom_text(aes(label = paste0(pct_missing, "%")),
               vjust = -0.3, size = 3) +
     labs(x = "Rank", y = "% missing", 
          title = paste0(type, " Annotation Coverage Gaps"),
          subtitle = paste0("Total Number of NCBI TaxIDs: ", total_n)) +
     theme_classic() +
     coord_cartesian(ylim = c(0, max(plot_dat$pct_missing) * 1.12))
p1

# NUMBER ACROSS GENES, BEFORE AND AFTER CUT-OFFS
library(scales)
library(flextable)

all_genes = unique(cutoffs$queryName)

temp <- hmmerHits_Clean %>%
     distinct(contig, gene) %>%
     count(gene, name = "w_Cuts") %>%
     right_join(tibble(gene = all_genes), by = "gene") %>%
     replace_na(list(w_Cuts = 0)) %>%
     bind_rows(
          summarise(., gene = "Total", w_Cuts = sum(w_Cuts))
     )
temp2 <- hmmerHits %>%
     rename(gene = QueryName) %>% 
     distinct(TargetName, gene) %>%
     count(gene, name = "No_Cuts") %>%
     right_join(tibble(gene = all_genes), by = "gene") %>%
     replace_na(list(No_Cuts = 0)) %>%
     bind_rows(
          summarise(., gene = "Total", No_Cuts = sum(No_Cuts))
     )

tab_data <- temp2 %>% left_join(., temp, by = "gene") %>% 
     arrange(desc(w_Cuts)) %>% 
     filter(No_Cuts != 0) %>% 
     left_join(., cutoffs %>% 
                    rename(gene = queryName, eValue = eValue_cutoffs) %>% 
                    mutate(eValue = scientific(eValue)),
               by = "gene")

t1 <- tab_data %>% flextable() %>% 
     bg(bg = "white", part = "all") %>% 
     bg(i = ~ gene == "Total", bg = "#DCEEFF", part = "body") %>% 
     bold(i = ~ gene == "Total") %>% 
     align(j = -1, align = "right", part = "all") %>%
     set_caption(paste0(cruise, " Transcript Counts (BEFORE and AFTER e-Value Cuts)")) %>% 
     autofit()
t1

# NUMBER ACROSS GENES, BEFORE AND AFTER CUT-OFFS
clean

clean %>%
     drop_na(phylum) %>% 
     count(gene, phylum) %>%
     ggplot(aes(x = n, y = gene, fill = phylum)) +
     geom_col() +
     labs(
          x = "Gene",
          y = "Number of Contigs",
          fill = "Phylum",
          title = "Phylum Diversity Across Genes"
     ) +
     theme_cowplot()

# GENES LOST BETWEEN HMMERHIT AND TAX DATAFRAMES
temp <- hmmerHits_Clean %>% distinct(contig, gene) %>% 
     count(gene, name = "hmmerHits")
temp2 <- tax_clean %>% 
     left_join(., 
               hmmerHits_Clean %>% distinct(contig, gene),
               by = "contig") %>% 
     distinct(contig, gene) %>% 
     count(gene, name = "Annotations")

tab_data <- temp %>% left_join(., temp2, by = "gene")

t2 <- tab_data %>% flextable() %>% 
     bg(bg = "white", part = "all") %>% 
     align(j = -1, align = "right", part = "all") %>%
     set_caption(paste0(cruise, " Annotation Coverage")) %>% 
     autofit()
t2

################################################################################
################################################################################

# Save Files

write.csv(unclean, paste0(res_dir, cruise, "_uncleanTAX_MM.csv"), row.names = F)
write.csv(clean, paste0(res_dir, cruise, "_cleanTAX_MM.csv"), row.names = F)

svg(filename = paste0(plo_dir, cruise, "_numberNA_MM.svg"), width = 7, height = 5, bg = "transparent")
p1; dev.off()

save_as_image(t1, path = paste0(tab_dir, cruise, "_counts.png"))
save_as_image(t2, path = paste0(tab_dir, cruise, "_annotLosses.png"))
