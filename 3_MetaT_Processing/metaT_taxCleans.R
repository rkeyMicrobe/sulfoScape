lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD SCRIPT PACKAGES
library("tidyverse")
library("cowplot")
library("feather")

dat_dir = "data/metaT/"
res_dir = "results/metaT/metaT_rentrez/"

#------------------------------------------------------------------------------

# Rentrez Tax Tables 
pa <- read.csv(paste0(res_dir, "G2PA_taxTable.csv")) %>% tibble
ns <- read.csv(paste0(res_dir, "G2NS_taxTable.csv")) %>% tibble

# Amplicon 2019. Pr2v5 and Silva 138.1
asv_euk <- read_feather(paste0(dat_dir, "g123_18S_taxonomy.feather")) %>% tibble %>% 
     select(-featureID, -specific) %>% distinct
asv_pro <- read_feather(paste0(dat_dir, "g123_16S_taxonomy.feather")) %>% tibble %>% 
     select(-featureID, -specific) %>% distinct

# Make changes, step-by-step...
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
#pa_clean <- cleanTax(data = pa, helperDF = zk) 
pa_clean <- cleanTax(data = pa, helperDF = asv_euk) 
pa_clean <- cleanTax(data = pa_clean, helperDF = asv_pro) 
sort(unique(pa_clean$phylum))
sort(unique(pa$phylum))


#ns_clean <- cleanTax(data = ns, helperDF = zk) 
ns_clean <- cleanTax(data = ns, helperDF = asv_euk) 
ns_clean <- cleanTax(data = ns_clean, helperDF = asv_pro) 

# Check corrections...
df1 <- pa_clean %>% filter(if_any(-c(strain, kingdom), is.na)) %>% select(-kingdom)
df2 <- ns_clean %>% filter(if_any(-c(strain, kingdom), is.na)) %>% select(-kingdom)

df1
df2

df <- rbind(df1, df2) %>% select(-cruise) %>% distinct

# Write out to do manual corrections...
write.csv(df, paste0(res_dir, "G2_taxCorrections.csv"), row.names = F)

# Table Regens...

checkGaps <- function(data = NULL, type = NULL){

     df <- data %>% select(-kingdom)
     ranks <- c("domain","phylum","class","order","family","genus","species","strain")
     N <- nrow(df)
     
     df2 <- df %>% mutate(across(all_of(ranks), ~na_if(., "")))
     
     na_summary <- df %>%
          summarise(across(all_of(ranks), ~sum(is.na(.)))) %>%
          pivot_longer(everything(), names_to = "rank", values_to = "n_missing") %>%
          mutate(
               n_present = N - n_missing,
               pct_missing = round(100 * n_missing / N, 1)
          )
     
     plot_dat <- na_summary %>%
          mutate(rank = factor(rank, levels = ranks)) %>% 
          mutate(total = n_missing + n_present)
     
     total_n = unique(plot_dat$total)
     
     ggplot(plot_dat, aes(rank, pct_missing)) +
          geom_col() +
          geom_text(aes(label = paste0(pct_missing, "%")),
                    vjust = -0.3, size = 3) +
          labs(x = "Rank", y = "% missing", 
               title = paste0(type, " Annotation Coverage Gaps"),
               subtitle = paste0("Total Number of NCBI TaxIDs: ", total_n)) +
          theme_classic() +
          coord_cartesian(ylim = c(0, max(plot_dat$pct_missing) * 1.12))
     
}
p1 <- checkGaps(data = pa_clean, type = "G2PA")
p2 <- checkGaps(data = ns_clean, type = "G2NS")
plot_grid(p1, p2, ncol = 1)

###############################################################################

# Save cleaned dataframes
write.csv(pa_clean, paste0(res_dir, "G2PA_taxTable_Clean.csv"), row.names = F)
write.csv(ns_clean, paste0(res_dir, "G2NS_taxTable_Clean.csv"), row.names = F)



