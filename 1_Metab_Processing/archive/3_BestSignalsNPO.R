lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())

library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v6.csv") %>% tibble 
iSTDs <- library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")

################################################################################
# Handle ATL_GOM datasets
################################################################################

hp <- read.csv("data/NPO_G2/QCd/QEQC_HILIC_Pos5.4_2017-2019_allSulfurStds.csv")
hn <- read.csv("data/NPO_G2/QCd/QEQC_HILIC_Neg5.4_2017-2019_allSulfurStds.csv")
rp <- read.csv("data/NPO_G2/QCd/QEQC_RP_Aq5.1_2017-2019_allSulfurStds3.csv")

cleanQC <- function(data = NULL){
     preData <- data
     colnames(preData) <- as.character(preData[1, ])
     preData <- preData[-1, ] 
     colnames(preData)[colnames(preData) == "Mass.Feature"] <- "Compound_Name"
     
     df <- preData %>% tibble %>% 
          select(Replicate.Name, Compound_Name, Area, Area.with.QC) %>% 
          mutate(sampleID = substr(Replicate.Name, 8, nchar(Replicate.Name)),
                 Area = as.numeric(Area),
                 Area.with.QC = as.numeric(Area.with.QC)) %>%
          rename(Compound_Name = Compound_Name, 
                 areaRAW = Area,
                 areaQCd = Area.with.QC) %>% 
          select(-Replicate.Name) %>% 
          mutate(sampleGroup = str_extract(sampleID, "(?<=_)[^_]+(?=_)")) %>% 
          tibble
     
     df <- df %>% filter(str_detect(sampleID, "Std_"))
     
     cat("Data processing complete.")
     return(df)
}

hp <- cleanQC(hp)
hn <- cleanQC(hn)
rp <- cleanQC(rp)

# Handle Names
cleanNames <- function(data = NULL, type = NULL, Mode = NULL){
     library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v5.csv") %>% 
          filter(Column == type & ionMode == Mode)
     nameKey <- read.csv("data/1_metaData/nameKey_update.csv") %>% 
          filter(Durham_Name %in% unique(library$Compound_Name)) %>% distinct
     
     data %>% distinct %>% 
          left_join(., nameKey, by = c("Compound_Name" = "name")) %>% 
          relocate("Durham_Name", .before = everything()) %>% tibble %>% 
          rename(oldName = Compound_Name,
                 Compound_Name = Durham_Name) %>% 
          tibble
}
hp <- cleanNames(data = hp, type = "HILIC", Mode = "Pos") %>% 
     mutate(Compound_Name = case_when(
          oldName == "(3-Carboxypropyl)trimethylammonium (TMAB)" ~ 
               "(3-Carboxypropyl)trimethylammonium (TMAB)",
          oldName == "Acetyl-L-carnitine" ~ "Acetyl-L-Carnitine",
          T ~ Compound_Name
     ))

hn <- cleanNames(data = hn, type = "HILIC", Mode = "Neg") 

rp <- cleanNames(data = rp, type = "RP", Mode = "Pos") %>% 
     mutate(Compound_Name = case_when(
          oldName == "Betaine" ~ "Betaine)",
          oldName == "Choline" ~ "Choline",
          T ~ Compound_Name
     ))


hp %>% filter(is.na(Compound_Name)) %>% select(Compound_Name, oldName) %>% distinct
hn %>% filter(is.na(Compound_Name)) %>% select(Compound_Name, oldName) %>% distinct
rp %>% filter(is.na(Compound_Name)) %>% select(Compound_Name, oldName) %>% distinct

################################################################################
# Consider MixSets for HILICs
################################################################################

mixClean <- function(data = NULL) {
     mixKey <- read.csv("data/1_metaData/metab_MixSets.csv") %>% tibble %>% 
          select(Compound_Name, mixSet_2019) %>% 
          drop_na(mixSet_2019) %>% distinct
     
     data %>%  drop_na(Compound_Name) %>% 
          left_join(mixKey, by = "Compound_Name") %>%
          filter(
               str_detect(sampleID, "Std_H2OinMatrix") |
                    (str_detect(sampleID, "Mix1") & mixSet_2019 == "Mix1") |
                    (str_detect(sampleID, "Mix2") & mixSet_2019 == "Mix2")
          )
}
hp <- mixClean(data = hp)
hn <- mixClean(data = hn)

################################################################################
# Handle ATL_GOM datasets
################################################################################

getSums <- function(data = NULL, sepType = NULL){
     x <- data %>% drop_na(Compound_Name) %>% 
          select(sampleID, Compound_Name, areaQCd) %>% 
          group_by(Compound_Name) %>% 
          reframe(total = mean(areaQCd, na.rm = T)) %>% 
          mutate(method = sepType) %>% 
          rename(metab = Compound_Name)
     return(x)
}

bests <- rbind(getSums(data = hn, sepType = "HN"),
               getSums(data = hp, sepType = "HP"),
               getSums(data = rp, sepType = "RP")) %>% 
     pivot_wider(names_from = method, values_from = total) %>% 
     mutate(across(everything(), ~replace(., is.infinite(.), 0))) %>% 
     filter(!metab %in% iSTDs) %>% 
     mutate(
          bestSignal = pmax(HN, HP, RP, na.rm = TRUE), # Highest signal
          detectionFreq = rowSums(!is.na(cbind(HN, HP, RP))), # Non-missing modes
          bestMode = case_when(
               bestSignal == HN ~ "HN",
               bestSignal == HP ~ "HP",
               bestSignal == RP ~ "RP",
               TRUE ~ "NA"
          )
     ) %>%
     arrange(desc(bestSignal))

# Consider Overrides, based on Ingall Library Notes and Heal Ratios from mSystems
applyOverrides <- function(data = NULL) {
     x <- data %>%
          mutate(bestMode = case_when(
               metab == "Cysteinylgylcine (Cys-Gly)" ~ "RP",
               metab == "Glutathione" ~ "HP",
               metab == "Glutathione Disulfide, Oxidized" ~ "HP",
               metab == "Choline" ~ "HP",
               metab == "Cysteinesulfinic Acid" ~ "HN",
               metab == "Butyryl-L-Carnitine" ~ "HP",
               metab == "Thiamine Monophosphate (Vit B1 Phosphate)" ~ "HP",
               metab == "Allopurinol" ~ "HN",
               metab == "Deoxyadenosine" ~ "HP", 
               metab == "Imidazoleacrylic Acid" ~ "HN",
               metab == "Inosine" ~ "HN",
               TRUE ~ bestMode
          ))
     return(x)
}

 applyOverrides <- function(data = NULL) {
     x <- data %>%
          mutate(bestMode = case_when(
               metab == "Cysteinylgylcine (Cys-Gly)" ~ "RP",
               metab == "Glutathione" ~ "HP",
               metab == "Glutathione Disulfide, Oxidized" ~ "HP",
               metab == "Choline" ~ "HP",
               metab == "Cysteinesulfinic Acid" ~ "HN",
               metab == "Butyryl-L-Carnitine" ~ "HP",
               metab == "Thiamine Monophosphate (Vit B1 Phosphate)" ~ "HP",
               TRUE ~ bestMode
          ))
     return(x)
}

bests_cleaned <- applyOverrides(data = bests) %>% 
     arrange(metab)

# Save info
dir = "data/ATL_GOM/BMIS_ATL_GOM_"
bests_cleaned
write.csv(bests_cleaned, paste0(dir, "signal_index.csv"), row.names = F)

################################################################################
# Sanity Check: Heal Comparisons
################################################################################
best_cleaned = bests_cleaned
hb <- read.csv("data/1_metaData/heal_ratios/NPO_RFs_Heal.csv") %>%
     select(Precursor.Ion.Name, ionMode) %>% 
     rename(name = Precursor.Ion.Name)

names <- read.csv("data/1_metaData/nameKey_update.csv") %>% 
     filter(name %in% unique(hb$name)) %>% distinct
rk <- best_cleaned %>% select(metab, bestMode) %>% distinct
cmpds <- rk$metab

hb <- hb %>% 
     left_join(., names, by = "name") %>% 
     relocate("Durham_Name", .before = everything()) %>% tibble %>% 
     rename(oldName = name,
            metab = Durham_Name) %>% 
     filter(metab %in% cmpds) %>% 
     tibble



 # mode comparisons

match <- hb %>% left_join(., rk, by = "metab") %>% 
     mutate(match = ifelse(is.na(ionMode) | is.na(bestMode), NA, bestMode == ionMode))

match %>% filter(match == FALSE)



