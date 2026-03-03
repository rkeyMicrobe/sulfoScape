lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")


library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble 
iSTDs <- library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)
sigIndex <- read.csv("data/1_metaData/bestSignals_index.csv") %>% 
     select(-IngallName1, -IngallName2) %>% tibble %>% distinct
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")  %>% tibble 

################################################################################
# Clean sample Meta data for 3 cruises
################################################################################

fetchCSV <- function(file, dataset){
     x <- read.csv(file) %>% tibble %>% filter(collectionSet == dataset) %>% 
          filter(str_detect(sampleID, "Smp")) %>% 
          select(sampleID, collectionSet, longitude, latitude, depth_m) %>% 
          rename(lat = latitude, lon = longitude, depth = depth_m) %>% 
          mutate(depth = case_when(depth == "UW" ~ 0,
                                   T ~ as.numeric(depth))) %>% 
          mutate(lat = as.numeric(lat), lon = as.numeric(lon*-1), depth = as.numeric(depth)) %>% 
          mutate(depth_cat = case_when(depth <= 15 ~ "0-15m",
                                       depth > 15 & depth <= 100 ~ "16-100m",
                                       depth > 100 & depth <= 200 ~ "101-200m",
                                       depth > 200 ~ ">200m")
          )
     return(x)
}
smpMet_gom <- fetchCSV(file = "data/1_metaData/sample_info/ATL_GOM_sampleMeta.csv", dataset = "GOM")
smpMet_atl <- fetchCSV(file = "data/1_metaData/sample_info/ATL_GOM_sampleMeta.csv", dataset = "ATL") 

smpMet_npo <- read.csv("data/1_metaData/sample_info/NPO_G2_sampleMeta.csv")  %>% tibble %>% 
     filter(str_detect(sampleID, "Smp")) %>% 
     select(sampleID, latitude, depth_m) %>% 
     rename(lat = latitude, depth = depth_m) %>% 
     mutate(depth = case_when(depth == "Underway" ~ 0, T ~ as.numeric(depth)),
            lat = as.numeric(lat), depth = as.numeric(depth)) %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "0-15m",
                                  depth > 15 & depth <= 100 ~ "16-100m",
                                  depth > 100 & depth <= 200 ~ "101-200m",
                                  depth > 200 ~ ">200m")
     )

#Save
#dir = "data/1_metaData/sample_info/"
#write.csv(smpMet_gom, paste0(dir, "GOM_Y2021_metaData_clean.csv"), row.names = F)
#write.csv(smpMet_atl, paste0(dir, "ATL_Y2021_metaData_clean.csv"), row.names = F)
#write.csv(smpMet_npo, paste0(dir, "NPO_G2.Y2017_metaData_clean.csv"), row.names = F)

################################################################################
# Apply Best Signals to Subset - ATLGOM
################################################################################

# Get ATLGOM Concentration dataframes
fetchInputs <- function(file, sepType){
     read.csv(file) %>% tibble %>% 
          select(Compound_Name, sampleID, conc_nM_L) %>%
          mutate(method = sepType)    
}
path = "data/ATL_GOM/ATL_GOM_"
data <- list(RP = fetchInputs(file = paste0(path, "RP_Metabs.csv"), sepType = "RP"),
             HP = fetchInputs(file = paste0(path, "HP_Metabs.csv"), sepType = "HP"),
             HN = fetchInputs(file = paste0(path, "HN_Metabs.csv"), sepType = "HN"))

# Pull relevant cmpds
bests <- sigIndex %>% select(Compound_Name, Ingall_Choice, ATLGOM_Choice)
bsHP = bests %>% filter(Ingall_Choice == "HP") %>% distinct %>% pull(Compound_Name) 
bsHN = bests %>% filter(Ingall_Choice == "HN") %>% distinct %>% pull(Compound_Name)
bsRP = bests %>% filter(Ingall_Choice == "RP") %>% distinct %>% pull(Compound_Name)

# Apply to each mode and bind together
rp <- data$RP %>% filter(Compound_Name %in% bsRP)
hp <- data$HP %>% filter(Compound_Name %in% bsHP)
hn <- data$HN %>% filter(Compound_Name %in% bsHN)

atlgom <- rbind(rp, hp, hn) %>% select(-method)

atlgom %>% select(Compound_Name) %>% distinct %>% 
     group_by(Compound_Name) %>% 
     reframe(count = n()) %>% 
     filter(count > 1)

# Correction for Homarine Stock Typo in Metabolite Library
atlgom  <- atlgom %>% tibble %>% 
     mutate(conc_nM_L = ifelse(Compound_Name == "Homarine", conc_nM_L * 0.04, conc_nM_L))

#-------------------------------------------------------------------------------
# Breaking ATLGOM to ATL and GOM
#-------------------------------------------------------------------------------

# ALT dataset Final
smps <- unique(smpMet_atl$sampleID)
atl <- atlgom %>% filter(sampleID %in% smps) 
write.csv(atl, "data/ATL_Metabs_Conc.csv", row.names = F)

# GOM dataset Final
smps <- unique(smpMet_gom$sampleID)
gom <- atlgom %>% filter(sampleID %in% smps) 
write.csv(gom, "data/GOM_Metabs_Conc.csv", row.names = F)

################################################################################
# Apply Best Signals to Subset - NPO
################################################################################

# Get NPO Concentration dataframes
fetchInputs <- function(file, sepType){
     read.csv(file) %>% tibble %>% 
          select(Compound_Name, sampleID, conc_nM_L) %>%
          mutate(method = sepType)    
}
path = "data/NPO_G2/NPO_G2_"
data <- list(RP = fetchInputs(file = paste0(path, "RP_Metabs.csv"), sepType = "RP"),
             HP = fetchInputs(file = paste0(path, "HP_Metabs_ProxyApplied.csv"), sepType = "HP"),
             HN = fetchInputs(file = paste0(path, "HN_Metabs_ProxyApplied.csv"), sepType = "HN"))

# Pull relevant cmpds
bests <- sigIndex %>% select(Compound_Name, Ingall_Choice, NPO_Choice)
bsHP = bests %>% filter(Ingall_Choice == "HP") %>% distinct %>% pull(Compound_Name) 
bsHN = bests %>% filter(Ingall_Choice == "HN") %>% distinct %>% pull(Compound_Name)
bsRP = bests %>% filter(Ingall_Choice == "RP") %>% distinct %>% pull(Compound_Name)

# Apply to each mode and bind together
rp <- data$RP %>% filter(Compound_Name %in% bsRP)
hp <- data$HP %>% filter(Compound_Name %in% bsHP)
hn <- data$HN %>% filter(Compound_Name %in% bsHN)

npog2 <- rbind(rp, hp, hn) %>% select(-method)

npog2 %>% select(Compound_Name) %>% distinct %>% 
     group_by(Compound_Name) %>% 
     reframe(count = n()) %>% 
     filter(count > 1)
length(unique(npog2$Compound_Name))
sort(unique(npog2$Compound_Name))
sort(sort(unique(data$HN$Compound_Name)))
sort(sort(unique(data$HP$Compound_Name)))
sort(sort(unique(data$RP$Compound_Name)))
#-------------------------------------------------------------------------------
# Add on Heal mSystem Dataset
#-------------------------------------------------------------------------------

# *** Appending Heal mSystem dataframe
# Any metabs shared, ours will be used instead of Heal

heal <- read.csv("data/NPO_G2/Heal_msystems_Conc2.csv") %>% tibble %>% 
     select(Compound_Name, Sample.ID, nM_Concentration) %>% 
     rename(name = Compound_Name, sampleID = Sample.ID, conc_nM_L = nM_Concentration)

names <- nameKey %>% filter(name %in% unique(heal$name)) %>% distinct

heal <- heal %>% left_join(., names, by = "name") %>% 
     rename(Compound_Name = Durham_Name) %>% 
     drop_na(Compound_Name) %>% 
     select(-name) %>% 
     mutate(sampleID = paste0("Smp_", sampleID))

# What Compounds are in Heal and not in NPOG2? 
list1 <- unique(npog2$Compound_Name)
list2 <- unique(heal$Compound_Name)
cmpds <- setdiff(list2, list1)

# Subset Heal to those not found in NPOG2. Then Combine dataframes
heal_subset <- heal %>% filter(Compound_Name %in% cmpds)
npog2_heal <- rbind(npog2, heal_subset) %>% distinct

# NPO dataset Final
smps <- unique(smpMet_npo$sampleID)
npo <- npog2_heal %>% filter(sampleID %in% smps)
write.csv(npo, "data/NPO_Metabs_Conc.csv", row.names = F)
