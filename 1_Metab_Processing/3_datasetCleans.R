lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())


library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v5.csv") %>% tibble 
iSTDs <- library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)

################################################################################
# Handle ATL_GOM datasets
################################################################################

# Load in datasets:
fetchCSV <- function(file, sepType){
read.csv(file) %>% tibble %>% 
          filter(str_detect(sampleID, "Smp")) %>% 
          select(sampleID, Compound_Name, BMISNormalizedArea) %>% 
          distinct %>% 
          pivot_wider(., names_from = Compound_Name, values_from = BMISNormalizedArea)
}
dir = "data/ATL_GOM/bmis/BMIS_ATL_GOM_"
HN <- fetchCSV(file = paste0(dir, "HN_normMetabs.csv"), sepType = "HN")
HP <- fetchCSV(file = paste0(dir, "HP_normMetabs.csv"), sepType = "HP")
RP <- fetchCSV(file = paste0(dir, "RP_normMetabs.csv"), sepType = "RP")

HN <- fetchCSV(file = paste0(dir, "HN_normMetabs.csv"), sepType = "HN")
dat <- read.csv("data/ATL_GOM/ATL_GOM_Metabs_Fin.csv") %>% tibble()
met <- read.csv("data/1_metaData/sample_info/ATL_Y2021_metaData_clean.csv") %>% tibble()



# File Checks
CheckMeans <- function(data){
     data %>% 
          pivot_longer(., cols = -sampleID, names_to = "Compound_Name", values_to = "value") %>% 
          group_by(Compound_Name) %>% reframe(value = mean(value)) %>% 
          arrange(value)   
}
x <- CheckMeans(data = HN) %>% arrange(Compound_Name) %>% filter(is.na(value) | value <= 0) %>% print
x <- CheckMeans(data = HP) %>% arrange(Compound_Name) %>% filter(is.na(value) | value <= 0) %>% print
x <- CheckMeans(data = RP) %>% arrange(Compound_Name) %>% filter(is.na(value) | value <= 0) %>% print

# Find best signal across separation modes:
getSums <- function(data = NULL, sepType = NULL){
     x <- data %>% distinct
     sums <- colSums(x[, -which(names(x) == "sampleID")], na.rm = TRUE)
     x <- data.frame(metab = names(sums), total = sums, method = sepType)
     return(x)
}
x <- getSums(data = HN, sepType = "HN"); dim(x)
x <- getSums(data = HP, sepType = "HP"); dim(x)
x <- getSums(data = RP, sepType = "RP"); dim(x)


 bests <- rbind(getSums(data = HN, sepType = "HN"),
               getSums(data = HP, sepType = "HP"),
               getSums(data = RP, sepType = "RP")) %>% 
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
 
# Consider Overrides, based on Ingall Library Notes
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

atlGom_index <- bests_cleaned
sort(unique(atlGom_index$metab))

################################################################################
# Handle NPO data
################################################################################

dir = "data/NPO_G2/bmis/BMIS_NPO_G2_"
HN <- fetchCSV(file = paste0(dir, "HN_normMetabs.csv"), sepType = "HN")
HP <- fetchCSV(file = paste0(dir, "HP_normMetabs.csv"), sepType = "HP")
RP <- fetchCSV(file = paste0(dir, "RP_normMetabs.csv"), sepType = "RP")

bests <- rbind(getSums(data = HN, sepType = "HN"),
               getSums(data = HP, sepType = "HP"),
               getSums(data = RP, sepType = "RP")) %>% 
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
bests_cleaned <- applyOverrides(data = bests)

bests_cleaned
dir = "data/NPO_G2/BMIS_NPO_G2_"
write.csv(bests_cleaned, paste0(dir, "signal_index.csv"), row.names = F)

npoG2_index <- bests_cleaned
sort(unique(npoG2_index$metab))

################################################################################
# Check Comparisons
################################################################################

# Metabolites
# Shared Compounds Checks
list1 = npoG2_index %>% pull(metab) %>% unique
list2 = atlGom_index %>% pull(metab) %>% unique
# Lengths of each list
length1 <- length(list1)
length2 <- length(list2)

#Checks
shared_list12 <- intersect(list1, list2)

if (setequal(shared_list12, list1) && setequal(shared_list12, list2)) {
     cat("All compounds are shared between NPO and ATL/GOM. \n")
} else {
     cat("Not all compounds are shared between NPO and ATL/GOM.\n")
     
     # Compounds in list1 but not in list2
     non_shared_from_list1 <- setdiff(list1, list2)
     if (length(non_shared_from_list1) > 0) {
          cat("Compounds in NPO (list1) but not in ATL/GOM (list2):\n")
          cat(paste(non_shared_from_list1, collapse = "\n"), "\n")
     }
     
     # Compounds in list2 but not in list1
     non_shared_from_list2 <- setdiff(list2, list1)
     if (length(non_shared_from_list2) > 0) {
          cat("Compounds in ATL/GOM (list2) but not in NPO (list1):\n")
          cat(paste(non_shared_from_list2, collapse = "\n"), "\n")
     }
}
shared <- intersect(list1, list2); length(shared)
cat("Compounds are shared between NPO and ATL/GOM. \n")
sort(shared)

# Best Signal Compares
dim(npoG2_index)
dim(atlGom_index)

match_index <- npoG2_index %>% select(metab, bestMode) %>% 
     left_join(., atlGom_index %>% select(metab, bestMode), 
               by = "metab" ) %>% 
     rename(bestSig_NPO = bestMode.x, bestSig_ATLGOM = bestMode.y) %>% 
     mutate(NPO_ATLGOM = case_when(bestSig_NPO == bestSig_ATLGOM ~ "Good", T ~ "Bad")) %>% 
     distinct

match_index %>% filter(NPO_ATLGOM == "Bad")

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

dir = "data/1_metaData/sample_info/"
write.csv(smpMet_gom, paste0(dir, "GOM_Y2021_metaData_clean.csv"), row.names = F)
write.csv(smpMet_atl, paste0(dir, "ATL_Y2021_metaData_clean.csv"), row.names = F)
write.csv(smpMet_npo, paste0(dir, "NPO_G2.Y2017_metaData_clean.csv"), row.names = F)

################################################################################
# Make final dataframe for metabolite analysis - ATLGOM
################################################################################

bestSignal = read.csv("data/ATL_GOM/BMIS_ATL_GOM_signal_index.csv") %>% tibble %>% 
     mutate(best = paste0(bestMode, "_", metab))
bestSignal

 fetchInputs <- function(file, sepType){
     read.csv(file) %>% tibble %>% 
          select(Compound_Name, sampleID, conc_nM_L) %>% 
          mutate(method = sepType)    
 }
data <- rbind(fetchInputs(file = "data/ATL_GOM/ATL_GOM_RP_Metabs.csv", sepType = "RP"),
              fetchInputs(file = "data/ATL_GOM/ATL_GOM_HP_Metabs.csv", sepType = "HP"),
              fetchInputs(file = "data/ATL_GOM/ATL_GOM_HN_Metabs.csv", sepType = "HN"))  %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSignal$best)
sort(unique(data$Compound_Name))
data <- data %>% mutate(method_origin = substr(best, 1, 2)) %>% select(-best, -method)
atlgom_env <- data


write.csv(atlgom_env, "data/ATL_GOM/Metabs_Fin.csv", row.names = F)

# Correction for Homarine Stock Typo in Metabolite Library
final_clean  <- atlgom_env %>% tibble %>% 
     mutate(conc_nM_L = ifelse(Compound_Name == "Homarine", conc_nM_L * 0.04, conc_nM_L))

write.csv(atlgom_env, "data/ATL_GOM/Metabs_Fin_Clean.csv", row.names = F)



