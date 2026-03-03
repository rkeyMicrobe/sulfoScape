lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)

# Function to validate filenames
check_filename <- function(filename, keywords) {
     missing_keywords <- keywords[!sapply(keywords, grepl, filename)]
     if (length(missing_keywords) > 0) {
          cat("Warning: The file", filename, "is missing the following keywords:", paste(missing_keywords, collapse = ", "), "\n")
     }
}

# PLACE YOUR INPUTS HERE....
method      = "HN" # Choices: HP, HN
Mode        = "Neg" # Choices: Pos, Neg
file_QCd = "data/ATL_GOM/QCd/QEQC_GOM_ATL_HILIC_Neg_exp8Clean.csv"
file_BMISd = "data/ATL_GOM/bmis/BMIS_ATL_GOM_HN_normMetabs.csv"

 # It will load in if you wrote mode and method correctly
check_filename(file_QCd, Mode)
preBMIS <- read.csv(file_QCd) 
check_filename(file_BMISd, method)
posBMIS <- read.csv(file_BMISd)

###################################################################################################
###################################################################################################
# CALLING GENERAL FILES. DOUBLE CHECK IF THESE ARE OKAY!
smp_MetaData <- read.csv("data/1_metaData/sample_info/ATL_GOM_sampleMeta.csv")
metab_library <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") 
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")
mixKey <- read.csv("data/1_metaData/metab_MixSets.csv") 

# CALLING GENERAL VARIABLES FOR THIS SCRIPT, DON'T TOUCH
region1     = "ATL" # Atlantic Ocean 
region2     = "GOM" # Gulf of Mexico
type        = "HILIC"
waterType = paste0(region1, "_", region2) 
dataType  = paste0(region1, "_", region2, "_", method); dataType 
data_path = paste0("data/", waterType, "/"); data_path
res_path  = paste0("results/RF_", waterType ,"/"); res_path
dataset_list = c(region1, region2)

###################################################################################################
# TIDYING UP FILE INPUTS
###################################################################################################

meta <- smp_MetaData %>% filter(collectionSet %in% dataset_list)
lib <- metab_library %>% filter(Column == type & ionMode == Mode)
istd <-  metab_library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name) %>% unique
nameKey <- nameKey %>% filter(Durham_Name %in% unique(lib$Compound_Name)) %>% distinct
mixKey <- mixKey %>% filter(ionMode == Mode)

data = preBMIS %>% tibble
# Tidy Pre-BMIS Input (these are Skyline integrations that were passed through Quality Control script)
cleanQC <- function(data = NULL, names = nameKey){
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
     cat("Data processing complete.")
     return(df)
}
preData <- preBMIS %>% cleanQC() %>% distinct %>% 
    # mutate(Compound_Name = str_squish(Compound_Name)) %>% 
     left_join(., nameKey, by = c("Compound_Name" = "name")) %>% 
     relocate("Durham_Name", .before = everything()) %>% tibble %>% 
     rename(oldName = Compound_Name,
            Compound_Name = Durham_Name) %>% 
     tibble %>% 
     mutate(Compound_Name = case_when(
          oldName == "(3-Carboxypropyl)trimethylammonium (TMAB)" ~ 
                     "(3-Carboxypropyl)trimethylammonium (TMAB)",
          oldName == "Acetyl-L-carnitine" ~ 
                     "Acetyl-L-Carnitine",
          T ~ Compound_Name
            ))
preData %>% filter(Compound_Name == "TMAB")
unique(sort(preData$Compound_Name))
preData %>% filter(is.na(Compound_Name)) %>% select(Compound_Name, oldName) %>% distinct

# Tidy Post-BMIS Input
posdata <- posBMIS %>%  tibble %>% 
     select(sampleID, Compound_Name, BMISNormalizedArea) %>% 
     rename(areaBMIS = BMISNormalizedArea) %>% 
     merge(., meta %>% tibble %>% select(sampleID, type, sampleGroup, replicate, full_vol_filtered_L), 
           by = "sampleID") %>% tibble
unique(sort(posdata$Compound_Name))

# What compounds are missing, Pre and Post BMIS?
pre_list <- preData %>% pull(Compound_Name) %>% unique
post_list <- posdata %>% pull(Compound_Name) %>% unique
cat("Compounds in Post-BMIS that are not in Pre-BMIS input: \n"); setdiff(post_list, pre_list)
cat("\nNumber of unique compounds in Data Input:\n"); length(post_list)


preData %>% pull(sampleID) %>% unique
posdata %>% pull(sampleID) %>% unique
posBMIS %>% pull(sampleID) %>% unique

###################################################################################################
# RF CALCULATIONS, PT 1
base::options(scipen = 999)
set.seed(1)
cat("Starting Part 1: Acquiring Compound RF Ratios in the Pre-BMIS Data input:")
 ###################################################################################################

# Select the 2021 mix set to pull concentrations and mix types
mixInfo <- mixKey %>% 
     select(Compound_Name, ionMode, conc_uM_2021, mixSet_2021) %>% 
     drop_na() %>% distinct %>% 
     rename(Concentration_uM = conc_uM_2021,
            HILIC_Mix = mixSet_2021)

# Acquiring means for each compound across sample types (this will average the replicates)
means <- preData %>% distinct %>% 
     left_join(., mixInfo, by = "Compound_Name") %>% 
     mutate(Area = areaQCd) %>% 
     group_by(sampleGroup, Compound_Name) %>%
     reframe(areaMean = mean(Area, na.rm = T)) %>% 
     drop_na(areaMean)

list <- sort(unique(means$Compound_Name))
cat("Compounds containing no Peak Areas across any samples: \n"); sort(setdiff(pre_list, list))

# Pull out only Standards from 'Means' object
stds <-  means %>% distinct %>%  
     filter(str_detect(sampleGroup, "Stds") |
                 str_detect(sampleGroup, "H2OinMatrix")) %>% 
     left_join(., mixInfo, by = "Compound_Name") %>% 
     select(Compound_Name, Compound_Name, sampleGroup, areaMean, Concentration_uM, HILIC_Mix) %>% tibble %>%
     filter(!Compound_Name %in% istd) %>% 
     distinct %>% tibble 
unique(sort(stds$Compound_Name))

list2 <- sort(unique(stds$Compound_Name))
cat("Compounds containing no Peak Areas across any 'Std' samples: \n"); sort(setdiff(list, list2))

# Prioritize Mix types
mix1_cmpds <- stds %>% filter(HILIC_Mix == "Mix1") %>% pull(Compound_Name) %>% unique() %>% sort
mix2_cmpds <- stds %>% filter(HILIC_Mix == "Mix2") %>% pull(Compound_Name) %>% unique() %>% sort

# Acquire areas for sample Group conditions
stds_h20 <- stds %>% 
     filter(
          (Compound_Name %in% mix1_cmpds & str_detect(sampleGroup, "Mix1")) |
               (Compound_Name %in% mix2_cmpds & str_detect(sampleGroup, "Mix2"))
     ) %>% 
     filter(str_detect(sampleGroup, "InH2O")) %>% 
     mutate(sampleGroup = "stdsH20") %>% 
     select(Compound_Name, Concentration_uM, areaMean) %>% 
     rename(stdsH20 = areaMean)
unique(sort(stds_h20$Compound_Name))


stds_mat <- stds %>% 
     filter(
          (Compound_Name %in% mix1_cmpds & str_detect(sampleGroup, "Mix1")) |
               (Compound_Name %in% mix2_cmpds & str_detect(sampleGroup, "Mix2")) 
     ) %>% 
     filter(str_detect(sampleGroup, "4uMStds")) %>% 
     filter(str_detect(sampleGroup, "InMatrix")) %>% 
     mutate(sampleGroup = "stdsMatrix") %>% 
     select(-HILIC_Mix, -sampleGroup, -Concentration_uM) %>% 
     rename(stdsMat = areaMean)
unique(sort(stds_mat$Compound_Name))

h20_mat <- stds %>% filter(sampleGroup == "H2OinMatrix") %>% 
     mutate(sampleGroup = "stdsMatrix") %>% 
     select(-HILIC_Mix, -sampleGroup, -Concentration_uM) %>% 
     rename(H2OMat = areaMean)
unique(sort(h20_mat$Compound_Name))

temp <- unique(c(stds_h20$Compound_Name, stds_mat$Compound_Name, h20_mat$Compound_Name))

rf_data <- 
     data.frame(Compound_Name = temp) %>% drop_na() %>% 
     left_join(., stds_h20, by = "Compound_Name") %>% 
     left_join(., stds_mat, by = "Compound_Name") %>% 
     left_join(., h20_mat, by = "Compound_Name") %>% tibble
unique(sort(rf_data$Compound_Name))

# Calculate Ratios
rf_ratios <- rf_data %>%
     mutate(MatrixInH2O = H2OMat  *10^-3,
            StdsInMatrix = stdsMat *10^-3,
            StdsinH2O = stdsH20 *10^-3,
            MatrixInH2O_adj = case_when(is.na(MatrixInH2O) ~ 0, T  ~ MatrixInH2O)
     ) %>% 
     mutate(RFs_mat = (StdsInMatrix - MatrixInH2O_adj) / Concentration_uM,
            RFs_h20 = StdsinH2O / Concentration_uM) %>% 
     select(-Concentration_uM, -stdsH20, -stdsMat, -H2OMat) %>% 
     
     mutate(RF_Ratio = RFs_mat / RFs_h20,
            RF_Ratio2 = case_when(RF_Ratio < 0 ~ 1, T  ~ RF_Ratio)
     )
unique(sort(rf_ratios$Compound_Name))

# Apply threshold, isolate the lower 1% of RF_Ratio_Final values, change them to 1
summary_rf_ratio <- rf_ratios %>%
     summarize(
          min_rf = min(RF_Ratio2, na.rm = TRUE),
          max_rf = max(RF_Ratio2, na.rm = TRUE),
          quantile_1 = quantile(RF_Ratio2, 0.01, na.rm = TRUE)  # 1st percentile
     )
threshold <- summary_rf_ratio$quantile_1  

final <- rf_ratios %>% mutate(
     RF_Ratio_Final = case_when(
          RF_Ratio2 <= threshold ~ 1,  
          TRUE ~ RF_Ratio2
     )
)
unique(sort(final$Compound_Name))

final %>% arrange(Compound_Name)
dir = paste0("data/" , waterType, "/ratios/")
write.csv(final, paste0(dir, "RF_", dataType, "_Ratios.csv"), row.names = F)

###################################################################################################
cat("Part 2: CALCULATE NMOL / L CONCENTRATIONS FOR ENVIRONMENTAL SAMPLES")
###################################################################################################

ratios <- final %>% 
     select(Compound_Name, RF_Ratio_Final, RFs_h20) %>% 
     mutate(injVol_uL =  400, tubeVol_mL = 1) 

# Append Ratio dataframe to post-BMIS dataframe containing adjusted peak areas
data <- posdata %>% distinct %>% 
     select(sampleID, sampleGroup, Compound_Name, areaBMIS, full_vol_filtered_L) %>% 
     rename(filtered_Full_L = full_vol_filtered_L) %>% 
     left_join(., ratios, by = "Compound_Name") %>% 
     filter(str_detect(sampleID,  "Smp_")) %>% distinct %>% 
     tibble %>% print 

final <- data %>% drop_na(RF_Ratio_Final) %>% 
     select(Compound_Name, sampleID, injVol_uL, filtered_Full_L, 
            RFs_h20, RF_Ratio_Final, areaBMIS
     ) %>% 
     mutate(areaBMIS = areaBMIS, 
            volume = (injVol_uL *10^-6) / filtered_Full_L
     ) %>% 
     mutate(calc1 = (areaBMIS) / (RFs_h20),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / RF_Ratio_Final)
     ) %>% 
     select(Compound_Name, sampleID, conc_nM_L)

write.csv(final, paste0(data_path, dataType, "_Metabs.csv"), row.names = F)
