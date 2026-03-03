lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())

# PLACE YOUR DATA INPUTS
posBMIS <- read.csv("data/NPO_G2/bmis/BMIS_NPO_G2_RP_normMetabs.csv")
preBMIS <- read.csv("data/NPO_G2/QCd/QEQC_RP_Aq5.1_2017-2019_allSulfurStds3.csv") 
smp_MetaData <- read.csv("data/1_metaData/sample_info/NPO_G2_sampleMeta.csv")
metab_library <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") 
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")

# CALLING GENERAL VARIABLES FOR THIS SCRIPT
region      = "NPO" # North Pacific Ocean 
cruise      = "G2" # Gradients 2 Cruise - 2017 survey
type        = "RP"
method      = "RP" 
Mode        = "Pos" 
waterType = paste0(region, "_", cruise); waterType
dataType  = paste0(region, "_", cruise, "_", method); dataType
data_path = paste0("data/", waterType, "/"); data_path
res_path  = paste0("results/RF_", waterType ,"/"); res_path

###################################################################################################
# TIDY UP FILE INPUTS
###################################################################################################

# Metabolite information
lib <- metab_library %>% filter(Group != "Internal Standard" & Column == type & ionMode == Mode)
istd <-  metab_library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name) %>% unique
nameKey <- nameKey %>% filter(Durham_Name %in% unique(lib$Compound_Name)) %>% distinct
cat("We are operating in: "); lib %>% select(Column, ionMode) %>% distinct

# Tidy Pre-BMIS Input (these are Skyline integrations that were passed through Quality Control script)
cleanQC <- function(data = NULL, names = nameKey){
     preData <- data
     colnames(preData) <- as.character(preData[1, ])
     preData <- preData[-1, ] 
     colnames(preData)[colnames(preData) == "Mass.Feature"] <- "Compound_Name"
     
     unique(preData$Replicate.Name) 
     
     df <- preData %>%
          select(Replicate.Name, Compound_Name, Area, Area.with.QC) %>%
          mutate(
               runDate = substr(Replicate.Name, 1, 6),  # New column for runDate
               sampleID = substr(Replicate.Name, 8, nchar(Replicate.Name)),
               Area = as.numeric(Area),
               Area.with.QC = as.numeric(Area.with.QC)
          ) %>%
          rename(
               Compound_Name = Compound_Name,
               areaRAW = Area,
               areaQCd = Area.with.QC
          ) %>%
          select(-Replicate.Name) %>%
          mutate(sampleGroup = str_extract(sampleID, "(?<=_)[^_]+(?=_)")) %>%
          tibble
     cat("Data processing complete.")
     return(df)
}
temp <- cleanQC(data = preBMIS)
preData <- temp %>% 
     rename(name = Compound_Name) %>% 
     filter(!name %in% istd) %>% 
     left_join(., nameKey, by = "name") %>% 
     relocate("Durham_Name", .before = everything()) %>% 
     rename(oldName = name,
            Compound_Name = Durham_Name) %>% 
     tibble 
unique(preData$runDate)
unique(preData$sampleID)
unique(preData$sampleGroup)
unique(preData$Compound_Name)
preData %>% filter(is.na(Compound_Name)) %>% select(oldName, Compound_Name) %>% distinct

 # Tidy Post-BMIS Input
posdata <- posBMIS %>%  tibble %>% 
     filter(!Compound_Name %in% istd) %>% 
     select(sampleID, Compound_Name, BMISNormalizedArea) %>% 
     rename(areaBMIS = BMISNormalizedArea) %>% 
     merge(., smp_MetaData %>% tibble %>% 
                select(sampleID, type, sampleGroup, replicate, full_vol_filtered_L), 
           by = "sampleID") %>% tibble
unique(posdata$Compound_Name)

# What compounds are missing, Pre and Post BMIS?
pre_list <- preData %>% pull(Compound_Name) %>% unique
post_list <- posdata %>% pull(Compound_Name) %>% unique
cat("Compounds in Post-BMIS that are not in Pre-BMIS input: \n"); setdiff(post_list, pre_list)
cat("\nNumber of unique compounds in Data Input:\n"); length(post_list)

###################################################################################################
# ACQUIRING RFs TO CALCULATE RF RATIOS
base::options(scipen = 999)
set.seed(1)
cat("Starting Part 1: Acquiring Compound RF Ratios in the Pre-BMIS Data input:")
###################################################################################################

unique(preData$runDate)
unique(preData$sampleGroup)
# Acquiring means for each compound across sample types (this will average the replicates)

means <- preData %>% distinct %>%
     filter(str_detect(sampleGroup, "Std") | str_detect(sampleGroup, "H2OinMatrix")) %>% 
     filter(runDate == "190725") %>% 
     mutate(Area = areaQCd) %>% 
     group_by(sampleGroup, Compound_Name) %>%
     reframe(MeanArea = ifelse(all(is.na(Area)), NA, max(Area, na.rm = TRUE)))
list <- sort(unique(means$Compound_Name))
cat("Compounds containing no Peak Areas across any samples: \n"); sort(setdiff(pre_list, list))
unique(means$sampleGroup)

temp <- lib %>% select(Compound_Name, Concentration_uM) %>% distinct %>% 
     filter(!Compound_Name == "Isobutyryl-L-Carnitine")

# Pull out only Standards from 'Means' object
stds <- means %>% distinct %>% 
     left_join(., temp, by = "Compound_Name") %>% 
     select(Compound_Name, sampleGroup, MeanArea, Concentration_uM) %>% tibble %>%
     filter(!Compound_Name %in% istd) %>% 
     distinct %>% tibble 
unique(stds$sampleGroup)


# Acquire areas for sample Group conditions
stds_h20 <- stds %>% 
     filter(sampleGroup == "250nMStdsInH2O") %>% 
     select(Compound_Name, Concentration_uM, MeanArea) %>% 
     rename(stdsH20 = MeanArea)
stds_mat <- stds %>% 
     filter(sampleGroup == "250nMStdsInMatrix") %>%
     select(-sampleGroup, -Concentration_uM) %>% 
     rename(stdsMat = MeanArea)

h20_mat <- stds %>% 
     filter(sampleGroup == "H2OinMatrix") %>% 
     select(-sampleGroup, -Concentration_uM) %>% 
     rename(H2OMat = MeanArea)

temp <- unique(c(stds_h20$Compound_Name, stds_mat$Compound_Name, h20_mat$Compound_Name))

rf_data <- 
     data.frame(Compound_Name = temp) %>% drop_na() %>% 
     left_join(., stds_h20, by = "Compound_Name") %>% 
     left_join(., stds_mat, by = "Compound_Name") %>% 
     left_join(., h20_mat, by = "Compound_Name") %>% tibble

# Calculate Ratios
rf_ratios <- rf_data %>%
     mutate(MatrixInH2O = H2OMat*10^-3,
            StdsInMatrix = stdsMat*10^-3,
            StdsinH2O = stdsH20*10^-3,
            MatrixInH2O_adj = case_when(is.na(MatrixInH2O) ~ 0, T  ~ MatrixInH2O)
     ) %>% 
     mutate(RFs_mat = (StdsInMatrix - MatrixInH2O_adj) / Concentration_uM,
            RFs_h20 = StdsinH2O / Concentration_uM) %>% 
     select(-Concentration_uM, -stdsH20, -stdsMat, -H2OMat) %>% 
     
     mutate(RF_Ratio = RFs_mat / RFs_h20,
            RF_Ratio2 = case_when(RF_Ratio < 0 ~ 1, T  ~ RF_Ratio)
     )

# Apply threshold, isolate the lower 1% of RF_Ratio_Final values, change them to 1 if under 0.1
summary_rf_ratio <- rf_ratios %>%
     summarize(
          min_rf = min(RF_Ratio2, na.rm = TRUE),
          max_rf = max(RF_Ratio2, na.rm = TRUE),
          quantile_1 = quantile(RF_Ratio2, 0.01, na.rm = TRUE)  # 1st percentile
     )
# Replace quantile_1 with 0.1 if it exceeds 0.1
if (summary_rf_ratio$quantile_1 > 0.1) {
     summary_rf_ratio$quantile_1 <- 0.1
}

threshold <- summary_rf_ratio$quantile_1  

final <- rf_ratios %>% mutate(
     RF_Ratio_Final = case_when(
          RF_Ratio2 <= threshold ~ 1,  
          TRUE ~ RF_Ratio2
     )
)

dir = paste0("data/" , waterType, "/ratios/")
write.csv(final, paste0(dir, "/RF_", dataType, "_Ratios.csv"), row.names = F)

###################################################################################################
cat("Part 2: CALCULATE NMOL / L CONCENTRATIONS FOR ENVIRONMENTAL SAMPLES")
###################################################################################################

stockTubes <- lib %>% select(Compound_Name, Concentration_uM) %>% 
     drop_na(Concentration_uM)
temp <- stockTubes %>% count(Compound_Name) %>% filter(n > 1) %>% pull(Compound_Name)
stockTubes <- stockTubes %>% filter(Compound_Name != temp)

temp <- final %>% count(Compound_Name) %>% filter(n > 1) %>% pull(Compound_Name)
ratios <- final %>%
     filter(!Compound_Name %in% temp) %>% 
     select(Compound_Name, RFs_h20, RF_Ratio_Final
     ) %>% 
     mutate(injVol_uL =  400, tubeVol_mL = 1) %>% 
     left_join(., stockTubes, by =  "Compound_Name")

# Append Ratio dataframe to post-BMIS dataframe containing adjusted peak areas
data <- posdata %>% distinct %>% 
     select(sampleID, sampleGroup, Compound_Name, areaBMIS, full_vol_filtered_L) %>% 
     rename(filtered_Full_L = full_vol_filtered_L) %>% 
     left_join(., ratios, by = "Compound_Name") %>% 
     filter(str_detect(sampleID,  "Smp_")) %>% distinct %>% 
     tibble %>% print 

final <- data %>% drop_na(RF_Ratio_Final) %>% 
     select(Compound_Name, sampleID, Concentration_uM, injVol_uL, filtered_Full_L, 
            RFs_h20, RF_Ratio_Final, areaBMIS
     ) %>% 
     mutate(areaBMIS = areaBMIS, 
            filtered_Half_L = filtered_Full_L /2,
            volume = (injVol_uL *10^-6) / filtered_Full_L,
            dilution_factor = 2
     ) %>% 
     mutate(calc1 = (areaBMIS) / (RFs_h20),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / RF_Ratio_Final) * dilution_factor
     ) %>% 
     select(Compound_Name, sampleID, conc_nM_L)

write.csv(final, paste0(data_path, dataType, "_Metabs.csv"), row.names = F)

###################################################################################################
cat("END OF SCRIPT")
###################################################################################################
