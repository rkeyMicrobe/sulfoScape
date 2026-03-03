lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())


# LOAD DATA INPUTS
posBMIS <- read.csv("data/NPO_G2/bmis/BMIS_NPO_G2_HN_normMetabs.csv")
preBMIS <- read.csv("data/NPO_G2/QCd/QEQC_HILIC_Neg5.4_2017-2019_allSulfurStds.csv") 
method      = "HN" # Choices: HP or HN
Mode        = "Neg" # Pos or Neg


###################################################################################################
###################################################################################################
smp_MetaData <- read.csv("data/1_metaData/sample_info/NPO_G2_sampleMeta.csv")
metab_library <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v3.csv") 
std_mixSets <- read.csv("data/1_metaData/metabolites_MixSetsAngie.csv")
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")

# CALLING VARIABLES FOR THIS SCRIPT
region      = "NPO" # North Pacific Ocean 
cruise      = "G2" # Gradients 2 Cruise - 2017 survey
type        = "HILIC"
waterType = paste0(region, "_", cruise) 
dataType  = paste0(region, "_", cruise, "_", method) 
data_path = paste0("data/", waterType, "/")
res_path  = paste0("results/RF_", waterType ,"/")

###################################################################################################
# TIDY UP FILE INPUTS
###################################################################################################

# Metabolite information
lib <- metab_library %>% filter(Group != "Internal Standard" & Column == type & ionMode == Mode)
istd <-  metab_library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name) %>% unique
nameKey <- nameKey %>% filter(Durham_Name %in% unique(lib$Compound_Name)) %>% distinct
mixKey <- std_mixSets %>% filter(ionMode == Mode) %>% 
     select(Compound_Name, ionMode, Concentration_uM, HILIC_Mix) %>% 
     distinct 
cat("We are operating in: "); lib %>% select(Column, ionMode) %>% distinct


# Tidy Pre-BMIS Input (these are Skyline integrations that were passed through Quality Control script)
cleanQC <- function(data = NULL, names = nameKey){
     preData <- data
     colnames(preData) <- as.character(preData[1, ])
     preData <- preData[-1, ] 
     colnames(preData)[colnames(preData) == "Mass.Feature"] <- "Compound_Name"
     
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
preData <- cleanQC(data = preBMIS) %>% 
     rename(name = Compound_Name) %>% 
     filter(!name %in% istd) %>% 
     left_join(., nameKey, by = "name") %>% 
     relocate("Durham_Name", .before = everything()) %>% 
     rename(oldName = name,
            Compound_Name = Durham_Name) %>% 
     tibble
unique(preData$sampleID)
unique(preData$Compound_Name)
preData %>% filter(is.na(Compound_Name)) %>% select(oldName, Compound_Name) %>% distinct

# Tidy Post-BMIS Input
posdata <- posBMIS %>%  tibble %>% 
     filter(!Compound_Name %in% istd) %>% 
     select(sampleID, Compound_Name, BMISNormalizedArea) %>% 
     rename(areaBMIS = BMISNormalizedArea) %>% 
     merge(., smp_MetaData %>% tibble %>% select(sampleID, type, sampleGroup, replicate, full_vol_filtered_L), 
           by = "sampleID") %>% tibble

# What compounds are missing, Pre and Post BMIS?
pre_list <- preData %>% pull(Compound_Name) %>% unique
post_list <- posdata %>% pull(Compound_Name) %>% unique
cat("Compounds in Post-BMIS that are not in Pre-BMIS input: \n"); setdiff(post_list, pre_list)
cat("\nNumber of unique compounds in Data Input:\n"); length(post_list)

###################################################################################################
###################################################################################################
# ACQUIRING RFs TO CALCULATE RF RATIOS
base::options(scipen = 999)
set.seed(1)
cat("Starting Part 1: Acquiring Compound RF Ratios in the Pre-BMIS Data input:")
###################################################################################################
###################################################################################################

# Pull out only Standards from 'Means' object
stds <-  preData %>% 
     filter(Compound_Name %in% post_list) %>% 
     distinct %>%  
     left_join(., mixKey, by = "Compound_Name") %>% 
     mutate(Area = areaQCd) %>% 
     filter(str_detect(sampleGroup, "Stds")) %>% 
     mutate(Col.Key = ifelse(grepl("InH2O",sampleID),"StdsinH2O",
                             ifelse(grepl("H2OInM",sampleID),"MaxtrixInH2O",
                                    ifelse(grepl("H2OinM",sampleID),"MaxtrixInH2O","StdsInMatrix")))) %>% 
     
     mutate(Mix.Key = ifelse(grepl("Mix1", sampleID),"Mix1",
                             ifelse(grepl("Mix2", sampleID),"Mix2","NoMix"))) %>% 
     select(Compound_Name, Compound_Name, sampleGroup, Area, HILIC_Mix, runDate, Col.Key, Mix.Key) %>% 
     filter(!Compound_Name %in% istd) %>% 
     distinct %>% tibble

# Choosing 2018 std set run, this is not 2017 or 2019!
rf_data <- stds %>%
     filter(runDate == "180709") %>%
     group_by(Compound_Name, Col.Key, Mix.Key) %>%
     reframe(MeanArea = ifelse(all(is.na(Area)), NA, max(Area, na.rm = TRUE))) %>% 
     left_join(., mixKey, by = "Compound_Name") %>%
     mutate(key.match = Mix.Key == HILIC_Mix,
            Col.Key = ifelse(!key.match & Col.Key =="StdsInMatrix", "MatrixInH2O", Col.Key)) %>%
     filter(!(!key.match & Col.Key == "StdsinH2O")) %>%
     select(-key.match, -Mix.Key) %>%
     unique() %>% 
     mutate(MeanArea = MeanArea)

rf_ratios <- rf_data %>%
     pivot_wider(., names_from = "Col.Key", values_from = "MeanArea") %>%
     mutate(MatrixInH2O = MatrixInH2O*10^-3,
            StdsInMatrix = StdsInMatrix*10^-3,
            StdsinH2O = StdsinH2O*10^-3,
            MatrixInH2O_adj = case_when(is.na(MatrixInH2O) ~ 0, T  ~ MatrixInH2O)
     ) %>% 
     mutate(RFs_mat = (StdsInMatrix - MatrixInH2O_adj) / Concentration_uM,
            RFs_h20 = StdsinH2O / Concentration_uM) %>% 
     select(Compound_Name, RFs_mat, RFs_h20) %>% 
     
     mutate(RF_Ratio = RFs_mat / RFs_h20,
            RF_Ratio2 = case_when(RF_Ratio < 0 ~ 1, T  ~ RF_Ratio)
     )


heal <- read.csv(paste0("data/1_metaData/heal_ratios/RFs_Heal_NPO_", method, ".csv")) %>%  tibble %>% 
     rename(name =  Precursor.Ion.Name,
            Heal_Ratio = RF.ratio) %>% 
     merge(., nameKey, by = "name", all.x = T) %>% 
     rename(Compound_Name = Durham_Name) %>% 
     select(name, Compound_Name, Heal_Ratio) %>% distinct %>% 
     drop_na(Compound_Name) %>% 
     select(-name)

rf_ratios <- rf_ratios %>% left_join(., heal, by = "Compound_Name")
rf_ratios

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

# Save
dir = paste0("data/" , waterType, "/ratios/")
write.csv(final, paste0(dir, "/RF_", dataType, "_Ratios.csv"), row.names = F)

# Find out what is missing of the list of metabolites
x <- final %>%
     mutate(
          RF_Ratio_Final = ifelse(is.na(RF_Ratio_Final), 0, RF_Ratio_Final),
          Heal_Ratio = ifelse(is.na(Heal_Ratio), 0, Heal_Ratio)
     ) %>% 
     mutate(
          RF_Ratio_Final_rounded = round(RF_Ratio_Final, 3),
          Heal_Ratio_rounded = round(Heal_Ratio, 3)
     ) %>%
     filter(RF_Ratio_Final_rounded != Heal_Ratio_rounded)
write.csv(x, "data/RatioCalcs_HP_Sulfurs.csv", row.names = F)

# Extract situations where final RF Ratio is zero to catalog why
df <- rbind(final %>% filter(is.na(RF_Ratio_Final)),
            final %>% filter(RF_Ratio_Final < 0))
cmpds <- lib %>% filter(Compound_Name %in% 
                             c(lib %>% filter(is.na(Concentration_uM)) %>% pull(Compound_Name), 
                               df$Compound_Name)) %>% 
     select(Compound_Name, Concentration_uM)
count <-  metab_library %>% 
     filter(Compound_Name %in% df$Compound_Name) %>% 
     mutate(mode = paste0(Column, "_", ionMode)) %>% 
     group_by(Compound_Name, mode) %>% 
     summarise(presence = ifelse(n() > 0, 1, 0), .groups = 'drop') %>% 
     pivot_wider(names_from = mode, values_from = presence, values_fill = 0) %>%
     replace(is.na(.), 0)

ratios_NA <- df %>% 
     merge(., cmpds, by = "Compound_Name", all.y = T) %>% 
     merge(., count, by = "Compound_Name", all.x = T) %>% tibble

write.csv(ratios_NA, paste0(dir, "/RF_", dataType, "_RatiosMissing.csv"), row.names = F)


###################################################################################################
cat("Part 2: CALCULATE NMOL / L CONCENTRATIONS FOR ENVIRONMENTAL SAMPLES")
###################################################################################################

stockTubes <- mixKey %>% select(Compound_Name, Concentration_uM)

ratios <- final %>% 
     select(Compound_Name, RFs_h20, RF_Ratio_Final, Heal_Ratio
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
            RFs_h20, RF_Ratio_Final, Heal_Ratio, areaBMIS
     ) %>% 
     mutate(areaBMIS = areaBMIS, 
            filtered_Half_L = filtered_Full_L /2,
            volume = (injVol_uL *10^-6) / filtered_Full_L,
            dilution_factor = 4
     )) %>% 
     mutate(calc1 = (areaBMIS) / (RFs_h20),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / RF_Ratio_Final) * dilution_factor
     ) %>% 
     select(Compound_Name, sampleID, conc_nM_L)

write.csv(final, paste0(data_path, dataType, "_Metabs.csv"), row.names = F)
