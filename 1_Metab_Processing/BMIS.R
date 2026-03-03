lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD SCRIPT PACKAGES
library("ggplot2")
library("tidyr")
library("cowplot")
library("readr")
library("stringr")
library("dplyr")

# Important Acronyms in this script
# NPO = North Pacific Ocean
# G2 = Gradient Cruise 2 survey taken in the year 2017
# ATL = Atlantic Ocean
# GOM = Gulf of Mexico
# HP = Hilic Positive
# HN = Hilic Negative
# RP = Reverse Phase Aqueous
# QEQC = Post-Quality Control
# BMIS = Best Matched Internal Standards

# Load sample data and metadata
waterRegion = "ATL" 
type = "GOM"
method = "RP"
method2 = "RP"
file_QCd = "data/ATL_GOM/QCd/QEQC_GOM_ATL_RP_expClean3.csv"

# Function to validate filenames
check_filename <- function(filename, keywords) {
     missing_keywords <- keywords[!sapply(keywords, grepl, filename)]
     if (length(missing_keywords) > 0) {
          cat("Warning: The file", filename, "is missing the following keywords:", paste(missing_keywords, collapse = ", "), "\n")
     }
}
# It will load in if you wrote mode and method correctly
check_filename(filename = file_QCd, keywords = c(waterRegion, type, method2))
preBMIS <- read.csv(file_QCd) 

#########################################################################################################
#########################################################################################################
dataType <- paste0(waterRegion, "_", type); dataType
data_path <- paste0("data/", waterRegion, "_", type, "/QCd/"); data_path
result_path <- paste0("data/", dataType, "/bmis/"); result_path
sample.meta <-  read_csv(paste0("data/1_metaData/sample_info/", dataType, "_sampleMeta.csv"))
metabolites <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")

# clean up
nameKey <- nameKey %>% filter(Durham_Name %in% unique(metabolites$Compound_Name)) %>% distinct

# Set output file paths
iSTD_bar_result <- paste0(result_path, "1a_", dataType, "_", method, "_iSTD_bar.svg")
iSTD_dot_result <- paste0(result_path, "1b_", dataType, "_", method, "_iSTD_dot.svg")
rsd_result <- paste0(result_path, "2a_", dataType, "_", method, "_rsdINJ.svg")
norm_rsd_data <- paste0(result_path, "2b_", dataType, "_", method, "_norm_rsd.csv")
new_normalized_data <- paste0(result_path, "BMIS_", dataType, "_", method, "_normMetabs.csv")

#########################################################################################################
cat("PART 1: Tidying and Checks of the QCd Skyline Output")
#########################################################################################################

# Adjust dataframe structure and retrieve relevant compound names from the nameKey file
# This links the 'Skyline_Name' to the 'Durham_Name' for uniformity across datasets

# Apply 'cleanQC' to the sample data and merge with the nameKey
names = nameKey %>% select(name, Durham_Name) %>% distinct
cleanQC <- function(data = NULL){
     preData <- data
     colnames(preData) <- as.character(preData[1, ])
     preData <- preData[-1, ] 
     colnames(preData)[colnames(preData) == "Mass.Feature"] <- "Compound_Name"

     pqc <- preData %>% 
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
     return(pqc)
}
data_QCd <- preBMIS %>% cleanQC() %>% 
     merge(., names, by.x = "Compound_Name", by.y = "name", all.x = TRUE) %>% 
     relocate("Durham_Name", .before = everything()) %>% 
     rename(oldName = Compound_Name, 
            Compound_Name = Durham_Name) %>% 
     tibble %>% print()
data_QCd %>% filter(is.na(Compound_Name)) %>% select(Compound_Name, oldName) %>% distinct

# Quick overview of the dataset:
cat("\nNumber of unique Smps in Data Input:\n"); length(unique(data_QCd$sampleID))
cat("\nNumber of unique compounds in Data Input:\n"); length(unique(data_QCd$Compound_Name))
unique(data_QCd$sampleID)

#--------------------------------------------------------
# Create the final dataframe structure by removing unwanted columns and preparing for analysis
#--------------------------------------------------------

sample.data = data_QCd %>% select(-oldName, -areaRAW)
unique(data_QCd$sampleGroup)
m = metabolites

# Merge data and metadata inputs
sample.info <- sample.meta %>% 
     select(sampleID, type, sampleGroup, replicate) %>% 
     mutate(Inj_vol = as.numeric(1)) %>%  
     mutate(Inj_vol = ifelse(grepl('Half', sampleID), 0.5, Inj_vol)) %>% 
     select(-sampleGroup)
unique(sample.info$sampleID)

# Merge the sample info with cleaned data...
 data <- sample.data %>%
     merge(., sample.info, by = "sampleID", all.x = TRUE) %>%
     filter(!type == "Std") %>% 
     mutate(Area = ifelse(!is.na(Inj_vol) & Inj_vol != 0.5, areaQCd, areaQCd * 2)) %>%
     rename(oldArea = areaQCd) %>%
     select(sampleID, type, sampleGroup, replicate, Inj_vol, 
            Compound_Name, oldArea, Area) %>%
     tibble
 unique(data$sampleGroup)
 
  # Filter the metabolite library based on the method (HN, HP, or RP)
filter_metabolites <- function(data, method) {
     if (method == "HN") {
          filtered_data <- data %>% 
               filter(Column == "HILIC" & ionMode == "Neg")
     } else if (method == "HP") {
          filtered_data <- data %>% 
               filter(Column == "HILIC" & ionMode == "Pos")
     } else if (method == "RP") {
          filtered_data <- data %>% 
               filter(Column == "RP" & ionMode == "Pos")
     } else {
          stop("Invalid method. Choose 'HN', 'HP', or 'RP'.")
     }
     return(filtered_data)
}
m <- filter_metabolites(m, method) %>% 
     drop_na(Concentration_uM) %>% 
     select(Compound_Name, Group, Concentration_uM, Column, ionMode)
m %>% select(Column, ionMode) %>% distinct
metabolites = m

#########################################################################################################
cat("PART 2: DETERMINE 'GOOD' INTERNAL STANDARDS TO USE")
#########################################################################################################

# Generate the list of internal standards (IS) from the metabolite library
IS.list <- metabolites %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)
# Exclude problematic IS
IS.list <- IS.list[!IS.list %in% 
                        c("Guanosine Monophosphate, 15N5 (GMP, 15N5)",
                          "Adenosine Monophosphate, 15N5 (AMP, 15N5)")]

cat("\nInternal Standards from your method choice, contained in the library:\n"); sort(IS.list)

# Filter sample data for those containing the selected internal standards
IS.dat <- data %>% filter(Compound_Name %in% IS.list) 
cat("\nInternal Standards from your method choice, contained in the data file:\n")
sort(unique(IS.dat$Compound_Name))

# Update IS list based on presence in the data and display the filtered list
IS.list <- IS.dat %>% pull(Compound_Name) %>% unique
cat("\nInternal Standards that will be used in BMIS script:\n"); IS.list

# Create a bar plot to visualize the areas of internal standards
bar <- ggplot(IS.dat, aes(x = sampleID, y = Area, fill = Compound_Name)) + 
     geom_bar(stat = "identity") + 
     facet_wrap( ~Compound_Name, scales = "free_y") +
     ggtitle("Internal Standards: Quality Controlled Areas") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5),
           axis.text.y = element_text(size = 7),
           legend.position = "bottom",
           strip.text = element_text(size = 7))
bar
svg(iSTD_bar_result, width = 15, height = 8); bar; dev.off()

# Function to identify 'bad' compounds based on missing or zero peak areas in over 50% of samples
identify_bad_compounds <- function(data, threshold = 0.50) {
     total_samples <- n_distinct(data$sampleID)
     bad_compounds <- data %>%
          group_by(Compound_Name) %>%
          mutate(Area = na_if(Area, 0)) %>% 
          summarise(
               num_NA_samples = sum(is.na(Area)),
               num_zero_samples = mean(Area == 0, na.rm = TRUE)) %>% 
          filter(num_NA_samples > total_samples * threshold | is.nan(num_zero_samples)) %>% 
          pull(Compound_Name)
     return(bad_compounds)
}
badIS <- identify_bad_compounds(data = bar$data)
cat("\nInternal Standards that will be dropped due to having little or no peak area:\n"); badIS

# Filter out the 'bad' internal standards from the dataset and update the IS list
IS.dat <- IS.dat %>% filter(!Compound_Name %in% badIS) 
IS.list <- IS.list[!unlist(IS.list) %in% unlist(badIS)]
cat("\nList of good standards that will be used:\n"); print(sort(IS.list))

# Output note for the user: Ensure that both pools (Poo) and transect samples (Smp) are available
cat(paste("NOTE!! You need pools (Poo), and transect samples (Smp) to continue this script:\n",
          " Sample Types that are present:\n")); unique(data$type)

#########################################################################################################
cat("PART 3: NORMALIZE AREAS USING INTERNAL STANDARDS")
#########################################################################################################

# Replace NA values in the 'Area' column with zero for downstream processing
data$Area[is.na(data$Area)] <- 0

# Generate a list of compounds that were detected across the dataset (non-blank samples)
Compds.detect <- data %>%      
     filter(type != "Blk") %>%
     select(sampleID, Compound_Name, Area) %>%
     group_by(Compound_Name) %>%
     summarise(percent.present = sum(!is.na(Area))/n()) %>%
     filter(percent.present > 0)
print(Compds.detect)

# Calculate the mean area values for each Internal Standard (IS)
IS.means <- IS.dat %>% 
     filter(type != "Blk") %>%
     mutate(Area = replace_na(Area, 0)) %>% 
     group_by(Compound_Name) %>%
     summarise(ave = mean(Area, rm.na = T)) %>% 
     rbind(., c("inj.vol", 1)) 
print(IS.means)

# Generate a dataframe that includes normalized areas for IS compounds
IS.dat_wArea <- data %>% 
     filter(!Compound_Name %in% badIS) %>%  
     select(sampleID, Compound_Name, Area, Inj_vol) %>% distinct %>% 
     group_by(sampleID, Compound_Name) %>% 
     reframe(Inj_vol = mean(Inj_vol, na.rm = T), 
             Area = mean(Area, na.rm = T)) %>% 
     filter(Compound_Name %in% sort(Compds.detect$Compound_Name)) %>% 
     pivot_wider(., names_from = "Compound_Name", values_from = "Area") %>% 
     mutate(inj.vol = Inj_vol) %>%
     select(-Inj_vol) %>% 
     relocate(inj.vol, .after = last_col()) %>% 
     as.data.frame()

# Create an extended list of internal standards, starting with 'inj.vol'
IS.listPlus <- c("inj.vol", IS.list)

# Test normalization with 'inj.vol'
this.IS <- IS.listPlus[1]
area.norm <- IS.dat_wArea[,-1] %>% 
     sapply(FUN = function(x) x / IS.dat_wArea[,grep(this.IS, names(IS.dat_wArea))]) %>%
     as.data.frame() %>% 
     mutate(sampleID = IS.dat_wArea$sampleID) %>%
     gather(Compound, Area_Norm, -sampleID) %>% tibble()
head(area.norm %>% tibble)

# Adjust normalized areas by the mean value of the current IS (injection volume)
this.mean <- IS.means %>% filter(Compound_Name == this.IS) %>% select(ave) %>% as.numeric
area.norm <- area.norm %>% mutate(Area_Norm = Area_Norm * this.mean)
key <- ncol(area.norm)
count <- length(which(!is.na(area.norm$Area_Norm))) / length(unique(area.norm$sampleID))
names(area.norm)[key] <- paste(this.IS, "Norm.Area", sep = ".")
print(area.norm)

# Loop through the remaining internal standards to normalize areas by each IS
for (i in 2:length(IS.listPlus)) {
     this.IS <- IS.listPlus[i]
     if(length(IS.dat_wArea[, grep(this.IS, names(IS.dat_wArea))]) != 0) {
          # Replace 0 with NA in the IS.dat_wArea column to avoid Inf results
          IS.dat_wArea_filtered <- IS.dat_wArea
          IS.dat_wArea_filtered[, grep(this.IS, names(IS.dat_wArea))] <- 
               ifelse(IS.dat_wArea_filtered[, grep(this.IS, names(IS.dat_wArea))] == 0, 
                      NA, 
                      IS.dat_wArea_filtered[, grep(this.IS, names(IS.dat_wArea))])
          
          this.norm <- IS.dat_wArea_filtered[,-1] %>%
               sapply(FUN = function(x) x / IS.dat_wArea_filtered[, grep(this.IS, names(IS.dat_wArea_filtered))]) %>%
               as.data.frame() %>%
               mutate(sampleID = IS.dat_wArea$sampleID) %>%
               gather(Compound, Area_Norm, -sampleID)
          # Handle NA values generated by division
          this.norm <- this.norm %>% mutate(Area_Norm = ifelse(is.infinite(Area_Norm) | is.nan(Area_Norm), NA, Area_Norm))
          
          this.mean <- IS.means %>% filter(Compound_Name == this.IS) %>% select(ave) %>% as.numeric
          this.norm <- this.norm %>% mutate(Area_Norm = Area_Norm * this.mean)
          
          key <- ncol(area.norm)
          area.norm[, key + 1] <- this.norm$Area_Norm
          names(area.norm)[key + 1] <- paste(this.IS, "Norm.Area", sep = ".")
          
          count <- length(which(!is.na(this.norm$Area_Norm))) / 
               length(unique(this.norm$sampleID))
     }
}
area.norm %>% tibble() %>% print() 
length(unique(area.norm$Compound))

# Merge the newly normalized data with the original dataset
mydata_new <- area.norm %>% mutate(Run.Cmpd = paste(sampleID, Compound))
mydata <- data %>% mutate(Run.Cmpd = paste(sampleID, Compound_Name))
dat <- full_join(mydata, mydata_new)
dat <- dat %>% select(-Compound, -oldArea) %>% arrange(sampleID, Compound_Name)

#########################################################################################################
#########################################################################################################
cat("PART 4: COMPARING RELATIVE STANDARD DEVIATIONS")
#########################################################################################################
#########################################################################################################

# Calculate the relative standard deviation (RSD) for each compound by sample type
dat2 <- dat %>%
     filter(Compound_Name %in% IS.means$Compound_Name) %>%
     select(-(sampleID:Area)) %>%
     gather(key = "MIS", value = "Adjusted_Area", factor_key = TRUE, -Run.Cmpd) %>% distinct %>% 
     left_join(dat %>% select(type, Compound_Name, Run.Cmpd) %>% distinct, 
               by = join_by(Run.Cmpd)) %>%
     mutate(Adjusted_Area = as.numeric(Adjusted_Area))

# Calculate RSD for sample type "Smp" (Samples)
smpdat <- dat2 %>%
     filter(type == "Smp") %>%
     dplyr::group_by(Compound_Name, MIS) %>%
     filter(!is.na(Adjusted_Area)) %>% 
     summarise(RSD_ofSmp = sd(Adjusted_Area) / mean(Adjusted_Area)) %>% 
     mutate(RSD_ofSmp = RSD_ofSmp) %>% 
     arrange(MIS)
smpdat
unique(mydata_new$sampleID)

# Calculate RSD for pool samples ("Poo")
alldat <- dat2 %>%
     filter(type == "Poo") %>%
     dplyr::group_by(Compound_Name, MIS) %>%
     filter(!is.na(Adjusted_Area)) %>% 
     summarise(RSD_ofPoo = sd(Adjusted_Area) / mean(Adjusted_Area)) %>% 
     mutate(RSD_ofPoo = RSD_ofPoo) %>% 
     left_join(smpdat) %>% 
     arrange(MIS)
alldat

# Extract data for injection volume normalization ("inj.vol.Norm.Area")
injectONlY <- alldat %>%
     filter(MIS == "inj.vol.Norm.Area" ) %>%
     rename(Orig_RSD_ofPoo = RSD_ofPoo, 
            Orig_RSD_ofSmp = RSD_ofSmp) %>%
     select(-MIS)

injectONlY_toPlot <- alldat %>% filter(MIS == "inj.vol.Norm.Area" ) 

# Set a threshold (cutoff) to identify significant changes in RSD
cut.off <- 0.4

# Calculate the difference between original RSD and new RSD, and determine acceptance based on cutoff
newalldat <- left_join(alldat, injectONlY, by = "Compound_Name") %>%
     mutate(del_RSD = (Orig_RSD_ofPoo - RSD_ofPoo), 
            percentDiff = del_RSD / Orig_RSD_ofPoo) %>%
     mutate(accept_MIS = (percentDiff > cut.off))

# Plot the RSD values, highlighting which compounds meet the cutoff for significant difference
plot <- ggplot() +
     geom_point(dat = newalldat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp, fill = accept_MIS))+ 
     scale_fill_manual(values=c("white","dark gray")) +
     geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
     facet_wrap(~ Compound_Name) +
     theme_bw()
 print(plot)
svg(rsd_result, width = 10, height = 8); plot; dev.off()

# Exclude blank samples and retain compounds that were detected
no.blank.dat <- dat %>% 
     filter(type != "Blk") %>%
     filter(Compound_Name %in% Compds.detect$Compound_Name) %>% 
     select(Compound_Name, type, sampleID, 'inj.vol.Norm.Area':ncol(.)) 
str(no.blank.dat)

# Gather the data for RSD calculations across all normalization methods (Normer)
rsd.stats <- no.blank.dat %>% 
     gather(Normer, Value, -Compound_Name, -type, -sampleID) %>% 
     mutate(sampleID = str_sub(sampleID, start = 1, end = -3)) %>% 
     mutate(sampleID = gsub("_Half", "", as.character(.$sampleID))) %>% 
     mutate(sampleID = gsub("_Full", "", as.character(.$sampleID))) %>%
     mutate(sampleID = gsub("_QC1to1", "", as.character(.$sampleID))) %>%
     mutate(sampleID = gsub("_dda", "", as.character(.$sampleID))) %>%
     group_by(Compound_Name, type, sampleID, Normer) %>%
     summarise(m = mean(Value, na.rm = T), 
               sd = sd(Value, na.rm = T), 
               rsd = sd / m) 
print(rsd.stats)

# Filter and clean the RSD statistics, keeping only meaningful values
rsd.clean <- rsd.stats %>% 
     filter(!is.na(m)) %>% 
     filter(Normer != "Area") %>% distinct()
unique(rsd.clean$Compound_Name)

# Calculate the variability in raw RSD values for "Poo" sample type using injection volume normalization
RawData.Variability <- rsd.clean %>%
     filter(Normer == "inj.vol.Norm.Area",
            type == "Poo",
            !is.na(rsd)) %>%
     rename(raw_rsd = rsd) %>% ungroup() %>%
     select(Compound_Name, sampleID, raw_rsd)
unique(RawData.Variability$Compound_Name)

# Set a low cutoff for filtering compounds with significant variability in raw RSD
low.cutoff <- 0.1

# Join RSD data with raw variability and calculate the difference in RSD for each compound
rsd.more <- full_join(rsd.clean, RawData.Variability, by = c("Compound_Name", "sampleID")) %>%
     mutate(del_poo = raw_rsd - rsd,
            percent_diff = del_poo / raw_rsd) %>%
     filter(raw_rsd > low.cutoff)
unique(rsd.more$Compound_Name)

#########################################################################################################
#########################################################################################################
cat("PART 5: GENERATING BMIS POOL MODELS")
#########################################################################################################
#########################################################################################################

# Set a cutoff for significant RSD differences
cutoff = 0.4

# Generate the BMIS model for pooled samples ("Poo")
PooModel <- rsd.more %>% filter(type == "Poo") %>%
     filter(percent_diff > cutoff) %>%
     select(-m, -sd, -del_poo, -percent_diff) %>%
     group_by(Compound_Name, Normer) %>%
     summarise(Mean.rsd = mean(rsd, na.rm = T)) %>%
     summarise(PooModelRSD = min(Mean.rsd),
               Poo.Picked.IS = unique(Normer)[which.min(Mean.rsd)][1])
print(PooModel)

# Combine the detected compounds and the BMIS model, filling missing IS with 'inj.vol.Norm.Area'
Models <- full_join(Compds.detect, PooModel, by = "Compound_Name") %>%
     select(-percent.present) %>%
     mutate(Poo.Picked.IS = ifelse(is.na(Poo.Picked.IS), "inj.vol.Norm.Area", Poo.Picked.IS))

# Merge RSD statistics with the BMIS model to get a complete dataset for RSD calculations
rsd.total <- full_join(rsd.stats, Models) %>% filter(!is.na(Compound_Name))

# Filter the RSD data for compounds with RSD greater than the low cutoff (0.1)
low.cutoff <- 0.1
rsd.more <- full_join(rsd.clean, RawData.Variability, by = c("Compound_Name", "sampleID")) %>%
     mutate(del_poo = raw_rsd - rsd, 
            percent_diff = del_poo / raw_rsd) %>%
     filter(raw_rsd > low.cutoff)

# Long and short formats for internal standard (IS) options based on RSD
AMIS.long <- rsd.more %>%
     ungroup() %>%
     filter(percent_diff > cutoff,
            type =="Poo") %>%
     select(Compound_Name, Normer, rsd) %>%
     arrange(Compound_Name , rsd) 
AMIS.short <- AMIS.long%>%
     group_by(Compound_Name) %>%
     summarise(IS.Options = paste(Normer, collapse=";"))
print(AMIS.short)

# Add the selected IS and model RSD to the complete RSD dataset
rsd.total <- rsd.total %>%
     mutate(PooPlus.IS = Poo.Picked.IS) %>%
     mutate(PooPlusModelRSD = PooModelRSD)

# Replace IS with the corresponding .Norm.Area where applicable
for (i in 1:nrow(rsd.total)){
     cmpd <- rsd.total$Compound_Name[i]
     if(length(grep(cmpd, IS.list))>0){
          newIS <- paste0(IS.list[grep(cmpd, IS.list)],".Norm.Area")
          rsd.total$PooPlus.IS[i] <- newIS
     }
}

# Remove unnecessary columns from the final RSD dataset, save it
rsd.total <- rsd.total %>% select(-PooModelRSD,-PooPlusModelRSD) 
write.csv(rsd.total, norm_rsd_data, row.names = F) 

# Prepare the final BMIS models dataset by selecting unique IS for each compound
models <- rsd.total %>% ungroup %>%
     select(Compound_Name, PooPlus.IS) %>%
     group_by(Compound_Name) %>%
     summarise(PooPlusModel.IS = unique(PooPlus.IS))

# Merge the final models with the original dataset
dat2 <- dat
dat.join <- as.data.frame(full_join(dat2, models)) %>%
     mutate(PooPlusModel = NA)

# Split the dataset by internal standard (IS) to calculate BMIS normalization
split.on.IS <- as.factor(dat.join$PooPlusModel.IS)
split.dat.join <- split(dat.join, split.on.IS)

# Apply BMIS normalization for each internal standard group
for (i in 1:length(split.dat.join)){
     col.key <-  which(names(split.dat.join[[i]]) == names(split.dat.join)[i])
     split.dat.join[[i]]$PooPlusModel <- split.dat.join[[i]][,col.key]
}

# Combine the split datasets and prepare the final BMIS-normalized dataset
unsplit.dat.join <- do.call(rbind, split.dat.join) %>% 
     rename(BMISNormalizedArea = PooPlusModel, BMIS.IS = PooPlusModel.IS) %>%
     tibble() %>% arrange(sampleID, Compound_Name) %>% distinct %>% 
     print() 
print(unsplit.dat.join)

# Save the BMIS-normalized data to a CSV file
df <- as.data.frame(unsplit.dat.join)
write.csv(df, new_normalized_data, row.names = F)

unique(df$sampleGroup)
unique(df$sampleID)
unique(df %>% filter(!Compound_Name %in% IS.list) %>% pull(Compound_Name))

