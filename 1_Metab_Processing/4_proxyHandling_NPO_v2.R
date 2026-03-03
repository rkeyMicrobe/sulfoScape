lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())

################################################################################
# Fetch Inputs
################################################################################

# NPO Post-BMIS datasets
fetchInputs <- function(dir, sepType){
     smpMet_npo <- read.csv("data/1_metaData/sample_info/NPO_G2_sampleMeta.csv")  %>% tibble %>% 
          filter(str_detect(sampleID, "Smp")) %>% 
          select(sampleID, full_vol_filtered_L)
     istd <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% 
          filter(Group == "Internal Standard") %>% pull(Compound_Name) %>% unique
     
     read.csv(paste0(dir, sepType, "_normMetabs.csv")) %>% tibble %>% 
          filter(str_detect(sampleID, "Smp")) %>% 
          filter(!Compound_Name %in% istd) %>% 
          left_join(., smpMet_npo, by = "sampleID") %>% 
          select(Compound_Name, sampleID, BMISNormalizedArea, full_vol_filtered_L) %>% 
          mutate(method = sepType) %>% 
          rename(filtered_vol = full_vol_filtered_L,
                 areaBMIS = BMISNormalizedArea)
}

npo_postBMIS <- rbind(fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "RP"),
                      fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HP"),
                      fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HN"))

atlgom_postBMIS <- rbind(fetchInputs(dir = "data/ATL_GOM/bmis/BMIS_ATL_GOM_", sepType = "RP"),
                         fetchInputs(dir = "data/ATL_GOM/bmis/BMIS_ATL_GOM_", sepType = "HP"),
                         fetchInputs(dir = "data/ATL_GOM/bmis/BMIS_ATL_GOM_", sepType = "HN"))
smps <- read.csv("data/1_metaData/sample_info/ATL_Y2021_metaData_clean.csv") %>% 
     pull(sampleID) %>% unique

atlgom_postBMIS <- atlgom_postBMIS %>% filter(sampleID %in% smps)

################################################################################
# Load in ATLGOM and NPO RF Ratios for all Compounds
################################################################################

# Get Ratio Info
fetchRatios <- function(dir1, dir2, method = NULL){
     x1 <- read.csv(paste0(dir1, method, "_Ratios.csv")) %>% tibble %>% 
          select(Compound_Name, RFs_h20, RF_Ratio_Final, Heal_Ratio) %>% 
          rename(G2_Ratio = RF_Ratio_Final)
     x2 <- read.csv(paste0(dir2, method, "_Ratios.csv")) %>% tibble %>% 
          select(Compound_Name, RFs_h20, RF_Ratio_Final) %>% 
          rename(Other_Ratio = RF_Ratio_Final,
                 Other_RFs_h20 = RFs_h20)
     
     x <- data.frame(Compound_Name = unique(sort(c(x1$Compound_Name, x2$Compound_Name)))) %>% 
          left_join(., x1, by = "Compound_Name") %>% 
          left_join(., x2, by = "Compound_Name") %>% tibble
     return(x)
}
RP_Ratios  <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                          dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "RP")
HP_Ratios <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                         dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "HP")
HN_Ratios <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                         dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "HN")

# Handle library concentration (also consider the NPO HILIC Mix Concentration)
mixKey <- read.csv("data/1_metaData/metabolites_MixSetsAngie.csv") %>% tibble %>% 
     select(Compound_Name, Concentration_uM, Column, ionMode) %>% 
     mutate(mode = paste0(Column, "_", ionMode)) %>% select(-Column, -ionMode) %>% 
     mutate(mode = case_when(mode == "HILIC_Pos" ~ "HP",
                             mode == "HILIC_Neg" ~ "HN",
                             T ~ mode)) %>% distinct 
lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v3.csv") %>% tibble %>% 
     select(Compound_Name, Column, ionMode, Concentration_uM) %>% 
     filter(Column == "RP") %>% mutate(mode = "RP") %>% 
     select(Compound_Name, Concentration_uM, mode)
concentrations <- rbind(mixKey, lib) %>% distinct

# Merge Ratio and Concentration Info (and any other pertinent variables)
process_ratios <- function(ratios, method, stocks) {
     remove <- ratios %>% count(Compound_Name) %>% filter(n > 1) %>% pull(Compound_Name)
     x <- stocks %>% filter(mode == method) 
     ratios %>% filter(!Compound_Name %in% remove) %>%
          left_join(x, by = "Compound_Name") %>% 
          # Injection info from LCM work
          mutate(injVol_uL =  400, tubeVol_mL = 1)
     
}

RP_Ratios <- process_ratios(ratios = RP_Ratios, method = "RP", stocks = concentrations)
HP_Ratios <- process_ratios(ratios = HP_Ratios, method = "HP", stocks = concentrations)
HN_Ratios <- process_ratios(ratios = HN_Ratios, method = "HN", stocks = concentrations)

################################################################################
# Clean the three inputs for the scan
################################################################################

# Load in ATLGOM Best Signal Index
signal_index <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, modeChoice) %>% 
     mutate(best = paste0(modeChoice, "_", Compound_Name)) 
bestSigs = signal_index %>% pull(best)

# Load in cmpds needing proxies
cmpds_help <- read.csv("data/proxyHandling_choices.csv") %>% 
     pull(Compound_Name) %>% sort()

# Clean Ratio dataframes, considering best sigs
ratios <- rbind(HP_Ratios, RP_Ratios, HN_Ratios) %>% 
     mutate(best = paste0(mode, "_", Compound_Name)) %>% 
     filter(best %in% bestSigs) %>% 
     mutate(mode = substr(best, 1, 2)) %>% select(-best) %>% distinct %>% 
     select(-Heal_Ratio)

# Clean NPO BMIS'd peak areas before scan
npo_candidates <- npo_postBMIS %>% filter(Compound_Name %in% cmpds_help)
unique(npo_candidates$Compound_Name)

# Clean the npo candidates to only consider the best signal chosen from ATLGOM
npo_canidates_clean <- npo_candidates %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSigs) %>%
     select(-best)

# Make summaries of all compounds in atlGOM to compare
atlgom_summaries <- atlgom_postBMIS %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSigs) %>% 
     group_by(Compound_Name) %>%
     reframe(
          mean = mean(areaBMIS, na.rm = TRUE),
          min = min(areaBMIS, na.rm = TRUE),
          max = max(areaBMIS, na.rm = TRUE)) %>% 
     filter(!mean <= 0)

# Make summaries of the compounds needing proxies...
npo_summaries <- npo_canidates_clean %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSigs) %>% 
     group_by(Compound_Name) %>%
     reframe(
          mean = mean(areaBMIS, na.rm = TRUE),
          min = min(areaBMIS, na.rm = TRUE),
          max = max(areaBMIS, na.rm = TRUE))  %>% 
     filter(!mean <= 0) %>% 
     select(-min, -max)

################################################################################
# SYSTEMATIC SCAN considering Absolute Differences
################################################################################

# Initialize an empty dataframe to store results
final_proxies <- data.frame()
tops = 75

# Loop through each compound in npo_summaries
for (i in 1:nrow(npo_summaries)) {
     # Get the current compound's mean and name
     candidate_mean <- npo_summaries$mean[i]
     candidate_name <- npo_summaries$Compound_Name[i]
     
     # Find the 3 closest proxies in atlgom_summaries
     proxies <- atlgom_summaries %>%
          mutate(candidate_mean = candidate_mean, 
                 diff_BMIS = abs(candidate_mean - mean)) %>%
          arrange(diff_BMIS) %>%
          slice_head(n = tops) %>%
          mutate(Candidate = candidate_name) %>%
          rename(potentialProxy = Compound_Name,
                 potentialProxy_mean = mean)
     
     # Append the results to the final dataframe
     final_proxies <- bind_rows(final_proxies, proxies)
}
scan_proxies <- final_proxies %>% tibble %>% 
     select(Candidate, potentialProxy, 
            candidate_mean, potentialProxy_mean, diff_BMIS)

# Append info on to BMIS_chosen and Subjective_Chosens potentials
temp <- read.csv("data/proxyHandling_choices.csv") %>% tibble
cmpds <- unique(c(temp$BMIS_RSD, temp$Chem1, temp$Chem2, temp$Chem3))
cmpds <- cmpds[!is.na(cmpds)]

 ################################################################################
# Perform NPO environment calculations given Proxy'd info
################################################################################

# Pull Scan findings and append choices from BMIS Relative Standard Dev and Chemical Structure Choices
cand_prox <- scan_proxies %>% 
     select(Candidate, potentialProxy) %>% 
     rename(BMIS_AD = potentialProxy) %>% 
     filter(!BMIS_AD %in% cmpds) %>% 
     left_join(., temp, by = c("Candidate" = "Compound_Name")) 

cmpds <- unique(c(cand_prox$BMIS_AD, cand_prox$BMIS_RSD, cand_prox$Chem1, cand_prox$Chem2, cand_prox$Chem3))

proxyInfo <- ratios %>% filter(Compound_Name %in% cmpds) %>% 
     select(-RFs_h20, -G2_Ratio) %>% arrange(Compound_Name)

# Combine Cand_prox and Proxy_Info dataframes to get all combinations
proxyInfo2 <- cand_prox %>% 
     pivot_longer(., cols = -Candidate, names_to = "Proxy_Type", values_to = "Proxy_Name") %>% 
     left_join(., proxyInfo, by = c("Proxy_Name" = "Compound_Name")) %>% 
     mutate(Proxy_Type = case_when(Proxy_Type == "BMIS_AD" ~ "AbsDiff", 
                                   Proxy_Type == "BMIS_RSD" ~ "RelStdDev",
                                   Proxy_Type == "Chem1" ~ "Human1",
                                   Proxy_Type == "Chem2" ~ "Human2",
                                   Proxy_Type == "Chem3" ~ "Human3"
                                   )) %>% drop_na(Other_Ratio) %>% 
     rename(RFs_h20 = Other_RFs_h20, Ratio = Other_Ratio) %>% 
     mutate(Combination = paste0(Candidate, " : ", Proxy_Name))

# Append this info onto transect data with BMIS areas, perform calculations
data <- npo_canidates_clean %>% 
     left_join(., proxyInfo2, by = c("Compound_Name" = "Candidate")) # warning expected given given  multiple proxies

results <- data %>%
     mutate(filtered_vol = filtered_vol / 2, 
            volume = (injVol_uL * 1e-6) / filtered_vol,
            dilution_factor = 4, 
            calc1 = areaBMIS / RFs_h20, 
            calc2 = calc1 * volume,
            conc_nM_L = calc2 * (1 / Ratio) * dilution_factor
     ) %>%
     select(Compound_Name, Proxy_Type, Proxy_Name, Combination, sampleID, conc_nM_L) %>%
     distinct() %>% 
     drop_na(Proxy_Type)


# Now, Pull the environment results from ATLGOM for the Candidates -- ATL samples only!
reference <- read.csv("data/ATL_GOM/ATL_GOM_Metabs_Fin.csv") %>% 
     filter(Compound_Name %in% cand_prox$Candidate) %>% 
     filter(sampleID %in% smps) %>% 
     mutate(Proxy_Type = "Reference",
            Proxy_Name = Compound_Name,
            Combination = paste0(Compound_Name, " : ", Proxy_Name)
            ) %>% 
     select(Compound_Name, Proxy_Type, Proxy_Name, Combination, sampleID, conc_nM_L)
     
temp <- rbind(results, reference) %>% ungroup
dir = "results/proxyHandling/"

category_Order = c("Reference", "RelStdDev", "AbsDiff", "Human1", "Human2", "Human3")
category_colors <- c(
     "Reference" = "#000000",  # Black
     "RelStdDev" = "#CC79A7",  # Yellow
     "AbsDiff" = "#009E73",    # Green
     "Human1" = "#E69F00",     # Orange
     "Human2" = "#E69F00",     # Blue
     "Human3" = "#E69F00"      # Purple
)

# Split the data by Compound_Name
plots <- split(temp, temp$Compound_Name)

#removes <- c("Eicosapentaenoic Acid (EPA)"
#unique(plots$Taurocyamine$Proxy_Name)
#plots$`Taurocyamine` <- plots$`Taurocyamine` %>% filter(!Proxy_Name %in% removes)

# Loop over each split dataframe and save the plot
lapply(names(plots), function(compound_name) {
     # Get the subset dataframe
     df <- plots[[compound_name]]
     
     referenceValues = df %>% 
          filter(Proxy_Type == "Reference") %>% 
          group_by(Compound_Name) %>% 
          reframe(
               mean = mean(conc_nM_L, na.rm = TRUE),
               min = min(conc_nM_L, na.rm = TRUE),
               max = max(conc_nM_L, na.rm = TRUE)) 
     
     # Create the plot
     p <- ggplot(df, aes(y = Proxy_Name, 
                         x = conc_nM_L, color = Proxy_Type)) +
          geom_point() +
          #geom_bar(stat = "identity", position = "dodge") +
          geom_vline(xintercept = referenceValues$min, color = "grey") +
          geom_vline(xintercept = referenceValues$mean, color = "black") +
          geom_vline(xintercept = referenceValues$max, color = "grey") +
          scale_color_manual(values = category_colors) +
          labs(
               x = "Mean Value",
               y = "ProxyType",
               fill = "Proxy_Type",
               title = paste("Mean Values by Combination for", compound_name),
               subtitle = "Bars Colored by Proxy Type"
          ) +
          theme(legend.position = "right") 
     
      # Save the plot
     file_name <- paste0("results/proxyHandling/", compound_name, ".png")
     ggsave(file_name, plot = p, width = 16, height = 15, dpi = 300, bg = "white")
})

################################################################################
# Perform NPO environment calculations given Proxy'd info
################################################################################
npo_postBMIS <- rbind(fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "RP"),
                      fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HP"),
                      fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HN"))

npo_postBMIS %>% select(Compound_Name) %>% distinct


# Load in ATLGOM Best Signal Index
signal_index <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, modeChoice) %>% 
     mutate(best = paste0(modeChoice, "_", Compound_Name))  
bestSigs = signal_index %>% pull(best)

npo_transect <- npo_postBMIS %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSigs)
unique(npo_transect$Compound_Name)


# Read in what was chosen based on plots above
picks <- read.csv("data/1_metaData/proxyPicks.csv") %>% tibble %>% drop_na() %>% 
     left_join(., proxyInfo, by = c("ProxyPick" = "Compound_Name"))

# Calculate environment concentrations with the picked proxies in NPO
npo_candidates <- npo_transect %>% 
     filter(Compound_Name %in% picks$Compound_Name) %>% 
     left_join(., picks, by = "Compound_Name") %>% 
     rename(waterRF = Other_RFs_h20, Ratio = Other_Ratio) %>% 
     mutate(dilution_factor = case_when(method %in% c("HP", "HN") ~ 4, method == "RP" ~ 2)) %>% 
     # Perform Calculations
     mutate(areaBMIS = areaBMIS, 
            filtered_Half_L = filtered_vol /2,
            volume = (injVol_uL *10^-6) / filtered_vol) %>% 
     mutate(calc1 = (areaBMIS) / (waterRF),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / Ratio) * dilution_factor) %>% 
     select(Compound_Name, sampleID, conc_nM_L)
unique(npo_candidates$Compound_Name)

npo_nonCandidates <- npo_transect %>% 
     filter(!Compound_Name %in% picks$Compound_Name) %>% 
     left_join(., ratios, by = "Compound_Name") %>% 
     rename(waterRF = Other_RFs_h20, Ratio = Other_Ratio) %>% 
     mutate(dilution_factor = case_when(method %in% c("HP", "HN") ~ 4, method == "RP" ~ 2)) %>% 
     # Perform Calculations
     mutate(areaBMIS = areaBMIS, 
            filtered_Half_L = filtered_vol /2,
            volume = (injVol_uL *10^-6) / filtered_vol) %>% 
     mutate(calc1 = (areaBMIS) / (waterRF),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / Ratio) * dilution_factor) %>% 
     select(Compound_Name, sampleID, conc_nM_L)
unique(npo_nonCandidates$Compound_Name)

combined <- rbind(npo_candidates, npo_nonCandidates) %>% arrange(Compound_Name)
unique(combined$Compound_Name)

cmpds_we_have <- unique(combined$Compound_Name)

# Add on Heal Information from mSystems
nameKey <- read_csv("data/1_metaData/nameKey_update.csv") %>% tibble()
msystems <- read_csv("data/Heal_msystems_Conc2.csv") %>% tibble() %>% 
     select(Compound_Name, `Sample ID`, nM_Concentration) %>% 
     rename(sampleID = `Sample ID`, conc_nM_L = nM_Concentration)

# Update msystems data with new names from nameKey
msystems_updated <- msystems %>% 
     left_join(nameKey, by = c("Compound_Name" = "name")) %>% 
     filter(!is.na(Durham_Name)) %>%      
     rename(oldName = Compound_Name, 
            Compound_Name = Durham_Name) %>% 
     select(-oldName)  %>% 
     mutate(sampleID = paste0("Smp_", sampleID))

# Ensure unique compounds are not duplicated in the final dataset
list1 <- unique(msystems_updated$Compound_Name)
list2 <- unique(combined$Compound_Name)  
msystems_updated <- msystems_updated %>% 
     filter(!Compound_Name %in% intersect(list1, list2)) %>% 
     select(Compound_Name, sampleID, conc_nM_L)
unique(msystems_updated$Compound_Name)

msystems_updated
combined
unique(combined$Compound_Name)


# Combine outs with Heal
combined2 <- rbind(combined, msystems_updated)

unique(combined2$sampleID)
directory = "data/NPO_G2/NPO_G2_"
write.csv(combined2, paste0(directory, "Metabs_Fin.csv"), row.names = F)


