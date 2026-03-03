lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")


library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble 
iSTDs <- library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)

bestModes <-library %>% select(Compound_Name, modeChoice) %>% drop_na()

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
#write.csv(smpMet_gom, paste0(dir, "GOM_Y2021_metaData_clean.csv"), row.names = F)
#write.csv(smpMet_atl, paste0(dir, "ATL_Y2021_metaData_clean.csv"), row.names = F)
#write.csv(smpMet_npo, paste0(dir, "NPO_G2.Y2017_metaData_clean.csv"), row.names = F)

################################################################################
# Make final dataframe for metabolite analysis - ATLGOM
################################################################################

bsHP = bestModes %>% filter(modeChoice == "HP") %>% distinct %>% pull(Compound_Name) 
bsHN = bestModes %>% filter(modeChoice == "HN") %>% distinct %>% pull(Compound_Name)
bsRP = bestModes %>% filter(modeChoice == "RP") %>% distinct %>% pull(Compound_Name)

fetchInputs <- function(file, list, sepType){
     read.csv(file) %>% tibble %>% 
          select(Compound_Name, sampleID, conc_nM_L) %>%
          filter(Compound_Name %in% list) %>% 
          mutate(method = sepType)    
}
path = "data/ATL_GOM/ATL_GOM_"
data <- rbind(fetchInputs(file = paste0(path, "RP_Metabs.csv"), list = bsRP, sepType = "RP"),
              fetchInputs(file = paste0(path, "HP_Metabs.csv"), list = bsHP, sepType = "HP"),
              fetchInputs(file = paste0(path, "HN_Metabs.csv"), list = bsHN, sepType = "HN"))
sort(unique(data$Compound_Name))
atlgom_env <- data

write.csv(atlgom_env, "data/ATL_GOM/Metabs_Fin2.csv", row.names = F)

# Correction for Homarine Stock Typo in Metabolite Library
final_clean  <- atlgom_env %>% tibble %>% 
     mutate(conc_nM_L = ifelse(Compound_Name == "Homarine", conc_nM_L * 0.04, conc_nM_L))

write.csv(final_clean, "data/ATL_GOM/Metabs_Fin_Clean2.csv", row.names = F)

################################################################################
# Check Alignment with Collaborator
################################################################################
baddies = c("Smp_S5UWT2-D-S16E39D5-E_A", "Smp_S15E35D5_D", "Smp_S2E13D11_D")

lib <- library %>% select(Compound_Name, Group) %>% distinct
smps <- smpMet_atl %>% 
     filter(depth <= 15) %>% 
     filter(!sampleID %in% baddies) %>% 
     pull(sampleID) %>% unique

plot_data <- final_clean %>% 
     filter(sampleID %in% smps,
     ) %>% 
     merge(., lib, by = "Compound_Name") %>%
     group_by(sampleID, Group) %>%
     reframe(concentration = sum(conc_nM_L, na.rm = TRUE)) 

ggplot(plot_data, aes(x = sampleID, y = concentration, fill = Group)) +
     geom_bar(stat = "identity") +
     ylab("concentration (nmol)") +
     xlab("") +
     ylim(0, 150) +
     theme_gray() +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

################################################################################
# Make final dataframe for metabolite analysis - NPO
################################################################################

bsHP = bestModes %>% filter(modeChoice == "HP") %>% distinct %>% pull(Compound_Name) 
bsHN = bestModes %>% filter(modeChoice == "HN") %>% distinct %>% pull(Compound_Name)
bsRP = bestModes %>% filter(modeChoice == "RP") %>% distinct %>% pull(Compound_Name)

# Load in environmental concentrations
path = "data/NPO_G2/NPO_G2_"
data <- rbind(fetchInputs(file = paste0(path, "RP_Metabs.csv"), list = bsRP, sepType = "RP"),
              fetchInputs(file = paste0(path, "HP_Metabs.csv"), list = bsHP, sepType = "HP"),
              fetchInputs(file = paste0(path, "HN_Metabs.csv"), list = bsHN, sepType = "HN")) 
sort(unique(data$Compound_Name))
npo_env <- data

# Append Heal msystem Concentrations
heal <- read.csv("data/NPO_G2/Heal_msystems_Conc2.csv") %>% tibble %>% 
     select(Compound_Name, Sample.ID, nM_Concentration) %>% 
     rename(sampleID = Sample.ID, name = Compound_Name)
cmpds <- unique(heal$name)
## Clean her names up
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")  %>% tibble %>% 
     filter(name %in% cmpds) %>% distinct
     
heal <- heal %>% 
     left_join(., nameKey, by = "name") %>% 
     rename(Compound_Name = Durham_Name) %>% 
     select(Compound_Name, sampleID, nM_Concentration)

# For those not contained in the heal dataset, we will take these from our data and append to hers
list <- unique(data$Compound_Name)
list2 <- unique(heal$Compound_Name)
cmpds <- setdiff(list, list2)

subdata <- data %>% 
     filter(Compound_Name %in% cmpds) %>% 
     select(-method) %>% 
     rename(nM_Concentration = conc_nM_L)

dataHeal <- rbind(subdata, heal) %>% drop_na(Compound_Name)

#------------------------------------------------------------------
# Handling Proxy situations... 
#------------------------------------------------------------------

# Load in Proxy Picks
picks <- read.csv("data/1_metaData/proxyPicks.csv") %>% tibble %>% drop_na() 
proxy <- unique(picks$ProxyPick)
cmpds <- unique(picks$Compound_Name)

# Fetch BMIS'd NPO info
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
data <- rbind(fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "RP"),
              fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HP"),
              fetchInputs(dir = "data/NPO_G2/bmis/BMIS_NPO_G2_", sepType = "HN"))  %>% 
     mutate(best = paste0(method, "_", Compound_Name)) %>% 
     filter(best %in% bestSignal$best)
sort(unique(data$Compound_Name))

# Fetch NPO Ratios
npo_ratios <- 
     rbind(read.csv("data/NPO_G2/ratios/RF_NPO_G2_HP_Ratios.csv") %>% mutate(mode = "HP"),
           read.csv("data/NPO_G2/ratios/RF_NPO_G2_HN_Ratios.csv") %>% mutate(mode = "HN"),
           read.csv("data/NPO_G2/ratios/RF_NPO_G2_RP_Ratios.csv") %>% mutate(mode = "RP") %>% 
     select(Compound_Name, RFs_mat, RFs_h20, RF_Ratio, RF_Ratio2, Heal_Ratio, RF_Ratio_Final, mode)
     ) %>% 
     
     
     filter(Compound_Name %in% cmpds) %>% 
     mutate(best = paste0(mode, "_", Compound_Name)) %>% 
     filter(best %in% bestSignal$best) %>% 
     mutate(mode = substr(best, 1, 2)) %>% select(-best) %>% distinct %>% 
     select(-Heal_Ratio) %>% 
     tibble

# ---------------------------------------------------------------------------

# Grab their ATLGOM info -- these will replace in NPO ratio dataframe
atlgom_ratios <- rbind(read.csv("data/ATL_GOM/ratios/RF_ATL_GOM_HP_Ratios.csv") %>% mutate(mode = "HP"),
                    read.csv("data/ATL_GOM/ratios/RF_ATL_GOM_HN_Ratios.csv") %>% mutate(mode = "HN"),
                    read.csv("data/ATL_GOM/ratios/RF_ATL_GOM_RP_Ratios.csv") %>% mutate(mode = "RP")
                    ) %>% 
     mutate(best = paste0(mode, "_", Compound_Name)) %>% 
     filter(best %in% bestSignal$best) %>% 
     mutate(mode = substr(best, 1, 2)) %>% select(-best) %>% distinct %>% 
     tibble %>% 
     filter(Compound_Name %in% cmpds) %>% 
     select(Compound_Name, RFs_mat, RFs_h20, RF_Ratio, RF_Ratio2, RF_Ratio_Final, mode)

# Replace with new content
proxy_values <- picks %>% left_join(atlgom_ratios, by = c("ProxyPick" = "Compound_Name"))

npo_ratios_updated <- npo_ratios %>%
     left_join(proxy_values, by = "Compound_Name") %>%
     mutate(
          RFs_mat   = coalesce(RFs_mat.y, RFs_mat.x),
          RFs_h20   = coalesce(RFs_h20.y, RFs_h20.x),
          RF_Ratio  = coalesce(RF_Ratio.y, RF_Ratio.x),
          RF_Ratio2 = coalesce(RF_Ratio2.y, RF_Ratio2.x),
          RF_Ratio_Final = coalesce(RF_Ratio_Final.y, RF_Ratio_Final.x),
          mode = coalesce(mode.y, mode.x)
     ) %>%
     select(Compound_Name, RFs_mat, RFs_h20, RF_Ratio, RF_Ratio2, RF_Ratio_Final, mode) %>% 
     arrange(Compound_Name) 
cmpds <- npo_ratios_updated$Compound_Name

# Handle library concentration (also consider the NPO HILIC Mix Concentration)
concMix <- read.csv("data/1_metaData/metab_MixSets.csv") %>% tibble %>% 
     select(Compound_Name, conc_uM_2019, Column, ionMode) %>% 
     rename(Concentration_uM = conc_uM_2019) %>% 
     mutate(mode = paste0(Column, "_", ionMode)) %>% select(-Column, -ionMode) %>% 
     mutate(mode = case_when(mode == "HILIC_Pos" ~ "HP",
                             mode == "HILIC_Neg" ~ "HN",
                             T ~ mode)) %>% distinct 
lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, Column, ionMode, Concentration_uM) %>% 
     filter(Column == "RP") %>% mutate(mode = "RP") %>% 
     select(Compound_Name, Concentration_uM, mode)
concentrations <- rbind(concMix, lib) %>% distinct %>% 
     filter(Compound_Name %in% cmpds) %>% 
     mutate(best = paste0(mode, "_", Compound_Name)) %>% 
     filter(best %in% bestSignal$best) 

#npo_ratios_updated %>%count(Compound_Name) %>%filter(n > 1)
#concentrations %>%count(Compound_Name) %>%filter(n > 1)

npo_ratios_dedup <- npo_ratios_updated %>%
     filter(!(Compound_Name %in% c("3OC10-L-Homoserine Lactone (HSL)", "p-Coumaric Acid")))
concentrations_dedup <- concentrations %>%
     filter(Compound_Name != "Isobutyryl-L-Carnitine")

# Append Info on
ratio_finals <- npo_ratios_dedup %>%
     left_join(concentrations_dedup, by = "Compound_Name") %>% 
     select(-mode.x, -best) %>% 
     rename(mode = mode.y) %>% 
     mutate(injVol_uL =  400, tubeVol_mL = 1)

# Perform peak area to concentration conversions...

final <- data %>% select(-method, -best) %>% 
     left_join(., ratio_finals) %>% 
     mutate(dilution_factor = case_when(method %in% c("HP", "HN") ~ 4, method == "RP" ~ 2)) %>% 
     mutate(areaBMIS = areaBMIS, 
            filtered_Half_L = filtered_vol /2,
            volume = (injVol_uL *10^-6) / filtered_vol) %>% 
     mutate(calc1 = (areaBMIS) / (RFs_h20),
            calc2 = calc1 * volume, 
            conc_nM_L = calc2 * (1 / RF_Ratio_Final) * dilution_factor) %>% 
     select(Compound_Name, sampleID, conc_nM_L)

#























