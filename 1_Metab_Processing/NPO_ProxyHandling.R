lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())

################################################################################
# 
################################################################################

proxyIndex <- read.csv("data/1_metaData/proxyChoice.csv") %>% tibble
choiceIndex <- read.csv("data/1_metaData/bestSignals_index.csv") %>% tibble
choices = unique(proxyIndex$ProxyPick) 
cmpds = unique(proxyIndex$Compound_Name) 
all <- c(choices, cmpds)

choiceIndex %>% filter(Compound_Name %in% cmpds)

################################################################################
# Fetch Ratios - ATLGOM and NPO
################################################################################

# Get Metab Ratio Info
fetchRatios <- function(dir1, dir2, method = NULL){
     x1 <- read.csv(paste0(dir1, method, "_Ratios.csv")) %>% tibble %>% 
          select(Compound_Name, RFs_h20, RF_Ratio_Final) %>% 
          rename(G2_Ratio = RF_Ratio_Final)
     x2 <- read.csv(paste0(dir2, method, "_Ratios.csv")) %>% tibble %>% 
          select(Compound_Name, RFs_h20, RF_Ratio_Final) %>% 
          rename(Other_Ratio = RF_Ratio_Final,
                 Other_RFs_h20 = RFs_h20)
     
     x <- list(NPO = x1, ATLGOM = x2)
     return(x)
}
RP_Ratios  <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                          dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "RP")
HP_Ratios <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                         dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "HP")
HN_Ratios <- fetchRatios(dir1 = "data/NPO_G2/ratios/RF_NPO_G2_", 
                         dir2 = "data/ATL_GOM/ratios/RF_ATL_GOM_", method = "HN")

# HILIC Positive
ratio_NPO = HP_Ratios$NPO
ratio_ATLGOM =  HP_Ratios$ATLGOM

info <- choiceIndex %>% filter(Compound_Name %in% cmpds) %>% filter(Ingall_Choice == "HP")
info_cmpds <- info$Compound_Name
missings <- setdiff(info_cmpds, ratio_NPO$Compound_Name)
filler <- tibble(
     Compound_Name = missings,
     RFs_h20 = NA_real_,
     G2_Ratio = NA_real_
)

dat1 <- rbind(ratio_NPO %>% filter(Compound_Name %in% cmpds),
              filler)
dat2 <- ratio_ATLGOM %>% filter(Compound_Name %in% choices) %>% rename(ProxyPick = Compound_Name)

list <- unique(dat1$Compound_Name)
temp <- proxyIndex %>% 
     filter(Compound_Name %in% list)

list <- unique(dat1$Compound_Name)

newInfo <- temp %>% 
     left_join(., dat1, by = "Compound_Name") %>% 
     left_join(., dat2, by = "ProxyPick")

HP_Ratios$NPO <- 
     rbind(HP_Ratios$NPO, filler) %>% 
     left_join(newInfo %>% 
                    select(Compound_Name, Other_RFs_h20, Other_Ratio), 
               by = "Compound_Name") %>%
     mutate(
          RFs_h20 = coalesce(Other_RFs_h20, RFs_h20),
          G2_Ratio = coalesce(Other_Ratio, G2_Ratio)
     ) %>%
     select(-Other_RFs_h20, -Other_Ratio) %>% 
     distinct

# HILIC negative
ratio_NPO = HN_Ratios$NPO
ratio_ATLGOM =  HN_Ratios$ATLGOM

# What all do we need?
info <- choiceIndex %>% filter(Compound_Name %in% cmpds) %>% filter(Ingall_Choice == "HN")
info_cmpds <- info$Compound_Name
missings <- setdiff(info_cmpds, ratio_NPO$Compound_Name)
filler <- tibble(
     Compound_Name = missings,
     RFs_h20 = NA_real_,
     G2_Ratio = NA_real_
)


dat1 <- rbind(ratio_NPO %>% filter(Compound_Name %in% cmpds), filler)
dat2 <- ratio_ATLGOM %>% filter(Compound_Name %in% choices) %>% rename(ProxyPick = Compound_Name)

list <- unique(dat1$Compound_Name)
temp <- proxyIndex %>% filter(Compound_Name %in% list)

dat3 <- HP_Ratios$ATLGOM %>% filter(Compound_Name %in% choices) %>% rename(ProxyPick = Compound_Name)
dat4 <- RP_Ratios$ATLGOM %>% filter(Compound_Name %in% choices) %>% rename(ProxyPick = Compound_Name)

newInfo <- temp %>% 
     left_join(., dat1, by = "Compound_Name") %>% 
     left_join(., dat3, by = "ProxyPick") %>% 
     left_join(., dat2, by = "ProxyPick") %>% 
     mutate(
          Other_RFs_h20 = coalesce(Other_RFs_h20.x, Other_RFs_h20.y),
          Other_Ratio = coalesce(Other_Ratio.x, Other_Ratio.y)
     ) %>% 
     select(Compound_Name, ProxyPick, Other_RFs_h20, Other_Ratio) %>% 
     distinct
     
HN_Ratios$NPO <- 
     rbind(HN_Ratios$NPO, filler) %>%
     left_join(newInfo %>% 
                    select(Compound_Name, Other_RFs_h20, Other_Ratio), 
               by = "Compound_Name") %>%
     mutate(
          RFs_h20 = coalesce(Other_RFs_h20, RFs_h20),
          G2_Ratio = coalesce(Other_Ratio, G2_Ratio)
     ) %>%
     select(-Other_RFs_h20, -Other_Ratio)

# Reverse Phase
dat1 <- RP_Ratios$NPO %>% filter(Compound_Name %in% cmpds)
temp <- HN_Ratios$NPO

################################################################################
# Fetch other pieces for peak to nmol conversions
# Add onto the Ratio dataframes
################################################################################

# Get Metab concentrations from 2019 
mixKey <- read.csv("data/1_metaData/metab_MixSets.csv") %>% tibble %>%
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
concentrations <- rbind(mixKey, lib) %>% distinct

# Merge Concentration Info (and any other pertinent variables) to NPO Ratios
process_ratios <- function(ratios, method, stocks) {
     
     x <- stocks %>% filter(mode == method) 
     ratios %>%
          left_join(x, by = "Compound_Name") %>% 
          arrange(Compound_Name) %>% 
          # Injection info from LCM work
          mutate(injVol_uL =  400, 
                 tubeVol_mL = 1)
     
}

NPO_info <- list(
     RP = process_ratios(RP_Ratios$NPO, method = "RP", stocks = concentrations),
     HP = process_ratios(HP_Ratios$NPO, method = "HP", stocks = concentrations),
     HN = process_ratios(ratios = HN_Ratios$NPO, method = "HN", stocks = concentrations)
)
temp <- NPO_info$HN
temp <- HN_Ratios$NPO
################################################################################
# Fetch BMIS PEAKS - NPO
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

sort(unique(npo_postBMIS$Compound_Name))
npo_postBMIS %>% filter(Compound_Name == "Cysteinolic Acid")
################################################################################
# Make new concentration dataframes
################################################################################

processDF <- function(ratioInfo = NULL, BMIS_df = NULL, mode = NULL){
     postBMIS = BMIS_df %>% 
          left_join(., ratioInfo, by = "Compound_Name") %>% 
          filter(str_detect(sampleID,  "Smp_")) %>% distinct %>% 
          tibble %>% print 
     
     final <- postBMIS %>% drop_na(G2_Ratio) %>% 
          select(Compound_Name, sampleID, Concentration_uM, injVol_uL, filtered_vol, 
                 RFs_h20, G2_Ratio, areaBMIS) %>% 
          mutate(areaBMIS = areaBMIS, 
                 filtered_Half_L = filtered_vol /2,
                 volume = (injVol_uL *10^-6) / filtered_vol,
                 dilution_factor = 4
          ) %>% 
          mutate(calc1 = (areaBMIS) / (RFs_h20),
                 calc2 = calc1 * volume, 
                 conc_nM_L = calc2 * (1 / G2_Ratio) * dilution_factor
          ) %>% 
          select(Compound_Name, sampleID, conc_nM_L)
     return(final)   
}


final <- processDF(ratioInfo =  NPO_info$HP, 
                   BMIS_df = npo_postBMIS %>% filter(method == "HP"), 
                   mode = "HP")
data_path = "data/NPO_G2/NPO_G2_HP_"
write.csv(final, paste0(data_path, "Metabs_ProxyApplied.csv"), row.names = F)

final <- processDF(ratioInfo =  NPO_info$HN, 
                  BMIS_df = npo_postBMIS %>% filter(method == "HN"), 
                  mode = "HN")
data_path = "data/NPO_G2/NPO_G2_HN_"
write.csv(final, paste0(data_path, "Metabs_ProxyApplied.csv"), row.names = F)

# No modifications for Reverse Phase


