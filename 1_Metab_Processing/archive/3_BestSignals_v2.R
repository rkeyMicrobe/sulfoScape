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
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")  %>% tibble 

################################################################################
# Best Sigs according to Ingalls Standard Libraries
################################################################################

ingall <- read.csv("data/1_metaData/std_libraries/Ingalls_Lab_Standards.csv") %>%
     tibble %>%
     select(Compound_Name, Compound_Name_Original, Column, z, Priority) %>% tibble %>% 
     filter(Priority == "TRUE") %>% 
     mutate(modeChoice = case_when(Column == "HILIC"& z < 0 ~ "HN",
                                   Column == "HILIC"& z > 0 ~ "HP",
                                   T ~ "RP")
     ) %>% 
     rename(name1 = Compound_Name,
            name2 = Compound_Name_Original)

# Those having multiple modes?     
drops <- ingall %>% group_by(name1) %>% 
     reframe(n = n()) %>% 
     filter(n > 1) %>% 
     pull(name1)
ingall <- ingall %>% filter(!name1 %in% drops)

# Acquiring Durham Names that will help link to our library
list <- unique(c(ingall$name1, ingall$name2))

temp <- ingall %>% 
     rename(name = name1) %>% 
     left_join(., nameKey %>% filter(name %in% ingall$name1) %>% distinct, by = "name") %>% 
     rename(name1 = name, Durham_Name1 = Durham_Name) %>%
     rename(name = name2) %>% 
     left_join(., nameKey %>% filter(name %in% ingall$name2) %>% distinct, by = "name") %>% 
     rename(name2 = name, Durham_Name2 = Durham_Name)

temp <- temp %>% 
     mutate(Compound_Name = case_when(is.na(Durham_Name1) ~ Durham_Name2,
                                      is.na(Durham_Name2) ~ Durham_Name1,
                                      T ~ Durham_Name1
     )) %>% 
     arrange(Compound_Name)

list <- unique(temp$Compound_Name)
library %>% filter(!Compound_Name %in% list) %>% pull(Compound_Name) %>% unique

temp <- temp %>% select(Compound_Name, name1, name2, modeChoice)
write.csv(temp, "data/1_metaData/bestSignals_Ingall.csv", row.names = F)

# Link to current Library
dim(library)
library2 <- library %>% distinct %>%  
     merge(., temp %>% select(Compound_Name, modeChoice) %>% distinct, 
           by = "Compound_Name", all.x = T) 
dim(library2)

################################################################################
# Handle ATL_GOM datasets
################################################################################

n <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_HN_normMetabs.csv") %>% tibble
p <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_HP_normMetabs.csv") %>% tibble
r <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_RP_normMetabs.csv") %>% tibble

processDF <- function(data = NULL){
     dat <- data %>% 
          filter(type == "Poo") %>% 
          select(sampleID, sampleGroup, Compound_Name, BMISNormalizedArea) %>% 
          rename(AreaNorm = BMISNormalizedArea) %>% 
          group_by(Compound_Name) %>% 
          reframe(area = sum(AreaNorm, na.rm = T))
     return(dat)
}

n <- processDF(data = n) %>% rename("HN" = area)     
p <- processDF(data = p) %>% rename("HP" = area)     
r <- processDF(data = r) %>% rename("RP" = area)     

getBests <- function(HN_data, HP_data, RP_data){
     all <- full_join(HN_data, HP_data, by = "Compound_Name") %>% 
          full_join(., RP_data, by = "Compound_Name") %>% 
          mutate(detectionFreq = rowSums(!is.na(cbind(HN, HP, RP)))) %>%  # Non-missing modes)
          mutate(across(everything(), ~replace(., is.na(.), 0))) %>% 
          mutate(
               bestSignal = pmax(HN, HP, RP, na.rm = TRUE), # Highest signal
               bestSignal = case_when(
                    bestSignal == HN ~ "HN",
                    bestSignal == HP ~ "HP",
                    bestSignal == RP ~ "RP",
                    TRUE ~ "NA"
               )
          ) %>%
          arrange(desc(bestSignal))
     return(all)
}
all <- getBests(HN_data = n, HP_data = p, RP_data = r)
write.csv(temp, "data/1_metaData/bestSignals_ATLGOM_poolSums.csv", row.names = F)

# Consider the situations not found in Ingalls, Append to newest library
list <- library2 %>% filter(is.na(modeChoice)) %>% pull(Compound_Name) %>% unique
missing <- all %>% filter(Compound_Name %in% list) %>% 
     select(Compound_Name, bestSignal)

temp1 <- library2 %>% filter(is.na(modeChoice)) %>% 
     left_join(., missing, by = "Compound_Name") %>% 
     select(-modeChoice) %>% rename(modeChoice = bestSignal)
temp2 <- library2 %>% filter(!is.na(modeChoice))     

final <- rbind(temp1, temp2)
dim(final)

write.csv(final, "data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv", row.names = F)

################################################################################
# Sanity Comparison
################################################################################
# 
# hb <- read.csv("data/1_metaData/heal_ratios/NPO_RFs_Heal.csv") %>%
#      select(Precursor.Ion.Name, ionMode) %>% 
#      rename(name = Precursor.Ion.Name)
# 
# names <- nameKey %>% filter(name %in% unique(hb$name)) %>% distinct
# 
# hb <- hb %>% 
#      left_join(., names, by = "name") %>% 
#      relocate("Durham_Name", .before = everything()) %>% tibble %>% 
#      rename(oldName = name,
#             Compound_Name = Durham_Name) %>% 
#      tibble
# 
# rsk <- final %>% select(Compound_Name, modeChoice) %>% distinct %>% 
#      filter(Compound_Name %in% hb$Compound_Name)
# 
# match <- hb %>% left_join(., rsk, by = "Compound_Name") %>% 
#      mutate(match = ifelse(is.na(ionMode) | is.na(modeChoice), NA, 
#                            modeChoice == ionMode
#                            )
#             )
# temp <- match %>% filter(match == FALSE)
# 
# temp2 <- temp %>% 
#      left_join(., all %>% select(Compound_Name, bestSignal) %>% distinct,
#                by = "Compound_Name") 
#      
# 



