lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())

library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble 
ingall <- read.csv("data/1_metaData/std_libraries/Ingalls_Lab_Standards.csv") %>% tibble 
iSTDs <- library %>% filter(Group == "Internal Standard") %>% pull(Compound_Name)
nameKey <- read.csv("data/1_metaData/nameKey_update.csv")  %>% tibble 

################################################################################
# Create the Best Signal Dataframe that will link to Durham Library Naming Format
################################################################################

#-------------------------------------------------------------------------------
# Ingall Choice
#-------------------------------------------------------------------------------

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

temp <- temp %>% select(Compound_Name, name1, name2, modeChoice) %>% 
     rename(IngallName1 = name1, IngallName2 = name2)
write.csv(temp, "data/1_metaData/BestSignals_Ingall.csv", row.names = F)

length(unique(library$Compound_Name))
length(unique(temp$Compound_Name))
library %>% filter(!Compound_Name %in% list) %>% pull(Compound_Name) %>% unique

ingallChoices = temp %>% rename(Ingall_Choice = modeChoice)

#-------------------------------------------------------------------------------
# ATLGOM - Highest Pool
#-------------------------------------------------------------------------------

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
pools <- getBests(HN_data = n, HP_data = p, RP_data = r)

poolsATLGOM <- pools %>% select(Compound_Name, bestSignal) %>% 
     rename(ATLGOM_Choice = bestSignal)

#-------------------------------------------------------------------------------
# NPO - Highest Pool
#-------------------------------------------------------------------------------

n <- read.csv("data/NPO_G2/bmis/BMIS_NPO_G2_HN_normMetabs.csv") %>% tibble
p <- read.csv("data/NPO_G2/bmis/BMIS_NPO_G2_HP_normMetabs.csv") %>% tibble
r <- read.csv("data/NPO_G2/bmis/BMIS_NPO_G2_RP_normMetabs.csv") %>% tibble

n <- processDF(data = n) %>% rename("HN" = area)     
p <- processDF(data = p) %>% rename("HP" = area)     
r <- processDF(data = r) %>% rename("RP" = area)     

pools <- getBests(HN_data = n, HP_data = p, RP_data = r)

poolsNPOG2 <- pools %>% select(Compound_Name, bestSignal) %>% 
     rename(NPO_Choice = bestSignal)

#-------------------------------------------------------------------------------
# Append all together to make the index
#-------------------------------------------------------------------------------

choices <- ingallChoices %>% 
     left_join(., poolsATLGOM, by = "Compound_Name") %>% 
     left_join(., poolsNPOG2, by = "Compound_Name") %>% 
     select(Compound_Name, Ingall_Choice, ATLGOM_Choice, NPO_Choice,
            IngallName1, IngallName2) %>% 
     arrange(NPO_Choice)

write.csv(choices, "data/1_metaData/bestSignals_index.csv", row.names = F)

