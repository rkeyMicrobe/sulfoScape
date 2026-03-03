lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
library("viridis")
library("flextable")


###############################################################################
# Bring MasterFiles
###############################################################################
source("./makeMasters.R")

theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

directory = "results/supplementFigs/"

###############################################################################
# Ocean Stoichiometry - Region
###############################################################################

# Function to extract element counts from the empirical formula
normalize_formula <- function(formula) {
     str_replace_all(formula, "([A-Z][a-z]?)(?!\\d)", "\\11")  
}
extract_element_count <- function(formula, element) {
     formula <- normalize_formula(formula)  # Normalize first
     match <- str_extract(formula, paste0(element, "\\d+"))  # Extract full element and number
     count <- as.numeric(str_replace(match, element, ""))  # Remove element letter, keep number
     ifelse(is.na(count), 0, count)  # Replace NA with 0
}

stoichiometryFetch <- function(dataframe = NULL){
     
     means <- dataframe %>% 
          group_by(Compound_Name) %>% 
          reframe(Compound_Name, meanVal = mean(conc_nM_L, na.rm = T)) %>% 
          distinct %>% 
          drop_na(meanVal)
     # Append Empirical Formula info
     lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
          select(Compound_Name, Empirical_Formula) %>% distinct
     dat <- means %>% left_join(., lib, by = "Compound_Name")
     # Apply extraction to dataset
     dat_stoich <- dat %>%
          mutate(
               C = sapply(Empirical_Formula, extract_element_count, element = "C"),
               N = sapply(Empirical_Formula, extract_element_count, element = "N"),
               S = sapply(Empirical_Formula, extract_element_count, element = "S")
          )
     # Calculate Stochiometrics
     stochs <- dat_stoich %>%
          mutate(
               C_contr = C * meanVal,
               N_contr = N * meanVal,
               S_contr = S * meanVal
          )  
}

atl <- stoichiometryFetch(dataframe = ATL_metab_results$metab_data_0.15) %>% mutate(Condition = "ATL")
gom <- stoichiometryFetch(dataframe = GOM_metab_results$metab_data_0.15) %>% mutate(Condition = "GOM")
npo <- stoichiometryFetch(dataframe = NPO_metab_results$metab_data_0.15) %>% mutate(Condition = "NPO")

# Apply redfield row for this table
hoRedfield <- tibble(
     Condition = "Ho_ et al. 2003",
     C = 95,
     N = 12,
     S = 1
)

bates <- tibble(
     Condition = "Bates et al. 1994",
     C = 182,
     N = 27,
     S = 1
)

matrai_susp <- tibble(
     Condition = "Matrai & Eppley 1989 (Suspended)",
     C = 224,
     N = 27,
     S = 1
)

matrai_part <- tibble(
     Condition = "Matrai & Eppley 1989 (Sinking)",
     C = 119,
     N = 17,
     S = 1
)

# Generate the region data table
dat <- rbind(atl, gom, npo) %>% 
     group_by(Condition) %>% 
     reframe(
          C = sum(C_contr, na.rm = TRUE),
          N = sum(N_contr, na.rm = TRUE),
          S = sum(S_contr, na.rm = TRUE)
     ) %>% 
     mutate(
          # Adjust S to 1
          C = C / S,
          N = N / S,
          S = S / S
     )

# Combine all rows, then calculate ratios
dat <- bind_rows(hoRedfield, bates, matrai_susp, matrai_part, dat) %>%
     mutate(
          C = round(C),
          N = round(N),
          S = round(S),
          `C:N (N=1)` = paste0(round(C / N), ":1"),
          `C:S (S=1)` = paste0(round(C / S), ":1")
     ) %>%
     select(Condition, C, N, S, `C:N (N=1)`, `C:S (S=1)`) %>% 
     mutate(Condition = case_when(
          Condition == "Matrai & Eppley 1989 (Suspended)" ~ "Matrai & Eppley 1989\n(Suspended)",
          Condition == "Matrai & Eppley 1989 (Sinking)" ~ "Matrai & Eppley 1989\n(Sinking)",
          TRUE ~ Condition
     ))

colnames(dat) <- c("Condition", "C", "N", "S", "C:N\n(N=1)", "C:S\n(S=1)")

# Generate FlexTable
ft <- flextable(dat) %>%
     theme_vanilla() %>%
     set_caption("Elemental Ratios in Different Oceans") %>%
     colformat_num(j = 2:4, digits = 0) %>%
     align(align = "left", part = "all") %>%
     autofit()
ft
save_as_image(ft,
              path = paste0(directory, "stochs_region_Table3.png"),
              webshot = "webshot2", zoom = 2, expand = 1
)

###############################################################################
# Ocean Stoichiometry - Region * Depth Category
###############################################################################

stoichiometryFetch_depths <- function(data = NULL, region = NULL){
     dataframe = data %>% 
          select(Compound_Name, sampleID, conc_nM_L, Group, depth) %>% 
          mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                       depth > 15 & depth <= 100 ~ "D2_16-100m",
                                       depth > 100 & depth <= 200 ~ "D3_101-200m",
                                       depth > 200 ~ "D4_>200m")) 
     means <- dataframe %>% 
          group_by(Compound_Name, depth_cat) %>% 
          reframe(Compound_Name, meanVal = mean(conc_nM_L, na.rm = T)) %>% 
          distinct %>% 
          drop_na(meanVal)
     
     # Append Empirical Formula info
     lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
          select(Compound_Name, Empirical_Formula) %>% distinct
     dat <- means %>% left_join(., lib, by = "Compound_Name")
     # Apply extraction to dataset
     dat_stoich <- dat %>%
          mutate(
               C = sapply(Empirical_Formula, extract_element_count, element = "C"),
               N = sapply(Empirical_Formula, extract_element_count, element = "N"),
               S = sapply(Empirical_Formula, extract_element_count, element = "S")
          )
     # Calculate Stochiometrics
     stochs <- dat_stoich %>%
          mutate(
               C_contr = C * meanVal,
               N_contr = N * meanVal,
               S_contr = S * meanVal
          )  
     
     dat <- stochs %>% mutate(Condition = region) %>% 
          group_by(Condition, depth_cat) %>% 
          reframe(C = sum(C_contr, na.rm = TRUE),
                  N = sum(N_contr, na.rm = TRUE),
                  S = sum(S_contr, na.rm = TRUE)) %>% 
          mutate( # Adjust S to 1 
               C = C / S * 1,
               N = N / S * 1,
               S = S / S * 1 
          ) %>% # Calculate proportions
          mutate('C:N' = C / N,
                 'C:S' = C / S)
     
     dat <- dat %>% 
          mutate(across(where(is.numeric), \(x) round(x, 1)))  
}
atl <- stoichiometryFetch_depths(data = ATL_metab_results$metab_data, region = "ATL")
gom <- stoichiometryFetch_depths(data = GOM_metab_results$metab_data, region = "GOM")
npo <- stoichiometryFetch_depths(data = NPO_metab_results$metab_data, region = "NPO")
data <- rbind(atl, gom, npo) %>% distinct

# Apply redfield row for this table
redfield <- tibble(Condition = "Redfield",
                   depth_cat = "N/A",
                   C = 81.5, N = 12.3, S = 1,
                   `C:N` = C / N, 
                   `C:S` = C / S,
)
depthOrder <- c("D1_0-15m", "D2_16-100m", "D3_101-200m", "D4_>200m")
conditionOrder <- c("NPO", "GOM", "ATL")

depth_expanded <- expand.grid(Condition = unique(data$Condition), 
                              depth_cat = depthOrder, ratio = c("C:N", "C:S"))
data_complete <- data %>% 
     select(Condition, depth_cat, `C:N`, `C:S`) %>%
     pivot_longer(cols = -c(Condition, depth_cat), names_to = "ratio", values_to = "value") %>%
     right_join(depth_expanded, by = c("Condition", "depth_cat", "ratio")) %>%
     mutate(value = ifelse(is.na(value), 0, value),  # Fill missing values with 0
            depth_cat = factor(depth_cat, levels = depthOrder),  
            Condition = factor(Condition, levels = conditionOrder)) 


data <- data %>% select(Condition, depth_cat, 'C:N', 'C:S') %>% 
     pivot_longer(., cols = -c(Condition, depth_cat), names_to = "ratio", values_to = "value")

redfield <- tibble(Condition = "Redfield",
                   depth_cat = "N/A",
                   C = 81.5, N = 12.3, S = 1,
                   `C:N` = C / N, 
                   `C:S` = C / S)
redfield_CN <- round(redfield$`C:N`, 2)
redfield_CS <- round(redfield$`C:S`, 2)



# ggplot
p <- ggplot(data_complete, aes(y = factor(depth_cat, levels = rev(depthOrder)), 
                               x = value, fill = Condition)) +
     geom_col(position = position_dodge(width = 0.6), color = "black") +
     facet_wrap(~ratio, scales = "free_x") + 
     labs(subtitle = paste0("Redfield Benchmarks \nC:N = ", redfield_CN, " | C:S = ", redfield_CS),
          y = "Depth Category", 
          x = "Ratio Value") +
     scale_fill_manual(values = c("#FFB86C", "#BD93F9", "#50FA7B")) +
     scale_y_discrete(expand = expansion(mult = 0.15)) +
     facet_theme +
     theme(axis.title.y = element_blank(),
           legend.position = "bottom"
     )
p
svg(paste0(directory, "regionDepths_stochs.svg"), width = 6, height = 4); p; dev.off()

#---------------------------------------------------------------------------------------------------------------

# All
 getNumbers <- function(data = NULL, ocean = NULL){
      dataframe = data %>% 
           select(Compound_Name, sampleID, conc_nM_L, Group, depth) %>% 
           mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                        depth > 15 & depth <= 100 ~ "D2_16-100m",
                                        depth > 100 & depth <= 200 ~ "D3_101-200m",
                                        depth > 200 ~ "D4_>200m")) 
      
      dataframe %>% 
           filter(!conc_nM_L <= 0) %>% 
           select(Compound_Name, depth_cat) %>% 
           distinct %>% 
           group_by(depth_cat) %>% 
           reframe(n_cmpds = n()) %>% 
           mutate(cruise = ocean)   
 }
 data = rbind(getNumbers(data = NPO_metab_results$metab_data, ocean = "NPO"),
              getNumbers(data = GOM_metab_results$metab_data, ocean = "GOM"),
              getNumbers(data = ATL_metab_results$metab_data, ocean = "ATL")
 ) %>% 
      pivot_wider(., names_from = "cruise", values_from = "n_cmpds")
 
# SULFUR GROUP DOM
 getNumbers <- function(data = NULL, ocean = NULL){
      dataframe = data %>% 
           select(Compound_Name, sampleID, conc_nM_L, Group, Sulfur_Group, depth) %>% 
           mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                        depth > 15 & depth <= 100 ~ "D2_16-100m",
                                        depth > 100 & depth <= 200 ~ "D3_101-200m",
                                        depth > 200 ~ "D4_>200m")) 
      
      dataframe %>% 
           filter(Sulfur_Group %in% c("Sulfide", "Sulfonate")) %>% 
           filter(!conc_nM_L <= 0) %>% 
           select(Compound_Name, Sulfur_Group, depth_cat) %>% distinct %>% 
           group_by(Sulfur_Group, depth_cat) %>% 
           reframe(n_cmpds = n()) %>% 
           mutate(cruise = ocean)   
 }
 data = rbind(getNumbers(data = NPO_metab_results$metab_data, ocean = "NPO"),
              getNumbers(data = GOM_metab_results$metab_data, ocean = "GOM"),
              getNumbers(data = ATL_metab_results$metab_data, ocean = "ATL")
 ) %>% 
      pivot_wider(., names_from = "depth_cat", values_from = "n_cmpds") %>% 
      arrange(Sulfur_Group, cruise)
 
 ###############################################################################
# Ocean Stoichiometry - Latitude, NPO
###############################################################################

dataframe = NPO_metab_results$metab_data

# Assign Depth Categories
dataframe <- dataframe %>% select(Compound_Name, sampleID, conc_nM_L, Group, latitude, depth) %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                  depth > 15 & depth <= 100 ~ "D2_16-100m")) %>% 
     mutate(lat = round(latitude)) %>% 
     select(Compound_Name, sampleID, conc_nM_L, lat, depth_cat) %>%
     drop_na(conc_nM_L)

# Append Empirical Formula info
lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, Empirical_Formula) %>% distinct
dat <- dataframe %>% left_join(., lib, by = "Compound_Name")

index <- dat %>% select(Compound_Name, Empirical_Formula) %>% distinct %>% 
     mutate(
          C = sapply(Empirical_Formula, extract_element_count, element = "C"),
          N = sapply(Empirical_Formula, extract_element_count, element = "N"),
          S = sapply(Empirical_Formula, extract_element_count, element = "S")
     ) %>% 
     select(-Empirical_Formula)

# Generate a table showing stiochiometric contributions
dat_stoich <- dat %>% 
     left_join(., index, by = "Compound_Name") %>% 
     mutate(
          C_contr = C * conc_nM_L,
          N_contr = N * conc_nM_L,
          S_contr = S * conc_nM_L
     ) %>% 
     group_by(sampleID) %>% 
     reframe(lat, depth_cat, 
             C = sum(C_contr, na.rm = TRUE),
             N = sum(N_contr, na.rm = TRUE),
             S = sum(S_contr, na.rm = TRUE)
     ) %>% distinct %>% 
     mutate(C = C / S * 1,
            N = N / S * 1,
            S = S / S * 1 
     ) %>%
     mutate('CN' = C / N,
            'CS' = C / S) %>% 
     mutate(lat = factor(lat, levels = sort(unique(lat), decreasing = F))) %>% 
     select(-C, -N, -S) %>% 
     mutate(lat = as.character(lat))

summary <- dat_stoich %>%
     group_by(lat, depth_cat) %>%
     summarise(
          CN_mean = mean(CN, na.rm = TRUE),
          CN_sd = sd(CN, na.rm = TRUE),
          CS_mean = mean(CS, na.rm = TRUE),
          CS_sd = sd(CS, na.rm = TRUE)
     ) %>%
     ungroup()

plotDat <- summary %>%
     pivot_longer(cols = c(CN_mean, CS_mean), names_to = "Ratio", values_to = "Value") %>%
     pivot_longer(cols = c(CN_sd, CS_sd), names_to = "Ratio_sd", values_to = "SD") %>%
     filter(substring(Ratio_sd, 1, 2) == substring(Ratio, 1, 2)) %>%  # Match SD to mean values
     select(-Ratio_sd)

p1 <- ggplot(plotDat %>% filter(Ratio == "CN_mean"), 
             aes(y = lat, x = Value, fill = depth_cat)) +
     geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.6) +  # Bars with black outline
     geom_errorbar(aes(xmin = Value - SD, xmax = Value + SD), 
                   width = 0.2, position = position_dodge(width = 0.8)) + 
     scale_fill_manual(values = c("D1_0-15m" = "#A7C7E7", "D2_16-100m" = "#44475A")) +
     labs(title = "C:N Trend", x = "C:N Ratio", y = "Latitude (°N)") +
     theme(legend.position = "none") +
     facet_theme +
     coord_flip()
p2 <- ggplot(plotDat %>% filter(Ratio == "CS_mean"), 
             aes(y = lat, x = Value, fill = depth_cat)) +
     geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.6) +  # Bars with black outline
     geom_errorbar(aes(xmin = Value - SD, xmax = Value + SD), 
                   width = 0.2, position = position_dodge(width = 0.8)) + 
     scale_fill_manual(values = c("D1_0-15m" = "#A7C7E7", "D2_16-100m" = "#44475A")) +
     labs(title = "C:S Trend", x = "C:S Ratio", y = "Latitude (°N)") +
     theme(legend.position = "none") +
     facet_theme+
     coord_flip()
p <- plot_grid(p1, p2, ncol = 1)
p
directory = "results/supplementFigs/"
svg(paste0(directory, "npoLats_stochs_stats.svg"), width = 6, height = 4); p; dev.off()

getNumbers <- function(data = NULL, ocean = NULL){
     dataframe = data %>% 
          select(Compound_Name, sampleID, conc_nM_L, Group, depth, latitude) %>% 
          mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                       depth > 15 & depth <= 100 ~ "D2_16-100m",
                                       depth > 100 & depth <= 200 ~ "D3_101-200m",
                                       depth > 200 ~ "D4_>200m")) %>% 
          mutate(lat = round(latitude)) %>% 
          filter(depth_cat %in% c("D1_0-15m", "D2_16-100m")) %>% 
          select(Compound_Name, sampleID, conc_nM_L, lat, depth_cat) %>%
          drop_na(conc_nM_L)
     
     dataframe %>% 
          filter(!conc_nM_L <= 0) %>% 
          select(Compound_Name, depth_cat, lat) %>% 
          distinct %>% 
          group_by(depth_cat, lat) %>% 
          reframe(n_cmpds = n()) %>% 
          mutate(cruise = ocean)   
}
data = getNumbers(data = NPO_metab_results$metab_data, ocean = "NPO") %>% 
     pivot_wider(., names_from = "cruise", values_from = "n_cmpds")

p <- ggplot(data, aes(x = NPO, y = as.character(lat), fill = depth_cat)) +
     geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.6) + 
     geom_text(aes(label = NPO), position = position_dodge(width = 0.8),
          vjust = 0.5, hjust = -0.5, angle = 90, size = 3) +
     xlim(0,100) +
     scale_fill_manual(values = c("D1_0-15m" = "#A7C7E7", "D2_16-100m" = "#44475A")) +
     labs(title = "Number of Metabs", x = "C:S Ratio", y = "Latitude (°N)") +
     theme(legend.position = "none") +
     facet_theme +
     coord_flip()
p
svg(paste0(directory, "npoLats_stochs_metaNumber.svg"), width = 6, height = 4); p; dev.off()


###############################################################################
# Ocean Stoichiometry - Depth, GOM
###############################################################################

dataframe = GOM_metab_results$metab_data

# Assign Depth Categories
dataframe <- dataframe %>% 
     select(Compound_Name, sampleID, conc_nM_L, depth) %>% 
     drop_na(conc_nM_L)

# Append Empirical Formula info
lib <- read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, Empirical_Formula) %>% distinct
dat <- dataframe %>% left_join(., lib, by = "Compound_Name")
index <- dat %>% select(Compound_Name, Empirical_Formula) %>% distinct %>% 
     mutate(
          C = sapply(Empirical_Formula, extract_element_count, element = "C"),
          N = sapply(Empirical_Formula, extract_element_count, element = "N"),
          S = sapply(Empirical_Formula, extract_element_count, element = "S")
     ) %>% 
     select(-Empirical_Formula)

# Generate a table showing stiochiometric contributions
dat_stoich <- dat %>% 
     left_join(., index, by = "Compound_Name") %>% 
     mutate(
          C_contr = C * conc_nM_L,
          N_contr = N * conc_nM_L,
          S_contr = S * conc_nM_L
     ) %>% 
     group_by(sampleID) %>% 
     reframe(depth, 
             C = sum(C_contr, na.rm = TRUE),
             N = sum(N_contr, na.rm = TRUE),
             S = sum(S_contr, na.rm = TRUE)
     ) %>% distinct %>% 
     mutate(C = C / S * 1,
            N = N / S * 1,
            S = S / S * 1 
     ) %>%
     mutate('CN' = C / N,
            'CS' = C / S) %>% 
     select(-C, -N, -S) 

summary <- dat_stoich %>%
     group_by(depth) %>%
     summarise(
          CN_mean = mean(CN, na.rm = TRUE),
          CN_sd = sd(CN, na.rm = TRUE),
          CS_mean = mean(CS, na.rm = TRUE),
          CS_sd = sd(CS, na.rm = TRUE)
     ) %>%
     ungroup()


plotDat <- summary %>%
     pivot_longer(cols = c(CN_mean, CS_mean), names_to = "Ratio", values_to = "Value") %>%
     pivot_longer(cols = c(CN_sd, CS_sd), names_to = "Ratio_sd", values_to = "SD") %>%
     filter(substring(Ratio_sd, 1, 2) == substring(Ratio, 1, 2)) %>%  # Match SD to mean values
     select(-Ratio_sd)

order <- rev(unique(plotDat$depth))

p <- ggplot(plotDat, aes(y = factor(as.character(depth), levels = order), 
                    x = Value, fill = Ratio)) +
     geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.6) +  # Bars with black outline
     geom_errorbar(aes(xmin = Value - SD, xmax = Value + SD), 
                   width = 0.2, position = position_dodge(width = 0.8)) + 
     facet_wrap(~Ratio, scales = "free_x") +
     scale_fill_manual(values = c("CN_mean" = "#ffb3dc", "CS_mean" = "#FFB86C")) +
     theme(legend.position = "none") +
     facet_theme +
     labs(title = "Ratio Trends", x = "Ratio Value", y = "Depth (m)")
p
directory = "results/supplementFigs/"
svg(paste0(directory, "gomDepths_stochs_stats.svg"), width = 6, height = 5); p; dev.off()

getNumbers <- function(data = NULL, ocean = NULL){
     dataframe = data %>% 
          select(Compound_Name, sampleID, conc_nM_L, Group, depth, latitude) %>% 
          mutate(lat = round(latitude)) %>% 
          select(Compound_Name, sampleID, conc_nM_L, lat, depth) %>%
          drop_na(conc_nM_L)
     
     dataframe %>% 
          filter(!conc_nM_L <= 0) %>% 
          select(Compound_Name, depth) %>% 
          distinct %>% 
          group_by(depth) %>% 
          reframe(n_cmpds = n()) %>% 
          mutate(cruise = ocean)   
}

data = getNumbers(data = GOM_metab_results$metab_data, ocean = "GOM") %>% 
     pivot_wider(., names_from = "cruise", values_from = "n_cmpds")

depth_levels <- rev(c("15", "40", "55", "64", "68", "200", "500"))  # as strings
depth_colors <- c(
     "#084594",  
     "#2171b5",  
     "#4292c6",  
     "#6baed6",  
     "#9ecae1",  
     "#c6dbef",  
     "#deebf7"
)

# Plot
p <- ggplot(data, aes(x = GOM, 
                      y = factor(as.character(depth), levels = depth_levels),
                      fill = factor(as.character(depth), levels = depth_levels))
            ) +
     geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.6) +
     geom_text(aes(label = GOM), position = position_dodge(width = 0.8),
               vjust = 0.5, hjust = -0.5, size = 3) +
     xlim(0, 220) +
     scale_fill_manual(values = depth_colors) +
     labs(title = "Number of Metabs", x = "Metabolite Number", y = "Depth (m)") +
     theme(legend.position = "none") +
     facet_theme
p
svg(paste0(directory, "gomDepth_stochs_metaNumber.svg"), width = 3.5, height = 5); p; dev.off()