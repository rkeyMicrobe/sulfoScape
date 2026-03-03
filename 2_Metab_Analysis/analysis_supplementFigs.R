lapply(names(sessionInfo()$otherPkgs), function(pkg) {
     detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
})
rm(list = ls()); gc(); cat("\014")

# LOAD PACKAGES, THEME
library("tidyverse")
library("cowplot")
library("flextable")
library("vegan")
library("viridis")

###############################################################################
# Bring MasterFiles and Functions
###############################################################################
source("./makeMasters.R")

applyDepthCats <- function(data = NULL){
     data %>% tibble() %>% 
          mutate(depth_cat = case_when(depth <= 15 ~ "0-15m",
                                       depth > 15 & depth <= 100 ~ "16-100m",
                                       depth > 100 & depth <= 200 ~ "101-200m",
                                       depth > 200 ~ ">200m"))
}

theme_set(theme_cowplot())
facet_theme <- theme(strip.background = element_rect(fill = "#44475A"),
                     strip.text = element_text(color = "white", face = "bold"))

dat_dir <- "data/"
dat_dir <- "results/analysisFigs/"
sup_dir <- "results/supplementFigs/"

#################################################################################
# Supplement Table - Metabolite MZ info    X
#################################################################################

s_library <-  read.csv("data/1_metaData/std_libraries/Durham_Metab_Library_v7.csv") %>% tibble %>% 
     select(Compound_Name, Sulfur_Group, Column, ionMode, Retention_Time) %>% distinct %>% 
     filter(!is.na(Sulfur_Group))

# Name Cleans
anionNames <- read.csv("data/1_metaData/nameKey_Anion.csv")

# Clean names
s_library <- s_library %>% 
     left_join(., anionNames, by = "Compound_Name") %>% 
     mutate(Compound_Name = coalesce(New_Name, Compound_Name)) %>% 
     select(-New_Name)


unique(s_library$Sulfur_Group)
length(unique(s_library$Compound_Name))

flex_table <- s_library %>%
     arrange(Sulfur_Group, Compound_Name, desc(Retention_Time)) %>%
     flextable() %>%
     set_header_labels(
          Compound_Name = "Compound\nName",
          Sulfur_Group = "Sulfur\nCategory",
          Column = "Column\nType",
          ionMode = "Ion\nMode",
          Retention_Time = "Retention\nTime (min)"
     ) %>%
     theme_vanilla() %>%
     autofit() %>%
     width(j = c(1, 1.5, 2, 2.5, 3), width = 2) %>%
     height_all(height = 0.3) %>%
     bg(bg = "white", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfonium", j = 2, bg = "#A6CEE3", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfonate", j = 2, bg = "#FDBF6F", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Disulfide", j = 2, bg = "#FFFF99", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfate ester", j = 2, bg = "#CAB2D6", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfide", j = 2, bg = "#FB9A99", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfoxide", j = 2, bg = "#C1E1C1", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Thiol", j = 2, bg = "#D2B48C", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfinic acid", j = 2, bg = "#CCCCCC", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Sulfoxonium", j = 2, bg = "#FADADD", part = "body") %>%
     bg(i = ~ Sulfur_Group == "Thione", j = 2, bg = "#8DD3C7", part = "body") %>%
     color(part = "header", color = "black") %>%
     bold(part = "header", bold = TRUE)
flex_table
save_as_image(
     flex_table,
     path = "results/supplementFigs/sulfur_RT_Table2.png",
     webshot = "webshot2",
     zoom = 2,
     expand = 1
)

###############################################################################
# Bring BMIS'd areas to find unique metabs - ATLGOM
###############################################################################
n <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_HN_normMetabs.csv") %>% tibble
p <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_HP_normMetabs.csv") %>% tibble
r <- read.csv("data/ATL_GOM/bmis/BMIS_ATL_GOM_RP_normMetabs.csv") %>% tibble

best <- read.csv("data/1_metaData/bestSignals_index.csv") %>% 
     select(Compound_Name, Ingall_Choice) %>% distinct
bsHP = best %>% filter(Ingall_Choice == "HP") %>% distinct %>% pull(Compound_Name) 
bsHN = best %>% filter(Ingall_Choice == "HN") %>% distinct %>% pull(Compound_Name)
bsRP = best %>% filter(Ingall_Choice == "RP") %>% distinct %>% pull(Compound_Name)

# Apply to each mode and bind together
r <- r %>% filter(Compound_Name %in% bsRP)
p <- p %>% filter(Compound_Name %in% bsHP)
n <- n %>% filter(Compound_Name %in% bsHN)


processDF <- function(data = NULL){
     dat <- data %>% 
          filter(BMISNormalizedArea > 0) %>% 
          select(sampleID, sampleGroup, Compound_Name, BMISNormalizedArea) %>% 
          rename(AreaNorm = BMISNormalizedArea) %>% 
          group_by(Compound_Name) %>% 
          reframe(area = sum(AreaNorm, na.rm = T))
     return(dat)
}

n <- processDF(data = n) %>% rename("HN" = area)   
p <- processDF(data = p) %>% rename("HP" = area)    
r <- processDF(data = r) %>% rename("RP" = area)

length(unique(c(n$Compound_Name, p$Compound_Name, r$Compound_Name)))

###############################################################################
###############################################################################
# Metabolites Detected Table and Number of Unique Metabolites  X
###############################################################################
###############################################################################

# Choose either all depths or 0-15m for this part. 
atl <- ATL_metab_results$metab_data %>% mutate(set = "ATL") %>% applyDepthCats()
gom <- GOM_metab_results$metab_data %>% mutate(set = "GOM") %>% applyDepthCats()
npo <- NPO_metab_results$metab_data %>% mutate(set = "NPO") %>% applyDepthCats() %>% 
     filter(!Compound_Name %in% c("Acetyl Coenzyme A (Acetyl-CoA)", "Arachidonate, D8"))


atl_s_count <- atl %>%
     filter(!is.na(conc_nM_L), conc_nM_L > 0, Sulfur_Group != "NA") %>%
     distinct(Compound_Name) %>%
     nrow()

gom_s_count <- gom %>%
     filter(!is.na(conc_nM_L), conc_nM_L > 0, Sulfur_Group != "NA") %>%
     distinct(Compound_Name) %>%
     nrow()

npo_s_count <- npo %>%
     filter(!is.na(conc_nM_L), conc_nM_L > 0, Sulfur_Group != "NA") %>%
     distinct(Compound_Name) %>%
     nrow()

getTally <- function(data = NULL){
     dat = data %>% 
          group_by(Compound_Name) %>% 
          reframe(Compound_Name, Group, Sulfur_Group, 
                  mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
          distinct %>% 
          filter(!mean_cmpd_conc <= 0 ) %>% 
          ungroup %>% group_by(Group) %>% 
          summarise(cruise = n())
     return(dat)
}
l1 = getTally(data = atl) %>% rename(ATL = cruise)
l1 = l1 %>% bind_rows(tibble(Group = "Total_Number", ATL = sum(l1$ATL)))

l2 = getTally(data = gom) %>% rename(GOM = cruise)
l2 = l2 %>% bind_rows(tibble(Group = "Total_Number", GOM = sum(l2$GOM)))

l3 = getTally(data = npo) %>% rename(NPO = cruise)
l3 = l3 %>% bind_rows(tibble(Group = "Total_Number", NPO = sum(l3$NPO)))

l4 = rbind(atl %>% select(Compound_Name, Group, conc_nM_L),
      gom %>% select(Compound_Name, Group, conc_nM_L),
      npo %>% select(Compound_Name, Group, conc_nM_L)) %>% 
     group_by(Compound_Name) %>% 
     reframe(Compound_Name, Group, 
             mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
     distinct %>% 
     filter(!mean_cmpd_conc <= 0 ) %>% 
     ungroup %>% group_by(Group) %>% 
     summarise('Combined' = n())
l4 = l4 %>% bind_rows(tibble(Group = "Total_Number", Combined = sum(l4$Combined)))


Combined = l1 %>% 
     left_join(., l2, by = "Group") %>% 
     left_join(., l3, by = "Group") %>% 
     left_join(., l4, by = "Group") %>% distinct
Combined[is.na(Combined)] <- 0

table = flextable(Combined) %>%
     set_caption("Detected Metabolite Group Counts Across Regions") %>% 
     autofit() %>%
     theme_vanilla() %>% 
     align(align = "center", part = "all") %>% 
     bg(j = "Combined", bg = "#F5F5F5") %>% 
     bg(i = ~ Group == "Total_Number", bg = "#FFFACD") 
     
table

save_as_image(table,
              path = paste0(sup_dir, "Supplements/TallyTable_detected_surface.png"),
              webshot = "webshot2", zoom = 2, expand = 1
)

# Number detected across all cruises, all depths
getCmpds <- function(dat1, dat2, dat3){
     rbind(dat1 %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           dat2 %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           dat3 %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name)
     ) %>% distinct
}
cmpds <- getCmpds(dat1 = ATL_metab_results$metab_data, dat2 = GOM_metab_results$metab_data, dat3 = NPO_metab_results$metab_data)
cmpds <- cmpds %>% left_join(., group_index, by = "Compound_Name") %>% select(-Sulfur_Group)
cmpds %>% group_by(Group) %>% reframe(count = n())

# Number detected across all cruises, 0-15m depths
cmpds <- getCmpds(dat1 = ATL_metab_results$metab_data_0.15, dat2 = GOM_metab_results$metab_data_0.15, dat3 = NPO_metab_results$metab_data_0.15)
cmpds <- cmpds %>% left_join(., group_index, by = "Compound_Name") %>% select(-Sulfur_Group)
cmpds %>% group_by(Group) %>% reframe(count = n())


#################################################################################
# ORINDATIONS - NMDS    X
#################################################################################
atl <- ATL_metab_results$metab_data %>% mutate(set = "ATL") %>% applyDepthCats()
gom <- GOM_metab_results$metab_data %>% mutate(set = "GOM") %>% applyDepthCats()
npo <- NPO_metab_results$metab_data %>% mutate(set = "NPO") %>% applyDepthCats() %>% 
     filter(!Compound_Name %in% c("Acetyl Coenzyme A (Acetyl-CoA)", "Arachidonate, D8"))

meta <- read.csv("data/1_metaData/sample_info/ATL_GOM_sampleMeta.csv") %>%
     filter(collectionSet == "ATL") %>% 
     select(sampleID, collectionDate) %>% drop_na()


getNMDS <- function(data = NULL) {
     nmds_input <- data  %>% 
          select(sampleID, Compound_Name, conc_nM_L) %>% distinct %>% 
          pivot_wider(names_from = Compound_Name, values_from = conc_nM_L) %>% 
          mutate(across(where(is.numeric), ~ replace_na(.x, 0)))    
     
     sample_ids <- nmds_input$sampleID
     compounds <- nmds_input %>%
          select(-sampleID) %>%
          select(where(~ sum(.x) > 0))
     nmds_input_filtered <- bind_cols(tibble(sampleID = sample_ids), compounds)
     
     # Extract matrix
     mat <- as.matrix(nmds_input_filtered[, -1])
     rownames(mat) <- nmds_input_filtered$sampleID
     
     # Perform NMDS
     set.seed(123)
     nmds <- metaMDS(mat, distance = "bray", k = 2, trymax = 300)
     return(nmds)
}

atl_nmds <- getNMDS(data = atl)
gom_nmds <- getNMDS(data = gom)
npo_nmds <- getNMDS(data = npo)

# Extract out NMDS1 and NMDS2 for plotting
get_NMDSscores <- function(data = NULL, nmds_result = NULL){
     nmds_sites <- as.data.frame(nmds_result$points) %>% rownames_to_column(var = "sampleID")
     samp_meta <- data %>% select(sampleID, set, depth, depth_cat) %>% distinct
     data_merged <- left_join(nmds_sites, samp_meta, by = "sampleID")
     
     return(data_merged)   
}
dat_atl <- get_NMDSscores(data = atl, nmds_result = atl_nmds) %>% 
     left_join(., meta, by = "sampleID") %>% 
     mutate(date_parsed = ymd(collectionDate),
            doy = day(date_parsed)) %>% 
     filter(!MDS1 > 1000)  %>% 
     distinct
dat_atl %>% filter(MDS1 > 1000)

dat_gom <- get_NMDSscores(data = gom, nmds_result = gom_nmds) %>% 
     left_join(., gom %>% select(sampleID, longitude), by = "sampleID")  %>% 
     distinct

dat_npo <- get_NMDSscores(data = npo, nmds_result = npo_nmds) %>% 
     left_join(., npo %>% select(sampleID, latitude), by = "sampleID") %>% 
     distinct

# Plot out
plot_nmds <- function(data, color_var, color_label, color_option) {
     ggplot(data, aes(x = MDS1, y = MDS2)) +
          geom_point(aes(shape = depth_cat, fill = !!sym(color_var)), 
                     size = 3, color = "black", stroke = 0.8) +
          scale_fill_viridis_c(name = color_label, option = color_option, direction = -1) +
          scale_shape_manual(
               values = c("0-15m" = 21, "16-100m" = 24, "101-200m" = 23, ">200m" = 22),
               name = "Depth Category"
          ) +
          labs(title = "NMDS 'Bray-Curtis' Distance", x = "NMDS1", y = "NMDS2") 
}
p_atl <- plot_nmds(dat_atl, "doy", "Day of the Month", "cividis"); p_atl
p_gom <- plot_nmds(dat_gom, "longitude", "Longitude (W)", "magma"); p_gom
p_npo <- plot_nmds(dat_npo, "latitude", "Latitude (N)", "viridis"); p_npo

svg(paste0(sup_dir, "atl_nmds.svg"), width = 5.25, height = 3); print(p_atl); dev.off()
svg(paste0(sup_dir, "gom_nmds.svg"), width = 5.25, height = 3); print(p_gom); dev.off()
svg(paste0(sup_dir, "npo_nmds.svg"), width = 5.25, height = 3); print(p_npo); dev.off()

#################################################################################
# Permanova - NMDS    X
#################################################################################

atl_extract <- as.data.frame(atl_nmds$points) %>% rownames_to_column(var = "sampleID")
gom_extract <- as.data.frame(gom_nmds$points) %>% rownames_to_column(var = "sampleID")
npo_extract <- as.data.frame(npo_nmds$points) %>% rownames_to_column(var = "sampleID")

#meta <- read.csv("data/1_metaData/sample_info/ATL_GOM_sampleMeta.csv") %>% 
#     filter(collectionSet == "ATL") %>% select(sampleID, collectionDate)

atl_env <- dat_atl %>% select(sampleID, depth, collectionDate) %>%
     column_to_rownames("sampleID")
gom_env <- dat_gom %>% 
     column_to_rownames("sampleID")
npo_env <- dat_npo %>% column_to_rownames("sampleID") %>% drop_na()
npo = npo %>% drop_na(latitude)
getDistanceMatrix <- function(data) {
     metab_input <- data %>% 
          select(sampleID, Compound_Name, conc_nM_L) %>%
          distinct() %>%
          pivot_wider(names_from = Compound_Name, values_from = conc_nM_L) %>%
          mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
     
     sample_ids <- metab_input$sampleID
     compounds <- metab_input %>%
          select(-sampleID) %>%
          select(where(~ sum(.x) > 0))
     
     metab_matrix <- as.matrix(compounds)
     rownames(metab_matrix) <- sample_ids
     
     vegdist(metab_matrix, method = "bray")
}
atl_perNov <- adonis2(getDistanceMatrix(atl) ~ depth + collectionDate, data = atl_env, permutations = 999, by = "margin")
gom_perNov <- adonis2(getDistanceMatrix(gom) ~ depth + longitude, data = gom_env, permutations = 999, by = "margin")
npo_perNov <- adonis2(getDistanceMatrix(npo) ~ depth + latitude, data = npo_env, permutations = 999, by = "margin")

###############################################################################
# Spatial Look from Culture Data - Metabs of Interest
###############################################################################
sort(unique(NPO_metab_results$metab_data$Compound_Name))

interest <- c("Cysteinolate", "Cysteinesulfinate", "Homotaurine")

# Lisa Culture data, quick view of spatials 
x <- NPO_metab_results$metab_data %>% 
     filter(Compound_Name %in% interest) %>% 
     select(-type, -station, -Group) %>% 
     filter(depth <= 100) %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "0-15m",
                                  depth > 15 & depth <= 100 ~ "16-100m")) 

plot_latitude_trend <- function(data, compound) {
     p <- data %>% 
          filter(Compound_Name == compound) %>% 
          ggplot(aes(x = latitude, y = conc_nM_L)) +
          geom_point(size = 4, fill = "black") +
          geom_point(size = 3, aes(color = latitude)) +
          geom_smooth(method = "loess", se = FALSE, linetype = "dashed", color = "black", size = 0.75) +
          facet_wrap(~depth_cat, ncol = 1) +
          scale_color_viridis(option = "D", name = "Latitude (°N)") +
          labs(x = "Latitude (°N)", y = "Concentration (nM)") +
          theme(
               legend.position = "right",
               axis.text.x = element_text(angle = 45, hjust = 1)
          ) + facet_theme
     return(p)
}
p <- plot_latitude_trend(x, "Cysteinesulfinate"); p
svg(paste0(sup_dir, "Cysteinesulfinate_latitude.svg"), width = 5, height = 5); print(p); dev.off()
p <- plot_latitude_trend(x, "Cysteinolate"); p
svg(paste0(sup_dir, "Cysteinolate_latitude.svg"), width = 5, height = 5); print(p); dev.off()
p <- plot_latitude_trend(x, "Homotaurine"); p
svg(paste0(sup_dir, "Homotaurine_latitude.svg"), width = 5, height = 5); print(p); dev.off()

###############################################################################
# SAMPLE SIZE DEPTH TABLE...
###############################################################################

get_smps <- function(data, cruise) {
     data$metab_data %>% applyDepthCats() %>%
          select(sampleID, depth_cat) %>% distinct() %>%
          mutate(Cruise = cruise)
}
combined <- bind_rows(
     get_smps(ATL_metab_results, "ATL"),
     get_smps(GOM_metab_results, "GOM"),
     get_smps(NPO_metab_results, "NPO")
)

table <- combined %>%
     drop_na(depth_cat) %>% 
     group_by(Cruise, depth_cat) %>%
     summarise(count = n(), .groups = "drop") %>%
     pivot_wider(names_from = depth_cat, values_from = count, values_fill = 0) %>%
     mutate(Total = rowSums(across(where(is.numeric)))) %>% 
     select(Cruise, '0-15m',  '16-100m', '101-200m', '>200m', Total)
table

cruise_colors <- c("ATL" = "#50fa7b", "GOM" = "#bd93f9", "NPO" = "#ffb86c")

flex_table <- table %>%
     flextable() %>%
     set_header_labels(
          Cruise = "Ocean \n Region",
          `0-15m` = "Surface \n (0-15m)",
          `16-100m` = "Subsurface \n (16–100m)",
          `101-200m` = "Deep Euphotic \n (101–200m)",
          `>200m` = "Mesopelagic \n (>200m)",
          Total = "Total # \n (n)"
     ) %>%
     theme_vanilla() %>%
     autofit() %>%
     bg(i = ~ Cruise == "ATL", j = 1, bg = cruise_colors["ATL"]) %>%
     bg(i = ~ Cruise == "GOM", j = 1, bg = cruise_colors["GOM"]) %>%
     bg(i = ~ Cruise == "NPO", j = 1, bg = cruise_colors["NPO"]) %>%
     color(part = "header", color = "black") %>%
     bold(part = "header", bold = TRUE)
flex_table

save_as_image(
     flex_table,
     path = "results/supplementFigs/sampleSize_Table2.png",
     webshot = "webshot2",
     zoom = 2,
     expand = 1
)

# Get Sampling Sites
combined <- bind_rows(
     get_smps(ATL_metab_results, "ATL"),
     get_smps(GOM_metab_results, "GOM"),
     get_smps(NPO_metab_results, "NPO")
) %>% 
     mutate(sampleID2 = paste0(Cruise, "_", sampleID))

NPO_metab_results$metab_data %>%
     mutate(sampleGroup = gsub("^S(\\d+)E\\d+D(\\d+)_?[A-Z]?$", "S\\1D\\2", gsub("^Smp_", "", sampleID))) %>%
     distinct(sampleGroup) %>%
     nrow()
GOM_metab_results$metab_data %>%
     mutate(sampleGroup = gsub("^S(\\d+)E\\d+D(\\d+)_?[A-Z]?$", "S\\1D\\2", gsub("^Smp_", "", sampleID))) %>%
     distinct(sampleGroup) %>%
     nrow()
ATL_metab_results$metab_data %>%
     mutate(sampleGroup = gsub("^S(\\d+)E\\d+D(\\d+)_?[A-Z]?$", "S\\1D\\2", gsub("^Smp_", "", sampleID))) %>%
     distinct(sampleGroup) %>%
     nrow()

length(unique(ATL_metab_results$metab_data$depth)) 
length(unique(GOM_metab_results$metab_data$depth)) 
length(unique(NPO_metab_results$metab_data$depth)) 

min(ATL_metab_results$metab_data$depth) ; max(ATL_metab_results$metab_data$depth) 
min(GOM_metab_results$metab_data$depth) ; max(GOM_metab_results$metab_data$depth) 
min(NPO_metab_results$metab_data$depth) ; max(NPO_metab_results$metab_data$depth) 

###############################################################################
# END
###############################################################################




