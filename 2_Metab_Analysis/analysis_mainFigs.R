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

df_res <- "out/mainFigs/"
fig_res <- "results/mainFigs/"

###############################################################################
# Group Totals - Regional - Main Figure
###############################################################################

datasets <- list(
     ATL = ATL_metab_results$metab_data_0.15,
     GOM = GOM_metab_results$metab_data_0.15,
     NPO = NPO_metab_results$metab_data_0.15
)

# Function to process
processDF <- function(metab_data, dataset_name) {
     totals <- metab_data %>%
          group_by(Compound_Name) %>%
          reframe(
               Compound_Name, Group, Sulfur_Group,
               med_cmpd_conc = median(conc_nM_L, na.rm = TRUE),
               mean_cmpd_conc = mean(conc_nM_L, na.rm = TRUE)
          ) %>% distinct()
     
     summary <- totals %>%
          group_by(Group) %>%
          summarise(
               sum = sum(mean_cmpd_conc, na.rm = TRUE),
               n_compounds = n()
          ) %>%
          arrange(desc(sum))
     
     summary %>%
          mutate(
               dataset = dataset_name,
               n_samples = length(unique(metab_data$sampleID))
          )
}
# Apply on surface samples across 3 dfs
data <- imap_dfr(datasets, ~ processDF(.x, .y))

# Filter shared groups
shared_groups <- data %>%
     group_by(Group) %>%
     summarise(n = n_distinct(dataset)) %>%
     filter(n == 3) %>%
     pull(Group)
filtered <- data %>% 
     filter(Group != "Internal Standard")

# Pull that for plot
npo_order <- filtered %>%
     filter(dataset == "NPO") %>%
     arrange(desc(sum)) %>%
     pull(Group)
gom_only <- filtered %>%
     filter(dataset == "GOM", !Group %in% npo_order) %>%
     arrange(desc(sum)) %>%
     pull(Group)
custom_order <- c(npo_order, gom_only)
filtered <- filtered %>%
     mutate(Group = factor(Group, levels = rev(custom_order)))

group_dataset_grid <- expand.grid(
     dataset = c("NPO", "GOM", "ATL"),
     Group = unique(filtered$Group),
     stringsAsFactors = FALSE
)

filtered_complete <- group_dataset_grid %>%
     left_join(filtered, by = c("dataset", "Group")) %>%
     mutate(
          sum = replace_na(sum, 0),
          n_compounds = replace_na(n_compounds, 0)
     )

# Calculate the total number of samples for each dataset
sample_counts <- filtered_complete %>%
     group_by(dataset) %>%
     summarise(n_samples = unique(n_samples)[1]) %>%
     mutate(facet_label = paste0(case_when(
          dataset == "NPO" ~ "North Pacific",
          dataset == "ATL" ~ "Atlantic",
          dataset == "GOM" ~ "Gulf of Mexico"
     ), "\n(n = ", n_samples, ")")
     )
# Create a named vector for the custom labeller
facet_labels <- setNames(sample_counts$facet_label, sample_counts$dataset)
# Plot out
filtered_complete <- filtered_complete %>% mutate(dataset = factor(dataset, levels = c("NPO", "GOM", "ATL")))
#p <- ggplot(filtered, aes(y = factor(Group, levels = rev(order)), 
p <- ggplot(filtered_complete, aes(y = Group, 
                          x = sum, fill = dataset)) +
     geom_bar(stat = "identity", position = "dodge", color = "black") + 
     geom_text(aes(label = n_compounds, x = sum + 0.1), 
               position = position_dodge(width = 0.9), hjust = 0) +
     
     scale_fill_manual(values = c("ATL" = "#50FA7B", "GOM" = "#BD93F9", "NPO" = "#FFB86C")) + 
     facet_wrap(~ dataset, ncol = 3, scales = "free",
                labeller = as_labeller(facet_labels)) +  # Use custom labels with sample counts
     labs(x = "Concentration (nM)",
          y = "Compound Group"
     )  +
     theme(
          legend.position = "none",  
          strip.text = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 14, face = "bold")
     ) + 
     facet_theme
p
svg(paste0(fig_res, "Cruise_BroadMetab_Means.svg"), width = 14, height = 5); p; dev.off()

#datasets$ATL %>% filter(Group == "Other") %>% group_by(Compound_Name) %>% reframe(sum = sum(conc_nM_L)) %>% arrange(sum)

###############################################################################
# Group Totals 2 - Regional - How many unique detecteds across all cruises?
###############################################################################

# What was detected?
cmpds <- 
     rbind(datasets$ATL %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           datasets$GOM %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           datasets$NPO %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name)
     ) %>% distinct

cmpds <- cmpds %>% left_join(., group_index, by = "Compound_Name") %>% select(-Sulfur_Group)

cmpds %>% 
     group_by(Group) %>% 
     reframe(count = n())

# Contributions of metab groups in each region
filtered %>% 
     group_by(dataset) %>% 
     mutate(inventory = sum(sum, na.rm = T),
            percent = (sum / inventory) *100) %>% 
     select(Group, dataset, percent) %>% 
     pivot_wider(., names_from = dataset, values_from = percent)
filtered %>% 
     group_by(dataset) %>% 
     select(Group, dataset, sum) %>% 
     pivot_wider(., names_from = dataset, values_from = sum)

###############################################################################
# Sulfur Totals - Regional - Main Figure
###############################################################################

# Surface sulfurs only
sulfur_data <- list(
     ATL = ATL_metab_results$metab_data_0.15 %>% filter(Group == "Sulfur"),
     GOM = GOM_metab_results$metab_data_0.15 %>% filter(Group == "Sulfur"),
     NPO = NPO_metab_results$metab_data_0.15 %>% filter(Group == "Sulfur")
)

getTotals <- function(data){
     data %>% 
          group_by(Compound_Name) %>% 
          reframe(Compound_Name, Group, Sulfur_Group,
                  mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
          distinct  
}
processDFs <- function(data, dataset_name) {
     n_samples <- data %>% pull(sampleID) %>% unique() %>% length()
     
     totals <- getTotals(data) %>%
          mutate(nSmps = n_samples)
     
     summary <- totals %>%
          group_by(Sulfur_Group) %>%
          summarise(
               sum = sum(mean_cmpd_conc, na.rm = TRUE),
               n_compounds = n()
          ) %>%
          arrange(desc(sum)) %>%
          mutate(
               dataset = dataset_name,
               n_samples = n_samples
          )
     
     return(summary)
}
sulfur_summary <- imap_dfr(sulfur_data, ~ processDFs (.x, .y))
# Calculate the total number of samples for each dataset
sample_counts <- sulfur_summary %>%
     group_by(dataset) %>%
     mutate(facet_label = paste0(case_when(
          dataset == "NPO" ~ "North Pacific",
          dataset == "ATL" ~ "Atlantic",
          dataset == "GOM" ~ "Gulf of Mexico"
     ), "\n(n = ", n_samples, ")"))

# Create a named vector for the custom labeller
facet_labels <- setNames(sample_counts$facet_label, sample_counts$dataset)
# Define the order based on NPO median concentrations
order <- sulfur_summary %>% filter(dataset == "NPO") %>%
     arrange(desc(sum)) %>%
     pull(Sulfur_Group) %>%  rev()
#order <- c("Sulfoxonium", "Sulfate ester", order )

#Plot
combined_summaries <- sulfur_summary %>%
     mutate(dataset = factor(dataset, levels = c("NPO", "GOM", "ATL")))
p <- ggplot(combined_summaries, aes(y = factor(Sulfur_Group, levels = order),
                                    x = sum, fill = dataset)) +
     geom_bar(stat = "identity", position = "dodge", color = "black") +  
     geom_text(aes(label = n_compounds, x = sum * 0.1),  # Move text to the left
               position = position_dodge(width = 0.9), hjust = 1.2,
               size = 5) +  # Adjust hjust to push left
     scale_fill_manual(values = c("ATL" = "#50FA7B", "GOM" = "#BD93F9", "NPO" = "#FFB86C")) + 
     facet_wrap(~ dataset, ncol = 1, scales = "free",
                labeller = as_labeller(facet_labels)) + 
     labs(
          y = "Concentration (nM)",
          x = "Sulfur Group",
          subtitle = "0-15m Depths"
     ) +
     theme(
          legend.position = "none",
          strip.text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 14, face = "bold")
     ) +
     facet_theme
p
# Save
svg(paste0(fig_res, "cruise_sulfurGroups_cumulatives.svg"), width = 5, height = 10); p; dev.off()

###############################################################################
# Sulfur Totals2 - Regional - How many unique detecteds across all cruises? How many Samples?
###############################################################################

cmpds <- 
     rbind(sulfur_data$ATL %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           sulfur_data$GOM %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name),
           sulfur_data$NPO %>% filter(!conc_nM_L <= 0) %>% select(Compound_Name)
     ) %>% distinct

cmpds <- cmpds %>% left_join(., group_index, by = "Compound_Name") %>% select(-Group)

cmpds %>% 
     group_by(Sulfur_Group) %>% 
     reframe(count = n())

sulfur_data$ATL %>% pull(sampleID) %>% unique %>% length
sulfur_data$GOM %>% pull(sampleID) %>% unique %>% length
sulfur_data$NPO %>% pull(sampleID) %>% unique %>% length

combined_summaries %>% group_by(dataset) %>% 
     mutate(total = sum(sum, na.rm = T),
            percent = sum / total) %>%
     select(Sulfur_Group, dataset, percent) %>% 
     pivot_wider(., names_from = Sulfur_Group, values_from = percent)
combined_summaries %>% group_by(dataset) %>% 
     select(Sulfur_Group, dataset, sum) %>% 
     pivot_wider(., names_from = Sulfur_Group, values_from = sum)

###############################################################################
# Sulfur Tops - Who are tops in the two abundant group categories? 
###############################################################################

groupsInterest <- c("Sulfonium", "Sulfonate")

sulfur_data <- list(
     ATL = ATL_metab_results$metab_data_0.15 %>% filter(Sulfur_Group %in% groupsInterest),
     GOM = GOM_metab_results$metab_data_0.15 %>% filter(Sulfur_Group %in% groupsInterest),
     NPO = NPO_metab_results$metab_data_0.15 %>% filter(Sulfur_Group %in% groupsInterest)
)

isolateTops <- function(data = NULL, top_n = NULL){
     data %>% 
          group_by(Compound_Name, Sulfur_Group) %>% 
          reframe(value = sum(conc_nM_L)) %>% 
          arrange(desc(value)) %>% 
          head(top_n)
     
}
list1 = isolateTops(data = sulfur_data$ATL, top_n = 8) %>% pull(Compound_Name) 
list2 = isolateTops(data = sulfur_data$GOM, top_n = 8) %>% pull(Compound_Name)
list3 = isolateTops(data = sulfur_data$NPO, top_n = 8) %>% pull(Compound_Name)

tops <- unique(c(list1, list2, list3))
tops

 ###############################################################################
# Sulfur Tops - Who are tops in the two abundant group categories? 
###############################################################################

getTotals <- function(data){
     data %>% 
          group_by(Compound_Name) %>% 
          reframe(Compound_Name, Group, Sulfur_Group,
                  mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
          distinct  
}


# Break Down most dominant
groupsInterest <- c("Sulfonium", "Sulfonate")
atl_dom = sulfur_data$ATL %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)
gom_dom = sulfur_data$GOM %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)
npo_dom = sulfur_data$NPO %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)


fetchPlot <- function(data = NULL){
     # Compute total mean concentration per Sulfur Group
     sulfur_totals <- data %>%
          group_by(Sulfur_Group) %>%
          summarize(total_conc = sum(mean_cmpd_conc, na.rm = TRUE))
     temp <- data %>%
          left_join(sulfur_totals, by = "Sulfur_Group") %>%
          mutate(percentage = (mean_cmpd_conc / total_conc) * 100)  # Compute % contribution
     # Group compounds making up <1% as "Other"
     temp <- temp %>%
          mutate(Compound_Name = ifelse(percentage < 1, "Other", Compound_Name)) %>%
          group_by(Sulfur_Group, Compound_Name) %>%
          summarize(mean_cmpd_conc = sum(mean_cmpd_conc), .groups = "drop")
     
     # Define separate color palettes for Sulfonium (warm) and Sulfonate (cool)
     pastel_palette <- c("#FFED6F", "#D48AC2", "#FFB347", "#A5B8FF",
                         "#FB8072", "#FCCDE5", "#6272A4", 
                         "#A8E6CF", "#B3DE69", "#A5B8FF")
     
     # Ensure we have enough colors for unique compounds
     num_compounds <- length(unique(temp$Compound_Name))
     
     if (num_compounds > length(pastel_palette)) {
          pastel_palette <- colorRampPalette(pastel_palette)(num_compounds)
     }
     # Assign black color to "Other" category
     compound_colors <- setNames(pastel_palette, unique(temp$Compound_Name))
     compound_colors["Other"] <- "#282A36"  # Black for "Other"
     # Reorder compounds within each Sulfur Group (largest to smallest)
     temp <- temp %>%
          group_by(Sulfur_Group) %>%
          mutate(Compound_Name = reorder(Compound_Name, -mean_cmpd_conc)) %>%
          ungroup()
     # Plot with stacked bars, ordered compounds, and designated colors
     ggplot(temp, aes(y = Sulfur_Group, x = mean_cmpd_conc, fill = Compound_Name)) +
          geom_bar(stat = "identity", position = "stack", color = "black") +
          scale_fill_manual(values = compound_colors) +
          theme(legend.position = "right")
}

p1 <- fetchPlot(data = npo_dom); p1
p2 <- fetchPlot(data = gom_dom); p2
p3 <- fetchPlot(data = atl_dom); p3

p <- plot_grid(p1, p2, p3, nrow = 3); p 
svg(paste0(fig_res, "cruise_sulfurGroups_doms_bars.svg"), width = 10, height = 10); p; dev.off()

###############################################################################
# Sulfur Tops - Percentages of the tops? 
###############################################################################
groupsInterest <- c("Sulfonium", "Sulfonate")
tops

getTotals <- function(data){
     data %>% 
          group_by(Compound_Name) %>% 
          reframe(Compound_Name, Group, Sulfur_Group,
                  mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
          distinct  
}


atl_dom = sulfur_data$ATL %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)
gom_dom = sulfur_data$GOM %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)
npo_dom = sulfur_data$NPO %>% getTotals %>% filter(Sulfur_Group %in% groupsInterest)

# Tops vs All
findPercents <- function(data = NULL){
     data %>%
          mutate(top = Compound_Name %in% tops) %>%
          group_by(Sulfur_Group, top) %>%
          summarise(total = sum(mean_cmpd_conc, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = top, values_from = total, names_prefix = "top_") %>%
          mutate(perc_top = 100 * top_TRUE / (top_TRUE + top_FALSE))
}
findPercents(data = atl_dom)
findPercents(data = gom_dom)
findPercents(data = npo_dom)

# Regional Composition
findDominants <- function(data = NULL){
     data %>%
          group_by(Sulfur_Group) %>%
          mutate(group_total = sum(mean_cmpd_conc, na.rm = TRUE)) %>%
          slice_max(mean_cmpd_conc, n = 2, with_ties = FALSE) %>%
          mutate(percent = 100 * mean_cmpd_conc / group_total) %>%
          select(Sulfur_Group, Compound_Name, percent, mean_cmpd_conc)   
}
findDominants(data = atl_dom)
findDominants(data = gom_dom)
findDominants(data = npo_dom)

# coefficient of variation (CV)
library(fmsb)
shortNames <- read.csv("data/1_metaData/std_libraries/sulfur_shortNames.csv")
atl_raw <- ATL_metab_results$metab_data %>%
     filter(Compound_Name %in% tops) %>%
     mutate(region = "ATL") %>% 
     select(Compound_Name, sampleID, latitude, depth, region, conc_nM_L)

gom_raw <- GOM_metab_results$metab_data %>%
     filter(Compound_Name %in% tops) %>%
     mutate(region = "GOM") %>% 
     select(Compound_Name, sampleID, latitude, depth, region, conc_nM_L)

npo_raw <- NPO_metab_results$metab_data %>%
     filter(Compound_Name %in% tops) %>%
     mutate(region = "NPO") %>% 
     select(Compound_Name, sampleID, latitude, depth, region, conc_nM_L)

sulf_all <- bind_rows(atl_raw, gom_raw, npo_raw) 

cv_table <- sulf_all %>%
     mutate(depth_bin = ifelse(depth <= 15, "Surf", "SubSurf")) %>%
     group_by(region, depth_bin, Compound_Name) %>%
     summarise(mean_nM = mean(conc_nM_L, na.rm = TRUE),
               sd_nM = sd(conc_nM_L, na.rm = TRUE),
               CV = sd_nM / mean_nM,
               .groups = "drop")

surf <- cv_table %>% filter(depth_bin == "Surf")
subsurf <- cv_table %>% filter(depth_bin == "SubSurf")

# Filter for depth group
table = surf
cv_subset <- table %>%
     select(region, Compound_Name, CV) %>%
     pivot_wider(names_from = Compound_Name, values_from = CV) %>%
     mutate(` ` = NA) %>%
     relocate(` `, .after = region) %>% 
     filter(region == "NPO")

# Make radar matrix
cv_radar <- cv_subset %>%
     column_to_rownames("region")
cv_radar <- rbind(
     rep(2, ncol(cv_radar)),  # max
     rep(0, ncol(cv_radar)),  # min
     cv_radar
)
radarchart(cv_radar,
           pcol = c("forestgreen", "darkorange", "steelblue"),
           pfcol = c("#228B2210", "#FF8C0010", "#4682B410"),
           plwd = 2,
           lty = 1,
           cglcol = "grey", 
           cglty = 3,
           axislabcol = "black",
           axistype = 1,
           caxislabels = c("0", "0.5", "1.0", "1.5", "2.0"),
           vlcex = 0.7,
           title = paste("CV Profiles of Sulfur Compounds - 0 to 15m. NPO"#, depth_group
                         ))

legend("topright", legend = rownames(cv_radar)[3:5],
       col = c("forestgreen", "darkorange", "steelblue"),
       lty = 1, lwd = 2, bty = "n")

# -------------------------------------------------------------

npo_raw <- NPO_metab_results$metab_data %>%
     filter(Compound_Name %in% tops) %>%
     mutate(region = "NPO") %>% 
     select(Compound_Name, sampleID, latitude, depth, region, conc_nM_L) %>% 
     mutate(NPO_zone = case_when(
          latitude < 32.5 ~ "NPSG",
          latitude >= 32.5 & latitude < 36.2 ~ "STZ",
          latitude >= 36.2 ~ "NTZ"
     ))

cv_latitudinal <- npo_raw %>%
     filter(depth <= 15) %>% 
     group_by(NPO_zone, Compound_Name) %>%
     summarise(mean_nM = mean(conc_nM_L, na.rm = TRUE),
               sd_nM = sd(conc_nM_L, na.rm = TRUE),
               CV = sd_nM / mean_nM,
               .groups = "drop")
plot_radar <- function(df) {
     mat <- df %>%
          select(NPO_zone, Compound_Name, CV) %>%
          pivot_wider(names_from = Compound_Name, values_from = CV) %>%
          column_to_rownames("NPO_zone")
     
     mat <- rbind(rep(2, ncol(mat)), rep(0, ncol(mat)), mat)
     
     radarchart(mat,
                pcol = c("forestgreen", "darkorange", "steelblue"),
                pfcol = c("#228B2210", "#FF8C0010", "#4682B410"),
                plwd = 2, lty = 1,
                cglcol = "grey", cglty = 3,
                axislabcol = "black", axistype = 1,
                caxislabels = c("0", "0.5", "1.0", "1.5", "2.0"),
                vlcex = 0.7,
                title = paste("NPO - Surface CV by Latitudinal Zone"))
     
     legend("topright", legend = rownames(mat)[3:5],
            col = c("forestgreen", "darkorange", "steelblue"),
            lty = 1, lwd = 2, bty = "n")
}
plot_radar(cv_latitudinal)

# -------------------------------------------------------------
cv_table <- sulf_all %>%
     drop_na(depth) %>%
     mutate(depth_bin = case_when(
          depth <= 15 ~ "D1",
          depth > 15 & depth <= 100 ~ "D2",
          depth > 100 & depth <= 200 ~ "D3",
          depth > 200 ~ "D4"
     )) %>%
     group_by(region, Compound_Name, depth_bin) %>%
     summarise(mean_nM = mean(conc_nM_L, na.rm = TRUE),
               sd_nM = sd(conc_nM_L, na.rm = TRUE),
               CV = sd_nM / mean_nM,
               .groups = "drop") %>%
     pivot_wider(names_from = depth_bin, values_from = CV)

region_radar <- function(region_input) {
     df <- cv_table %>%
          filter(region == region_input) %>%
          group_by(Compound_Name) %>%
          summarise(across(D1:D4, ~mean(.x, na.rm = TRUE))) %>%
          column_to_rownames("Compound_Name") %>%
          t() %>%
          as.data.frame()
     
     df <- df[rowSums(!is.na(df)) > 0, ]
     
     print("Filtered df (after dropping NA rows and columns):")
     print(df)
     
     if (nrow(df) < 2 || ncol(df) < 2) {
          warning("Not enough data to plot radar chart for ", region_input)
          return(NULL)
     }
     
     df <- rbind(rep(2, ncol(df)), rep(0, ncol(df)), df)
     
     radarchart(df,
                pcol = rep_len(c("skyblue", "tomato", "goldenrod", "purple"), nrow(df)-2),
                pfcol = rep_len(c("#87CEEB30", "#FF634730", "#DAA52030", "#80008030"), nrow(df)-2),
                plwd = 2,
                lty = 1,
                axistype = 1,
                caxislabels = c("0", "0.5", "1.0", "1.5", "2.0"),
                cglcol = "grey", 
                cglty = 3,
                vlcex = 0.7,
                title = paste("CV Profiles by Depth -", region_input)
     )
     
     legend("topright", legend = rownames(df)[3:nrow(df)],
            col = rep_len(c("skyblue", "tomato", "goldenrod", "purple"), nrow(df)-2),
            lty = 1, lwd = 2, bty = "n")
}


region_radar("NPO")
region_radar("GOM")
region_radar("ATL")


###############################################################################
# NPO Latitude Dots, Dominant Sulfurs
###############################################################################

x <- NPO_metab_results$metab_data_0.15 %>% filter(Compound_Name %in% tops) %>% 
     mutate(latitude = round(latitude)) %>% 
     group_by(Compound_Name, latitude) %>% 
     reframe(Compound_Name, Sulfur_Group,
             mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
     ungroup %>% distinct %>% 
     complete(Compound_Name, fill = list(mean_cmpd_conc = NA)) %>% 
     mutate(mean_cmpd_conc = replace_na(mean_cmpd_conc, 0)) %>% 
     drop_na(latitude)
sort(unique(x$Compound_Name))
x %>% filter(is.na(Compound_Name))
cmpd_order <- c("Sulfolactate", "DHPS", "Taurine",
                "Isethionate", "Cysteinolate", "DMS-Ac",
                "DMSP", "Gonyol")
salinityFraont = 32.5
chlaFront = 36.2

p <- ggplot(x, aes(x = factor(Compound_Name, levels = cmpd_order), 
                   y = latitude, size = mean_cmpd_conc, 
                   fill = latitude)) +
     geom_point(shape = 21, color = "black", stroke = 0.8) +  
     scale_size(range = c(2, 10), name = "Mean Concentration (nM)") +  
     scale_fill_viridis(option = "D", name = "Mean Concentration (nM)") +
     labs(y = "Latitude (N)") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
           axis.text.y = element_text(size = 15),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size = 20)) +
     geom_hline(yintercept = salinityFraont, linetype = "dashed") +
     geom_hline(yintercept = chlaFront, linetype = "solid")
p
svg(paste0(fig_res, "npo_latitudeDots.svg"), width = 8, height = 4.5); p; dev.off()

# Reporting Scheme in Text Document
x <- NPO_metab_results$metab_data %>% filter(Compound_Name %in% tops) %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "D1",
                                  depth > 15 & depth <= 100 ~ "D2",
                                  depth > 100 & depth <= 200 ~ "D3",
                                  depth > 200 ~ "D4")) %>% 
     mutate(latitude = round(latitude)) %>% 
     group_by(Compound_Name, depth_cat, latitude) %>% 
     reframe(Compound_Name, Sulfur_Group, depth_cat,
             mean_cmpd_conc = mean(conc_nM_L, na.rm = T)) %>% 
     ungroup %>% distinct %>% 
     complete(Compound_Name, fill = list(mean_cmpd_conc = NA)) %>% 
     mutate(mean_cmpd_conc = replace_na(mean_cmpd_conc, 0)) %>% 
     drop_na(latitude)
unique(x$latitude)

dat = x %>% pivot_wider(., values_from = mean_cmpd_conc, names_from = latitude)

dat_long <- dat %>%
     pivot_longer(
          cols = `26`:`41`,
          names_to = "Latitude",
          values_to = "Concentration"
     ) %>%
     mutate(
          Latitude = as.numeric(Latitude)
     )
p <- ggplot(dat_long, aes(x = as.character(Latitude),
                          y = Concentration, 
                          color = depth_cat,
                          group = interaction(Compound_Name, depth_cat))) +
     geom_line() +
     geom_point() +
     facet_wrap(~ Compound_Name, scales = "free_y") +
     labs(x = "Latitude (°N)",
          y = "Concentration (nM)",
          color = "Depth Group") +
     theme_bw() +facet_theme
p
svg(paste0(fig_res, "npo_latitudeTracking2.svg"), width = 8, height = 4.5); p; dev.off()

###############################################################################
# Sulfur Totals of Top 8 - Spatials NPO, GOM - Main Figure, Supplement Figure
###############################################################################

# Surface and SubSurface depths
npo = NPO_metab_results$metab_data %>% 
     filter(Compound_Name %in% tops) %>% 
     select(-Group, -Sulfur_Group, -station, -type, -depth_m)
npo %>% pull(Compound_Name) %>%  unique %>%  sort

gom = GOM_metab_results$metab_data %>% filter(Compound_Name %in% tops) %>% 
     select(-Group, -Sulfur_Group, -station, -type, -depth_m)
gom %>% pull(Compound_Name) %>%  unique %>%  sort

# Depth look, geomTile
npo <- npo %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                  depth > 15 & depth <= 100 ~ "D2_16-100m")
     ) %>% 
     group_by(Compound_Name, depth_cat) %>% 
     reframe(value = mean(conc_nM_L, na.rm = T))
gom <- gom %>% 
     mutate(depth_cat = case_when(depth <= 15 ~ "D1_0-15m",
                                  depth > 15 & depth <= 100 ~ "D2_16-100m",
                                  depth > 100 & depth <= 200 ~ "D3_101-200m",
                                  depth > 200 ~ "D4_>200m")
     ) %>% 
     group_by(Compound_Name, depth_cat) %>% 
     reframe(value = mean(conc_nM_L, na.rm = T))

# Plot
max_val <- 5.5
min_val <- 0
color_palette <- c("#282A36", "#4299E1", "#72C6DA", "#A8E6CF", "#D4F1BE", "#F4E285", "#FAEDB8")
depth_order <- npo %>% pull(depth_cat) %>% unique
cmpd_order <- c("Sulfolactate", "DHPS", "Taurine", "Isethionate", "Cysteinolate", "DMS-Ac", "DMSP", "Gonyol")
npo2 <- npo %>% drop_na(depth_cat)

p <- ggplot(npo2, aes(x = factor(Compound_Name, levels = cmpd_order), 
                     y = factor(depth_cat, levels = rev(depth_order)), 
                     fill = value)) +
     geom_tile(color = "black", size = .5) +
     geom_text(aes(label = round(value, 2)), color = "black", size = 5) +
     scale_fill_gradientn(colors = color_palette, limits = c(min_val, max_val)) + 
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           panel.grid = element_blank(),
           axis.title = element_blank())
p
suppDir = "results/supplementFigs/"
svg(paste0(suppDir, "npo_mtbs_depths.svg"), width = 6, height = 3); p; dev.off()

#
max_val <- 1.5
min_val <- 0
color_palette <- c("#282A36", "#4299E1", "#72C6DA", "#A8E6CF", "#D4F1BE", "#F4E285", "#FAEDB8")
depth_order <- gom %>% pull(depth_cat) %>% unique
p <- ggplot(gom, aes(x = factor(Compound_Name, levels = cmpd_order), 
                     y = factor(depth_cat, levels = rev(depth_order)), 
                     fill = value)) +
     geom_tile(color = "black", size = .5) +
     geom_text(aes(label = round(value, 3)), color = "black", size = 5) +
     scale_fill_gradientn(colors = color_palette, limits = c(min_val, max_val)) + 
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           panel.grid = element_blank(),
           axis.title = element_blank())
p
svg(paste0(suppDir, "gom_mtbs_depths.svg"), width = 7, height = 4); p; dev.off()


# Find Fold change
gom %>%
     filter(depth_cat %in% c('D1_0-15m', 'D2_16-100m')) %>% 
     group_by(Compound_Name) %>%
     summarise(
          max_value = max(value, na.rm = TRUE),
          min_value = min(value, na.rm = TRUE),
          fold_change = max_value / min_value,
          range = max_value - min_value
     ) %>%
     arrange(desc(fold_change))

npo %>%
     filter(depth_cat %in% c('D1_0-15m', 'D2_16-100m')) %>% 
     group_by(Compound_Name) %>%
     summarise(
          max_value = max(value, na.rm = TRUE),
          min_value = min(value, na.rm = TRUE),
          fold_change = max_value / min_value,
          range = max_value - min_value
     ) %>%
     arrange(desc(fold_change))


#################################################################################
# LISA'S Culture work
#################################################################################

culture <- read.csv("data/CultureWork/rpomTests_RK.csv") %>% tibble
anionNames <- read.csv("data/1_metaData/nameKey_Anion.csv")
shortNames <- read.csv("data/1_metaData/std_libraries/sulfur_shortNames.csv")

# Clean names
data <- culture %>% 
     rename(Compound_Name = COMPOUND) %>% 
     left_join(., anionNames, by = "Compound_Name") %>% 
     mutate(Compound_Name = coalesce(New_Name, Compound_Name)) %>% 
     left_join(., shortNames, by = "Compound_Name") %>% 
     mutate(Compound_Name = coalesce(shortName, Compound_Name)) %>% 
     select(-New_Name, -shortName)

     
data_summary <- data %>%
     group_by(Compound_Name, DAYS) %>%
     summarise(mean_OD = mean(OD, na.rm = TRUE),
               sd_OD = sd(OD, na.rm = TRUE), .groups = "drop") %>% 
     filter(Compound_Name != "Ergothione")

unique(data$Compound_Name)
compound_colors <- c(
     "Cysteinolate" = "#D48AC2",  # Specific color
     "No Carbon Source" = "#000000"  # Black for No Carbon
)
unique_compounds <- setdiff(unique(data_summary$Compound_Name), names(compound_colors))
other_colors <- viridis::viridis(length(unique_compounds), option = "D", begin = 0.2, end = 0.8)
final_colors <- setNames(c(compound_colors, other_colors), c(names(compound_colors), unique_compounds))


p <- ggplot(data_summary, aes(x = DAYS, y = mean_OD, color = Compound_Name, group = Compound_Name)) +
     geom_line(size = 1) +  
     geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), width = 0.2, position = position_dodge(width = .2)) + 
     geom_point(size = 3, alpha = 0.9) +
     scale_color_manual(values = final_colors) +  
     labs(x = "Days", y = "Optical Density (OD)", title = "OD Changes Over Time for Compounds") +
     theme(
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold")
     )
p
svg(paste0(fig_res, "sulfurSelects_culture.svg"), width = 6, height = 3.5); p; dev.off()




















