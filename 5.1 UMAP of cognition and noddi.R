# Add UMAP library
library(umap)

# ───────────────────────────────────────────────────────────────
# UMAP Visualization of NODDI Data Colored by Age
# ───────────────────────────────────────────────────────────────

# Prepare NODDI data for UMAP (using the same imaging variables as in CCA)
noddi_umap_data <- NoddiCog_imputed %>%
  select(all_of(img_set)) %>%
  # Remove any remaining missing values
  drop_na() %>%
  # Scale the data for UMAP
  scale()

# Get the corresponding metadata for the filtered data
umap_metadata <- NoddiCog_imputed %>%
  select(UnQID, Age, sex, PIDN, DCDate) %>%
  # Filter to match the rows in noddi_umap_data
  slice(which(complete.cases(select(NoddiCog_imputed, all_of(img_set)))))

# Run UMAP
set.seed(123)  # For reproducibility
umap_result <- umap(noddi_umap_data, 
                   n_neighbors = 15,
                   min_dist = 0.5,
                   metric = "euclidean")

# Create UMAP dataframe
umap_df <- data.frame(
  UnQID = umap_metadata$UnQID,
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Age = umap_metadata$Age,
  sex = umap_metadata$sex,
  PIDN = umap_metadata$PIDN
) %>%
  # Join with CDR data if available
  left_join(patient %>% select(UnQID, cdr_global), by = "UnQID") %>%
  mutate(
    cdr_category = case_when(
      cdr_global == 0 ~ "CDR 0",
      cdr_global == 0.5 ~ "CDR 0.5", 
      cdr_global > 0.5 ~ "CDR >0.5",
      TRUE ~ "Missing"
    )
  )

# ───────────────────────────────────────────────────────────────
# Plot UMAP colored by Age
# ───────────────────────────────────────────────────────────────

# Basic UMAP colored by age
p_umap_age <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Age)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_gradientn(
    colours = c("#2c7bb6", "#ffff8c", "#d7191c"),
    name = "Age (years)"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP Visualization of NODDI Data",
    subtitle = "Colored by Age",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

print(p_umap_age)

# UMAP colored by sex for comparison
p_umap_sex <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sex)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = "tomato")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP Visualization of NODDI Data",
    subtitle = "Colored by Sex",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

print(p_umap_sex)

# UMAP colored by CDR category
p_umap_cdr <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cdr_category)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(
    values = c("CDR 0" = "#2c7bb6", "CDR 0.5" = "#ffff8c", "CDR >0.5" = "#d7191c", "Missing" = "gray"),
    name = "CDR Category"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP Visualization of NODDI Data",
    subtitle = "Colored by CDR Category",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

print(p_umap_cdr)

# ───────────────────────────────────────────────────────────────
# Enhanced UMAP with density contours
# ───────────────────────────────────────────────────────────────

p_umap_age_enhanced <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Age)) +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.3) +
  scale_color_gradientn(
    colours = c("#2c7bb6", "#ffff8c", "#d7191c"),
    name = "Age (years)"
  ) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#ffff8c", "#d7191c"),
    name = "Density"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP Visualization with Density Contours",
    subtitle = "NODDI data colored by Age",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

print(p_umap_age_enhanced)

# ───────────────────────────────────────────────────────────────
# Faceted UMAP by age groups
# ───────────────────────────────────────────────────────────────

umap_df <- umap_df %>%
  mutate(
    age_group = cut(Age, 
                   breaks = c(0, 65, 75, 85, 100),
                   labels = c("<65", "65-75", "76-85", ">85"))
  )

p_umap_faceted <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = age_group)) +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ age_group) +
  scale_color_viridis_d(name = "Age Group") +
  theme_minimal(base_size = 12) +
  labs(
    title = "UMAP Visualization by Age Groups",
    subtitle = "NODDI data faceted by age categories",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

print(p_umap_faceted)

# ───────────────────────────────────────────────────────────────
# Save UMAP plots
# ───────────────────────────────────────────────────────────────

# Combine main UMAP plots
umap_combined <- (p_umap_age | p_umap_sex) / p_umap_cdr +
  plot_annotation(
    title = "UMAP Analysis of NODDI White Matter Metrics",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

print(umap_combined)

ggsave("NTM2025_UMAP_NODDI_Analysis.pdf", umap_combined,
       width = 16, height = 12, units = "in", dpi = 300)

# Save individual plots
ggsave("NTM2025_UMAP_Age.pdf", p_umap_age, width = 10, height = 8, dpi = 300)
ggsave("NTM2025_UMAP_Age_Enhanced.pdf", p_umap_age_enhanced, width = 10, height = 8, dpi = 300)

# ───────────────────────────────────────────────────────────────
# Optional: UMAP parameter exploration
# ───────────────────────────────────────────────────────────────

# Function to try different UMAP parameters
explore_umap_parameters <- function(data, neighbors = c(5, 15, 30), min_dists = c(0.01, 0.1, 0.5)) {
  plots <- list()
  
  for (n in neighbors) {
    for (md in min_dists) {
      set.seed(123)
      umap_temp <- umap(data, n_neighbors = n, min_dist = md)
      
      temp_df <- data.frame(
        UMAP1 = umap_temp$layout[, 1],
        UMAP2 = umap_temp$layout[, 2],
        Age = umap_metadata$Age
      )
      
      p <- ggplot(temp_df, aes(x = UMAP1, y = UMAP2, color = Age)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_gradientn(colours = c("#2c7bb6", "#ffff8c", "#d7191c")) +
        theme_minimal() +
        labs(
          title = paste0("n_neighbors = ", n, ", min_dist = ", md),
          x = "UMAP1", y = "UMAP2"
        ) +
        theme(legend.position = "none")
      
      plots[[paste0("n", n, "_md", md)]] <- p
    }
  }
  
  return(plots)
}

# Uncomment to explore different UMAP parameters
# umap_param_plots <- explore_umap_parameters(noddi_umap_data)
# param_grid <- wrap_plots(umap_param_plots, ncol = 3) + 
#   plot_annotation(title = "UMAP Parameter Exploration")
# print(param_grid)