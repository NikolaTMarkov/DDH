#######
#NTM 2025 
#Calculates correlation between different regions in fiso ficv odi or fa MD



 library ( tidyverse )
 library ( pheatmap ) 
 library ( gridExtra )
 library ( ggplotify )   # for saving grobs with ggsave
 
 load ( '/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/LongWideNODDI_DWI.RData' )

# Define earth palette reversed so that green is the low correlation and brown is high
earth_palette <- colorRampPalette (
  						rev ( c ("#543005", "#bf812d", "#f6e8c3", "#c7eae5", "#01665e" ) ) ) ( 100 )

# Helper: return a pheatmap object (not drawn) for one metric
corr_heatmap_obj <- function ( df , metric = "ficv" , main_prefix = "Pearson Correlations" ) {
					cols <- names ( df ) %>% keep ( ~ str_ends ( .x , paste0 ( "\\." , metric ) ) )
								if ( length ( cols ) < 2 ) stop ( "Not enough columns for metric: " , metric )

					mat <- df %>% select ( all_of ( cols ) )
  					corr_mat <- cor ( mat , use = "complete.obs" , method = "pearson" ) # "pairwise.complete.obs", "pearson" , "spearman" "kendall"

  # --- Clean labels: remove ".metric" suffix from col/row names ---
  clean_labels <- str_remove(cols, paste0("\\.", metric))
  rownames(corr_mat) <- clean_labels
  colnames(corr_mat) <- clean_labels

  pheatmap(
    corr_mat,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color = earth_palette,
    main = paste(main_prefix, "—", toupper(metric)),
    fontsize_row = 6, fontsize_col = 6,
    silent = TRUE
  )
}


# 1) Build objects
hm_ficv <- corr_heatmap_obj(noddi_wide, "ficv")
hm_fiso <- corr_heatmap_obj(noddi_wide, "fiso")
hm_odi  <- corr_heatmap_obj(noddi_wide, "odi")

# 2) Arrange in a grid for display
grid_arranged <- grid.arrange(hm_ficv$gtable, hm_fiso$gtable, hm_odi$gtable, ncol = 3)

# 3) Save each heatmap to its own file
ggsave("heatmap_ficv.png", as.ggplot(hm_ficv$gtable), width = 10, height = 8, dpi = 300)
ggsave("heatmap_fiso.png", as.ggplot(hm_fiso$gtable), width = 10, height = 8, dpi = 300)
ggsave("heatmap_odi.png",  as.ggplot(hm_odi$gtable),  width = 10, height = 8, dpi = 300)

# 4) Save the combined grid as well
ggsave("heatmaps_grid.png", as.ggplot(grid_arranged), width = 18, height = 7, dpi = 300)

################
# Repeat for DTI
hm_fa <- corr_heatmap_obj(dti_wide, "fa")
hm_md <- corr_heatmap_obj(dti_wide, "md")

grid_arranged <- grid.arrange(hm_fa$gtable, hm_md$gtable, ncol = 2)

ggsave("heatmaps_dti.png", as.ggplot(grid_arranged), width = 14, height = 7, dpi = 300)
