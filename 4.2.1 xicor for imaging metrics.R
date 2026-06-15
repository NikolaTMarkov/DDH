# NTM — 2025 August
# Matches DWI (NODDI/DTI) to cognition-ready matrices
# Replaces Pearson with Xi correlation (XICOR) for tract×tract heatmaps

──────────────────────────────────────────────────────────────────────

 library ( tidyverse )    
 library ( pheatmap )
 library ( gridExtra )
 library ( ggplotify )
 library ( XICOR )

load ( '/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/LongWideNODDI_DWI.RData' )

earth_palette <- colorRampPalette( rev( c("#543005", "#bf812d", "#f6e8c3", "#c7eae5", "#01665e")))(100)


# Minimal patient table (PIDN join helper)
patient <- UCSF_db$general %>%
  transmute(UnQID, PIDN, DCDate_record = as.Date(DCDate), Age = age_at_DCDate)

# Helper: attach PIDN by UnQID
attach_pidn <- function(df){
  df %>%
    left_join(select(patient, "UnQID", "PIDN"), by = "UnQID") %>%
    relocate(PIDN, UnQID, .before = 1)
}



# ──────────────────────────────────────────────────────────────────────────────
# 6) Xi correlation helpers
# ──────────────────────────────────────────────────────────────────────────────
# Pairwise Xi that handles NAs and small n
xi_pair <- function(x, y, min_n = 3){
  keep <- stats::complete.cases(x, y)
  if (sum(keep) < min_n) return(NA_real_)
  as.numeric(XICOR::xicor(x[keep], y[keep] ,  ties='average' ))  # add ties='average' if desired
}

# Build a symmetric Xi matrix for a numeric data.frame/matrix
xi_corr_mat <- function(mat){
  mat <- as.data.frame(mat)
  p <- ncol(mat)
  out <- matrix(NA_real_, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(p)){
    out[i, i] <- 1
    if (i < p){
      for (j in (i+1):p){
        val <- xi_pair(mat[[i]], mat[[j]])
        out[i, j] <- val
        out[j, i] <- val
      }
    }
  }
  out
}

# Xi-based heatmap factory
corr_heatmap_obj <- function(df, metric = "ficv", main_prefix = "Xi Correlations") {
  cols <- names(df) %>% keep(~ stringr::str_ends(.x, paste0("\\.", metric)))
  if (length(cols) < 2) stop("Not enough columns for metric: ", metric)

  mat <- df %>% dplyr::select(dplyr::all_of(cols))
  corr_mat <- xi_corr_mat(mat)

  # --- Clean labels: remove ".metric" suffix from col/row names ---
  clean_labels <- stringr::str_remove(cols, paste0("\\.", metric))
  rownames(corr_mat) <- clean_labels
  colnames(corr_mat) <- clean_labels

  pheatmap::pheatmap(
    corr_mat,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color = earth_palette,
    breaks = seq(0, 1, length.out = 101),   # Xi values are usually 0–1
    main = paste(main_prefix, "—", toupper(metric)),
    fontsize_row = 6, fontsize_col = 6,
    silent = TRUE
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# 7) Build heatmaps (Xi) and save
# ──────────────────────────────────────────────────────────────────────────────
hm_ficv <- corr_heatmap_obj(noddi_wide, "ficv")
hm_fiso <- corr_heatmap_obj(noddi_wide, "fiso")
hm_odi  <- corr_heatmap_obj(noddi_wide, "odi")

# Arrange a combined view for quick inspection
grid_arranged <- grid.arrange(hm_ficv$gtable, hm_fiso$gtable, hm_odi$gtable, ncol = 3)



#################
hm_fa <- corr_heatmap_obj(dti_wide, "fa")
hm_md <- corr_heatmap_obj(dti_wide, "md")

dev.new()
# Arrange a combined view for quick inspection
grid_arranged <- grid.arrange(hm_fa$gtable, hm_md$gtable, ncol = 2)


# Save individual and combined figures
#ggsave("heatmap_ficv_xi.png", as.ggplot(hm_ficv$gtable), width = 10, height = 8, dpi = 300)
#ggsave("heatmap_fiso_xi.png", as.ggplot(hm_fiso$gtable), width = 10, height = 8, dpi = 300)
#ggsave("heatmap_odi_xi.png",  as.ggplot(hm_odi$gtable),  width = 10, height = 8, dpi = 300)
#ggsave("heatmaps_grid_xi.png", as.ggplot(grid_arranged),  width = 18, height = 7, dpi = 300)

# ──────────────────────────────────────────────────────────────────────────────
# 8) Optional notes:
# - If your Xi should treat ties explicitly, use: XICOR::xicor(x, y, ties = "average")
#   inside xi_pair().
# - If you *know* Xi will be nonnegative in your context, you can switch to
#   breaks = seq(0, 1, length.out = 101) and use a 0→1 color palette.
# - To extend to FA/MD matrices, build fa_wide/md_wide similarly and pass to
#   corr_heatmap_obj(df = fa_wide, metric = "fa") etc.
# ──────────────────────────────────────────────────────────────────────────────
