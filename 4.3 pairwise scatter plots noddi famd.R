

 load ( '/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/LongWideNODDI_DWI.RData' )
# --- List available regions (cleaned names) for a metric ----------------------
list_metric_regions <- function(df, metric = "ficv") {
  cols <- names(df)[stringr::str_ends(names(df), paste0("\\.", metric))]
  sort(stringr::str_remove(cols, paste0("\\.", metric)))
}

# --- Base-R pairs() helper (fixed + regex region matching) --------------------
pair_scatter_matrix_base <- function(df,
                                     metric = "ficv",
                                     regions = NULL,      # exact names OR regex patterns (case-insensitive)
                                     n_max = 12,          # set NULL for all
                                     file = NULL,         # e.g., "pairs_ficv.png"
                                     width = NULL, height = NULL,  # pixels (PNG)
                                     point_alpha = 0.25,
                                     point_cex = 0.5,
                                     main_prefix = "Pairs — ") {
  stopifnot(is.data.frame(df))

  # 1) columns for this metric
  cols_all <- names(df)[stringr::str_ends(names(df), paste0("\\.", metric))]
  if (length(cols_all) < 2) stop("Not enough columns for metric: ", metric)

  clean_all <- stringr::str_remove(cols_all, paste0("\\.", metric))

  # 2) optional selection: exact names or regex/patterns
  cols <- cols_all
  if (!is.null(regions)) {
    # expand patterns to matching clean names
    idx <- logical(length(clean_all))
    for (pat in regions) {
      idx <- idx | grepl(pat, clean_all, ignore.case = TRUE)
    }
    cols <- cols_all[idx]
    if (length(cols) < 2) {
      stop("Fewer than 2 matching regions for metric: ", metric,
           "\nRequested patterns: ", paste(regions, collapse = ", "),
           "\nAvailable (use list_metric_regions(df, '", metric, "'))")
    }
  }

  # 3) keep top-variance regions if requested
  if (!is.null(n_max) && length(cols) > n_max) {
    vars <- vapply(cols, function(z) stats::var(df[[z]], na.rm = TRUE), numeric(1))
    cols <- cols[order(vars, decreasing = TRUE)][seq_len(n_max)]
  }

  # 4) working data + clean labels
  dat <- df[, cols, drop = FALSE]
  clean_labels <- stringr::str_remove(cols, paste0("\\.", metric))
  # ensure numeric
  dat[] <- lapply(dat, function(v) suppressWarnings(as.numeric(v)))
  dat <- stats::na.omit(dat)
  colnames(dat) <- clean_labels

  # --- custom panels (no '...' forwarded to avoid pch duplication) ------------
  panel_lower_points <- function(x, y) {
    points(x, y, pch = 16, cex = point_cex, col = rgb(0, 0, 0, alpha = point_alpha))
  }

  panel_upper_cor <- function(x, y, digits = 2, cex.cor = 1.2) {
    oldusr <- par("usr"); on.exit(par(usr = oldusr), add = TRUE)
    par(usr = c(0, 1, 0, 1))
    r <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
    txt <- if (is.finite(r)) formatC(r, digits = digits, format = "f") else "NA"
    text(0.5, 0.5, txt, cex = cex.cor)
  }

  panel_diag_hist <- function(x) {
    oldusr <- par("usr"); on.exit(par(usr = oldusr), add = TRUE)
    par(usr = c(oldusr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE, breaks = "FD")
    if (length(h$counts)) {
      y <- h$counts / max(h$counts, na.rm = TRUE)
      rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y, col = "grey80", border = "white")
    }
  }

  title_txt <- paste0(main_prefix, toupper(metric),
                      " (", ncol(dat), " regions; n=", nrow(dat), ")")

  # 5) optional file output (PNG)
  if (!is.null(file)) {
    p <- ncol(dat)
    if (is.null(width))  width  <- max(1600, 140 * p)
    if (is.null(height)) height <- width
    png(filename = file, width = width, height = height, res = 300)
    on.exit(dev.off(), add = TRUE)
  }

  # 6) draw (NOTE: no extra args that would be forwarded to panels)
  pairs(dat,
        lower.panel = panel_lower_points,
        upper.panel = panel_upper_cor,
        diag.panel  = panel_diag_hist,
        main = title_txt)
}





###############
pair_scatter_matrix_base(noddi_wide, metric = "ficv", n_max = 12,
                         file = "pairs_ficv.png")
pair_scatter_matrix_base(noddi_wide, metric = "fiso", n_max = 12,
                         file = "pairs_fiso.png")
pair_scatter_matrix_base(noddi_wide, metric = "odi",  n_max = 12,
                         file = "pairs_odi.png")

pair_scatter_matrix_base(dti_wide, metric = "fa", n_max = 12,
                         file = "pairs_fa.png")
pair_scatter_matrix_base(dti_wide, metric = "md", n_max = 12,
                         file = "pairs_md.png")


pair_scatter_matrix_base(
  noddi_wide,
  metric  = "ficv",
  regions = c( "fusiform"    ,
                "superiortemporal" ,
              "supramarginal"),
  n_max   = NULL,
  file    = "pairs_ficv_selected.png"
)
