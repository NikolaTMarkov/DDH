###############################################################################
#NTM Aug 2025 
# NODDI + cognition + outliers + sex-aware analyses
# fits Robust linear model to NODDI values to explore effect of age
###############################################################################

 library ( tidyverse )
 library ( broom )

# ---- load your data ----
# load("/mnt/data/UCSF_db_Branch_Summer2024.RData")
 load ( "~/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/UCSF_db_Branch_Summer2024.RData" )

# ---- switches ----
 METRICS_TO_USE  <- c ( "ficv" , "fiso" , "odi" )



 # ---- patient + CDR ----
 patient <- UCSF_db$general %>%
						transmute (
								PIDN , UnQID,
								DCDate = as.Date ( DCDate ) ,
								Age    = age_at_DCDate ) %>%
						left_join (
								UCSF_db$demographics %>%
								transmute (
								UnQID ,
								sex = case_when (
          									gender == 1 ~ "Male" ,
          									gender == 2 ~ "Female" ,
          									TRUE ~ NA_character_ ) ) ,
    							by = "UnQID" ) %>%
  						left_join ( UCSF_db$cdr %>%
							transmute (
        						UnQID ,
        						cdr_global = as.numeric ( cdr_global ) ,  # be explicit
        						CDR_class = factor (
          									case_when (
            									is.na ( cdr_global )	~ "Unknown" ,
            									cdr_global == 0			~ "CU" ,
            									cdr_global == 0.5		~ "MCI"	,
            									cdr_global  > 0.5		~ "SCI" ,
            									TRUE					~ "Unknown" ) ,
          									levels = c ( "Unknown" , "CU" , "MCI" , "SCI" ) ) ) ,
    						by = "UnQID" ) %>%
  						mutate ( sex = factor ( sex , levels = c ( "Male","Female" ) ) )

# Colors for plots
 palette_cdr <- c ( "Unknown" = "gray70" ,
                 "CU"      = "chartreuse4" ,
                 "MCI"     = "orange" ,
                 "SCI"     = "red" )



# ---- helper: longify a NODDI table (expects UnQID, tiv, parcel cols) ----
longify_noddi <- function ( df , metric_name ) {
 				 # 1) Join patient, keep DCDate from patient only
  				d <- df %>%
    					dplyr::select ( -DCDate ) %>% # ignore DCDate from noddi
    					left_join ( patient , by = "UnQID" ) %>%
    					dplyr::select (
      							PIDN , UnQID , DCDate , Age , sex , CDR_class ,
      							#starts_with ( "global" ) ,
      							starts_with ( "wm_lh_" ) , starts_with ( "wm_rh_" ) )

  # 2) Coerce parcel columns to numeric and keep rows with any NODDI value present
  parcel_cols <- setdiff( 
  						names ( d ) , 
  						c ( "PIDN" , "UnQID" , "DCDate" , "Age" , "sex" , "CDR_class" ) ) #returns everything else than c
  d <- d %>%
    		mutate ( across ( all_of ( parcel_cols ) , ~ suppressWarnings ( as.numeric ( .x ) ) ) ) %>%
    		filter ( if_any ( all_of ( parcel_cols ) , ~ !is.na ( .x ) ) )

  # 3) Now choose earliest visit **among rows with NODDI data**
  d <- d %>%
    		group_by ( PIDN ) %>%
    		{ if ( any ( !is.na ( .$DCDate ) ) ) slice_min ( . , DCDate , with_ties = FALSE ) else slice_head ( . , 1 ) } %>%
    		ungroup ()

  # 4) Long format, collapse L/R, average duplicates, tag metric
  d %>%
    	pivot_longer (
      		cols = all_of ( parcel_cols ) ,
			names_to = "parcel", values_to = "value"  ) %>%
    		mutate ( parcel = stringr :: str_remove ( parcel , "^wm_(lh|rh)_" ) ) %>%
    		group_by ( PIDN , UnQID , Age , sex , CDR_class , parcel ) %>% 
    		summarise ( value = mean ( value , na.rm = TRUE ) , .groups = "drop" ) %>%
    		filter ( !is.na ( value ) & !is.nan ( value ) ) %>%   # drop empty cells
    		mutate ( metric = metric_name )
 }#End of longify noddi function



# ---- build long data ----
noddi_sources <- list (
		ficv = UCSF_db$imaging_noddi_ficv ,
  		fiso = UCSF_db$imaging_noddi_fiso ,
  		odi  = UCSF_db$imaging_noddi_odi )

noddi_long <- purrr::imap_dfr ( noddi_sources [ METRICS_TO_USE ] , longify_noddi )

# ---- outliers: loess residual ± 3*IQR within (metric, parcel) ----
flag_outliers <- function ( d ) {
  					d %>%
						group_by ( metric, parcel ) %>%
    					arrange ( Age ) %>%
    						group_modify ( ~ {
      							x <- .x
      							ok <- !is.na( x$value ) & !is.na ( x$Age )
      					pred <- rep ( NA_real_ , nrow ( x ) )
      						if ( sum ( ok ) >= 8 ) {
        						pred [ ok ] <- predict ( loess ( value ~ Age , data = x [ ok , ] , span = 0.75 , na.action = na.exclude ) ,
                            	newdata = x [ ok , ] ) }
      					x$loess_fit <- pred
      					x } ) %>%
    					mutate (
      						resid  = value - loess_fit ,
      						med    = median ( resid , na.rm = TRUE ) ,
      						iqr    = IQR ( resid , na.rm = TRUE ) ,
      						out_hi = med + 3*iqr ,
      						out_lo = med - 3*iqr ,
      						is_outlier = resid > out_hi | resid < out_lo ) %>%
    					ungroup () }#end outlier flag

 noddi_flagged <- flag_outliers ( noddi_long )
 bad_ids <- noddi_flagged %>% filter ( is_outlier ) %>% distinct ( UnQID ) %>% pull ( )
 noddi_clean <- noddi_long %>% filter( !UnQID %in% bad_ids )



# ---- plots: one per metric, sex as rows, BIG RED STARS if GLM is significant ----
plot_metric_by_sex <- function(metric_name, sexes = c("Male","Female")) {
  d <- dplyr::filter(noddi_clean, metric == metric_name)
  if (nrow(d) == 0) return(invisible(NULL))

  for (s in sexes) {
    ds <- dplyr::filter(d, sex == s)
    if (nrow(ds) == 0) next

    # Per-parcel GLM and star code
    star_df <- ds %>%
      dplyr::group_by(parcel) %>%
      dplyr::group_modify(~{
        x <- dplyr::filter(.x, !is.na(Age) & !is.na(value))
        if (nrow(x) < 3) return(tibble::tibble(star = ""))
        fit <- try(glm(value ~ Age, data = x, family = gaussian()), silent = TRUE)
        if (inherits(fit, "try-error")) return(tibble::tibble(star = ""))
        sm <- summary(fit)$coefficients
        p  <- if ("Age" %in% rownames(sm)) sm["Age","Pr(>|t|)"] else NA_real_
        star <- dplyr::case_when(
          is.na(p)       ~ "",
          p < 0.001      ~ "***",
          p < 0.01       ~ "**",
          p < 0.05       ~ "*",
          TRUE           ~ ""
        )
        tibble::tibble(star = star)
      }) %>%
      dplyr::ungroup()

    ds_lab <- ds %>% dplyr::left_join(star_df, by = "parcel")

    p <- ggplot(ds_lab, aes(x = Age, y = value)) +
      geom_point(alpha = 0.35, size = 0.8, colour = "black") +
      geom_smooth(method = "glm", method.args = list(family = gaussian()),
                  se = TRUE, linewidth = 0.9, colour = "black") +
      facet_wrap(~ parcel, scales = "free_y") +
      # Big red stars in the top-right of each facet (only if star != "")
      geom_text(
        data = dplyr::filter(ds_lab, star != ""),
        aes(x = Inf, y = Inf, label = star),
        inherit.aes = FALSE,
        hjust = 1.1, vjust = 1.1,  # nudge inside the panel
        colour = "red", size = 8    # BIG red asterisks
      ) +
      labs(
        title = paste0(toupper(metric_name), " vs Age (outliers removed) — ", s),
        x = "Age (years)", y = "Value"
      ) +
      theme_bw(9) +
      theme(strip.text = element_text(size = 7),
            axis.text  = element_text(size = 7))

    dev.new(); print(p)
  }
  invisible(NULL)
}

# run for each metric
invisible(lapply(METRICS_TO_USE, plot_metric_by_sex))

########################
###############################################################################
# RLM per tract: value ~ Age + sex  (for NODDI metrics)
# Output: forest-style plots of Age beta (±95% CI) by parcel, faceted by metric
###############################################################################

library(tidyverse)
library(MASS)      # rlm
library(broom)     # tidy summaries

# -------------------------------------------------------------------
# ASSUMPTIONS:
# You already built `patient` and `noddi_long` using your pipeline.
# `noddi_long` has columns: UnQID, PIDN (optional), DCDate, Age, sex,
#                           CDR_class (optional), metric, parcel, value
# Metrics are one of: "ficv","fiso","odi"
# -------------------------------------------------------------------

# --- 0) Minimal hygiene -------------------------------------------------------
noddi_clean <- noddi_long %>%
  filter(metric %in% c("ficv","fiso","odi")) %>%
  filter(!is.na(Age), !is.na(value), !is.na(sex)) %>%
  mutate(sex = factor(sex, levels = c("Male","Female")))  # keep consistent

# --- 1) Robust betas: value ~ Age + sex --------------------------------------
compute_rlm_betas <- function(data = noddi_clean, metrics = c("ficv","fiso","odi"), min_n = 8){

  data %>%
    filter(metric %in% metrics) %>%
    group_by(metric, parcel) %>%
    group_modify(~{
      x <- .x

      # guard against tiny groups
      if (nrow(x) < min_n || length(unique(x$Age)) < 3) {
        return(tibble(term = "Age", estimate = NA_real_, std.error = NA_real_,
                      conf.low = NA_real_, conf.high = NA_real_, n = nrow(x)))
      }

      # fit RLM with Huber psi (robust to outliers)
      fit <- tryCatch(
        rlm(value ~ Age + sex , data = x, psi = psi.huber, maxit = 100),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(tibble(term = "Age", estimate = NA_real_, std.error = NA_real_,
                      conf.low = NA_real_, conf.high = NA_real_, n = nrow(x)))
      }

      # broom::tidy for rlm returns estimate & std.error; p-values are not standard
      tt <- broom::tidy(fit) %>% filter(term == "Age")

      if (nrow(tt) == 0) {
        return(tibble(term = "Age", estimate = NA_real_, std.error = NA_real_,
                      conf.low = NA_real_, conf.high = NA_real_, n = nrow(x)))
      }

      # Wald-style 95% CI (approx) since rlm doesn’t give exact CIs by default
      est  <- tt$estimate[1]
      se   <- tt$std.error[1]
      ci_l <- est - 1.96 * se
      ci_h <- est + 1.96 * se

      tibble(
        term      = "Age",
        estimate  = est,
        std.error = se,
        conf.low  = ci_l,
        conf.high = ci_h,
        n         = nrow(x)
      )
    }) %>%
    ungroup()
}

rlm_betas <- compute_rlm_betas()

# --- 2) Order parcels by ascending Age beta on fICV --------------------------
ficv_order <- rlm_betas %>%
  filter(metric == "ficv", term == "Age") %>%
  arrange(-estimate) %>%
  pull(parcel)

# if some parcels only exist in fiso/odi, keep them after ficv list
all_parcels <- rlm_betas %>% distinct(parcel) %>% pull(parcel)
parcel_levels <- c(ficv_order, setdiff(all_parcels, ficv_order))

rlm_betas <- rlm_betas %>%
  mutate(parcel = factor(parcel, levels = parcel_levels))

# --- 3) Plot: β(Age) ± 95% CI per parcel, faceted by metric ------------------
# Tip: Use coord_flip() for long parcel lists.
library(patchwork)

p1 <- rlm_betas %>% filter(metric == "ficv", term=="Age") %>%
  ggplot(aes(x=parcel, y=estimate)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "red", width = 0.3) +
  geom_point(size = 1.8, color = "black") +
  labs(y="β_Age", x=NULL, title="ficv") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_scientific()) +
  theme_bw(11) +
  theme(axis.text.y=element_text(size=7) ,
  axis.text.x = element_text(angle=45, hjust=1),
  	plot.margin = margin(t = 10, r = 0, b = 50, l = 50)) # top, right, bottom, left))

p2 <- rlm_betas %>% filter(metric == "fiso", term=="Age") %>%
  ggplot(aes(x=parcel, y=estimate)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "red", width = 0.3) +
  geom_point(size = 1.8, color = "black") +
  labs(y=NULL, x=NULL, title="fiso") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_scientific()) +
  theme_bw(11) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1) ,
       plot.margin = margin(t = 10, r = 0, b = 50, l = 0)) # top, right, bottom, left)

p3 <- rlm_betas %>% filter(metric == "odi", term=="Age") %>%
  ggplot(aes(x=parcel, y=estimate)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), color = "red", width = 0.3) +
  geom_point(size = 1.8, color = "black") +
  labs(y=NULL, x=NULL, title="odi") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_scientific()) +
  theme_bw(11) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1) ,
        plot.margin = margin(t = 10, r = 0, b = 50, l = 0)) # top, right, bottom, left)

beta_plot <- p1 + p2 + p3 + plot_layout(ncol=3)
beta_plot

# --- 4) (Optional) Save results table ----------------------------------------
# write_csv(select(rlm_betas, metric, parcel, estimate, std.error, conf.low, conf.high, n),
#           "noddi_rlm_age_betas.csv")


















