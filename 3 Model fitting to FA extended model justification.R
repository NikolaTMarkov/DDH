# NTM — August 2025
# Goal: choose between none/linear/quadratic age effects (gender-adjusted),
# and attach transparent decision-support statistics for each tract.


 library ( tidyverse )


# --- Load data -----------------------------------------------------------
 load ( '~/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/FA_UCSF_db_082025.RData' )

# Confirm that gender is a factor and define age groups
 firstVisit_fa_long_combined <- firstVisit_fa_long_combined %>%
		mutate ( gender = as.factor ( gender ) ) %>%
  		mutate(
    		age_group = cut ( Age , breaks = seq ( floor ( min ( Age , na.rm = TRUE ) ) , ceiling ( max ( Age , na.rm = TRUE ) ) , by = 5 ) , include.lowest = TRUE ) ,
    		age_bin_midpoint = map_dbl ( age_group , ~ mean ( as.numeric ( stringr::str_extract_all ( .x , "[0-9\\.]+" ) [[ 1 ]] ) ) ) )

# Model selection ------------------------------------------
# Decision rule:
#  1) If curvature is significant (quadratic vs linear) AND has evidence ( ΔAIC <= -2 OR w_nonlin >= 0.70 ): choose "non_linear"
#  2) Else, if overall Age (linear vs null) is significant: choose "linear"
#  3) Else: "none"
select_best_model <- function ( df , alpha = 0.05 , aic_cut = -2 , w_cut = 0.70 ) {

 dat <- df %>%
		select ( FA , Age , gender ) %>%
		filter ( is.finite ( FA ) , is.finite ( Age ) , !is.na ( gender ) )

  		# Safety: if too few rows or no variation, return NAs
 if ( nrow ( dat ) < 10 || dplyr::n_distinct ( dat$Age ) < 5 ) {
		return ( tibble (
			best_model	= "none",
			n_obs		= nrow(dat),
			age_overall_p	= NA_real_ , f_age = NA_real_ , df_age_num = NA_real_ , df_age_den = NA_real_ ,
			curvature_p   = NA_real_ , f_curv = NA_real_ , df_curv_num = NA_real_ , df_curv_den = NA_real_ ,
			aic_lin	= NA_real_ , aic_quad = NA_real_ , delta_aic = NA_real_ ,
			w_linear = NA_real_ , w_nonlin = NA_real_ ,
			bic_lin = NA_real_ , bic_quad = NA_real_ , delta_bic = NA_real_ ,
			adjr2_lin = NA_real_ , adjr2_quad = NA_real_ , delta_adjr2 = NA_real_ ,
			linear_slope = NA_real_ , linear_slope_se = NA_real_ , linear_p = NA_real_ ,
			quad_p	= NA_real_ )) }# end table for results collection

  # Fit: Null, Linear, Quadratic (all gender-adjusted)
  null_fit   <- lm (FA ~ gender ,  data = dat )
  linear_fit <- lm ( FA ~ Age + gender ,  data = dat )
  quad_fit   <- lm ( FA ~ poly ( Age , 2 ) + gender , data = dat )

  # Nested F-tests
  a_age  <- anova ( null_fit ,  linear_fit )
  a_curv <- anova(linear_fit, quad_fit)

  age_overall_p <- a_age$`Pr(>F)` [ 2 ]
  f_age         <- a_age$F [ 2 ]
  df_age_num    <- a_age$Df [ 2 ]
  df_age_den    <- a_age$Res.Df [ 2 ]

  curvature_p   <- a_curv$`Pr(>F)` [ 2 ]
  f_curv        <- a_curv$F [ 2 ]
  df_curv_num   <- a_curv$Df [ 2 ]
  df_curv_den   <- a_curv$Res.Df [ 2 ]

  # Information criteria + weights (linear vs quadratic set)
  aic_lin  <- AIC ( linear_fit ) ;  aic_quad  <- AIC ( quad_fit )
  bic_lin  <- BIC ( linear_fit ) ;  bic_quad  <- BIC ( quad_fit )
  delta_aic <- aic_quad - aic_lin
  delta_bic <- bic_quad - bic_lin

  aics <- c ( aic_lin, aic_quad )
  dAIC <- aics - min ( aics )
  ww   <- exp ( -0.5 * dAIC ) / sum ( exp ( -0.5 * dAIC ) )
  names ( ww ) <- c ( "linear" , "quadratic" )
  w_linear <- unname ( ww [ "linear" ] )
  w_nonlin <- unname ( ww [ "quadratic" ] )

  # Adjusted R²
  adjr2_lin  <- summary ( linear_fit )$adj.r.squared
  adjr2_quad <- summary ( quad_fit )$adj.r.squared
  delta_adjr2 <- adjr2_quad - adjr2_lin

  # Linear slope
  lin_tab <- summary ( linear_fit )$coefficients
  linear_slope    <- unname ( lin_tab [ "Age" , "Estimate" ] )
  linear_slope_se <- unname ( lin_tab [ "Age" , "Std. Error" ] )
  linear_p        <- unname ( lin_tab [ "Age" , "Pr(>|t|)" ] )

  # Quadratic term p-value
  coef_tab <- summary ( quad_fit )$coefficients
  quad_row <- grep ( "poly\\(Age, 2\\).*2$" , rownames ( coef_tab ) , value = TRUE )
  quad_p   <- if ( length ( quad_row )) coef_tab[quad_row , "Pr(>|t|)" ] else NA_real_

  # Decision rule
  choose_quad <- 	( !is.na ( age_overall_p ) && age_overall_p < alpha ) && #Has linear effect of age at least
  					( !is.na ( curvature_p ) && curvature_p < alpha ) && 	#Has significant polynomial
                 	( ( !is.na ( delta_aic ) && delta_aic <= aic_cut ) ||  #Evidence of non trivial curvature
                   	( !is.na ( w_nonlin )   && w_nonlin >= w_cut ) )		#Evidence of non trivial curvature

  best_model <- case_when (
  				choose_quad ~ "non_linear" ,
				!is.na ( age_overall_p ) && age_overall_p < alpha ~ "linear" ,
				TRUE ~ "none" )
	tibble (
		best_model = best_model ,
		n_obs = nrow ( dat ) ,
    # Overall Age effect (Linear vs Null)
    	age_overall_p , f_age , df_age_num , df_age_den ,
    # Curvature beyond linear (Quadratic vs Linear)
   		curvature_p ,  f_curv , df_curv_num , df_curv_den ,
    # ICs and evidence
    	aic_lin , aic_quad , delta_aic , w_linear , w_nonlin ,
    	bic_lin , bic_quad , delta_bic,
    # Fit quality
    	adjr2_lin , adjr2_quad , delta_adjr2 ,
    # Coeff diagnostics
		linear_slope , linear_slope_se , linear_p , quad_p ) }

# --- Run per tract ------------------------------------------------------------
 model_results <- firstVisit_fa_long_combined %>%
			group_by ( tract ) %>%
			nest ( ) %>%
			mutate ( model_info = purrr::map ( data , select_best_model ) ) %>%
			unnest ( model_info ) %>%
  			select ( -data ) %>%
			ungroup	(	)

# Multiple-testing control across tracts (optional but recommended)
 model_results <- model_results %>%
  mutate ( 	age_overall_p_fdr = p.adjust ( age_overall_p , method = "fdr" ) ,
    		curvature_p_fdr   = p.adjust ( curvature_p ,   method = "fdr" ) )

# Order table: non_linear first ( by curvature_p ) , then linear ( by age_overall_p ) , then none
 model_results_ordered <- model_results %>%
		mutate ( best_model = factor ( best_model , levels = c ( "non_linear" , "linear" , "none" ) ) ) %>%
		arrange ( best_model ,
    				case_when ( best_model == "non_linear" ~ curvature_p ,
              					best_model == "linear"     ~ age_overall_p ,
								TRUE ~ Inf ) )

# --- (Optional) predictions and plot ------------------------------------------
# For visualization we predict using the most frequent gender within each tract.
	predict_fa_bins <- function ( df , model_type ) {
			dat <- df %>% arrange ( age_bin_midpoint )
			ref_gender <- dat %>%
			count ( gender , sort = TRUE ) %>% slice( 1 ) %>% pull( gender )

			if ( model_type == "linear" ) {
				fit <- lm ( FA ~ Age + gender , data = dat )
			} else if ( model_type == "non_linear" ) {
				fit <- lm ( FA ~ poly ( Age , 2 ) + gender, data = dat )
			} else {
				return ( tibble (
				age_bin_midpoint = unique ( dat$age_bin_midpoint ) ,
      			predicted_FA = mean ( dat$FA , na.rm = TRUE ) 
      		)) }

  tibble ( age_bin_midpoint = sort ( unique ( dat$age_bin_midpoint ) ) ) %>%
		mutate ( predicted_FA = predict ( fit , newdata = tibble ( Age = age_bin_midpoint , gender = ref_gender))) }

# Build predictions per tract using the chosen model
 fa_predictions <- firstVisit_fa_long_combined %>%
		left_join ( model_results_ordered %>% select ( tract , best_model ) , by = "tract" ) %>%
		group_by ( tract , best_model ) %>%
		nest ( ) %>%
		mutate ( predictions = purrr::map2 ( data , best_model , predict_fa_bins ) ) %>%
		select ( -data ) %>%
		unnest ( predictions ) %>%
		ungroup ( )

# Preserve tract facet order
 tract_order <- model_results_ordered$tract
 firstVisit_fa_long_combined <- firstVisit_fa_long_combined %>%
			mutate ( tract = factor ( tract , levels = tract_order ) )
 fa_predictions <- fa_predictions %>%
			mutate ( tract = factor ( tract , levels = tract_order ) )

# Plot 
	ggplot ( ) +
	geom_point ( data = firstVisit_fa_long_combined ,
             aes ( x = Age , y = FA ) ,
             color = "grey70" , alpha = 0.4 , size = 0.5 ) +
  	stat_summary ( data = firstVisit_fa_long_combined,
       		aes ( x = age_bin_midpoint , y = FA ) ,
       		fun = mean , geom = "point" , color = "black" , size = 1.5 ) +
	geom_line ( data = fa_predictions ,
			aes ( x = age_bin_midpoint , y = predicted_FA , color = best_model ) , linewidth = 1 ) +
  	facet_wrap	(	~ tract , scales = "free_y" ) +
  	scale_color_manual ( values = c ( 	"linear" = "#d95f02" ,
                                		"non_linear" = "#1b9e77" ,
                                		"none" = "#7570b3" ) ) +
	labs ( title = "Tract FA vs Age (5-year bins)" ,
	subtitle = "Gender-adjusted Linear vs Quadratic model choice with transparent statistics",
       	x = "Age (years)" , y = "Fractional Anisotropy (FA)" , color = "Chosen model" ) +
  	theme_bw ( base_size = 13 ) +
	theme (	strip.text = element_text ( size = 7 ) ,
        	legend.position = "top" ,
        	axis.text = element_text ( size = 8 ) ,
        	plot.title = element_text ( face = "bold" ) )

# Inspect results
 print ( model_results_ordered , n = nrow ( model_results_ordered ) )
 # path <- '/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Figures/Figure sources/FAMD evolution with age sex corrected'
 # write.csv (model_results_ordered, file = paste ( path , "Statistics of FA non linear vs linear model choice.csv" , sep = "/" ) )
 
 
 
 dev.new ( width = 7 , height = 4)
# ================================
# Extra plot: only 3 selected tracts
# ================================
selected_tracts <- c ( "uncinate_fasciculus" ,
                     "anterior_corona_radiata" ,
                     "posterior_limb_of_internal_capsule" )

ggplot() +
  geom_point ( data = filter ( firstVisit_fa_long_combined , tract %in% selected_tracts ) ,
             aes ( x = Age , y = FA ) ,
             color = "grey70" , alpha = 0.4 , size = 0.75 ) +
  stat_summary ( data = filter ( firstVisit_fa_long_combined , tract %in% selected_tracts ) ,
                aes ( x = age_bin_midpoint , y = FA ) ,
               fun = mean , geom = "point", color = "black" , size = 1.5 ) +
  geom_line ( data = filter ( fa_predictions , tract %in% selected_tracts ) ,
            aes ( x = age_bin_midpoint , y = predicted_FA , color = best_model ) ,
            linewidth = 1 ) +
  facet_wrap ( ~ tract , scales = "free_y" ) +
  scale_color_manual ( values = c ( "linear" = "#d95f02" ,
                                "non_linear" = "#1b9e77" ,
                                "none" = "#7570b3" ) ) +
  labs ( title = "Selected Pathways: FA vs Age (Sex-adjusted)" ,
       x = "Age" ,
       y = "Fractional Anisotropy (FA)" ,
       color = "Model Match" ) +
  theme_bw ( base_size = 13 ) +
  theme ( strip.text = element_text ( size = 9 ) ,
        legend.position  = "top" ,
        axis.text = element_text (size = 12 ) ,
        plot.title = element_text ( face = "bold" ) )

