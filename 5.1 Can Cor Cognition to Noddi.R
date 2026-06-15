######################
#NTM Oct 2025
# Calculates canonical correlation between white matter NODDI maps and composite cognitive variables
# Remove high missingness rows, impute when only few variables not estimated
# Add in test that could be up to year in time separate from MRI scan
#Compute canonical correlation plot and explore if AGE, Inflammation, Sex are  critical variables
##########################################################################################################
######################################################################################################################
###########################
# NTM Oct 2025
# Match NODDI scans to nearest cognitive assessments (±1 year)
###########################

 library ( tidyverse )
 library ( yacca )
 library ( ggpubr )
 library ( patchwork )
 library ( viridis )
 library ( reshape2 )



tolerance_days <- 365 * 2


# --- Load data --------------------------------------------------
load ( "~/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/UCSF_db_Branch_Summer2024.RData" )
load ( "~/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/LongWideNODDI_DWI.RData" )

patient <- patient %>%
			left_join ( UCSF_db$cdr %>% select ( "UnQID" , "cdr_global") , by = "UnQID"   )

#List cognitive variables 
 bedside_vars <- c ( "gds_tot" , "mmse_tot" , "memoryzscore" , "bsexzscore" , "an_corr" )
 infoprocessingspeed_vars <- c ( "spatial", "verbal" ) #, "verbal"
 cog_vars <- c ( infoprocessingspeed_vars , bedside_vars )


 bedside <- UCSF_db$bedside [ , c ( "UnQID" , "DCDate" , bedside_vars ) ]
 infoprocessingspeed <- UCSF_db$infoprocessingspeed [ , c ( "UnQID" , "DCDate" , infoprocessingspeed_vars ) ]


 Proc_speed <- infoprocessingspeed %>%
					left_join ( patient %>% select ( UnQID , PIDN ) , by = "UnQID" ) %>%
					relocate ( UnQID , PIDN , DCDate , .before = everything ( ) )

 Proc_exec <- bedside %>%
					left_join ( patient %>% select ( UnQID , PIDN ) , by = "UnQID" ) %>%
					relocate ( UnQID , PIDN , DCDate , .before = everything ( ) )

 noddi_wide <- noddi_wide %>%
					left_join ( patient %>% select ( UnQID , Age , sex ) , by = "UnQID" ) %>%
					relocate ( UnQID , Age , sex , DCDate , .before = everything ( ) )

 NoddiCog <- noddi_wide
			NoddiCog [ cog_vars ] <- NA_real_


###################################################
# --- Matching function for a single row ---
####################################################################
match_one_row <- function ( a ) {
  pid  <- a [[ "PIDN" ]]
  date <- as.Date ( a [[ "DCDate" ]] )

  out <- a  # start with original row

  # --- Infoprocessing speed vars ---
  if ( !is.na ( pid ) && !is.na ( date ) ) {
    ps_rows <- Proc_speed [ Proc_speed$PIDN == pid ,  ,  drop = FALSE ]
    if ( nrow ( ps_rows ) > 0 ) {
      for ( v in infoprocessingspeed_vars ) {
        ok <- which ( abs ( as.numeric ( difftime ( ps_rows$DCDate , date , units = "days" ) ) ) <= tolerance_days )
        if ( length ( ok ) > 0 ) {
          idx <- ok [ which.min ( abs ( as.numeric ( difftime ( ps_rows$DCDate [ ok ] , date , units = "days" ) ) ) ) ]
          out [[ v ]] <- ps_rows [[ v ]] [ idx ]
        }
      }
    }
  }

  # --- Bedside vars ---
  if ( !is.na ( pid ) && !is.na ( date ) ) {
    be_rows <- Proc_exec [ Proc_exec$PIDN == pid , , drop = FALSE ]
    if ( nrow ( be_rows ) > 0 ) {
      for ( v in bedside_vars ) {
        ok <- which ( abs ( as.numeric ( difftime ( be_rows$DCDate , date , units = "days" ) ) ) <= tolerance_days )
        if ( length ( ok ) > 0 ) {
          idx <- ok [ which.min ( abs ( as.numeric ( difftime ( be_rows$DCDate [ ok ] , date , units = "days" ) ) ) ) ]
          out [[ v ]] <- be_rows [[ v ]] [ idx ]
        }
      }
    }
  }

  return ( out )
}

# --- Apply row-wise data fill in over NoddiCog ---
NoddiCog_filled <- NoddiCog %>%
  						mutate ( row_id = row_number ( ) ) %>%
						split ( .$row_id ) %>%
  						map_dfr (match_one_row ) %>% 
  						select ( -row_id ) %>%
						mutate ( across ( all_of ( cog_vars ) , as.numeric ) )
				
NoddiCog_filtered <- NoddiCog_filled %>%
							filter ( rowMeans ( is.na ( select ( . , all_of ( cog_vars ) ) ) ) <= 0.35 ) #keep only rows with less than x% missing



  
 impute_age_window <- function ( df , vars , age_col = "Age" , window = 5 ) {
							df_out <- df
  
  							# keep original data reference (so we only use real observed values)
  							df_ref <- df
  
  								for ( v in vars ) {
    								for ( i in seq_len ( nrow ( df_out ) ) ) {
      										if ( is.na ( df_out [[ v ]][ i ] )) {
												current_age <- df_out [[ age_col ]][ i ]
        										if ( is.na ( current_age ) ) next
        
        								# subset reference data to same-age window and non-missing var
        				neighbors <- df_ref %>%
          								filter ( !is.na ( .data [[ v ]] ),
                								 abs ( .data [[ age_col ]] - current_age ) <= window )
        
        				if (nrow ( neighbors ) > 0 ) {
          					df_out [[ v ]][ i ] <- median ( neighbors [[ v ]] , na.rm = TRUE )
 						       }
    						  }
    						}
  						}
  				df_out
				}#End impute age window function

age_window <- 5

NoddiCog_imputed <- impute_age_window (
						df   	=	 NoddiCog_filtered ,
  						vars 	=  	cog_vars ,
  						age_col	= 	"Age",
  						window	=	age_window )

 NoddiCog_imputed <-  NoddiCog_imputed [ complete.cases ( NoddiCog_imputed  ) , ] 
  
 #───────────────────────────────────────────────────────────────
# Keep only the first visit per patient (earliest DCDate)
#───────────────────────────────────────────────────────────────
NoddiCog_imputed <- NoddiCog_imputed %>%
						group_by(PIDN) %>%
						arrange(DCDate, .by_group = TRUE) %>%
						slice_tail(n = 1) %>% ##For last visit use slice_tail
						ungroup()
 
  
  
##########################################################################################################
##########################################################################################################
#───────────────────────────────────────────────────────────────
# 1. Define variable groups
#───────────────────────────────────────────────────────────────
 cog_set   <- cog_vars                                     # already defined earlier
 img_set   <- grep ( "\\.(ficv|fiso|odi)$" , names ( NoddiCog_imputed )  , value = TRUE ) # all noddi metrics

#───────────────────────────────────────────────────────────────
# 2. Prepare matrices
#───────────────────────────────────────────────────────────────
 X_img <- NoddiCog_imputed %>%
					select ( all_of ( img_set ) ) %>%
					#mutate ( across ( everything ( ) , scale ) ) %>% 
					as.matrix ()

 Y_cog <- NoddiCog_imputed %>%
					select ( all_of ( cog_set ) ) %>%
					#mutate ( across ( everything ( ) , scale ) ) %>% 
					as.matrix ( )

#───────────────────────────────────────────────────────────────
# 3. Run canonical correlation analysis
#───────────────────────────────────────────────────────────────
cca_res <- yacca::cca(X_img, Y_cog, 
                     xscale = TRUE,    # Scale X to unit variance
                     yscale = TRUE,    # Scale Y to unit variance  
                     standardize.scores = TRUE)
#yacca provides the canonical variates directly when standardize.scores = TRUE

U <- cca_res$canvarx   # Canonical variates for imaging data
V <- cca_res$canvary   # Canonical variates for cognitive data


#───────────────────────────────────────────────────────────────
# 2. Compute structure correlations (variable loadings)
#───────────────────────────────────────────────────────────────
xstruct <- cor(X_img, U)  # Use original data (scaling handled by yacca)
ystruct <- cor(Y_cog, V)

struct <- list(xstruct = xstruct, ystruct = ystruct)


#───────────────────────────────────────────────────────────────
# 3. Assemble tidy data for scatterplots
#───────────────────────────────────────────────────────────────


###########
#scree plot
#Represents the shared variance between canonical variates. It cumulates above 100 because 
#1 and 2 already explain 103% of what would have been perfect correlation addin 3 is not necessary

#Calculate significance of canonical variates
 cca_ftest <- yacca::F.test.cca ( cca_res )
 print ( cca_ftest )

 
 corrsq_df <- tibble(
		CV = factor ( names ( cca_res$corrsq ) , levels = names ( cca_res$corrsq ) ) ,
  		R2 = as.numeric ( cca_res$corrsq ) ,
		CV_num = 1 : length ( cca_res$corrsq ) ,
		p_value = cca_ftest$p.value  )  %>%
		  			mutate (
						significance = case_when (
									p_value < 0.001 ~ "***" ,
      								p_value < 0.01 ~ "**" , 
									p_value < 0.05 ~ "*" ,
									TRUE ~ "" ) ,
    						label_with_sig = paste0 ( round ( R2 , 2) , significance ) )


 corrsq_df %>%
 		ggplot( aes ( x = CV_num , y = R2 ) ) +
		geom_line ( color = "darkgreen" , linewidth = 1 ) +
		geom_point ( size = 3 , color = "firebrick" ) +
			scale_x_continuous (
					breaks = corrsq_df$CV_num ,
    					labels = corrsq_df$CV ,
    					limits = c ( 0.5 , max ( corrsq_df$CV_num ) + 0.1 ) ) +
    		scale_y_continuous(
  					limits = c ( min ( corrsq_df$R2 ) , max ( corrsq_df$R2 ) * 1.1 ) )+
  		theme_minimal ( base_size = 14 ) +
  		labs (
    				title = "Canonical Correlations (R²)" ,
    				subtitle = "Shared variance per canonical dimension\n*p < 0.05, **p < 0.01, ***p < 0.001" ,
    				x = "Canonical Variate Pair", 
    				y = expression ( R^2 ) ) +
  		geom_text ( aes ( label = label_with_sig  ), vjust = - 0.8 , size = 5 )

################
#Plot correlation matrices
 img_cor <- cor( X_img  )
 img_cor_melted <- melt ( img_cor )

 ggplot ( img_cor_melted , aes ( x = Var1 , y = Var2 , fill = value ) ) +
 		geom_tile ( ) +
		scale_fill_gradient2 ( low = "blue" , high = "red" , mid = "white" , 
                       midpoint = 0 , limit = c ( -1 , 1 ) , space = "Lab" ,
                       name = "Correlation") +
		theme_minimal ( ) +
		theme ( 	axis.text.x = element_blank () ,
  				axis.text.y = element_text ( size = 7 ) ) +
  		labs ( title = "Imaging Variables Correlation Matrix" , x =  "" , y = "" ) +
		coord_fixed ( )

# Compute correlation matrix for cognitive data
 cog_cor <- cor ( Y_cog )
 cog_cor_melted <- melt ( cog_cor )

 ggplot ( cog_cor_melted , aes ( x = Var1 , y = Var2 , fill = value ) ) +
  			geom_tile ( ) +
			scale_fill_gradient2 ( low = "blue" , high = "red" , mid = "white" , 
                       midpoint = 0 , limit = c ( -1 , 1 ) , space = "Lab" ,
                       name = "Correlation" ) +
  			theme_minimal ( ) +
  			theme	( axis.text.x = element_text ( angle = 45 , vjust = 1 , hjust = 1 ) ) +
  			labs ( title = "Cognitive Variables Correlation Matrix" ,
       			x = "", y = "") +
  			coord_fixed()

# Compute cross-correlations between imaging and cognitive variables
 cross_cor <- cor ( X_img , Y_cog )
 cross_cor_melted <- melt ( cross_cor )

 ggplot ( cross_cor_melted , aes ( x = Var1 , y = Var2 , fill = value ) ) +
			geom_tile ( ) +
  			scale_fill_gradient2 ( low = "blue" , high = "red" , mid = "white" , 
                       				midpoint = 0 , limit = c ( -0.25 , 0.25 ) , space = "Lab" ,
                      					 name = "Correlation" ) +
  			theme_minimal ( ) +
  			theme ( axis.text.x = element_text ( angle = 45 , vjust = 1 , hjust = 1 , size = 7 )  ,
  					 legend.position = "bottom") +
  			labs ( title = "Imaging vs Cognitive Variables Correlation Matrix" ,
       				x = "Imaging Variables" , y = "Cognitive Variables" ) +
  			coord_fixed ( )


#───────────────────────────────────────────────────────────────
# 4. Canonical scatterplots (Age-colored)
#───────────────────────────────────────────────────────────────
cv_df <- tibble(
	UnQID	= NoddiCog_imputed$UnQID ,
	DCDate	= NoddiCog_imputed$DCDate ,
	PIDN	= NoddiCog_imputed$PIDN ,
	Age     = NoddiCog_imputed$Age ,
	sex	  = NoddiCog_imputed$sex ,
	CV1_img = as.numeric ( U [ , "CV 1" ] ) ,
	CV1_cog = as.numeric ( V [ , "CV 1" ] ) ,
	CV2_img = as.numeric ( U [ , "CV 2" ] ) ,
	CV2_cog = as.numeric ( V [ , "CV 2" ] )  )
cv_df <- cv_df %>%
  	left_join ( patient %>% select ( UnQID , cdr_global ) , by = "UnQID" )

# Create CDR category for coloring
cv_df <- cv_df %>%
		mutate ( cdr_category = case_when (
		cdr_global == 0 ~ "CDR 0" ,
		cdr_global == 0.5 ~ "CDR 0.5" , 
		cdr_global > 0.5 ~ "CDR >0.5" ,
		TRUE ~ "Missing" ) )
	
	
	
age_range <- range ( cv_df$Age , na.rm = TRUE )

 make_cv_scatter <- function ( df , xvar , yvar , cv_label , colorize = "Age" , AddCDR = FALSE )  {
						df <- df [ order ( df [[ colorize ]] , decreasing = TRUE ) ,  ]

		# Start building the plot
	p <- ggplot ( df , aes ( x = .data [[ xvar ]] , y = .data [[ yvar ]] , color = .data [[ colorize ]] ) ) 
  
		# Add CDR points first (if requested)
  		if ( AddCDR == TRUE ) {
    		cdr_data <- df %>% filter ( cdr_category %in% c ( "CDR 0.5" , "CDR >0.5" ) )
    
	p <- p + 
			geom_point (
					data = cdr_data ,
					aes ( x = .data [[ xvar ]] , y = .data [[ yvar ]] ) ,
					color = case_when (
							cdr_data$cdr_category == "CDR 0.5" ~ "black" ,
							cdr_data$cdr_category == "CDR >0.5" ~ "purple" ) ,
					shape = 16,
					size = 5 ) }
  
  		# Add the main points and other elements
	p <- p + 
    		geom_point ( alpha = 1 , size = 2 ) +
			geom_smooth ( method = "lm" , se = FALSE , color = "black" , size = 0.9 ) +
    		scale_color_gradientn (
      				colours = c ( "#2c7bb6" , "#ffff8c" , "#d7191c" ) ,
      				limits = age_range , 
      				name = "Age (years)" ) +
    		ggpubr :: stat_cor (
      				method = "pearson" ,
      				label.x.npc = 0.1 , label.y.npc = 0.9 ,
      				color = "black" , size = 4 , fontface = "bold" ) +

    		theme_minimal ( base_size = 24) +
    		labs (
				title = paste0 ( "Canonical Variate " , cv_label ) ,
				subtitle = "Imaging vs Cognition" ,
				x = paste0 ( "U" , cv_label , " (Imaging variate)" ) ,
				y = paste0 ( "V" , cv_label , " (Cognitive variate)" ) )
  
  	return ( p )
	}# End Canonical variates colored by age


 make_cv_scatter (cv_df, "CV1_img", "CV1_cog", "1", "Age", FALSE ) 
 make_cv_scatter ( cv_df , "CV2_img" , "CV2_cog" , "2" ,  "sex", FALSE ) +
  			scale_color_manual ( values = c ( "Male" = "steelblue" , "Female" = "tomato" ) )


#───────────────────────────────────────────────────────────────
# 5. Top contributors by structure correlation (CV1)
#───────────────────────────────────────────────────────────────
top_bar_from_struct <- function( df , cv , set , n_top = Inf ) {
  tibble(variable = rownames(df), value = df[, cv]) %>%
    arrange(desc(abs(value))) %>%
    slice_head(n = n_top) %>%
    mutate(
      metric_type = case_when(
        grepl("fiso", variable, ignore.case = TRUE) ~ "fiso",
        grepl("odi", variable, ignore.case = TRUE) ~ "odi",
        grepl("ficv", variable, ignore.case = TRUE) ~ "ficv",
        TRUE ~ "other"
      ),
      fill_color = case_when(
        metric_type == "fiso" ~ "#2ca02c",  # Green
        metric_type == "odi" ~ "#d62728",   # Red
        metric_type == "ficv" ~ "#1f77b4",  # Blue
        TRUE ~ "#212121"    # Charcoal
      ),  variable = gsub("\\.", "•", variable) 
    ) %>%
    ggplot(aes(x = reorder(variable, -value), y = value)) +
    geom_col(aes(fill = fill_color ) , alpha = 0.9 ) + #, width = 1
    geom_hline(yintercept = c ( -0.1, 0.1 ), color = "gray", linetype = "dashed", alpha = 1)+
    scale_fill_identity() +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme( axis.text.y = element_text (size = 12 ) ) +
    labs(
      title = paste0(set, " contributors to ", cv),
      x = NULL, y = "Structure correlation"
    )
}

top_bar_from_struct( xstruct, "CV 1", "Imaging", n_top = 35 )


###
top_bar_from_struct <- function(df, cv, set, n_top = Inf) {
  tibble(variable = rownames(df), value = df[, cv]) %>%
    arrange(desc(abs(value))) %>%
    slice_head(n = n_top) %>%
    mutate(
      # Use the actual variable names as the fill categories
      metric_type = variable
    ) %>%
    ggplot(aes(x = reorder(variable, -value), y = value, fill = metric_type)) +
    geom_col(alpha = 0.9) +
    geom_hline(yintercept = c(-0.5, 0.5), color = "gray", linetype = "dashed", alpha = 1) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 12),
  					 legend.position = "none" ) +
    labs(
      title = paste0(set, " contributors to ", cv),
      x = NULL, y = "Structure correlation",
      fill = "Cognitive Test"
    )
}

top_bar_from_struct(ystruct, "CV 1", "Cognition", n_top = Inf) +
  scale_fill_manual(
    values = c(
      "spatial" = "#8B4513",        # Saddle brown
      "bsexzscore" = "#CD853F",     # Peru
      "verbal" = "#D2691E",         # Chocolate
      "mmse_tot" = "#BC8F8F",       # Rosy brown
      "an_corr" = "#A0522D",        # Sienna
      "gds_tot" = "#6D4C41",        # Coffee brown
      "memoryzscore" = "#8D6E63"    # Taupe brown
    )
  )


#───────────────────────────────────────────────────────────────
# 6. Combine and export
#───────────────────────────────────────────────────────────────
full_plot <- (p_cv1 | p_cv2) / (p_img_cv1 | p_cog_cv1) +
  plot_annotation(
    title = "Canonical Correlation Analysis — White-Matter NODDI × Cognition",
    subtitle = sprintf("Canonical correlations: %s",
                       paste(round(cca_res$corr[1:4], 2), collapse = ", ")),
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

print(full_plot)

ggsave("NTM2025_CCA_full_plot_byAge.pdf", full_plot,
       width = 14, height = 10, units = "in", dpi = 300)



