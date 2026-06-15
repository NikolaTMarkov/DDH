#:NTM 2025 apr analysis DTI Data
 library ( MASS ) #for the RLM
 library ( rstatix ) #For wilcoxon test of FA ~ CDR
 library ( boot ) #for bootstrapping the CI estimates
 library ( forcats ) #Just to reorder the Y axis of the beta estimates
 library ( tidyverse )
 library ( ggpmisc ) #to add statistics on plot
 library ( broom ) #for the stats table
 library (zoo) #for rolling and moving calculations

#Load
load ( "/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/UCSF_db_Branch_Summer2024.RData" )

 dti_fa <- as_tibble ( UCSF_db$imaging_dti_combat_fa )
 cdr <- as_tibble ( UCSF_db$cdr )
 demographics <- as_tibble ( UCSF_db$demographics )
 
 #Sex 
 demographics <- demographics %>%
 	mutate ( gender = factor ( gender , levels = c ( 1 , 2  ) , labels = c ( "male" , "female"  ) ) ) 
 
## #Create a cognition factor and colors for plot of cognition
 cdr <- cdr %>%
	#it is important to match the level order for the two
			mutate(
    			CDR_class = factor (
      							case_when (
        						is.na ( cdr_global ) ~ "Unknown" ,    # Gray for NA values
        						cdr_global == 0 ~ "CU" ,          # chartreuse4 for unimpaired
        						cdr_global == 0.5 ~ "MCI" ,       # Orange for MCI
        						TRUE ~ "SCI"                     # Red for severe cognitive impairment
      								) , levels = c ( "CU", "Unknown", "MCI", "SCI" ) ) ,
    			CDR_palette = factor (
      							case_when (
        						is.na ( cdr_global ) ~ "gray70" ,
        						cdr_global == 0 ~ "chartreuse4" ,
        						cdr_global == 0.5 ~ "orange" ,
        						cdr_global > 0.5  ~ "red" 
        							) , levels = c  (  "chartreuse4"  ,  "gray70"  ,  "orange" ,  "red"  )  ) ) #end cdr class and palette
 color_mapping  <-  c (  "Unknown" = "gray70" ,
 				 "CU" = "chartreuse4" , 
 				 	"MCI" = "orange" , 
 				 		"SCI" = "red"  )

#Associate with age and select only the first recorded visit

  firstVisit_fa <- UCSF_db$general %>% 
  					select ( UnQID, PIDN, Age = age_at_DCDate ) %>%
  					left_join ( dplyr:: select (  demographics , UnQID , gender ) , by  = "UnQID"  ) %>%
  					left_join ( UCSF_db$imaging_dti_combat_fa , by = "UnQID" ) %>%
  					left_join ( dplyr::select ( cdr , UnQID, CDR_global = cdr_global, CDR_class , CDR_palette   ) , by  = "UnQID" ) %>%
  					filter ( complete.cases ( . ) ) %>%
  						group_by ( PIDN ) %>%
  						arrange ( DCDate ) %>%
  						slice ( 1 ) %>% #Select the first visit of each patient based on the DCDate
  						ungroup ( )
  						
 firstVisit_fa_long <- firstVisit_fa %>%
					select ( -date_diff , -DCDate , -PIDN , - CDR_global) %>%  # remove unwanted columns
 					pivot_longer (
						cols = - c ( UnQID ,  Age , gender , CDR_class , CDR_palette ) ,
						names_to = "tract" ,
						values_to = "FA" )
						
#Combine left and right as we find them to be simiar in pattern to each other
# First, identify left/right tracts and create a unified name
firstVisit_fa_long_combined <- firstVisit_fa_long %>%
					 mutate(
						tract_base = str_remove(tract, "_[lr]$"),  # removes the _l or _r suffix
						hemisphere = case_when(
							str_detect(tract, "_l$") ~ "left",
							str_detect(tract, "_r$") ~ "right",
							TRUE ~ "midline" )  ) %>%
					group_by (UnQID , Age , gender , CDR_class , CDR_palette , tract_base ) %>%
					summarise ( FA = mean (FA , na.rm = TRUE ) , .groups = 'drop' ) %>% 
					rename ( tract = tract_base )

##########################
####################
#Outlier detection and plot
#(need to remove them before the long combined first visit)
 
 
# Compute residuals and robust outlier detection per tract
outlier_detection_robust <- firstVisit_fa_long_combined %>%
  group_by(tract , gender) %>%
  arrange(Age , .by_group = TRUE ) %>%
  mutate(
    # Loess fit per tract
    loess_fit = predict(loess(FA ~ Age, span = 0.75 , data = pick( everything () ) ) ),
    
    # Residuals
    residuals = FA - loess_fit,
    
    # Robust spread estimation per tract
    med_resid = median(residuals, na.rm = TRUE),
    iqr_resid = IQR(residuals, na.rm = TRUE),
    
    # Calculate fences using standard 1.5 IQR rule
    upper_fence = med_resid + 3 * iqr_resid,
    lower_fence = med_resid - 3 * iqr_resid,
    
    #mad_resid = mad(residuals, constant = 1.4826, na.rm = TRUE),
	#upper_fence = med_resid + 3 * mad_resid ,
	#lower_fence = med_resid - 3 * mad_resid ,
	
    # Identify outliers clearly
    is_outlier = residuals > upper_fence | residuals < lower_fence
  ) %>%
  ungroup()

# Visualization to verify results clearly
# Split data by gender
outlier_male <- filter(outlier_detection_robust, gender == "male")
outlier_female <- filter(outlier_detection_robust, gender == "female")

# Plot for males
p_male <- ggplot(outlier_male, aes(x = Age, y = FA)) +
  geom_point(aes(color = is_outlier), size = 2, alpha = 0.7) +
  geom_line(aes(y = loess_fit), color = "blue", linewidth = 1) +
  facet_wrap(~tract, scales = "free_y") +
  scale_color_manual(values = c("chartreuse4", "red")) +
  labs(
    title = "Robust Outlier Detection – Males",
    color = "Outlier"
  ) +
  theme_bw()

# Plot for females
p_female <- ggplot(outlier_female, aes(x = Age, y = FA)) +
  geom_point(aes(color = is_outlier), size = 2, alpha = 0.7) +
  geom_line(aes(y = loess_fit), color = "blue", linewidth = 1) +
  facet_wrap(~tract, scales = "free_y") +
  scale_color_manual(values = c("chartreuse4", "red")) +
  labs(
    title = "Robust Outlier Detection – Females",
    color = "Outlier"
  ) +
  theme_bw()

# Print plots separately
print(p_male)
dev.new()
print(p_female)

 
 
#####Remove the outliers from the Dataset
#outilier in one tract removes the whole scan 
 OutliersFA <- outlier_detection_robust %>% filter ( is_outlier == TRUE ) %>% 
 											distinct( UnQID) %>%  
 											pull( UnQID )
 firstVisit_fa_long_combined <- firstVisit_fa_long_combined %>% 
 								filter ( !UnQID %in% OutliersFA ) 
 
 # SOM Figure 1 plot of the fa against age
 
# Plot: LOESS by gender and CDR class
 firstVisit_fa_long_combined %>%
	group_split ( gender ) %>%
	
	walk( ~ {
		df_gender  <-  .x  %>%  arrange ( CDR_class )
		g <- unique ( df_gender$gender )
		p <- ggplot ( df_gender , aes ( x = Age , y = FA , color = CDR_class ) ) +
			geom_point ( alpha = 0.6 , size = 1.5 ) +
			geom_smooth ( method = "loess" , se = TRUE , color = "blue" ) +
			scale_color_manual ( values = color_mapping ) +
			facet_wrap ( ~tract , scales = "free_y" ) +
				labs (
				title = paste ( "Fractional Anisotropy by Age –" , g ) ,
				x = "Age (years)",
				y = "Average FA per tract",
				color = "Cognition" ) +
			theme_bw ( ) +
      		theme (
        		strip.text = element_text ( size = 8 ) ,
        		axis.text.x = element_text ( size = 7 ) ,
        		axis.text.y = element_text ( size = 7 ) )
			print ( p )
			dev.new ( ) } )

					
#Resistent linear model to identify if AGE affects FA
#RLM does not assume normal distribution therefore it has not P-value
# instead we estimate significance of the effect by calculating bootstraped confidence intervals (CI)
# If CI does include zero the effect is likelly not significant if it does not inculde zero it is



# Bootstrap function to extract Age coefficient
 boot_age_coef <- function ( data , indices ) {
	fit <- rlm ( FA ~ Age + gender , data = data [ indices , ] )
	coef ( fit ) [ "Age" ] }

# Run rlm and bootstrap CI per tract
 boot_results <- firstVisit_fa_long_combined %>%
	group_split ( tract ) %>%
		map_df ( function ( df ) {
    tract_name <- unique ( df$tract )

    	# Fit RLM
    	model <- rlm ( FA ~ Age + gender , data = df )
    	est <- coef ( model ) [ "Age" ]

    	# Bootstrap CI
    	boot_out <- boot ( data = df , statistic = boot_age_coef , R = 1000)
    	ci <- boot.ci ( boot_out , type = "perc" )

    tibble (
      tract = tract_name ,
      estimate = est ,
      ci_lower = ci$perc [ 4 ] ,
      ci_upper = ci$perc [ 5 ] )
 } ) #End boot function




 ggplot ( boot_results , aes ( x = fct_reorder ( tract , -estimate ) , y = estimate ) ) +
	geom_point ( size = 2 ) +
	geom_errorbar ( aes ( ymin = ci_lower , ymax = ci_upper ) , width = 0.2 , color = "red" ) +
	geom_hline ( yintercept = 0 , linetype = "dashed" , color = "gray50" ) +
	coord_flip ( ) +
	labs (
    title = "Gender adjusted RLM effect of Age on FA" ,
    x = "Tract (sorted by effect size)" ,
    y = "Estimated Beta (Age Effect)" ) +
	theme_minimal ( base_size = 12 )




###############Boxplots
# Add age group as ordered factor
 firstVisit_fa_long_combined <- firstVisit_fa_long_combined %>%
	mutate(
		age_group = cut (
		Age ,
		breaks = c ( -Inf  , 65 ,  Inf ) , 
		labels = c ( "Pre 65" , "Post 65" ) ,
		right = FALSE ) )

# Pairwise Wilcoxon tests within each age group




within_age_stats <- firstVisit_fa_long_combined %>%
  group_by(age_group, gender) %>%
  filter(n_distinct(CDR_class) >= 2) %>%
  rstatix::wilcox_test(FA ~ CDR_class, p.adjust.method = "fdr") %>%
  mutate(comparison = paste(group1, "vs", group2)) %>%
  select(age_group, gender, comparison, p.adj) %>%
  arrange(age_group, gender, p.adj)

print(within_age_stats, n = nrow(within_age_stats))




#PLOT firstVisit_fa_long_combined %>%
 firstVisit_fa_long_combined %>%
	group_split ( gender ) %>%
		walk ( ~ {
		df_gender <- .x
		g <- unique ( df_gender$gender )

		p <- ggplot ( df_gender , aes ( x = age_group , y = FA , fill = CDR_class ) ) +
      	geom_jitter (
        aes ( color = CDR_class ) ,
        alpha = 0.4 , size = 0.5,
        position = position_jitterdodge ( jitter.width = 0.2 , dodge.width = 0.8 ) ) +
        geom_boxplot ( outlier.shape = NA , alpha = 0.8 , position = position_dodge ( width = 0.8 ) ) +
		facet_wrap ( ~tract , scales = "free_y" ) +
		scale_fill_manual ( values = color_mapping ) +
      	scale_color_manual ( values = color_mapping ) +
      	labs (
        title = paste ( "FA by Age Group, Tract, and CDR —" , g ) ,
        x = "Age Group" ,
        y = "Fractional Anisotropy (FA)" ,
        fill = "Cognition" ,
        color = "Cognition" ) +
      	theme_bw ( ) +
      	theme (
        	strip.text = element_text ( size = 9 ) ,
        	axis.text.x = element_text ( angle = 30 , hjust = 1 ) )

	print ( p )
	dev.new ( ) })





#save ( list = c ( "color_mapping" , "firstVisit_fa_long_combined" , "mta" , "UCSF_db" ) , file = "FA_UCSF_db_082025.RData" )                      


