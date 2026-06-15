#NTM 2025 August
# Generates long and wide forms of noddi and dwi images
#Calculates correlations ( pearson , spearmann?)
#Generates the dataset( "LongWideNODDIDWI.R")

####
#NOTE this script removes all the outliers as identified via IQR by the scripts for figure 2 and 3

 library ( tidyverse )


 load ( "~/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/UCSF_db_Branch_Summer2024.RData" )
 
 load("/Users/nikolamarkov/Library/CloudStorage/GoogleDrive-markovokean@gmail.com/My Drive/7 Iosif brain synchrony and pathways aging/Data/UnQID_outliers.RData")
 
 
 patient <- UCSF_db$general %>%
					transmute ( UnQID ,
            		PIDN ,
            		DCDate_record = as.Date ( DCDate ) ,
            		Age = age_at_DCDate ) %>%
  					left_join ( UCSF_db$demographics %>%
      							transmute (
        							UnQID,
        							sex = factor ( gender ,
                     					levels = c ( 1 , 2 ) ,
                     					labels = c ( "Male", "Female" ) )
      										) , by = "UnQID" )

#attach PIDN to any imaging table by UnQID
 attach_pidn <- function ( df ) { df %>%
								left_join ( select ( patient ,"UnQID" , "PIDN" ) , by = "UnQID" ) %>%
    							relocate ( PIDN , UnQID , .before = 1 ) }

#NODDI: keep only WM parcels (wm_lh_* / wm_rh_*), strip prefixes, average L/R
 collapse_noddi_wm_lr <- function ( df , metric_name ) {
  				# Attach PIDN
  				d <- attach_pidn ( df )

  				# Identify wm_lh_* and wm_rh_* columns
				wm_lh_cols <- grep ( "^wm_lh_" , names ( d ) , value = TRUE )
  				wm_rh_cols <- grep ( "^wm_rh_" , names ( d ) , value = TRUE )
  				wm_cols    <- c ( wm_lh_cols , wm_rh_cols )
  				wm_base    <- union ( sub ( "^wm_lh_" , "" , wm_lh_cols ) ,
                      					sub ( "^wm_rh_" ,  "" ,  wm_rh_cols ) )

  				# Keep rows with any WM value and coerce numerics
  				d <- d %>%
  						dplyr :: select ( PIDN , UnQID , DCDate , all_of ( wm_cols ) ) %>%
    				mutate ( across ( all_of (    wm_cols  ) ,
                  		~ suppressWarnings ( as.numeric ( .x ) ) ) ) %>%
					filter ( if_any (  all_of ( wm_cols )  , ~ !is.na ( .x ) ) )

  # Average L/R for each shared base parcel; if only one hemi present, use it
  	for ( b in wm_base ) {
			lh <- paste0 ( "wm_lh_" , b )
    		rh <- paste0 ( "wm_rh_" , b )
    		lh_exists <- lh %in% names ( d )
    		rh_exists <- rh %in% names ( d )
    	if ( lh_exists && rh_exists ) {
      d [[ b ]] <- rowMeans ( cbind ( d [[ lh ]] , d [[ rh ]] ) , na.rm = TRUE )
    	} else if ( lh_exists ) {
      d [[ b ]] <- d [[ lh ]]
    	} else if ( rh_exists ) {
      d [[ b ]] <- d [[ rh ]] } }

  # Longify the averaged parcels only
  d_long <- d %>%
    dplyr :: select ( PIDN , UnQID , DCDate ,  all_of ( wm_base ) ) %>%
    tidyr::pivot_longer( cols = all_of ( wm_base ) ,
                        names_to = "parcel", values_to = "value" ) %>%
    dplyr::filter ( !is.na ( value ) & !is.nan ( value ) ) %>%
    dplyr::mutate ( metric = metric_name )

  d_long } #end function noddi collapse
##########################################################################



# 4) DTI (ComBat): longify FA/MD white‑matter tracts
#    We keep UnQID + tract columns, attach PIDN, and longify.
 longify_dti_combat <- function ( df , metric_name ) {
  # Attach IDs and drop only date_diff (if present)
		d <- attach_pidn ( df ) %>%
    					dplyr::select ( - any_of ( "date_diff" ) )

  # Identify lateralized columns and their base names
  lat_cols <- grep ( "_(l|r)$" , names( d ) , value = TRUE )
  bases    <- unique ( sub ( "_(l|r)$" , "" , lat_cols ) )

  # Average L/R per base; if only one side exists, use it
	for ( b in bases ) {
    	l <- paste0 ( b , "_l" )
    	r <- paste0 ( b , "_r" )
    	has_l <- l %in% names ( d )
    	has_r <- r %in% names ( d )
    if ( has_l && has_r ) {
      	d [[ b ]] <- rowMeans ( cbind ( d [[ l ]] , d [[ r ]] ) , na.rm = TRUE )
    } else if (	has_l ) {
      	d [[ b ]] <- d [[ l ]]
    } else if ( has_r ) {
      d [[ b ]] <- d [[ r ]]
    } }#end for loop

  # IDs to keep if present
  id_cols <- intersect ( c ( "PIDN" , "UnQID" , "DCDate" ) , names ( d ) )
  # Final tract columns = everything except IDs and raw L/R columns
  tract_cols <- setdiff ( names ( d ) , c ( id_cols , lat_cols ) )

  d %>%
    dplyr::select ( all_of ( id_cols ) , all_of ( tract_cols ) ) %>%
    pivot_longer ( cols = all_of ( tract_cols ) ,
                        names_to = "parcel" , values_to = "value" ) %>%
    dplyr::filter ( !is.na ( value ) & !is.nan ( value ) ) %>%
    dplyr::mutate ( metric = metric_name )
 }#end longify dti combat function


# ── Build modality‑specific long tables ───────────────────────────────────────

# NODDI
 noddi_ficv_long <- collapse_noddi_wm_lr ( UCSF_db$imaging_noddi_ficv , "ficv" ) %>%
 						filter ( ! as.character ( UnQID ) %in% OutliersNODDI ) #Prune out outliers
 noddi_fiso_long <- collapse_noddi_wm_lr ( UCSF_db$imaging_noddi_fiso , "fiso" ) %>%
 						filter ( ! as.character ( UnQID ) %in% OutliersNODDI )
 noddi_odi_long  <- collapse_noddi_wm_lr ( UCSF_db$imaging_noddi_odi , "odi" ) %>%
 						filter ( ! as.character ( UnQID ) %in% OutliersNODDI )

 noddi_long <- bind_rows ( noddi_ficv_long , noddi_fiso_long , noddi_odi_long )

# DTI (ComBat)
fa_long <- longify_dti_combat ( UCSF_db$imaging_dti_combat_fa, "fa" ) %>%
						filter ( ! as.character ( UnQID ) %in% OutliersFA )
md_long <- longify_dti_combat ( UCSF_db$imaging_dti_combat_md, "md" ) %>%
						filter ( ! as.character ( UnQID ) %in% OutliersFA )
dti_long <- bind_rows ( fa_long , md_long )
# ── Combine if desired ───────────────────────────────────────────────────────
imaging_long_all <- bind_rows ( noddi_long , fa_long , md_long )

###W iden
library(tidyverse)

noddi_wide <- noddi_long %>%
 			 mutate ( parcel_metric = paste ( parcel , metric , sep = "." ) ) %>%  # make composite name
  					select ( PIDN , UnQID , DCDate , parcel_metric , value ) %>%
  				pivot_wider (
    						names_from = parcel_metric ,
   							values_from = value )

dti_wide <- dti_long  %>%
 			 mutate ( parcel_metric = paste ( parcel , metric , sep = "." ) ) %>%  # make composite name
  					select ( PIDN , UnQID , DCDate , parcel_metric , value ) %>%
  				pivot_wider (
    						names_from = parcel_metric ,
   							values_from = value )

#save(list = ls() , file = "LongWideNODDI_DWI.RData")