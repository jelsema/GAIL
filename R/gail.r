


#' Geo-Assignment of Irregular Locations
#'
#' @description
#' Applies random allocation of irregular spatial units to nearby regular spatial units. Both sets of 
#' of spatial units are assummed to be represented by their centroid (hence they should be point objects,
#' and not shape objects). Several methods are available for allocating irregular units, including equal 
#' probability, inverse centroid distance, and use of an index variable. Alternatively the user may 
#' specify a custom function to determine allocation probabilities.
#' 
#' @param sp_units Data frame of regular spatial units, must have columns named `longitude` and `latitude`.
#' @param cases Data frame containing the the cases, see details.
#' @param suid Column name in both `sp_units` and `cases` which contains the spatial unit identification / name.
#' @param num_cases Column name from `cases` containing the number of cases, see details.
#' @param max_dist The maximum distance at which two locations can be considered neighbors.
#' @param group Name of variable in `sp_units` to denote regular and irregular spatial units.
#' @param labels Vector of values for the regular and irregular units (the values from `group`), with the first.
#'         element being the regular units, and the second element being the irregular units.
#' @param RAP Function to be used for calculating assignment probabilities, see details.
#' @param seed If given, sets the seed for the RNG.
#' @param convert Logical, convert returned object to `Spatial*` type objects (from **sp**).
#' @param ... Space for additional arguments (e.g., for `RAP`).
#'  
#' 
#' @details 
#'  - `cases` : The dataset of cases should be aggregated and have columns for the spatial unit name 
#'  (given by `suid`) and the number of cases. If the column given by `num_cases` does not exist in `cases`,
#'   then `cases` will be aggregated to create the column.
#'  
#'  - `RAP` : The argument controls how irregular spatial units are allocated to regular spatial units.
#'  There are several formats this can take. If `RAP` is a function, it will be used. Otherwise `gail_rap` 
#'  is used. By default the internal function `gail_rap` is used with arguments `method` and `index_var`
#'  (see documentation of [gail_rap]) determined by:
#'  
#'    - If `RAP` is a column name from `sp_units` then `method="index"` and `index_var = RAP` .
#'    - If `RAP="icd"` then `method="icd"` and `index_var = NULL` .
#'    - Otherwise, `method="equal"` and `index_var = NULL` .
#'  
#'  The user may also specify only `method` and/or `index_val` to control the behavior of `gail_rap` appropriately.
#'  
#'  - `labels` : The general (but not default) format of this is `labels=c("Standard", "Irregular")`. If `labels`
#'   is missing, `gail` finds the unique values of `group` and takes the first element as the regular spatial 
#'   units, and the second element as the irregular spatial units.
#' 
#' 
#' 
#' @return
#' Returns a list containing DESCRIBE THE RESULTS
#' 
#' 
#' @import dplyr 
#' @import sf
#' @import units
#' @importFrom purrr map_dfr
#' 
#' @export
#' 
gail <- function( sp_units, cases, suid, num_cases, max_dist, group, labels, RAP=gail_rap, seed=NULL, convert=FALSE, ... ){
  
  ## Some error-checking on arguments
  if(   missing(sp_units) | missing(cases)   ){
    stop("Missing at least one of: zips, cases")
  }
  if(   missing(max_dist)  |  !is.numeric(max_dist)   ){
    stop("max_dist must be numeric")
  }
  if( is.null(group) ){
    stop("group must be a variable name in sp_units identifying regular / irregular units")
  }
  
  ##
  ## Extract some information from the arguments
  ##
  
  cc <- match.call()
  
  if( missing(labels) ){
    labels <- unique( sp_units[[ suid ]] )[1:2]
  }
  
  if( !is.function(RAP) ){
    ## if RAP is not a function, set RAP to 
    ## be gail_rap and set a default method
    if(   any(colnames(sp_units)==RAP)  ){
      method    <- "index"
      index_val <- colnames(sp_units)[ which(colnames(sp_units)==RAP) ]
    } else if( RAP == "icd") {
      method    <- "icd"
      index_val <- NULL
    } else{
      method    <- "equal"
      index_val <- NULL
    }
    RAP    <- gail_rap
  } else{
    method    <- cc$method
    index_val <- cc$index_val
    if( is.null(method) & any(colnames(sp_units)==index_val) ){
      method <- "index"
    }
    
  }
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  
  ## Filter and separate spatial units
  ## Remove any units with NA in them
  
  if( !("sf" %in% class(sp_units))  ){
    sp_units <- sp_units %>%  dplyr::filter( complete.cases(.) ) %>%
      st_as_sf( coords = c("longitude", "latitude"), agr = "aggregate", ...  )
  } else{
    sp_units <- sp_units %>%  dplyr::filter( complete.cases(.) )
  }
  
  units_reg <- sp_units %>% dplyr::filter( !!sym( group ) == labels[1] )
  units_irr <- sp_units %>% dplyr::filter( !!sym( group ) == labels[2] )
  
  if( !is.null( cc$crs) ){
    units( max_dist ) <- units( st_distance( units_reg[1,], units_reg[2,] ) )
  }
  
  
  
  #units_reg <- sp_units %>% dplyr::filter( complete.cases(.) & !!sym( group ) == labels[1] ) %>%
  #  st_as_sf( coords = c("longitude", "latitude"), agr = "aggregate", ... )
  #
  #units_irr <- sp_units %>% dplyr::filter( complete.cases(.) & !!sym( group ) == labels[2] ) %>%
  #  st_as_sf( coords = c("longitude", "latitude"), agr = "aggregate", ... )
  #
  # units_oth <- sp_units %>% dplyr::filter( !(!!sym(group) %in% labels) & !is.na(!!sym(group)) ) %>%
  #   st_as_sf( coords = c("longitude", "latitude"), agr = "aggregate", crs=crs )
  
  ## Number of neighbors for the irregular regions
  num_neighbors <- map_dfr( units_irr[[ suid ]] ,
                            ~gail_nn( .x, df1=units_reg, df2=units_irr, max_dist=max_dist, suid=suid ) )
  
  if(  min(num_neighbors) < 3  ){
    stop("Some irregular sp_units have less than 3 regular unit neighbors. Try increasing max_dist")
  }
  
  ## Get number of cases for each suid / unit and join to main dataframes
  
  if( !any( colnames(cases) == num_cases ) ){
    
    cases <- cases %>% 
      dplyr::group_by( !!sym( suid ) ) %>% 
      dplyr::summarize( !!num_cases := n() )
    
  }
  
  
  
  #cases_sum <- cases %>% 
  #  dplyr::group_by( !!sym( suid ) ) %>% 
  #  dplyr::summarize( ncase = n() )
  
  ## The final sets of regular, irregular, and "other" regions
  units_reg <- units_reg %>% dplyr::left_join( cases , by=suid )
  units_irr <- units_irr %>% dplyr::left_join( cases , by=suid ) %>% dplyr::filter( !!sym( num_cases ) > 0  )
  units_oth <- cases     %>% dplyr::anti_join( units_reg , by=suid ) %>% dplyr::anti_join( units_irr , by=suid )
  
  ## Loop through all of the IRREGULAR sp_units and apply RAP()
  assigned_units <- list()
  for( ii in 1:nrow(units_irr) ){
    
    n_to_assign <- units_irr[[ num_cases ]][ii]
    units_reg   <- units_reg %>% dplyr::mutate(
      dist_to_irr = st_distance( units_reg , units_irr[ii,] )
    )
    
    ## The set of regular sp_units within max_dist of the irregular
    reg_kk <- units_reg %>% dplyr::filter( dist_to_irr <= max_dist )
    
    allocation_probs <- RAP( rUnits=units_reg, iUnits=units_irr[ii,]   , 
                             rUclose=reg_kk  , max_dist=max_dist  ,
                             method=method   , index_val=index_val, ... )
    
    assigned_units[[ii]] <- sample( reg_kk[[ suid ]], n_to_assign,
                                    replace=TRUE, prob=allocation_probs )
  }
  
  ## Convert to vector and summarize into counts
  assigned_units <- tibble( !!suid := unlist( assigned_units ) ) %>%
    dplyr::group_by( !!sym(suid) ) %>%
    dplyr::summarize( newcases = n() )
  
  units_reg <- units_reg %>% dplyr::select( -dist_to_irr )
  
  ## Merge the counts into the data
  final_assignment <- units_reg %>% 
    dplyr::left_join( assigned_units , by=suid  ) %>%
    dplyr::mutate( !!num_cases :=  !!sym( num_cases ) + newcases ) %>%
    dplyr::select( -newcases )
  
  
  ## Convert from sf to sp if requested
  if( convert==TRUE ){
    units_reg <- as( units_reg , "Spatial" )
    units_irr <- as( units_irr , "Spatial" )
    units_oth <- as( units_oth , "Spatial" )
    final_assignment <- as( final_assignment , "Spatial" )
  }
  
  
  ## Make and return the results
  return_obj <- list(
    rap_method       = list( method=method, index_val=index_val ),
    units_reg        = units_reg ,
    units_irr        = units_irr ,
    units_oth        = units_oth ,
    final_assignment = final_assignment
  )
  
  return( return_obj )
  
  
  if( FALSE ){
    # !!sym( num_cases ) 
    # - from -
    # ncase
  }
  
  
}










