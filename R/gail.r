


#' Geo-Assignment of Irregular Locations
#'
#' @description
#' Stochastic allocation of irregular spatial units to nearby regular spatial units. Both sets of 
#' of spatial units are assummed to be represented by their centroid (hence they should be POINT objects,
#' and not POLYGON objects). Several methods are available for allocating irregular units, including equal 
#' probability, inverse centroid distance, and use of an index variable. Alternatively the user may 
#' specify a custom function to determine allocation probabilities.
#' 
#' @param sp_units Data frame of regular spatial units, see details.
#' @param cases Data frame containing the the cases, see details.
#' @param suid Column name in both `cases` and `sp_units` which contains the spatial unit identification / name.
#' @param num_cases Column name from `cases` containing the number of cases, see details.
#' @param max_dist The maximum distance at which two locations can be considered neighbors.
#' @param RAP Function to be used for calculating assignment probabilities, see details.
#' @param seed If given, sets the seed for the RNG.
#' @param unit_value The units of distance, default to meters. See [units::set_units].
#' @param convert Logical, convert returned object to `Spatial*` type objects (from **sp**).
#' @param ... Space for additional arguments (e.g., for `RAP`).
#'  
#' 
#' @details 
#'  - `sp_units` is the set of spatial units that are the "destination", the cases will be aggregated to
#'    this set of units. It should not contain duplicate spatial units, and should not have
#'    any spatial units on a different scale. It must either be an object of class `sf` with POINT or 
#'    POLYGON geometry, or it must have variables `longitude` and `latitude` to be convertable into 
#'    an object of class `sf`.
#'  
#'  - `cases` contains the units (and associated number of cases) to be allocated. 
#'  As with `sp_units`, this must either already have class `sf` or have `latitude` 
#'  and `longitude` to be converted to an `sf` object.
#'  
#'  The dataset of cases should be aggregated and have a column for the number of cases. 
#'  If the column given by `num_cases` does not exist in `cases`, then `cases` will be 
#'  aggregated by `suid` to create the column.
#'  
#'  - Cases will first be allocated deterministically by exact matches with `suid`.
#'  Any remaining cases will be allocated stochastically.  
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
#' 
#' @seealso 
#' [gail_mc]
#' 
#' @return
#' Returns a list containing four elements:
#' - `rap_method`: The method of stochastic allocation
#' - `units_reg`: A copy of the input `sp_units` with the number of cases detrministically allocated.
#' - `units_irr`: The irregular spatial units and the number of cases which were stochastically allocated 
#'                from these units to the regular spatial unit.
#' - `sp_units`: A copy of the input `sp_units` along with the number of cases (following allocation).
#' 
#' 
#' @importFrom data.table :=
#' @importFrom dplyr group_by summarize sym
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom sf st_as_sf st_crs st_distance st_geometry_type
#' @import units
#' 
#' @export
#' 

gail <- function( sp_units, cases, suid, num_cases, max_dist, RAP=gail_rap, seed=NULL, unit_value="m", convert=FALSE, ... ){
  
  
  ## Some error-checking on arguments
  if(   missing(sp_units) | missing(cases)   ){
    stop("Missing at least one of: sp_units, cases")
  }
  if(   missing(max_dist)  |  !is.numeric(max_dist)   ){
    stop("max_dist must be numeric")
  }
  
  ## Attach units to max_dist
  ## Eventually this can probably be removed
  ## This line doesn't work for some reason
  #   max_dist <- units::as_units( max_dist , value=unit_value )
  
  ## Should first check if max_dist has a unit attached.
  ## If yes, then use that. Otherwise, set it
  units( max_dist ) <- unit_value
  
  
  ##
  ## Extract some information from the arguments
  ##
  
  cc <- match.call()
  
  if( !is.function(RAP) ){
    ## if RAP is not a function, set RAP to 
    ## be gail_rap and set a default method
    if( any(colnames(sp_units)==RAP)  ){
      method    <- "index"
      index_val <- colnames(sp_units)[ which(colnames(sp_units)==RAP) ]
    } else if( RAP == "icd") {
      method    <- "icd"
      index_val <- NULL
    } else{
      method    <- "equal"
      index_val <- NULL
    }
    RAP <- gail_rap
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
    sp_units <- sp_units %>%
      sf::st_as_sf( coords = c("longitude", "latitude"), agr = "aggregate", ...  )
  } else{
    if(  !any( sf::st_geometry_type( sp_units ) %in% c("POINT", "POLYGON") ) ){
      mssg1 <- "Features of sp_units must have POINT or POLYGON geometry."
      mssg2 <- "\n       See ?sf::st_set_geometry and ?sf::st_as_sf to set appropriate geometry."
      stop( paste0(mssg1,mssg2) )
    }
    # Recommend this code?
    # sp_units <- st_set_geometry( sp_units, value=NULL )
    # sp_units <- st_as_sf( sp_units , coords = c("longitude", "latitude"), agr = "aggregate")
    # sp_units <- sp_units %>%  dplyr::filter( complete.cases(.) )
  }
  
  shapefile_crs_format <- sf::st_crs( sp_units )
  
  
  ## Detect if cases is simple features object
  ## If NOT, set as sf
  ## If YES, check to make sure it is POINT geometry
  if( !("sf" %in% class(cases))  ){
    
    ## Aggregate cases if needed
    if( !any( colnames(cases) == num_cases ) ){
      cases <- cases %>% 
        dplyr::group_by( !!sym( suid ) ) %>% 
        dplyr::summarize( !!num_cases := n() )
    }
    
    cases <- sf::st_as_sf( cases , coords = c("longitude", "latitude"),
                       crs = shapefile_crs_format  , agr = "aggregate" )
    
  } else{
    if(  !any( sf::st_geometry_type( sp_units ) %in% c("POINT", "POLYGON") ) ){
      mssg1 <- "Features of cases must have POINT or POLYGON geometry."
      mssg2 <- "\n       See ?sf::st_set_geometry and ?sf::st_as_sf to set appropriate geometry."
      stop( paste0(mssg1,mssg2) )
    }
    sf::st_crs(cases) <- shapefile_crs_format
    
  }
  
  ## Is num_cases an existing column in sp_units
  ## If NOT, then initialize it.
  if(  !(any(colnames(sp_units) == num_cases))  ){
    sp_units[[ num_cases ]] <- 0
  }
  
  ##
  ## Deterministic Allocation
  ##
  
  ## Is suid a column in BOTH datasets?
  ## If YES, then add to num_cases for those spatial units.
  ## Also REMOVE from cases those values which are de
  if(  any(colnames(sp_units)==suid)  ){
    
    allocated_cases <- rep( 0 , nrow(cases) )
    
    for( ii in 1:nrow(cases) ){
      idx_suid <- which( sp_units[[suid]] == cases[[suid]][ii]  )
      
      if( length(idx_suid)==1 ){
        sp_units[[num_cases]][idx_suid] <- sp_units[[num_cases]][idx_suid] + cases[[num_cases]][ii]
        allocated_cases[ii] <- 1
      }
    }
    
    ## Generate the sets of cases allocated deterministically
    ## and those needed to be allocated stochastically
    cases_detrm <- sp_units
    cases_stoch <- cases[ allocated_cases==0 , ]
  } else{
    
    ## Maybe delete these? Is it really necessary to report this?
    ## Would need to rework object names below, but would save memory for large datasets
    cases_detrm <- sp_units
    cases_stoch <- cases
    
  }
  
  ##
  ## Stochastic Allocation
  ##
  
  ## Only stochastically allocate if there are any spatial units left to allocate
  if( nrow(cases_stoch)>0 ){
    
    ## Number of neighbors for the irregular regions
    num_neighbors <- map_dfr( cases_stoch[[ suid ]] ,
                              ~gail_nn( .x, df1=cases_detrm, df2=cases_stoch, 
                                        max_dist=max_dist, suid=suid, ... ) )
    
    if(  min(num_neighbors) < 3  ){
      stop("Some irregular sp_units have less than 3 regular unit neighbors. Try increasing max_dist")
    }
    
    ## Loop through all of the IRREGULAR sp_units and apply RAP()
    for( ii in 1:nrow(cases_stoch) ){
      
      n_to_assign <- cases_stoch[[ num_cases ]][ii]
      
      ## The set of regular sp_units within max_dist of the cases to be allocated
      # dist_kk <- units::set_units(  sf::st_distance( sp_units , cases_stoch[ii,] ),   value=unit_value )
      dist_kk <- sf::st_distance( sp_units , cases_stoch[ii,] )
      units( dist_kk ) <- unit_value
      
      
      idx_kk  <- which( dist_kk <= max_dist )
      reg_kk  <- sp_units[ idx_kk , ]
      
      
      allocation_probs <- RAP( rUnits=sp_units, iUnits=cases_stoch[ii,]   , 
                               rUclose=reg_kk  , rUdist=dist_kk[idx_kk],  max_dist=max_dist  ,
                               method=method   , index_val=index_val  )
      
      allocations <- sample( 1:length(idx_kk) , n_to_assign, replace=TRUE, prob=allocation_probs )
      n_assigned  <- table( c( 1:length(idx_kk) , allocations )  ) - 1
      sp_units[[ num_cases ]][ idx_kk ] <- sp_units[[ num_cases ]][ idx_kk ] + n_assigned
      
    }
    
  }

  ## Convert from sf to sp if requested
  if( convert==TRUE ){
    cases_detrm <- as( cases_detrm , "Spatial" )
    cases_stoch <- as( cases_stoch , "Spatial" )
    sp_units    <- as( sp_units    , "Spatial" )
  }
  
  
  ## Make and return the results
  return_obj <- list(
    rap_method = list( method=method, index_val=index_val ),
    units_reg  = cases_detrm ,
    units_irr  = cases_stoch ,
    sp_units   = sp_units
  )
  
  return( return_obj )
  
  
}





















