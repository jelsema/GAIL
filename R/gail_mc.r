


#' Monte Carlo Estimation of Variability from GAIL
#'
#' @description
#' For a specified statistic, performs a monte carlo simulation to obtain a distribution based on
#' repeated allocation using [gail]. 
#' 
#' 
#' @param sp_units Data frame of regular spatial units.
#' @param cases Data frame containing the the cases.
#' @param suid Column name in both `cases` which contains the spatial unit identification / name.
#' @param num_cases Column name from `cases` containing the number of cases, see details.
#' @param max_dist The maximum distance at which two locations can be considered neighbors.
#' @param RAP Function to be used for calculating assignment probabilities, see details.
#' @param unit_value The units of distance, default to meters. See [units::set_units].
#' @param spat_fun A function to compute the statistic of interest over the spatial units. Can also accept
#'        "moran", "localmoran", and "localmoran_full".
#' @param nsim The number of Monte Carlo replications to run.
#' @param seed If given, sets the seed for the RNG.
#' @param ... Space for additional arguments (e.g., for [gail]).
#' 
#' @details 
#' Most of the arguments are passed to the [gail] arguments of the same names. They are included 
#' (instead of using `...`) to facilitate understanding this function.
#' 
#' The `spat_fun` argument performs some computation on the dataset. Currently built-in options are:
#' - `spat_fun = "moran"` will compute Moran's I using `spdep::moran`.
#' - `spat_fun = "localmoran"` will compute the local Moran's I for each spatial unit using 
#'               `spdep::localmoran` and return the vector of local Moran's I.
#' - `spat_fun = "localmoran_full"` will return the full results of `spdep::localmoran`.
#' 
#' The monte carlo simulation repeatedly re-allocates the cases to the regular 
#' spatial units. For each allocation, the result from the argument `spat_fun` is saved. This could
#' be a scaler (e.g., Moran's I), a vector (e.g. local Moran's I), or something more complex. 
#' 
#' A custom function for `spat_fun` should take the arguments:
#'  - `dat_values` is a vector of numeric data representing the value (number of cases) for each spatial unit
#'  - `nghb_list` is the neighborhood matrix, e.g. as returned by `spdep::nb2listw`
#'  - `...` for extra arguments.
#' For example, the code of `spat_fun` for Moran's I is: yyy
#' 
#' ```
#'    spat_fun <- function( dat_values, nghb_list, ... ){
#'        spdep::moran( x     = dat_values,
#'                      listw = nghb_list,
#'                      n     = length(dat_values),
#'                      S0    = spdep::Szero(nghb_list), ...)$I
#'     }
#' ```
#' 
#' @return
#' The return format depends on the nature of the output of `spat_fun`.
#' - If `spat_fun` returns a scaler, `gail_mc` returns a vector.
#' - If `spat_fun` returns a vector, `gail_mc` returns an array.
#' - If `spat_fun` returns an array, `gail_mc` returns a list.
#' 
#' @seealso
#' [gail]
#' 
#' 
#' @importFrom magrittr %>%
#' @import units
#' @importFrom spdep localmoran moran nb2listw poly2nb Szero
#' 
#' 
#' @export
#' 

gail_mc <- function( sp_units, cases, suid, num_cases, max_dist, 
                     RAP=gail_rap, unit_value="m", spat_fun, nsim=10, seed=NULL, ... ){
  
  ########################################
  ## Check and/or make the function from some choices
  
  if( missing(spat_fun) ){
    spat_fun <- "moran"
  }
  
  if( is.character(spat_fun) ){
    
    if( spat_fun=="moran" ){
      spat_fun <- function( dat_values, nghb_list, ... ){
        spdep::moran( x     = dat_values,
                      listw = nghb_list, 
                      n     = length(dat_values),
                      S0    = spdep::Szero(nghb_list), ...)$I
      }
    } else if( spat_fun=="localmoran" ){
      spat_fun <- function( dat_values, nghb_list, ... ){
        as.data.frame( spdep::localmoran( 
          x     = dat_values,
          listw = nghb_list 
        ))$Ii
      }
      
    } else if( spat_fun=="localmoran_full" ){
      spat_fun <- function( dat_values, nghb_list, ... ){
        as.data.frame( spdep::localmoran( 
          x     = dat_values,
          listw = nghb_list 
        ))
      }
      
    }
    
    
  }
  
  ########################################
  ## Extract information from arguments
  nhbd       <- spdep::poly2nb( sp_units )
  Wlist      <- spdep::nb2listw( nhbd , style="W" )
  return_obj <- list()
  
  ########################################
  ## Run GAIL repeatedly
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  vector_of_seeds <- round( runif( nsim, 1, 1e9 ) )
  
  ## Potentially replace with parallel: foreach() %dopar% { ... }
  ## Add an argument ncores
  for( ii in 1:nsim ){
    
    gail01a <- gail( sp_units   = sp_units  ,
                     cases      = cases     , 
                     suid       = suid      ,
                     num_cases  = num_cases ,
                     max_dist   = max_dist  , ## distance in meters
                     seed       = vector_of_seeds[ii],
                     ...
    )
    
    data_vec <- gail01a[["sp_units"]][[num_cases]]
    return_obj[[ii]] <- spat_fun( data_vec, Wlist )
    
  }
  
  
  ########################################
  ## Collapse the results storage appropriately
  
  dim_sfo <- dim( return_obj[[1]] )
  if( is.null(dim_sfo) ){
    if( length( return_obj[[1]] )==1 ){
      return_obj <- unlist( return_obj )
    } else{
      spat_fun_output <- t( sapply( return_obj, FUN="cbind" ) )
    }
  } else{
    return_obj <- list()
  }
  
  
  return( return_obj )
  
}