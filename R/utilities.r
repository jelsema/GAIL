
################################################################################
##  
##  Utility functions for the Geo-assignment of Irregular Locations project
##  
################################################################################




#' Number of Neighbors
#'
#' @description
#' Calculates the the number of neighbors for a spatial unit.
#' 
#' @param val Name of spatial unit from df2.
#' @param df1 Data frame of regular spatial units, must have columns `latitude` and `longitude`.
#' @param df2 Data frame of irregular spatial units, must have columns named `latitude` and `longitude`.
#' @param max_dist The maximum distance at which two locations can be considered neighbors.
#' @param suid The column name which contains the ID of the spatial unit.
#' @param ... Space for additional arguments, passed to [sf::st_distance].
#' 
#' @details
#' This is an internal function designed to be used within [gail]
#' as part of a purrr type operation.
#' 
#' 
#' @return
#' Returns the number of spatial units from `df1` that are within
#' distance `max_dist` of the unit `val` from `df2`.
#' 
#' @seealso
#' [gail]
#' 
#' @import units
#' @importFrom sf st_distance
#' 
gail_nn <- function( val, df1, df2, max_dist, suid, ... ){
  
  # df2b     <- df2 %>% filter( !!sym( suid ) == val )
  ridx <- which( df2[[suid]] == val )
  df2b <- df2[ ridx , ]
  
  dist_vec <- sf::st_distance( df1 , df2b, ... )
  units(dist_vec) <- units(max_dist)
  
  return( list(  nneigh = sum(dist_vec <= max_dist) )  )
  
}





#' Region Assignment Probability
#'
#' @description
#' Function to calculate the (relative) probability of assigning an irregular unit to a
#' neighboring regular unit.
#' 
#' @param rUnits Data frame of regular spatial units.
#' @param iUnits Data frame of irregular spatial unit being allocated.
#' @param rUclose The set of regular units that are considered close (within `max_dist`)
#' @param rUdist Vector of distances of `rUclose` to `iUnits`
#' @param max_dist The maximum distance at which two locations can be considered neighbors.
#' @param method The method of allocation, see details.
#' @param index_val If `method="index"`, the name of the variable (in `rUnits`) to be used as an index variable.
#' 
#' @details 
#'  - `method = "index"` : the variable given by `index_val` (in `rUnits`) is used in the form `X / sum(X)`. 
#'  - `method = "equal"` : all neighboring spatial units (in `rUclose`) are given equal probability.
#'  - `method = "icd"`   : inverse centroid distance, the closer the centroid of a regular spatial unit,
#'  the higher the probability. If `D` is the distance to the centroid, the probabilities are given by:
#'  `(1/D) / sum(1/D) = 1 / ( D * sum(1/D) )`. Similar to the index method, but uses a computed variable.
#' 
#' @note
#' This is an internal function designed to be used within [gail]. This function
#' provides several default options for the `RAP` argument to `gail`. Users can alternatively 
#' supply a custom function as the argument for the `RAP` argument in `gail`, provided it accepts
#' the same arguments and the returns is the same type and dimension.
#' 
#' The returned value is passed to `sample` as the argument `prob`.
#' 
#' 
#' @return
#' Returns a set of weights which has length equal to the number of spatial units in `rUclose`.
#' 
#' @seealso
#' [gail]
#' 
#' @export
#' 
gail_rap <- function( rUnits, iUnits, rUclose, rUdist, max_dist, method, index_val ){
  
  if( method=="index" ){
    allocation_probs <- rUclose[[index_val]] / sum(rUclose[[index_val]])
  } else if( method=="equal" ){
    allocation_probs <- rep( 1 / nrow(rUclose), nrow(rUclose) )
  } else if( method=="icd" ){
    allocation_probs <- 1 / c( as.numeric(rUdist) )
  }
  
  return( allocation_probs )
  
}




