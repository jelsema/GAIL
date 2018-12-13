



#' Geo-Assignment of Irregular Locations to Shapes
#'
#' @description
#' Allocates the number of cases/events from irregular spatial units onto to regular set of spatial units.
#' The regular spatial units are assumed to have POLYGON geometry, while the irregular spatial units are
#' assumed to have POINT geometry. 
#' 
#' Note that this function is largely a wrapper for using `sf::st_contains` with basic
#' aggregation before and after. It is provided primarily for (relative) consistency with [gail].
#' 
#' 
#' @param sp_units Dataset of regular spatial units, see details.
#' @param cases Dataset containing the the cases, see details.
#' @param suid Column name in both `sp_units` and `cases` which contains the spatial unit identification / name.
#' @param num_cases Column name from `cases` containing the number of cases, see details.
#' @param group Optional name of variable in `cases` to use for filtering. See details.
#' @param label Optional value of `group` by which to filter. See details.
#' @param convert Logical, convert returned object to `Spatial*` type objects (from **sp**).
#' @param ... Space for additional arguments.
#' 
#' 
#' 
#' @details 
#'  - `sp_units` should be an object of class `sf` (from the **sf** package). The 
#'  argument `convert` can be used to convert the results into a `Spatial*` object (from 
#'  the **sp** package).
#'  
#'  - `cases` contains the units (and associated number of cases) to be allocated. 
#'  This should be an object of class `sf`. If not, it must have coordinates `latitude` 
#'  and `longitude` and will be converted to an `sf` object.
#'  The dataset of cases should be aggregated and have a column for the number 
#'     of cases. If the column given by `num_cases` does not exist in `cases`,
#'     then `cases` will be aggregated by `suid` to create the column.
#'  
#'  - `group` and `label` : If provided, these arguments are used to filter `cases`
#'    to rows where `group == label`.
#' 
#' 
#' @return
#' Returns a copy of `sp_units` after having added a column named `num_cases`, which contains 
#' the number of cases which fall within the boundry of each spatial unit.
#' 
#' @import sf
#' @importFrom data.table :=
#' 
#' @seealso
#' [gail]
#' 
#' @export
#' 
gail_shp <- function( sp_units, cases, suid, num_cases, group, label , convert=FALSE , ... ){
  
  
  ## Filter the cases if needed
  if( !missing(group) & !missing(label) ){
    if( !is.null(group) & !is.null(label)  ){
      cases <- cases %>% filter( !!sym(group) == label )
    }
  }
  

  ## Check to see if the data are sf objects
  ## If not, set them to be such
  if(  !("sf" %in% class(sp_units))  ){
    # sp_units <- st_as_sf( sp_units ,  coords = c("longitude", "latitude"),  agr = "aggregate", ... )
    stop( "sp_units does not have class sf and gail_shp() cannot detect boundaries.")
    #  Add boundaries or use gail() for POINT (non-polygon) geometries.
  } else{
    if(  any( sf::st_geometry_type( sp_units ) != "POLYGON" ) ){
      stop("gail_shp() expects features of sp_units to have POLYGON geometry.")
      # See ?sf::st_set_geometry and ?sf::st_as_sf to set appropriate geometry.
      # gail() expects features of sp_units to have POINT geometry.
    }
  }
  shapefile_crs_format <- st_crs( sp_units )
  
  
  if(  !("sf" %in% class(cases))  ){
    
    ## Aggregate cases if needed
    if( !any( colnames(cases) == num_cases ) ){
      cases <- cases %>% 
        dplyr::group_by( !!sym( suid ) ) %>% 
        dplyr::summarize( !!num_cases := n() )
    }
    
    cases <- st_as_sf( cases , coords = c("longitude", "latitude"),
                          crs = shapefile_crs_format  , agr = "aggregate" )
  } else{
    st_crs(cases) <- shapefile_crs_format
  }
  
  
  ## Sum up the cases for each of the spatial units
  idx_list <- st_contains( sp_units, cases )
  
  sp_units[[ num_cases ]] <- 0
  for( ii in 1:length(idx_list ) ){
    tidx <- idx_list[[ii]] 
    if( length(tidx) > 0 ){
      sp_units[[ num_cases ]][ii] <- sum( cases[[ num_cases ]][tidx] , na.rm=TRUE )
    }
  }
  
  
  ## Convert from sf to sp if requested
  if( convert==TRUE ){
    sp_units <- as( sp_units , "Spatial" )
  }
  
  ## Return results
  return( sp_units )
  
}










