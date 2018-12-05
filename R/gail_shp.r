



#' Geo-Assignment of Irregular Locations to Shapes
#'
#' @description
#' DESCRIPTION OF FUNCTION
#' 
#' 
#' @param sp_units Data frame of regular spatial units, see details.
#' @param cases Data frame containing the the cases, see details.
#' @param num_cases Column name from `cases` containing the number of cases, see details.
#' @param group Optional name of variable in `cases` to use for filtering. See details.
#' @param label Optional value of `group` by which to filter. See details.
#' @param convert Logical, convert returned object to `Spatial*` type objects (from **sp**).
#' @param ... Space for additional arguments.
#' 
#' 
#' 
#' @details 
#'  - `cases` : The dataset of cases should be aggregated and have a column for the number 
#'     of cases. If the column given by `num_cases` does not exist in `cases`,
#'     then `cases` will be aggregated to create the column.
#'  
#'  - `group` and `label` : If provided, these arguments are used to filter `cases`
#'    to rows where `group == label`.
#' 
#' @note
#' MAYBE SOME COMMENTS.
#' 
#' 
#' @return
#' Returns a list containing DESCRIBE THE RESULTS
#' 
#' 
#' @seealso
#' [gail]
#' 
#' @export
#' 
gail_shp <- function( sp_units, cases, num_cases, group, label , convert=FALSE , ... ){
  
  
  ## Filter the cases if needed
  if( !missing(group) & !missing(label) ){
    cases <- cases %>% filter( !!sym(group) == label )
  }
  
  ## Aggregate cases if needed
  if( !any( colnames(cases) == num_cases ) ){
    cases <- cases %>% 
      dplyr::group_by( !!sym( suid ) ) %>% 
      dplyr::summarize( !!num_cases := n() )
  }
  
  
  ## Check to see if the data are sf objects
  ## If not, set them to be such
  if(  !("sf" %in% class(sp_units))  ){
    sp_units <- st_as_sf( sp_units ,  coords = c("longitude", "latitude"), 
                          agr = "aggregate", ... )
  }
  shapefile_crs_format <- st_crs( sp_units )
  
  
  if(  !("sf" %in% class(cases))  ){
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










