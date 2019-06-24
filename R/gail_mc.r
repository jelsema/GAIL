


#' Monte Carlo Estimation of Variability from GAIL
#'
#' @description
#' 
#' 
#' @param data the data ... replace me  
#' 
#' @details 
#' 
#' 
#' @return
#' Returns a list containing DESCRIBE THE RESULTS
#' 
#' 
#' 
#' @export
#' 
#' 
#' 
#' 
#' 

gail_mc <- function( data ){
  
  ## GAIL should NOT REQUIRE either set of spatial units to be POINT.
  
  ## Some error-checking on arguments
  if(   missing(sp_units) | missing(cases)   ){
    stop("Missing at least one of: sp_units, cases")
  }
  if(   missing(max_dist)  |  !is.numeric(max_dist)   ){
    stop("max_dist must be numeric")
  }
  
  
  
}