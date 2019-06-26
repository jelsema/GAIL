##
## This is the PACKAGE documentation
##

#' Geo-Assignment of Irregular Locations
#' 
#' @docType package
#' @name GAIL-package
#' @rdname GAIL-package
#' 
#' @description
#' A package for handling allocating a set of irregular spatial units (e.g. P.O. Box ZIP codes) 
#' onto a regular set of spatial units (e.g. Standard ZIP codes, Census Tract, County, etc).
#' 
#' 
#' @details
#' The motivating use-case for this package wasto allocate cases of some disease (e.g., influenza) 
#' from PO Box ZIP codes to nearby Standard ZIP codes. Often data from PO Box ZIP codes are discarded, 
#' or all of the data are aggregated to a larger spatial unit such as county. To preserve the cases 
#' associated with PO Box ZIP codes at the more fine spatial resolution (Standard ZIP code), we 
#' stochastically allocate cases from PO Box ZIP codes to 
#' 
#' @import methods
#' @importFrom utils globalVariables
#' 
"_PACKAGE"

globalVariables( strsplit( "ax ay case case_rate case_status data geometry index iregion mx my ncases npop nrural region rregion rural_rate rural_status rx ry sample_n sx sy type", split=" ")[[1]] )












