#' Generate Artificial Polygon Regions
#'
#' @description
#' Create a set of spatial polygon regions 
#' 
#' @param npoints Number of internal points to select
#' @param type Type of spatial units to generate, accepts "regular" and "irregular". See details.
#' @param nedge Data frame containing the the cases, see details.
#' @param seed If given, sets the seed for the RNG.
#' @param suid Optional spatial unit ID, used as prefix for generated spatial units. 
#' @param ... Space for additional arguments (e.g., for `fields::cover.design`).
#'  
#' 
#' @details 
#'  The `type` argument works as follows:
#'  - `type="regular"` produces a regular grid of units with `nedge` units along each dimension.
#'  - `type="irregular"` randomly draws 1000 points uniformly across a 100x100 region (excluding a band of 
#'  width 5 around the edge). Then uses `fields::cover.design` to select `npoints` of them.
#'  Creates regions using Voronoi tesselation via `deldir::deldir` and `deldir::tile.list`.
#'  
#' 
#' @examples
#' \dontrun{
#'   loca_reg <- gail_gen_regions( npoints=40, type="regular"  , nedge=10, suid="reg" )
#'   loca_irr <- gail_gen_regions( npoints=40, type="irregular", nedge=6 , seed=42 , P=-20, Q=20 )
#'   
#'   ggplot( loca_reg ) + 
#'     geom_sf( aes(fill=region) )
#'   ggplot( loca_irr ) +
#'     geom_sf( aes(fill=region) )
#' }
#' 
#' 
#' @importFrom dplyr bind_rows
#' @import sf
#' @importFrom fields cover.design
#' @import deldir
#' @importFrom stringr str_pad
#' @importFrom stringr str_length
#' 
#' @export
gail_gen_regions <- function( npoints, type="irregular", nedge=5, seed=NULL, suid=NULL, ... ){
  
  ## Reference to create SF polygons:
  ## https://stackoverflow.com/questions/44335246/polygons-from-coordinates
  
  if( is.numeric(seed) ){ set.seed( seed ) }
  if( is.null(suid) ){ suid <- "suid" }
  
  cc <- match.call(expand.dots = TRUE)
  
  if( type == "regular" ){
    
    make_row <- function(x){ 
      st_polygon( list(
        x + cbind( rep(0,5), (ii-1)*rep(100/nedge,5) ) 
      ))
    } 
    
    edge_seq <- seq(0, 100, 100/nedge)
    nlp      <- length(edge_seq) - 1
    base_box <- as.matrix(
      cbind( edge_seq[ c(1, 2, 2, 1, 1)  ],
             edge_seq[ c(1, 1, 2, 2, 1)  ])
    )
    
    base_row <- list()
    for( ii in 1:nedge ){
      base_row[[ii]] <- base_box + cbind( (ii-1)*rep(100/nedge,5), rep(0,5) )
    }
    
    sim_shp00 <- list()
    for( ii in 1:nedge ){
      sim_shp00[ (1:nedge) + (ii-1)*nedge ] <- lapply( base_row, FUN=make_row )
    }
    
    
    
    
  } else if( type == "irregular" ){
    
    longitude <- runif( 1000, 5, 95 )
    latitude  <- runif( 1000, 5, 95 )
    
    if( !is.null( cc$P ) ){
      P <- as.numeric( paste0( as.character( cc$P ), collapse="" ) )
      Q <- as.numeric( paste0( as.character( cc$Q ), collapse="" ) )
    }
    
    outpoints <- fields::cover.design( as.matrix( cbind(longitude,latitude) ) , npoints, P=P, Q=Q )
    #outpoints <- fields::cover.design( as.matrix( cbind(longitude,latitude) ) , npoints  )
    edge_seq  <- seq(0, 100, 100/nedge)
    out_des   <- dplyr::bind_rows( 
      as.data.frame( outpoints$design ), 
      data.frame( longitude =   0,      latitude = edge_seq    ),
      data.frame( longitude = 100,      latitude = edge_seq    ),
      data.frame( longitude = edge_seq, latitude =   0         ),
      data.frame( longitude = edge_seq, latitude = 100         )
    )
    
    del_out   <- deldir::deldir( out_des[["longitude"]], out_des[["latitude"]] , rw=c(0,100,0,100) )
    del_out02 <- deldir::tile.list(del_out)
    
    sim_shp00 <- list()
    for( ii in 1:length(del_out02) ){
      sim_shp00[[ii]] <- st_polygon(
        list( as.matrix( cbind(
          longitude = c( del_out02[[ii]]$x, del_out02[[ii]]$x[1] ) ,
          latitude  = c( del_out02[[ii]]$y, del_out02[[ii]]$y[1] )
        )))
      )
    }
    
  }
  
  ## Finalize and return the sf POLYGON object
  suid_n  <- length(sim_shp00)
  suid_c  <- str_pad( 1:suid_n , width=str_length(suid_n), side="left", pad="0"  )
  sim_shp <- st_sf( region   = paste0( suid, suid_c ) , 
                    geometry = st_sfc(sim_shp00)     )
  
  return( sim_shp )
  
}




#' Sets Incidence Rates for Simulated Regions
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Create the underlying rate of cases for each regular spatial unit. This function can be bypassed
#' if the user creates a variable names 'true_rate' in the set of regular spatial units. 
#' 
#' @param units_reg Set of regular spatial units
#' @param rate_base Vector of length 2 giving upper and lower bounds for uniform distribution. 
#'                  The base rate is randomly allocated between these two values. 
#'                  Default is `c(0.03, 0.07)`.
#' @param rate_spec A `data.frame` describing coordinates and change to base rate. See details.
#' 
#' @param seed If given, sets the seed for the RNG. 
#' @param ... Space for additional arguments (e.g., for `fields::cover.design`).
#'  
#' 
#' @details 
#' Describe the data.frame input.
#'  


gail_set_rate <- function( units_reg, rate_base=c(0.03,0.07), rate_spec=NULL ){
  
  if( is.null(rate_spec) ){
    stop("Must provide argument 'rate_spec'")
  }
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  #rate_base <- c(0.03, 0.07)
  nregion <- length( unique(units_reg[["region"]]) )
  
  ## Get the neighborhood matrix
  #nhbd <- spdep::poly2nb( units_reg )
  #Wmat <- spdep::nb2mat( nhbd, style="B" )
  
  # mvec <- as.matrix( rep(0,nregion) )
  n_rate_gen <- 2000
  loca_rate  <- data.frame( 
    longitude = runif(n_rate_gen,0,100) , 
    latitude  = runif(n_rate_gen,0,100)
  )
  
  xmat <- matrix( NA, nrow=n_rate_gen, ncol=nrow(rate_spec) )
  for( ii in 1:nrow(rate_spec) ){
    dx  <- abs( loca_rate[["longitude"]] - rate_spec[["mx"]][ii] )
    dy  <- abs( loca_rate[["latitude"]]  - rate_spec[["my"]][ii] )
    dxy <- ((dx^2) / rate_spec[["ax"]][ii]^2) + ((dy^2) / rate_spec[["ay"]][ii]^2) 
    xmat[,ii] <- (1 - 0.25*dxy)^2 * ( dxy^2 <= 4 )
    xmat[,ii] <- xmat[,ii]/max(xmat[,ii])
  }
  
  ## Initial Method: also averaged the uniform variable over the init points
  # mvec1 <- runif( nregion, min(rate_base), max(rate_base) )
  # mvec2 <- (xmat %*% rate_spec[["efc"]] ) 
  # mvec2 <- ifelse( mvec2==0, 1, mvec2 )
  # mvec  <- mvec1 * mvec2
  # 
  # loca_rate_st <- sf::st_as_sf( loca_rate , 
  #                               coords = c("longitude", "latitude") )
  # 
  # loca_rate_ruc <- st_contains( units_reg, loca_rate_st )
  # units_reg[["rate"]] <- NA
  # for( ii in 1:nregion ){
  #   units_reg[["rate"]][ii] <- mean( mvec[ loca_rate_ruc[[ii]] ]  )
  # }
  
  mvec <- (xmat %*% rate_spec[["efc"]] )
  mvec <- ifelse( mvec==0, 1, mvec )
  
  loca_rate_st <- sf::st_as_sf( loca_rate ,
                                coords = c("longitude", "latitude") )
  
  loca_rate_ruc <- st_contains( units_reg, loca_rate_st )
  mvec_base <- runif( nregion, min(rate_base), max(rate_base) )
  
  true_rate <- rep(NA,nregion)
  for( ii in 1:nregion ){
    true_rate[ii] <- mvec_base[ii] * mean( mvec[ loca_rate_ruc[[ii]] ]  )
  }
  
  return( true_rate )
  
}
















