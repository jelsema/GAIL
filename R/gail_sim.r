#' Generate Artificial Polygon Regions
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Create a set of spatial polygon regions.
#' For a more comprehensive description of how to use the simulation framework, 
#' see the vignette "GAIL Simulation."
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
#' @return
#' Returns an object of class `sf` with POLYGON geometry. This represents a set of areal 
#' 
#' 
#' @seealso
#' [gail_sim_rate], [gail_sim_pop], [gail_sim_index], [gail_sim_assign]
#' 
#' 
#' @examples
#' \dontrun{
#'  loca_reg <- gail_gen_regions( npoints=40, type="regular"  , nedge=10, suid="reg" )
#'  loca_irr <- gail_gen_regions( npoints=40, type="irregular", nedge=6 , seed=42 , P=-20, Q=20 )
#'   
#'  ggplot( loca_reg ) + 
#'    geom_sf( aes(fill=region) )
#'  ggplot( loca_irr ) +
#'    geom_sf( aes(fill=region) )
#' }
#' 
#' 
#' @importFrom deldir deldir tile.list
#' @importFrom dplyr bind_rows
#' @importFrom fields cover.design
#' @importFrom sf st_sf st_sfc st_polygon
#' @importFrom stats runif
#' @importFrom stringr str_length str_pad
#' 
#' @export
gail_sim_regions <- function( npoints, type="irregular", nedge=5, seed=NULL, suid=NULL, ... ){
  
  ## Reference to create SF polygons:
  ## https://stackoverflow.com/questions/44335246/polygons-from-coordinates
  
  if( is.numeric(seed) ){ set.seed( seed ) }
  if( is.null(suid) ){ suid <- "suid" }
  
  cc <- match.call(expand.dots = TRUE)
  
  if( type == "regular" ){
    
    make_row <- function(x){ 
      sf::st_polygon( list(
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
    
    longitude <- stats::runif( 1000, 5, 95 )
    latitude  <- stats::runif( 1000, 5, 95 )
    
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
      sim_shp00[[ii]] <- sf::st_polygon(
        list( as.matrix( cbind(
          longitude = c( del_out02[[ii]]$x, del_out02[[ii]]$x[1] ) ,
          latitude  = c( del_out02[[ii]]$y, del_out02[[ii]]$y[1] )
        )))
      )
    }
    
  }
  
  ## Finalize and return the sf POLYGON object
  suid_n  <- length(sim_shp00)
  suid_c  <- stringr::str_pad( 1:suid_n , width=str_length(suid_n), side="left", pad="0"  )
  sim_shp <- sf::st_sf( region   = paste0( suid, suid_c ) , 
                    geometry = sf::st_sfc(sim_shp00)     )
  
  return( sim_shp )
  
}




#' Sets Incidence Rates for Simulated Regions
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Create the underlying rate of cases for each regular spatial unit which can follow a
#' user-specified pattern. This can be used to generate the rate of cases, as well as the
#' rate of individuals being in the irregular spatial unit.
#' 
#' For creating simulated data, this function can be bypassed if the user creates variables 
#' named 'case_rate' and 'rural_rate' in the set of regular spatial units. 
#' 
#' @param units_reg Set of regular spatial units
#' @param rate_base Vector of length 2 giving upper and lower bounds for uniform distribution. 
#'                  The base rate is randomly allocated between these two values. 
#'                  Default is `c(0.03, 0.07)`.
#' @param rate_spec A `data.frame` describing coordinates and change to base rate. See details.
#' @param seed If given, sets the seed for the RNG. 
#' @param ... Space for additional arguments (e.g., for `fields::cover.design`).
#' 
#' @details 
#' Each row of `rate_spec` creates a 'hotspot' (or 'coldspot') in terms of the incidence 
#' rate of cases. That is: areas in the spatial domain which have an increased or decreased rate. This
#' should be a `data.frame` (or comprable object) with columns: `mx`, `my`, `ax`, `ay`, `efc`. The
#' hotspot is centered at the point (`mx`, `my`), while `ax` and `ay` control the size in the x and y 
#' directions, respectively, with larger values corresponding to larger range in that dimension. The `efc`
#' value is the effect size, which acts as a multiplier for the base rate of indidence.
#' 
#' The incidence rate is generated by first drawing a base rate for each spatial unit 
#' from a uniform distribution with bounds given by `rate_base`.
#' 
#' Then n=2000 points are drawn uniformly across the spatial domain (100x100 square). These points are
#' given a weight on the interval `[0, 1]`, which decreases from 1 at the center of the hotspot down to 0.
#' The weight is multiplied by the hotspot effect size (`efc`), and shifted so that each of the 2000
#' points have a value on the interval `[1, efc]`. The mean effect is taken for all individuals contained
#' within a spatial region, and that mean is used as a multiplier for the base rate of that region.
#' 
#' 
#' 
#' @return
#' Returns a vector of length `nrow(units_reg)` which contains 
#' 
#' 
#' @seealso
#' [gail_sim_regions], [gail_sim_pop], [gail_sim_index], [gail_sim_assign]
#' 
#' @examples 
#' \dontrun{
#'  ## Generate Regions
#'  loca_reg <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  
#'  ## Generate incidence rate
#'  rate_spec <- data.frame(
#'    mx  = c(25, 60), 
#'    my  = c(25, 80), 
#'    ax  = c(10, 40), 
#'    ay  = c(25, 20),
#'    efc = c( 0.15 , 4.0 )
#'  )
#'  loca_reg[["case_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'  ggplot( loca_reg ) +
#'    geom_sf( aes(fill=case_rate) )
#'  
#'  ## Generate rate of being in irregular locations
#'  irr_spec <- data.frame(
#'    mx  = c(85, 20, 25, 60), 
#'    my  = c(15, 80, 25, 80), 
#'    ax  = c(20, 10, 10, 40), 
#'    ay  = c(20, 20, 25, 20),
#'    efc = c(4.0, 4.0, 0.15 , 0.15 )
#'  )
#'   
#'  loca_reg[["rural_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'     
#'  ggplot( loca_reg ) +
#'    geom_sf( aes(fill=rural_rate) )
#'     
#' }
#' 
#' 
#' @importFrom sf st_as_sf st_contains
#' @importFrom stats runif
#' 
#' 
#' @export
gail_sim_rate <- function( units_reg, rate_base=c(0.03,0.07), rate_spec=NULL, seed=NULL ){
  
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
    longitude = stats::runif(n_rate_gen,0,100) , 
    latitude  = stats::runif(n_rate_gen,0,100)
  )
  
  xmat <- matrix( NA, nrow=n_rate_gen, ncol=nrow(rate_spec) )
  for( ii in 1:nrow(rate_spec) ){
    dx  <- abs( loca_rate[["longitude"]] - rate_spec[["mx"]][ii] )
    dy  <- abs( loca_rate[["latitude"]]  - rate_spec[["my"]][ii] )
    dxy <- ((dx^2) / rate_spec[["ax"]][ii]^2) + ((dy^2) / rate_spec[["ay"]][ii]^2) 
    xmat[,ii] <- (1 - 0.25*dxy)^2 * ( dxy^2 <= 4 )
  }
  
  loca_rate_st <- sf::st_as_sf( loca_rate ,
                                coords = c("longitude", "latitude") )
  
  loca_rate_ruc <- sf::st_contains( units_reg, loca_rate_st )
  mvec_base <- stats::runif( nregion, min(rate_base), max(rate_base) )
  
  case_rate <- rep(NA,nregion)
  for( ii in 1:nregion ){
    idx_region    <- loca_rate_ruc[[ii]]
    rate_change   <- sum( colMeans( xmat[ idx_region , ] ) * rate_spec[["efc"]] )
    case_rate[ii] <- ilogit( logit(mvec_base[ii]) + rate_change )
  }
  
  return( case_rate )
  
}


#' Simulate Population for GAIL
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Create a simulated population distributed across the spatial domain. This population can be
#' sampled to create cases, and to create individuals reporting the irregular spatial unit using
#' [gail_sim_assign].
#' 
#' @param units_reg Set of regular spatial units (cases get allocated to this set).
#' @param units_irr Set of irregular spatial units (cases get allocated from this set).
#' @param method Method of simulating population: 'uniform' or 'beta'. See details.
#' @param npop Size of population to simulate.
#' @param beta_setup If `method='beta'`, a `data.frame` describing how to distribute the locations. See details.
#' 
#' @param seed If given, sets the seed for the RNG. 
#' @param ... Space for additional arguments (e.g., for `fields::cover.design`).
#'  
#' 
#' @details 
#' For `method='uniform'` points are simulated uniformly across the 100x100 spatial domain.
#' For `method='irregular'` then three additional parameters are necessary: A list of centers, 
#' and a list of values for alpha and beta. For each center, points are simulated from a beta 
#' distribution. 
#' 
#' The argument `beta_setup` sets the parameters for the beta method of distributing the population.
#' This can be used to generate a population which is clustered in certain areas. Thise should be a
#' `data.frame` (or comprable object) with columns: `nn`, `mx`, `my`, `sx`, `sy`. These columns
#' represent the number of individuals in the cluster (`nn`). The cluster is centered at the point 
#' (`mx`, `my`), while `sx` and `sy` are the standard deviation in the x and y dimensions, respectively.
#' In one dimension, the cluster will be drawn from a beta distribution (scaled to `[0, 100]`) with mean
#' of `mx` and standard deviation of `sx`. The \eqn{\alpha} and \eqn{\beta} parameters of the beta distribution
#' are derived from the mean and standard deviation.
#' 
#' 
#' @seealso
#' [gail_sim_regions], [gail_sim_rate], [gail_sim_index], [gail_sim_assign]
#' 
#' 
#' @examples 
#' \dontrun{
#'  ## Generate Regions
#'  loca_reg <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  loca_irr <- gail_gen_regions( npoints=40, type="irregular", nedge=6 , seed=42 , P=-20, Q=20 )
#'  
#'  ## Generate incidence rate
#'  rate_spec <- data.frame(
#'    mx  = c(25, 60), 
#'    my  = c(25, 80), 
#'    ax  = c(10, 40), 
#'    ay  = c(25, 20),
#'    efc = c( 0.15 , 4.0 )
#'  )
#'  loca_reg[["case_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'  
#'  ## Generate rate of being in irregular locations
#'  irr_spec <- data.frame(
#'    mx  = c(85, 20, 25, 60), 
#'    my  = c(15, 80, 25, 80), 
#'    ax  = c(20, 10, 10, 40), 
#'    ay  = c(20, 20, 25, 20),
#'    efc = c(4.0, 4.0, 0.15 , 0.15 )
#'  )
#'   
#'  loca_reg[["rural_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'     
#'  ## Generate population
#'  beta_setup <- data.frame(
#'    nn=c(5000, 1000, 500),
#'    mx=c(50, 25, 60), 
#'    my=c(50, 25, 80), 
#'    sx=c(30, 10, 10), 
#'    sy=c(30, 10, 5)
#'  )
#'  loca_pop <- gail_sim_pop( loca_irr, loca_irr, method="beta", 
#'                            beta_setup=beta_setup, seed=42 )
#'  
#' }
#' 
#' 
#' @importFrom dplyr mutate
#' @importFrom purrr pmap
#' @importFrom stats rbeta runif
#' @importFrom sf st_as_sf
#' 
#' @export
gail_sim_pop <- function( units_reg, units_irr, method="uniform", npop=100000, beta_setup=NULL, seed=NULL ){
  
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  ## ##################################################
  ## May want to add functionality to allow for non-square regions (e.g., a state)
  ## This would need:
  ##  - Extracting the bounds for the region of interest
  ##  - Simulating within these bounds
  ##  - Checking that the simulated points are all within the region 
  ##   
  ##  - Alternatively: Generate the sample sizes from a poisson distribution
  ##  - Simulate that number of points within each region individually.
  ## 
  
  
  ## ##################################################
  ## Uniform distributed population 
  ##
  
  if( method=="uniform" ){
    
    locas01 <- data.frame( 
      longitude = stats::runif(npop,0,100) , 
      latitude  = stats::runif(npop,0,100)
    )
    
    locas02 <- sf::st_as_sf( locas01 , coords = c("longitude", "latitude") )
    
  }
  
  
  ## ##################################################
  ## Beta distributed population 
  ## 
  
  if( method=="beta" ){
    
    if( is.null(beta_setup) ){
      beta_setup <- data.frame(
        nn=npop, mx=50, my=50, sx=30, sy=30
      )
    }
    
    beta_setup <- dplyr::mutate( beta_setup,
                                 rx = (100/mx - 1),
                                 ax = (rx / ( ((sx/100)^2) * ( (rx+1)^3 ) )) - (1/(rx+1)),
                                 bx = rx*ax,
                                 ry = (100/my - 1),
                                 ay = (ry / ( ((sy/100)^2) * ( (ry+1)^3 ) )) - (1/(ry+1)),
                                 by = ry*ay
    )
    
    beta_parms_x <- list(
      n      = beta_setup[["nn"]],
      shape1 = beta_setup[["ax"]],
      shape2 = beta_setup[["bx"]]
    )
    beta_parms_y <- list(
      n      = beta_setup[["nn"]],
      shape1 = beta_setup[["ay"]],
      shape2 = beta_setup[["by"]]
    )
    
    longitude <- unlist( purrr::pmap( beta_parms_x , .f=stats::rbeta ) )
    latitude  <- unlist( purrr::pmap( beta_parms_y , .f=stats::rbeta ) )
    locas02 <- sf::st_as_sf( data.frame(longitude, latitude)*100 , 
                             coords = c("longitude", "latitude") )
    
  }
  
  ## ##################################################
  ## Add region ID to the population  
  ## ruc = "regular unit contains"
  ## iuc = "irregular unit contains"
  
  ruc <- sf::st_contains( units_reg, locas02 )
  units_reg[["pop"]] <- sapply( ruc, FUN=length )
  
  locas02[["region"]] <- factor( NA, levels=unique( units_reg[["region"]]) )
  for( ii in 1:nrow(units_reg) ){
    locas02[["region"]][ ruc[[ii]]  ] <- units_reg[["region"]][ii]
  }
  
  
  iuc <- sf::st_contains( units_irr, locas02 )
  units_irr[["pop"]] <- sapply( iuc, FUN=length )
  
  locas02[["iregion"]] <- factor( NA, levels=unique( units_irr[["region"]]) )
  for( ii in 1:nrow(units_irr) ){
    locas02[["iregion"]][ iuc[[ii]]  ] <- units_irr[["region"]][ii]
  }
  
  id_n  <- nrow(locas02)
  id_c  <- str_pad( 1:id_n , width=str_length(id_n), side="left", pad="0"  )
  locas02[["pop_id"]] <- paste0( "lpid_", id_c )
  
  ## Make return object
  
  return( locas02 )
  
}



#' Generate a Simulated Similarity Index
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Creates a spatially dependent similarity index across the regions. One of the method of 
#' geo-allocation with [gail] is to use such an index.
#' 
#' @param units_sp Set of spatial units
#' @param tau Nugget variance parameter, see details.
#' @param phi Spatial dependence parameter, see details.
#' @param seed If given, sets the seed for the RNG. 
#'  
#' 
#' @details 
#' Generates a spatially-dependent similarity index *y* by first simulating 
#' from a multivariate normal distribution:
#' 
#' \eqn{ Y \sim \mbox{MVN}( 0, \tau^{2}( I - \phi H) ) }
#' 
#' where H is the neighborhood matrix, and I is the identity matrix. Then
#' the index is converted to a `[0, 15]` uniform random variable.
#'  
#'  
#' @seealso
#' [gail_sim_regions], [gail_sim_rate], [gail_sim_pop], [gail_sim_assign]
#' 
#' @examples 
#' \dontrun{
#'  
#'  ## Generate Regions
#'  loca_reg <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  loca_irr <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  
#'  ## Generate incidence rate
#'  rate_spec <- data.frame(
#'    mx  = c(25, 60), 
#'    my  = c(25, 80), 
#'    ax  = c(10, 40), 
#'    ay  = c(25, 20),
#'    efc = c( 0.15 , 4.0 )
#'  )
#'  loca_reg[["case_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'  
#'  ## Generate rate of being in irregular locations
#'  irr_spec <- data.frame(
#'    mx  = c(85, 20, 25, 60), 
#'    my  = c(15, 80, 25, 80), 
#'    ax  = c(20, 10, 10, 40), 
#'    ay  = c(20, 20, 25, 20),
#'    efc = c(4.0, 4.0, 0.15 , 0.15 )
#'  )
#'   
#'  loca_reg[["rural_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'     
#'  ## Generate population
#'  beta_setup <- data.frame(
#'    nn=c(5000, 1000, 500),
#'    mx=c(50, 25, 60), 
#'    my=c(50, 25, 80), 
#'    sx=c(30, 10, 10), 
#'    sy=c(30, 10, 5)
#'  )
#'  loca_pop <- gail_sim_pop( loca_reg, loca_irr, method="beta", 
#'                            beta_setup=beta_setup, seed=42 )
#'  
#'  loca_reg[["index"]] <- gail_sim_index( loca_reg, tau=1.5, phi=0.05 )
#'  loca_irr[["index"]] <- gail_sim_index( loca_irr, tau=1.5, phi=0.05 )
#'  
#' }
#' 
#' 
#' @importFrom MASS mvrnorm
#' @importFrom spdep poly2nb
#' @importFrom spdep nb2mat
#' @importFrom stats pnorm uniroot
#'  
#' @export
gail_sim_index <- function( units_sp, tau, phi, seed=NULL ){
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  #rate_base <- c(0.03, 0.07)
  nregion <- length( unique(units_sp[["region"]]) )
  
  ## Get the neighborhood matrix
  nhbd <- spdep::poly2nb( units_sp )
  Wmat <- spdep::nb2mat( nhbd, style="B" )
  Imat <- diag(nregion)
  
  ## May need to add a check here
  ## Revert to manually compute generalized Sigma^{-1/2} if needed?
  
  ## Create small helper function
  min_eig <- function( x ){
    mm <- (Imat - x*Wmat)
    min( eigen(mm)$values )
  }
  
  max_phi <- stats::uniroot( min_eig, interval=c(0,20) )$root
  
  if( phi >= max_phi ){
    stop( paste0("Argument 'phi' is too large, must be < ", max_phi) )
  }
  
  ## May want to add in parameterization to control mean?
  ## Is it needed, if there is spatial correlation?
  # YES -- to be able to align regular and irregular units.
  mvec <- as.matrix( rep(0,nregion) )
  smat <- (tau^2) * solve( Imat - phi * Wmat )
  # smat <- tau * MASS::ginv( Imat - pp*Wmat )
  
  idx_norm <- MASS::mvrnorm( 1, mvec, smat )
  idx_unif <- round( stats::pnorm( idx_norm ) * 15 )
  
  return( idx_unif )
  
}






#' Allocate Simulated Population
#'
#' @description
#' Function for the simulation framework in GAIL.
#' Determined the spatial unit (both regular and irregular) into which a population individual
#' falls, and randomly samples whether an individual is case vs non-case and whether an individual
#' is in regular or non-regular spatial unit.
#' 
#' @param units_reg Set of regular spatial units
#' @param units_irr Set of irregular spatial units
#' @param loca_pop True locations of the population
#' @param seed If given, sets the seed for the RNG. 
#'  
#' 
#' @details 
#' The argument `units_reg` must have columns `case_rate` and `rural_rate`. These will
#' ordinarily be generated by other simulation functions in the process of simulating
#' the data.
#' 
#' This function assigns case/non and regular/irregular by:
#' 1. Assigning entire population to the regular spatial units
#' 2. For each regular spatial unit, independently sampling based on case rate (to generate cases)
#'    and rural rate (to generate individuals reporting the irregular spatial unit).
#' 
#' @seealso
#' [gail_sim_regions], [gail_sim_rate], [gail_sim_pop], [gail_sim_index]
#' 
#' @examples 
#' \dontrun{
#'  
#'  ## Generate Regions
#'  loca_reg <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  loca_irr <- gail_gen_regions( npoints=40, type="regular", nedge=10, suid="reg" )
#'  
#'  ## Generate incidence rate
#'  rate_spec <- data.frame(
#'    mx  = c(25, 60), 
#'    my  = c(25, 80), 
#'    ax  = c(10, 40), 
#'    ay  = c(25, 20),
#'    efc = c( 0.15 , 4.0 )
#'  )
#'  loca_reg[["case_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'  
#'  ## Generate rate of being in irregular locations
#'  irr_spec <- data.frame(
#'    mx  = c(85, 20, 25, 60), 
#'    my  = c(15, 80, 25, 80), 
#'    ax  = c(20, 10, 10, 40), 
#'    ay  = c(20, 20, 25, 20),
#'    efc = c(4.0, 4.0, 0.15 , 0.15 )
#'  )
#'   
#'  loca_reg[["rural_rate"]] <- gail_sim_rate( loca_reg,  rate_base=c(0.03,0.07), 
#'                                             rate_spec=rate_spec, seed=42 )
#'     
#'  ## Generate population
#'  beta_setup <- data.frame(
#'    nn=c(5000, 1000, 500),
#'    mx=c(50, 25, 60), 
#'    my=c(50, 25, 80), 
#'    sx=c(30, 10, 10), 
#'    sy=c(30, 10, 5)
#'  )
#'  loca_pop <- gail_sim_pop( loca_reg, loca_irr, method="beta", 
#'                            beta_setup=beta_setup, seed=42 )
#'  
#'  ## Simulate index
#'  loca_reg[["index"]] <- gail_sim_index( loca_reg, tau=1.5, phi=0.05 )
#'  loca_irr[["index"]] <- gail_sim_index( loca_irr, tau=1.5, phi=0.05 )
#'  
#'  
#'  ## Assign cases and spatial unit
#'  gsa01 <- gail_sim_assign( loca_reg, loca_irr, loca_pop, seed=42 )
#'  
#'  
#' }
#' 
#' @importFrom dplyr group_by filter left_join mutate rename sample_n select summarize
#' @importFrom purrr map2
#' @importFrom magrittr %>%
#' @importFrom sf st_coordinates st_geometry
#' @importFrom tidyr nest unnest
#'  
#' @export
gail_sim_assign <- function( units_reg, units_irr, loca_pop, seed=NULL ){
  
  if( is.numeric(seed) ){
    set.seed( seed )
  }
  
  if( !("case_rate" %in% colnames( units_reg )) ){
    stop( "Argument 'units_reg' must have column 'case_rate'")
  }
  if( !("rural_rate" %in% colnames( units_reg )) ){
    stop( "Argument 'units_reg' must have column 'rural_rate'")
  }
  
  ## Step 1: Assign population as case/non case and reg/irr
  
  ## Maybe replace with a loop 
  ##   - Would enable foreach() 
  ##   - Also might be able to switch to base R
  ## Could replace independent sampling with (possibility):
  ##   - Generate bivariate normal
  ##   - Compute Phi(Z)
  ##   - Set Z[1] <= X as case and Z[2] <= Y as rural
  ##   - Would need fairly WEAK correlation, just to get some shape
  
  pop_table <- table( loca_pop[["region"]] )
  
  full_pop <- data.frame(  region = names(pop_table), 
                           npop   = as.numeric(pop_table) )
  
  units_reg <- units_reg %>%
    dplyr::left_join( full_pop, by="region" ) %>%
    dplyr::mutate( 
      ncases = round( npop * case_rate  ),
      nrural = round( npop * rural_rate )
    )
  
  loca_pop_nest <- loca_pop %>%
    dplyr::group_by( region ) %>%
    tidyr::nest( ) %>%
    dplyr::left_join( units_reg , by="region" )
  
  loca_pop_cases <- loca_pop_nest %>%
    dplyr::mutate( case_status = purrr::map2( data, ncases, dplyr::sample_n ) ) %>%
    dplyr::select( region, case_status ) %>%
    tidyr::unnest() %>%
    dplyr::select( -geometry, -region, -iregion ) %>%
    dplyr::mutate( case = "case" )
  
  loca_pop_types <- loca_pop_nest %>%
    dplyr::mutate( rural_status = purrr::map2( data, nrural, dplyr::sample_n ) ) %>%
    dplyr::select( region, rural_status ) %>%
    tidyr::unnest() %>%
    dplyr::select( -geometry, -region, -iregion) %>%
    dplyr::mutate( type = "irregular" )
  
  loca_pop_assigned <- loca_pop %>%
    dplyr::left_join( loca_pop_cases , by="pop_id" ) %>%
    dplyr::left_join( loca_pop_types , by="pop_id" ) %>%
    dplyr::rename( rregion = region ) %>%
    dplyr::mutate(
      rregion = as.character(rregion),
      iregion = as.character(iregion),
      case    = ifelse( is.na(case), "noncase", case ),
      type    = ifelse( is.na(type), "regular", type ),
      region  = ifelse( type=="regular", rregion, iregion )
    )
  
  loca_coords <- dplyr::as_tibble( sf::st_coordinates( loca_pop_assigned ) )
  colnames(loca_coords) <- c("longitude", "latitude")
  sf::st_geometry( loca_pop_assigned ) <- NULL
  loca_pop_assigned <- data.frame( loca_pop_assigned, loca_coords )
  
  ## Step 2: Compute summaries for units_reg and units_irr
  
  ## Geometry == POLYGON
  
  if( "index" %in% colnames(units_reg) ){
    sp_units <- units_reg %>%
      dplyr::select( region, index, npop ) 
    
    prep_cases01 <- units_reg %>% 
      dplyr::select( region, index, geometry  ) %>%
      rbind( units_irr ) 
  } else{
    sp_units <- units_reg %>%
      dplyr::select( region, npop ) 
    
    prep_cases01 <- units_reg %>% 
      dplyr::select( region, geometry  ) %>%
      rbind( units_irr )
  }
  
  prep_cases02 <- loca_pop_assigned %>%
    dplyr::filter( case=="case" ) %>%
    dplyr::select( region, case  ) %>%
    dplyr::group_by( region ) %>%
    dplyr::summarize( ncases = sum( case=="case" ) ) %>%
    dplyr::rename( region = region )
  
  cases <- prep_cases01 %>%
    dplyr::left_join( prep_cases02 , by="region" ) %>%
    dplyr::filter( !is.na(ncases) )
  
  
  ## Make output object
  # units_reg_truth   = Ground truth, how many cases there actually are in a spatial unit
  # sp_units          = set of regular spatial units with population and index
  # cases             = set of cases with spatial unit
  # loca_pop_assigned = population with case and rurality status, for arbitrary use
  
  return_obj <- list(
    units_reg_truth = units_reg,
    sp_units        = sp_units,
    cases           = cases,
    loca_pop        = loca_pop_assigned
  )
  
  return( return_obj )
  
}








