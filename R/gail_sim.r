




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
#' @import dplyr bind_rows
#' @import sf
#' @import fields cover.design
#' @import deldir
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'   tmp_out_01 <- gail_gen_regions( npoints=40, type="regular"  , nedge=10, suid="reg" )
#'   tmp_out_02 <- gail_gen_regions( npoints=40, type="irregular", nedge=6 , seed=42 , P=-20, Q=20 )
#'   
#'   ggplot( tmp_out_01 ) + 
#'     geom_sf( aes(fill=region) )
#'   ggplot( tmp_out_02 ) +
#'     geom_sf( aes(fill=region) )
#' }
#' 



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
    
    xdim <- runif( 1000, 5, 95 )
    ydim <- runif( 1000, 5, 95 )
    
    if( !is.null( cc$P ) ){
      P <- as.numeric( paste0( as.character( tmp_out_02$P ), collapse="" ) )
      Q <- as.numeric( paste0( as.character( tmp_out_02$Q ), collapse="" ) )
    }
    
    outpoints <- fields::cover.design( as.matrix( cbind(xdim,ydim) ) , npoints, P=P, Q=Q )
    #outpoints <- fields::cover.design( as.matrix( cbind(xdim,ydim) ) , npoints  )
    edge_seq  <- seq(0, 100, 100/nedge)
    out_des   <- dplyr::bind_rows( 
      as.tibble( outpoints$design ), 
      data.frame( xdim =   0,      ydim = edge_seq    ),
      data.frame( xdim = 100,      ydim = edge_seq    ),
      data.frame( xdim = edge_seq, ydim =   0         ),
      data.frame( xdim = edge_seq, ydim = 100         )
    )
    
    del_out   <- deldir::deldir( out_des[["xdim"]], out_des[["ydim"]] , rw=c(0,100,0,100) )
    del_out02 <- deldir::tile.list(del_out)
    
    sim_shp00 <- list()
    for( ii in 1:length(del_out02) ){
      sim_shp00[[ii]] <- st_polygon(
        list( as.matrix( cbind(
          xdim = c( del_out02[[ii]]$x, del_out02[[ii]]$x[1] ) ,
          ydim = c( del_out02[[ii]]$y, del_out02[[ii]]$y[1] )
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




















