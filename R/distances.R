#' @title Euclidean distance calculator in 2D or 3D
#'
#' @author Gabriel Gaona
#' @param destmat  matrix, datat.frame or tibble of target points 2D or 3D
#' @param centre Numeric. Coordinates 2D or 3D of central point.
#' @details Given a matriz of xy or xyz coordinates for a set of point locations
#'   compute the euclidean distances from a centre coordinate. This can be used
#'   in a planar coordinate system (e.g. projected coordinates in UTM)
#' @return Depends on input, returns a Numeric vector of distances to each
#'   location from the specified centre.
#' @importFrom splancs gridpts
#'
#' @export
distances <- function(destmat, centre = NULL) {
  # Control for few centre coordiantes
  if(!is.null(centre) & length(centre) < 2) {
    stop("Centre must be a c(x, y) or c(x, y, z) vector")
  }

  # Control for to many centre coordiantes
  if(!is.null(centre) & length(centre) > 3) {
    warning("'centre' more than 3 columns. Only first three columns are being used as xyz coordinates")
    centre <- unlist(centre)[1:3]
  }

  # Control for few destmat coordiantes
  if(ncol(destmat) < 2) {
    stop("'destmat' must be of class matrix, data.frame or tibble with x,y or x, y, z columns")
  }

  # Control for to many destmat coordiantes
  if(ncol(destmat) > 3) {
    warning("'desmat' more than 3 columns. Only first three columns are being used as xyz coordinates")
    destmat <- destmat[,1:3]
  }

  # Control for centre == null
  if(is.null(centre)) {
    warning("Centre  == NULL. Calculating mean_centre from input locations")
    centre <- mean_centre(destmat)
  }

  # Control for same dimensions in centre and destmat
  if(ncol(destmat) != length(centre)) {
    if(ncol(destmat) < length(centre)) {
      centre <- centre[1:2]
      warning("More 'centre' dimensions than 'destmat'. Extra dimensions droped...")
    } else {
      destmat <- destmat[,1:2]
      warning("More 'destmat' dimensions than 'centre'. Extra dimensions droped...")
    }
  }

  centre <- unlist(centre)
  # Compute euclidean distance
  sq = destmat
  i = 1
  while(i <= length(centre)){
    sq[[i]] <- (sq[[i]] - centre[i]) ^ 2
    i = i + 1
  }
  di <- sqrt(rowSums(sq))
  names(di) <- NULL
  # PROVIDE THE VECTOR OF DISTANCES, di, AS A RETURN PARAMETER TO THE CALLING FUNCTION
  return(di)
}
