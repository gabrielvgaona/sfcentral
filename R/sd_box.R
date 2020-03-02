#' @title Standard deviation box calculator in 2D or 3D
#'
#' @author Gabriel Gaona
#' @param points  matrix, datat.frame or tibble of points 2D or 3D
#' @param centre  Numeric. Coordinates 2D or 3D of central point. Default NULL,
#'   performs a calculation of mean_centre() from point localities
#' @param weights Numeric. Same length of number of points.
#' @param output Character. "coord" to return coordinates or sd box in 2D or 3D.
#'   "param" to return sd box attributes 2D or 3D. Default c("coord", "param"),
#'   to return both as list.
#' @return Depends on input, "coords" returns a data.frame of 2 or 3 columns and
#'   4 or 8 point coordinates. "param" returns a data.frame with centre
#'   coordinates, standard deviation in each axis, space(area for 2D, volume for
#'   3D) and number of dimensions in coordinates.
#' @importFrom Hmisc wtd.var
#' @importFrom stats sd
#' @export
#' @rdname sd_box
sd_box <- function(points, centre=NULL, weights=NULL,
                     output = c("coord", "param")) {

  if(!is.null(centre) & length(centre) < 2) {
    stop("Centre must be a c(x, y) or c(x, y, z) vector")
  }

  if(ncol(points) < 2) {
    # ERROR: TOO FEW COLUMNS IN DESTINATION COORDINATE MATRIX
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    stop("'points' must be of class matrix, data.frame or tibble with xy or xyz columns")
  }

  if(ncol(points) > 3) {
    warning("'points' more than 3 columns. Only first three columns are being used as xyz coordinates")
    points <- points[,1:3]
  }

  if(is.null(centre))  centre <- mean_centre(points, weights = weights)

  if(ncol(points) != length(centre)) {
    if(ncol(points) < length(centre)) {
      centre <- centre[1:2]
      warning("More 'centre' dimensions than 'points'. Extra dimensions droped...")
    } else {
      points <- points[,1:2]
      warning("More 'points' dimensions than 'centre'. Extra dimensions droped...")
    }
  }

  # STORE THE COUNT OF POINTS/CASES IN THE SOURCE DATASET
  n <- dim(points)[1]


  # INITIALIZE FUNCTION VARIABLE WITH PARAMETER VALUE
  dist <- distances(points, centre)

  # TEST WHETHER A SUFFICIENT NUMBER OF POINTS WERE SUPPLIED
  if(length(dist) < 3) stop("Not enough values to compute SDD")

  SD <- vector("numeric", ncol(points))
  names(SD) <- names(points)
  if(!is.null(weights)) {
    #PERFORM THE WEIGHTED STANDARD DEVIATION DISTANCE COMPUTATION (WEIGHTED SDD)
    SDD <- sqrt(sum((weights*dist^2)/((sum(weights)) - 2) ) )

    # COMPUTE AND ADD THE STANDARD DEVIATION OF THE X AND Y COORDINATES
    i = 1
    while(i <= ncol(points)){
      SD[i] <- sqrt(Hmisc::wtd.var(points[[i]], weights))
      i = i + 1
    }
  } else {
    # PERFORM THE STANDARD DEVIATION DISTANCE COMPUTATION (SDD)
    SDD <- sqrt(sum(dist^2/(length(dist) - 2) ) )

    # COMPUTE AND ADD THE STANDARD DEVIATION OF THE X AND Y COORDINATES
    i = 1
    while(i <= ncol(points)){
      SD[i] <- sd(points[[i]])
      i = i + 1
    }
  }

  # COMPUTE THE AREA OF THE SD BOX
  space <- prod(2*SD)

  # STORE THE COORDINATES OF EACH CORNER OF THE SD BOX IN SEPARATE OBJECTS


  dirs <- t(expand.grid(rep(list(c(1, -1)), length(SD))))
  box.points <- t(centre + dirs * SD)
  box.points <- as.data.frame(box.points)
  colnames(box.points) <- names(points)
  class(box.points) <- c("box_coords", class(box.points))

  param <- as.data.frame(rbind(c(Centre = centre,
          sd = SD,
          box_space = space,
          dim = ncol(points))))
  class(param) <- c("box_param", class(param))


  out <- list("coord" = box.points,
              "param" = param)[output]

  class(out) <- c("sd_box", class(out))

  if(length(out) == 1) out <- out[[1]]
  return(out)
}

