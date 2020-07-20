#' @title Standard deviation distance calculator (2D)
#' @author Gabriel Gaona
#' @param points  matrix, datat.frame or tibble of points 2D
#' @param x,y  Optional. numeric vectors for `x` and `y` coordinates of points.
#'   Default NULL perform extraction from point localities (`x = points[,1]` and
#'   `y <- points[,2]`)
#' @param centre  Numeric. Coordinates 2D of central point. Default NULL,
#'   performs a calculation of `mean_centre()` from point localities
#' @param weights Numeric. Same length of number of points.
#' @param output Character. `"coord"` to return coordinates for sd_distance
#'   points in 2D. "param" to return sd_distance attributes for 2D. Default
#'   c("coord", "param") to return both as list.
#' @return Depends on input, "coords" returns a data.frame of 2 and 360 point
#'   coordinates. "param" returns a data.frame with weigted indicator (`1 = TRUE`,
#'   or `0 = FALSE`), centre coordinates, standard deviation distance radius, and
#'   area of sd_distance.
#'
#' @export
#' @rdname sd_distance

sd_distance <- function (points = NULL,
                         x = NULL,
                         y = NULL,
                         centre = NULL,
                         weights = NULL,
                         output = c("coord", "param")) {

  if (all(is.null(points), is.null(x), is.null(y))) {
    stop("Points data.frame or x, y vectors must be provided")
  }

  # Control for few points coordiantes
  if(ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with x,y columns")
  }

  # Control for to many points coordiantes
  if(ncol(points) > 2) {
    warning("'points' more than 2 columns. Only first two columns are being used as xy coordinates")
    points <- points[,1:2]
  }
  if (is.null(points)) points <- cbind(x, y)
  if (!is.null(points)) x <- points[[1]]
  if (!is.null(points)) y <- points[[2]]

  weighted <- FALSE
  if (!is.null(weights)) weighted <- TRUE

  if (!is.null(centre)) {
    if (length(centre) == 2) {
      stop("Invalid center: must be a vector of length = 2")
    }
  }

  if (is.null(centre)) centre <- mean_centre(points = points)

  n <- dim(points)[1]
  dist <- sqrt((x - centre[1]) ^ 2 + (y - centre[2]) ^ 2)

  if (weighted) {
    SDD <- sqrt(sum((weights * dist ^ 2) / (sum(weights) - 2)))
  } else {
    SDD <- sqrt(sum(dist ^ 2 / (length(dist) - 2)))
  }

  sddarea <- pi * SDD^2
  sddperim <- pi * 2 * SDD
  B <- min(SDD, SDD)
  A <- max(SDD, SDD)
  d2 <- (A - B) * (A + B)
  phi <- 2 * pi * seq(0, 1, len = 360)
  sp <- sin(phi)
  cp <- cos(phi)
  r <- SDD * SDD / sqrt(B ^ 2 + d2 * sp ^ 2)
  xy <- r * cbind(cp, sp)
  al <- 0 * pi / 180
  ca <- cos(al)
  sa <- sin(al)

  coordsSDD <- t(t(xy %*% rbind(c(ca, sa), c(-sa, ca))) + centre)
  coordsSDD <- as.data.frame(coordsSDD)
  names(coordsSDD) <- colnames(points)[1:2]
  class(coordsSDD) <- c("sdd_coords", class(coordsSDD))

  result.sdd <- as.data.frame(
    rbind(
      c(weighted = weighted,
        centre = centre,
        radius = SDD,
        area = sddarea,
        perimeter = sddperim
      )
    )
  )
  class(result.sdd) <- c("sdd_param", class(result.sdd))

  out <- list("coord" = coordsSDD, "param" = result.sdd)[output]
  class(out) <- c("sd_distance", class(out))

  if (length(out) == 1)
    out <- out[[1]]
  return(out
  )
}

