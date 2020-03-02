#' @title Standard deviation ellipse calculator (2D)
#'
#' @author Gabriel Gaona
#' @param points  matrix, datat.frame or tibble of points 2D
#' @param x,y  Optional. numeric vectors for `x` and `y` coordinates of points.
#'   Default NULL perform extraction from point localities (`x = points[,1]` and
#'   `y <- points[,2]`)
#' @param centre  Numeric. Coordinates 2D of central point. Default NULL,
#'   performs a calculation of `mean_centre()` from point localities
#' @param weights Numeric. Same length of number of points.
#' @param output Character. `"coord"` to return coordinates for sd_ellipse
#'   points in 2D. "param" to return sd_ellipse attributes for 2D. Default
#'   c("coord", "param") to return both as list.
#' @return Depends on input, "coords" returns a data.frame of 2 and 360 point
#' coordinates. "param" returns a data.frame with weigted indicator (`1 = TRUE`,
#' or `0 = FALSE`), centre coordinates, standard deviation distance radius, and
#' area of sd_distance.
#' @importFrom purrr map_dbl
#' @export
#' @rdname sd_ellipse
sd_ellipse <- function (points = NULL, x = NULL, y = NULL,
                       centre = NULL, weights = NULL,
                       output = c("coord", "param")) {

  if (all(is.null(points), is.null(x), is.null(y))) {
    stop("Points data frame or x, y vectors must be provided")
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
  if (dim(points)[2] != 2) {
    stop("Points must have two columns 'x' and 'y'")
  }
  if (!is.null(points)) x <- points[[1]]
  if (!is.null(points)) y <- points[[2]]

  n <- dim(points)[1]
  if (!is.null(centre)) {
    if (length(centre) != 2) {
      stop("Invalid center: must be a vector of length = 2")
    }
  }

  if (is.null(centre)) centre <- mean_centre(points = points,
                                                   weights = weights)


  tmpA <- points
  tmpA2 <- tmpA ^ 2
  tmpAc <- t(t(tmpA) - centre)
  tmpAc2 <- tmpAc ^ 2
  tmpAcp <- purrr::map_dbl(seq_len(nrow(tmpA)), ~ prod(tmpAc[.x, ]))

  if (!is.null(weights)) {
    top1 <- sum(rowSums(t(t(weights * tmpAc2) * c(1,-1))))

    top2 <- sqrt(sum(top1) ^ 2 + 4 * sum(weights * tmpAcp) ^ 2)
    bottom <- (2 * sum(weights * tmpAcp))
    tantheta <- (top1 + top2) / bottom
  } else {
    top1 <- sum(rowSums(t(t(tmpAc2) * c(1,-1))))
    top2 <- sqrt(top1 ^ 2 +  4 * sum(tmpAcp) ^ 2)
    bottom <- (2 * sum(tmpAcp))
    tantheta <- (top1 + top2) / bottom
  }

  if (tantheta < 0) {
    theta <- 180 + atan_d(tantheta)
  } else {
    theta <- atan_d(tantheta)
  }
  sintheta <- sin_d(theta)
  costheta <- cos_d(theta)
  sin2theta <- sintheta ^ 2
  cos2theta <- costheta ^ 2
  sinthetacostheta <- sintheta * costheta
  if (!is.null(weights)) {
    sigmax <-
      sqrt(2) * sqrt((sum(colSums(weights * tmpAc2) * c(cos2theta, sin2theta)) -
        2 * sum(weights * tmpAcp) * sinthetacostheta) / sum(weights))

    sigmay <-
      sqrt(2) * sqrt((sum(colSums(weights * tmpAc2) * c(sin2theta, cos2theta)) +
        2 * (sum(weights * tmpAcp)) *  sinthetacostheta) / sum(weights))
  } else {
    sigmax <-
      sqrt(2) * sqrt((sum(colSums(tmpAc2) * c(cos2theta, sin2theta)) -
        2 * sum(tmpAcp) * sinthetacostheta) / (n - 2))
    sigmay <-
      sqrt(2) * sqrt((sum(colSums(tmpAc2) * c(sin2theta, cos2theta)) +
        2 * (sum(tmpAcp)) *  sinthetacostheta) /
        (n - 2))
  }

  if (sigmax > sigmay) {
    Major <- "SigmaX"
    Minor <- "SigmaY"
  }  else {
    Major <- "SigmaY"
    Minor <- "SigmaX"
  }
  lengthsigmax <- 2 * sigmax
  lengthsigmay <- 2 * sigmay
  areaSDE <- pi * sigmax * sigmay
  eccentricity <- sqrt(1 - ((min(sigmax, sigmay) ^ 2) / (max(sigmax, sigmay) ^ 2)))
  B <- min(sigmax, sigmay)
  A <- max(sigmax, sigmay)
  d2 <- (A - B) * (A + B)
  phi <- 2 * pi * seq(0, 1, len = 360)
  sp <- sin(phi)
  cp <- cos(phi)
  r <- sigmax * sigmay / sqrt(B ^ 2 + d2 * sp ^ 2)
  xy <- r * cbind(cp, sp)
  al <- (90 - theta) * pi / 180
  ca <- cos(al)
  sa <- sin(al)
  coordsSDE <- xy %*% rbind(c(ca, sa), c(-sa, ca)) +
    matrix(centre, nrow = 360, ncol = 2, byrow = TRUE)
  coordsSDE <- as.data.frame(coordsSDE)
  names(coordsSDE) <- names(points)[1:2]
  class(coordsSDE) <- c("sde_coords", class(coordsSDE))

  if (sigmax < sigmay) {
    Theta.Corr <- theta
  }  else {
    Theta.Corr <- theta + 90
  }

  result.sde <- as.data.frame(rbind(
    c(weighted = !is.null(weights), centre = centre, Sigma.x = sigmax,
      Sigma.y = sigmay, Theta = theta, Eccentricity = eccentricity,
      area.sde = areaSDE, tan_theta = tantheta, SinTheta = sintheta,
      CosTheta = costheta, SinThetaCosTheta = sinthetacostheta,
      Sin2Theta = sin2theta, Cos2Theta = cos2theta, ThetaCorr = Theta.Corr)
    ))
  class(result.sde) <- c("sde_param", class(result.sde))

  out <- list("coord" = coordsSDE,
              "param" = result.sde)[output]
  class(out) <- c("calc_sde", class(out))

  if(length(out) == 1) out <- out[[1]]
  return(out)
}


