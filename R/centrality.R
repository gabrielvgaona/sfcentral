#' @title Centrality calculator in 2D or 3D
#' @author Gabriel Gaona
#' @note Mean, Median and Central features calculators are modified versions
#'   of `aspace::*()` from Ron Buliung & Randy Bui (2012)
#' @param points,points2  matrix, datat.frame or tibble of points 2D or 3D
#' @param weights Numeric. Same length of number of points.
#' @param dist Atomic numeric, Default 100. Hold distance value between i
#'   and ith iterations.
#' @details Central minimum distance can only performed in 2D. Other centrality
#'   calculations could be used in 2D or 3D points
#' @return Depends on input, returns a Numeric vector of 2 or 3 values
#' @importFrom splancs gridpts
#' @importFrom grDevices chull
#' @importFrom stats dist median sd

#' @export
#' @rdname centrality
median_centre <- function(points) {
  # Control for few columns in points
  if (ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with 2 or 3 numeric columns")
  }

  # Control for to many columns in points
  if (ncol(points) > 3) {
    warning(
      "'points' more than 3 columns.\nOnly first three columns are being used as xyz coordinates"
    )
    points <- points[, 1:3]
  }

  # Compute median centre
  mc <- vector("numeric", ncol(points))
  i = 1
  while (i <= ncol(points)) {
    mc[i] <- stats::median(points[[i]])
    i = i + 1
  }
  names(mc) <- colnames(points)
  mc
}

#' @export
#' @rdname centrality
mean_centre <- function(points, weights = NULL) {
  # Control for few columns in points
  if (ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with 2 or 3 numeric columns")
  }

  #  Control for to many columns in points
  if (ncol(points) > 3) {
    warning(
      "'points' more than 3 columns.\nOnly first three columns are being used as xyz coordinates"
    )
    points <- points[, 1:3]
  }

  # Mean centre calculation
  n <- dim(points)[1]
  weights <- if (is.null(weights)) 1/n else weights / sum(weights)
  i = 1
  while (i <= ncol(points)) {
      points[,i] <- points[,i] * weights
      i = i + 1
    }

  centre <- colSums(points)
  centre
}


#' @export
#' @rdname centrality
central_feature <- function(points) {
  # Control for few columns in points
  if (ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with 2 or 3 numeric columns")
  }

  #  Control for to many columns in points
  if (ncol(points) > 3) {
    warning(
      "'points' more than 3 columns.\nOnly first three columns are being used as xyz coordinates"
    )
    points <- points[, 1:3]
  }

  # Central feature calculation
  d <- stats::dist(points, diag = TRUE, upper = TRUE)
  ds <- rowSums(as.matrix(d))
  unlist(points[which.min(ds), ])
}

#' @importFrom splancs gridpts
#' @export
#' @rdname centrality
centre_min_dis <- function(points, dist = 100) {
  # Control for few columns in points
  if (ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with 2 numeric columns")
  }

  #  Control for to many columns in points
  if (ncol(points) > 3) {
    warning(
      "'points' more than 2 columns.\nOnly first two columns are being used as 2D coordinates"
    )
    points <- points[, 1:2]
  }

  # Initializing
  x <- c() #HOLD X-COORD OF CMD FOR EACH ITERATION
  y <- c() #HOLD Y-COORD OF CMD FOR EACH ITERATION
  d <- c() #HOLD DISTANCE VALUE BETWEEN I AND ITH ITERATIONS
  n <- c() #HOLD ITERATION NUMBER
  cells <-
    c() #HOLD NUMBER OF CELLS IN EACH GRID CREATED FOR EACH ITERATION

  # INITIALIZE OBJECTS, COUNTERS, GRID SIZE
  i <- 1
  xy <- points[1, ]
  xy[] <- 0
  d[i] <- dist
  n <- 0
  cells[i] <- 0
  dx <- 1 #INITIALIZE GRID SPACING IN X
  dy <-
    1 #INITIALIZE GRID SPACING IN Y, LARGER NUMBER = MORE CELLS, HERE IT IS SET TO 1 CELLS IN X,Y

  # GENERATE MCP
  hpts <- chull(points)

  MCP <- cbind(points[hpts, 1], points[hpts, 2])

  while (d[i] >= dist) {
    grid <- gridpts(MCP, dx, dy)

    sumdist.CMD <- vector("numeric", nrow(grid))
    j = 1
    while (j <= nrow(grid)) {
      sumdist.CMD[j] <- sum(distances(points, centre = grid[j, ]))
      j = j + 1
    }
    CMD <- grid[which.min(sumdist.CMD), ]
    #Dump CMD for each interation
    xy[i + 1, ] <- CMD
    #distance between current and previous
    d[i + 1] <- distances(CMD, xy[i, ])
    n[i + 1] <- dx
    cells[i + 1] <- nrow(grid)
    i <- i + 1
    dx <- dx + 1
    dy <- dy + 1
  }
  names(xy) <- names(points)
  #Results
  result <- data.frame(centre = xy,
                       dist = d,
                       cells = cells)
  return(result[nrow(result), ])
}

#' @export
#' @rdname centrality
central_feature2pts <- function(points, points2) {
  # Control for few columns in points
  if(ncol(points) < 2) {
    stop("'points' must be of class matrix, data.frame or tibble with 2 or 3 numeric columns")
  }

  # Control for few columns in points2
  if(ncol(points2) < 2) {
    stop("'points2' must be of class matrix, data.frame or tibble with 2 or 3 numeric columns")
  }

  # Control for to many columns in points
  if(ncol(points) > 3) {
    warning("'points' more than 3 columns. Only first three columns are being used as xyz coordinates")
    points <- points[,1:3]
  }

  # Control for to many columns in points2
  if(ncol(points2) > 3) {
    warning("'points2' more than 3 columns. Only first three columns are being used as xyz coordinates")
    points2 <- points2[,1:3]
  }

  # Control for same size in points and points2
  if(ncol(points) != ncol(points2)) {
    if(ncol(points) > ncol(points2)) {
      points <- points[,1:2]
      warning("More dimensions in 'points' than 'points2'. Extra dimensions droped...")
    } else {
      points2 <- points2[,1:2]
      warning("More dimensions in 'points2' than 'points'. Extra dimensions droped...")
    }
  }
  # Central feature between two point patterns
  d <- vector("numeric", nrow(points2))
  i = 1
  while(i <= nrow(points2)){
    d[i] <- sum(distances(destmat = points, centre = points2[i,]))
    i = i + 1
  }

  unlist(points2[which.min(d),])
}
