#' @title Trigonometric functions to compute in degrees
#' @author Gabriel Gaona
#' @description To simplify the use of radians in R
#' @param x Numeric vector of angles in degrees
#' @return Numeric vector of trigonometric computation
#' @keywords trigonometry, degrees, as_radians

#' @export
#' @rdname trigonometric
sin_d <- function(x) return(sin(x * pi / 180))

#' @export
#' @rdname trigonometric
cos_d <- function(x) return(cos(x * pi / 180))

#' @export
#' @rdname trigonometric
tan_d <- function(x) return(tan(x * pi / 180))

#' @export
#' @rdname trigonometric
asin_d <- function(x) return(asin(x * 180 / pi))

#' @export
#' @rdname trigonometric
acos_d <- function(x) return(acos(x * 180 / pi))

#' @export
#' @rdname trigonometric
atan_d <- function(x) return(atan(x * 180 / pi))

#' @export
#' @rdname trigonometric
as_radians <- function(x) return(x * pi / 180)
