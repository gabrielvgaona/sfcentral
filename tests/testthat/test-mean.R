test_that("mean center works", {
  x <- sf::st_sf(geometry = sf::st_sfc(list(sf::st_point(c(1, 1)),
                                            sf::st_point(c(2, 2)),
                                            sf::st_point(c(3, 3)))))
  y <- st_coordinates(x)
  expect_equal(st_coordinates(st_central_point(x)),
               matrix(apply(y, 2, mean), ncol = 2,
                      dimnames = list(1, c("X", "Y"))))
  expect_equal(nrow(st_central_point(x)), 1)
  expect_error(st_central_point(y))
  expect_s3_class(st_central_point(x), c("sf", "data.frame"))
})
