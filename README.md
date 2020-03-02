# R package aspace2

Centrographic statistics and computational geometries for spatial point patterns in 2D/3D

## Description:

`aspace2` package is a fork of aspace package from Randy Bui, Ron N. 
Buliung and Tarmo K. Remmel (2012). Like these, aspace2 package contains the
same collection of functions for computing centrographic statistics (e.g., 
standard distance, standard deviation ellipse, standard deviation box) for 
observations taken at point locations in 2D or 3D. The aspace2 library was 
conceived to avoid the default output assigns and file writings. This is 
primarilly intended to be used in a tidy* context, (e.g in `dplyr::groupy_by`,
`dplyr::summarise`, `dplyr::mutate`, etc, contexts).

## Installation:

`aspace2` package is not in CRAN, but cand be installed by `devtools` package
functions.

```
# install.packages(devtools)

# from GitLab (main repository)
devtools::install_gitlab('gavg712/aspace2.git')
```

## Troubleshooting

Any issue found? Please let me know. Report bus, issues, etc. at:

* At GitLab (main repository): https://gitlab.com/gavg712/aspace2/issues
