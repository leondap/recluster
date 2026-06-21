.recluster_clean_geom <- function(x, crs_ref) {
  
  x <- sf::st_make_valid(x)
  sf::st_crs(x) <- crs_ref
  
  gt <- unique(as.character(sf::st_geometry_type(x)))
  
  if (any(gt == "GEOMETRYCOLLECTION")) {
    x <- sf::st_collection_extract(x, "POLYGON")
  }
  
  x <- x[!sf::st_is_empty(x), ]
  
  x
}





.recluster_area_safe <- function(x, crs_ref, crs_area = NULL) {
  
  x <- .recluster_clean_geom(x, crs_ref = crs_ref)
  
  if (!is.null(crs_area)) {
    x <- sf::st_transform(x, crs_area)
  }
  
  x <- sf::st_make_valid(x)
  x <- sf::st_union(x)
  
  as.numeric(sf::st_area(x))
}




.recluster_intersection_area_safe <- function(x, y, crs_ref, crs_area = NULL) {
  
  x <- .recluster_clean_geom(x, crs_ref = crs_ref)
  y <- .recluster_clean_geom(y, crs_ref = crs_ref)
  
  if (!is.null(crs_area)) {
    x <- sf::st_transform(x, crs_area)
    y <- sf::st_transform(y, crs_area)
  }
  
  x <- sf::st_union(sf::st_make_valid(x))
  y <- sf::st_union(sf::st_make_valid(y))
  
  z <- try(sf::st_intersection(x, y), silent = TRUE)
  
  if (inherits(z, "try-error") || is.null(z) || all(sf::st_is_empty(z))) {
    return(0)
  }
  
  z <- sf::st_make_valid(z)
  
  gt <- unique(as.character(sf::st_geometry_type(z)))
  
  if (any(gt == "GEOMETRYCOLLECTION")) {
    z <- sf::st_collection_extract(z, "POLYGON")
  }
  
  z <- z[!sf::st_is_empty(z)]
  
  if (length(z) == 0 || all(sf::st_is_empty(z))) {
    return(0)
  }
  
  z <- sf::st_union(z)
  
  as.numeric(sf::st_area(z))
}
