  clean_geom <- function(x) {
  x <- sf::st_make_valid(x)
  sf::st_crs(x) <- sf::st_crs(matpunti)
  
  gt <- unique(as.character(sf::st_geometry_type(x)))
  
  if (any(gt == "GEOMETRYCOLLECTION")) {
    x <- sf::st_collection_extract(x, "POLYGON")
  }
  
  x <- x[!sf::st_is_empty(x), ]
  x
}






area_safe <- function(x) {
  x <- clean_geom(x)
  
  if (!is.null(crs_area)) {
    x <- sf::st_transform(x, crs_area)
  }
  
  x <- sf::st_make_valid(x)
  x <- sf::st_union(x)
  
  as.numeric(sf::st_area(x))
}





intersection_area_safe <- function(x, y) {
  x <- clean_geom(x)
  y <- clean_geom(y)
  
  if (!is.null(crs_area)) {
    x <- sf::st_transform(x, crs_area)
    y <- sf::st_transform(y, crs_area)
  }
  
  x <- sf::st_union(sf::st_make_valid(x))
  y <- sf::st_union(sf::st_make_valid(y))
  
  z <- try(sf::st_intersection(x, y), silent = F)
  
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