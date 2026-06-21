biodecrypt <- function(mat, id, alpha = NULL, ratio = 2.5, buffer = 90000, 
                       polygon = NULL, checkdist = TRUE, minimum = 7, 
                       plot = FALSE, map = NULL, xlim = NULL, ylim = NULL, 
                       main = NULL, crs_area = NULL) {
  
  res <- list()
  res$type <- "sep"
  
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  
  oldw <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = oldw), add = TRUE)
  
  colnames(mat) <- c("Long", "Lat")
  
  usenondupl <- which(!(duplicated(cbind(mat, id))))
  mat <- mat[usenondupl, , drop = FALSE]
  id <- id[usenondupl]
  
  unknown <- is.na(id) | id == 0 | id == "0"
  taxa_ids <- sort(unique(id[!unknown]))
  taxa <- length(taxa_ids)
  
  if (taxa < 1) stop("No known taxa in id.")
  
  id_internal <- rep(0, length(id))
  for (j in seq_along(taxa_ids)) {
    id_internal[id == taxa_ids[j]] <- j
  }
  
  names(taxa_ids) <- seq_along(taxa_ids)
  
  matpunti <- sf::st_as_sf(
    as.data.frame(mat),
    coords = c("Long", "Lat"),
    crs = 4326
  )
	
  crs_ref <- sf::st_crs(matpunti)
	
  if (!is.null(polygon)) {
    polygon <- sf::st_transform(polygon, sf::st_crs(matpunti))
    polygon <- sf::st_make_valid(polygon)
  }
  

  distances <- matrix(0, nrow(mat), taxa)
  
  if (is.null(alpha)) {
    alpha <- rep(10, taxa)
  }
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, taxa)
  }
  
  if (length(alpha) != taxa) {
    stop("alpha must have length 1 or the same number of known taxa.")
  }
  
  if (is.null(xlim)) xlim <- c(min(mat[, 1]), max(mat[, 1]))
  if (is.null(ylim)) ylim <- c(min(mat[, 2]), max(mat[, 2]))
  
  if (plot) {
    plot(cbind(xlim, ylim), type = "n", main = main)
    if (!is.null(map)) plot(map, add = TRUE)
  }
  
  vectab <- matrix(0, nrow(mat), taxa + 1)
  hulls <- list()
  hullpl <- list()
  areas <- rep(NA, taxa)
  alphaused <- rep(NA, taxa)
  
  for (spec in seq_len(taxa)) {
    
    taxsp <- which(id_internal == spec)
    hullas <- matpunti[taxsp, ]
    
    if (nrow(hullas) >= minimum) {
      
      hullasp <- unique(mat[taxsp, , drop = FALSE])
      hullspat <- NULL
      hull <- NULL
      
      for (increasea in seq(alpha[spec], 20, by = 1)) {
        
        hulltry <- try({
          hull_tmp <- alphahull::ahull(hullasp, alpha = increasea)
		hull2 <- recluster:::ah2sf(hull_tmp)
        
          
          hull_sf <- sf::st_as_sf(hull2)
          sf::st_crs(hull_sf) <- sf::st_crs(matpunti)
          hull_sf <- .recluster_clean_geom(hull_sf, crs_ref = crs_ref)
      
          list(
            hull = hull_tmp,
            hullspat = hull_sf
          )
        }, silent = F)
        
        if (!inherits(hulltry, "try-error")) {
          hull <- hulltry$hull
          hullspat <- hulltry$hullspat
          alphaused[spec] <- increasea
          break
        }
      }
      
      if (is.null(hullspat)) {
        hull <- sf::st_convex_hull(sf::st_make_valid(sf::st_union(hullas)))
        hullspat <- clean_geom(hull)
        alphaused[spec] <- NA
      }
      
    } else {
      
      hull <- sf::st_convex_hull(sf::st_make_valid(sf::st_union(hullas)))
      hullspat <- clean_geom(hull)
      alphaused[spec] <- NA
    }
    
    if (!is.null(polygon)) {
      hullspat <- sf::st_intersection(hullspat, polygon)
      hullspat <- clean_geom(hullspat)
    }
    
    hulls[[spec]] <- hullspat
    hullpl[[spec]] <- hull
    
    if (plot) {
      plot(hullspat, add = TRUE)
      points(hullas)
    }
    
    areas[spec] <- .recluster_area_safe(
  hulls[[spec]],
  crs_ref = crs_ref,
  crs_area = crs_area
   )
    
    pointin <- unique(unlist(sf::st_contains(hulls[[spec]], matpunti)))
    if (length(pointin) > 0) {
      vectab[pointin, spec] <- 1
    }
    
    dist_matrix <- sf::st_distance(matpunti, hullspat)
    distances[, spec] <- apply(dist_matrix, 1, min)
  }
  
  vectab[, ncol(vectab)] <- rowSums(vectab[, seq_len(taxa), drop = FALSE])
  
  id2_internal <- id_internal
  
  uncertain2 <- which(vectab[, ncol(vectab)] == 0 & id_internal == 0)
  inside <- which(vectab[, ncol(vectab)] == 1 & id_internal == 0)
  
  if (length(uncertain2) > 0) {
    
    distancesunc <- distances[uncertain2, , drop = FALSE]
    
    for (k in seq_along(uncertain2)) {
      
      dist <- distancesunc[k, ]
      ord <- order(dist)
      attribution <- ord[1]
      ordereddist <- dist[ord]
      
      if (length(ordereddist) > 1 &&
          ordereddist[2] > buffer &&
          (ordereddist[2] / ordereddist[1]) > ratio) {
        id2_internal[uncertain2[k]] <- attribution
      }
    }
  }
  
  if (length(inside) > 0) {
    
    distancesunc <- distances[inside, , drop = FALSE]
    
    for (k in seq_along(inside)) {
      
      attr <- which(vectab[inside[k], seq_len(taxa)] == 1)
      diste <- distancesunc[k, ]
      diste <- diste[-attr]
      
      if (length(diste) > 0 && min(diste) > buffer) {
        id2_internal[inside[k]] <- attr
      }
    }
  }
  
  if (checkdist) {
    
    check <- which(id2_internal > 0 & id_internal == 0)
    
    if (length(check) > 0) {
      
      known <- which(id_internal > 0)
      use <- id_internal[known]
      
      for (ch in seq_along(check)) {
        
        attrp <- id2_internal[check[ch]]
        
        dist1 <- sf::st_distance(
          matpunti[check[ch], ],
          matpunti[known, ]
        )
        
        mini <- which.min(dist1)
        
        if (attrp != use[mini]) {
          id2_internal[check[ch]] <- 0
        }
      }
    }
  }
  
  intersect <- matrix(NA, taxa, taxa)
  sympatry <- matrix(NA, taxa, taxa)
  
  if (taxa > 1) {
    
    for (k in 1:(taxa - 1)) {
      for (c in (k + 1):taxa) {
        
        inter_area <- .recluster_intersection_area_safe(
  hulls[[k]],
  hulls[[c]],
  crs_ref = crs_ref,
  crs_area = crs_area
)
        
        intersect[k, c] <- inter_area
        intersect[c, k] <- inter_area
        
        denom <- areas[c] + areas[k] - inter_area
        
        if (is.na(denom) || denom <= 0) {
          sympatry[k, c] <- NA
        } else {
          sympatry[k, c] <- inter_area / denom
        }
        
        sympatry[c, k] <- sympatry[k, c]
      }
    }
  }
  
  rownames(intersect) <- colnames(intersect) <- as.character(taxa_ids)
  rownames(sympatry) <- colnames(sympatry) <- as.character(taxa_ids)
  names(areas) <- as.character(taxa_ids)
  names(alphaused) <- as.character(taxa_ids)
  names(hulls) <- as.character(taxa_ids)
  names(hullpl) <- as.character(taxa_ids)
  
  id2 <- rep(NA, length(id2_internal))
  id2[id2_internal == 0] <- if (is.character(id)) "0" else 0
  id2[id2_internal > 0] <- as.character(taxa_ids[id2_internal[id2_internal > 0]])
  
  if (!is.character(id)) {
    suppressWarnings(id2 <- as.numeric(id2))
  }
  
  if (plot) {
    points(matpunti, col = as.numeric(as.factor(id2)), cex = 0.5)
  }
  
  n_unknown <- length(which(id_internal == 0))
  
  res$areas <- areas
  res$intersections <- intersect
  res$sympatry <- sympatry
  res$NUR <- if (n_unknown > 0) {
    length(which(id2_internal == 0)) / n_unknown * 100
  } else {
    NA
  }
  res$table <- cbind(mat, id2 = id2, id = id)
  res$hulls <- hulls
  res$hullpl <- hullpl
  res$alphaused <- alphaused
  res$taxa_ids <- taxa_ids
  res$crs_area <- crs_area
  res$table_internal <- cbind(mat, id2_internal = id2_internal, id_internal = id_internal)
  return(res)
}
