
biodecrypt.cross <- function(
  mat, id, alpha = NULL, ratio = 2.5, buffer = 90,
  checkdist = TRUE, minimum = 7, polygon = NULL,
  map = NULL, xlim = NULL, ylim = NULL, main = NULL,
  runs = 10, test = TRUE, crs_area = NULL
) {
  
  res <- list()
  res$type <- "cross"
  
  colnames(mat) <- c("Long", "Lat")
  
  unknown <- is.na(id) | id == 0 | id == "0"
  taxa_ids <- sort(unique(id[!unknown]))
  taxa <- length(taxa_ids)
  
  if (taxa < 2) {
    stop("Cross validation requires at least two known taxa.")
  }
  
  id_internal <- rep(0, length(id))
  for (j in seq_along(taxa_ids)) {
    id_internal[id == taxa_ids[j]] <- j
  }
  
  use_known <- which(id_internal > 0)
  idr <- id_internal[use_known]
  matrix <- mat[use_known, , drop = FALSE]
  
  if (is.null(alpha)) {
    alpha <- rep(8, taxa)
  }
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, taxa)
  }
  
  if (length(alpha) != taxa) {
    stop("alpha must have length 1 or the same number of known taxa.")
  }
  
  q <- aggregate(rep(1, length(idr)) ~ idr, FUN = "sum")
  
  if (min(q[, 2]) < 4) {
    stop("The cross validation procedure requires a minimum of 4 distinct points per taxon.")
  }
  
  if (test) {
    test_run <- biodecrypt(
      mat, id,
      alpha = alpha,
      ratio = ratio,
      buffer = buffer,
      map = map,
      polygon = polygon,
      checkdist = checkdist,
      minimum = minimum,
      xlim = xlim,
      ylim = ylim,
      main = main,
      crs_area = crs_area
    )
    
    res$NUR <- test_run$NUR
    res$areas <- test_run$areas
    res$intersections <- test_run$intersections
    res$sympatry <- test_run$sympatry
    res$table <- test_run$table
  }
  
  species <- sort(unique(idr))
  
  which_sp <- vector("list", length(species))
  how_many <- numeric(length(species))
  ord <- NULL
  
  for (sp in seq_along(species)) {
    idx <- which(idr == species[sp])
    which_sp[[sp]] <- sample(idx)
    how_many[sp] <- length(which_sp[[sp]]) / runs
    ord <- c(ord, which_sp[[sp]])
  }
  
  matrixnewbs <- cbind(matrix[ord, , drop = FALSE], idr[ord])
  colnames(matrixnewbs) <- c("Long", "Lat", "id")
  
  attr <- rep(NA, nrow(matrixnewbs))
  
  for (giro in seq_len(runs)) {
    
    memberbs <- matrixnewbs[, "id"]
    
    for (spe in seq_along(species)) {
      
      qualitot <- which(matrixnewbs[, "id"] == species[spe])
      
      first <- floor((giro - 1) * length(qualitot) / runs) + 1
      last <- floor(giro * length(qualitot) / runs)
      
      if (first <= last) {
        memberbs[qualitot[first:last]] <- 0
      }
    }
    
    finalebs <- biodecrypt(
      matrixnewbs[, c("Long", "Lat"), drop = FALSE],
      memberbs,
      alpha = alpha,
      minimum = minimum,
      ratio = ratio,
      buffer = buffer,
      checkdist = checkdist,
      polygon = polygon,
      map = map,
      xlim = xlim,
      ylim = ylim,
      main = main,
      crs_area = crs_area
    )
    
    attri <- which(memberbs == 0)
    
    if (!is.null(finalebs$table_internal)) {
      attr[attri] <- finalebs$table_internal[attri, "id2_internal"]
    } else {
      attr[attri] <- finalebs$table[attri, 3]
    }
  }
  
  res$cross <- cbind(
    original = matrixnewbs[, "id"],
    predicted = attr,
    MIR = rep(0, nrow(matrixnewbs)),
    NIR = rep(0, nrow(matrixnewbs)),
    Long = matrixnewbs[, "Long"],
    Lat = matrixnewbs[, "Lat"]
  )
  
  res$cross[which(res$cross[, "predicted"] > 0 &
                    res$cross[, "original"] != res$cross[, "predicted"]), "MIR"] <- 1
  
  res$cross[which(res$cross[, "predicted"] == 0), "NIR"] <- 1
  
  res$MIR <- (sum(res$cross[, "MIR"]) / nrow(res$cross)) * 100
  res$NIR <- (sum(res$cross[, "NIR"]) / nrow(res$cross)) * 100
  
  res$taxa_ids <- taxa_ids
  
  return(res)
}
