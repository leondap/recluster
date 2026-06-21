biodecrypt.wrap <- function(
  mat, id,
  alpha = c(1, 5, 10, 15),
  alphamat = NULL,
  ratio = c(2, 3, 4, 5),
  buffer = c(0, 40, 80, 120, 160),
  checkdist = TRUE,
  polygon = NULL,
  minimum = 7,
  map = NULL,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  save = TRUE,
  name = "res_cross.txt",
  runs = 10,
  crs_area = NULL
) {
  
  res <- list()
  
  unknown <- is.na(id) | id == 0 | id == "0"
  taxa_ids <- sort(unique(id[!unknown]))
  taxa <- length(taxa_ids)
  
  if (taxa < 2) {
    stop("At least two known taxa are required.")
  }
  
  if (!is.null(alphamat)) {
    alphamat <- as.matrix(alphamat)
    
    if (nrow(alphamat) != taxa) {
      stop("alphamat must have one row for each known taxon.")
    }
    
    al <- ncol(alphamat)
  } else {
    al <- length(alpha)
  }
  
  rat <- length(ratio)
  buf <- length(buffer)
  
  n_alpha_cols <- if (is.null(alphamat)) 1 else taxa
  
  res_cross <- matrix(
    NA,
    nrow = al * rat * buf,
    ncol = 5 + n_alpha_cols
  )
  
  riga <- 1
  
  for (alphav in seq_len(al)) {
    for (ratiov in seq_len(rat)) {
      for (bufferv in seq_len(buf)) {
        
        if (is.null(alphamat)) {
          alphause <- rep(alpha[alphav], taxa)
        } else {
          alphause <- alphamat[, alphav]
        }
        
        print(c(alphav, ratiov, bufferv))
        
        cross <- biodecrypt.cross(
          mat, id,
          ratio = ratio[ratiov],
          buffer = buffer[bufferv],
          alpha = alphause,
          checkdist = checkdist,
          minimum = minimum,
          map = map,
          polygon = polygon,
          xlim = xlim,
          ylim = ylim,
          main = main,
          runs = runs,
          test = TRUE,
          crs_area = crs_area
        )
        
        res_cross[riga, 1] <- cross$MIR
        res_cross[riga, 2] <- cross$NIR
        res_cross[riga, 3] <- cross$NUR
        res_cross[riga, 4] <- ratio[ratiov]
        res_cross[riga, 5] <- buffer[bufferv]
        if (is.null(alphamat)) {
  res_cross[riga, 6] <- alpha[alphav]
} else {
  res_cross[riga, 6:ncol(res_cross)] <- alphause
}
        
        if (save) {
          write.table(
            res_cross,
            file = name,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
          )
        }
        
        riga <- riga + 1
      }
    }
  }
  
  colnames(res_cross) <- c(
    "MIR", "NIR", "NUR", "ratio", "buffer",
    paste0("alpha_", seq_len(n_alpha_cols))
  )
  
  res$table <- res_cross
  res$taxa_ids <- taxa_ids
  res$crs_area <- crs_area
  
  return(res)
}
