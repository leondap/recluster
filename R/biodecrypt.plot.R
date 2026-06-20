biodecrypt.plot <- function(
  x, minsize = 0.3, pchid = 1, cexid = 0.1, square = 0.001,
  col = c("red", "darkgreen", "blue", "purple"),
  attributed = NULL, hull = TRUE, NUR = "black", fading = 50, ...
) {
  
  recycle_col <- function(cols, n) {
    rep(cols, length.out = n)
  }
  
  if (x$type == "sep") {
    
    data <- as.data.frame(x$table)
    
    if (!is.null(x$table_internal)) {
      data_int <- as.data.frame(x$table_internal)
      id2_plot <- as.numeric(data_int[, "id2_internal"])
      id_plot <- as.numeric(data_int[, "id_internal"])
    } else {
      id2_plot <- suppressWarnings(as.numeric(data[, 3]))
      id_plot <- suppressWarnings(as.numeric(data[, 4]))
    }
    
    newcol <- as.character(paste(data[, 1], data[, 2], sep = "_"))
    data <- cbind(data, newcol)
    data$id2_plot <- id2_plot
    data$id_plot <- id_plot
    
    quanti <- aggregate(rep(1, nrow(data)) ~ data$newcol, FUN = "sum")
    quanti <- quanti[quanti[, 2] > 1, ]
    
    zeri <- which(data$newcol %in% quanti[, 1] & data$id_plot == 0)
    
    if (length(zeri) > 0) {
      data2 <- data[-zeri, ]
    } else {
      data2 <- data
    }
    
    maxid <- max(data2$id2_plot, na.rm = TRUE)
    cols_use <- recycle_col(col, maxid)
    
    matcol <- cbind(
      matrix(1, nrow(data2), 2),
      matrix(NA, nrow(data2), 3)
    )
    
    ids <- sort(unique(data2$id2_plot))
    ids <- ids[!is.na(ids) & ids > 0]
    
    for (sp in ids) {
      wh <- which(data2$id2_plot == sp)
      matcol[wh, 3] <- as.vector(col2rgb(cols_use[sp]))[1]
      matcol[wh, 4] <- as.vector(col2rgb(cols_use[sp]))[2]
      matcol[wh, 5] <- as.vector(col2rgb(cols_use[sp]))[3]
    }
    
    if (is.null(attributed)) {
      attributed <- "fade"
    }
    
    if (attributed == "fade") {
      lower <- which(data2$id_plot == 0 & data2$id2_plot > 0)
      
      if (length(lower) > 0) {
        for (cc in seq_along(lower)) {
          for (colo in 3:5) {
            matcol[lower[cc], colo] <- matcol[lower[cc], colo] +
              ((255 - matcol[lower[cc], colo]) / (100 / fading))
          }
        }
      }
    }
    
    if (hull) {
      hull_obj <- if (!is.null(x$hulls)) x$hulls else x$hullpl
      
      for (sphull in seq_along(hull_obj)) {
        plot(
          hull_obj[[sphull]],
          col = NA,
          border = cols_use[sphull],
          add = TRUE,
          ...
        )
      }
    }
    
    matcol[which(is.na(matcol[, 3])), 3:5] <- col2rgb(NUR)
    
    recluster.plot.pie(
      data2[, 1],
      data2[, 2],
      mat = as.matrix(matcol),
      square = square,
      minsize = minsize,
      proportional = FALSE,
      add = TRUE
    )
    
    if (attributed == "points") {
      
      lat <- as.numeric(data2[, 2])
      long <- as.numeric(data2[, 1])
      
      latsqo <- floor(lat / square) * square
      longsqo <- floor(long / square) * square
      
      newcoord0 <- cbind(
        aggregate(long ~ longsqo + latsqo, FUN = "mean"),
        aggregate(lat ~ longsqo + latsqo, FUN = "mean")[, 3]
      )
      
      newcoord0[, 5] <- paste(newcoord0[, 1], newcoord0[, 2])
      
      lat <- as.numeric(data2[data2$id_plot > 0, 2])
      long <- as.numeric(data2[data2$id_plot > 0, 1])
      
      latsqo <- floor(lat / square) * square
      longsqo <- floor(long / square) * square
      
      newcoord1 <- cbind(
        aggregate(long ~ longsqo + latsqo, FUN = "mean"),
        aggregate(lat ~ longsqo + latsqo, FUN = "mean")[, 3]
      )
      
      newcoord1[, 5] <- paste(newcoord1[, 1], newcoord1[, 2])
      
      newcoord3 <- newcoord0[newcoord0[, 5] %in% newcoord1[, 5], ]
      
      points(newcoord3[, 3:4], cex = cexid, pch = pchid)
    }
  }
  
  if (x$type == "cross") {
    
    data2 <- x$cross
    matcol <- cbind(data2[, c(5, 6)], matrix(NA, nrow(data2), 3))
    
    maxid <- max(data2[, 1], na.rm = TRUE)
    cols_use <- recycle_col(col, maxid)
    
    for (sp in 1:maxid) {
      wh <- which(data2[, 2] == sp)
      matcol[wh, 3] <- as.vector(col2rgb(cols_use[sp]))[1]
      matcol[wh, 4] <- as.vector(col2rgb(cols_use[sp]))[2]
      matcol[wh, 5] <- as.vector(col2rgb(cols_use[sp]))[3]
    }
    
    matcol[which(data2[, 3] == 1), 3:5] <- col2rgb("black")
    matcol[which(data2[, 4] == 1), 3:5] <- col2rgb("white")
    matcol <- matcol[complete.cases(matcol), ]
    
    recluster.plot.pie(
      matcol[, 1],
      matcol[, 2],
      mat = as.matrix(matcol),
      square = square,
      minsize = minsize,
      proportional = FALSE,
      add = TRUE
    )
  }
}
