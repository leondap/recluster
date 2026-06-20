biodecrypt.optimise <- function(tab, coef = c(2, 1, 1), penalty = 10) {
  
  res <- list()
  
  if (ncol(tab) < 6) {
    stop("tab must contain at least MIR, NIR, NUR, ratio, buffer and one alpha column.")
  }
  
  tab <- as.data.frame(tab)
  
  MIR <- as.numeric(tab[, 1])
  NIR <- as.numeric(tab[, 2])
  NUR <- as.numeric(tab[, 3])
  
  score <- MIR^coef[1] + NIR^coef[2] + NUR^coef[3]
  
  mini <- min(score, na.rm = TRUE)
  quali <- which((score - mini) < penalty)
  
  if (length(quali) == 0) {
    stop("No valid optimisation rows found.")
  }
  
  if (length(quali) == 1) {
    
    res$ratio <- as.numeric(tab[quali, 4])
    res$buffer <- as.numeric(tab[quali, 5])
    res$MIR <- MIR[quali]
    res$NIR <- NIR[quali]
    res$NUR <- NUR[quali]
    res$alpha <- as.numeric(tab[quali, 6:ncol(tab)])
    
  } else {
    
    score_sel <- score[quali]
    w <- 1 / score_sel
    w <- w / sum(w, na.rm = TRUE)
    
    res$ratio <- sum(as.numeric(tab[quali, 4]) * w, na.rm = TRUE)
    res$buffer <- sum(as.numeric(tab[quali, 5]) * w, na.rm = TRUE)
    
    res$MIR <- mean(MIR[quali], na.rm = TRUE)
    res$NIR <- mean(NIR[quali], na.rm = TRUE)
    res$NUR <- mean(NUR[quali], na.rm = TRUE)
    
    alpha_mat <- as.matrix(tab[quali, 6:ncol(tab), drop = FALSE])
    storage.mode(alpha_mat) <- "numeric"
    
    res$alpha <- colSums(alpha_mat * w, na.rm = TRUE)
  }
  
  res$score <- mini
  res$selected_rows <- quali
  res$penalty <- penalty
  res$coef <- coef
  
  return(res)
}
