recluster.region <- function(mat, tr = 50, dist = "simpson", method = "ward.D2",
                              members = NULL, phylo = NULL,
                              mincl = 2, maxcl = 10,
                              rettree = FALSE, retmat = FALSE, retmemb = FALSE) {
  rows <- nrow(mat)
  clusters <- maxcl - mincl + 1
  mat2 <- as.matrix(mat)

  orig_names <- rownames(mat2)
  if (is.null(orig_names)) orig_names <- as.character(seq_len(rows))

  tab2 <- array(NA_integer_, dim = c(rows, tr, clusters),
                dimnames = list(orig_names, NULL, as.character(mincl:maxcl)))

  if (data.class(mat) == "dist") {
    distam <- mat
  } else {
    distam <- recluster.dist(mat2, phylo = phylo, dist = dist)
  }
  dista <- as.matrix(distam)

  for (i in seq_len(tr)) {
    samp_idx <- sample(seq_len(rows))
    dista2 <- dista[samp_idx, samp_idx, drop = FALSE]

    if (method == "pam") {
      tree <- NULL
    } else {
      if (method == "diana") {
        tree <- cluster::diana(as.dist(dista2))
      } else {
        tree <- hclust(as.dist(dista2), members = members, method = method)
      }
    }

    for (cut in mincl:maxcl) {
      if (method == "pam") {
        cuttr <- cluster::pam(as.dist(dista2), k = cut)$clustering
      } else if (method == "diana") {
        cuttr <- cutree(as.hclust(tree), k = cut)
      } else {
        cuttr <- cutree(tree, k = cut)
      }
      names(cuttr) <- orig_names[samp_idx]
      tab2[, i, (cut - mincl + 1)] <- cuttr[orig_names]
    }
  }

  if (any(is.na(tab2))) {
    bad_rows <- unique(which(apply(is.na(tab2), 1, any)))
    warning("NAs produced in membership table for rows: ",
            paste(orig_names[bad_rows], collapse = ", "),
            ". Check rownames consistency and sampling.")
        stop("NAs found in tab2 — aborting. See warning for problematic rows.")
  }

matrices <- array(NA_real_, dim = c(rows, rows, clusters),
                    dimnames = list(orig_names, orig_names, as.character(mincl:maxcl)))
  for (sel in seq_len(clusters)) {
    tabsel <- tab2[, , sel, drop = FALSE]
    tabsel_df <- as.data.frame(tabsel, stringsAsFactors = TRUE)
    # conversione in factor colonna per colonna (già factors ma ripeto in sicurezza)
    tabsel_df[] <- lapply(tabsel_df, function(x) factor(x, levels = unique(x)))
    matrices[, , sel] <- as.matrix(cluster::daisy(tabsel_df, metric = "gower"))
  }

  res <- list()
  if (retmemb) res$memb <- tab2
  if (retmat) res$matrices <- matrices

  pamsol <- matrix(NA_integer_, nrow = rows, ncol = clusters,
                   dimnames = list(orig_names, as.character(mincl:maxcl)))
  res$solutions <- matrix(NA_real_, nrow = clusters, ncol = 3,
                          dimnames = list(NULL, c("k", "silh", "ex.diss")))
  res$solutions[, "k"] <- seq(mincl, maxcl)

  for (pamr in 1:clusters) {

  if (method == "pam") {
    pamsol[,pamr] <- pam(as.dist(matrices[,,pamr]), k = mincl - 1 + pamr)$clustering

  } else if (method == "diana") {
    pami <- diana(as.dist(matrices[,,pamr]))
    pamsol[,pamr] <- cutree(as.hclust(pami), k = mincl - 1 + pamr)
    if (rettree) res$tree[[pamr]] <- pami

  } else {
    pami <- hclust(as.dist(matrices[,,pamr]), method = method)
    pamsol[,pamr] <- cutree(pami, k = mincl - 1 + pamr)
    if (rettree) res$tree[[pamr]] <- pami
  }

  res$solutions[pamr,3] <- recluster.expl(distam, pamsol[,pamr])
  res$solutions[pamr,2] <- mean(cluster::silhouette(pamsol[,pamr], dist = distam)[,3])
}
res$grouping <- pamsol
  return(res)
}
