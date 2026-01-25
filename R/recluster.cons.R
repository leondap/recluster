recluster.cons <- function(mat, phylo = NULL, tr = 100, p = 0.5,
                           dist = "simpson", method = "average",
                           blenghts = TRUE, select = FALSE,
                           seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # distanza
  if (inherits(mat, "dist")) {
    distance <- mat
  } else {
    distance <- recluster.dist(mat, phylo, dist)
  }

  dist_mat <- as.matrix(distance)
  n <- nrow(dist_mat)
  idx <- seq_len(n)

  trees <- vector("list", tr)
  RSS <- if (select) numeric(tr) else NULL

  for (i in seq_len(tr)) {
    samp <- sample(idx)
    dist1 <- dist_mat[samp, samp]

    tree <- as.phylo(hclust(as.dist(dist1), method = method))
    trees[[i]] <- tree

    if (select) {
      RSS[i] <- attr(
        nnls.tree(as.dist(dist1), tree, rooted = TRUE),
        "RSS"
      )
    }
  }

  if (select) {
    keep <- RSS < median(RSS, na.rm = TRUE)
    trees <- trees[keep]
    RSS <- RSS[keep]
  }

  cons <- compute.brlen(
    consensus(trees, p = p, check.labels = TRUE),
    method = "Grafen"
  )

  if (blenghts) {
    cons <- tryCatch(
      nnls.tree(distance, cons, rooted = TRUE, trace = FALSE),
      error = function(e) cons
    )
  }

  res <- list(
    cons = multi2di(cons, random = TRUE),
    trees = trees,
    RSS = RSS
  )

  return(res)
}
