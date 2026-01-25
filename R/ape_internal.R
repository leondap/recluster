# Copied from ape internal helpers 

.recluster_compressTipLabel<-function (x, ref = NULL) 
{
    if (!is.null(attr(x, "TipLabel"))) 
        return(x)
    if (is.null(ref)) 
        ref <- x[[1]]$tip.label
    n <- length(ref)
    if (length(unique(ref)) != n) 
        stop("some tip labels are duplicated in tree no. 1")
    relabel <- function(y) {
        label <- y$tip.label
        if (!identical(label, ref)) {
            if (length(label) != length(ref)) 
                stop("one tree has a different number of tips")
            ilab <- match(label, ref)
            if (any(is.na(ilab))) 
                stop("one tree has different tip labels")
            ie <- match(1:n, y$edge[, 2])
            y$edge[ie, 2] <- ilab
        }
        y$tip.label <- NULL
        y
    }
    x <- unclass(x)
    x <- lapply(x, relabel)
    attr(x, "TipLabel") <- ref
    class(x) <- "multiPhylo"
    x
}

.recluster_uncompressTipLabel<-function (x) 
{
    Lab <- attr(x, "TipLabel")
    if (is.null(Lab)) 
        return(x)
    clx <- class(x)
    class(x) <- NULL
    for (i in seq_along(x)) x[[i]]$tip.label <- Lab
    class(x) <- clx
    attr(x, "TipLabel") <- NULL
    x
}
