recluster.expl<-function (dist, clust) {
    mat <- as.vector(dist)
    same <- as.vector(dist(clust))
    beta <- sum(mat, na.rm=T)
    partial <- sum(mat[which(same > 0)], na.rm=T)
    return(partial/beta)
}
