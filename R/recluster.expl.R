recluster.expl<-function(dist,clust){
	mat<-as.vector(dist)
	same<-as.vector(dist(clust))
	beta <- sum(mat)
	partial<-sum(mat[which(same>0)])
	return(partial/beta)
}
