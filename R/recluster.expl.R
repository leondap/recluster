recluster.expl<-function(mat,clust){
mat<-as.matrix(mat)
beta <- sum(mat, na.rm=T)
betapartial <- 0
        for (row in 1:nrow(mat)) {
            for (col in 1:ncol(mat)) {
			if (!(is.na(mat[row,col])) & !(is.na(clust[row])) & !(is.na(clust[col]))) {
               	 	if (clust[row] != clust[col]) {
                  betapartial <- betapartial + mat[row, col]
                }
            }
        }
	}
return(betapartial/beta)
}
