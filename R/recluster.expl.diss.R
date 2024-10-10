recluster.expl.diss<-function (tree, dist, maxcl = NULL, mincl=NULL, maxnode=NULL, expld=TRUE) { 
	nclust <- NULL
	result <- NULL
	res <- NULL
	tree <- reorder(tree)
	mat <- nodeHeights(tree)
	if(is.null(mincl)){mincl<-1}
	if(is.null(maxcl)){maxcl<-nrow(mat)-1}
	mat2 <- rbind(c(0,1),mat[order(mat[, 1], mat[, 2]), ])
	mat2 <- mat2[mat2[, 1] + mat2[, 2] != 0, ]
	mat2 <- mat2[!duplicated(round(mat2[, 1],5)), ]
	if(is.null(maxnode)){maxnode<-(nrow(mat2)-1)}
	if(maxnode>nrow(mat2)-1){maxnode<-(nrow(mat2)-1)}
	mat2<-mat2[2:(maxnode+1),]
	matrix <- matrix(data = NA, ncol = nrow(mat2), nrow = length(tree$tip.label))
	comp <- rownames(as.matrix(dist))
	for (cl in 1:nrow(mat2)) {
		res <- treeSlice(tree, mat2[cl, 1] - 1e-06, trivial = TRUE)
		sub <- length(res)
		nclust[cl] <- sub
		for (subtrees in 1:sub) {
			tips<-res[[subtrees]]$tip.label
			matrix[match(tips, comp), cl] <- subtrees                           
        	}
    	}
    	cluster <- NULL
	if(expld){
		for (loops in 1:(nrow(mat2))) {
			cluster[loops] <-recluster.expl(dist, as.numeric(matrix[,loops]))
		}
      }
	select<-which(nclust>=mincl & nclust<= maxcl)   
	rownames(matrix) <- comp
	result$matrix <- matrix[,select]
	result$expl.div <- cluster[select]
	result$nclust <- nclust[select]
	return(result)
}

