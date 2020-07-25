\name{recluster.cons}
\alias{recluster.cons}
\title{
Consensus tree among re-sampled trees 
}
\description{
This function creates a series of trees by resampling the order of sites in the original dissimilarity matrix. Then, it computes a consensus among them. The resulting tree is independent of the original row order.
}
\usage{
recluster.cons(mat, phylo = NULL, tr = 100, p = 0.5, 
dist = "simpson", method = "average", blenghts=TRUE, select=FALSE)
}
\arguments{
  \item{mat}{
A matrix containing sites (rows) and species (columns) or any dissimilarity matrix.
}
  \item{phylo}{
An ultrametric and rooted tree for species phylogeny having the same labels as in mat columns. Only required to compute phylogenitic beta-diversity indexes.
}
  \item{tr}{
The number of trees to be used for the consensus.
}
  \item{p}{
A numeric value between 0.5 and 1 giving the proportion for a clade to be represented in the consensus tree.
}
  \item{dist}{
A beta-diversity index (the Simpson index by default) included in recluster.dist or any custom binary dissimilarity to be specified according to the syntax of designdist function of the vegan package.
}
   \item{method}{
Any clustering method allowed by hclust.
}
   \item{blenghts}{
A logical indicating if non-negative least squares branch lengths should be computed.
}
   \item{select}{
A logical indicating if only trees having a fit higher than the median value in the least squares regression should be included in the consensus analysis.
}
}
\details{
According to the primitive "consensus" function from the "ape" package, p must range between 0.5 and 1. Select = TRUE can allow lowering polytomies by removing trees with topology showing particularly low correlation with the distance matrix. Row names are required.
}
\value{
	\item{cons}{
The consensus tree, an object of class phylo.
	}
	\item{trees}{
The trees used to construct the final consensus tree.
	}
	\item{RSS}{
The Residual Sum of Squares for the trees resulting if select=TRUE.
	}
}
\references{
Dapporto L., Ramazzotti M., Fattorini S., Talavera G., Vila R., Dennis R. 
"recluster: an unbiased clustering procedure for beta-diversity turnover"
Ecography (2013), 36:1070-1075.
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\examples{
#Faunistic beta diversity
data(datamod,treemod)
tree<-recluster.cons(datamod,tr=10)
plot(tree$cons,direction="downwards")

#Phylogenetic beta diversity
tree_p<-recluster.cons(datamod,treemod,dist="phylosort",tr=10)
plot(tree_p$cons, direction="downwards")
}
\keyword{ cluster }
