\name{recluster.expl.diss}
\alias{recluster.expl.diss}
\title{
Cuts a phylogenetic tree and provides cluster membership of areas for custom of all possible clustering solutions and their explained dissimilarity.
}
\description{
This function cuts a phylogenetic tree at all its nodes, and provides membership for each element in the series of resulting clusters and computes the fraction of dissimilarity explained by each solution.
}
\usage{
recluster.expl.diss(tree, dist, maxcl=NULL, mincl=NULL, maxnode=NULL, expld=TRUE) 
}
\arguments{
  \item{tree}{
A phylo tree
}
 \item{dist}{
A dissimilarity matrix.
}
 \item{maxcl}{
A custom number indicating the solution with the minimum number of clusters. If NULL the minimum number of clusters is returned.
}
 \item{mincl}{
A custom number indicated the solution with the maximum number of clusters. If NULL the maximum number of clusters is returned
}
 \item{maxnode}{
A custom number indicated the most external node for the cut. If NULL all the nodes will be cut
}
\item{expld}{
A logical. If TRUE then the matrix for explained dissimilarity is computed.
}
}
\details{
When polytomic nodes are involved in a cut the number of clusters at that cut could increase more than one unit. It is also possible that at the first cut more than two cluster are identified, it is thus possible to obatin a first solution showing a higher number of clusters then the miminum number included in mincl. Holt at al. (2013) identified levels of explained dissimilarity to be used as a reliable threshold to assess a tree cut. When cases are highly numerous maxnode can be set in order to avoid a very long computation keeping in mind that a cut at node 6 can produce solutions with >6 clusters
}
\value{
 \item{matrix}{
A matrix indicating cluster membership of each site in each cut of the tree.
 }
 \item{expl.div}{
A vector indicating the explained dissimilarity for each cut.
 }
  \item{nclust}{
A vector indicating the number of clusters resulting from each cut.
 }
}	
\references{
Dapporto L., Ciolli G., Dennis R.L.H., Fox R., Shreeve, T.G. "A new procedure for extrapolating turnover regionalization at mid?small spatial scales, tested on B ritish butterflies."
Methods in Ecology and Evolution (2015), 6:1287-1297.
}
\author{
Leonardo Dapporto
}
\examples{
data(datamod)
sor_tree<- recluster.cons(datamod, dist="sorensen")
sor_diss <- recluster.dist (datamod, dist="sorensen")
expl_diss <- recluster.expl.diss (sor_tree$cons,sor_diss)
expl_diss 
}
\keyword{ cluster }
