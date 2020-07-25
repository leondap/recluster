\name{recluster.region}
\alias{recluster.region}
\title{
A new clustering method based on continuous consensus.
}
\description{
This function is specifically designed to allow regionalization analysis when occurrence of zero and tied values are particularly numerous. This often occurs when using turnover indices at small or mid-spatial scales where huge barriers do not occur. The function requires as an input a matrix having areas in rows and species occurrence (1,0) in columns. It also allows for inclusion of a phylogenetic tree to compute phylogenetic beta-diversity. The indices to be used are those allowed by recluster.dist, but custom indices can be introduced (see recluster.dist). Alternatively a dissimilarity metrix obtained with any function can be provided. The function requires the input of a custom number of trees (n=50 by default) and of an interval of mincl-maxcl values (by default 2-3) indicating the number of regions to be identified. Clustering methods implemented by hclust are allowed as well as Partition Around Medoids and DIANA. The ward.2D default method is generally the best performing method, but ward.D, complete clustering, Partition Around Medoids and DIANA can also perform well. The function produces n trees by randomly re-ordering the original row order. Then, the trees are cut to different nodes (from the mincl-1th to the maxcl-1th node), thus producing an increasing number of clusters. Then, the function compares clustering solutions among the same cuts in different re-sampled trees. This produces a dissimilarity matrix among cells represented by the times each pair of areas is located in different clusters over the different solutions extracted for different trees for the same cut. This dissimilarity is standardised by the number of re-sampled trees to produce values ranging from 0 (for pairs of cells always belonging to the same cluster) to 1 (pairs never belonging to the same cluster). A final hierarchical clustering is applied generating an interval of maxcl-mincl. Since the number of clusters requested by the user interval cannot precisely match the mean number of clusters obtained by the tree cuts, the clustering solution for each k value is obtained by selecting the dissimilarity matrix obtained from the nearest mean number of clustering solutions. 
}
\usage{
recluster.region (mat,tr=50,dist="simpson",method="ward.D2",phylo=NULL, mincl=2,maxcl=3,
rettree=FALSE,retmat=FALSE,retmemb=FALSE)}
\arguments{
  \item{mat}{
A binary presence-absence community matrix or any dissimilarity matrix.
}
    \item{tr}{
The number of trees to be included in the consensus.
}
   \item{dist}{
One among the beta-diversity indexes allowed by recluster.dist or a custom binary dissimilarity specified according to the syntax of designdist function of the vegan package. Not required when the input is a dissimilarity matrix.
}
   \item{method}{
Any clustering method allowed by hclust but also "pam" and "diana".
}
   \item{phylo}{
An ultrametric and rooted phylogenetic tree for species having the same labels as in mat columns. Only required for phylogenetic beta-diversity indices.
}
   \item{mincl}{
The minimum number of regions requested
}
   \item{maxcl}{
The maximum number of regions requested
}
   \item{rettree}{
Logical, if TRUE the final trees are returned. 
}
   \item{retmat}{
Logical, if TRUE the new dissimilarity matrices are returned. 
}
   \item{retmemb}{
Logical, if TRUE the memberships for areas in different random trees is returned. 
}
}
\details{
Like other evaluators for goodness of clustering solutions, the funtion provides silhouette values and the explained dissimilarity. The explained dissimilarity (sensu Holt et al. 2013) is represented by the ratio between sums of mean dissimilarities among members of different clusters and the sum of all dissimilarities of the matrix. This value clearly tends to 1 when all areas are considered as independent groups. Silhouette width measures the strength of any partition of objects from a dissimilarity matrix by comparing the minimum distance between each cell and the most similar cell belonging to any other cluster and the mean distance between that cell and the others belonging to the same cluster (see silhouette function in the cluster package). Silhouette values range between -1 and +1, with a negative value suggesting that most cells are probably located in an incorrect cluster. 
}
\value{
\item{memb}{
An array with different matrices indicating for each area (rows) the membership in each random tree (columns) in each cut (matrix).
 }
\item{matrices}{
The new dissimilarity matrices. Up-right cells provided as NAs.
 }
\item{nclust}{
Mean number of clusters among random trees obtained by different cuts.
 }
\item{solutions}{
A matrix providing number of clusters for each solution (k), the associated mean number of clusters obtained by cuts (clust), the silhouette (silh) value and the explained dissimilarity (ex.diss).
 }
\item{grouping}{
A matrix indicating cluster membership of each site in each solution for different numbers of clusters.
 }
}
\references{
Dapporto L. et al.
A new procedure for extrapolating turnover regionalization at mid-small spatial scales, tested on  British butterflies.
Methods Ecol Evol (2015), 6, 1287-1297
}
\author{
Leonardo Dapporto
}
\examples{
data(dataisl)
simpson<-recluster.dist(dataisl)
turn_cl<-recluster.region(simpson,tr=10,rettree=TRUE)
#plot the three for three clusters
plot(turn_cl$tree[[2]])
#inspect cluster membership
turn_cl$grouping
}
\keyword{ cluster }
