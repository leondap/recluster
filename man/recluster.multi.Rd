\name{recluster.multi}
\alias{recluster.multi}
\title{
Multiscale bootstrap based on a consensus tree
}
\description{
Given an initial tree and a data matrix, this function computes bootstrap for nodes as done by recluster.boot. Different levels of bootstrap can be computed by varying the proportions of species sampled from the original matrix.
}
\usage{
recluster.multi(tree, mat, phylo = NULL, tr = 100, p = 0.5, 
dist = "simpson", method = "average", boot = 1000, levels = 2, step = 1)
}
\arguments{
  \item{tree}{
A reference phylo tree for sites presumably constructed with recluster.cons function.
}
  \item{mat}{
The matrix used to construct the tree.
}
  \item{phylo}{
An ultrametric and rooted phylo tree for species having the same labels as in mat columns. Only required for phylogenetic beta-diversity indexes.
}
  \item{tr}{
The number of trees to be included in the consensus.
}
  \item{p}{
A numeric value between 0.5 and 1 giving the proportion for a clade to be represented in the consensus tree.
}
  \item{dist}{
One among the twelve beta-diversity indexes \code{"simpson" "sorensen" "nestedness" "beta3" "richness" "jaccard" "phylosor" "phylosort" "phylosorpd" "unifrac" "unifract" "unifractpd"}. Any custom binary dissimilarity can also be specified according to the syntax of designdist function of the vegan package.
}
  \item{method}{
Any clustering method allowed by hclust.
}
  \item{boot}{
The number of trees used for bootstrap computation.
}
  \item{levels}{
The number of levels to be used in multiscale bootstrap.
}
  \item{step}{
The increase in ratio between the first level (x1) and the next ones.
}
}
\details{
Computation can be time consuming. It is suggested to assess the degree of row bias by recluster.hist and recluster.node.strength to optimize the number of consensus trees before starting the analysis.
}
\value{
A matrix indicating the percentage of bootstrap trees replicating each node for each level.
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
data(datamod)
tree<-recluster.cons(datamod,tr=10)
multiboot<-recluster.multi(tree$cons,tr=10,boot=50,datamod,levels=2,step=1)
recluster.plot(tree$cons,multiboot,1,2,direction="downwards")
}

\keyword{ cluster }
