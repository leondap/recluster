\name{recluster.node.strength}
\alias{recluster.node.strength}
\title{
Evaluating order row bias in a cluster
}
\description{
This function helps to understand the magnitude of row bias by computing a first tree with the original order of areas. Then it creates a default series of six trees by recluster.cons with increasing consensus rule from 50% to 100% by a default 10% step. Then, it scores the percentage of times each node is replicated among the six runs of consensus thresholds. Low values indicate high row bias.
}
\usage{
recluster.node.strength(mat, phylo = NULL, dist = "simpson", 
nodelab.cex=0.8, tr = 100, levels=6, method = "average", ...)
}
\arguments{
  \item{mat}{
A matrix containing sites (rows) and species (columns).
}
  \item{phylo}{
An ultrametric and rooted phylogenetic tree for species having the same labels as in mat columns. Only required for phylogenitic beta-diversity indexes.
}
  \item{tr}{
The number of trees to be used for the consensus.
}
   \item{dist}{
A beta-diversity index (the Simpson index by default) included in recluster.dist or any custom binary dissimilarity to be specified according to the syntax of designdist function of the vegan package.
}
   \item{nodelab.cex}{
the cex() parameter for controlling the size of the labels on the nodes (see \code{?nodelabels}).
}
   \item{levels}{
The number of levels of different consensus threshold to be used.
}
    \item{method}{
Any clustering method allowed by hclust.
}
  \item{...}{
Arguments to be passed to plot.phylo methods, see the ape package manual and \code{?plot.phylo}.
}
}
\details{
It has to be noted that values obtained by this function are not bootstrap supports for nodes but a crude indication of the magnitude of the row bias. Nodes with low value in this analysis can have strong bootstrap support and vice versa. This preliminary analysis can avoid that the use of a strict consensus (100%), when  the row bias is strong, collapse nodes potentially highly supported by bootstrap. Morever, when the bias is low, it is reasonable to reduce the number of consensus trees in successive analyses, thus fastening the bootstrap procedure. 
}
\value{
A cluster with percentages of recurrence over different consensus runs for each node.
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
recluster.node.strength(datamod, tr=10)
}
\keyword{ cluster }
