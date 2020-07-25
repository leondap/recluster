\name{recluster.dist}
\alias{recluster.dist}
\title{
Compute a dissimilarity matrix using a battery of beta-diversity indices
}
\description{
This function computes dissimilarity matrices based on the two most popular partitions of faunistic and phylogenitic beta-diversity. In particular Jaccard = beta3 + richness (Carvalho et al. 2012), Jaccard = Jturnover + Jnestedness (Baselga, 2012) and Sorensen = Simpson + nestedness (Baselga 2010) for faunistic indexes and Unifrac = Unifrac_turn + Unifrac_PD and PhyloSor = PhyloSor_turn + Phylosor_PD (Leprieur et al. 2012). Any other binary index can be included in brackets by using the syntax of designdist function of the vegan package.
}
\usage{
recluster.dist(mat, phylo=NULL, dist="simpson")
}
\arguments{
  \item{mat}{
A matrix containing sites (rows) and species (columns).
}
  \item{phylo}{
An ultrametric and rooted phylogenetic tree for species having the same labels as in mat columns. Only required for phylogenitic beta-diversity indexes.
}
  \item{dist}{
One among the 14 beta-diversity indexes \code{"simpson" "sorensen" "nestedness" "beta3" "richness" "jaccard" "jturnover" "jnestedness" "phylosor" "phylosort" "phylosorpd" "unifrac" "unifract" "unifracpd"}. Any custom binary dissimilarity can also be specified according to the syntax of designdist function of the vegan package.
}
}
\details{
Syntax for binary indices in vegdist: J, number of common species; A and B, number of species exclusive of the first and of the second site. 
}
\value{
An object of class dist (see vegan:designdist for further details)
}
\references{
Baselga A. 
"Partitioning the turnover and nestedness components of beta diversity." 
Global Ecol Biogeogr (2010), 19: 134-143.

Carvalho J. C., Cardoso P., Gomes P. 
"Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns." 
Global Ecol Biogeogr (2012), 21: 760-771. 

Leprieur F., Albouy C., De Bortoli J., Cowman P.F., Bellwood D.R., Mouillot D.  
"Quantifying Phylogenetic Beta Diversity: Distinguishing between 'True' Turnover of Lineages and Phylogenetic Diversity Gradients."
Plos One (2012), 7
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\keyword{ array }
