\name{recluster.test.dist}
\alias{recluster.test.dist}
\title{
Test variation lost by a bidimensional configuration when the coordinates of the elements are reduced to the configuration of the barycentres of a given series of groups.
}
\description{
This function evaluates the amount of variation maintained by a bidimensional configuration after the elements are reduced to the barycentres according to a grouping variable. If elements of different groups are randomly scattered in the configuration, almost all barycentres are expected to attain a rather central position with respect to the original elements, which would result in a small mean distance between barycentres. Conversely, if the elements of different groups are strictly clustered in the representation, the distances among barycentres are expected to be similar to the distances among original elements.
}
\usage{
recluster.test.dist(mat1,mat2,member,perm=1000,elev=2)
}
\arguments{
  \item{mat1}{
The bidimensional configuration before computing barycentres for groups.
}
 \item{mat2}{
The bidimensional configuration after computing barycentres for groups.
}
 \item{member}{
A vector indicating group membership for each element.
}
 \item{perm}{
The number of permutations.
}
\item{elev}{
The power of distances (by default 2:squared distances).
}
}
\details{
The function produces a ratio between the mean squared pairwise distance for all elements and the mean squared pairwise distance for barycentres. This ratio is calculated for the overall configuration and for the two axes separately. The function also provides a test for the significance of the variation preserved by barycentres by creating a custom number of matrices (1000 by default) by randomly sampling the original vector defining groups. Then it computes the frequency of mean squared distance ratios in random configurations higher than the observed ratio. 
}
\value{
 \item{ratio}{
The ratio between mean distances among original elements and barycentres over the overall configuration.
 }
 \item{ratioX}{
The ratio between mean distances among original elements and barycentres on the X axis.
 }
  \item{ratioY}{
The ratio between mean distances among original elements and barycentres on the Y axis.
 }
 \item{test}{
The permutation test for variation maintained over the overall configuration.
 }
 \item{testX}{
The permutation test for variation maintained along the X axis.
 }
  \item{testY}{
The permutation test for variation maintained along the Y axis.
 }
}	
\references{
Dapporto L., Voda R., Dinca V., Vila R.
"Comparing population patterns for genetic and morphological markers with uneven sample sizes. An example for the butterfly Maniola jurtina"
Methods Ecol Evol (2014), 5, 834-843.
}
\author{
Leonardo Dapporto
}
\examples{
data(dataisl)
#Define groups of islands
memb<-c(2,3,5,7,5,3,1,1,2,5,1,3,1,1,5,2,2,1,2,4,1,3,1,5,2,1,7,6,1,1,1) 
#Compute bidimensional representation for elements
pcoa<-cmdscale(recluster.dist(dataisl))
bar<-aggregate(pcoa~memb,FUN="mean")[,2:3]
# test if the variation has been significantly lost
recluster.test.dist(pcoa,bar,memb,perm=100)
}
\keyword{ cluster }
