\name{recluster.procrustes}
\alias{recluster.procrustes}
\title{
Computes a procrustes analysis between two matrices even if only a subset of cases are shared.
}
\description{
This function computes a procrustes analysis (as done by the vegan procrustes function) but it also allows including a subset of cases shared between the two matrices and some unshared cases. The shared cases must be listed first and in the same order in the two matrices. Moreover, the number of shared cases must be indicated. The function applies a procrustes analysis by scaling, mirroring ad rotating the second matrix to minimizing its dissimilarity from the first on the basis of shared cases. Then, the same transformation is applied to the unshared cases of the second matrix. Finally, it allows including the matrices of coordinates for variables as obtained, for example, by PCA.
}
\usage{
recluster.procrustes(X, Y, Yv=FALSE, num=nrow(X), scale = TRUE, ...)
}
\arguments{
  \item{X}{
Target matrix.
}
  \item{Y}{
Matrix to be rotated.
}
  \item{Yv}{
Matrix of variables for the matrix to be rotated.
}
  \item{num}{
number of shared cases between the target matrix and the matrix to be rotated (by default all).
}
  \item{scale}{
number of shared cases between the target matrix and the matrix to be rotated (by default all).
}
 \item{...}{
See procrustes() for other parameters
}
}
\details{
recluster.procrustes uses the vegan function procrustes to rotate a configuration (Y) to maximum similarity with another target matrix configuration (X) on the basis of a series of shared objects (rows). These objects must be in the same order in the two X and Y matrices. In case of additional cases (rows) in both the X and Y matrices, the same transformation is applied to the case of the Y matrices which are not shared with X. Moreover, the same transformation can be applied to an additional Yv matrix likely representing the coordinates of variables as obtained for example by PCA or other ordination methods. The functions returns an object of the class "procrustes" as implemented in vegan.
}
\value{
\item{Yrot}{
Rotated matrix Y.
 }
\item{X}{
Target matrix.
 }
\item{Yvrot}{
Rotated matrix of variables Yv.
 }
\item{ss}{
Sum of squared differences between X and Yrot on the basis of shared objects.
 }
\item{rotation}{
Orthogonal rotation matrix on the basis of shared objects.
 }
\item{translation}{
Translation of the origin on the basis of shared objects.
 }
\item{scale}{
Scaling factor on the basis of shared objects.
 }
\item{xmean}{
The centroid of the target on the basis of shared objects.
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
#Create and plot a target matrix
ex1 <-rbind(c(1,5),c(5,5),c(3,4),c(3,6))
plot(ex1,col=c(1:4),pch=19,xlim=c(0,6),ylim=c(0,6),cex=2)
#Create and plot a matrix to be rotated. Only the points 1-4 are shared
ex2<-rbind(c(3,1),c(3,3),c(2.5,2),c(3.5,2),c(3,4))
plot(ex2,col=c(1:5),pch=19,xlim=c(0,6),ylim=c(0,6),cex=2)

#Perform the procrustes and plot the matrices
procr1<-recluster.procrustes(ex1,ex2,num=4)
plot(procr1$X,col=c(1:4),pch=19,xlim=c(-4,4),ylim=c(-4,4),cex=2)
plot(procr1$Yrot,col=c(1:5),pch=19,xlim=c(-4,4),ylim=c(-4,4),cex=2)
}
\keyword{ cluster }
