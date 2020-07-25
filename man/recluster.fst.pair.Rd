\name{recluster.fst.pair}
\alias{recluster.fst.pair}
\title{
Compute pairwise indexes of genetic differentiation among populations
}
\description{
This function computes pairwise indexes of genetic differentiation among populations based on a distance matrix and on a vector for populations.
}
\usage{
recluster.fst.pair(dist,vect,setzero=F,setnazero=F)
}
\arguments{
  \item{dist}{
A distance matrix.
}
  \item{vect}{
A vector indicating population membership. Cases must be in the some order of the distance matrix.
}
  \item{setzero}{
A logical indicating if negative values should be set to zero
}
\item{setnazero}{
A logical indicating if NA values should be set to zero
}
}
\details{
The formulas used for pairwise calculations between i and j populations are 
Dstij = Htij - Hsij
Gstij = Dstij/Ht
Dij = (Dstij/(1-Hsij))*2
G'stij = Gstij/((1-Hsij)/(1+Hsij))
see also recluster.fst for a discussion of indexes
}
\value{
\item{Dstm}{
The Dst distance matrix.
}
\item{Gstm}{
The Gst distance matrix.
}
\item{Dm}{
The D distance matrix.
}
\item{G1stm}{
The G'st distance matrix.
}
}
\references{

Jost L. 
"GST and its relatives do not measure differentiation."
Mol Ecol (2008), 17:4015-4026.

Meirmans P. G., Hedrick P. W. 
"Assessing population structure: FST and related measures: Invited Technical Reviwev."
Mol Ecol Res (2011), 11: 5-18.

Nei M. 
Molecular evolutionary genetics (1987), Columbia University Press.

Whitlock M.C. 
"G'ST and D do not replace FST." Mol Ecol (2011), 20: 1083-1091.

}
\author{
Leonardo Dapporto
}
\examples{
datavirtual<-data.frame(replicate(20,sample(0:1,60,rep=TRUE)))
dist<-recluster.dist(datavirtual)
population<-c(rep(1,20),rep(2,20),rep(3,20))
recluster.fst.pair(dist,population)
}
\keyword{differentiation}
