\name{recluster.fst}
\alias{recluster.fst}
\title{
Compute some indexes of genetic differentiation
}
\description{
This function computes some indexes of genetic differentiation based on a distance matrix and on a vector for populations.
}
\usage{
recluster.fst(dist,vect,setzero=F,setnazero=F)
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
There has been a large dabate around FST like indexes. Two main indexes are culcalated by this function: the absolute differentiation (Dst) and the standardized differentiation (Gst) (Nei, 1987) .Dst is calculated as: 
Dst = Ht - Hs
where Ht represents the average distances among all the specimens in the sample, and Hs is the average of the intra-area (or intra-sub-area) distances. Thus, Dst represents the average genetic differentiation among areas in p-distance units. Gst is a standardized index defined as:
Gst = Dst/Ht
representing the fraction of the total genetic differentiation encompassed by the differentiation among areas (Nei, 1987). This index ranges from negative values to 1 (complete differentiation). 
Negative values in Gst and Dst (intra-area differentiation higher than inter-area differentiation) can have different subtle meanings, but are most often generated as bias due to relatively small sample sizes; usually they are set to zero (Meirmans & Hedrick, 2011) and we applied this solution. 
In the species showing no mutations in the sample, Gst returns a NA value (while Dst equals to zero). These cases can be also set to zero
The use of Dst and Gst has been debated as a measure of population diversification for extremely variable markers (as micro-satellites) as it tends to underestimate differentiation among populations and to strongly depend on intra-population variability (Jost, 2008; Whitlock, 2011). D and G-st indices are less affected by high values of Hs 
}
\value{
\item{Ht}{
The average distances among all the specimens in the sample.
}
\item{lengthHt}{
The number of distances among all the specimens in the sample.
}
\item{Hs}{
The average distances among the specimens of the same populations.
}
\item{lengthHs}{
The number of distances among the specimens of the same populations.
}
\item{Dst}{
The Dst value.
}
\item{Gst}{
The Gst value.
}
\item{D}{
The D value.
}
\item{G1st}{
The G'st value.
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
datavirtual<-data.frame(replicate(10,sample(0:1,30,rep=TRUE)))
dist<-recluster.dist(datavirtual)
population<-c(rep(1,20),rep(2,20),rep(3,20))
recluster.fst(dist,population)
}
\keyword{differentiation}
