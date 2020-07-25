\name{recluster.plot.col}
\alias{recluster.plot.col}
\title{
Plotting data in RGB space
}
\description{
This function plots a matrix obtained by recluster.col in the RGB space. 
}
\usage{
recluster.plot.col(mat,cext=0.3,cex=1,cex.axis=0.7,cex.lab=0.8,pch=16,text=TRUE,
add=F,xlim=NULL,ylim=NULL,ylab="Axis 2",xlab="Axis 1",...)
}
\arguments{
  \item{mat}{
A matrix inherited by recluster.col.
}
 \item{cext}{
Dimension for labels of row names.
}
 \item{cex}{
Dimension of dots.
}
 \item{cex.axis}{
Dimension of axis labels.
}
 \item{cex.lab}{
Dimension of labels.
}
 \item{text}{
A logical indicating if row names should be plotted.
}
 \item{pch}{
The shape of the dots (See par()).
}
 \item{add}{
A logical indicating if the plot should be added to a precedent graph.
}
 \item{xlim}{
The limit values for x-axis, if NULL the values in the orignal matrix is used.
}
 \item{ylim}{
The limit values for y-axis, if NULL the values in the orignal matrix is used.
}
 \item{ylab}{
The label of the y-axis
}
 \item{xlab}{
The label of the x-axis
}
 \item{...}{
See par() for other graphical parameters
}
}
\value{
A colour plot. 
}	
\references{
Kreft H., Jetz, W. 2010. 
"A framework for delineating biogeographic regions based on species distributions" 
J Biogeogr (2010),37: 2029-2053.

Dapporto, L., Fattorini, S., Voda, R., Dinca, V., Vila, R. "Biogeography of western Mediterranean butterflies: combining turnover and nestedness components of faunal dissimilarity." 
J Biogeogr (2014), 41: 1639-1650. 
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\examples{
data(datamod)
sordiss<- recluster.dist(datamod,dist="sorensen")
points<-cmdscale(sordiss)
col<-recluster.col(points)
recluster.plot.col(col)
}
\keyword{ plot }
