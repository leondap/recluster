\name{recluster.col}
\alias{recluster.col}
\title{
Projecting a two dimensional plot in RGB space
}
\description{
This function projects a two dimensional matrix into a RGB space with red, green, yellow and blue at its four corners. RGB combination for each case corresponding to its position in this space is provided together with new coordinates.
}
\usage{
recluster.col(mat,st=TRUE,rot=TRUE)
}
\arguments{
  \item{mat}{
A matrix containing two dimensional coordinates for cases.
}
  \item{st}{
Logical, if TRUE then values in axes are standardized between 0 and 1, if FALSE then original values are maintained.
}
  \item{rot}{
Logical, if TRUE then the axis with highest variance is oriented on the x-axis.
}
}
\value{
A matrix with the first two columns representing the coordinates and the third, fourth and fifth representing the red, green and blue components, respectively. 	
}
\references{
Kreft H., Jetz, W. 
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
col
}
\keyword{ plot }
