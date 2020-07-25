\name{recluster.plot.matrix}
\alias{recluster.plot.matrix}
\title{
Plot the values of the cells of a matrix in grey scale
}
\description{
This function plots the values of the cells of a matrix in grey scale.
}
\usage{
recluster.plot.matrix(mat)
}
\arguments{
  \item{mat}{
A dissimilarity matrix.
}
}
\value{
A plot of cell values.
}
\references{
Dapporto, L., Fattorini, S., Voda, R., Dinca, V., Vila, R. "Biogeography of western Mediterranean butterflies: combining turnover and nestedness components of faunal dissimilarity." 
J Biogeogr (2014), 41: 1639-1650. 
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\examples{
data(datamod)
simpdiss<- recluster.dist(datamod)
recluster.plot.matrix(simpdiss)
}
\keyword{ dissimilarity }
