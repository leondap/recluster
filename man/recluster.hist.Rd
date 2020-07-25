\name{recluster.hist}
\alias{recluster.hist}
\title{
Histogram of dissimilarity with tied and zero values
}
\description{
This function creates a histogram with the values of a dissimilarity matrix where the number of cells with zero value are explicitely showed in the first bar. Moreover, it provides the percentage of cells having equal values in the matrix.
}
\usage{
recluster.hist(x)
}
\arguments{
  \item{x}{
A dissimilarity matrix.
}
}
\value{
An histogram with supplementary information. The first bar only shows the zero values.
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
simpdiss<- recluster.dist(datamod)
recluster.hist(simpdiss)
}
\keyword{ dissimilarity }
