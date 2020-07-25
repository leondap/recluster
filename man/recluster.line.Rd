\name{recluster.line}
\alias{recluster.line}
\title{
Identifies a line in a configuration and computes its intercept and angular coefficient}
\description{
This function identifies a line in a configuration based on different criteria and produces its slope and intercept values. I can be used together with recluster.rotate to rotate a configuration based on a custom line.
}
\usage{
recluster.line(mat,type="maxd",X1=NULL,X2=NULL)
}
\arguments{
  \item{mat}{
The bidimensional configuration.
}
  \item{type}{
The type of line to be computed: "maxd" is the line connecting the most distant points, "regression" is the regression line between X and Y values, "points" is the line connecting two custom points of the configuration (X1 and X2).
}
  \item{X1}{
The row number in mat of the first custom point.
}
  \item{X2}{
The row number in mat of the second custom point.
}
}
\value{
\item{m}{
The slope of the line.
 }
\item{q}{
The intercept of the line.
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
#Compute bidimensional representation for islands
pcoa<-cmdscale(recluster.dist(dataisl))
#Compute the line
lin<-recluster.line(pcoa)
}
\keyword{configurations}
