\name{recluster.rotate}
\alias{recluster.rotate}
\title{
Rotates a bidimensional configuration according to a line
}
\description{
This function rotates the points of a configuration to a new configuration where a line identified by its intercept and its angular coefficient is rotated to become horizontal. The function can also flip or centre a configuration
}
\usage{
recluster.rotate(table,m=FALSE,q=FALSE,flip="none",centre=TRUE)
}
\arguments{
  \item{table}{
The bidimensional configuration.
}
  \item{m}{
The line slope.
}
  \item{q}{
The line intercept
}
  \item{flip}{
The kind of flip, no flip, "none"; "hor", flip horizontally; "ver", flip vertically; "both", flip vertically and horizontally.
}
\item{centre}{
A logical. If TRUE the configuration, after transformation is centered to the mean X and Y values.
}
}
\value{
\item{table2}{
The transformed bidimensional configuration.
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
plot (pcoa)
#Compute the line
lin<-recluster.line(pcoa)
transf<-recluster.rotate(pcoa,m=lin$m,q=lin$q)
plot(transf)
}
\keyword{configurations}
