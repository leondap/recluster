\name{recluster.check}
\alias{recluster.check}
\title{
An internal workaround for the new ape package
}
\description{
This function solve a label problem producing unreliable low values of bootstrap
}
\usage{
recluster.check(tree, tip)
}

\arguments{
  \item{tree}{
An object of class phylo
}
\item{tip}{
A tip name
}
}
\value{
An oject of the class phylo
}
\references{
Dapporto L., Ramazzotti M., Fattorini S., Talavera G., Vila R., Dennis R. 
"recluster: an unbiased clustering procedure for beta-diversity turnover"
Ecography (2013), 36:1070-1075.
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\keyword{ internal }
