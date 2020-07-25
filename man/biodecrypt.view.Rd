\name{biodecrypt.view}
\alias{biodecrypt.view}
\title{
A function to view the hulls and manually adjust alpha values before applying biodecrypt
}
\description{
Biodecrypt is absed on hulls chracterized by a concavity determined by an initial alpha 
value which is increased until a given fraction of samples are included in the resulting 
hulls (default 95%). Because different alpha values can better fit the distribution of a 
given cryptic species, the optimal alpha value can be evaluated by series of 
cross-validation analyses using biodecrypt.wrap but also according to the researcher's 
perception. In this case biodecrypt.view provides a rapid representation of the hulls 
based on custom parameters.
}
\usage{
biodecrypt.view(mat,id,alpha=NULL, fraction=0.95, partCount=10, map=NULL,main=NULL, 
xlim=NULL,ylim=NULL,clipToCoast="terrestrial",cex0=0.2,cex1=0.25)
}
\arguments{
  \item{mat}{
A matrix for longitude and latitude (in decimal degrees) of known and unknown cases.
}
   \item{id}{
A vector indicating species membership of each record (in the same order of mat). 
Known records are indicated with 1,2..n representing species membership, unidentified 
records with 0.
}
   \item{alpha}{
A vector indicating an initial alpha value for each species. If NULL a default values of 
8 for all species is used.
}
    \item{fraction}{
The minimum fraction of occurrences that must be included in each hull.
}
\item{partCount}{
The maximum number of disjunct polygons that are allowed..
}
   \item{clipToCoast}{
Either "no" (no clipping), "terrestrial" (only terrestrial part of range is kept) or 
"aquatic" (only non-terrestrial part is clipped).
}
\item{map}{
A map to be plotted during the procedure to show the separation progress.
}
\item{xlim}{
Longitude boudaries for the map. 
}
\item{ylim}{
latitude boudaries for the map. 
}
\item{main}{
The title to be plotted on the graph
}
  \item{cex0}{
The size of dots for non-identified records. 
}
   \item{cex1}{
The size of dots for identified records. 
}
}
\value{
\item{areas}{
The hull areas for all the species (in km squares).
 }
\item{intersections}{
The areas of intersections among hulls for each pair of species.
 }
\item{sympatry}{
The fraction of the overlap area over the total area of the two hulls.
 }
}
\references{
Platania L. et al.
Assigning occurrence data to cryptic taxa improves climatic niche assessments: 
biodecrypt, a new tool tested on European butterflies. Glocal Ecology and Biogeography 
(2020)
}
\author{
Leonardo Dapporto
}
\examples{
# Create an example for a dataset
mat<-rbind(cbind(rnorm(n = 20, mean = 1, sd = 4),rnorm(n = 20, mean = 40, sd = 3)),
cbind(rnorm(n = 20, mean = 7, sd = 5),rnorm(n = 20, mean = 45, sd = 2)))

id<-c(rep(1,20),rep(2,20))
id[sample(c(1:40))[1:10]]<-0

# Make the separation with default parameters 
# alpha gets high to include 95% of defaults cases and hulls are almost convex 
\dontrun{biodecrypt.view(mat,id, clipToCoast="no")}


# Make the separation with custom parameters
# With a lower fraction hulls becomes more concave. Excluded dots works as a punctiform 
# sub-hull in the attribution. 
\dontrun{biodecrypt.view(mat,id, alpha=c(1,5), fraction=0.80, clipToCoast="no")}
}
\keyword{record attribution}
