\name{biodecrypt.cross}
\alias{biodecrypt.cross}
\title{
Perform a cross validation analysis to test the attribution of biodecrypt on attributed
 records
}
\description{
The function biodecrypt.cross wraps the biodecrypt function to carry out 
cross-validation of known cases thus verifying the robustness of the attribution of 
unknown cases. This function requires the same input of biodecrypt (coordinates and 
vector with attribution together with values of distance ratio, buffer and alpha). 
Moreover this function requires a "runs" value defining the number of different runs, 
thus the fraction of test records included in each run. In each run, randomly selected 
group of test records (actually identified to a given species) are regarded as 
unidentified (0 value) and the biodecrypt function is carried out to attribute them.
The analysis is repeated as often as defined in runs (a runs value of 10 will perform 
a ten-fold cross-validation based on the initial selection of ten randomly distributed 
subsets).
}
\usage{
biodecrypt.cross(mat,id,alpha=NULL,ratio=2.5,buffer=90,fraction=0.95, partCount = 10, 
checkdist=T, clipToCoast="terrestrial", proj = "+proj=longlat +datum=WGS84",minimum=7,
map=NULL,xlim=NULL,ylim=NULL,main=NULL,runs=10,test=T)
}
\arguments{
  \item{mat}{
A matrix for longitude and latitude (in decimal degrees) for all records.
}
   \item{id}{
A vector indicating species membership of each record (in the same order of mat). 
Identified records are indicated with 1,2..n, unidentified records with 0.
}
   \item{alpha}{
A vector indicating an initial alpha value for each species. If NULL, the default value
 of 8 for all species is used.
}
   \item{ratio}{
The minimum ratio between the distance from the second distant hull compared to the 
closest hull to allow attribution. Default 2.5 times.
}
   \item{buffer}{
A distance buffer from hulls(in km).
}
   \item{fraction}{
The minimum fraction of occurrences that must be included in polygon.
}
\item{partCount}{
The maximum number of disjunct polygons that are allowed..
}
   \item{clipToCoast}{
Either "no" (no clipping), "terrestrial" (only terrestrial part of range is kept) or 
"aquatic" (only non-terrestrial part is clipped).
}
   \item{checkdist}{
Logical, if TRUE cases attributed to a given species based on relative distance from 
hulls but closer to an identified record of another species are not attributed. 
}
   \item{proj}{
the projection information for mat. In this version, the default is the only supported 
option.
}
   \item{minimum}{
The minimum number of specimens required to build alpha hulls. If the number of 
identified specimens is lower, convex hulls are calculated to improve procedure 
stability.
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
The name to be plotted on the graph
}
\item{runs}{
The number of runs among which the cases are randomly assigned as non-attributed 
records
}
\item{test}{
A logical, if TRUE, a biodecrypt analysis is also carried out to compute NUR.
}
}
\details{
 The procedure attributes the subsets of identified records to the test group (unknown 
cases) as evenly as possible among runs both in terms of total number of test records and records belonging to the same original species. If the number of runs equates the number of records, then each identified record is individually attributed in a jackknife procedure. Subsequently, the attribution vector obtained is provided and compared with the original membership and two values are provided: the percentages of identified cases attributed to a wrong species (Mis-Identified Records, MIR) and the percentage of known cases not attributed to any species (Non-attributed Identified Records, NIR). The function also has an option to calculate the percentage of Non-attributed Unidentified Records (NUR) representing the fraction of unknown records that could not be attributed to a species after a typical biodecrypt analysis using the parameters provided by the user and the complete set of records. 
}
\value{
\item{type}{
"cross" an argument to be passed to biodecrypt.plot.
 }
\item{NUR}{
The percentage of Non-attributed Unidentified Records.
 }
\item{areas}{
The hull areas for all the species (in km squares).
 }
\item{intersections}{
The areas of intersections among hulls for each pair of species.
 }
\item{sympatry}{
The fraction of the overlap area compared to the total area of the two hulls.
 }
\item{table}{
The result table of the test (if test=TRUE) with Longitude and Latitude for each 
occurrence datum, its id after the biodecrypt procedure (id2) and its initial id (id).
 }
\item{cross}{
The result table with the original attribution (original), the attribution obtained 
after cross validation (predicted) and the classification as MIR or NIR. Longitude and 
Latitude are also provided.
 }
\item{MIR}{
The percentage of Mis-Identified Records.
 }
\item{NIR}{
The percentage of Non-Identified Records.
 }
}
\references{
Platania L. et al.
Assigning occurrence data to cryptic taxa improves climatic niche assessments: 
biodecrypt, a new tool tested on European butterflies. Glocal Ecology and Biogeography 
(2020).
}
\author{
Leonardo Dapporto
}
\examples{
## Not Run
## Create an example for a dataset

#mat<-rbind(cbind(rnorm(n = 20, mean = 1, sd = 4),rnorm(n = 20, mean = 40, sd = 3)),
#cbind(rnorm(n = 20, mean = 7, sd = 5),rnorm(n = 20, mean = 45, sd = 2))

#id<-c(rep(1,20),rep(2,20))
#id[sample(c(1:40))[1:10]]<-0

#cross<-biodecrypt.cross(mat,id)
#plot(mat,type="n")
#biodecrypt.plot(cross)
}
\keyword{record attribution}
