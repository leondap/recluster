\name{biodecrypt.wrap}
\alias{biodecrypt.wrap}
\title{
Wraps the biodecrypt.cross analysis to compare the performance of biodecrypt among different parameters.
}
\description{
The function biodecrypt.wrap wraps the biodecrypt.cross analysis by using all possible combinations of a series of distance ratio, alpha and buffer values to compare their resulting MIR, NIR and NUR.
}
\usage{
biodecrypt.wrap(mat,id,alpha=c(1,5,10,15),alphamat=NULL,ratio=c(2,3,4,5),
buffer=c(0,40,80,120,160),fraction=0.95, partCount=10, checkdist=T, 
clipToCoast="terrestrial", proj="+proj=longlat +datum=WGS84", minimum=7, 
map=NULL,xlim=NULL,ylim=NULL,main=NULL,save=T,name="res_cross.txt",runs=10)
}
\arguments{
  \item{mat}{
A matrix for longitude and latitude (in decimal degrees) for all records.
}
  \item{id}{
A vector indicating species membership of each record (in the same order of mat). Identified records are indicated with 1,2..n, unidentified records with 0.
}
  \item{alpha}{
A vector indicating the initial alpha values. It will be the same for all species
}
  \item{alphamat}{
A matrix indicating different alpha values for different species (optional).
}
  \item{ratio}{
The values of ratio.
}
  \item{buffer}{
The values of buffer.
}
  \item{fraction}{
The minimum fraction of occurrences that must be included in polygon.
}
\item{partCount}{
The maximum number of disjunct polygons that are allowed..
}
\item{checkdist}{
Logical, if TRUE cases attributed to a given species based on relative distance from hulls but closer to an identified record of another species are not attributed. 
}
   \item{clipToCoast}{
Either "no" (no clipping), "terrestrial" (only terrestrial part of range is kept) or "aquatic" (only non-terrestrial part is clipped).
}
      \item{proj}{
the projection information for mat. In this version, the default is the only supported option.
}
   \item{minimum}{
The minimum number of specimens required to build alpha hulls. If the number of identified specimens is lower, convex hulls are calculated to improve procedure stability.  
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
\item{save}{
Logical, if TRUE a result table is saved after each biodecrypt.cross run
}
\item{name}{
The name of the saved file
}
\item{runs}{
The number of runs among which the cases are randomly assigned as non-attributed records
}
}
\details{
 The resulting table can be passed to biodecrypt.optimise to compute the best combination of alpha, buffer and ratio.
}
\value{
\item{table}{
The result table indicating for each cross validation test the MIR, NIR and NUR values together with the used ratio, buffer and alpha values.
}
}
\references{
Platania L. et al.
Assigning occurrence data to cryptic taxa improves climatic niche assessments: biodecrypt, a new tool tested on European butterflies. Glocal Ecology and Biogeography (2020).
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

\dontrun{wrap_data_fast<-biodecrypt.wrap(mat,id, alpha=c(1,4), ratio=2, 
		buffer=20, runs=2, clipToCoast="no")}
\dontrun{parameters<-biodecrypt.optimise(wrap_data_fast$table,penalty=10)}

#Make the example with default 10 runs and more values
\dontrun{wrap_data<-biodecrypt.wrap(mat,id, alpha=c(1,4), ratio=c(2,4), 
		buffer=c(20,50), clipToCoast="no")}
\dontrun{parameters<-biodecrypt.optimise(wrap_data$table)}

#inspect the optimised parameters
\dontrun{parameters}

#Use different alpha for the two species
#alpha for first

\dontrun{alpha1<-c(1,3)}

#alpha for second
\dontrun{alpha2<-c(1,5)}

\dontrun{alphamat<-cbind(alpha1,alpha2)}

\dontrun{wrap_data<-biodecrypt.wrap(mat,id, alphamat=alphamat, ratio=c(2,4),  
                           buffer=c(20,50), clipToCoast="no")}

\dontrun{parameters<-biodecrypt.optimise(wrap_data$table, penalty=20)}

#inspect the optimised parameters

\dontrun{parameters}

}
\keyword{record attribution}
