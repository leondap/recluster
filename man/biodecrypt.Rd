\name{biodecrypt}
\alias{biodecrypt}
\title{
An algorithm to attribute unidentified occurrence data based on a subset of identified 
records
}
\description{
biodecrypt uses the function ahull from package alphahull to construct concave 
hulls with custom concavity (alpha) for each taxon. This function can also remove sea or 
ground areas from the analysis based on a SpatialPolygonsDataFrame representing the area of interest. 
The main input is represented by: i) a matrix of longitude and 
latitude (decimal degrees of longitude and latitude, WGS84) for all occurrence records, ii) a vector indicating species membership of 
each record in the same order of the matrix (1,2..n for known species and 0 for cases 
to be attributed). By using spatial coordinates, the list of identified records and 
alpha values, biodecrypt computes a concave hull for each species based on known 
records. Then, the function attempts to attribute unknown cases to their most likely 
species based on the comparison of hull localtion, geometry and the location of 
occurrence data (see details). 
}
\usage{
biodecrypt(mat, id, alpha = NULL, ratio = 2.5, buffer = 90, polygon=NULL, checkdist = T, minimum = 7, plot=T, map = NULL, 
    xlim = NULL, ylim = NULL, main = NULL)
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
A distance buffer from hulls (metres).
}
   \item{polygon}{
A SpatialPolygonsDataFrame with area of interests (ground or sea). Typically obtained from Natural Earth (https://www.naturalearthdata.com/). 
If NULL no removal is applied
}
   \item{checkdist}{
Logical, if TRUE cases attributed to a given species based on relative distance from 
hulls but closer to an identified record of another species are not attributed. 
}
   \item{minimum}{
The minimum number of specimens required to build alpha hulls. If the number of 
identified specimens is lower, convex hulls are calculated to improve procedure stability. 
}
\item{map}{
A map to be plotted during the procedure to show the separation progress. 
}
\item{plot}{
Flag to FALSE is plotting the result is not required
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
}
\details{
Once the hulls for species are drawn according to the distribution of known records, 
each unidentified record could be either: i) inside more than one hull, ii) inside a single hull, or iii) outside all hulls. 
- Cases inside more than one hull
 In this case, the function cannot attribute the unidentified records to a species and 
 only the a priori identified records belonging to intersection areas are passed to the final vector as identified. 
- Cases inside a single hull
 The unidentified records falling inside a single hull are attributed to that species 
 if their distance to any other hull is higher than a buffer value provided by the user. Unidentified records inside the buffer of another hull are not attributed.
- Cases outside all hulls
 The unidentified records which do not fall inside any hull are attributed to the 
 closest hull if: i) the distance from the second nearest hull is higher than the 
 buffer and if ii) the ratio between the minimum distance to the second closest hull and to the closest hull is more than a value indicated by the user (ratio).
- Check for distances from the nearest identified record
 As described above, the attribution of unknown records is strictly determined by the 
 distance from the hulls. The biodecrypt function also contains an option (checkdist=T) to check if cases attributed to a given species based on relative distance from hulls are closer to an identified record of another species, which may occasionally occur. If this option is selected (as in the default settings) these cases are not attributed to any species.
In general, alpha values <5 return small and fragmented hulls with low predictive power; moreover, depending on the topology of occurrence data, some low alpha values could return no hulls. In this cases, biodecrypt tries to increase alpha up the the lowest value for which a hull is obtained. The alpha values used is stored in the alphaused vector.
}
\value{
\item{type}{
"sep" an argument to be passed to biodecrypt.plot.
 }
\item{areas}{
The areas of hulls for all the species (in km squares).
 }
\item{intersections}{
The areas of intersections among hulls for each pair of species.
 }
\item{sympatry}{
The fraction of the overlap area compared to the total area of the two hulls.
 }
\item{NUR}{
The percentage of Non-attributed Unidentified Records (NUR).
 }
\item{table}{
The result table with Longitude and Latitude for each occurrence datum, its id after 
the biodecrypt procedure (id2, the result of the procedure) and its initial id (id).
 }
\item{hulls}{
The hulls in sf format.
 }
\item{hullspl}{
The hulls in alphahull format.
 }
\item{alphaused}{
The values of alpha used fopr each species (see details).
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

n1<-7
n2<-7
mat<-rbind(cbind(rnorm(n = n1, mean = 1, sd = 2),rnorm(n = n1, mean = 40, sd = 2)),cbind(rnorm(n = n2, mean = 7, sd = 2),rnorm(n = n2, mean = 45, sd = 2)))

id<-c(rep(1,n1),rep(2,n2))
id[sample(c(1:(n1+n2)))[1:round((n1+n2)/4,0)]]<-0


# Make the separation with custom parameters
attribution<-biodecrypt(mat,id, alpha=c(10,10))

#plot the results
plot(mat,type="n")
biodecrypt.plot(attribution)

#Group plots into pies
biodecrypt.plot(attribution, square=2, minsize=0.5)


# Make the separation with custom parameters
# With a lower fraction values the first hull (alpha equal to 1) can become more 
#concave. Excluded dots works as a punctiform sub-hull in the attribution. 
attribution<-biodecrypt(mat,id, alpha=c(1,5), buffer=20, ratio=2, minimum=5)

#plot the results
plot(mat,type="n")
biodecrypt.plot(attribution)
	
\keyword{record attribution}
