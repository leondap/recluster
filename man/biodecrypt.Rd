\name{biodecrypt}
\alias{biodecrypt}
\title{
An algorithm to attribute unidentified occurrence data based on a subset of identified 
records
}
\description{
biodecrypt uses the function getDynamicAlphaHull of the rangeBuilder R package 
(https://github.com/ptitle/rangeBuilder, accessed 2020/02/05) to construct concave 
hulls. This function supports an initial alpha value which is increased until a given 
fraction of samples are included in the resulting hulls (default 95%) and a maximum 
number of partition if satistified (default 10). This function can also remove sea or 
ground areas from the polygons when terrestrial or aquatic organisms, respectively, 
are under study. The main input is represented by: i) a matrix of longitude and 
latitude for all occurrence records, ii) a vector indicating species membership of 
each record in the same order of the matrix (1,2..n for known species and 0 for cases 
to be attributed). By using spatial coordinates, the list of identified records and 
alpha values, biodecrypt computes a concave hull for each species based on known 
records. Then, the function attempts to attribute unknown cases to their most likely 
species based on the comparison of hull localtion, geometry and the location of 
occurrence data (see details). 
}
\usage{
biodecrypt(mat,id,alpha=NULL,ratio=2.5,buffer=90,fraction=0.95, partCount = 10, 
clipToCoast="terrestrial", checkdist=T, proj = "+proj=longlat +datum=WGS84", minimum=7,
map=NULL,xlim=NULL,ylim=NULL,main=NULL)
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
A distance buffer from hulls (in km).
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
identified specimens is lower, convex hulls are calculated to improve procedure stability. 
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
\dontrun{attribution<-biodecrypt(mat,id, clipToCoast="no")}

#plot the results
\dontrun{plot(mat,type="n")}
\dontrun{biodecrypt.plot(attribution)}
#Group plots into pies
\dontrun{biodecrypt.plot(attribution, square=5, minsize=0.8)}


# Make the separation with custom parameters
# With a lower fraction values the first hull (alpha equal to 1) can become more 
#concave. Excluded dots works as a punctiform sub-hull in the attribution. 
\dontrun{attribution<-biodecrypt(mat,id, alpha=c(1,5), buffer=20, ratio=2, fraction=0.80, 
	minimum=3, clipToCoast="no")}

#plot the results
\dontrun{plot(mat,type="n")}
\dontrun{biodecrypt.plot(attribution)}
}
\keyword{record attribution}
