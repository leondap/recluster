\name{recluster.plot.pie}
\alias{recluster.plot.pie}
\title{
Plotting pies with RGB colours on a custom coordinate space
}
\description{
This function groups cases based on a space grid in a user defined set of coordinates 
(usually longitude and latitude) and plot them in pies using RGB colours. The function 
can either use an output from recluster.col function or compute colours based on any 
distance matrix where the cases are in the same order as in the latitude and longitude 
data.
}
\usage{
recluster.plot.pie(long, lat, mat=NULL, distance=NULL, loc=NULL, areas=NULL, square=2,
                   map=NULL,add=FALSE,minsize=NULL,proportional=T,xlim=NULL,ylim=NULL,
                   main=NULL,xlab=NULL,ylab=NULL,...)
}
\arguments{
	\item{long}{
A vector indicating longitude for cases.
}
	\item{lat}{
A vector indicating latitude for cases.
}
  \item{mat}{
A matrix inherited by recluster.col.
}
  \item{distance}{
A dissimilarity matrix for cases.
}
\item{loc}{
A list of localities to group cases, if available.
}
\item{square}{
The grid to be used to divide cases into groups (2 degrees latitude and longitude by 
default).
}
\item{areas}{
An additional vector to divide groups (e.g. islands versus continents).
}
\item{map}{
A map to be plotted.
}
 \item{add}{
A logical. If TRUE then the points are added to an existing graph.
}
 \item{minsize}{
Dimension for the dimension of a single-case pie.
}
 \item{proportional}{
A logical. If TRUE then the point area is proportional to the number of cases.
}
 \item{xlim}{
Limits of the plot in the x-axis.
}
 \item{ylim}{
Limits of the plot in the y-axis.
}
\item{main}{
The title of the graph.
}
 \item{xlab}{
The label of x-axis
}
 \item{ylab}{
The label of y-axis
}
 \item{...}{
See par() for other graphical parameters
  }
}
\value{
A colour plot. 
}
\references{
Hernandez Roldan J.L., Dapporto L., Dinca V, Vicente J.C., Hornett E.A., Sichova J., 
Lukhtanov V.L., Talavera G. & Vila, R. Integrative analyses unveil speciation linked to host plant shift in Spialia butterflies. Molecular Ecology (2016) 25: 4267-4284.
}
\author{
Leonardo Dapporto
}
\examples{
# create a virtual dataset and a corresponding distance matrix
lat<-runif(50,min=20,max=40)
long<-runif(50,min=20,max=40)
datavirtual<-data.frame(replicate(20,sample(0:1,50,rep=TRUE)))
dist<-recluster.dist(datavirtual)

# Make a plot using a custom distance
recluster.plot.pie(long,lat,distance=dist,xlab="Longitude",ylab="Latitude")

# Make a plot using a recluster.col matrix
colours<-recluster.col(cmdscale(dist))
recluster.plot.pie(long,lat,mat=colours,xlab="Longitude",ylab="Latitude")

# Make points of equal size
recluster.plot.pie(long,lat,mat=colours,xlab="Longitude", proportional=FALSE,
ylab="Latitude")

# Reduce the grid
recluster.plot.pie(long,lat,distance=dist,square=1, xlab="Longitude",ylab="Latitude")

# Reduce the size of the plots
recluster.plot.pie(long,lat,distance=dist,xlab="Longitude",ylab="Latitude", minsize=0.5)
# Use a custom colour matrix
pcoa<-cmdscale(dist)
colour<-recluster.col(pcoa)
recluster.plot.col(colour)
recluster.plot.pie(long,lat,mat=colour,xlab="Longitude",ylab="Latitude")

# Include an additional factor for separating dots in groups(e.g. two continents)
continent<-rep(1,50)
continent[which(long>25)]<-2
recluster.plot.pie(long,lat,distance=dist,xlab="Longitude",ylab="Latitude", 
areas=continent)
}
\keyword{ plot }
