\name{biodecrypt}
\alias{biodecrypt}
\title{
Attribute unidentified occurrence records using spatial hulls built from identified records
}
\description{
\code{biodecrypt} uses occurrence coordinates from identified records to build spatial hulls for each taxon. 
The function then attempts to attribute unidentified records to the most likely taxon based on their position 
relative to these hulls.

Concave hulls are built using \code{ahull} from the \pkg{alphahull} package. If alpha hull construction fails, 
or if the number of identified records is below the chosen minimum, convex hulls are used as a fallback.
}
\usage{
biodecrypt(mat, id, alpha = NULL, ratio = 2.5, buffer = 90000,
           polygon = NULL, checkdist = TRUE, minimum = 7,
           plot = FALSE, map = NULL, xlim = NULL, ylim = NULL,
           main = NULL, crs_area = NULL)
}
\arguments{
  \item{mat}{
A two-column matrix or data frame containing longitude and latitude in decimal degrees. 
Columns are interpreted as longitude and latitude and are assigned the names \code{Long} and \code{Lat}.
}
  \item{id}{
A vector indicating taxon membership for each record, in the same order as \code{mat}. 
Unidentified records must be coded as \code{0}. Identified records can be coded as consecutive integers 
(\code{1, 2, ...}), non-consecutive integers, or character labels.
}
  \item{alpha}{
Initial alpha value for each taxon. If \code{NULL}, a default value is used for all taxa. 
If a single value is provided, it is recycled for all taxa.
}
  \item{ratio}{
Minimum ratio between the distance to the second closest hull and the distance to the closest hull 
required to attribute an unidentified record outside all hulls. Default is \code{2.5}.
}
  \item{buffer}{
Distance buffer from hulls. Records are not attributed if the relevant alternative hull is within this distance.
}
  \item{polygon}{
Optional \code{sf} polygon object representing the area of interest. Hulls are intersected with this polygon 
before attribution and area calculations. If \code{NULL}, no clipping is applied.
}
  \item{checkdist}{
Logical. If \code{TRUE}, records attributed on the basis of hull distances are checked against the nearest 
identified record. If the nearest identified record belongs to another taxon, the record is left unattributed.
}
  \item{minimum}{
Minimum number of identified records required to attempt alpha hull construction. If fewer records are available, 
a convex hull is used.
}
  \item{plot}{
Logical. If \code{TRUE}, hulls and points are plotted during the procedure.
}
  \item{map}{
Optional map object to be added to the plot.
}
  \item{xlim}{
Longitude limits for plotting.
}
  \item{ylim}{
Latitude limits for plotting.
}
  \item{main}{
Main title for the plot.
}
  \item{crs_area}{
Optional coordinate reference system used for area and intersection calculations. 
If \code{NULL}, areas are calculated in the original coordinate reference system. 
For large areas or accurate comparisons of hull overlap, users should provide a suitable projected 
equal-area CRS, for example \code{3035} for Europe.
}
}
\details{
The input coordinates are assumed to be longitude and latitude in decimal degrees, WGS84 
(\code{EPSG:4326}). For each taxon represented by identified records, \code{biodecrypt} builds a spatial hull. 
Alpha hulls are attempted first. If the alpha hull cannot be converted to a valid polygon, the function 
tries progressively larger alpha values. If no valid alpha hull is obtained, a convex hull is used.

Each unidentified record can fall into one of three cases:

\enumerate{
  \item Inside more than one hull. In this case, the record is not attributed.
  \item Inside a single hull. The record is attributed to that taxon if its distance to all other hulls is 
  greater than \code{buffer}.
  \item Outside all hulls. The record is attributed to the nearest hull if the distance to the second nearest 
  hull is greater than \code{buffer} and if the ratio between the second nearest and nearest hull distances 
  is greater than \code{ratio}.
}

If \code{checkdist = TRUE}, records attributed to a taxon are further checked against the nearest identified 
record. If the nearest identified record belongs to another taxon, attribution is removed.

The argument \code{crs_area} affects only area and intersection calculations. It does not change the input 
coordinates or the coordinate reference system used for point-in-polygon and distance operations. When 
\code{crs_area = NULL}, area calculations are performed in the original CRS. This may be sufficient for small 
or exploratory analyses, but for large geographic extents a suitable projected CRS is recommended.
}
\value{
A list with the following components:

\item{type}{
The string \code{"sep"}, used by \code{biodecrypt.plot}.
}
\item{areas}{
Areas of the hulls for all taxa. Units depend on the CRS used for area calculation. If \code{crs_area} is a 
metric CRS, areas are in square metres.
}
\item{intersections}{
Matrix of pairwise intersection areas among hulls.
}
\item{sympatry}{
Matrix of pairwise overlap values, calculated as intersection area divided by union area.
}
\item{NUR}{
Percentage of Non-attributed Unidentified Records.
}
\item{table}{
Table with longitude, latitude, attributed id after the procedure (\code{id2}), and original id (\code{id}).
}
\item{table_internal}{
Internal numeric version of the result table. This is used for plotting and for handling non-consecutive 
or character taxon labels.
}
\item{hulls}{
Hulls in \code{sf} format.
}
\item{hullpl}{
Original hull objects used during hull construction.
}
\item{alphaused}{
Alpha values effectively used for each taxon. \code{NA} indicates that the convex hull fallback was used.
}
\item{taxa_ids}{
Original taxon identifiers corresponding to the internal numeric taxon codes.
}
\item{crs_area}{
The CRS used for area and intersection calculations.
}
}
\references{
Platania L. et al. Assigning occurrence data to cryptic taxa improves climatic niche assessments: 
biodecrypt, a new tool tested on European butterflies. Global Ecology and Biogeography (2020).
}
\author{
Leonardo Dapporto
}
\examples{
n1 <- 7
n2 <- 7

mat <- rbind(
  cbind(rnorm(n = n1, mean = 1, sd = 2),
        rnorm(n = n1, mean = 40, sd = 2)),
  cbind(rnorm(n = n2, mean = 7, sd = 2),
        rnorm(n = n2, mean = 45, sd = 2))
)

id <- c(rep(1, n1), rep(2, n2))
id[sample(seq_len(n1 + n2))[1:round((n1 + n2) / 4, 0)]] <- 0

attribution <- biodecrypt(mat, id, alpha = c(10, 10))

plot(mat, type = "n")
biodecrypt.plot(attribution)

biodecrypt.plot(attribution, square = 2, minsize = 0.5)

attribution <- biodecrypt(mat, id, alpha = c(1, 5), buffer = 20, ratio = 2, minimum = 5)

plot(mat, type = "n")
biodecrypt.plot(attribution)

id_char <- as.character(id)
id_char[id == 1] <- "pippo"
id_char[id == 2] <- "minnie"

attribution_char <- biodecrypt(mat, id_char, alpha = c(10, 10))

plot(mat, type = "n")
biodecrypt.plot(attribution_char)
}
\keyword{record attribution}
\keyword{spatial}
\keyword{alpha hull}
