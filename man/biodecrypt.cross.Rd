\name{biodecrypt.cross}
\alias{biodecrypt.cross}
\title{
Cross-validation of biodecrypt attribution
}
\description{
\code{biodecrypt.cross} wraps \code{biodecrypt} to perform cross-validation on known records. In each run, a subset of records with known taxon membership is temporarily treated as unidentified and then attributed using \code{biodecrypt}. Predicted and original memberships are compared to estimate misidentification and non-identification rates.
}
\usage{
biodecrypt.cross(mat, id, alpha = NULL, ratio = 2.5, buffer = 90,
                 checkdist = TRUE, minimum = 7, polygon = NULL,
                 map = NULL, xlim = NULL, ylim = NULL, main = NULL,
                 runs = 10, test = TRUE, crs_area = NULL)
}
\arguments{
  \item{mat}{
A two-column matrix or data frame containing longitude and latitude in decimal degrees.
}
  \item{id}{
A vector indicating taxon membership for each record, in the same order as \code{mat}. Unidentified records must be coded as \code{0}. Identified records can be coded as consecutive integers, non-consecutive integers, or character labels.
}
  \item{alpha}{
Initial alpha value for each taxon. If \code{NULL}, a default value of \code{8} is used for all taxa. If a single value is provided, it is recycled for all taxa.
}
  \item{ratio}{
Minimum ratio between the distance to the second closest hull and the distance to the closest hull required for attribution. Default is \code{2.5}.
}
  \item{buffer}{
Distance buffer from hulls.
}
  \item{checkdist}{
Logical. If \code{TRUE}, records attributed to a taxon but closer to an identified record of another taxon are left unattributed.
}
  \item{minimum}{
Minimum number of identified records required to attempt alpha hull construction. If fewer records are available, convex hulls are used.
}
  \item{polygon}{
Optional \code{sf} polygon object representing the area of interest. Hulls are intersected with this polygon before attribution and area calculations. If \code{NULL}, no clipping is applied.
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
  \item{runs}{
Number of cross-validation runs. A value of \code{10} performs a ten-fold cross-validation. If \code{runs} equals the number of known records, the procedure approximates a jackknife.
}
  \item{test}{
Logical. If \code{TRUE}, a full \code{biodecrypt} analysis is also performed to compute NUR and hull overlap summaries.
}
  \item{crs_area}{
Optional coordinate reference system used for area and intersection calculations. If \code{NULL}, areas are calculated in the original coordinate reference system. For large areas or accurate overlap estimates, users should provide a suitable projected equal-area CRS.
}
}
\details{
The function first keeps only records with known membership. For each taxon, known records are randomly ordered and split across \code{runs}. At each run, one subset per taxon is temporarily coded as unidentified (\code{0}), and \code{biodecrypt} is used to reassign those records.

The resulting predicted membership is compared with the original membership. Two summary statistics are returned: Mis-Identified Records (MIR), the percentage of known records attributed to the wrong taxon, and Non-Identified Records (NIR), the percentage of known records left unattributed.

If \code{test = TRUE}, a standard \code{biodecrypt} run is also performed on the full input dataset, allowing the calculation of Non-attributed Unidentified Records (NUR), hull areas, intersections and sympatry.
}
\value{
A list with the following components:

\item{type}{
The string \code{"cross"}, used by \code{biodecrypt.plot}.
}
\item{NUR}{
Percentage of Non-attributed Unidentified Records, returned if \code{test = TRUE}.
}
\item{areas}{
Hull areas for all taxa, returned if \code{test = TRUE}. Units depend on the CRS used for area calculation.
}
\item{intersections}{
Matrix of pairwise intersection areas among hulls, returned if \code{test = TRUE}.
}
\item{sympatry}{
Matrix of pairwise overlap values, returned if \code{test = TRUE}.
}
\item{table}{
Result table from the full \code{biodecrypt} analysis, returned if \code{test = TRUE}.
}
\item{cross}{
Cross-validation table with original membership, predicted membership, MIR/NIR flags, longitude and latitude.
}
\item{MIR}{
Percentage of Mis-Identified Records.
}
\item{NIR}{
Percentage of Non-Identified Records.
}
\item{taxa_ids}{
Original taxon identifiers corresponding to the internal numeric taxon codes.
}
}
\references{
Platania L. et al. Assigning occurrence data to cryptic taxa improves climatic niche assessments: biodecrypt, a new tool tested on European butterflies. Global Ecology and Biogeography (2020).
}
\author{
Leonardo Dapporto
}
\examples{
mat <- rbind(
  cbind(rnorm(n = 20, mean = 1, sd = 4),
        rnorm(n = 20, mean = 40, sd = 3)),
  cbind(rnorm(n = 20, mean = 7, sd = 5),
        rnorm(n = 20, mean = 45, sd = 2))
)

id <- c(rep(1, 20), rep(2, 20))
id[sample(seq_len(40))[1:10]] <- 0

cross <- biodecrypt.cross(mat, id)

plot(mat, type = "n")
biodecrypt.plot(cross)

id_char <- as.character(id)
id_char[id == 1] <- "pippo"
id_char[id == 2] <- "minnie"

cross_char <- biodecrypt.cross(mat, id_char)

plot(mat, type = "n")
biodecrypt.plot(cross_char)
}
\keyword{record attribution}
\keyword{cross-validation}
\keyword{spatial}
