\name{biodecrypt.wrap}
\alias{biodecrypt.wrap}
\title{
Compare biodecrypt parameter combinations by cross-validation
}
\description{
\code{biodecrypt.wrap} runs \code{biodecrypt.cross} over all combinations of ratio, buffer and alpha values, returning MIR, NIR and NUR for each parameter set.
}
\usage{
biodecrypt.wrap(mat, id, alpha = c(1, 5, 10, 15), alphamat = NULL,
                ratio = c(2, 3, 4, 5),
                buffer = c(0, 40, 80, 120, 160),
                checkdist = TRUE, polygon = NULL, minimum = 7,
                map = NULL, xlim = NULL, ylim = NULL, main = NULL,
                save = TRUE, name = "res_cross.txt", runs = 10,
                crs_area = NULL)
}
\arguments{
  \item{mat}{
A two-column matrix or data frame containing longitude and latitude in decimal degrees.
}
  \item{id}{
A vector indicating taxon membership for each record, in the same order as \code{mat}. Unidentified records must be coded as \code{0}. Identified records can be coded as consecutive integers, non-consecutive integers, or character labels.
}
  \item{alpha}{
Vector of initial alpha values to test. Each value is applied to all taxa, unless \code{alphamat} is provided.
}
  \item{alphamat}{
Optional matrix of alpha values. Rows correspond to taxa and columns correspond to alternative alpha combinations to test.
}
  \item{ratio}{
Vector of ratio values to test.
}
  \item{buffer}{
Vector of buffer values to test.
}
  \item{checkdist}{
Logical. If \code{TRUE}, records attributed to a taxon but closer to an identified record of another taxon are left unattributed.
}
  \item{polygon}{
Optional \code{sf} polygon object representing the area of interest. If \code{NULL}, no clipping is applied.
}
  \item{minimum}{
Minimum number of identified records required to attempt alpha hull construction. If fewer records are available, convex hulls are used.
}
  \item{map}{
Optional map object to be added to plots.
}
  \item{xlim}{
Longitude limits for plotting.
}
  \item{ylim}{
Latitude limits for plotting.
}
  \item{main}{
Main title for plots.
}
  \item{save}{
Logical. If \code{TRUE}, the result table is saved after each \code{biodecrypt.cross} run.
}
  \item{name}{
Name of the output file saved when \code{save = TRUE}.
}
  \item{runs}{
Number of cross-validation runs.
}
  \item{crs_area}{
Optional coordinate reference system used for area and intersection calculations. If \code{NULL}, areas are calculated in the original coordinate reference system.
}
}
\details{
The resulting table can be passed to \code{biodecrypt.optimise} to select an optimal combination of alpha, buffer and ratio values.

If \code{alphamat} is provided, it must have one row for each known taxon and one column for each alpha combination to test.
}
\value{
A list with the following components:

\item{table}{
A table reporting, for each parameter combination, MIR, NIR, NUR, ratio, buffer and alpha values.
}
\item{taxa_ids}{
Original taxon identifiers corresponding to the internal numeric taxon codes.
}
\item{crs_area}{
The CRS used for area and intersection calculations.
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

\dontrun{
wrap_data_fast <- biodecrypt.wrap(
  mat, id,
  alpha = c(1, 4),
  ratio = 2,
  buffer = 20,
  runs = 2
)

parameters <- biodecrypt.optimise(wrap_data_fast$table, penalty = 10)
parameters
}

\dontrun{
wrap_data <- biodecrypt.wrap(
  mat, id,
  alpha = c(1, 4),
  ratio = c(2, 4),
  buffer = c(20, 50)
)

parameters <- biodecrypt.optimise(wrap_data$table)
parameters
}

\dontrun{
alpha1 <- c(1, 3)
alpha2 <- c(1, 5)

alphamat <- cbind(alpha1, alpha2)

wrap_data <- biodecrypt.wrap(
  mat, id,
  alphamat = alphamat,
  ratio = c(2, 4),
  buffer = c(20, 50)
)

parameters <- biodecrypt.optimise(wrap_data$table, penalty = 20)
parameters
}

\dontrun{
id_char <- as.character(id)
id_char[id == 1] <- "pippo"
id_char[id == 2] <- "minnie"

wrap_char <- biodecrypt.wrap(
  mat, id_char,
  alpha = c(1, 4),
  ratio = 2,
  buffer = 20,
  runs = 2
)
}
}
\keyword{record attribution}
\keyword{cross-validation}
\keyword{spatial}
