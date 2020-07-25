\name{biodecrypt.optimise}
\alias{biodecrypt.optimise}
\title{
Comparing the values obtained by biodecrypt.wrap, it optimises the combination of alpha, 
buffer and ratio values to be used with biodecrypt function.
}
\description{
The function biodecrypt.optimise analyses the output of biodecrypt.wrap. By default, a 
combination of MIR^2+NIR+NUR is used as a penalty value for the different combinations 
of the parameters (providing a  higher importance to MIR). The exponents can be changed 
by the user. Since the method showing the lowest penalty in cross-validation might not 
necessarily be the optimal value for the final analysis, all the combinations showing a 
penalty value not higher than a certain threshold compared with the analysis showing the 
lowest penalty should be considered as similarly good. We provided a value of 10 as a default, representing a variation of about 3 for each addendum of the penalty. The optimal parameters can then be calculated as mean values of distance ratio, alpha and buffer among those used in these cross-validation analyses, weighted by 1/penalty in order to provide an increasing contribution to the solutions showing the lowest penalty values. 
}
\usage{
biodecrypt.optimise(tab,coef=c(2,1,1), penalty=10)
}
\arguments{
  \item{tab}{
A matrix obtained with biodecrypt.wrap.
}
   \item{coef}{
The three exponents to be applied to MIR, NIR and NUR, respectively, to calculate the 
penalties.
}
   \item{penalty}{
The penalty threshold for inclusion in the calculation.
}
}
\value{
\item{ratio}{
The optimized ratio value.
}
\item{buffer}{
The optimized buffer value.
}
\item{alpha}{
The optimized alpha value.
}
\item{MIR}{
The weighted average MIR among selected combinations.
}
\item{NIR}{
The weighted average MIR among selected combinations.
}
\item{NUR}{
The weighted average MIR among selected combinations.
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
#See the example provided in biodecrypt.wrap
}
\keyword{record attribution}
