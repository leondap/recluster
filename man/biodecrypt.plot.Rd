\name{biodecrypt.plot}
\alias{biodecrypt.plot}
\title{
Plotting biodecrypt and biodecrypt.cross results.
}
\description{
The function plots the results of biodecrypt and biodecrypt.cross analyses. It provides 
plot with circles with different colours to identify different kinds of records. Records known a priori can be distinguished in the plot from records attributed by biodecrypt as likely belonging to a given species or as NUR (or MIR and NIR in biodecrypt.cross).  
}
\usage{
biodecrypt.plot(x,minsize=0.3,pchid=1,cexid=0.1,square=0.001,col=c("red","darkgreen",
"blue","purple"), attributed=c("fade","points"), NUR="black", fading=50, ... )
}
\arguments{
  \item{x}{
An object obtained by biodecrypt or biodecrypt.cross
}
   \item{minsize}{
The size of the dots to be plotted.
}
   \item{pchid}{
The pch of the points marking known cases in case when attributed="points".
}
\item{cexid}{
The size of the points marking known cases in case when attributed="points".
}
\item{square}{
The size of square grid to which occurrence are collapsed and organized in pies. If the 
value is lower than data resolution then records are not grouped in pies.
}
\item{col}{
The colours to be attributed to species: 1...n.
}
\item{attributed}{
The method to plot known records.
}
\item{NUR}{
The colour for NUR records after biodecrypt.
}
\item{fading}{
The degree of fading for the colours of records attributed by biodecrypt if 
attributed="fading" (100 makes the points white).
}
\item{...}{
other parameters of the default plot
}

}
\details{
 The function adds dots to a previous plot (usually a map). The records with a priori 
known attribution (1...n in id) are marked with a point inside the dots (
attributed="points") or by fading the colour of the dots for the records that have been 
attributed by biodecrypt (attributed="fading"). In the results of biodecrypt.cross, MIR 
are represented as black dots and NIR as white dots. For biodecrypt black default colour 
for NUR can be changed. 
}
\value{
a plot
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
	#See examples in biodecrypt and biodecrypt.cross
}
\keyword{record attribution}
