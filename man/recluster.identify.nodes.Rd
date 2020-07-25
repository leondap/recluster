\name{recluster.identify.nodes}
\alias{recluster.identify.nodes}
\title{
Evaluating solutions in multiscale bootstrap
}
\description{
This function helps to understand different behaviours of node supports in multiscale bootstrap by i) plotting trends of support values in different bootstrap scales, ii) identifying the boostrap scale with highest diversification between two groups of nodes and iii) identifying nodes into two classes according to the best bootstrap level identified in (ii) and ploting their mean support values. 
}
\usage{
recluster.identify.nodes(mat, low=TRUE)
}
\arguments{
  \item{mat}{
A matrix containing nodes (rows) and bootstrap levels (columns) as obtained by recluster.multi.
}
  \item{low}{
A logical value indicating if lower scales should be favoured in the selection.
}
  }
\details{
This function recognizes nodes showing different trends of support in multiscale bootstrap. In the analysis of turnover in biogeography some nodes may show a substantial increase in support in a multiscale bootstrap. Areas connected by these nodes may host a few species responsible for turnover, but the biogeographic pattern with respect is clear. Other nodes may show a slow (or no) increase in support. In this case, the links among areas can be considered as uncertain. Partitioning Around Medioids is used to identify two classes of nodes at each level, then the bootstrap scale showing the best diversification in two classes is identified by silhouette scores weighted by differences in mean values between classes. If "low" is set to TRUE the function favours low scales. 
}
\value{
A plot with bootstrap supports and their means (diamonds) for the best combination of two groups of nodes (black and red).
	\item{scale}{
The best bootstrap scale to identify two groups of nodes.
}
	\item{nodes}{
A vector containing classification for nodes in the best bootstrap scale.
}
	}
\references{
Dapporto L., Ramazzotti M., Fattorini S., Talavera G., Vila R., Dennis R. 
"recluster: an unbiased clustering procedure for beta-diversity turnover"
Ecography (2013), 36:1070-1075.
}
\author{
Leonardo Dapporto and Matteo Ramazzotti
}
\examples{
data(multiboot)
recluster.identify.nodes(multiboot)
}
\keyword{ bootstrap }
