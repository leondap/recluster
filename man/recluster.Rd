\name{recluster}
\alias{recluster}
\alias{reCluster}
\docType{package}
\title{
Ordination methods for biodiversity patterns.
}
\description{
The analysis of different aspects of biodiversity requires specific algorithms. 
For example, in regionalisation analyses, the high frequency of ties and zero values in 
dissimilarity matrices produced by Beta-diversity turnover produces hierarchical 
cluster dendrograms whose topology and bootstrap supports are affected by the order of 
rows in the original matrix. Moreover, visualisation of biogeographical regionalisation 
can be facilitated by a combination of hierarchical clustering and multi-dimensional 
scaling. The recluster package provides robust techniques to visualise and analyse 
pattern of biodiversity and to improve occurrence data for cryptic taxa.
}
\details{
\tabular{ll}{
Package: \tab recluster\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2020-05-09\cr
License: \tab GPL (>= 2.0)\cr
}
}
\author{
Leonardo Dapporto, Matteo Ramazzotti, Simone Fattorini, Roger Vila, Gerard Talavera, 
Roger H.L. Dennis
Maintainer: Leonardo Dapporto <leondap@gmail.com>
}
\references{
Dapporto L., Ramazzotti M., Fattorini S., Talavera G., Vila R., Dennis R. 
"recluster: an unbiased clustering procedure for beta-diversity turnover"
Ecography (2013), 36:1070-1075.

Dapporto, L., Fattorini, S., Voda, R., Dinca, V., Vila, R. "Biogeography of western 
Mediterranean butterflies: combining turnover and nestedness components of faunal 
dissimilarity." 
J Biogeogr (2014), 41: 1639-1650.

Dapporto L., Voda R., Dinca V., Vila R.
"Comparing population patterns for genetic and morphological markers with uneven sample 
sizes. An example for the butterfly Maniola jurtina"
Methods Ecol Evol (2014), 5, 834-843.

Platania L., Menchetti M. Dinca V., Corbella C., Kay-Lavelle I., Vila R., Wiemers M., 
Schweiger O., Dapporto L. "Assigning occurrence data to cryptic taxa improves climatic 
niche assessments: biodecrypt, a new tool tested on European butterflies". Glocal 
Ecology and Biogeography (2020).

\url{
https://github.com/leondap/recluster
}
}
\keyword{ package }
\examples{
#load model data provided with the package 
\dontrun{
data(datamod)

#explore zero and tied values in the data set
simpdiss<- recluster.dist(datamod)
recluster.hist(simpdiss)

#create and view unbiased consensus tree (100% rule)
constree_full<-recluster.cons(datamod, tr=10, p=1)
plot(constree_full$cons,direction="downwards")

#compute and view node strength
recluster.node.strength(datamod, tr=10)

#create and view unbiased consensus tree (50% rule)
constree_half<-recluster.cons(datamod, tr=10, p=0.5)
plot(constree_half$cons, direction="downwards")

#the latter is the correct tree
tree<-constree_half$cons

#perform and view bootstrap on nodes
boot<-recluster.boot(tree, datamod, tr=10, p=0.5, boot=50)
recluster.plot(tree,boot)

#perform and view multiscale bootstrap on nodes
multiboot<- recluster.multi(tree, datamod, tr=10, boot=50, levels=2, step=1)
recluster.plot(tree,multiboot,low=1,high=2, direction="downwards")

#project and plot a bi-dimensional plot in the RGB colour space
sordiss<- recluster.dist(datamod,dist="sorensen")
points<-cmdscale(sordiss)
col<-recluster.col(points)
recluster.plot.col(col)

#inspect explained diversity for different cuts of a tree
tree<-recluster.cons(datamod, dist="sorensen",tr=10, p=0.5)
expl_div<-recluster.expl.diss(tree$cons,sordiss)
expl_div

#Select cut #4 and group data in RGB space
ncol<-recluster.group.col(col,expl_div$matrix[,4])

#Plot mean values for clusters
recluster.plot.col(ncol$aggr)

#Plot mean colours for sites in the geographic space
lat<-c(2,2,2,1,3,1,1,3,3)
long<-c(1,5,3,3,3,1,5,1,5)
recluster.plot.sites.col(long, lat, ncol$all,text=TRUE)

#Use recluster.region procedure on island butterflies
data(dataisl)
simpson<-recluster.dist(dataisl)
turn_cl<-recluster.region(simpson,tr=10,rettree=TRUE)
turn_cl

#Select solution with three cluster and plot the tree.
plot(turn_cl$tree[[2]])
turn_cl$grouping

#Perform a procrustes with uneven sample size
#Create and plot a target matrix
ex1 <-rbind(c(1,5),c(5,5),c(3,4),c(3,6))
plot(ex1,col=c(1:4),pch=19,xlim=c(0,6),ylim=c(0,6),cex=2)
#Create and plot a matrix to be rotated. Only the points 1-4 are shared
ex2<-rbind(c(3,1),c(3,3),c(2.5,2),c(3.5,2),c(3,4))
plot(ex2,col=c(1:5),pch=19,xlim=c(0,6),ylim=c(0,6),cex=2)

#Perform the procrustes on points 1-4
#Apply the transformation to point 5 of ex2 and plot the matrices
procr1<-recluster.procrustes(ex1,ex2,num=4)
plot(procr1$X,col=c(1:4),pch=19,xlim=c(-4,4),ylim=c(-4,4),cex=2)
plot(procr1$Yrot,col=c(1:5),pch=19,xlim=c(-4,4),ylim=c(-4,4),cex=2)

# Create an example for biodecrypt
mat<-rbind(cbind(rnorm(n = 20, mean = 1, sd = 4),rnorm(n = 20, mean = 40, sd = 3)),
cbind(rnorm(n = 20, mean = 7, sd = 5),rnorm(n = 20, mean = 45, sd = 2)))

id<-c(rep(1,20),rep(2,20))
id[sample(c(1:40))[1:10]]<-0

# Perform biodecrypt with default parameters 
# alpha gets high to include 95% of defaults cases and hulls are almost convex 
attribution<-biodecrypt(mat,id, clipToCoast="no")
#plot the results
plot(mat,type="n")
biodecrypt.plot(attribution)
}
}
