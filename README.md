![](logo.png)

## <b>Ordination Methods for the Analysis of Beta-Diversity Indices</b>


The analysis of different aspects of biodiversity requires specific algorithms. For example, in regionalisation analyses, the high frequency of ties and zero values in dissimilarity matrices produced by Beta-diversity turnover produces hierarchical cluster dendrograms whose topology and bootstrap supports are affected by the order of rows in the original matrix. Moreover, visualisation of biogeographical regionalisation can be facilitated by a combination of hierarchical clustering and multi-dimensional scaling. The recluster package provides robust techniques to visualise and analyse pattern of biodiversity and to improve occurrence data for cryptic taxa.


<b>Sources:</b>

Dapporto L., Ramazzotti M., Fattorini S., Talavera G., Vila R., & Dennis R. L. (2013). recluster: an unbiased clustering procedure for beta‐diversity turnover. Ecography, 36(10), 1070-1075.

Platania L., Menchetti M., Dinca V., Corbella C., Kay-Lavelle I., Vila R., Wiemers M., Schweiger O., Dapporto L. (2020). Assigning occurrence data to cryptic taxa improves climatic niche assessments: biodecrypt, a new tool tested on European butterflies. Global Ecology and Biogeography, DOI:10.1111/geb.13154.

To install recluster use:
```
install.packages("remotes")

remotes::install_github("leondap/recluster")
```

### The zero problem

The dissimilarity indexes of turnover can provide fundamental information in the analysis of beta-diversity patterns but the distance matrices produced on occurrence data show peculiar features:

1) the triangular relationship is rarely met
2) in case of highly nested pattern, the dissimilarity matrix contains many pairs having zero and tied dissimilarity producing problematic cases in clustering

We inspect the data of butterflies of Western Mediterranean islands

```
library(recluster)
data(dataisl)
data(treebut)
```

Calculate the Simpson dissimilatity matrix based on occurrence data

```
simpdiss_b <- recluster.dist(dataisl)
recluster.hist(simpdiss_b)
```

The recluster package includes a function to inspect the fraction of zero and tied values

```
recluster.hist(simpdiss_b)
```
https://github.com/leondap/images/blob/main/cluster%201.png?raw=true(https://github.com/leondap/images/blob/main/histogram.png?raw=true)

So there are 27 tied 0 values and no other tied cells. This potentially creates problems in hierarchical clustering because at the first step the algorithm has to choose between several minimum dissimilarity to create the first link. Actually, when tied values occur, several solutions are possible at each agglomeration step and many different trees may be generated. Statistical packages typically do not cope with this flaw but use arbitrary linking rules (e.g. select the first pair with dissimilarity = 0 in a matrix). Bootstrapping systematically re-samples species but maintains the original order of sites. Consequently, the pairs linked in the first reference tree are more likely to be linked during the entire bootstrap procedure resulting in false strong supports.

The effect of the bias can be easily realised by launching several times a hierarchical clustering by randomly reordering the order of the areas in the matrix which should not affect the topology

```
dataislrnd<-dataisl[sample(1:nrow(dataisl)),]
simpdiss_rnd <- recluster.dist(dataislrnd)
plot(hclust(simpdiss_rnd))
```
Two random examples here show the possible differences among solutions:

![](https://github.com/leondap/images/blob/main/cluster%201.png?raw=true)
![](https://github.com/leondap/images/blob/main/cluster%202.png?raw=true)

Any biogeographical interpretations based on these spurious relationships might be inconsistent

### To solve this problem, recluster.cons function produces a series of trees after randomly re-ordering the row order and creates a consensus tree (0.5 consensus rule by default). By replicating these procedure on subsets of species a bootstrap analysis among the consensus trees can be also computed by recluster.boot.

![](https://github.com/leondap/images/blob/main/flow%20chart.png?raw=true)

```
tree_bf <- recluster.cons (dataisl, p=0.5,method="average")$cons
plot(tree_bf, direction="downwards")
```
![](https://github.com/leondap/images/blob/main/consensus%20tree.png?raw=true)

This tree is resistant to changing the order of the row and shows that the internal division among some groups of areas have no meaning

The function ‘recluster.boot’ allows bootstrapping of nodes in the original consensus tree by applying a user- defined number of consensus trees with user-defined numbers of sampled species (level = 1 means same number of the original matrix i.e species number x 1). 
```
boot_bf <- recluster.boot (tree_bf, dataisl, tr=20, method="average", boot=100, level=1)
recluster.plot (tree_bf, boot_bf)
```
![](https://github.com/leondap/images/blob/main/bootstrap.png?raw=true)

Most nodes received weak support when the number of species randomly sampled with replacement was the same as in the original dataset (level =1). 
In a set of highly nested assemblages, as displayed by most island assemblages, turnover is encompassed by a substantially reduced percentage of species which can univocally link an islands to another area (first case). All bootstrap iterations excluding these species resulted in these islands missing any turnover signal. By applying a multiscale bootstrap, the increase in the number of species randomly selected with repetition can provide greater opportunities for these special taxa to enter the bootstrap matrices, thus increasing the support for these nodes. On the other hand (second case), when a node has a weak (× 1) support because it equally links-up intermediate areas, the increase in the number of species is expected to produce a slower increase in support. 

the ‘recluster.multi’ function to perform multiscale bootstrap analysis. This function requires the same inputs as ‘recluster.boot’ and a number of different scales to be applied as a multiplier for the species sampled at each step. The results are stored in a matrix providing bootstrap values for each node (rows) for each bootstrap scale (columns). 
Try with a muliscale bootstrap with 10 levels starting for x1 to x10 level
```
multiboot_bf <- recluster.multi (tree_bf, dataisl, tr=20, method="average", boot=100, levels=10, step=1)
```
The results can be inspected on the tree by indicating whatever pair of levels (1 and 10 in the example)
```
recluster.plot (tree_bf, multiboot_bf, 1, 10)
```
![](https://github.com/leondap/images/blob/main/bootstrapmulti.png?raw=true)

According to our expectations of two kinds of nodes, some nodes obtained a high values at level = 10, while others did not substantially change. It must be noted that by indefinitely multiplying the number of species, all nodes would attain 100% support at some point. 
Identifying the two kinds of nodes permits recognition as to which links among areas are actually supported by data, even on the basis of a restricted set of species, and which links are uncertain. The ‘recluster.identify.nodes’ function helps in a proper selection of the parameters to ascertain which nodes belong to each class creating two groups of nodes.
```
id_bf<-recluster.identify.nodes(multiboot_bf)
id_bf
```
![](https://github.com/leondap/images/blob/main/identify.png?raw=true)

The function produces a plot where fast growing nodes are marked in black and slow growing nodes in red. The id_bf$scale value also indicates the best scale to be used to identify the two kinds of nodes. Now the first and the third level bootstrap can be plotted on the tree marking with black and red colours strongly and weakly supported nodes

```
recluster.plot (tree_bf, multiboot_bf, 1, 3, id=id_bf$nodes)
```
![](https://github.com/leondap/images/blob/main/multiscale2.png?raw=true)

### Making maps for zooregionalisation

First open the dataset from a previous paper (Dapporto et al 2014)
```
databut <- read.csv("https://raw.githubusercontent.com/leondap/files/main/jbi12315-sup.csv")
```
Extraxct the information about the areas and retain data on butterfly occurrence only
```
latitude<-databut[,6]
longitude<-databut[,7]
richness<-databut[,5]
names<-databut[,3]
databut<-databut[,8:ncol(databut)]
rownames(databut)<-names
head(databut)
```

[UNDER CONSTRUCTION]



References

