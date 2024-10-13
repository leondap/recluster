![](logo.png)

## <b>Ordination Methods for the Analysis of Species Diversity </b>


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

```

The recluster package includes a function to inspect the fraction of zero and tied values

```
recluster.hist(simpdiss_b)
```
![](https://github.com/leondap/images/blob/main/histogram.png?raw=true)

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

make the recluster.cons tree
```
tree_simp<-recluster.cons(databut,dist="simpson",tr=100)
plot(tree_simp$cons, direction="downwards",cex=0.5)
```
![](https://github.com/leondap/images/blob/main/cluster_but.png?raw=true)

The presence of several polytomies indicates the need for the recluster.cons producere instead of a typical hierarchical clustering

Now we have to decide which cut to apply to identify the regions. According to Holt et al (2013) recognition of regions can be facilitated by identifying a cut in the tree explaining at least 90% of the dissimilarity. To do this it is enough to compute the ratio between the sum of cell values representing distances among sites belonging to different clusters and the sum of all the values included in the dissimilarity matrix. The function recluster.expl.diss computes this value for all possible cuts of a dendrogram.

```
dist_simp<-recluster.dist(databut)
expl_diss<-recluster.expl.diss(tree_simp$cons, dist_simp)
expl_diss
```
which returns a large matrix expl_diss$matrix representing cluster attribution of area for each cut and two vectors $expl.div and $nclust including the explained dissimilarity at each cut and teh number of resulting clusters, respectively:
```

$expl.div
 [1] 0.5504787 0.7565929 0.7917113 0.8400110 0.8659223 0.8674571 0.9362333 0.9392667 0.9435704 0.9497930 0.9510247
[12] 0.9519164 0.9521848 0.9524516 0.9887708 0.9938220 0.9942584 0.9944479 0.9945420 0.9950057 0.9951586 0.9952174
[23] 0.9967261 0.9967790 0.9984160 0.9984552 0.9984908 0.9985263 0.9997911 0.9999752 1.0000000

$nclust
 [1]  2  3  4  5  6  7  8 10 11 12 13 15 16 17 29 35 36 38 39 42 43 44 50 51 58 59 60 61 70 72 73


```
Explaining that the seventh cut (creating 8 clusters) explains more than 93% of dissimilairity and that successive cuts only slightly increase this value

At this stage the hierarchical cluster procedure can be paired with Principal Coordinates Analysis.  We used the cmdscale function. 

```
pcoa_simp<-cmdscale(dist_simp)
```

According to Kreft & Jetz (2010), the two-dimensional configuration from PCoA can be plotted into a two-dimensional space, where the four corners are represented by the colours pure red, yellow, green and blue. The position of each area in this space can be represented by a combination of RGB colours. The function recluster.col projects a two-dimensional configuration into a new RGB space. The flag st=F permits to maintain PCoA axes in the same unit of the original dissimilarity matrix (e.g. turnover dissimilarity)
```
colours_simp<-recluster.col(pcoa_simp, st=F)
recluster.plot.col(colours_simp, cex=1.5, cext=0.6)
```
![](https://github.com/leondap/images/blob/main/recluster.plot.col.png?raw=true)

Kreft & Jetz (2010) and Holt et al. (2013) aggregated sites in the colour space according to the results obtained by cluster analysis by computing mean coordinates for sites belonging to the same group. This is made with the recluster.group.col function, requiring a table as obtained by the recluster.col function and a vector describing group membership for each case. This vector can be obtained directly from the table provided by recluster.expl.diss$matrix by selecting the column corresponding to the selected cut. In our case the first cut producing >90% of explained dissimilarity is #4. So the group membership is obtain as:
```
membership<-expl_diss$matrix[,7]
```
and the colours for the resulting clusters as:
```
new_colours_sim<-recluster.group.col (colours_simp,membership)
recluster.plot.col(new_colours_sim$aggr,text=F,cex=2)
```
![](https://github.com/leondap/images/blob/main/recluster.plot.aggr.png?raw=true)

Now the areas with their colours can be plotted in a map by using recluster.plot.pie function together with a map. The very low square value avoid that areas are merged in pies.
```
library(rworldmap)
library(rworldxtra)
map <- getMap(resolution = "high")
recluster.plot.pie(longitude,latitude,mat=new_colours_sim$all,square=0.001,minsize=0.3,xlab="Longitude", ylab="Latitude")
plot(map, add=T)
```
![](https://github.com/leondap/images/blob/main/map%20final%203.png?raw=true)

Where the relationships among areas appear very clear in a sight.

### Zooregionalisation at mid-small scale using geographic cells
When spatial units are close to one another and not separated by conspicuous barriers (like sea straits or mountain chains) the existence of many tied values, especially zero value, affects the possibility to create  dendrograms even by the recluster.cons algorithm because equivalent solutions are significantly different from one another and in practice makes impossible to obtain any reliable 50% consensus among trees

https://methodsblog.com/2016/08/11/biogeographic-regions/

The recluster.region algorithm overcomes the row-order bias and the restrictive 50% consensus rule.
The example is based on a dataset of British butterflies provided in the original paper describing recluster.region (Dapporto et al. 2015)
```
data <- read.csv("https://raw.githubusercontent.com/leondap/files/main/mee312415-sup-0005-2001-9-selected.csv")
```
Extract coordinates and occurrence data
```
coordin<-data[,3:4]
table_09<-data[,5:ncol(data)]
```
Compute the Simpson dissimilarity index and inspect its values showing a high incidence of zero values
```
simpson<-recluster.dist(table_09)
recluster.hist(simpson)
```
![](https://github.com/leondap/images/blob/main/recluster.hist.region.png?raw=true)

Compute the PCoA on dissimilarities and project the configuration in RGB space
```
pcoa<-cmdscale(simpson)
rgbcol<-recluster.col(pcoa)
recluster.plot.col(rgbcol)
```
![](https://github.com/leondap/images/blob/main/RGBregion.png?raw=true)

Perform the recluster.region analysis using the Ward method (all methods implemented in hclust plus "pam" and "diana" are allowed) with a number of cluster solutions ranging from 2 (minclust=2) and 4 (maxcl=4), based on 100 random trees (takes some minutes)
```
solution<-recluster.region(simpson,method="ward.D",tr=100,mincl=2,maxcl=4)
```
Select the number of clusters to plot (2 to 4)
```
clusters<-4
```
Compute cluster colours and plot the solution
```
grp<-recluster.group.col(rgbcol,solution$grouping[,clusters-1])
plot(coordin[,1], coordin[,2],  col = rgb(grp$all[, 3], grp$all[,4], grp$all[, 5], maxColorValue = 255), cex = 0.7, pch=15)
```
![](https://github.com/leondap/images/blob/main/region4x.png?raw=true)

To enhance contrasts among colours to a maximum variability by keeping the relative poisition of cluster it is enough to compute a new recluster.col
```
newcol<-recluster.col(grp$all[,1:2])
plot(coordin[,1], coordin[,2],  col = rgb(newcol[, 3], newcol[,4], newcol[, 5], maxColorValue = 255), cex = 0.7, pch=15)
```
![](https://github.com/leondap/images/blob/main/region4xcolours.png?raw=true)

### The biodecrypt functions
Occurrence data are fundamental to macroecology, but accuracy is often compromised when multiple units are lumped together (e.g., in recently separated cryptic species, in genetic lineages or in citizen science records). When such units are at least in part allopatric unidentified occurrences can be objectively attributed to the most probable unit based on a subset of identified records. The objective of the algorithm is to reliably attribute species membership to a set of ambiguous records belonging to two (or more) cryptic entities based on the distribution of a subset of accurately determined records. The main idea is that records from an area where only one taxon occurs can be attributed with confidence, while records from the areas of sympatry or too far from any ascertained record cannot be reliably attributed.
The main inputs for the functions are a matrix with longitude and latitude (decimal longitude and latitude, WGS84) for all the occurrence data and a vector (in the same order) providing their identification. The  identified records must be indicated in the vector with a sequential numeric value (1, 2, …, n), which represents the verified membership to the nth unit. The occurrence data with unknown identification (unidentified records) are marked with a 0. Based on this vector and on the geographic coordinates of identified records, biodecrypt builds concave hulls of distribution for each species (ahull algorithm).
After the construction of the alpha-hulls, biodecrypt attempts the attribution of unidentified records to the most likely unit. For this aim, biodecrypt also requires a buffer and a ratio value (explained below). Based on hull geometry and their relative position, each unidentified record could be either: (a) inside more than one hull, (b) inside a single hull, or (c) outside all hulls. The three cases are treated separately and referred to the figure below which represents an example of the assignment procedure based on the overlapping distributions of Polyommatus icarus (red) and Polyommatus celina (blue) in Iberia. Red and blue dots represent the sites from where specimens of the respective species have been sequenced, empty circles represent sites with unidentified records. The continuous red and blue lines represent the hulls obtained for the two species based on occurrence of sequenced specimens, the dotted blue line represents the buffer of the P. celina hull (for clarity, the buffer of the P. icarus hull is not represented). 

![](https://github.com/leondap/images/blob/main/biodecrypt_fig1.jpg?raw=true)

1 Cases inside more than one hull
In this case, the function cannot attribute the unidentified records to a species (record 1 in the Figure) and only the a priori identified records belonging to intersection areas are passed to the final vector as identified.

2 Cases inside a single hull
The unidentified records falling inside a single hull are attributed to that species if their distance to any other hull is higher than the buffer value (in km) provided by the user (record 2 in the Figure). Unidentified records inside the buffer of another hull are not attributed (record 3 in the Figure).

3 Cases outside all hulls
The unidentified records that do not fall inside any hull are attributed to the closest hull if: (a) the distance from the second nearest hull is higher than the buffer and if (b) the ratio between the minimum distance to the second closest hull and to the closest hull is more than the ratio value indicated by the user. For example, in Figure 2 record 4 is not attributed while record 5 is attributed to Polyommatus celina.

4 Check for distances from the nearest identified record
As described above, the attribution of unattributed records is strictly determined by the distance from the hulls. The biodecrypt function also contains an option (“checkdist”) to check if records attributed to a given species based on relative distance from hulls are closer to an identified record of another species, which may occasionally occur. If this option is selected (default) these records are not attributed to any species (record 6 in the Figure).

Biodecrypt also computes the area of overlap among hulls and can exclude sea areas for terrestrial organisms and land area for marine ones based on a polygon defining them.
Open the libraries, the map and the land polygon from natural earth
```
library(recluster)
library(rworldmap)
library(rworldxtra)
map <- getMap(resolution = "low")
library(rnaturalearth)
polygon <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf")
```
Open the data for the Polyommatus icarus and celina taxa and obtain the matrix of coordinates (mat) and the id vector (id)
```
data <- read.csv("https://raw.githubusercontent.com/leondap/files/refs/heads/main/Polyommatus.csv")
mat<-data[,c(7:6)]
id<-data[,3]
```
Run a biodecrypt analysis with alpha=5 for both species and buffer=50000 metres
```
biodecrypt1<-biodecrypt(mat, id,alpha=c(5,5),map=map, buffer=50000, polygon=polygon)
plot(map, xlim=range(mat[,1]),ylim=range(mat[,2]))
biodecrypt.plot(biodecrypt1, col=c("red","blue"))

```
![](https://github.com/leondap/images/blob/main/Biodecrypr_fig2.jpg?raw=true)



###References


Dapporto, L., Ciolli, G., Dennis, R. L., Fox, R., & Shreeve, T. G. (2015). A new procedure for extrapolating turnover regionalization at mid‐small spatial scales, tested on British butterflies. Methods in Ecology and Evolution, 6(11), 1287-1297.

Holt, B.G., Lessard, J.-P., Borregaard, M.K., Fritz S.A., Araújo, M.B., Dimitrov, D., Fabre, P.-H. Graham, C.H., Graves, G.R., Jønsson, K.A., Nogués-Bravo, D., Wang Z., Whittaker, R.J., Fjeldså, J. & Rahbek, C. (2013) An update of Wallace’s zoogeographic regions of the world Science, 339, 74–78.

Kreft, H. & Jetz, W. (2010) A framework for delineating biogeographic regions based on species distributions. Journal of Biogeography, 37, 2029–2053.

