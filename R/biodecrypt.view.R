biodecrypt.view<-function(mat,id,alpha=NULL, fraction=0.95, partCount=10, map=NULL,main=NULL, xlim=NULL,ylim=NULL,clipToCoast="terrestrial",cex0=0.2,cex1=0.25){
	res<-NULL
	taxa<-max(id)
	if(is.null(alpha)){
		alpha=rep(8,taxa)
	}
	if(is.null(xlim)){
		xlim<-c(min(mat[,1]),max(mat[,1]))
	}
	if(is.null(ylim)){
		ylim<-c(min(mat[,2]),max(mat[,2]))
	}
	plot(cbind(xlim,ylim),type="n")
	if(!is.null(map)){
		plot(map,add=T)
	}
	points(mat,col="grey",cex=cex0)
	points(mat,col=id,cex=cex1)
	hulls<-list()
	areas<-NULL
	oldw <- getOption("warn")
		options(warn = -1)
		for (sp in 1:taxa){
			taxsp<-which(id==sp)
			hulla<-mat[taxsp,]
			hullas<-hulla[!duplicated(hulla), ]
			hull<-rangeBuilder::getDynamicAlphaHull(hullas[,c(1,2)],fraction = fraction, partCount = partCount, buff = 0, clipToCoast=clipToCoast,initialAlpha=alpha[sp])[[1]]
			plot(hull,border=sp,add=T)
			hulls[sp]<-hull 
			areas[sp]<-area(hulls[[sp]])/1000000
		}
		intersect<-matrix(NA,taxa,taxa)
		sympatry<-intersect
		for(k in 1:(taxa-1)){
			for(c in 2:taxa){
				inter<-raster::intersect(hulls[[k]],hulls[[c]])
				if(!is.null(inter)){
					intersect[k,c]<-(raster::area(inter)/1000000)
					sympatry[k,c]<-(raster::area(inter)/1000000)/(areas[c]+areas[k]-area(inter)/1000000)
					intersect[c,k]<-intersect[k,c]
					sympatry[c,k]<-sympatry[k,c]
				}
				if(is.null(inter)){
					intersect[k,c]<-0
					sympatry[k,c]<-0
					intersect[c,k]<-intersect[k,c]
					sympatry[c,k]<-sympatry[k,c]
				}		
			}
		}	
		options(warn = oldw)
		res$areas<-areas
		res$intersections<-intersect
		res$sympatry<-sympatry
		return(res)
}