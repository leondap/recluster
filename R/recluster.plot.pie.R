recluster.plot.pie<-function(long, lat, mat=NULL, distance=NULL, loc=NULL, areas=NULL, square=2,map=NULL,add=FALSE,minsize=NULL,proportional=T,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,ylab=NULL,...){
	if(is.null(mat) & is.null(distance)){
 		stop("A distance matrix or a colour matrix from recluster.col must be provided")
	}
	if(is.null(loc)){
		if(is.null(areas)){
			areas<-rep(1,length(lat))
		}
		latsq<-floor(lat/square)*square
		longsq<-floor(long/square)*square
		newcoord<-cbind(aggregate(long~longsq+latsq+areas, FUN="mean"),aggregate(lat~longsq+latsq+areas, FUN="mean")[,4])
		for (k in 1:nrow(newcoord)){
			quali1<-c(1:length(lat))[which(latsq==newcoord[k,2])]
			quali2<-quali1[which(longsq[quali1]==newcoord[k,1])]
			quali3<-quali2[which(areas[quali2]==newcoord[k,3])]
				if(length(quali3)>0){
				loc[quali3]<-k
			}
		}
	}
	if(!is.null(distance)){
		pcoall<-cmdscale(distance)
		mat<-recluster.col(pcoall,st=F,rot=F)
	}
	if(is.null(xlim)){
		xlim<-range(long)
		}
	if(is.null(ylim)){
		ylim<-range(lat)
		}
	xylim<-cbind(xlim,ylim)
	if(!(add)){
		plot(cbind(xylim[1:2],xylim[3:4]),type="n",main=main,xlab=xlab,ylab=ylab)
	}
	if(!is.null(map)){
		plot(map, asp = 2,add=T,xlab=xlab,ylab=ylab)
	}
	if(is.null(minsize)){
		minsize<-min(abs(range(long)[1]-range(long)[2]),abs(range(lat)[1]-range(lat)[2]))/30
	}
	for (i in 1:max(loc)){
		specim<-which(loc==i)
		if(length(specim)==1){
			specimens<-c(long[specim], lat[specim],mat[specim,])
			plotrix::floating.pie(specimens[1],specimens[2],1,radius=minsize,border=NA,col=rgb(specimens[5], specimens[6], specimens[7], maxColorValue = 255))
			plotrix::draw.circle(specimens[1],specimens[2],radius=minsize)
		}
	if(length(specim)>1){
		specimens<-cbind(long[specim], lat[specim],mat[specim,])
		if(length(specim)>3){
			dista<-dist(mat[specim,1:2])
			if(sum(dista)>0){
				mds<-cbind(c(1:length(specim)),cmdscale(dist(mat[specim,1:2]),k=1))
				specimens<-specimens [order(mds[,2]),]
				}
			}
		if(proportional){
			rad<-minsize*(length(specim))^0.25
		}else{
			rad<-minsize
		}
		plotrix::floating.pie(mean(specimens[,1]),mean(specimens[,2]),rep(1,length(specim)),border=NA,radius=rad,col=rgb(specimens[, 5], specimens[, 6], specimens[, 7], maxColorValue = 255))
		plotrix::draw.circle(mean(specimens[,1]),mean(specimens[,2]),radius=rad)
		}
	}
}



