biodecrypt.plot<-function(x,minsize=0.3,pchid=1,cexid=0.1,square=0.001,col=c("red","darkgreen","blue","purple"), attributed=NULL, hull=T, NUR="black", fading=50, ...){
if(x$type=="sep"){
		data<-as.data.frame(x$table)
		newcol<-as.character(paste(data[,1],data[,2],sep="_"))
		data<-cbind(data,newcol)
		quanti<-aggregate(rep(1,nrow(data))~ data[,5],FUN="sum")
		quanti<-quanti[which(quanti[,2]>1),]
		zeri<-which(data[,5]%in%quanti[,1] & data[,4]==0)
		if(length(zeri)>0){
			data2<-data[-zeri,]}else{
			data2<-data}
		matcol<-cbind(matrix(1,nrow(data2),2),matrix(NA, nrow(data2),3))
		for(sp in 1:max(as.numeric(data2[,3]))){
			#sp<-2
			which<-which(as.numeric(data2[,3])==sp)
			matcol[which,3]<-as.vector(col2rgb(col[sp]))[1]
			matcol[which,4]<-as.vector(col2rgb(col[sp]))[2]
			matcol[which,5]<-as.vector(col2rgb(col[sp]))[3]
		}
		if(is.null(attributed)){
			attributed<-"fade"
		}
		if(attributed=="fade"){
			lower<-which(data2[,4]==0 & data2[,3]>0)
			for (c in 1:length(lower)){
				for(colo in 3:5){
					matcol[lower[c],colo]<-matcol[lower[c],colo]+((255-matcol[lower[c],colo])/(100/fading))
				}
			}
		}
		if(hull){
			for(sphull in 1: length(x$hullpl)){
				plot(x$hullpl[[sphull]],col=col[sphull], add=T)
			}
		
		}
		matcol[which(is.na(matcol[,3])),3:5]<-col2rgb(NUR)
		recluster.plot.pie(data2[,1],data2[,2],mat=as.matrix(matcol),square=square,minsize=minsize,proportional=F,add=T)
		if(attributed=="points"){
			lat<-data2[,2]
			long<-data2[,1]
			latsqo<-floor(lat/square)*square
			longsqo<-floor(long/square)*square
			newcoord0<-cbind(aggregate(long~longsqo+latsqo, FUN="mean"),aggregate(lat~longsqo+latsqo, FUN="mean")[,3]) 
			newcoord0[,5]<-paste(newcoord0[,1],newcoord0[,2])
			lat<-data2[which(data2[,4]>0),2]
			long<-data2[which(data2[,4]>0),1]
			latsqo<-floor(lat/square)*square
			longsqo<-floor(long/square)*square
			newcoord1<-cbind(aggregate(long~longsqo+latsqo, FUN="mean"),aggregate(lat~longsqo+latsqo, FUN="mean")[,3])
			newcoord1[,5]<-paste(newcoord1[,1],newcoord1[,2])
			newcoord3<-newcoord0[which(newcoord0[,5] %in% newcoord1[,5]),]
			points(newcoord3[,3:4],cex=cexid,pch=pchid)
		}
	}
	if(x$type=="cross"){
		data2<-x$cross
		matcol<-cbind(data2[,c(5,6)],matrix(NA,nrow(data2),3))
		for(sp in 1:max(data2[,1])){
			#sp<-2
			which<-which(data2[,2]==sp)
			matcol[which,3]<-as.vector(col2rgb(col[sp]))[1]
			matcol[which,4]<-as.vector(col2rgb(col[sp]))[2]
			matcol[which,5]<-as.vector(col2rgb(col[sp]))[3]
		}
		matcol[which(data2[,3]==1),3:5]<-col2rgb("black")
		matcol[which(data2[,4]==1),3:5]<-col2rgb("white")
		matcol<-matcol[complete.cases(matcol),]
		recluster.plot.pie(matcol[,1],matcol[,2],mat=as.matrix(matcol),square=square,minsize=minsize,proportional=F,add=T)	
	}
}
