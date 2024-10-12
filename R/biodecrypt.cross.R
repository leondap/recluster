
biodecrypt.cross<-function(mat,id,alpha=NULL,ratio=2.5,buffer=90, checkdist=T, minimum=7, polygon=NULL, map=NULL,xlim=NULL,ylim=NULL,main=NULL,runs=10,test=T){
	res<-NULL
	res$type<-"cross"
	taxa<-max(id)
	idr<-id[which(id>0)]
	matrix<-mat[which(id>0),]
	if(is.null(alpha)){
		alpha=rep(8,taxa)
	}
	q<-aggregate(rep(1,length(idr))~idr,FUN="sum")
	min<-which(q[,2]==min(q[,2]))
	if(min(q[,2])<4){
		stop("The cross validation procedure requires a minimum of 4 distinct points")
	}
	if(test){
	test_run<-biodecrypt(mat,id,alpha=alpha,ratio=ratio,buffer=buffer,map=map,polygon=polygon, checkdist=checkdist, xlim=xlim,ylim=ylim,main=main)
	res$NUR<-test_run$NUR
	res$areas<-test_run$areas
	res$intersections<-test_run$intersections
	res$sympatry<-test_run$sympatry
	res$table<-test_run$table
	}
	species<-unique(idr)
	which_sp<-NULL
	how_may<-NULL
	order<-NULL
	for (sp in 1:length(species)){
		which_sp[[sp]]<-which(idr==species[sp])[sample(1:length(which(idr==species[sp])))]
		how_may[[sp]]<-length(which_sp[[sp]])/runs
		order<-c(order,which_sp[[sp]])
		}
	matrixnewbs<-cbind(matrix[order,],idr[order])
	attr<-rep(NA, nrow(matrix))
	for (giro in 1:runs){
		first<-NULL
		last<-NULL
		memberbs<-matrixnewbs[,3]
		for (spe in 1:length(species)){
			#spe<-2
			first[spe]<-round(1+((giro-1)*as.numeric(how_may[spe])))
			last[spe]<-round(giro*as.numeric(how_may[spe]),0)
			qualitot<-which(matrixnewbs[,3]==species[spe])
			if(first[spe]<=last[spe]){
			memberbs[qualitot[first[spe]:last[spe]]]<-0
			if(giro==runs){
			  memberbs[qualitot[first[spe]:qualitot[length(qualitot)]]]<-0
			}
			}
		}
		finalebs<-biodecrypt(matrixnewbs[,c(1,2)], memberbs, alpha=alpha, minimum=minimum, ratio=ratio,buffer=buffer, checkdist=checkdist,polygon=polygon, map=map,xlim=xlim,ylim=ylim,main=main) 
		attri<-which(memberbs==0)
		attr[attri]<-finalebs$table[attri,3]
	}
	res$cross<-cbind(matrixnewbs[,3],attr,rep(0,nrow(matrix)),rep(0,nrow(matrix)), matrixnewbs[,1:2])
	colnames(res$cross)<-c("original","predicted","MIR","NIR","Long","Lat")
	res$cross[which(res$cross[,2]>0 & res$cross[,1]!=res$cross[,2]),3]<-1
	res$cross[which(res$cross[,2]==0),4]<-1
	res$MIR<-(sum(res$cross[,3])/nrow(res$cross))*100
	res$NIR<-(sum(res$cross[,4])/nrow(res$cross))*100
	return(res)
}
