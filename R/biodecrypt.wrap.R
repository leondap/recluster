biodecrypt.wrap<-function(mat,id,alpha=c(1,5,10,15),alphamat=NULL,ratio=c(2,3,4,5),buffer=c(0,40,80,120,160), checkdist=T, polygon=NULL, minimum=7, map=NULL,xlim=NULL,ylim=NULL,main=NULL,save=T,name="res_cross.txt",runs=10){
	res<-NULL
	taxa<-max(id)
	al<-length(alpha)
	if(!is.null(alphamat)){
		al<-ncol(alphamat)
	}
	rat<-length(ratio)
	buf<-length(buffer)
	res_cross<-matrix(NA,al*rat*buf,6)
	if(!is.null(alphamat)){
		res_cross<-matrix(NA,al*rat*buf,(5+nrow(alphamat)))
	}

	riga<-1
	for(alphav in 1:al){
		for(ratiov in 1:rat){
			for (bufferv in 1:buf){
				alphause<-rep(alpha[alphav],taxa)
				if(!is.null(alphamat)){
					alphause<-alphamat[,alphav]
					addcol<-ncol(alphamat)
				}
				print(c(alphav,ratiov,bufferv))
				cross<-biodecrypt.cross(mat, id, ratio=ratio[ratiov],buffer=buffer[bufferv],alpha=alphause, checkdist=checkdist,map=map, polygon=polygon, runs=runs, test=T) 
				res_cross[riga,4]<-ratio[ratiov]
				res_cross[riga,5]<-buffer[bufferv]
				if(is.null(alphamat)){
					res_cross[riga,6]<-alphause[1]
					addcol<-1
				}
				if(!is.null(alphamat)){
					res_cross[riga,6:ncol(res_cross)]<-alphause
				}
				res_cross[riga,1]<-cross$MIR
				res_cross[riga,2]<-cross$NIR
				res_cross[riga,3]<-cross$NUR
				if(save){
					write.table(res_cross,name)
				}
				riga<-riga+1
			}
		}
	}
colnames(res_cross)<-c("MIR","NIR","NUR","ratio","buffer",rep("alpha",addcol))
res$table<-res_cross
return(res)
}
