recluster.fst<-function(dist,vect,setzero=F,setnazero=F){
	res<-NULL
	populations<-vect[!duplicated(vect)]
	spec<-length(populations)
	vectmat<-as.vector(dist)
	Ht<-mean(dist)
	res$Ht<-Ht
	res$lengthHt<-length(dist)
	omol<-as.vector(dist(vect))
	omol[which(omol>0)]<-1
	res$lengthHs<-length(vectmat[which(omol==0)])
	Hs<-mean(vectmat[which(omol==0)],na.rm=T)
	res$Hs<-Hs
	res$Gst<-((Ht-Hs)/(Ht))
	res$G1st<-((Ht-Hs)/(Ht))/(((spec-1)*(1-Hs))/(spec-1+Hs))
	res$Dst<-Ht-Hs
	res$D<-((Ht-Hs)/(1-Hs))*(spec/(spec-1))
	if (setzero){
		res$Gst[res$Gst<0]<-0
		res$Dst[res$Dst<0]<-0
		res$D[res$D<0]<-0
		res$G1st[res$G1st<0]<-0
		}
	if(setnazero) {
		if(Ht==0){
			res$Gst<-0
			res$G1st<-0
		}
	}
	return(res)
}

