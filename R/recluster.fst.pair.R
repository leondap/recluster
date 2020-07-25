recluster.fst.pair<-function(dist,vect,setzero=F,setnazero=F){
	res<-NULL
	populations<-vect[!duplicated(vect)]
	matrix<-as.matrix(dist)
	matrix2<-matrix(NA, length(populations),length(populations))
	colnames(matrix2)<-populations
	rownames(matrix2)<-populations
	matrix3<-matrix2
	matrix4<-matrix2
	matrix5<-matrix2
	phist<-recluster.fst(dist,vect,setzero=setzero, setnazero=setnazero)		
	for(col in 1:(length(populations)-1)){
		for(row in (col+1) :length(populations)){
			sum<-NULL
			vectsub<-vect[which(vect==populations[row]|vect==populations[col])]
			subm<-matrix[which(vect==populations[row]|vect==populations[col]),which(vect==populations[row]|vect==populations[col])]
			submfst<-recluster.fst(as.dist(subm), vectsub,setzero=setzero, setnazero=setnazero)
			matrix2[row,col]<-submfst$Dst/phist$Ht
			matrix3[row,col]<-submfst$Dst
			matrix4[row,col]<-(submfst$Dst/phist$Ht)/(((1-submfst$Hs))/(1+submfst$Hs))
			matrix5[row,col]<-((submfst$Dst)/(1-submfst$Hs))*2
			if (setnazero & phist$Ht==0){
				matrix2[row,col]<-0
				matrix4[row,col]<-0
			}
		}
	}		
	res$Gstm<-as.dist(matrix2)
	res$Dstm<-as.dist(matrix3)
	res$G1stm<-as.dist(matrix4)
	res$Dm<-as.dist(matrix5)
	return (res)
}