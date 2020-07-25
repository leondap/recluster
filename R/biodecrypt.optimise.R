biodecrypt.optimise<-function(tab,coef=c(2,1,1), penalty=10){
  	res<-NULL
	val<-tab[,1]^coef[1]+tab[,2]^coef[2]+tab[,3]^coef[3]
	mini<-min(val)
	quali<-which((val-mini)<penalty)
	if(length(quali)>1){
		tot<-sum(1/val[quali])
		weight<-tab[quali,(4:ncol(tab))]/val[quali]
		res$ratio<-(sum(weight[,1]))/tot
		res$buffer<-(sum(weight[,2]))/tot
		res$MIR<-mean(tab[quali,1])
		res$NIR<-mean(tab[quali,2])
		res$NUR<-mean(tab[quali,3])
		quantialpha<-ncol(weight)-2
		res$alpha<-NULL
			for(q in 1:quantialpha){
			res$alpha[q]<-(sum(weight[,(2+q)]))/tot
		}
	}
	if(length(quali)==1){
		res$ratio<-tab[quali,4]
		res$buffer<-tab[quali,5]
		res$MIR<-tab[quali,1]
		res$NIR<-tab[quali,2]
		res$NUR<-tab[quali,3]
		quantialpha<-ncol(tab)-5
		res$alpha<-NULL
			for(q in 1:quantialpha){
			res$alpha[q]<-tab[quali,(5+q)]	
		}
	}
	return(res)
}