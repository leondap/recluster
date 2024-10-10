recluster.region<-function(mat,tr=50,dist="simpson",method="ward.D2",members=NULL, phylo=NULL, mincl=2,maxcl=10,rettree=FALSE,retmat=FALSE,retmemb=FALSE){
 	res<-NULL
	clusters<-maxcl-mincl+1
	mat2<-as.matrix(mat)    
	rows<-nrow(mat2)    
	tab2<-array(NA,dim=c(rows,tr,clusters))
	rownames(tab2)<-rownames(mat2)        
	rownames(mat2)<-c(1:rows)
	if(data.class(mat)=="dist"){distam<-mat
		}else{
		distam<-recluster.dist(mat2, phylo=phylo, dist=dist)
	}
	dista<-as.matrix(distam)                
	for(i in 1:tr){ 
		samp<-sample(as.numeric(as.character(rownames(mat2))))
		dista2<-dista[samp,samp]
		if(method=="pam"){
 			tree<-NULL}else{
			if (method=="diana"){tree<-diana(as.dist(dista2))}else
				{tree<-hclust(as.dist(dista2),members=members,method=method)}
			}
			for (cut in mincl:maxcl){
				if(method=="pam"){
					cuttr<-pam(as.dist(dista2),k=cut)$clustering}else
				if(method=="diana"){
					cuttr<-cutree (as.hclust(tree), k=cut)}
					else{                
					cuttr<-cutree (tree, k=cut)}
					tab2[,i,(cut-mincl+1)]<-cuttr[match(rownames(tab2),names(cuttr))]
			}
		}
		tree<-NULL
		cuttr<-NULL
		matrices<-array(NA,dim=c(rows,rows,clusters))
		for(sel in 1:clusters){
			tabsel<-tab2[,,sel]
			tabsel<-as.data.frame(tabsel)
			tabsel[] <- lapply(tabsel, factor)
			matrices[,,sel]<-as.matrix(daisy(tabsel,metric="gower"))
		}
		tabsel<-NULL
		if(retmemb){res$memb<-tab2}
		tab2<-NULL
		if(retmat){res$matrices<-matrices}
		pamsol<-matrix(data=NA, nrow=rows,ncol=clusters)
		colnames(pamsol)<-c(mincl:maxcl)
		rownames(pamsol)<-rownames(mat2)
		res$solutions<-matrix(data=NA, nrow=clusters,ncol=3)
		colnames(res$solutions)<-c("k","silh","ex.diss")
		res$solutions[,1]<-c(mincl:maxcl)
		for (pamr in 1:clusters){
			if(method=="pam"){pamsol[,pamr]<-pam(as.dist(matrices[,,pamr]),k=mincl-1+pamr)$clustering}
			else{
			if (method=="diana"){
				pami<-diana(as.dist(matrices[,,pamr]))
				pamsol[,pamr]<-cutree(as.hclust(pami),k=mincl-1+pamr)
				}else{
				pami<-hclust(as.dist(matrices[,,pamr]),method=method)}
				pamsol[,pamr]<-cutree(pami,k=mincl-1+pamr)
				if(rettree){res$tree[[pamr]]<-pami}                                                                    
			}
			pami<-NULL
			res$solutions[pamr,3]<-recluster.expl(distam,pamsol[,pamr])
			res$solutions[pamr,2]<-mean(silhouette(pamsol[,pamr],dist=distam)[,3])
	}
	res$grouping<-pamsol
	return(res)
}
