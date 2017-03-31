library("impute") #Availabe from bioconductor
library("glasso")
library("QUIC")
library("RColorBrewer")
library("igraph")


########## Functions:

#Calculation of custom robust correlation. This tries to discriminate between technical outliers (artefacts) and true biological outliers.
#See supplemental methods of Papasaikas, Tejedor et al for details.
#Arguments are M -> a n x p matrix (n features -e.g ASEs- in the rows, p variables -e.g KDs- in the columns)
#and ct a scalar defining the threshold for finding close neighbors of variables
CRobCor <- function (M,ct=0.5) {
Mcor<-cov((M),use="all.obs")
Mcor<-cov2cor(Mcor)
DELTA<-array(data=0,dim=c(ncol(M),ncol(M),nrow(M) ),dimnames=list(colnames(M),colnames(M),rownames(M)) )	# p x p x n  dimensional array (obviously ofr large numbers of n this becomes impossible to compute)
MRobCor<-Mcor;

for (ev in 1:nrow(M)) {
  marg_M<-M[-c(ev),];
  marg_M<-scale(marg_M) 
  marg_Mcor<-cov(marg_M,use="all.obs")
  marg_Mcor<-cov2cor(marg_Mcor)
  DELTA[,,ev]<-abs(marg_Mcor-Mcor)		#marginal Delta Covariance
}

INFL<-apply(DELTA,3,sum)/(ncol(M)^2-ncol(M)); #Calculate "influence" of outliers

for (p_i in 1:(ncol(M)-1)){
  for (p_j in (p_i+1):ncol(M)){
    if (abs(Mcor[p_i,p_j])<0.25) next               
    cta=ct-0.15 /(   0.5*length(which(abs(Mcor[p_i,])>ct))+ 0.5*length(which(abs(Mcor[p_j,])>ct))          )
    WT=apply(DELTA[p_i,  c(p_j,which(abs(Mcor[p_i,])>cta)),],  2,sum)/(length(which(abs(Mcor[p_i,])>cta)) )    #Calculate Weights
    WT=WT+apply(DELTA[p_j,  c(p_i,which(abs(Mcor[p_j,])>cta)),],  2,sum)/(length(which(abs(Mcor[p_j,])>cta)) )    #Calculate Weights 
    CRij=cov.wt(M[,c(p_i,p_j)], ((1/WT^(0.8))),cor=TRUE)$cor  #Weighted correlation
    MRobCor[p_i,p_j]=CRij[1,2]
    MRobCor[p_j,p_i]=CRij[2,1]
  }
}

return(MRobCor)
}


#  CENTRALITY/AUTHORITY MEASURES. Return a matrix of authority scores and corresponding node rankings
CentralityRanking<-function(g){
  PR<-page.rank(g)$vector  		#Pagerank Score
  OPR<-order(PR,decreasing=TRUE)
  
  BC<-betweenness(g)		#Betweeness Centrality
  OBC<-order(BC,decreasing=TRUE)	
  
  CC<-closeness(g)		#Closeness Centrality
  OCC<-order(CC,decreasing=TRUE)	
  
  RANK<-rank(PR)+rank(BC)+rank(CC)
  ORANK<-order(RANK,decreasing=TRUE)	#Ordering of Nodes based on aggregate score.
  return(rbind(PR,OPR,BC,OBC,CC,OCC,RANK,ORANK))
}
###############################################


###### Scale vector in 0-1 range
Vscale <- function(v, a, b) {
  v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a
}


##### Network Reconstruction, Plotting and Analysis. Input is a covariance matrix C and a regularization paramether rho
CreateGraph<-function(C,rho=0.5){
  
  RHO<-matrix(data=rho,nrow=nrow(C),ncol=ncol(C)) #Create reqularization matrix
  #InvCov<-glasso(s=C,rho=RHO,maxit=20000,penalize.diagonal=TRUE)$wi  #Graphical Lasso estimation of sparse inverse covariance matrix
  InvCov<-QUIC(C,rho=RHO,tol=1e-02,maxIter=100,msg=0)$X  #0.1 for C4/ 0.3 for C3 / 0.5 for C2
  
  ADJ<-abs(InvCov) #Weighted adjacency matrix (taking absolute values)
  TADJ<-InvCov     #"True" Weighted Adjacency matrix (keeping the sign)
  GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
  TGRAO<-graph.adjacency(TADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
  
  labels=rownames(C)
  GRAO<-set.vertex.attribute(GRAO,"labels",value=labels)
  
  ### Remove isolated nodes:
  g<-induced.subgraph(GRAO,which(degree(GRAO)>0)) 
  Tg<-induced.subgraph(TGRAO,which(degree(GRAO)>0)) 
  
  ### Calculate Centrality scores and node rankings
  Centr_G<-CentralityRanking(GRAO)
  
  ######## COMMUNITY DETECTION #########
  fg_C<-fastgreedy.community(g, modularity=TRUE) #infomap now available in igraph is better...
  memb=fg_C
  names(memb$membership)<-get.vertex.attribute(g,"labels")
  csize=table(memb$membership)
  
  ##################################### Plotting Parameters for the Graph Object
  comps <-memb$membership
  colbar <- c(brewer.pal(8,"Dark2")[2:3],brewer.pal(8,"Set1"),brewer.pal(8,"Accent")[1],brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(8,"Set2"))
  V(g)$color <- colbar[comps]
  
  ##### ADD TRANSPARENCIES #####
  E(g)$width<-Vscale(E(g)$weight,1.5,9)
  V(g)$color <- colbar[comps]
  V(g)$color<-paste(V(g)$color,"44",sep="")
  E(g)$color<-ifelse(E(Tg)$weight<0,"#55885599","#BB111199")
  V(g)$label_color="#222222"
  l<-layout.fruchterman.reingold(g)
  
  plot(g, layout=l, vertex.size=50,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=V(g)$labels,vertex.label.dist=3,vertex.label.cex=0.7,vertex.label.font=4,
       vertex.label.color=V(g)$label_color,asp=0.72)
  ##########
  return(list(g, memb, Centr_G))
}





############ DATA FILTERING AND PREPROCESSING:
#Network Data. These are already transformed to Z-scores + Pvalues. This file is generated from the prepare_z_extra_data_2016.pl script.
Ztable<-read.table("LABCHIPS/labchip37_26_06_13_bare_Z_P.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""); 
#Splicing Annotation Data:
Annot<-read.table("main_annotation.txt",as.is=TRUE,na.strings="NAn",header=FALSE,row.names=1,sep="\t",quote=""); 

#Split data to a Z-values and a P-values matrix:
Zmatrix<-as.matrix(Ztable[,seq(1,(ncol(Ztable)-1),2)]);  
Pmatrix<-as.matrix(Ztable[,seq(2,ncol(Ztable),2)]);		

#Remove Uninformative genes. These are genes whose inclusion does not change in the HeLa context .
Zmatrix<-(Zmatrix[,-c(20,24,28,36)]);             #Do not include VEGFA (20th), SYK (24th), CND1 (28th), ???GADD45A (31st), LMNA (36)
Pmatrix<-as.matrix(Pmatrix[,-c(20,24,28,36)]);		#Do not include VEGFA (20th), SYK (24th), CCND1 (28th), ???GADD45A (31st), LMNA (36)

# Remove Sparse KD data and impute remaining missing values:
Zvals<-Zmatrix
SparseRows<-which(rowMeans(is.na(Zvals))*100>30)#Row numbers with >30% missing values
if (length(SparseRows)>0){
  Z<-Z[-c(SparseRows),]
}
Zvals<-impute.knn(as.matrix(Zvals),k=round(sqrt(nrow(Zvals))*0.25),rowmax=0.5)$data

### Transpose matrix to ASEs (rows) x KDs (columns)
t_Zvals=t(Zvals)
#Remove conditions not corresponding to KDs (e.g untransfected, mocks...). This is based on a collection of suffixes/prefixes, since data kept being added to the Labchip this kept expanding...
NonKDs=c("Untransfected","untransfected","Mock","CN ","CN_","CL_","GFP_","Control","_CN","DMSO","UNT_","CTRL_","ctr_4h","ctr_24h","ctr_37h","Contr_","Scrambled","0h","rep")
for (i in 1: length(NonKDs)){
word=NonKDs[i]
if (length(grep(word,rownames(Zvals) ) )>0  )  {   t_Zvals<-t_Zvals[1:ncol(Zvals),-c(grep(word,colnames(t_Zvals)))]       }   
}

#Knock Down (Column) Scaling. 
t_Zvals<-scale(t_Zvals)    #  Column Scaling -> Only Shape not Scale of KD effect is important. This is essential since for example not all KDs have the same efficiecny 

#Calculate Robust Correlation Matrix
C <- CRobCor(t_Zvals)

dev.off()
#Create graph, calculate centrality measures, identifu communities and plot
CreateGraph(C,rho=0.43)   # Rho was set to achieve FDR < 5%







































