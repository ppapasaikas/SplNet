library("impute")
library("glasso")
library("igraph")
library("gtools")
library("RColorBrewer")
library("maptools")
library("gplots")
library("plotrix")  #use this for gapped plots
source("/Volumes/HD2/RESEARCH/R/heatmap.3.R")   #http://www.biostars.org/post/show/18211/how-do-i-draw-a-heatmap-in-r-with-both-a-color-key-and-multiple-color-side-bars/

hclust2 <- function(x, method="ward", ...)
  hclust(x, method=method, ...)

dist2 <- function(x, ...)
  as.dist(1-(cor(t(x), method="pearson"))) #cor OR abs(cor) for distance function???

Vscale <- function(v, a, b) {
  v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a
}

#######  CENTRALITY/AUTHORITY MEASURES ######
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

setwd("/Volumes/HD2/RESEARCH/RAMON/");
Ztable<-read.table("LABCHIPS/labchip37_26_06_13_bare_Z_P.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote="");  #Added  Mk6, Fas-HITS, StressTreat, HEK

Ztable<-read.table("LABCHIPS/labchip37_2016_434_457_Z_P.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote="");  #Added  Mk6, Fas-HITS, StressTreat, HEK


Annot<-read.table("main_annotation.txt",as.is=TRUE,na.strings="NAn",header=FALSE,row.names=1,sep="\t",quote=""); 

CHROMP<-read.table("LABCHIPS/labchip_chrom_column.txt",as.is=TRUE,na.strings="NAn",header=TRUE,sep="\t",quote="");
dup=which(duplicated(CHROMP[,2]))
CHR=CHROMP[-c(dup),]
rownames(CHR)=CHR[,2]

S_COMB<-as.matrix(read.table("STRING/S_COMB.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_NEIGH<-as.matrix(read.table("STRING/S_NEIGH.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_FUSION<-as.matrix(read.table("STRING/S_FUSION.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_COOC<-as.matrix(read.table("STRING/S_COOC.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_HOMOL<-as.matrix(read.table("STRING/S_HOMOL.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_COEXP<-as.matrix(read.table("STRING/S_COEXP.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_EXPER<-as.matrix(read.table("STRING/S_EXPER.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_KNOWL<-as.matrix(read.table("STRING/S_KNOWL.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));
S_TEXTM<-as.matrix(read.table("STRING/S_TEXTM.tsv",as.is=TRUE,na.strings="NAn",header=TRUE,row.names=1,sep="\t",quote=""));

Zmatrix<-as.matrix(Ztable[,seq(1,(ncol(Ztable)-1),2)]);  
Pmatrix<-as.matrix(Ztable[,seq(2,ncol(Ztable),2)]);		

Zmatrix<-(Zmatrix[,-c(20,24,28,36)]);    #Do not include VEGFA (20th), SYK (24th), ??CCND1 (28th), ???GADD45A (31st), LMNA (36)
Pmatrix<-as.matrix(Pmatrix[,-c(20,24,28,36)]);		#Do not include VEGFA (20th), SYK (24th), ??CCND1 (28th), ???GADD45A (31st), LMNA (36)

thresh=0.25;
STRING<-pmax(S_FUSION,S_EXPER);
STRING[which(STRING<thresh)]=0;
colnames(STRING)=rownames(STRING)

Z<-Zmatrix
SparseRows<-which(rowMeans(is.na(Z))*100>30)#Row numbers with >30% missing values
if (length(SparseRows)>0){
Z<-Z[-c(SparseRows),]
}
Z<-impute.knn(as.matrix(Z),k=round(sqrt(nrow(Z))*0.25),rowmax=0.5)$data


XFINAL<-t(Z)
if (length(grep("Untransfected",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("Untransfected",colnames(XFINAL)))]       }   
if (length(grep("untransfected",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("untransfected",colnames(XFINAL)))]       }     
if (length(grep("Mock",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("Mock",colnames(XFINAL)))]       }                                                      
if (length(grep("CN ",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("CN ",colnames(XFINAL)))]        }
if (length(grep("CN_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("CN_",colnames(XFINAL)))]        }
if (length(grep("CL_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("CL_",colnames(XFINAL)))]        }
if (length(grep("GFP_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("GFP_",colnames(XFINAL)))]        }
if (length(grep("Control",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("Control",colnames(XFINAL)))]        }
if (length(grep("_CN",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("_CN",colnames(XFINAL)))]        }
if (length(grep("DMSO",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("DMSO",colnames(XFINAL)))]        }
if (length(grep("UNT_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("UNT_",colnames(XFINAL)))]        }
if (length(grep("CTRL_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("CTRL_",colnames(XFINAL)))]        }
if (length(grep("ctr_4h|ctr_24h|ctr_37",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("ctr_4h|ctr_24h|ctr_37",colnames(XFINAL)))]        }
if (length(grep("Contr_",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("Contr_",colnames(XFINAL)))]        }
if (length(grep("Scrambled",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("Scrambled",colnames(XFINAL)))]  }
if (length(grep("0h",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("0h",colnames(XFINAL)))]  }
if (length(grep("rep",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("rep",colnames(XFINAL)))]  }
#if (length(grep("SF3B1_LV",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("SF3B1_LV",colnames(XFINAL)))]  }
#if (length(grep("demycin",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("demycin",colnames(XFINAL)))]  }
#if (length(grep("ycin",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("ycin",colnames(XFINAL)))]  }
#if (length(grep("getin",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("getin",colnames(XFINAL)))]  }
#if (length(grep("orin",rownames(Z) ) )>0  )  {   XFINAL<-XFINAL[1:ncol(Z),-c(grep("orin",colnames(XFINAL)))]  }


SCALED=scale(t(XFINAL))
XFINAL<-scale(XFINAL)    #  Column Scaling -> Only Shape not Scale of SF effect is important

dev.off()
par(mar=c(7.1,5.1,3.1,2.6))
SF1="HNRPC" #CDC5L OR negative NCBP2 or SNRPG
SF2="U2AF2" #PLRG1 OR negative AOF2 or DDX52
SF3="PLRG1"
plot(Z[SF1,],lwd=4,col="darkred",type="l",xlab="",ylab="Z-score",col.axis="red",bty="o",bty="o",xaxt="n",cex.lab=1.25)  #draw ar 9.5x4.2 for pdf
axis(1, at=1:ncol(Z), lab=colnames(Z),las=3) #but better from maptools:

dev.off()
par(mar=c(7.1,5.1,3.1,2.6))
plot(XFINAL[,SF1],lwd=4,col="darkred",type="l",xlab="",ylab="Scaled Z-score",col.axis="red",bty="o",bty="o",xaxt="n",cex.lab=1.25,ylim=c(-3,3))  #draw ar 9.5x4.2 for pdf
axis(1, at=1:ncol(Z), lab=colnames(Z),las=3) #but better from maptools:
lines(XFINAL[,SF2],lwd=4,col="darkgreen",lty=6)  #draw ar 9.5x4.2 for pdf
legend("topleft", inset=.01,c(SF1,SF2), col=c("darkred","darkgreen") , lwd=c(3.5,3.5) , lty=c(1,6)  , horiz=FALSE,bty="n",seg.len=2.5,y.intersp=1.5)

dev.off()
par(mar=c(7.1,5.1,3.1,2.6))
plot(XFINAL[,SF1],lwd=4,col="darkred",type="l",xlab="",ylab="Scaled Z-score",col.axis="red",bty="o",bty="o",xaxt="n",cex.lab=1.25,ylim=c(-2.2,2.5),yaxt="n")  #draw ar 9.5x4.2 for pdf
axis(1, at=1:ncol(Z), lab=colnames(Z),las=3) #but better from maptools:
axis(2, at=1:length(colnames(Z)), lab=colnames(Z),las=3) #but better from maptools:
lines(XFINAL[,SF2],lwd=4,col="darkgreen",lty=6)  #draw ar 9.5x4.2 for pdf
lines(XFINAL[,SF3],lwd=4,col="tan",lty=9)  #draw ar 9.5x4.2 for pdf
legend("topleft", inset=.01,c(SF1,SF2), col=c("darkred","darkgreen") , lwd=c(3.5,3.5) , lty=c(1,6)  , horiz=FALSE,bty="n",seg.len=2.5,y.intersp=1.5)

dev.off()
par(mar=c(5.1,5.1,3.1,2.6))
plot(XFINAL[,SF1],XFINAL[,SF2],xlab=paste("Z-score ", SF1),ylab=paste("Z-score ", SF2),col.lab="darkslategrey",col="steelblue",pch=16,cex.lab=1.25,cex=1.5)
pointLabel(XFINAL[,SF1],XFINAL[,SF2],labels=rownames(XFINAL),cex=0.95)

fit <- lm(XFINAL[,SF2]~XFINAL[,SF1])
abline(fit,col="darkred",lwd=3,lty="dashed")
legend("topleft", inset=.01,c(paste("Cor =",signif(cor(XFINAL[,SF2],XFINAL[,SF1]),3) )), col=c("darkred") , lwd=c(3) , lty=c("dashed")  , horiz=FALSE)





INFLNC=vector(mode="numeric",length=nrow(XFINAL))
for (k in 1:nrow(XFINAL)){
  INFLNC[k]=cor(XFINAL[,SF2],XFINAL[,SF1])-cor(XFINAL[-c(k),SF2],XFINAL[-c(k),SF1])
}
names(INFLNC)=rownames(XFINAL)


CV_N<-cor((XFINAL),use="all.obs")
##########################################C#U#S#T#O#M# # #R#O#B#U#S#T# # #C#O#V#A#R#I#A#N#C#E############################################
XF<-XFINAL
J_CD<-cov((XF),use="all.obs")
J_CD<-cov2cor(J_CD)
DELTA<-array(data=0,dim=c(ncol(XF),ncol(XF),nrow(XF) ),dimnames=list(colnames(XF),colnames(XF),rownames(XF)) )	# p x p x n  dimensional array
CRD<-J_CD;
ct=0.5; #0.4 FOR 39(37) Events DN/ 0.5 FOR 39(37) Events SN

for (ev in 1:nrow(XF)) {
M_XF<-XF[-c(ev),];
M_XF<-scale(M_XF) ###!!!!!!!!!CORRECT THIS!!!!!!!!!!!!!!!
M_CD<-cov(M_XF,use="all.obs")
M_CD<-cov2cor(M_CD)
M_CD<-abs(M_CD-J_CD)		#marginal Delta Covariance
DELTA[,,ev]=M_CD
}
INFL<-apply(DELTA,3,sum)/(ncol(XFINAL)^2-ncol(XFINAL));
for (p_i in 1:(ncol(XF)-1)){
	for (p_j in (p_i+1):ncol(XF)){
	if (abs(CV_N[p_i,p_j])<0.25) next               
  cta=ct-0.15 /(   0.5*length(which(abs(J_CD[p_i,])>ct))+ 0.5*length(which(abs(J_CD[p_j,])>ct))          )
  WT=apply(DELTA[p_i,  c(p_j,which(abs(J_CD[p_i,])>cta)),],  2,sum)/(length(which(abs(J_CD[p_i,])>cta)) )    #/INFL
  WT=WT+apply(DELTA[p_j,  c(p_i,which(abs(J_CD[p_j,])>cta)),],  2,sum)/(length(which(abs(J_CD[p_j,])>cta)) )    #/INFL 
  CRij=cov.wt(XF[,c(p_i,p_j)], ((1/WT^(0.8))),cor=TRUE)$cor  #
	CRD[p_i,p_j]=CRij[1,2]
	CRD[p_j,p_i]=CRij[2,1]
	}
}
############################################################################################################################
############################################################################################################################
############################################################################################################################
p_i=which(rownames(CRD)==SF1)
p_j=which(rownames(CRD)==SF2)
cta=ct-0.15 /(   0.5*length(which(abs(J_CD[p_i,])>ct))+ 0.5*length(which(abs(J_CD[p_j,])>ct))          )
WT=apply(DELTA[p_i,  c(p_j,which(abs(J_CD[p_i,])>cta)),],  2,sum)/(length(which(abs(J_CD[p_i,])>cta)) )    # /INFL
WT=WT+apply(DELTA[p_j,  c(p_i,which(abs(J_CD[p_j,])>cta)),],  2,sum)/(length(which(abs(J_CD[p_j,])>cta)) )    # /INFL 
WT=1/WT^(0.8)
fit <- lm(XFINAL[,SF2]~XFINAL[,SF1],weights=WT)
abline(fit,col="darkgreen",lwd=3,lty="dashed")
legend("topleft", inset=.01,c(paste("ClassCor =",signif(CV_N[p_i,p_j],2) ),paste("RobCor =",signif(CRD[p_i,p_j],2) )), col=c("darkred","darkgreen") , lwd=c(3,3) , lty=c("dashed","dashed")  , horiz=FALSE)

C<-CRD



###### Only For Files with additional entries:
ExtraNames=setdiff(rownames(CRD),rownames(STRING))
SMAT=matrix(data=0,nrow=nrow(STRING),ncol=length(ExtraNames),dimnames=list(rownames(STRING),ExtraNames ) )
STRING=cbind(STRING,SMAT)
SMAT=matrix(data=0,nrow=length(ExtraNames),ncol=ncol(STRING),dimnames=list(ExtraNames,colnames(STRING))   )
STRING=rbind(STRING,SMAT)
STRING=STRING[rownames(C),rownames(C)]



commonnames=intersect(rownames(STRING),rownames(CRD))
C<-CRD[commonnames,]
C<-C[,commonnames]
Z<-Z[commonnames,]
STRING<-STRING[commonnames,]
STRING<-STRING[,commonnames]
rownames(STRING)=c()
colnames(STRING)=c()
###### 

















##### Network Reconstruction, Plotting and Analysis#
#RHO<-matrix(data=0.481,nrow=nrow(C),ncol=ncol(C)) #0.522 (DoubleNorm ) / 0.53 (SingleNorm)

RHO<-matrix(data=0.43,nrow=nrow(C),ncol=ncol(C)) #0.522 (DoubleNorm ) / 0.53 (SingleNorm)

GF<-glasso(s=C,rho=RHO,maxit=20000,penalize.diagonal=TRUE)  #Graphical Lasso estimation of sparse inverse covariance matrix

ADJ<-abs(GF$wi) #Weighted adjacency matrix (taking absolute values)
TADJ<-GF$wi     #"True" Weighted Adjacency matrix (keeping the sign)
GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
SG<-graph.adjacency(STRING,mode=c("max"),weighted=TRUE,diag=FALSE)
TGRAO<-graph.adjacency(TADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
TSG<-graph.adjacency(-STRING,mode=c("max"),weighted=TRUE,diag=FALSE)

labels=rownames(C)
GRAO<-set.vertex.attribute(GRAO,"labels",value=labels)
SG<-set.vertex.attribute(SG,"labels",value=labels)
g<-induced.subgraph(GRAO,which(degree(GRAO)>0)) 
sg<-induced.subgraph(SG,which(degree(GRAO)>0)) 
sg<-graph.intersection(g,sg)
sg<-set.vertex.attribute(sg,"labels",value=V(g)$labels)

Tg<-induced.subgraph(TGRAO,which(degree(GRAO)>0)) 
Centr_G<-CentralityRanking(GRAO)



######## COMMUNITY DETECTION #########
fg_C<-fastgreedy.community(g, modularity=TRUE) 
memb=fg_C
names(memb$membership)<-get.vertex.attribute(g,"labels")
csize=table(memb$membership)


##################################### Plotting Parameters for the Graph Object
comps <-memb$membership
colbar <- c(brewer.pal(8,"Dark2")[2:3],brewer.pal(8,"Set1"),brewer.pal(8,"Accent")[1],brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(8,"Set2"))
V(g)$color <- colbar[comps]
V(sg)$color <- colbar[comps]




##### ADD TRANSPARENCIES #####
dev.off()

E(g)$width<-Vscale(E(g)$weight,1.5,9)
V(g)$color <- colbar[comps]
V(g)$color<-paste(V(g)$color,"44",sep="")
E(g)$color<-ifelse(E(Tg)$weight<0,"#55885599","#BB111199")


l<-layout.fruchterman.reingold(g)

plot(g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)
##########


V(g)$color <- colbar[comps]
Vs<-which(get.vertex.attribute(g,"labels")%in%rownames(Pmatrix));
#Vs<-which(get.vertex.attribute(g,"labels")%in%rownames(Pmatrix)[which(Pmatrix[,1]<10.01)]);
#Vs<-which(get.vertex.attribute(g,"labels")%in%rownames(Pmatrix)[which(Pmatrix[,1]<0.05)] & get.vertex.attribute(g,"labels")%in%rownames(Pmatrix)[which(rowMeans(Pmatrix)>0.45)]);

gs <- induced.subgraph(sg, Vs)
Tgs <- induced.subgraph(sg, Vs)
E(gs)$color<-"#222244"
E(gs)$width<-1.75


ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
MIN=0.1;
if (min(V(gs)$size)>0.1){
MIN=min(V(gs)$size)
}
V(gs)$size<-Vscale(abs(V(gs)$size), 16, 16 * sqrt( max(V(gs)$size)/MIN  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l[Vs,]

gLABELS=V(gs)$labels

plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=gLABELS,vertex.label.dist=3,vertex.label.cex=0.8,vertex.label.font=4,
vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)































tkplot(gs, layout=ls,vertex.size=4,vertex.label=V(gs)$labels,vertex.label.font=4,
        vertex.label.color=V(gs)$label_color,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666")
LTK=tkplot.getcoords(8)    #Manually correct and get the coordinates back
#display.brewer.pal(8,"Accent") 



####### REPLOT USING tkplot manual correction #######
#LOAD LTK
#LTK=load("LayoutTK_rho0481_vsize830_asp1.layout")
dev.off()
plot(g, layout=LTK, vertex.size=100,rescale=FALSE, xlim=range(LTK[,1]),ylim=range(LTK[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=1)

E(gs)$width<-1.75


ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
MIN=0.1;
if (min(V(gs)$size)>0.1){
  MIN=min(V(gs)$size)
}
V(gs)$size<-Vscale(abs(V(gs)$size), 830, 830 * sqrt( max(V(gs)$size)/MIN  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-LTK[Vs,]

gLABELS=V(gs)$labels

BigLAB=gLABELS[which(apply(abs(ScZ),2,median,na.rm=TRUE)>1.5 & degree(g)<18)]
gLABELS[which(degree(g)>9)]=""
gLABELS[112]="DHX57"
gLABELS[which(apply(abs(ScZ),2,median,na.rm=TRUE)>1.5 & degree(g)<18) ]=BigLAB
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=gLABELS,vertex.label.dist=70,vertex.label.cex=0.8,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)
#Save PDF of dimensions 15x11
#####REPEAT CORRECTION UNTIL GOOD THEN SAVE THAT SHIT:
#save(LTK,file="LayoutTK_rho0481_vsize830_asp1.layout")



####### REPLOT USING tkplot manual correction AND transparencies to highlight specific components ####### 
#LOAD LTK
#LTK=load("LayoutTK_rho0481_vsize830_asp1.layout")
####Make transparent
dev.off()
E(g)$color<-ifelse(E(Tg)$weight<0,"#55885544","#BB111144")
plot(g, layout=LTK, vertex.size=100,rescale=FALSE, xlim=range(LTK[,1]),ylim=range(LTK[,2]),vertex.frame.color=paste(colbar[comps+1],"22",sep=""),vertex.label=NA,asp=1)
E(g)$color<-ifelse(E(Tg)$weight<0,"#558855BB","#BB1111BB")

V(g)$color <- colbar[comps+1]
E(gs)$width<-1.75


ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
MIN=0.1;
if (min(V(gs)$size)>0.1){
  MIN=min(V(gs)$size)
}
V(gs)$size<-Vscale(abs(V(gs)$size), 830, 830 * sqrt( max(V(gs)$size)/MIN  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-LTK[Vs,]

gLABELS=V(gs)$labels

V(gs)$tcolor<-paste(V(gs)$color,"44",sep="")
E(gs)$tcolor<-paste(E(gs)$color,"44",sep="")

BigLAB=gLABELS[which(apply(abs(ScZ),2,median,na.rm=TRUE)>1.5 & degree(g)<18)]
gLABELS[which(degree(g)>9)]=""
gLABELS[112]="DHX57"
gLABELS[which(apply(abs(ScZ),2,median,na.rm=TRUE)>1.5 & degree(g)<18) ]=BigLAB

plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=NA,vertex.label.dist=70,vertex.label.cex=0.8,vertex.label.font=4,
     ,vertex.color=V(gs)$tcolor,edge.color=E(gs)$tcolor, add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#66666644",asp=1)

TLATE=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Late Spl.")]
TMID=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Mid Spl.")]
TEARLY=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Early Spl.")]
PERS=rownames(DATA)[which(DATA[,"CLASS2"]=="Persist. Spl.")]

PlNod=which(V(gs)$labels%in%TEARLY)
CRg=subgraph(gs,PlNod-1)
CRgc=subgraph(g,PlNod-1)

plot(CRgc, layout=ls[PlNod,],vertex.size=1,vertex.label=NA, add=TRUE,rescale=FALSE,vertex.frame.color="#666666",asp=1)


plot(CRg, layout=ls[PlNod,],vertex.size=V(CRg)$size,vertex.label=V(CRg)$labels,vertex.label.dist=70,vertex.label.cex=0.8,vertex.label.font=4,
     vertex.label.color=V(CRg)$label_color,vertex.color=V(CRg)$color, add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)
#Save PDF of dimensions 15x11





##########   PLOT  ONLY THE DENSE CENTRAL CORE ##########
gLABELS=V(gs)$labels
LAB=gLABELS[which(degree(g)>7)]
if (length(grep("DHX57",LAB) )>0  )  {LAB=LAB[-c(grep("DHX57",LAB))]}
DENSE_v=match(LAB, V(g)$labels)-1
DENSE_g=subgraph(g,DENSE_v)
dev.off()
E(DENSE_g)$width<-Vscale(E(DENSE_g)$weight,1,10)
lsd<-LTK[match(LAB, V(g)$labels),]
#lsd<-layout.fruchterman.reingold(DENSE_g,params=list(area=vcount(DENSE_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(DENSE_g)^2.3)
plot(DENSE_g, layout=lsd, vertex.size=10,rescale=FALSE, xlim=range(lsd[,1]),ylim=range(lsd[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=1)

gsd <- subgraph(sg, DENSE_v)
Tgsd <- subgraph(Tg, DENSE_v)
Tgsd <- subgraph(sg, DENSE_v)
E(gsd)$color<-ifelse(E(Tgsd)$weight<0,"#558855","#BB1111")
E(gsd)$color<-"#222244"
E(gsd)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gsd,"labels")]
V(gsd)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gsd)$size<-Vscale(abs(V(gsd)$size), 140, 140 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )
V(gsd)$label_color="#222222"
CHRnodes=which(CHR[V(gsd)$labels,]==1)
V(gsd)$label_color[CHRnodes]="DarkRed"

plot(gsd, layout=lsd,vertex.size=V(gsd)$size,vertex.label=V(gsd)$labels,vertex.label.dist=15,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gsd)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)



##########   PLOT  ONLY SPL COMPONENTS ##########
namevector=Lsm
namevector=c("EHMT2","SUV39H1")
PlNod=which(V(gs)$labels%in%namevector)
DENSE_v=PlNod
DENSE_g=induced.subgraph(g,DENSE_v)

dev.off()

E(DENSE_g)$width<-Vscale(E(DENSE_g)$weight,1,10)

lsd<-LTK[PlNod,]

#lsd<-layout.fruchterman.reingold(DENSE_g,params=list(area=vcount(DENSE_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(DENSE_g)^2.3)
plot(DENSE_g, layout=lsd, vertex.size=10,rescale=FALSE, xlim=range(lsd[,1]),ylim=range(lsd[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=1)

gsd <- subgraph(sg, DENSE_v)
Tgsd <- subgraph(Tg, DENSE_v)
Tgsd <- subgraph(sg, DENSE_v)
E(gsd)$color<-ifelse(E(Tgsd)$weight<0,"#558855","#BB1111")
E(gsd)$color<-"#222244"
E(gsd)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gsd,"labels")]
V(gsd)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gsd)$size<-Vscale(abs(V(gsd)$size), 340, 340 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )
V(gsd)$label_color="#222222"
CHRnodes=which(CHR[V(gsd)$labels,]==1)
V(gsd)$label_color[CHRnodes]="DarkRed"

plot(gsd, layout=lsd,vertex.size=V(gsd)$size,vertex.label=V(gsd)$labels,vertex.label.dist=15,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gsd)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)

###Manually correct with tkplot
tkplot(gsd, layout=lsd,vertex.size=4,vertex.label=V(gsd)$labels,vertex.label.font=4,
       vertex.label.color=V(gsd)$label_color,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666")
LTKd=tkplot.getcoords(4)    #Manually correct and get the coordinates back

####Replot using corrected coordinates:
dev.off()
E(DENSE_g)$width<-Vscale(E(DENSE_g)$weight,1,10)
lsd<-LTKd
plot(DENSE_g, layout=lsd, vertex.size=10,rescale=FALSE, xlim=range(lsd[,1]),ylim=range(lsd[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=1)

V(gsd)$size<-Vscale(abs(V(gsd)$size), 1000, 1000 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

plot(gsd, layout=lsd,vertex.size=V(gsd)$size,vertex.label=V(gsd)$labels,vertex.label.dist=60,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gsd)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)

#Save PDF of dimensions ??x??
#####REPEAT CORRECTION UNTIL GOOD THEN SAVE THAT SHIT:
#save(LTKd,file="U2_LayoutTK_vsize1000_asp1.layout")




























####Plot subgraph consisiting of neighborhood of specific nodes
k=1   # consider k-th order neighbors
N_node1=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="21AS_37"))   ) )
N_node2=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="21AS_43"))   ) )
N_node3=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="H2O2_4h"))   ) )
N_node4=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_24h"))   ) )
N_node5=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_4h"))   ) )
N_node6=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_43"))   ) )
nodevector=c(N_node1,N_node2,N_node3,N_node4,N_node5,N_node6)
plotnodes=unique(nodevector)

#add repeated k-th order neighbors:
k=2   # consider k-th order neighbors
N_node1=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="21AS_37"))   ) )
N_node2=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="21AS_43"))   ) )
N_node3=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="H2O2_4h"))   ) )
N_node4=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_24h"))   ) )
N_node5=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_4h"))   ) )
N_node6=unlist(neighborhood(g, k ,nodes=c( which(V(g)$labels=="ctr_43"))   ) )
nodevector=c(N_node1,N_node2,N_node3,N_node4,N_node5,N_node6)
repeated=nodevector[which(duplicated(nodevector))]
plotnodes=unique(c(plotnodes,repeated))


dev.off()

##### ADD TRANSPARENCIES #####
dev.off()
E(g)$width<-Vscale(E(g)$weight,1,8)
V(g)$color <- colbar[comps]
V(g)$color<-paste(V(g)$color,"44",sep="")
E(g)$color<-ifelse(E(Tg)$weight<0,"#55885599","#BB111199")

subg=induced_subgraph(g,plotnodes)
l<-layout.fruchterman.reingold(subg)
plot(subg, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)
##########

subgs=induced_subgraph(gs,plotnodes)

plot(subgs, layout=l,vertex.size=V(subgs)$size/2,vertex.label=V(subgs)$labels,vertex.label.dist=1,vertex.label.cex=0.8,vertex.label.font=4,
     vertex.label.color=V(subgs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)
#E(gs)$width<-1.75














N_DFOA=unlist(neighborhood(g,1,nodes=c( which(V(g)$labels=="DFOA")-1)   ) )
N_Hemin=unlist(neighborhood(g,1,nodes=c( which(V(g)$labels=="Hemin")-1)   ))
N_CPX=unlist(neighborhood(g,1,nodes=c( which(V(g)$labels=="CPX")-1)   ))
IRONnodes=unique(c(N_DFOA,N_Hemin,N_CPX))

#LAB1=labels(CRD["Hemin",which( abs(CRD["Hemin",])+abs(CRD["DFOA",])>0.625 )])  #with 0.35 in glasso regularization
LAB1=labels(CRD["Hemin",which( abs(CRD["Hemin",])+abs(CRD["DFOA",]+CRD["CPX",])>0.85 )])  #with 0.35 in glasso regularization
LAB2=labels(CRD["Hemin",which( abs(CRD["Hemin",])>0.4 )])
LAB3=labels(CRD["Hemin",which( abs(CRD["DFOA",])> 0.4 )])
LAB4=labels(CRD["Hemin",which( abs(CRD["CPX",])> 0.4 )])

LAB=unique(c(LAB1,LAB2,LAB3,LAB4,"Hemin","DFOA"))
IRONnodes=match(LAB, V(g)$labels)-1
IRON_g=subgraph(g,IRONnodes)

dev.off()

E(IRON_g)$width<-Vscale(E(IRON_g)$weight,1,10)

l<-layout.fruchterman.reingold(IRON_g,params=list(area=vcount(IRON_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(IRON_g)^2.3)
plot(IRON_g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, IRONnodes)
Tgs <- subgraph(Tg, IRONnodes)
Tgs <- subgraph(sg, IRONnodes)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 40, 40 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=5,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)







LAB1=labels(CRD["Hemin",which( abs(CRD["Hemin",])+abs(CRD["DFOA",])+abs(CRD["CPX",])+abs(CRD["FTL",])+abs(CRD[,"ACO1"]) > 1.3 )])  #with 0.3 in glasso regularization
LAB1=LAB1[!LAB1=="FTL"]
LAB1=LAB1[!LAB1=="DFOA"]
LAB1=LAB1[!LAB1=="CPX"]
LAB1=LAB1[!LAB1=="Hemin"]
LAB1=LAB1[!LAB1=="ACO1"]
LAB2=labels(CRD["Hemin",which( abs(CRD["Hemin",])>0.41 )])
LAB3=labels(CRD["Hemin",which( abs(CRD["DFOA",])> 0.41 )])
LAB4=labels(CRD["Hemin",which( abs(CRD["CPX",])> 0.41 )])
LAB5=labels(CRD["Hemin",which( abs(CRD["FTL",])>0.41 )])
LAB6=labels(CRD["Hemin",which( abs(CRD["ACO1",])> 0.41 )])

LAB=unique(c(LAB1,LAB2,LAB3,LAB4,LAB5,LAB6))
LAB=LAB[!LAB==1]
if (length(grep("_\\d+h",LAB) )>0  )  {LAB=LAB[-c(grep("_\\d+h",LAB))]}
if (length(grep("TOP",LAB) )>0  )  {LAB=LAB[-c(grep("TOP",LAB))]}

LAB=unique(rownames(which( abs(CRD[,LAB])> 0.5,arr.ind=TRUE )))

IRONnodes=match(LAB, V(g)$labels)-1
IRONnodes=IRONnodes[!is.na(IRONnodes)]

IRON_g=subgraph(g,IRONnodes)

dev.off()

E(IRON_g)$width<-Vscale(E(IRON_g)$weight,1,10)

l<-layout.fruchterman.reingold(IRON_g,params=list(area=vcount(IRON_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(IRON_g)^2.3)
plot(IRON_g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, IRONnodes)
Tgs <- subgraph(Tg, IRONnodes)
Tgs <- subgraph(sg, IRONnodes)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 100, 100 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=15,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)







TOPO=c("CPT_3h","CPT_8h","CPT_24h","Etopo_3h","Etopo_8h","Etopo_24h","Doxor_3h","Doxor_8h","Doxor_24h","TOP2ci_3h","TOP2ci_8h","TOP2ci_24h")
#labels(CRD["Hemin",which( colSums(abs(CRD[TOPO,])) > 3.1 ) ]) 
LAB1=labels(CRD["Hemin",which( colSums(abs(CRD[TOPO,])) > 2.9 ) ])  #with 0.35 in glasso regularization
LAB1=LAB1[-c(which(LAB1%in%TOPO))]

LAB2=c("1");
for (i in 1:length(TOPO)){
LAB2=c(LAB2,   labels(CRD["Hemin",which( abs(CRD[TOPO[i],])>0.45 )]) )
}
LAB=unique(c(LAB1,LAB2))
LAB=LAB[!LAB==1]

LAB=unique(rownames(which( abs(CRD[,LAB])> 0.5,arr.ind=TRUE )))

IRONnodes=match(LAB, V(g)$labels)-1
IRONnodes=IRONnodes[!is.na(IRONnodes)]

IRON_g=subgraph(g,IRONnodes)

dev.off()

E(IRON_g)$width<-Vscale(E(IRON_g)$weight,1,10)

l<-layout.fruchterman.reingold(IRON_g,params=list(area=vcount(IRON_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(IRON_g)^2.3)
plot(IRON_g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, IRONnodes)
Tgs <- subgraph(Tg, IRONnodes)
Tgs <- subgraph(sg, IRONnodes)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 80, 80 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=10,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)




C["Spliceostatin",]=C["Spliceostatin",]+0.08
C[,"Spliceostatin"]=C[,"Spliceostatin"]+0.08
C["Meayamycin",]=C["Meayamycin",]-0.15
C[,"Meayamycin"]=C[,"Meayamycin"]-0.15
C["TG003",]=C["TG003",]+0.03
C[,"TG003"]=C[,"TG003"]+0.03
C["Cyclosporin",]=C["Cyclosporin",]+0.01
C[,"Cyclosporin"]=C[,"Cyclosporin"]+0.01
C["HNRPA2B1","SKIV2L2"]=C["HNRPA2B1","SKIV2L2"]+0.01
C["SKIV2L2","HNRPA2B1"]=C["SKIV2L2","HNRPA2B1"]+0.01

TOPO=V(g)$labels[grep("[Go][r0][0i][n3]",V(g)$labels)] #Cyclosporin/TG003
#TOPO=V(g)$labels[grep("[ai][tn][_i][Hn]",V(g)$labels)] #Meayamycin/Spliceostatin
#TOPO=V(g)$labels[grep("[OK][S6]_",V(g)$labels)] #OsmoticStress+MK6
#TOPO=V(g)$labels[grep("[ayeoG][tcr0][i0][n3]",V(g)$labels)] #DRUGS
#TOPO=V(g)$labels[grep("OS_\\d+",V(g)$labels)] #OsmoticStress
#TOPO=V(g)$labels[grep("[Da]_\\d+",V(g)$labels)] #ActD/Sta
#TOPO=c("CPT_3h","CPT_8h","CPT_24h","Etopo_3h","Etopo_8h","Etopo_24h","Doxor_3h","Doxor_8h","Doxor_24h","TOP2ci_3h","TOP2ci_8h","TOP2ci_24h","TOP2A","TOP1","TOPBP1")
#labels(CRD["Hemin",which( colSums(abs(CRD[TOPO,])) > 3.1 ) ]) 
LAB1=labels(CRD["SF3B4",which( colSums(abs(C[TOPO,])) > 1.4) ])  #3.5 for TOPO with 0.35 in glasso regularization, 1.2 for Sta/Act
LAB1=LAB1[-c(which(LAB1%in%TOPO))]

LAB2=c("1");
for (i in 1:length(TOPO)){
  LAB2=c(LAB2,   labels(C["SF3B4",which( abs(C[TOPO[i],])>0.48 )]) ) #0.425 for TOPO, 0.4 for Act/Sta, 0.75 for Meaya/Spliceost
}

LAB=unique(c(LAB1,LAB2,TOPO))
LAB=LAB[!LAB==1]

if (length(grep("OS_",LAB) )>0  )  {LAB=LAB[-c(grep("OS_",LAB))]}
if (length(grep("[Da]_\\d+",LAB) )>0  )  {LAB=LAB[-c(grep("[Da]_\\d+",LAB))]}
if (length(grep("TOP",LAB) )>0  )  {LAB=LAB[-c(grep("TOP",LAB))]}
#if (length(grep("[ay][tc]in",LAB) )>0  )  {LAB=LAB[-c(grep("[ay][tc]in",LAB))]}
#if (length(grep("in_H",LAB) )>0  )  {LAB=LAB[-c(grep("in_H",LAB))]}
if (length(grep("getin",LAB) )>0  )  {LAB=LAB[-c(grep("getin",LAB))]}
#if (length(grep("003",LAB) )>0  )  {LAB=LAB[-c(grep("003",LAB))]}
if (length(grep("_LV",LAB) )>0  )  {LAB=LAB[-c(grep("_LV",LAB))]}

LAB=unique(rownames(which( abs(C[,LAB])> 0.48,arr.ind=TRUE )))  #0.5/0.48 Normally, 0.75 for Meayam/Spliceostatin

if (length(grep("OS_",LAB) )>0  )  {LAB=LAB[-c(grep("OS_",LAB))]}
if (length(grep("[Da]_\\d+",LAB) )>0  )  {LAB=LAB[-c(grep("[Da]_\\d+",LAB))]}
if (length(grep("TOP",LAB) )>0  )  {LAB=LAB[-c(grep("TOP",LAB))]}
#if (length(grep("[ay][tc]in",LAB) )>0  )  {LAB=LAB[-c(grep("[ay][tc]in",LAB))]}
if (length(grep("in_L",LAB) )>0  )  {LAB=LAB[-c(grep("in_L",LAB))]}
if (length(grep("getin",LAB) )>0  )  {LAB=LAB[-c(grep("getin",LAB))]}
#if (length(grep("003",LAB) )>0  )  {LAB=LAB[-c(grep("003",LAB))]}
if (length(grep("Sudemycin",LAB) )>0  )  {LAB=LAB[-c(grep("Sudemycin",LAB))]}
if (length(grep("_LV",LAB) )>0  )  {LAB=LAB[-c(grep("_LV",LAB))]}



IRONnodes=match(LAB, V(g)$labels)-1
IRONnodes=IRONnodes[!is.na(IRONnodes)]

IRON_g=subgraph(g,IRONnodes)

dev.off()

E(IRON_g)$width<-Vscale(E(IRON_g)$weight,1,10)

l<-layout.fruchterman.reingold(IRON_g,params=list(area=vcount(IRON_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(IRON_g)^2.3)
plot(IRON_g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, IRONnodes)
Tgs <- subgraph(Tg, IRONnodes)
Tgs <- subgraph(sg, IRONnodes)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75


ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 80, 80 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
V(gs)$color=c(rep(brewer.pal(7,"Accent")[1],2),brewer.pal(7,"Set1")[3],rep(brewer.pal(7,"Set2")[7],2) ) #TG003, Cyclosporin Colors
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=10,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)


###Manually correct with tkplot
tkplot(gs, layout=ls,vertex.size=4,vertex.label=V(gs)$labels,vertex.label.font=4,
       vertex.label.color=V(gs)$label_color,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666")
LTKd=tkplot.getcoords(2)    #Manually correct and get the coordinates back

####### REPLOT USING tkplot manual correction #######
dev.off()
lsd=LTKd
plot(IRON_g, layout=lsd, vertex.size=100,rescale=FALSE, xlim=range(lsd[,1]),ylim=range(lsd[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=1)
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 80, 80 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$size<-V(gs)$size*30
#V(gs)$color=c(rep(brewer.pal(8,"Dark2")[2],12),rep(brewer.pal(7,"Set2")[7],2) ) #Spliceostatin, Meayamycin Colors
plot(gs, layout=lsd,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=120,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=1)






FACTOR="ZNF326"
LABZ=labels(CRD[FACTOR,which( abs(CRD[FACTOR,])>0.4 )])
IRONnodes=match(LABZ, V(g)$labels)-1

LABZ=c("1",LABZ);
LABZ=LABZ[!LABZ==1]
LABZ=unique(rownames(which( abs(CRD[,LABZ])> 0.5,arr.ind=TRUE )))
IRONnodes=match(LABZ, V(g)$labels)-1

IRONnodes=IRONnodes[!is.na(IRONnodes)]
IRON_g=subgraph(g,IRONnodes)
dev.off()

E(IRON_g)$width<-Vscale(E(IRON_g)$weight,1,10)

l<-layout.fruchterman.reingold(IRON_g,params=list(area=vcount(IRON_g)^2.3),niter=5000,coolexp=3,repulserad=vcount(IRON_g)^2.3)
plot(IRON_g, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, IRONnodes)
Tgs <- subgraph(Tg, IRONnodes)
Tgs <- subgraph(sg, IRONnodes)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 20, 20 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=5,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)















#####################################################################
GG=GRAO
#GG=g
#### Degree Distribution ####
plot(degree.distribution(GG), xlab="node degree")
lines(degree.distribution(GG))


CORE=rownames(DATA)[which(DATA[,"CLASS1"]=="Core Spl.")]
CORE=c("1",CORE);
CORE=CORE[!CORE==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
CORE_n=match(CORE, V(GG)$labels)-1
CORE_n=CORE_n[!is.na(CORE_n)]
CORE_g=subgraph(GG,CORE_n)

SF=rownames(DATA)[which(DATA[,"CLASS1"]=="Spl. Factors")]
SF=c("1",SF);
SF=SF[!SF==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
SF_n=match(SF, V(GG)$labels)-1
SF_n=SF_n[!is.na(SF_n)]
SF_g=subgraph(GG,SF_n)

CHROM=rownames(DATA)[which(DATA[,"CLASS1"]=="Chrom. Factors")]
CHROM=c("1",CHROM);
CHROM=CHROM[!CHROM==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
CHROM_n=match(CHROM, V(GG)$labels)-1
CHROM_n=CHROM_n[!is.na(CHROM_n)]
CHROM_g=subgraph(GG,CHROM_n)

OTHER=rownames(DATA)[which(DATA[,"CLASS1"]=="Other")]
OTHER=c("1",OTHER);
OTHER=OTHER[!OTHER==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
OTHER_n=match(OTHER, V(GG)$labels)-1
OTHER_n=OTHER_n[!is.na(OTHER_n)]
OTHER_g=subgraph(GG,OTHER_n)




PERS=rownames(DATA)[which(DATA[,"CLASS2"]=="Persist. Spl.")]
PERS=c("1",PERS);
PERS=PERS[!PERS==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
PERS_n=match(PERS, V(GG)$labels)-1
PERS_n=PERS_n[!is.na(PERS_n)]
PERS_g=subgraph(GG,PERS_n)

TMID=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Mid Spl.")]
TMID=c("1",TMID);
TMID=TMID[!TMID==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
TMID_n=match(TMID, V(GG)$labels)-1
TMID_n=TMID_n[!is.na(TMID_n)]
TMID_g=subgraph(GG,TMID_n)

TEARLY=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Early Spl.")]
TEARLY=c("1",TEARLY);
TEARLY=TEARLY[!TEARLY==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
TEARLY_n=match(TEARLY, V(GG)$labels)-1
TEARLY_n=TEARLY_n[!is.na(TEARLY_n)]
TEARLY_g=subgraph(GG,TEARLY_n)

TLATE=rownames(DATA)[which(DATA[,"CLASS2"]=="Trans. Late Spl.")]
TLATE=c("1",TLATE);
TLATE=TLATE[!TLATE==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
TLATE_n=match(TLATE, V(GG)$labels)-1
TLATE_n=TLATE_n[!is.na(TLATE_n)]
TLATE_g=subgraph(GG,TLATE_n)

EJC=rownames(DATA)[which(DATA[,"CLASS2"]=="EJC")]
EJC=c("1",EJC);
EJC=EJC[!EJC==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
EJC_n=match(EJC, V(GG)$labels)-1
EJC_n=EJC_n[!is.na(EJC_n)]
EJC_g=subgraph(GG,EJC_n)


SR=rownames(DATA)[which(DATA[,"CLASS2"]=="SR proteins")]
SR=c("1",SR);
SR=SR[!SR==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
SR_n=match(SR, V(GG)$labels)-1
SR_n=SR_n[!is.na(SR_n)]
SR_g=subgraph(GG,SR_n)

hnRNP=rownames(DATA)[which(DATA[,"CLASS2"]=="hnRNP proteins")]
hnRNP=c("1",hnRNP);
hnRNP=hnRNP[!hnRNP==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
hnRNP_n=match(hnRNP, V(GG)$labels)-1
hnRNP_n=hnRNP_n[!is.na(hnRNP_n)]
hnRNP_g=subgraph(GG,hnRNP_n)

dev.off()
#plot(degree.distribution(CORE_g), xlab="node degree")


plot(prop.table(table(degree(GG)))*100,type="l",col="black",ylim=c(0,80))
points(prop.table(table(degree(GG)))*100,type="o",col="black")
#lines(prop.table(table(degree(GG,CORE_n)))*100,type="l",col="red")
points(prop.table(table(degree(GG,CORE_n)))*100,type="o",col="red")
points(prop.table(table(degree(GG,SF_n)))*100,type="o",col="green")
points(prop.table(table(degree(GG,CHROM_n)))*100,type="o",col="blue")




dev.off()
LCOL=c("black","darkred","darkgreen","darkorange","darkblue")
plot(c(0,as.numeric(names(table(degree(GG))))),c(0, cumsum(  prop.table(table(degree(GG)) ) *100    ) )     ,type="l",col=LCOL[1],lty=2,ylim=c(0,100),xlab="Degree",ylab="Cumulative Frequency")
points(c(0,as.numeric(names(table(degree(GG))))),c(0, cumsum(  prop.table(table(degree(GG)) ) *100    ) )   ,type="o",lty="dashed",pch=20,lwd=2,col=LCOL[1])

points(c(0,as.numeric(names(table(degree(GG,CORE_n))))), c(0, cumsum(  prop.table(table(degree(GG,CORE_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[2])

points(c(0,as.numeric(names(table(degree(GG,SF_n))))), c(0, cumsum(  prop.table(table(degree(GG,SF_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[3])

points(c(0,as.numeric(names(table(degree(GG,OTHER_n))))), c(0, cumsum(  prop.table(table(degree(GG,OTHER_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[4])

points(c(0,as.numeric(names(table(degree(GG,CHROM_n))))), c(0, cumsum(  prop.table(table(degree(GG,CHROM_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[5])


legend("bottomright", inset=.1,c("All","Spliceosome","Splic. Factors","Other RNApr.","Chrom. Factors"), col=c(LCOL[1],LCOL[2],LCOL[3],LCOL[4],LCOL[5]) , lwd=c(3,3) , lty=c("dashed","solid","solid","solid","solid")  , horiz=FALSE,bty="n",y.intersp=2,seg.len=3)
#Save PDF at 10:9




dev.off()
LCOL=c("black","darkred","darkgreen","darkblue","darkorange","purple")
plot(c(0,as.numeric(names(table(degree(GG,CORE_n))))),c(0, cumsum(  prop.table(table(degree(GG,CORE_n)) ) *100    ) )     ,type="l",col=LCOL[1],lty=2,ylim=c(0,100),xlab="Degree",ylab="Cumulative Frequency")
points(c(0,as.numeric(names(table(degree(GG,CORE_n))))),c(0, cumsum(  prop.table(table(degree(GG,CORE_n)) ) *100    ) )   ,type="o",lty="dashed",pch=20,lwd=2,col=LCOL[1])

points(c(0,as.numeric(names(table(degree(GG,PERS_n))))), c(0, cumsum(  prop.table(table(degree(GG,PERS_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[2])

points(c(0,as.numeric(names(table(degree(GG,TEARLY_n))))), c(0, cumsum(  prop.table(table(degree(GG,TEARLY_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[3])

points(c(0,as.numeric(names(table(degree(GG,TMID_n))))), c(0, cumsum(  prop.table(table(degree(GG,TMID_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[4])

points(c(0,as.numeric(names(table(degree(GG,TLATE_n))))), c(0, cumsum(  prop.table(table(degree(GG,TLATE_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[5])

#points(c(0,as.numeric(names(table(degree(GG,EJC_n))))), c(0, cumsum(  prop.table(table(degree(GG,EJC_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[6])


legend("bottomright", inset=.1,c("Spliceosome","Persistent","Trans. Early","Trans. Mid","Trans. Late","EJC"), col=c(LCOL[1],LCOL[2],LCOL[3],LCOL[4],LCOL[5],LCOL[6]) , lwd=c(3,3) , lty=c("dashed","solid","solid","solid","solid","solid")  , horiz=FALSE,bty="n",y.intersp=2,seg.len=3)
#Save PDF at 10:9


dev.off()
LCOL=c("black","darkred","darkgreen","darkblue","darkorange","purple")
plot(c(0,as.numeric(names(table(degree(GG))))),c(0, cumsum(  prop.table(table(degree(GG)) ) *100    ) )     ,type="l",col=LCOL[1],lty=2,ylim=c(0,100),xlab="Degree",ylab="Cumulative Frequency")
points(c(0,as.numeric(names(table(degree(GG))))),c(0, cumsum(  prop.table(table(degree(GG)) ) *100    ) )   ,type="o",lty="dashed",pch=20,lwd=2,col=LCOL[1])

points(c(0,as.numeric(names(table(degree(GG,SF_n))))), c(0, cumsum(  prop.table(table(degree(GG,SF_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[2])

points(c(0,as.numeric(names(table(degree(GG,SR_n))))), c(0, cumsum(  prop.table(table(degree(GG,SR_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[3])

points(c(0,as.numeric(names(table(degree(GG,hnRNP_n))))), c(0, cumsum(  prop.table(table(degree(GG,hnRNP_n)) ) *100    ) )   ,type="o",pch=20,lwd=2,col=LCOL[4])

legend("bottomright", inset=.1,c("All","SFactors","SR proteins","hnRNP proteins"), col=c(LCOL[1],LCOL[2],LCOL[3],LCOL[4]) , lwd=c(3,3) , lty=c("dashed","solid","solid","solid")  , horiz=FALSE,bty="n",y.intersp=2,seg.len=3)
















LABs=Sm;
LABs=rownames(DATA)[which(DATA[,"CLASS1"]=="Core Spl.")]
LABs=c("1",LABs);
LABs=LABs[!LABs==1]
#LABs=unique(rownames(which( abs(CRD[,LABs])> 0.5,arr.ind=TRUE )))
NODESs=match(LABs, V(g)$labels)-1

NODESs=NODESs[!is.na(NODESs)]
SUBg=subgraph(g,NODESs)
dev.off()



E(SUBg)$width<-Vscale(E(SUBg)$weight,1,10)

l<-layout.fruchterman.reingold(SUBg,params=list(area=vcount(SUBg)^2.3),niter=5000,coolexp=3,repulserad=vcount(SUBg)^2.3)
plot(SUBg, layout=l, vertex.size=10,rescale=FALSE, xlim=range(l[,1]),ylim=range(l[,2]),vertex.frame.color=paste(colbar[comps+1],"77",sep=""),vertex.label=NA,asp=0.72)

gs <- subgraph(sg, NODESs)
Tgs <- subgraph(Tg, NODESs)
Tgs <- subgraph(sg, NODESs)
E(gs)$color<-ifelse(E(Tgs)$weight<0,"#558855","#BB1111")
E(gs)$color<-"#222244"
E(gs)$width<-1.75

ScZ=(t(Z))[,get.vertex.attribute(gs,"labels")]
V(gs)$size<-apply(abs(ScZ),2,median,na.rm=TRUE)  
V(gs)$size<-Vscale(abs(V(gs)$size), 20, 20 * sqrt( max(V(gs)$size)/min(V(gs)$size)  )     )

V(gs)$label_color="#222222"
CHRnodes=which(CHR[V(gs)$labels,]==1)
V(gs)$label_color[CHRnodes]="DarkRed"
ls<-l
plot(gs, layout=ls,vertex.size=V(gs)$size,vertex.label=V(gs)$labels,vertex.label.dist=2,vertex.label.cex=1,vertex.label.font=4,
     vertex.label.color=V(gs)$label_color,add=TRUE,rescale=FALSE,edge.lty="dotted",vertex.frame.color="#666666",asp=0.72)




TMID_n=match(TMID, V(GG)$labels)-1
TMID_n=TMID_n[!is.na(TMID_n)]
TMID_g=subgraph(GG,TMID_n)


NLIST=list(PERS,TEARLY,TMID,TLATE)
#DRAW SIMPLIFIED NETWORK TOPOLOGY:
#First find number of connections between modules:
modules=LMods
ModAdj=matrix(data=0,nrow=4,ncol=4)

for(i in 1:length(NLIST)){
  for(j in 1:length(NLIST)){
NodesA=match(NLIST[[i]], V(GG)$labels)-1
NodesA=NodesA[!is.na(NodesA)]
NodesB=match(NLIST[[j]], V(GG)$labels)-1
NodesB=NodesB[!is.na(NodesB)]
NodeUN=unique(c(NodesA,NodesB))
SUBA=subgraph(GG,NodesA)
SUBB=subgraph(GG,NodesB)
SUBC=subgraph(GG,NodeUN)
ModAdj[i,j]=abs(length(E(SUBC))-length(E(SUBA))-length(E(SUBB)))
NormF=length(V(SUBA))*length(V(SUBB))
if (i==j) {NormF=0.5*length(V(SUBA))*(length(V(SUBA))-1) }
print ( c(i,j,NormF,ModAdj[i,j]) )
ModAdj[i,j]=ModAdj[i,j]/NormF  #Normalize over product of community sizes 
  }
}
NModAdj=100*ModAdj 
NModAdj=NModAdj/max(max(NModAdj)) #Normalize over max
Ncsize=sqrt(sapply(NLIST,length)/max(sapply(NLIST,length))) #Normalized community size


CGRAPH<-graph.adjacency(ModAdj,mode=c("max"),weighted=TRUE,diag=TRUE)
V(CGRAPH)$shape<-"circle"
V(CGRAPH)$color <- c("darkred","darkgreen","darkblue","darkorange")
V(CGRAPH)$size  <-Ncsize*40
E(CGRAPH)$labels=round(ModAdj[upper.tri(ModAdj,diag=TRUE)][c(1,2,4,7,3,5,8,6,9,10)],digits=2)

NModAdj[which(NModAdj<0.02)]=0
for(i in 0:(length(NLIST)-1)){
  for(j in 0:(length(NLIST)-1)){
    if (ModAdj[i+1,j+1]>0){
      E(CGRAPH,c(i,j))$width <- NModAdj[i+1,j+1]*30
    }
  }
}

V(CGRAPH)$label=c("P","E","M","L")
sl<-layout.circle(CGRAPH,params=list(area=40000))

plot(CGRAPH, layout=sl,vertex.label.cex=2.5,vertex.label.color="white",vertex.label.font=2,edge.label=E(CGRAPH)$labels,edge.label.cex=1.5,edge.label.color="black")












