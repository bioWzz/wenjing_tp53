source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicFeatures","INSPEcT","Rsamtools","deSolve","BiocParallel","compiler","sqldf"))
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

source("http://bioconductor.org/biocLite.R")
biocLite("diptest")
library(TCGAbiolinks)
library(diptest)


setwd("D:/ivan/ivan6-cell/Results/结果数据")
load("D:/2016_11_24/labAandC.RData")
ls()
library(ggplot2)

###一、red
p <- ggplot()
layer1 <- geom_density(data=data,colour = "white",size = 1,aes(x = x, y = ..density..))
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer2 <- geom_density(data=data1,colour = "white",linetype = 2,size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "red",size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "red",linetype = 2,size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

##blue
p <- ggplot()
layer1 <- geom_density(data=data,colour = "blue",size = 1,aes(x = x, y = ..density..))
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer2 <- geom_density(data=data1,colour = "blue",linetype = 2,size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer2 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer3 <- geom_density(data=data3,colour = "white",size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer3 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=data4,colour = "white",linetype = 2,size = 1,aes(x = x, y = ..density..))
(p1 <- p1 + layer4 + xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


###二、eRNA
con <- file("D:/2016_11_24/2_FigureS6-文静/eRNA_A549.fasta", "r")
line=readLines(con,n=1)
while( length(line) != 0 ) {
  #no <- strsplit(line, '_')
  #print(strsplit(line, '_')[[1]][2])
  mo<-strsplit(strsplit(line, '_')[[1]][2],",")
  #print(as.vector(mo))
  write.table(as.vector(mo), file = "D:/2016_11_24/2_FigureS6-文静/eRNA1.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  
     line=readLines(con,n=1)
     line=readLines(con,n=1)
}
close(con)

f <- read.csv("D:/2016_11_24/2_FigureS6-文静/eRNA1.txt")
class(f)
index<-which(f[,1]=="NA")
rm(f[index,1])
f<-na.omit(f)
unique(f)
write.table(unique(f), file = "D:/2016_11_24/2_FigureS6-文静/eRNA.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE)


##pRNA
con <- read.csv("D:/2016_11_24/2_FigureS6-文静/pRNA_GENE.csv")
class(con)
dim(con)
head(con)
index<-which(con[,2]=="pRNA")
c<-con[index,1]
write.table(c, file = "D:/2016_11_24/2_FigureS6-文静/pRNA.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE) 

##lincRNA\protein_coding
f <- read.csv("D:/2016_11_24/2_FigureS6-文静/gentype.txt",sep = "\t")
class(f)
index<-which(f[,2]=="lincRNA")
write.table(f[index,1], file = "D:/2016_11_24/2_FigureS6-文静/lincRNA.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE) 
index1<-which(f[,2]=="protein_coding")
write.table(f[index1,1], file = "D:/2016_11_24/2_FigureS6-文静/protein_coding.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE) 

##画图前的数据处理
load("D:/2016_11_24/labAandC.RData")
#ls()
#class(GeneLabAREP2)
#head(GeneLabAREP2)
GeneLabCREP1=cbind(ENSG=row.names(GeneLabCREP1),GeneLabCREP1)
row.names(GeneLabCREP1)=NULL  #把row.names变为第一列ENSG

f <- file("D:/2016_11_24/2_FigureS6-文静/eRNA.txt","r")
line=readLines(f,n=1)
while( length(line) != 0 ) {
  result=subset(GeneLabCREP1,ENSG==line)
  write.table(result, file = "D:/2016_11_24/GeneLabCREP1/eRNA.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  line=readLines(f,n=1)
}
close(f)

f1<- file("D:/2016_11_24/2_FigureS6-文静/lincRNA.txt","r")
line=readLines(f1,n=1)
while( length(line) != 0 ) {
  result=subset(GeneLabCREP1,ENSG==line)
  write.table(result, file = "D:/2016_11_24/GeneLabCREP1/lincRNA.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  line=readLines(f1,n=1)
}
close(f1) 

f2 <- file("D:/2016_11_24/2_FigureS6-文静/pRNA.txt","r")
line=readLines(f2,n=1)
while( length(line) != 0 ) {
  result=subset(GeneLabCREP1,ENSG==line)
  write.table(result, file = "D:/2016_11_24/GeneLabCREP1/pRNA.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  line=readLines(f2,n=1)
}
close(f2)

f3 <- file("D:/2016_11_24/2_FigureS6-文静/protein_coding.txt","r")
line=readLines(f3,n=1)
while( length(line) != 0 ) {
  result=subset(GeneLabCREP1,ENSG==line)
  write.table(result, file = "D:/2016_11_24/GeneLabCREP1/protein_coding.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  line=readLines(f3,n=1)
}
close(f3)


###三、画total\4su图
library(ggplot2)
Express<- read.csv("D:/2016_11_24/GeneLabCREP1/eRNA.txt",sep = "\t",header=F)
# 数据进行初筛
Express[apply(Express[,2:4]<=1, FUN = any, 1), ] = NA
Express[apply(Express[,5:7]<=0, FUN = any, 1), ] = NA
Express[apply(Express[,8:16]<=0.01, FUN = any, 1), ] = NA
Express=na.omit(Express)

# total
total=Express[,2:4]
total=as.data.frame(as.numeric(unlist(total)))
t=median(as.numeric(unlist(total)))
total=log2(total)

# 4su
su=Express[,5:7]
su=as.data.frame(as.numeric(unlist(su)))
SU=median(as.numeric(unlist(su)))
su=log2(su)

##画概率分布图 
graphics.off()
p <- ggplot()
layer1 <- geom_density(data=total,colour = "red",size = 1,aes(x = total, y = ..density..))
(p1 <- p + layer1+xlim(-10,10))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer2 <- geom_density(data=su,colour = "red",linetype = 2,size = 1,aes(x = su, y = ..density..))
(p1 <- p1 + layer2 +xlim(-10,10)+ xlab("log(Expression)") + ylab("fraction of junctions"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))


###四、画3率
f <- read.csv("D:/2016_11_24/GeneLabAREP1/eRNA.txt",sep = "\t",header=F)
f1 <- read.csv("D:/2016_11_24/GeneLabAREP1/lincRNA.txt",sep = "\t",header=F)
f2 <- read.csv("D:/2016_11_24/GeneLabAREP1/pRNA.txt",sep = "\t",header=F)
f3 <- read.csv("D:/2016_11_24/GeneLabAREP1/protein_coding.txt",sep = "\t",header=F)

# 数据进行初筛
f[apply(f[,2:4]<=1, FUN = any, 1), ] = NA
f[apply(f[,5:7]<=0, FUN = any, 1), ] = NA
f[apply(f[,8:16]<=0.01, FUN = any, 1), ] = NA
f=na.omit(f)

f1[apply(f1[,2:4]<=1, FUN = any, 1), ] = NA
f1[apply(f1[,5:7]<=0, FUN = any, 1), ] = NA
f1[apply(f1[,8:16]<=0.01, FUN = any, 1), ] = NA
f1=na.omit(f1)

f2[apply(f2[,2:4]<=1, FUN = any, 1), ] = NA
f2[apply(f2[,5:7]<=0, FUN = any, 1), ] = NA
f2[apply(f2[,8:16]<=0.01, FUN = any, 1), ] = NA
f2=na.omit(f2)

f3[apply(f3[,2:4]<=1, FUN = any, 1), ] = NA
f3[apply(f3[,5:7]<=0, FUN = any, 1), ] = NA
f3[apply(f3[,8:16]<=0.01, FUN = any, 1), ] = NA
f3=na.omit(f3)

# synthsis
synthsisf=f[,8:10]
synthsisf=as.data.frame(as.numeric(unlist(synthsisf)))
sf=median(as.numeric(unlist(synthsisf)))
synthsisf=log2(synthsisf)

synthsisf1=f1[,8:10]
synthsisf1=as.data.frame(as.numeric(unlist(synthsisf1)))
sf1=median(as.numeric(unlist(synthsisf1)))
synthsisf1=log2(synthsisf1)

synthsisf2=f2[,8:10]
synthsisf2=as.data.frame(as.numeric(unlist(synthsisf2)))
sf2=median(as.numeric(unlist(synthsisf2)))
synthsisf2=log2(synthsisf2)

synthsisf3=f3[,8:10]
synthsisf3=as.data.frame(as.numeric(unlist(synthsisf3)))
sf3=median(as.numeric(unlist(synthsisf3)))
synthsisf3=log2(synthsisf3)


graphics.off()
library(ggplot2)
p <- ggplot()
layer1 <- geom_density(data=synthsisf,colour = "lightgreen",size = 1,aes(x = synthsisf, y = ..density..))
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer2 <- geom_density(data=synthsisf1,colour = "purple",size = 1,aes(x = synthsisf1, y = ..density..))
(p1 <- p1 + layer2 )
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

#layer3 <- geom_density(data=synthsisf2,colour = "yellow",size = 1,aes(x = synthsisf2, y = ..density..))
#(p1 <- p1 + layer3)
#p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer4 <- geom_density(data=synthsisf3,colour = "black",size = 1,aes(x = synthsisf3, y = ..density..))
(p1 <- p1 + layer4 + xlab("Gene's transcription rate") + ylab("fraction of genes"))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))



p2=ggplot()
layer2 <- geom_density(data=processing,colour = "lightgreen",size = 1,aes(x = processing,y = ..density..))+
(p<-p2+layer2+labs(x = "Gene's processing rate ", y = "fraction of genes", title = ""))
p+theme(panel.grid =element_blank(),panel.background = element_blank())  ## 删去网格线

p3=ggplot()
layer3 <- geom_density(data=degration,colour = "lightgreen",size = 1, aes(x = degration,y = ..density..)) 
(p<-p3+layer3+labs(x = "Gene's degration rate ", y = "fraction of genes", title = ""))
p+theme(panel.grid =element_blank(),panel.background = element_blank())  ## 删去网格线

# processing
processing=Express[,14:16]
processing=as.data.frame(as.numeric(unlist(processing)))
p=median(as.numeric(unlist(processing)))
processing=log2(processing)

# degration
degration=Express[,11:13]
degration=as.data.frame(as.numeric(unlist(degration)))
d=median(as.numeric(unlist(degration)))
degration=log2(degration)
 


###五、以行、列为单位画箱式图

graphics.off()
c<-colplot(f)
r<-rowplot(f,1)

colplot <- function(dataframe){
  #dataframe <- subset(dataframe, select = -class)
  dataframe<-dataframe[,-dim(dataframe)[2]]
  return(boxplot(dataframe,las=1))
}

rowplot <- function(dataframe,classnum){
  #index<-which(dataframe[,"class"]==classnum)
  #dataframe <- subset(dataframe[index,], select = -class )
  index<-which(dataframe[,dim(dataframe)[2]]==classnum)
  dataframe<-dataframe[index,][,-dim(dataframe)[2]]
  return(boxplot(dataframe,las=1))
}

##以行为单位画箱式图
index<-which(f[,"class"]==1)
f <- subset(f[index,], select = -class )
boxplot(f,las=1)

##以列为单位画箱式图
f <- subset(f, select = -class ) #删数据框最后一列（分类列）
boxplot(f,las=1)  #画全部列
boxplot(f["total_0"],las=1)+axis(1, at=1,lab="total_0")#画指定一列

#画指定列>2
library(ggplot2)
p <- ggplot()
layer1 <- geom_boxplot(data=f["total_0"],colour = "black",size = 1,aes(x="total_0",y =total_0))
(p1 <- p + layer1)
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

layer2 <- geom_boxplot(data=f["total_4"],colour = "black",size = 1,aes(x ="total_4", y =total_4))
(p1 <- p1 + layer2 )
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))



###六、画a散点图图gene-TX 3率+拟合直线
load("D:/2016_11_24/labAandC.RData")

a1<-a_gene_TX_process(GeneLabAREP1,TXLabAREP1,13)
##输入第一个是基因的数据框，第二个是转录本的数据框，第三个是i
##i=7是数据框中synthsis所在的第一个列数,输出为gene_TX合成的synthsis数据框，根据此数据框画图
##i=10是数据框中processing所在的第一个列数,输出为gene_TX合成的processing数据框，根据此数据框画图
##i=13是数据框中degration所在的第一个列数,输出为gene_TX合成的degration数据框，根据此数据框画图
a_gene_TX_process <- function(dataframe1,dataframe2,i){
  f <- read.table("D:/2016_11_24/散点图/mart_export.txt",sep = "\t",header=T)##基因和转录本的对应关系文件

  f1=cbind(Gene.ID=row.names(dataframe1),dataframe1[,i])
  row.names(f1)=NULL  
  f2=cbind(Transcript.ID=row.names(dataframe2),dataframe2[,i])
  row.names(f2)=NULL   
  total1 <- merge(f,f1,by="Gene.ID",all=T)
  total2 <- merge(total1,f2,by="Transcript.ID",all=T)
  total2=na.omit(total2)

  f3=cbind(Gene.ID=row.names(dataframe1),dataframe1[,i+1])
  row.names(f3)=NULL  
  f4=cbind(Transcript.ID=row.names(dataframe2),dataframe2[,i+1])
  row.names(f4)=NULL   
  total3 <- merge(f,f3,by="Gene.ID",all=T)
  total4 <- merge(total3,f4,by="Transcript.ID",all=T)
  total4=na.omit(total4)

  f5=cbind(Gene.ID=row.names(dataframe1),dataframe1[,i+2])
  row.names(f5)=NULL  
  f6=cbind(Transcript.ID=row.names(dataframe2),dataframe2[,i+2])
  row.names(f6)=NULL   
  total5 <- merge(f,f5,by="Gene.ID",all=T)
  total6 <- merge(total5,f6,by="Transcript.ID",all=T)
  total6=na.omit(total6)

  total<-rbind(total2,total4)
  data<-rbind(total,total6)
  write.table(data, file = paste("D:/2016_11_24/散点图/gene_TX",i, ".txt", sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
}


a2<-a_gene_TX_plot(13,"black")
##根据gene_TX合成的synthsis数据框画图，输入为i，线的颜色
a_gene_TX_plot <- function(i,linecolour){
  file = paste("D:/2016_11_24/散点图/gene_TX",i, ".txt", sep="")
  data  <- read.table(file,sep = "\t")
  r<-cor(data[,3],data[,4])
  ##lab:i=7,production rate\i=10,processing rate\i=13,degration rate
  graphics.off()
  z<-lm(data[,4]~data[,3])
  plot(data[,3:4],pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", r),xlab="gene production rate",ylab="TX production rate")
  #plot(data[,3:4],pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", r),xlab="gene processing rate",ylab="TX processing rate")
  #plot(data[,3:4],pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", r),xlab="gene degration rate",ylab="TX degration rate")
  lines(data[,3],fitted(z),lty=2,col=linecolour)
}


###六、画b散点图gene交叉3率+拟合直线
load("D:/2016_11_24/labAandC.RData")
#dataframe<-GeneLabAREP1

b<-b_plot(GeneLabAREP1,"black")

b_plot<- function(dataframe,linecolour){
  dataframe=na.omit(dataframe[,7:15])#7-15是数据框synthsis到degration的下标
  
  synthsis=dataframe[,1:3]
  synthsis=as.data.frame(as.numeric(unlist(synthsis)))

  processing=dataframe[,4:6]
  processing=as.data.frame(as.numeric(unlist(processing)))

  degration=dataframe[,7:9]
  degration=as.data.frame(as.numeric(unlist(degration)))

  x<-synthsis[,1]
  y<-processing[,1]
  z<-degration[,1]

  graphics.off()
  par(mfrow=c(1,3))
  m1<-lm(y~x)
  plot(x,y,pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", cor(x,y)),xlab="gene production rate",ylab="gene processing rate")
  lines(x,fitted(m1),lty=2,col=linecolour)

  m2<-lm(z~x)
  plot(x,z,pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", cor(x,z)),xlab="gene production rate",ylab="gene degration rate")
  lines(x,fitted(m2),lty=2,col=linecolour)

  m3<-lm(z~y)
  plot(y,z,pch=20,col= "darkblue",cex = 1,las=1,main =paste("r=", cor(y,z)),xlab="gene processing rate",ylab="gene degration rate")
  lines(y,fitted(m3),lty=2,col=linecolour)
}

####################################
#用ggplot2画散点图（同上）
graphics.off()
library(ggplot2)
p <- ggplot()
layer1 <-geom_point(data=data[,3:4],colour = "darkblue",size = 1.5,aes(x = data[,3], y = data[,4]))
(p1 <- p + layer1+labs(x = "gene production rate", y = "TX production rate", title = paste("r=", r))
+geom_smooth(aes(x = data[,3], y = data[,4]), method="lm", fullrange=TRUE,colour ="black",linetype = 2))
p1+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
####################################

###七、用DEseq2做基因差异表达
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

normal_exp <- read.table("D:/2016_11_24/COAD_readCount_normal.txt",sep = "\t",header=T)
tumor_exp <- read.table("D:/2016_11_24/COAD_readCount_tumor.txt",sep = "\t",header=T)


normal_exp <- subset(normal_exp, select = -GeneID.ID.WioutCli.Rep)
row.names(normal_exp)<-normal_exp[,1]
normal_exp <- subset(normal_exp, select = -Symbol.ID.WioutCli.Rep )

tumor_exp <- subset(tumor_exp, select = -GeneID.ID.WioutCli.Rep )
tumor_exp <- subset(tumor_exp, select = -Symbol.ID.WioutCli.Rep )

colnames(normal_exp) <- paste0(colnames(normal_exp), "_normal")
colnames(tumor_exp) <- paste(colnames(tumor_exp), 1:ncol(tumor_exp), sep="%")
countData <- cbind(normal_exp, tumor_exp)
countData <- round(countData)
pData <- data.frame(phenotype=factor(c(rep("normal", ncol(normal_exp)), 
                                       rep("tumor", ncol(tumor_exp))), levels=c("normal", "tumor")))
rownames(pData) <- colnames(countData)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = pData, design = ~ phenotype)
dds <- DESeq(dds)
dds <- replaceOutliersWithTrimmedMean(dds)
res <- results(dds, cooksCutoff=FALSE)
res <- res[order(res$padj),]
as.data.frame(res)

DiffExp(normal_exp, tumor_exp)


DiffExp <- function(normal_exp, tumor_exp) {
  library(DESeq2)
  colnames(normal_exp) <- paste0(colnames(normal_exp), "_normal")
  colnames(tumor_exp) <- paste(colnames(tumor_exp), 1:ncol(tumor_exp), sep="%")
  countData <- cbind(normal_exp, tumor_exp)
  countData <- round(countData)
  pData <- data.frame(phenotype=factor(c(rep("normal", ncol(normal_exp)), 
                                         rep("tumor", ncol(tumor_exp))), levels=c("normal", "tumor")))
  rownames(pData) <- colnames(countData)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = pData, design = ~ phenotype)
  dds <- DESeq(dds)
  dds <- replaceOutliersWithTrimmedMean(dds)
  res <- results(dds, cooksCutoff=FALSE)
  res <- res[order(res$padj),]
  as.data.frame(res)
}

###八、2个dataframe以行画红绿点,页面布局为mXn
load("D:/2016_11_24/labAandC.RData")
f1<-head(GeneLabAREP1)
f2<-head(TXLabAREP1)

r<-rowdataframe(f1,f2,2)

rowdataframe<- function(dataframe1,dataframe2,m){
  graphics.off()
  par(mfrow=c(m,ceiling(dim(dataframe1)[1]/m)))
  for (i in 1:dim(dataframe1)[1]){
    plot(1:dim(dataframe1)[2],round(dataframe1[i,]), pch=20,col="red",las=1,xlab="line",ylab="value")
    points(1:dim(dataframe1)[2],round(dataframe2[i,]), pch=20,col="blue")
  }
}

###九、gene匹配read数
f<- read.table("D:/2016_12_24/hg18.hg18",sep = "\t",header=F,stringsAsFactors = F)
f1 <- read.table("D:/2016_12_24/GRO-Seq_0h_GSM874647_GroSeq_sample1_rmdup_merged.bam.bed/GSM874647_GroSeq_sample1_rmdup_merged.bam.bed",sep = "\t",header=F,stringsAsFactors = F)
#f2 <- read.table("D:/2016_12_24/GRO-Seq_12h_WT_GSM874648_GroSeq_sample4_rmdup_merged.bam.bed/GSM874648_GroSeq_sample4_rmdup_merged.bam.bed",sep = "\t",header=F,stringsAsFactors = F)
#f3 <- read.table("D:/2016_12_24/GSM874631_Ivan_H3K4me3_T0_fastq.bowtie.bam.bed/GSM874631_Ivan_H3K4me3_T0_fastq.bowtie.bam.bed",sep = "\t",header=F,stringsAsFactors = F)
#f4 <- read.table("D:/2016_12_24/GSM874632_Ivan_H3K4me3_T12_fastq.bowtie.bam.bed/GSM874632_Ivan_H3K4me3_T12_fastq.bowtie.bam.bed",sep = "\t",header=F,stringsAsFactors = F)

f<- read.table("D:/2016_12_24/hg18.hg18",sep = "\t",header=F,stringsAsFactors = F)
f1 <- read.table("D:/2016_12_24/GRO-Seq_0h_GSM874647_GroSeq_sample1_rmdup_merged.bam.bed/GSM874647_GroSeq_sample1_rmdup_merged.bam.bed",sep = "\t",header=F,stringsAsFactors = F)
f1[7]=1

count=function(a){
  data=f1[f1[,1]==as.character(a[1]),]   
  data=data[data[,6]==as.character(a[7]),] 
  data=data[as.numeric(data[,2])>as.numeric(a[4]),]
  data=data[as.numeric(data[,3])<as.numeric(a[5]),]
  #return(sum(data[,7]))
  write.table(paste(strsplit(strsplit(a[9],";")[[1]][1]," ")[[1]][2], sum(data[,7])), file ="D:/2016_12_24/GRO-Seq_0H.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
}
c<-apply(f,1,count)


ff<-read.table("D:/2016_12_24/GRO-Seq_0H.txt",sep = "\t",header=F)
a<-substr(as.character(ff[,1]),1,15)
b<-as.numeric(substring(as.character(ff[,1]),17))
d<-data.frame(tapply(b,a,sum))
d<-cbind(ENSG=row.names(d),d)
row.names(d)=NULL
write.table(d, file ="D:/2016_12_24/GRO-Seq_0H_result.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 


###十、用R包TCGAbiolinks下载数据

##load the TCGAbiolinks package\Clinical query
library(TCGAbiolinks)
clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical")
json<-tryCatch(GDCdownload(clin.query),error = function(e) GDCdownload(clin.query, method = "client", directory = "F:/2017_03_07"))

clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
write.table(clinical.patient, file = "F:/2017_03_07/TCGA-BRCA/brca_patient.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
clinical.drug <- GDCprepare_clinic(clin.query, clinical.info = "drug")
clinical.admin <- GDCprepare_clinic(clin.query, clinical.info = "admin")
clinical.radiation <- GDCprepare_clinic(clin.query, clinical.info = "radiation")
clinical.stage_event <- GDCprepare_clinic(clin.query, clinical.info = "stage_event")
clinical.new_tumor_event <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
clinical.all <- GDCprepare_clinic(clin.query)


##Check with subtypes from TCGAprepare and update examples
BRCA_path_subtypes <- TCGAquery_subtype(tumor = "brca")
write.table(BRCA_path_subtypes, file = "F:/2017_03_07/TCGA-BRCA/brca_subtype.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


##Mutation data processing
acc.muse.maf <- GDCquery_Maf("BRCA", pipelines = "muse")
acc.varscan2.maf <- GDCquery_Maf("BRCA", pipelines = "varscan2")
acc.somaticsniper.maf <- GDCquery_Maf("BRCA", pipelines = "somaticsniper")
acc.mutect.maf <- GDCquery_Maf("BRCA", pipelines = "mutect")

write.table(acc.muse.maf, file = "F:/2017_03_07/TCGA-BRCA/brca_Mutation.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
 
##CNV data query
query <- GDCquery(project = "TCGA-BRCA",data.category = "Copy Number Variation",data.type = "Copy Number Segment")
GDCdownload(query, directory = "F:/2017_03_07")


##mRNA data query
acc.gbm <- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(acc.gbm, method = "api", directory = "F:/2017_03_07", chunks.per.download = 50)

##miRNA data query
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query, method = "client", directory = "F:/2017_03_07")
GDCdownload(query, method = "client")
GDCdownload(query, method = "api")
GDCdownload(query)
data<-GDCprepare(query)            
#miRNA data query
miRNA.query <- GDCquery(project = "TCGA-BRCA", data.category ="Gene expression",data.type="miRNA gene quantification",legacy=TRUE)
GDCdownload(miRNA.query, method = "client", directory = "F:/2017_03_07")
GDCdownload(miRNA.query, method = "client")
data<-GDCprepare(miRNA.query)
data <- GDCprepare(query = miRNA.query,
                     save = TRUE, 
                     save.filename = "exp.rda")



library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TARGET-AML", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query, method = "client")
GDCdownload(query, chunks.per.download = 50)
data <- GDCprepare(query)

###十一、4亚型对应每个表型画堆叠柱状图
##处理数据，选出对应亚型的样本级表型数据
f<-read.table("F:/2017上半年个人xi/BRCA/CNV_KMeans_BRCA_k=4.txt",header=F)
class(f)
index<-which(f[,2]=="1")
write.table(f[index,1], file = "F:/2017上半年个人xi/BRCA/subtype1.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE)

f1<-read.table("F:/2017上半年个人xi/BRCA/clinical_data.data",sep="\t",header=T)
f2 <- file("F:/2017上半年个人xi/BRCA/subtype4.txt","r")
line=readLines(f2,n=1)
while( length(line) != 0 ) {
  result=subset(f1,sampleID==line)
  write.table(result, file = "F:/2017上半年个人xi/BRCA/subtype44.txt", row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE) 
  line=readLines(f2,n=1)
}
close(f2)

##判断数据类型，factor和integer为一类（离散型a），numeric为另一类（连续型b），根据表型数据分别画4亚型堆叠柱状图
library(RColorBrewer)
f1<-read.table("F:/2017上半年个人xi/BRCA/clinical_data.data",sep="\t",header=T)
for (m in 2:length(colnames(f1))){ #除去97、125、136
  if(class(f1[,m])=="numeric"){
    g<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段
    
    NO1<-(max(sort(unique(f1[,m])))-min(sort(unique(f1[,m]))))%/%3+min(sort(unique(f1[,m]))) #min~NO1
    NO2<-(max(sort(unique(f1[,m])))-min(sort(unique(f1[,m]))))%/%3*2+min(sort(unique(f1[,m])))#NO1~NO2#NO2~max
    row.names(g)=c(paste(min(sort(unique(f1[,m]))),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(sort(unique(f1[,m])))))
    
    f1[216]=1
    da={}
    for(i in row.names(g)){
      data=f1[as.numeric(na.omit(f1[,m]))>=as.numeric(strsplit(i,"~")[[1]][1]) & as.numeric(na.omit(f1[,m]))<as.numeric(strsplit(i,"~")[[1]][2]),] 
      da=append(da,sum(data[,216]))
    }
    names(da)=row.names(g)
    da<-as.data.frame(da)
    index=f1[as.numeric(na.omit(f1[,m]))== max(sort(unique(f1[,m]))),]
    da[3,]<-da[3,]+sum(index[,216])
    
    c1<-read.table("F:/2017上半年个人xi/BRCA/subtype11.txt",sep="\t",header=F)
    c11<-as.data.frame(t(rbind(na.omit(c1[,m]),1)))
    index1=c11[as.numeric(c11[,1])< NO1,]
    g[1,1]<-sum(index1[,2])/dim(c11)[1]    
    index2=c11[as.numeric(c11[,1])>=NO1&as.numeric(c11[,1])<NO2,]
    g[2,1]<-sum(index2[,2])/dim(c11)[1]
    index3=c11[as.numeric(c11[,1])>=NO2,]
    g[3,1]<-sum(index3[,2])/dim(c11)[1]
    colnames(g)[1]=dim(c11)[1]
    
    c2<-read.table("F:/2017上半年个人xi/BRCA/subtype22.txt",sep="\t",header=F)
    c22<-as.data.frame(t(rbind(na.omit(c2[,m]),1)))
    index21=c22[as.numeric(c22[,1])< NO1,]
    g[1,2]<-sum(index21[,2])/dim(c22)[1]
    index22=c22[as.numeric(c22[,1])>=NO1&as.numeric(c22[,1])<NO2,]
    g[2,2]<-sum(index22[,2])/dim(c22)[1]
    index23=c22[as.numeric(c22[,1])>=NO2,]
    g[3,2]<-sum(index23[,2])/dim(c22)[1]
    colnames(g)[2]=dim(c22)[1]
    
    c3<-read.table("F:/2017上半年个人xi/BRCA/subtype33.txt",sep="\t",header=F)
    c33<-as.data.frame(t(rbind(na.omit(c3[,m]),1)))
    index31=c33[as.numeric(c33[,1])< NO1,]
    g[1,3]<-sum(index31[,2])/dim(c33)[1]
    index32=c33[as.numeric(c33[,1])>=NO1&as.numeric(c33[,1])<NO2,]
    g[2,3]<-sum(index32[,2])/dim(c33)[1]
    index33=c33[as.numeric(c33[,1])>=NO2,]
    g[3,3]<-sum(index33[,2])/dim(c33)[1]
    colnames(g)[3]=dim(c33)[1]
    
    c4<-read.table("F:/2017上半年个人xi/BRCA/subtype44.txt",sep="\t",header=F)
    c44<-as.data.frame(t(rbind(na.omit(c4[,m]),1)))
    index41=c44[as.numeric(c44[,1])<NO1,]
    g[1,4]<-sum(index41[,2])/dim(c44)[1]
    index42=c44[as.numeric(c44[,1])>=NO1&as.numeric(c44[,1])<NO2,]
    g[2,4]<-sum(index42[,2])/dim(c44)[1]
    index43=c44[as.numeric(c44[,1])>=NO2,]
    g[3,4]<-sum(index43[,2])/dim(c44)[1]
    colnames(g)[4]=dim(c44)[1]
    
    pdf(paste("F:/2017上半年个人xi/BRCA/barplot/样本/",colnames(f1)[m],".pdf",sep=""))
    barplot(as.matrix(g[,1:4]),main=colnames(f1)[m],col=brewer.pal(24,"Set3"),beside=FALSE,legend=paste(row.names(g),da[row.names(g),1]),args.legend=list(x="bottomleft",cex=0.8,text.width=0.15),xlab="subtype",ylab="frequency")
    dev.off()
  }
  else{
    g<-as.data.frame(matrix(numeric(0),nrow =length(as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))),ncol=4))
    row.names(g)=as.factor(sort(na.omit(gsub("^$",NA,unique(f1[,m])))))
    
    f1[216]=1
    da={}
    for(i in row.names(g)){
      data=f1[f1[,m]==i,] 
      da=append(da,sum(data[,216]))
    }
    names(da)=row.names(g)
    da<-as.data.frame(da)
    
    c1<-read.table("F:/2017上半年个人xi/BRCA/subtype11.txt",sep="\t",header=F)
    c1[216]=1
    for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
      index=c1[c1[,m]==i,]
      g[i,1]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c1[,m]))))
    }
    colnames(g)[1]=length(as.factor(na.omit(gsub("^$",NA,c1[,m]))))
    
    c2<-read.table("F:/2017上半年个人xi/BRCA/subtype22.txt",sep="\t",header=F)
    c2[216]=1
    for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
      index=c2[c2[,m]==i,]
      g[i,2]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c2[,m]))))
    }
    colnames(g)[2]=length(as.factor(na.omit(gsub("^$",NA,c2[,m]))))
    
    c3<-read.table("F:/2017上半年个人xi/BRCA/subtype33.txt",sep="\t",header=F)
    c3[216]=1
    for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
      index=c3[c3[,m]==i,]
      g[i,3]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c3[,m]))))
    }
    colnames(g)[3]=length(as.factor(na.omit(gsub("^$",NA,c3[,m]))))
    
    c4<-read.table("F:/2017上半年个人xi/BRCA/subtype44.txt",sep="\t",header=F)
    c4[216]=1
    for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
      index=c4[c4[,m]==i,]
      g[i,4]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c4[,m]))))
    }
    colnames(g)[4]=length(as.factor(na.omit(gsub("^$",NA,c4[,m]))))
    
    pdf(paste("F:/2017上半年个人xi/BRCA/barplot/样本/",colnames(f1)[m],".pdf",sep=""))
    barplot(as.matrix(g[,1:4]),main=colnames(f1)[m],col=brewer.pal(24,"Set3"),beside=FALSE,legend=paste(row.names(g),da[row.names(g),1]),args.legend=list(x="bottomleft",cex=0.8,text.width=0.15),xlab="subtype",ylab="frequency")
    dev.off()
    }
}
  
##a、表型数据是已知的几个分组（离散型），画堆叠柱状图
#Converted_Stage_nature2012---5
m=2
f1<-read.table("F:/2017上半年个人xi/BRCA/clinical_data.data",sep="\t",header=T)
g<-as.data.frame(matrix(numeric(0),nrow =length(as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))),ncol=4))
row.names(g)=as.factor(sort(na.omit(gsub("^$",NA,unique(f1[,m])))))

f1[216]=1
da={}
for(i in row.names(g)){
  data=f1[f1[,m]==i,] 
  da=append(da,sum(data[,216]))
}
names(da)=row.names(g)
da<-as.data.frame(da)

c1<-read.table("F:/2017上半年个人xi/BRCA/subtype11.txt",sep="\t",header=F)
c1[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c1[c1[,m]==i,]
  g[i,1]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c1[,m]))))
}
colnames(g)[1]=length(as.factor(na.omit(gsub("^$",NA,c1[,m]))))

c2<-read.table("F:/2017上半年个人xi/BRCA/subtype22.txt",sep="\t",header=F)
c2[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c2[c2[,m]==i,]
  g[i,2]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c2[,m]))))
}
colnames(g)[2]=length(as.factor(na.omit(gsub("^$",NA,c2[,m]))))
  
c3<-read.table("F:/2017上半年个人xi/BRCA/subtype33.txt",sep="\t",header=F)
c3[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c3[c3[,m]==i,]
  g[i,3]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c3[,m]))))
}
colnames(g)[3]=length(as.factor(na.omit(gsub("^$",NA,c3[,m]))))
  
c4<-read.table("F:/2017上半年个人xi/BRCA/subtype44.txt",sep="\t",header=F)
c4[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c4[c4[,m]==i,]
  g[i,4]<-sum(index[,216])/length(as.factor(na.omit(gsub("^$",NA,c4[,m]))))
}
colnames(g)[4]=length(as.factor(na.omit(gsub("^$",NA,c4[,m]))))
  
library(RColorBrewer)
pdf(paste("F:/2017上半年个人xi/BRCA/barplot/",colnames(f1)[m],".pdf",sep=""))
barplot(as.matrix(g[,1:4]),main=colnames(f1)[m],col=brewer.pal(24,"Set3"),beside=FALSE,legend=paste(row.names(g),da[row.names(g),1]),args.legend=list(x="bottomleft",cex=0.8,text.width=0.15),xlab="subtype",ylab="frequency")
dev.off()

##b、表型数据是一系列的数字（连续型），按区间段分3组，画堆叠柱状图
#Age_at_Initial_Pathologic_Diagnosis_nature2012----3
m=6
f1<-read.table("F:/2017上半年个人xi/BRCA/clinical_data.data",sep="\t",header=T)
g<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段

NO1<-(max(sort(unique(f1[,m])))-min(sort(unique(f1[,m]))))%/%3+min(sort(unique(f1[,m]))) #min~NO1
NO2<-(max(sort(unique(f1[,m])))-min(sort(unique(f1[,m]))))%/%3*2+min(sort(unique(f1[,m])))#NO1~NO2#NO2~max
row.names(g)=c(paste(min(sort(unique(f1[,m]))),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(sort(unique(f1[,m])))))

f1[216]=1
da={}
for(i in row.names(g)){
  data=f1[as.numeric(na.omit(f1[,m]))>=as.numeric(strsplit(i,"~")[[1]][1]) & as.numeric(na.omit(f1[,m]))<as.numeric(strsplit(i,"~")[[1]][2]),] 
  da=append(da,sum(data[,216]))
}
names(da)=row.names(g)
da<-as.data.frame(da)
index=f1[as.numeric(na.omit(f1[,m]))== max(sort(unique(f1[,m]))),]
da[3,]<-da[3,]+sum(index[,216])

c1<-read.table("F:/2017上半年个人xi/BRCA/subtype11.txt",sep="\t",header=F)
c11<-as.data.frame(t(rbind(na.omit(c1[,m]),1)))
index1=c11[as.numeric(c11[,1])< NO1,]
g[1,1]<-sum(index1[,2])/dim(c11)[1]    
index2=c11[as.numeric(c11[,1])>=NO1&as.numeric(c11[,1])<NO2,]
g[2,1]<-sum(index2[,2])/dim(c11)[1]
index3=c11[as.numeric(c11[,1])>=NO2,]
g[3,1]<-sum(index3[,2])/dim(c11)[1]
colnames(g)[1]=dim(c11)[1]

c2<-read.table("F:/2017上半年个人xi/BRCA/subtype22.txt",sep="\t",header=F)
c22<-as.data.frame(t(rbind(na.omit(c2[,m]),1)))
index21=c22[as.numeric(c22[,1])< NO1,]
g[1,2]<-sum(index21[,2])/dim(c22)[1]
index22=c22[as.numeric(c22[,1])>=NO1&as.numeric(c22[,1])<NO2,]
g[2,2]<-sum(index22[,2])/dim(c22)[1]
index23=c22[as.numeric(c22[,1])>=NO2,]
g[3,2]<-sum(index23[,2])/dim(c22)[1]
colnames(g)[2]=dim(c22)[1]
 
c3<-read.table("F:/2017上半年个人xi/BRCA/subtype33.txt",sep="\t",header=F)
c33<-as.data.frame(t(rbind(na.omit(c3[,m]),1)))
index31=c33[as.numeric(c33[,1])< NO1,]
g[1,3]<-sum(index31[,2])/dim(c33)[1]
index32=c33[as.numeric(c33[,1])>=NO1&as.numeric(c33[,1])<NO2,]
g[2,3]<-sum(index32[,2])/dim(c33)[1]
index33=c33[as.numeric(c33[,1])>=NO2,]
g[3,3]<-sum(index33[,2])/dim(c33)[1]
colnames(g)[3]=dim(c33)[1]
 
c4<-read.table("F:/2017上半年个人xi/BRCA/subtype44.txt",sep="\t",header=F)
c44<-as.data.frame(t(rbind(na.omit(c4[,m]),1)))
index41=c44[as.numeric(c44[,1])<NO1,]
g[1,4]<-sum(index41[,2])/dim(c44)[1]
index42=c44[as.numeric(c44[,1])>=NO1&as.numeric(c44[,1])<NO2,]
g[2,4]<-sum(index42[,2])/dim(c44)[1]
index43=c44[as.numeric(c44[,1])>=NO2,]
g[3,4]<-sum(index43[,2])/dim(c44)[1]
colnames(g)[4]=dim(c44)[1]

library(RColorBrewer)
pdf(paste("F:/2017上半年个人xi/BRCA/barplot/",colnames(f1)[m],".pdf",sep=""))
barplot(as.matrix(g[,1:4]),main=colnames(f1)[m],col=brewer.pal(24,"Set3"),beside=FALSE,legend=paste(row.names(g),da[row.names(g),1]),args.legend=list(x="bottomleft",cex=0.8,text.width=0.15),xlab="subtype",ylab="frequency")
dev.off()

###十二、批量读取文件名、fisher检验
all.files <-list.files("F:/2017上半年个人xi/BRCA/barplot/样本/显著/大十")
gsub(".pdf","",all.files)
write.table(gsub(".pdf","",all.files), file = "F:/2017上半年个人xi/BRCA/barplot/样本/显著/大十表型信息.txt",row.names = F,col.names = F,quote = F,sep = "\t") 

m=8
f1<-read.table("F:/2017上半年个人xi/BRCA/clinical_data.data",sep="\t",header=T)
g<-as.data.frame(matrix(numeric(0),nrow =length(as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))),ncol=4))
row.names(g)=as.factor(sort(na.omit(gsub("^$",NA,unique(f1[,m])))))

c1<-read.table("F:/2017上半年个人xi/BRCA/subtype11.txt",sep="\t",header=F)
c1[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c1[c1[,m]==i,]
  g[i,1]<-sum(index[,216])
}

c2<-read.table("F:/2017上半年个人xi/BRCA/subtype22.txt",sep="\t",header=F)
c2[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c2[c2[,m]==i,]
  g[i,2]<-sum(index[,216])
}

c3<-read.table("F:/2017上半年个人xi/BRCA/subtype33.txt",sep="\t",header=F)
c3[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c3[c3[,m]==i,]
  g[i,3]<-sum(index[,216])
}

c4<-read.table("F:/2017上半年个人xi/BRCA/subtype44.txt",sep="\t",header=F)
c4[216]=1
for (i in as.factor(na.omit(gsub("^$",NA,unique(f1[,m]))))){
  index=c4[c4[,m]==i,]
  g[i,4]<-sum(index[,216])
}

d<-t(rbind(g[2,],g[1,]+g[3,])) #Negative\non-Negative
colnames(d)[2]="Non"

g12<-d[1:2,1:2]
g13<-d[1:3,1:2]
g14<-d[1:4,1:2]
g23<-d[2:3,1:2]
g24<-d[2:4,1:2]
g34<-d[3:4,1:2]

fisher.test(g14, workspace = 200000, hybrid = FALSE,control = list(), or = 1, alternative = "two.sided",conf.int = TRUE, conf.level = 0.95)

#test <- fisher.test(g12, workspace = 200000, hybrid = FALSE,control = list(), or = 1, alternative = "two.sided",conf.int = TRUE, conf.level = 0.95)
#p<-test$p.value
#write.table(test$p.value, file = paste("F:/2017上半年个人xi/BRCA/barplot/",colnames(f1)[m],".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "", append = TRUE)


###十三、NMF非负矩阵分解、PCA主成分分析、一致聚类（分为2组：单细胞和bulk）
f<-read.table("F:/2017上半年个人xi/GSE75688/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt",sep="\t",header=T)

b={}
for (i in colnames(f)){
  if(strsplit(i,"_")[[1]][-1]=="Tumor"|strsplit(i,"_")[[1]][-1]=="Pooled"){
    b=append(b,i)
  }
}

bulk<-subset(f,select=b)
rownames(bulk)=f[,1]
bulk<-t(bulk)

sc<-subset(f,select=setdiff(colnames(f),b))   #差集setdiff
rownames(sc)=sc[,1]
sc<-t(sc[,-1:-3])

##主成分分析
pca=prcomp(sc, scale=F)
t(pca$ro)%*%pca$ro
round(t(pca$ro)%*%pca$ro)
pca$rota    #成分负荷
summary(pca)
screeplot(pca,type="line",main="碎石图",lwd=2)
pca$x       #成分得分
biplot(pca,c(1,2))  #呈现1、2成分的负荷，因子得分
#write.table(pca$rota, file="F:/2017上半年个人xi/GSE75688/rota.txt", row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)

#f[,4:6]*pca$rota[,1]
#colSums(f[,4:6]*pca$rota[,1])

z1<-colSums(t(sc)*pca$rota[,1])
z2<-colSums(t(sc)*pca$rota[,2])
z3<-colSums(t(sc)*pca$rota[,3])
z<-cbind(z1,z2,z3)

library(ConsensusClusterPlus)
title="F:/2017上半年个人xi/GSE75688/"   #目录
results = ConsensusClusterPlus(t(z),maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]
icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]
 
######lijiayao#######
f1<-read.table("F:/GSE.txt",sep="\t",header=T)
f2<-read.table("F:/miRNA标准名字.txt",sep="\t",header=T)
i<-intersect(f1[,1],f2[,1])
length(intersect(f1[,1],f2[,1]))
row.names(f1)=f1[,1]
row.names(f2)=f2[,1]
f1[i,128]=f2[i,2]
data<-cbind(f1[,128],f1[,1:ncol(f1)-1])
write.table(f1, file="F:/result.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
write.table(data, file="F:/data.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


###十四、tp53并排堆叠柱状图
#f<-read.table("F:/tp53/integrate.txt",sep="\t",header=T,stringsAsFactors = TRUE)
f<-read.table("F:/tp53/integrate.txt",sep="\t",header=T)
m<-cbind(f[,9],f[,11],f[,13],f[,15],f[,17])
#rownames(m)<-paste(f[,2],"/",f[,3])
rownames(m)<-paste(substr(as.character(f[,2]), 1,2),"/",substr(as.character(f[,3]), 1,1))

pdf("F:/tp53/plot.pdf")
barplot(t(m),main="tp53",cex.names=0.3,col=c("red","orange","green","yellow","blue"),beside=TRUE,
        legend=c(colnames(f)[9],colnames(f)[11],colnames(f)[13],colnames(f)[15],colnames(f)[17]),args.legend=list(x="topleft",cex=0.6,text.width=0.2),
        xlab=paste(colnames(f)[2],"/",colnames(f)[3]),ylab="frequency")
dev.off()

##计算tp53出现次数及个数画折线图
filenames<-list.files("F:/tp53/TCGA/",pattern = '*_tp53.csv') 
filenames<-"thym_tcga_Provisional_tp53.csv"
for(i in filenames){
  f<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T)
  f1<-read.csv(paste("F:/tp53/TCGA/",gsub("tp53.csv","clinical_data.tsv",i),sep=""),header=T,sep=" ")
  pdf(paste("F:/tp53/折线图1/",strsplit(i,".csv"),".pdf",sep=""))
  plot(as.matrix(table(table(f[,1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste(strsplit(i,"tcga")[[1]][1],"tcga(",dim(f1)[1]-1,"/",length(unique(f[,1])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
  text(as.data.frame(table(table(f[,1])))[,1],as.data.frame(table(table(f[,1])))[,2],round(as.data.frame(table(table(f[,1])))[,2]/sum(as.data.frame(table(table(f[,1])))[,2]),2))
  dev.off()
}

i="blca_tcga_Nature_tp53.csv"
paste(strsplit(i,"tcga")[[1]][1],"tcga",sep="")
f<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T)
plot(as.matrix(table(table(f[,1]))),type="o",pch=20,col= "darkblue",cex = 1,las=1,main =paste(strsplit(i,"tcga")[[1]][1],"tcga","( /",length(unique(f[,1])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(table(f[,1])))[,1],as.data.frame(table(table(f[,1])))[,2],as.data.frame(table(table(f[,1])))[,2]/sum(as.data.frame(table(table(f[,1])))[,2]))
f1<-read.csv(paste("F:/tp53/TCGA/",gsub("tp53.csv","clinical_data.tsv",i),sep=""),header=T,sep=" ")
dim(f1)[1]-1

##clinic sample\tp53 sample intersect
filenames<-list.files("F:/tp53/TCGA/",pattern = '*.tsv') 
da={}
for(m in filenames){
  da=append(da,strsplit(m,"_tcga")[[1]][1])
}

for(n in unique(da)){
  cancer=strsplit(n,"_tcga")[[1]][1]
  name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
  cat(paste0(cancer,"\t","clinic sample","\t"), file="F:/tp53/样本交集.txt", append=T)
  s1={}
  for(i in name1){
    f1<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep="\t")
    s1=append(s1,as.character(f1[,1]))
    cat(paste0(dim(f1)[1],"\t"), file="F:/tp53/样本交集.txt",append=TRUE)
  }
  cat("intersect\t",table(table(s1)),"\n", file="F:/tp53/样本交集.txt",append=TRUE)
  
  name2<-list.files("F:/tp53/TCGA/",pattern = paste(cancer,'.*.csv',sep="")) #tp53 sample
  cat(paste0(cancer,"\t","tp53 sample","\t"), file="F:/tp53/样本交集.txt", append=T)
  s2={}
  for(i in name2){
    f2<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep=",")
    s2=append(s2,unique(as.character(f2[,1])))
    cat(paste0(length(unique(as.character(f2[,1]))),"\t"), file="F:/tp53/样本交集.txt",append=TRUE)
  }
  cat("intersect\t",table(table(s2)),"\n", file="F:/tp53/样本交集.txt",append=TRUE)
}


cancer="gbm" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
cat(paste0(cancer,"\t","clinic sample","\t"), file="F:/tp53/交集.txt", append=T)
s1={}
for(i in name1[-length(name1)]){
  f1<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep="\t")
  s1=append(s1,unique(as.character(f1[,1])))
  cat(paste0(length(unique(as.character(f1[,1]))),"\t"), file="F:/tp53/交集.txt",append=TRUE)
}
cat("intersect\t",table(table(s1)),"\n", file="F:/tp53/交集.txt",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern = paste(cancer,'.*.csv',sep="")) #tp53 sample
cat(paste0(cancer,"\t","tp53 sample","\t"), file="F:/tp53/交集.txt", append=T)
s2={}
for(i in name2[-length(name2)]){
  f2<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep=",")
  s2=append(s2,unique(as.character(f2[,1])))
  cat(paste0(length(unique(as.character(f2[,1]))),"\t"), file="F:/tp53/交集.txt",append=TRUE)
}
cat("intersect\t",table(table(s2)),"\n", file="F:/tp53/交集.txt",append=TRUE)


##画韦恩图
install.packages("VennDiagram")
library(grid)
library(VennDiagram)
cancer="gbm" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
f1<-read.csv(paste("F:/tp53/TCGA/",name1[1],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,1]))
f2<-read.csv(paste("F:/tp53/TCGA/",name1[2],sep=""),header=T,sep="\t")
B=unique(as.character(f2[,1]))                    
f3<-read.csv(paste("F:/tp53/TCGA/",name1[3],sep=""),header=T,sep="\t")
C=unique(as.character(f3[,1]))                     
D<-venn.diagram(list(A=A,B=B,C=C),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
grid.draw(D)

##数据的重新整理
cancer="stad" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
s1={}
for(i in name1){
  f1<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep="\t")
  s1=append(s1,unique(as.character(f1[,1])))
}
write.table(unique(s1), file=paste("F:/tp53/new data/",cancer,"_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.csv',sep="")) #tp53 sample
s2={}
for(i in name2){
  f2<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep=",")
  s2=append(s2,unique(as.character(f2[,1])))
}
write.table(unique(s2), file=paste("F:/tp53/new data/",cancer,"_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


cancer="gbm" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
s1={}
for(i in name1[-length(name1)]){
  f1<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep="\t")
  s1=append(s1,unique(as.character(f1[,1])))
}
write.table(unique(s1), file=paste("F:/tp53/new data/",cancer,"_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.csv',sep="")) #tp53 sample
s2={}
for(i in name2[-length(name2)]){
  f2<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T,sep=",")
  s2=append(s2,unique(as.character(f2[,1])))
}
write.table(unique(s2), file=paste("F:/tp53/new data/",cancer,"_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


cancer="ov" 
f<-read.table("F:/tp53/ov.txt",header=T,sep="\t")
write.table(intersect(f[,1],f[,2]), file=paste("F:/tp53/new data/",cancer,"_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(intersect(f[,3],f[,4]), file=paste("F:/tp53/new data/",cancer,"_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



#红的包含\blca\coadread\luad\lusc\prad?\thca\ucec
cancer="ucec" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
f1<-read.csv(paste("F:/tp53/TCGA/",name1[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste("F:/tp53/TCGA/",name1[2],sep=""),header=T,sep="\t")
#setdiff(f2[,1],f1[,1])    #求属于x而不属于y的所有元素
write.table(f1[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f2[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f2[,1],f1[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.csv',sep="")) #tp53 sample
f3<-read.csv(paste("F:/tp53/TCGA/",name2[1],sep=""),header=T,sep=",")
f4<-read.csv(paste("F:/tp53/TCGA/",name2[2],sep=""),header=T,sep=",")
write.table(f3[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f4[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f4[,1],f3[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

#紫的\hnsc\kirc\ov????\stad
cancer="stad" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
f1<-read.csv(paste("F:/tp53/TCGA/",name1[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste("F:/tp53/TCGA/",name1[2],sep=""),header=T,sep="\t")
union(f1[,1],f2[,1])    #求并集
intersect(f1[,1],f2[,1])    #求交集
setdiff(f2[,1],intersect(f1[,1],f2[,1]))
write.table(setdiff(f2[,1],intersect(f1[,1],f2[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f1[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(f1[,1],f2[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.csv',sep="")) #tp53 sample
f3<-read.csv(paste("F:/tp53/TCGA/",name2[1],sep=""),header=T,sep=",")
f4<-read.csv(paste("F:/tp53/TCGA/",name2[2],sep=""),header=T,sep=",")
write.table(setdiff(f4[,1],intersect(f3[,1],f4[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f3[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(f3[,1],f4[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


cancer="ov" ##数据有问题，自己整理的
f<-read.table("F:/tp53/ov.txt",header=T,sep="\t")
write.table(setdiff(f[,2],intersect(f[,1],f[,2])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f[,1], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(f[,1],f[,2]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

write.table(setdiff(f[,4],intersect(f[,3],f[,4])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(f[,3], file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(f[,3],f[,4]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


#绿的\brca\gbm
cancer="gbm" 
name1<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.tsv',sep="")) #clinic sample
f1<-read.csv(paste("F:/tp53/TCGA/",name1[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste("F:/tp53/TCGA/",name1[2],sep=""),header=T,sep="\t")
f3<-read.csv(paste("F:/tp53/TCGA/",name1[3],sep=""),header=T,sep="\t")
write.table(setdiff(f1[,1],intersect(f1[,1],f2[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f2[,1],setdiff(intersect(f2[,1],f3[,1]),intersect(intersect(f1[,1],f2[,1]),f3[,1]))), file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f3[,1],intersect(f1[,1],f3[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(union(f1[,1],f2[,1]),f3[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_4_clinic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

name2<-list.files("F:/tp53/TCGA/",pattern =paste(cancer,'.*.csv',sep="")) #tp53 sample
f4<-read.csv(paste("F:/tp53/TCGA/",name2[1],sep=""),header=T,sep=",")
f5<-read.csv(paste("F:/tp53/TCGA/",name2[2],sep=""),header=T,sep=",")
f6<-read.csv(paste("F:/tp53/TCGA/",name2[3],sep=""),header=T,sep=",")
write.table(setdiff(f4[,1],intersect(f4[,1],f5[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_1_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f5[,1],setdiff(intersect(f5[,1],f6[,1]),intersect(intersect(f1[,1],f2[,1]),f3[,1]))), file=paste("F:/tp53/TCGAintegrate/",cancer,"_2_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(setdiff(f6[,1],intersect(f4[,1],f6[,1])), file=paste("F:/tp53/TCGAintegrate/",cancer,"_3_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(union(union(f4[,1],f5[,1]),f6[,1]), file=paste("F:/tp53/TCGAintegrate/",cancer,"_4_tp53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


#个别cancer(blca\brca\gbm\hnsc\kirc\prad\stad\ov)
filenames<-list.files("F:/tp53/TCGAintegrate/",pattern = paste("ov",'.*_tp53.txt',sep="")) 
for(i in filenames){
  f<-read.csv(paste("F:/tp53/TCGAintegrate/",i,sep=""),header=F)
  f1<-read.csv(paste("F:/tp53/TCGAintegrate/",gsub("tp53.txt","clinic.txt",i),sep=""),header=F,sep=" ")
  pdf(paste("F:/tp53/TCGAintegrate折线图/",strsplit(i,".txt"),".pdf",sep=""))
  plot(as.matrix(table(table(f[,1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste(strsplit(i,"tp53")[[1]][1],"tcga(",dim(f1)[1],"/",length(unique(f[,1])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
  text(as.data.frame(table(table(f[,1])))[,1],as.data.frame(table(table(f[,1])))[,2],round(as.data.frame(table(table(f[,1])))[,2]/sum(as.data.frame(table(table(f[,1])))[,2]),2))
  dev.off()
}


###十五、cBioPortal TCGA cancer clinic data process
###（一、总.txt 二、M.txt 三、M_FEMALE.txt 四、韦恩图 五、折线图 六、柱状图 七、散点图）
##一、总.txt
PATH="F:/tp53/clinic_data/"
path="F:/tp53/clinic_data_process/"
path1="F:/tp53/clinic_data_process_result/"
cancer="pcpg"                    ##注意注意要改！！！！！！！！
fname<-list.files(PATH,pattern = paste(cancer,'.*.txt',sep="")) 
for(name in fname){
  f<-read.csv(paste(PATH,name,sep=""),header=T,sep="\t")
  rep<-table(f[,4])[table(f[,4])>1]
  #as.data.frame(rep)[,1] #7个
  for(i in 1:dim(f)[1]+1){
    if(f[i,4] %in% as.data.frame(rep)[,1]){
      write.table(f[i,], file = paste(path,strsplit(name,"study_view_clinical_data.txt")[[1]][1],"1vs2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
    }
  }
  cancer_type_na<-is.na(f[,2])
  #f[cancer_type_na,4] #9个

  index<-which(f[,2]=="Pheochromocytoma")##注意注意要改！！！！！！！！
  #f[index,] #1093行
  write.table(f[index,], file = paste(path,strsplit(name,"study_view_clinical_data.txt")[[1]][1],"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  write.table(f[index,], file = paste(path,strsplit(name,"study_view_clinical_data.txt")[[1]][1],"总.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  
  for(i in 1:dim(f)[1]+1){
    if(is.na(f[i,2]) & f[i,4] %in% as.data.frame(rep)[,1]){#7行
      write.table(f[i,], file = paste(path,strsplit(name,"study_view_clinical_data.txt")[[1]][1],"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
      write.table(f[i,], file =paste(path,strsplit(name,"study_view_clinical_data.txt")[[1]][1],"总.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
    }
  }
}

##二、总.txt->M.txt
path="F:/tp53/clinic_data_process/"
path1="F:/tp53/clinic_data_process_result/"
cancer="pcpg" 
fname<-list.files(path,pattern = paste(cancer,'.*总.txt',sep="")) 
for(name in fname){
  if (is.na(file.info(paste(path,name,sep=""))[1])){
    next
  }
  f<-read.csv(paste(path,name,sep=""),header=F,sep="\t")
  index<-which(f[,8]=="YES")
  write.table(f[index,], file = paste(path1,strsplit(name,"总.txt")[[1]][1],"M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
}

##三、M.txt->M_FEMALE.txt
path1="F:/tp53/clinic_data_process_result/"
cancer="pcpg" 
fname1<-list.files(path1,pattern = paste(cancer,'.*M.txt',sep="")) 
for(name1 in fname1){
  f<-read.csv(paste(path1,name1,sep=""),header=F,sep="\t")
  index1<-which(f[,34]=="FEMALE")
  write.table(f[index1,], file = paste(path1,strsplit(name1,"M.txt")[[1]][1],"M_FEMALE.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  index2<-which(f[,34]=="MALE")
  write.table(f[index2,], file = paste(path1,strsplit(name1,"M.txt")[[1]][1],"M_MALE.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  
  rep<-table(f[,4])[table(f[,4])>1]
  for(i in 1:dim(f)[1]+1){
    if(f[i,4] %in% as.data.frame(rep)[,1]){
      write.table(f[i,], file = paste(path1,strsplit(name1,"M.txt")[[1]][1],"Primary_Metastatic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
    }
  }
  index3<-which(f[,34]!="MALE" & f[,12]!="Metastatic")
  write.table(f[index3,], file = paste(path1,strsplit(name1,"M.txt")[[1]][1],"M_del_MALE_Metastatic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
}

##四、韦恩图
#install.packages("VennDiagram")
library(grid)
library(VennDiagram)
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="pcpg" 
fname<-list.files(path1,pattern = paste(cancer,'.*.txt',sep="")) 
for(name in fname){
  f<-read.csv(paste(path1,name,sep=""),header=F,sep="\t")
  f1<-read.csv(paste(path2,strsplit(name,"_")[[1]][1],"_tcga_",strsplit(name,"_")[[1]][2],"_tp53.csv",sep=""),header=T,sep=",")
  A=unique(as.character(f[,5]))
  B=unique(as.character(f1[,1]))
  pdf(paste(path1,"韦恩图/",strsplit(name,".txt")[[1]][1],"_tp53.pdf",sep=""))
  D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
  grid.draw(D)
  dev.off()
}

##其余16种cancer
library(grid)
library(VennDiagram)
path1="F:/tp53/TCGA16cancer/process_result/" 
cancername<-list.files(path1,pattern = paste('.*_M.txt',sep=""))
for(i in cancername){
  cancer=gsub("_M.txt","",i)
  name<-list.files(path1,pattern = paste(cancer,'.*.txt',sep=""))
  f1<-read.csv(paste(path1,name[1],sep=""),header=F,sep="\t")
  f2<-read.csv(paste(path1,name[2],sep=""),header=F,sep="\t")
  A=unique(as.character(f1[,5]))
  B=unique(as.character(f2[,1]))
  pdf(paste(path1,"韦恩图/",strsplit(name[1],".txt")[[1]][1],"_tp53.pdf",sep=""))
  D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
  grid.draw(D)
  dev.off()
}

cancer="laml" 
name<-list.files(path1,pattern = paste(cancer,'.*.txt',sep=""))
f1<-read.csv(paste(path1,name[1],sep=""),header=F,sep="\t")
f2<-read.csv(paste(path1,name[2],sep=""),header=F,sep="\t")
A=unique(as.character(f1[,5]))
B=unique(as.character(f2[,1]))
pdf(paste(path1,"韦恩图/",strsplit(name[1],".txt")[[1]][1],"_tp53.pdf",sep=""))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
dev.off()
 

############################################################
f<-read.csv("F:/tp53/clinic_data/brca_Provisional_study_view_clinical_data.txt",header=T,sep="\t")
#f<-read.csv("F:/tp53/clinic_data/brca_Cell_study_view_clinical_data.txt",header=T,sep="\t")

rep<-table(f[,4])[table(f[,4])>1]
#as.data.frame(rep)[,1] #7个
for(i in 1:dim(f)[1]+1){
  if(f[i,4] %in% as.data.frame(rep)[,1]){
    write.table(f[i,], file = paste(path,"1vs2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  }
}

cancer_type_na<-is.na(f[,2])
#f[cancer_type_na,4] #9个

nouse<-setdiff(f[cancer_type_na,4],as.data.frame(rep)[,1]) #2个
for(i in 1:dim(f)[1]+1){
  if(is.na(f[i,2]) & f[i,4] %in% nouse){
    write.table(f[i,], file = paste(path,"nouse.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  }
}

index<-which(f[,2]=="Breast Cancer")
#f[index,] #1093行
write.table(f[index,], file = paste(path,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
write.table(f[index,], file = paste(path,"总.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)

for(i in 1:dim(f)[1]+1){
  if(is.na(f[i,2]) & f[i,4] %in% as.data.frame(rep)[,1]){#7行
    write.table(f[i,], file = paste(path,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
    write.table(f[i,], file =paste(path,"总.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  }
}

f<-read.csv("F:/tp53/clinic_data/brca_Nature_study_view_clinical_data.txt",header=T,sep="\t")
#f<-read.csv("F:/tp53/clinic_data/brca_Cell_study_view_clinical_data.txt",header=T,sep="\t")
#f<-read.csv(paste(path,"brca_Provisional_总.txt",sep=""),header=F,sep="\t")
index<-which(f[,8]=="YES")
#write.table(f[index,], file = paste(path,"M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
write.table(f[index,], file = paste(path1,"brca_Nature_","M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)


f<-read.csv(paste(path1,"brca_Provisional_","M.txt",sep=""),header=F,sep="\t")
index1<-which(f[,34]=="FEMALE")
write.table(f[index1,], file = paste(path1,"brca_Provisional_","M_FEMALE.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
index2<-which(f[,34]=="MALE")
#write.table(f[index2,], file = paste(path,"brca_Provisional_","M_MALE.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
write.table(f[index2,], file = paste(path1,"brca_Nature_","M_MALE.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)


f<-read.csv(paste(path1,"brca_Nature_","M.txt",sep=""),header=F,sep="\t")
rep<-table(f[,4])[table(f[,4])>1]
#as.data.frame(rep)[,1] #5个
for(i in 1:dim(f)[1]+1){
  if(f[i,4] %in% as.data.frame(rep)[,1]){
    write.table(f[i,], file = paste(path1,"brca_Nature_","Primary_Metastatic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)
  }
}

index3<-which(f[,34]!="MALE" & f[,12]!="Metastatic")
write.table(f[index3,], file = paste(path1,"brca_Nature_","M_del_MALE_Metastatic.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)

################################################

#install.packages("VennDiagram")
library(grid)
library(VennDiagram)
path1="F:/tp53/clinic_data_process_result/"
f<-read.csv(paste(path1,"brca_Provisional_M_FEMALE.txt",sep=""),header=F,sep="\t")
f1<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
A=unique(as.character(f[,5]))
B=unique(as.character(f1[,1]))
pdf(paste(path1,"韦恩图/","brca_Provisional_M_FEMALE_tp53.pdf",sep=""))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
dev.off()

install.packages("VennDiagram") 
library(grid)
library(VennDiagram)
path1="F:/tp53/clinic_data_process_result/"
f1<-read.csv(paste(path1,"brca_Cell_M.txt",sep=""),header=F,sep="\t")
f2<-read.csv(paste(path1,"brca_Nature_M.txt",sep=""),header=F,sep="\t")
f3<-read.csv(paste(path1,"brca_Provisional_M.txt",sep=""),header=F,sep="\t")
A=unique(as.character(f1[,5]))
B=unique(as.character(f2[,5]))
C=unique(as.character(f3[,5]))
pdf(paste(path1,"韦恩图/","brca_3_M.pdf",sep=""))
D<-venn.diagram(list(A=A,B=B,C=C),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
grid.draw(D)
dev.off()

f4<-read.csv("F:/tp53/TCGA/brca_tcga_Cell_tp53.csv",header=T,sep=",")
f5<-read.csv("F:/tp53/TCGA/brca_tcga_Nature_tp53.csv",header=T,sep=",")
f6<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
E=unique(as.character(f4[,1]))
#G=unique(as.character(f5[,1]))
H=unique(as.character(f6[,1]))
pdf(paste(path1,"韦恩图/","2.pdf",sep=""))
D<-venn.diagram(list(A=E,C=H),filename=NULL,lwd=1,lty=2,col=c('red','blue'),fill=c('red','blue'),cat.col=c('red','blue'),reverse=TRUE)
grid.draw(D)
dev.off()

path1="F:/tp53/clinic_data_process_result/"
f1<-read.csv(paste(path1,"brca_Cell_M.txt",sep=""),header=F,sep="\t")
f2<-read.csv(paste(path1,"brca_Nature_M.txt",sep=""),header=F,sep="\t")
f3<-read.csv(paste(path1,"brca_Provisional_M.txt",sep=""),header=F,sep="\t")
A=unique(as.character(f1[,5]))
#B=unique(as.character(f2[,5]))
C=unique(as.character(f3[,5]))
pdf(paste(path1,"韦恩图/","1.pdf",sep=""))
D<-venn.diagram(list(A=A,C=C),filename=NULL,lwd=1,lty=2,col=c('red','blue'),fill=c('red','blue'),cat.col=c('red','blue'),reverse=TRUE)
grid.draw(D)
dev.off()



f4<-read.csv("F:/tp53/TCGA/brca_tcga_Cell_tp53.csv",header=T,sep=",")
f6<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
write.table(setdiff(f4[,1],f6[,1]), file = paste(path1,"韦恩图/","TCGA.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)

f2<-read.csv(paste(path1,"brca_Provisional_Primary_Metastatic.txt",sep=""),header=F,sep="\t")
#f3<-read.csv(paste(path1,"brca_Provisional_M.txt",sep=""),header=F,sep="\t")
f6<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
write.table(intersect(f2[,5],f6[,1]), file = paste(path1,"韦恩图/","TCGA2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t", append = TRUE)

###################################################################

##五、折线图
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="pcpg" 
fname<-list.files(path1,pattern = paste(cancer,'.*.txt',sep="")) 
for(name in fname){
  f<-read.csv(paste(path1,name,sep=""),header=F,sep="\t")
  f1<-read.csv(paste(path2,strsplit(name,"_")[[1]][1],"_tcga_",strsplit(name,"_")[[1]][2],"_tp53.csv",sep=""),header=T,sep=",")
  i<-intersect(f1[,1],f[,5])
  table(table(f1[(f1[,1] %in% i),1]))
  pdf(paste(path1,"折线图/",strsplit(name,".txt")[[1]][1],"_tp53.pdf",sep=""))
  plot(as.matrix(table(table(f1[(f1[,1] %in% i),1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste(strsplit(name,".txt")[[1]][1],"_tp53","(",dim(f)[1],"/",length(intersect(f1[,1],f[,5])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
  text(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,1],as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2],round(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2]/sum(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2]),2))
  if (0 %in% as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,1]){
    mtext("comutation number-1 in tp53",side =3)
  }
  dev.off()
}

##其余16种cancer
path1="F:/tp53/TCGA16cancer/process_result/" 
cancername<-list.files(path1,pattern = paste('.*_M.txt',sep=""))
for(i1 in cancername){
  cancer=gsub("_M.txt","",i1)
  name<-list.files(path1,pattern = paste(cancer,'.*.txt',sep="")) 
  f<-read.csv(paste(path1,name[1],sep=""),header=F,sep="\t")
  f1<-read.csv(paste(path1,name[2],sep=""),header=F,sep="\t")
  i<-intersect(f1[,1],f[,5])
  table(table(f1[(f1[,1] %in% i),1]))
  pdf(paste(path1,"折线图/",strsplit(name[1],".txt")[[1]][1],"_tp53.pdf",sep=""))
  plot(as.matrix(table(table(f1[(f1[,1] %in% i),1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste(strsplit(name[1],".txt")[[1]][1],"_tp53","(",dim(f)[1],"/",length(intersect(f1[,1],f[,5])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
  text(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,1],as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2],round(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2]/sum(as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,2]),2))
  if (0 %in% as.data.frame(table(table(f1[(f1[,1] %in% i),1])))[,1]){
    mtext("comutation number-1 in tp53",side =3)
  }
  dev.off()
}
  
  
path1="F:/tp53/clinic_data_process_result/"
f<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
f1<-read.csv(paste(path1,"brca_Provisional_M.txt",sep=""),header=F,sep="\t")
pdf(paste(path1,"折线图/","brca_Provisional_M_tp53.pdf",sep=""))
plot(as.matrix(table(table(f[,1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste("brca_Provisional_M_tp53","(",dim(f1)[1],"/",length(unique(f[,1])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(table(f[,1])))[,1],as.data.frame(table(table(f[,1])))[,2],round(as.data.frame(table(table(f[,1])))[,2]/sum(as.data.frame(table(table(f[,1])))[,2]),2))
dev.off()
 
path1="F:/tp53/clinic_data_process_result/"
f<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
f1<-read.csv(paste(path1,"brca_Provisional_M_MALE.txt",sep=""),header=F,sep="\t")
i<-intersect(f[,1],f1[,5])
table(table(f[(f[,1] %in% i),1]))
pdf(paste(path1,"折线图/","brca_Provisional_M_MALE_tp53.pdf",sep=""))
plot(as.matrix(table(table(f[(f[,1] %in% i),1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste("brca_Provisional_M_MALE_tp53","(",dim(f1)[1],"/",length(intersect(f[,1],f1[,5])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1],as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2],round(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]/sum(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]),2))
if (0 %in% as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1]){
  #p=as.numeric(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1])-1
  mtext("comutation number-1 in tp53",side =3)
}
dev.off()


path1="F:/tp53/clinic_data_process_result/"
f<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
f1<-read.csv(paste(path1,"brca_Provisional_M_del_MALE_Metastatic.txt",sep=""),header=F,sep="\t")
i<-intersect(f[,1],f1[,5])
table(table(f[(f[,1] %in% i),1]))
pdf(paste(path1,"折线图/","brca_Provisional_M_del_MALE_Metastatic_tp53.pdf",sep=""))
plot(as.matrix(table(table(f[(f[,1] %in% i),1]))),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste("brca_Provisional_M_del_MALE_Metastatic_tp53","(",dim(f1)[1],"/",length(intersect(f[,1],f1[,5])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1],as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2],round(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]/sum(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]),2))
if (0 %in% as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1]){
  #p=as.numeric(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1])-1
  mtext("comutation number-1 in tp53",side =3)
}
dev.off()


path1="F:/tp53/clinic_data_process_result/"
f<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
f1<-read.csv(paste(path1,"brca_Provisional_Primary_Metastatic.txt",sep=""),header=F,sep="\t")
i<-intersect(f[,1],f1[,5])
table(table(f[(f[,1] %in% i),1]))
pdf(paste(path1,"折线图/","brca_Provisional_Primary_Metastatic_tp53.pdf",sep=""))
plot(as.matrix(table(table(f[(f[,1] %in% i),1]))),type="b",pch=20,col= "darkblue",cex = 1,las=1,main =paste("brca_Provisional_Primary_Metastatic_tp53","(",dim(f1)[1],"/",length(intersect(f[,1],f1[,5])),")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1],as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2],round(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]/sum(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,2]),2))
if (0 %in% as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1]){
  #p=as.numeric(as.data.frame(table(table(f[(f[,1] %in% i),1])))[,1])-1
  mtext("comutation number-1 in tp53",side =3)
}
dev.off()


##################################################################

##六、柱状图
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="pcpg" 
name<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=T,sep=",")
a<-as.data.frame(table(f[,5]))
rownames(a)<-a[,1]
Missense<-a["Missense_Mutation",2]
Truncating<-sum(a["Frame_Shift_Ins",2],a["Frame_Shift_Del",2],a["Nonsense_Mutation",2],a["Splice_Region",2],a["Splice_Site",2],na.rm=TRUE)
Inframe<-sum(a["In_Frame_Del",2],a["In_Frame_Ins",2],na.rm=TRUE)
b<-a[!a$Var1 %in% c("In_Frame_Del","Missense_Mutation","In_Frame_Ins"),]
#t(subset(b,select=-Var1))
b1<-a[a$Var1 %in% c("In_Frame_Del","In_Frame_Ins"),]
#t(subset(b1,select=-Var1))
m<-cbind(Missense,Truncating,t(subset(b,select=-Var1)),Inframe,t(subset(b1,select=-Var1)))
pdf(paste(path1,"柱状图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_tp53.pdf",sep=""))
barplot(m,main=paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_tp53 mutation",sep=""),cex.names=0.3,beside=TRUE, space=c(0,1),width=0.1,las=1,
        col=c("pink",rep("lightgreen",dim(m)[2]-2-dim(t(subset(b1,select=-Var1)))[2]),rep("orange",dim(t(subset(b1,select=-Var1)))[2]+1)),xlab="TP53 mutation type",ylab="number")
for (i in 1:dim(m)[2]){
  text(i/16+0.1*(i+1.5),m[1,i]-1,m[1,i],cex=0.7)
}
dev.off()
m<-as.data.frame(m)
m[is.na(m)] <- 0
m[2,]<-round(m[1,]/sum(Missense,Truncating,Inframe,na.rm=TRUE),3)
pdf(paste(path1,"柱状图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_tp53_1.pdf",sep=""))
barplot(as.matrix(m[2,]),main=paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_tp53 mutation",sep=""),cex.names=0.3,beside=TRUE, space=c(0,1),width=0.1,las=1,
        col=c("pink",rep("lightgreen",dim(m)[2]-2-dim(t(subset(b1,select=-Var1)))[2]),rep("orange",dim(t(subset(b1,select=-Var1)))[2]+1)),xlab="TP53 mutation type",ylab="mutation frequency")
for (i in 1:dim(m)[2]){
  text(i/16+0.1*(i+1.5),round((m[1,i]-1)/sum(Missense,Truncating,Inframe,na.rm=TRUE),3),round(m[1,i]/sum(Missense,Truncating,Inframe,na.rm=TRUE),3),cex=0.7)
}
dev.off()

##其余16种cancer
path1="F:/tp53/TCGA16cancer/process_result/" 
cancername<-list.files(path1,pattern = paste('.*_TP53.txt',sep=""))
for(i1 in cancername){
  cancer=gsub("_TP53.txt","",i1)
  name<-list.files(path1,pattern = paste(cancer,'.*_TP53.txt',sep="")) 
  f<-read.csv(paste(path1,name,sep=""),header=F,sep="\t")
  a<-as.data.frame(table(f[,5]))
  rownames(a)<-a[,1]
  Missense<-a["Missense_Mutation",2]
  Truncating<-sum(a["Frame_Shift_Ins",2],a["Frame_Shift_Del",2],a["Nonsense_Mutation",2],a["Splice_Region",2],a["Splice_Site",2],na.rm=TRUE)
  Inframe<-sum(a["In_Frame_Del",2],a["In_Frame_Ins",2],na.rm=TRUE)
  b<-a[!a$Var1 %in% c("In_Frame_Del","Missense_Mutation","In_Frame_Ins"),]
  #t(subset(b,select=-Var1))
  b1<-a[a$Var1 %in% c("In_Frame_Del","In_Frame_Ins"),]
  #t(subset(b1,select=-Var1))
  m<-cbind(Missense,Truncating,t(subset(b,select=-Var1)),Inframe,t(subset(b1,select=-Var1)))
  pdf(paste(path1,"柱状图/",strsplit(name,".txt")[[1]][1],".pdf",sep=""))
  barplot(m,main=paste(strsplit(name,"_")[[1]][1],"_tp53 mutation",sep=""),cex.names=0.3,beside=TRUE, space=c(0,1),width=0.1,las=1,
          col=c("pink",rep("lightgreen",dim(m)[2]-2-dim(t(subset(b1,select=-Var1)))[2]),rep("orange",dim(t(subset(b1,select=-Var1)))[2]+1)),xlab="TP53 mutation type",ylab="number")
  for (i in 1:dim(m)[2]){
    text(i/16+0.1*(i+1.5),m[1,i]-1,m[1,i],cex=0.7)
  }
  dev.off()
  m<-as.data.frame(m)
  m[is.na(m)] <- 0
  m[2,]<-round(m[1,]/sum(Missense,Truncating,Inframe,na.rm=TRUE),3)
  pdf(paste(path1,"柱状图/",strsplit(name,".txt")[[1]][1],"_1.pdf",sep=""))
  barplot(as.matrix(m[2,]),main=paste(strsplit(name,"_")[[1]][1],"_tp53 mutation",sep=""),cex.names=0.3,beside=TRUE, space=c(0,1),width=0.1,las=1,
          col=c("pink",rep("lightgreen",dim(m)[2]-2-dim(t(subset(b1,select=-Var1)))[2]),rep("orange",dim(t(subset(b1,select=-Var1)))[2]+1)),xlab="TP53 mutation type",ylab="mutation frequency")
  for (i in 1:dim(m)[2]){
    text(i/16+0.1*(i+1.5),round((m[1,i]-1)/sum(Missense,Truncating,Inframe,na.rm=TRUE),3),round(m[1,i]/sum(Missense,Truncating,Inframe,na.rm=TRUE),3),cex=0.7)
  }
  dev.off()
}



path1="F:/tp53/clinic_data_process_result/"
f<-read.csv("F:/tp53/TCGA/brca_tcga_Provisional_tp53.csv",header=T,sep=",")
a<-as.data.frame(table(f[,5]))
rownames(a)<-a[,1]
Missense<-a["Missense_Mutation",2]
Truncating<-sum(a["Frame_Shift_Ins",2],a["Frame_Shift_Del",2],a["Nonsense_Mutation",2],a["Splice_Region",2],a["Splice_Site",2],na.rm=TRUE)
Inframe<-a["In_Frame_Del",2] 
b<-a[!a$Var1 %in% c("In_Frame_Del","Missense_Mutation"),]
#t(subset(b,select=-Var1))
m<-cbind(Missense,Truncating,t(subset(b,select=-Var1)),Inframe)
pdf(paste(path1,"柱状图/","brca_Provisional_tp53.pdf",sep=""))
barplot(m,main="brca_Provisional_tp53 mutation",cex.names=0.4,col="red",beside=TRUE, space=c(0,1),width=0.1,las=1,ylim=c(0,max(m[1,])+1),
        xlab="TP53 mutation type",ylab="number")
for (i in 1:dim(m)[2]){
  text(i/16+0.1*(i+1.5),m[1,i]-0.5,m[1,i],cex=0.7)
}
dev.off()



##七、散点图
path2="F:/tp53/TCGA16cancer/process_result/" 
cancer="laml" 
name<-list.files(path2,pattern = paste(cancer,'.*_TP53.txt',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=F,sep="\t")
#path1="F:/tp53/clinic_data_process_result/"
#path2="F:/tp53/TCGA/"
#cancer="pcpg" 
#name<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
#f<-read.csv(paste(path2,name,sep=""),header=T,sep=",")
aa<-as.data.frame(cbind(f[,13],as.character(f[,15])))
aa[,3]<-0
for(i in 1:dim(aa)[1]){
  if(!"-" %in% aa[i,2]){
    aa[i,3]<-nchar(as.character(aa[i,2]))
  }
}
b<-data.frame(c1=as.character(unique(aa[,1])))
row.names(b)<-b[,1]
for(i in as.character(b[,1])){
  index=aa[as.character(aa[,1])==i,]
  b[i,2]=sum(as.numeric(index[,3]))
}
n<-as.data.frame(matrix(0,max(as.numeric(as.character(b[,1])))-min(as.numeric(as.character(b[,1])))+1,2))
rownames(n)<-min(as.numeric(as.character(b[,1]))):max(as.numeric(as.character(b[,1])))
n[,1]<-min(as.numeric(as.character(b[,1]))):max(as.numeric(as.character(b[,1])))
for(i in n$V1){
  if(is.na(b[as.character(i),"V2"])){
    next
  }
  n[as.character(i),2]=b[as.character(i),2]
}
#pdf(paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base.pdf",sep=""))
#plot(n,pch=20,col= "purple",cex = 1,las=1,main =paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""),xlab="Start Pos",ylab="base mutation number")
#dev.off()
#pdf(paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_1.pdf",sep=""))
#plot(n[,1],n[,2]/sum(n[,2]),pch=20,col="purple",cex = 1,las=1,main =paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""),xlab="Start Pos",ylab="base mutation frequency")
#dev.off()

pdf(paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base.pdf",sep=""))
plot(n,pch=20,col= "purple",cex = 1,las=1,main =paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""),xlab="Start Pos",ylab="base mutation number")
dev.off()
pdf(paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base_1.pdf",sep=""))
plot(n[,1],n[,2]/sum(n[,2]),pch=20,col="purple",cex = 1,las=1,main =paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""),xlab="Start Pos",ylab="base mutation frequency")
dev.off()


#a\散点图标颜色有点突变、有段突变的
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="ucs" 
name<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=T,sep=",")
aa<-as.data.frame(cbind(paste(as.character(f[,5]),as.character(f[,13]),sep="/"),as.character(f[,13]),as.character(f[,15])))
aa[,4]<-0
for(i in 1:dim(aa)[1]){
  if(!"-" %in% aa[i,3]){
    aa[i,4]<-nchar(as.character(aa[i,3]))
  }
}
e1<-subset(aa,aa[,4]<=1) #base<=1
#e2<-subset(aa,aa[,4]>1)  #base>1
b1<-data.frame(c1=as.character(unique(e1[,1])))
row.names(b1)<-b1[,1]
b1$c2 <- gsub("^.*/","",b1$c1)
for(i in as.character(b1[,1])){
  index=e1[as.character(e1[,1])==i,]
  b1[i,3]=sum(as.numeric(index[,4]))
}
n1<-as.data.frame(matrix(0,max(as.numeric(as.character(f[,13])))-min(as.numeric(as.character(f[,13])))+1,2))
rownames(n1)<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
n1[,1]<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
m1<-merge(n1,b1,by.x = "V1",by.y = "c2",all=T)
m1[is.na(m1[,"V3"]),"V3"]=0
m1$V2 <- gsub("/.*$*","",m1$c1)

e2<-subset(aa,aa[,4]>1)  #base>1
e2[,5]<-1
b2<-data.frame(c1=as.character(unique(e2[,1])))
row.names(b2)<-b2[,1]
b2$c2 <- gsub("^.*/","",b2$c1)
for(i in as.character(b2[,1])){
  index=e2[as.character(e2[,1])==i,]
  b2[i,3]=sum(as.numeric(index[,5]))
}
m2<-merge(n1,b2,by.x = "V1",by.y = "c2",all=T)
m2[is.na(m2[,"V3"]),"V3"]=0
m2$V2 <- gsub("/.*$*","",m2$c1)

#library(ggplot2)
#p <- ggplot()
#layer1 <- geom_point(data=m1,size =0.2,aes(x = V1,y = V3,colour = V2),shape="|")
#(p1 <- p  +layer1 + xlab("Start Pos") + ylab("base mutation number")
#+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
#+theme(plot.title=element_text(hjust = 0.5))
#+annotate("segment",x=m2$V1, xend=m2$V1+1, y=m2$V3, yend=m2$V3,colour="red",size=0.5)
#+annotate("text",x=m2$V1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
#p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
#          panel.background = element_blank(),axis.line = element_line(colour = "black"))
#ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_3.pdf",sep=""),width =8,height = 5)

library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1-min(m1[,1])+1,y = V3,colour = V2),shape="|")
(p1 <- p  +layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5))
+annotate("segment",x=m2$V1-min(m2[,1])+1, xend=m2$V1+1-min(m2[,1])+1, y=m2$V3, yend=m2$V3,colour="red",size=0.5)
+annotate("text",x=m2$V1-min(m2[,1])+1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_4.pdf",sep=""),width =8,height = 5)



##其余16种cancer#散点图标颜色有点突变、有段突变的\有点突变、没有段突变的
path2="F:/tp53/TCGA16cancer/process_result/" 
cancer="read" 
name<-list.files(path2,pattern = paste(cancer,'.*_TP53.txt',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=F,sep="\t")
aa<-as.data.frame(cbind(paste(as.character(f[,5]),as.character(f[,13]),sep="/"),as.character(f[,13]),as.character(f[,15])))
aa[,4]<-0
for(i in 1:dim(aa)[1]){
  if(!"-" %in% aa[i,3]){
    aa[i,4]<-nchar(as.character(aa[i,3]))
  }
}
e1<-subset(aa,aa[,4]<=1) #base<=1
e2<-subset(aa,aa[,4]>1)  #base>1
b1<-data.frame(c1=as.character(unique(e1[,1])))
row.names(b1)<-b1[,1]
b1$c2 <- gsub("^.*/","",b1$c1)
for(i in as.character(b1[,1])){
  index=e1[as.character(e1[,1])==i,]
  b1[i,3]=sum(as.numeric(index[,4]))
}
n1<-as.data.frame(matrix(0,max(as.numeric(as.character(f[,13])))-min(as.numeric(as.character(f[,13])))+1,2))
rownames(n1)<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
n1[,1]<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
m1<-merge(n1,b1,by.x = "V1",by.y = "c2",all=T)
m1[is.na(m1[,"V3"]),"V3"]=0
m1$V2 <- gsub("/.*$*","",m1$c1)

library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1,y = V3,colour = V2),shape="|")
(p1 <- p + layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5)))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base_3.pdf",sep=""),width =8,height = 5)

library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1-min(m1[,1])+1,y = V3,colour = V2),shape="|")
(p1 <- p + layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5)))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base_4.pdf",sep=""),width =8,height = 5)


library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1,y = V3,colour = V2),shape="|")
(p1 <- p  +layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5))
+annotate("segment",x=m2$V1, xend=m2$V1+1, y=m2$V3, yend=m2$V3,colour="red",size=0.5)
+annotate("text",x=m2$V1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base_3.pdf",sep=""),width =8,height = 5)

library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1-min(m1[,1])+1,y = V3,colour = V2),shape="|")
(p1 <- p  +layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5))
+annotate("segment",x=m2$V1-min(m2[,1])+1, xend=m2$V1+1-min(m2[,1])+1, y=m2$V3, yend=m2$V3,colour="red",size=0.5)
+annotate("text",x=m2$V1-min(m2[,1])+1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path2,"散点图/",strsplit(name,"_")[[1]][1],"_base_4.pdf",sep=""),width =8,height = 5)




 
#b\散点图标颜色没有段突变，只有点突变<=1cesc\dlbc\meso\tgct\thym\
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="thym" 
name<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=T,sep=",")
aa<-as.data.frame(cbind(paste(as.character(f[,5]),as.character(f[,13]),sep="/"),as.character(f[,13]),as.character(f[,15])))
aa[,4]<-0
for(i in 1:dim(aa)[1]){
  if(!"-" %in% aa[i,3]){
    aa[i,4]<-nchar(as.character(aa[i,3]))
  }
}
e1<-subset(aa,aa[,4]<=1) #base<=1
e2<-subset(aa,aa[,4]>1)  #base>1
b1<-data.frame(c1=as.character(unique(e1[,1])))
row.names(b1)<-b1[,1]
b1$c2 <- gsub("^.*/","",b1$c1)
for(i in as.character(b1[,1])){
  index=e1[as.character(e1[,1])==i,]
  b1[i,3]=sum(as.numeric(index[,4]))
}
n1<-as.data.frame(matrix(0,max(as.numeric(as.character(f[,13])))-min(as.numeric(as.character(f[,13])))+1,2))
rownames(n1)<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
n1[,1]<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
m1<-merge(n1,b1,by.x = "V1",by.y = "c2",all=T)
m1[is.na(m1[,"V3"]),"V3"]=0
m1$V2 <- gsub("/.*$*","",m1$c1)

#library(ggplot2)
#p <- ggplot()
#layer1 <- geom_point(data=m1,size =0.2,aes(x = V1,y = V3,colour = V2),shape="|")
#(p1 <- p + layer1 + xlab("Start Pos") + ylab("base mutation number")
#+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
#+theme(plot.title=element_text(hjust = 0.5)))
#p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
#          panel.background = element_blank(),axis.line = element_line(colour = "black"))
#ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_3.pdf",sep=""),width =8,height = 5)

library(ggplot2)
p <- ggplot()
layer1 <- geom_point(data=m1,size =0.2,aes(x = V1-min(m1[,1])+1,y = V3,colour = V2),shape="|")
(p1 <- p + layer1 + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5)))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_4.pdf",sep=""),width =8,height = 5)




#c\散点图标颜色没有点突变，只有段突变>1\\pcpg\\
path1="F:/tp53/clinic_data_process_result/"
path2="F:/tp53/TCGA/"
cancer="pcpg" 
name<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f<-read.csv(paste(path2,name,sep=""),header=T,sep=",")
aa<-as.data.frame(cbind(paste(as.character(f[,5]),as.character(f[,13]),sep="/"),as.character(f[,13]),as.character(f[,15])))
aa[,4]<-0
for(i in 1:dim(aa)[1]){
  if(!"-" %in% aa[i,3]){
    aa[i,4]<-nchar(as.character(aa[i,3]))
  }
}
#e1<-subset(aa,aa[,4]<=1) #base<=1
#e2<-subset(aa,aa[,4]>1)  #base>1
n1<-as.data.frame(matrix(0,max(as.numeric(as.character(f[,13])))-min(as.numeric(as.character(f[,13])))+1,2))
rownames(n1)<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
n1[,1]<-min(as.numeric(as.character(f[,13]))):max(as.numeric(as.character(f[,13])))
e2<-subset(aa,aa[,4]>1)  #base>1
e2[,5]<-1
b2<-data.frame(c1=as.character(unique(e2[,1])))
row.names(b2)<-b2[,1]
b2$c2 <- gsub("^.*/","",b2$c1)
for(i in as.character(b2[,1])){
  index=e2[as.character(e2[,1])==i,]
  b2[i,3]=sum(as.numeric(index[,5]))
}
m2<-merge(n1,b2,by.x = "V1",by.y = "c2",all=T)
m2[is.na(m2[,"V3"]),"V3"]=0
m2$V2 <- gsub("/.*$*","",m2$c1)

#library(ggplot2)
#p <- ggplot()
#(p1 <- p + xlab("Start Pos") + ylab("base mutation number")
#+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
#+theme(plot.title=element_text(hjust = 0.5))
#+annotate("segment",x=m2$V1, xend=m2$V1+1 , y=m2$V3, yend=m2$V3,colour="red",size=0.5)
#+annotate("text",x=m2$V1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
#p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
#          panel.background = element_blank(),axis.line = element_line(colour = "black"))
#ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_3.pdf",sep=""),width =8,height = 5)

library(ggplot2)
p <- ggplot()
(p1 <- p + xlab("Start Pos") + ylab("base mutation number")
+ggtitle(paste(strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base mutation",sep=""))
+theme(plot.title=element_text(hjust = 0.5))
+annotate("segment",x=m2$V1-min(m2[,1])+1, xend=m2$V1+1-min(m2[,1])+1 , y=m2$V3, yend=m2$V3,colour="red",size=0.5)
+annotate("text",x=m2$V1-min(m2[,1])+1, y=m2$V3+0.5, colour="red",size=1,label=paste(substr(gsub("_.*$*","",m2$V2), 1,1),substr(gsub("^.*_","",m2$V2), 1,1),sep="")))
p1+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave(file=paste(path1,"散点图/",strsplit(name,"_")[[1]][1],"_",strsplit(name,"_")[[1]][3],"_base_4.pdf",sep=""),width =8,height = 5)

 


  
  

library(grid)
library(VennDiagram)
#f1<-read.csv("F:/tp53/clinic_data/brca_Provisional_study_view_clinical_data.txt",header=T,sep="\t")#1105
f1<-read.csv("F:/tp53/TCGA/lgg_tcga_Provisional_clinical_data.tsv",header=T,sep="\t") 
f2<-read.csv("F:/tp53/TCGA/lgggbm_tcga_Cell_clinical_data.tsv",header=T,sep="\t") 
#f2<-read.csv("F:/tp53/brca320.txt",header=T,sep="\t")
#cc<-gsub("-0.*$*","",as.character(f1[,1]))
#i<-intersect(f1[,1],f2[,1])
A=unique(as.character(f1[,1]))
B=unique(as.character(f2[,1]))
pdf(paste("F:/tp53/","lgggbm_lgg.pdf",sep=""))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
dev.off()

f1<-read.csv("F:/tp53/TCGA/stes_tcga_Nature_clinical_data.tsv",header=T,sep="\t") 
f2<-read.csv("F:/tp53/TCGA/stad_tcga_Nature_clinical_data.tsv",header=T,sep="\t") 
f3<-read.csv("F:/tp53/TCGA/stad_tcga_Provisional_clinical_data.tsv",header=T,sep="\t") 
A=unique(as.character(f1[,1]))
B=unique(as.character(f2[,1]))
C=unique(as.character(f3[,1]))
pdf(paste("F:/tp53/","stes_stad.pdf",sep=""))
D<-venn.diagram(list(A=A,B=B,C=C),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
grid.draw(D)
dev.off()



###十六、从TCGA下载其余16种cancer的mutation及tp53 mutation对应的样本号的json文件，进行提取样本信息

#install.packages("rjson")
library("rjson")
path3="F:/tp53/TCGA16cancer/"
fname<-list.files(paste(path3,"data/",sep=""),pattern = '.*.json')
for(name in fname){
  f <- fromJSON(file = paste(path3,"data/",name,sep=""))
  df <- t(as.data.frame(f))
  write.table(df[grep("^submitter_id.*",row.names(df)),1],file =paste(path3,"process/",gsub("json","txt",name),sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}


#a、mutation有2个研究m_tp53都是大包小\blca\thca\lusc\luad\ucec\prad\laml\
library(grid)
library(VennDiagram)
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="laml" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
#f1<-read.csv(paste(PATH,"blca_Nature_study_view_clinical_data.txt",sep=""),header=T,sep="\t")
#f2<-read.csv(paste(PATH,"blca_Provisional_study_view_clinical_data.txt",sep=""),header=T,sep="\t")
name<-list.files(PATH,pattern = paste(cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste(PATH,name[2],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,4]))
B=unique(as.character(f2[,4]))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
i=intersect(f2[,4],f[,1])
a=f2[f2[,4] %in% i & as.character(f2[,8])=="YES",]
write.table(a,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
f5<-read.csv(paste(path2,name1[2],sep=""),header=T,sep=",")
E=unique(as.character(f4[,1]))
G=unique(as.character(f5[,1]))
D<-venn.diagram(list(A=E,B=G),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
j=intersect(gsub("-01|-03","",as.character(f5[,1])),f3[,1])
b=f5[gsub("-01|-03","",as.character(f5[,1])) %in% j,]
write.table(b,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


#b、mutation有2个研究m是大包小tp53相交取并\kirc\hnsc\stad\
library(grid)
library(VennDiagram)
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="stad" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
name<-list.files(PATH,pattern = paste("^",cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste(PATH,name[2],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,4]))
B=unique(as.character(f2[,4]))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
i=intersect(f2[,4],f[,1])
a=f2[f2[,4] %in% i & as.character(f2[,8])=="YES",]
write.table(a,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste("^",cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
f5<-read.csv(paste(path2,name1[2],sep=""),header=T,sep=",")
E=unique(as.character(f4[,1]))
G=unique(as.character(f5[,1]))
D<-venn.diagram(list(A=E,B=G),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
df=rbind(f4,f5)
j=intersect(gsub("-01","",as.character(df[,1])),f3[,1])
b=df[gsub("-01","",as.character(df[,1])) %in% j,]
write.table(b,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


#c、mutation有2个研究m相交取并tp53大包小\coad\read\
library(grid)
library(VennDiagram)
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="read" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
name<-list.files(PATH,pattern = paste("coad",cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste(PATH,name[2],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,4]))
B=unique(as.character(f2[,4]))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
i1=unique(f2[,4])
i2=setdiff(f1[,4],f2[,4])
i3=intersect(i1,f[,1])
i4=intersect(i2,f[,1])
a1=f2[f2[,4] %in% i3 & as.character(f2[,8])=="YES",]
a2=f1[f1[,4] %in% i4 & as.character(f1[,8])=="YES",]
write.table(a1,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(a2,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste("coad",cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
f5<-read.csv(paste(path2,name1[2],sep=""),header=T,sep=",")
E=unique(as.character(f4[,1]))
G=unique(as.character(f5[,1]))
D<-venn.diagram(list(A=E,B=G),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
j=intersect(gsub("-01","",as.character(f4[,1])),f3[,1])
b=f4[gsub("-01","",as.character(f4[,1])) %in% j,]
write.table(b,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

#d、mutation有2个研究m_tp53都是相交取并\ov\
library(grid)
library(VennDiagram)
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="ov" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
name<-list.files(PATH,pattern = paste("^",cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste(PATH,name[2],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,4]))
B=unique(as.character(f2[,4]))
D<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
i1=unique(f2[,4])
i2=setdiff(f1[,4],f2[,4])
i3=intersect(i1,f[,1])
i4=intersect(i2,f[,1])
a1=f2[f2[,4] %in% i3 & as.character(f2[,8])=="YES",]
a2=f1[f1[,4] %in% i4 & as.character(f1[,8])=="YES",]
write.table(a1,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(a2,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste("^",cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
f5<-read.csv(paste(path2,name1[2],sep=""),header=T,sep=",")
E=unique(as.character(f4[,1]))
G=unique(as.character(f5[,1]))
D<-venn.diagram(list(A=E,B=G),filename=NULL,lwd=1,lty=2,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),reverse=TRUE)
grid.draw(D)
df=rbind(f4,f5)
j=intersect(gsub("-01","",as.character(df[,1])),f3[,1])
b=df[gsub("-01","",as.character(df[,1])) %in% j,]
write.table(b,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


#e、mutation只有一个研究\sarc\lihc\
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="lihc" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
name<-list.files(PATH,pattern = paste(cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
i=intersect(f1[,4],f[,1])
a=f1[f1[,4] %in% i & as.character(f1[,8])=="YES",]
write.table(a,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
j=intersect(gsub("-01","",as.character(f4[,1])),f3[,1])
b=f4[gsub("-01","",as.character(f4[,1])) %in% j,]
write.table(b,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

 
#f、mutation有3个研究m相交取并tp53大包小\gbm\
library(grid)
library(VennDiagram)
PATH="F:/tp53/clinic_data/"
path2="F:/tp53/TCGA/"
path3="F:/tp53/TCGA16cancer/" 
cancer="gbm" 
f<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_m.txt",sep=""),header=F,sep="\t")
name<-list.files(PATH,pattern = paste(cancer,'.*.txt',sep="")) 
f1<-read.csv(paste(PATH,name[1],sep=""),header=T,sep="\t")
f2<-read.csv(paste(PATH,name[2],sep=""),header=T,sep="\t")
f22<-read.csv(paste(PATH,name[3],sep=""),header=T,sep="\t")
A=unique(as.character(f1[,4]))
B=unique(as.character(f2[,4]))
H=unique(as.character(f22[,4]))
D<-venn.diagram(list(A=A,B=B,C=H),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
grid.draw(D)
i1=setdiff(f1[,4],intersect(f1[,4],f2[,4])) 
i2=setdiff(f2[,4],setdiff(intersect(f2[,4],f22[,4]),intersect(intersect(f1[,4],f2[,4]),f22[,4]))) 
i3=setdiff(f22[,4],intersect(f1[,4],f22[,4]))
i4=intersect(i1,f[,1])
i5=intersect(i2,f[,1])
i6=intersect(i3,f[,1])
a1=f1[f1[,4] %in% i4 & as.character(f1[,8])=="YES",]
a2=f2[f2[,4] %in% i5 & as.character(f2[,8])=="YES",]
a3=f22[f22[,4] %in% i6 & as.character(f22[,8])=="YES",]
write.table(a1,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(a2,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(a3,file =paste(path3,"process_result/",cancer,"_M.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#tp53
f3<-read.csv(paste(path3,"process/","TCGA-",toupper(cancer),"_tp53.txt",sep=""),header=F,sep="\t") 
name1<-list.files(path2,pattern = paste(cancer,'.*_tp53.csv',sep="")) 
f4<-read.csv(paste(path2,name1[1],sep=""),header=T,sep=",")
f5<-read.csv(paste(path2,name1[2],sep=""),header=T,sep=",")
f55<-read.csv(paste(path2,name1[3],sep=""),header=T,sep=",")
E=unique(as.character(f4[,1]))
G=unique(as.character(f5[,1]))
M=unique(as.character(f55[,1]))
D<-venn.diagram(list(A=E,B=G,C=M),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
grid.draw(D)
j1=unique(gsub("-01","",as.character(f5[,1])))
j2=setdiff(gsub("-01","",as.character(f55[,1])),gsub("-01","",as.character(f5[,1])))
j3=intersect(j1,f3[,1])
j4=intersect(j2,f3[,1])
b1=f5[gsub("-01","",as.character(f5[,1])) %in% j3,]
b2=f55[gsub("-01","",as.character(f55[,1])) %in% j4,]
write.table(b1,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
write.table(b2,file =paste(path3,"process_result/",cancer,"_TP53.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



###十七、整合文件得到总表
##前11种cancer
a={}
path1="F:/tp53/clinic_data_process_result/cancer_M/"
#cancer="brca"
#fname<-list.files(path1,pattern = paste(cancer,'.*_M.txt',sep="")) 
fname<-list.files(path1,pattern = paste('.*_M.txt',sep="")) 
for(name in fname){
  cancer=gsub("_M.txt","",name)
  f<-read.csv(paste(path1,name,sep=""),header=F,sep="\t")
  a=append(a,paste(unique(as.character(f[,4])),cancer))
}
write.table(a,file ="F:/tp53/clinic_data_process_result/result1.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


f1<-read.csv(file ="F:/tp53/clinic_data_process_result/result1.txt",header=F,sep=" ")
f1[,3]=0
fcancer=unique(gsub("_","_tcga_",as.character(f1[,2])))
for(cancer in fcancer){
  f2<-read.csv(paste("F:/tp53/TCGA/",cancer,'_tp53.csv',sep=""),header=T,sep=",")
  f1[as.character(f1[,1]) %in% gsub("-01","",f2[,1]) & as.character(f1[,2])==gsub("_tcga_","_",cancer),3]=1
}
write.table(f1,file ="F:/tp53/clinic_data_process_result/result2.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

##后16种cancer
b={}
path2="F:/tp53/TCGA16cancer/process_result/"
fname<-list.files(path2,pattern = paste('.*_M.txt',sep="")) 
for(name in fname){
  cancer=gsub("_M.txt","",name)
  f<-read.csv(paste(path2,name,sep=""),header=F,sep="\t")
  b=append(b,paste(unique(as.character(f[,4])),cancer))
}
write.table(b,file ="F:/tp53/TCGA16cancer/result3.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

f1<-read.csv(file ="F:/tp53/TCGA16cancer/result3.txt",header=F,sep=" ")
f1[,3]=0
fcancer=unique(as.character(f1[,2]))
for(cancer in fcancer){
  f2<-read.csv(paste(path2,cancer,'_TP53.txt',sep=""),header=F,sep="\t")
  f1[as.character(f1[,1]) %in% gsub("-01","",f2[,1]) & as.character(f1[,2])==cancer,3]=1
}
write.table(f1,file ="F:/tp53/TCGA16cancer/result4.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

##共27种cancer
f1<-read.csv("F:/tp53/clinic_data_process_result/result2.txt",header=F,sep="\t")
f2<-read.csv("F:/tp53/TCGA16cancer/result4.txt",header=F,sep="\t")
#rbind(f1,f2)
write.table(rbind(f1,f2),file ="F:/tp53/TCGA16cancer/Result.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


f<-read.csv("F:/tp53/TCGA16cancer/Result.txt",header=F,sep="\t")
a0<-subset(f,as.character(f[,3])==0)
a1<-subset(f,as.character(f[,3])==1)
for(i in 1:nrow(a0)){
  if("_" %in% a0[i,2])
    f1<-read.csv(paste("F:/tp53/clinic_data_process_result/cancer_M/",a0[i,2],"_M.txt",sep=""),header=F,sep="\t")
  if(!"_" %in% a0[i,2])
    f1<-read.csv(paste("F:/tp53/TCGA16cancer/process_result/",a0[i,2],"_M.txt",sep=""),header=F,sep="\t")
}


 
###十八、TCGA10188样本
library("rjson")
f <- fromJSON(file="F:/tp53/TCGA10188/TCGA10188.json")
df <- t(as.matrix(f))
write.table(df,file ="F:/tp53/TCGA10188/TCGA10188.txt", row.names = F,col.names=F, quote = F,sep = "\n",append=TRUE)

library("rjson")
f <- fromJSON(file="F:/tp53/TCGA10188/TP53_3956.json")
df <- t(as.matrix(f))
write.table(df,file ="F:/tp53/TCGA10188/TP53_3956.txt", row.names = F,col.names=F, quote = F,sep = "\n",append=TRUE)


##mutation的样本是否发生tp53突变
f1<-read.csv(file ="F:/tp53/TCGA10188/TCGA10188_1.txt",header=F,sep="\t")
f2<-read.csv(file ="F:/tp53/TCGA10188/TP53_3956_1.txt",header=F,sep="\t")
f1[,12]=0
a=gsub("\"","",gsub("submitter_id=","",as.character(f1[,4])))
b=gsub("\"","",gsub("submitter_id=","",as.character(f2[,4])))
#a %in% b
c=gsub("\"","",gsub("project_id=","",as.character(f1[,1])))
d=gsub("\"","",gsub("project_id=","",as.character(f2[,1])))
#c==d
f1[gsub("\"","",gsub("submitter_id=","",as.character(f1[,4]))) %in% gsub("\"","",gsub("submitter_id=","",as.character(f2[,4]))) ,12]=1
write.table(f1,file ="F:/tp53/TCGA10188/Result1.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


f<-read.csv(file ="F:/tp53/TCGA10188/Result1.txt",header=T,sep="\t")
#cancer=unique(f[,1])
#race=na.omit(gsub("^$",NA,unique(f[,10])))
#ethnicity=c("nothispanicorlatino","hispanicorlatino")
#a=subset(f,as.character(f[,10])=="white"& as.character(f[,11])=="nothispanicorlatino")
#nrow(a[as.character(a[,1])=="HNSC",]) 
#sum(a[as.character(a[,1])=="HNSC",12])
for(i in na.omit(gsub("^$",NA,unique(f[,10])))){
  for(j in c("nothispanicorlatino","hispanicorlatino")){
    a=subset(f,as.character(f[,10])==i & as.character(f[,11])==j)
    for(k in unique(a[,1])){
      m=nrow(a[as.character(a[,1])==k,]) 
      n=sum(a[as.character(a[,1])==k,12])
      write.table(paste(i,j,k,m,n,sep="\t"),file ="F:/tp53/TCGA10188/race_ethnicity_cancer_M_TP53.txt", row.names = F,col.names=F, quote = F,sep = "\n",append=TRUE)
    }
  }
}


library("survival")
f<-read.csv(file ="F:/tp53/TCGA10188/Result1.txt",header=T,sep="\t",stringsAsFactors = F)
f$days_to_death<-as.numeric(gsub("NULL",NA,f$days_to_death))
f$vital_status<-gsub("alive",0,f$vital_status)
f$vital_status<-as.numeric(gsub("dead",1,f$vital_status))
dif <- survdiff(Surv(na.omit(f)$days_to_death,na.omit(f)$vital_status)~na.omit(f)$TP53_mutation)#求生存时间
#p_value <- pchisq(dif$chisq,1,lower.tail=F) #p= 0.133 
kmsurvival<-survfit(Surv(na.omit(f)$days_to_death,na.omit(f)$vital_status)~na.omit(f)$TP53_mutation,conf.type = "log-log")
#summary(kmsurvival)
plot(kmsurvival, lty = c('solid', 'dashed'), col=c('black','red'),
     xlab='survival time in days',ylab='survival probabilities')
legend('topright', c('no TP53 mutation','TP53 mutation'), lty=c('solid','dashed'),
       col=c('black','red'))


###十九、TCGA11160样本
library("rjson")
f <- fromJSON(file="F:/tp53/TCGA10188/metadata10before.json")
df <- t(as.matrix(f))
write.table(df,file ="F:/tp53/TCGA10188/metadata10before.txt", row.names = F,col.names=F, quote = F,sep = "\n",append=TRUE)
f1 <- fromJSON(file="F:/tp53/TCGA10188/metadata16after.json")
df1 <- t(as.matrix(f1))
write.table(df1,file ="F:/tp53/TCGA10188/metadata16after.txt", row.names = F,col.names=F, quote = F,sep = "\n",append=TRUE)


f1<-read.csv(file ="F:/tp53/TCGA10188/metadata10before_1.txt",header=F,sep="\t",stringsAsFactors = F)
for (i in 1:nrow(f1)){
  for (j in c(" entity_submitter_id = "," vital_status = "," days_to_death = "," days_to_birth = "," days_to_last_follow_up = ")){
    if(length(grep(j,f1[i,]))==1){
      k = f1[i,grep(j,f1[i,])]
      cat(k,"\t",file ="F:/tp53/TCGA10188/TCGA11160_1.txt",append=TRUE)
    }
  }
  if(length(grep("revision = ",f1[i,]))==1){
    cat(file ="F:/tp53/TCGA10188/TCGA11160_1.txt",sep = "\n",append=TRUE)
  }
}


f2<-read.csv(file ="F:/tp53/TCGA10188/metadata16after_1.txt",header=F,sep="\t",stringsAsFactors = F)
for (i in 1:nrow(f2)){
  for (j in c(" entity_submitter_id = "," vital_status = "," days_to_death = "," days_to_birth = "," days_to_last_follow_up = ")){
    if(length(grep(j,f2[i,]))==1){
      k = f2[i,grep(j,f2[i,])]
      cat(k,"\t",file ="F:/tp53/TCGA10188/TCGA11160_2.txt",append=TRUE)
    }
  }
  if(length(grep("revision = ",f2[i,]))==1){
    cat(file ="F:/tp53/TCGA10188/TCGA11160_2.txt",sep = "\n",append=TRUE)
  }
}


f1<-read.csv(file ="F:/tp53/TCGA10188/TCGA11151.txt",header=T,sep="\t")
f2<-read.csv(file ="F:/tp53/TCGA10188/Result1.txt",header=T,sep="\t")
a=subset(f1,f1$entity_submitter_id  %in% f2$submitter_id)
f=merge(f2,a,by.x = "submitter_id",by.y = "entity_submitter_id",all=T)
write.table(f,file ="F:/tp53/TCGA10188/Result2.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


library("survival")
##TP53 mutation or not 的生存曲线all cancer
pdf("F:/tp53/TCGA10188/all_survivalplot.pdf")
opar <- par(no.readonly=T) # 保存画图环境
par(mfrow=c(1,3))
par(oma=c(2,0,3,0))
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
a=subset(f,f$vital_status.y=="alive")
b=subset(f,f$vital_status.y=="dead")
alive=cbind(a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
dead=cbind(b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
e=rbind(alive,dead)
dif <- survdiff(Surv(e[,3],e[,2])~e[,1])#求生存时间
#p_value <- pchisq(dif$chisq,1,lower.tail=F)  
kmsurvival<-survfit(Surv(e[,3],e[,2])~e[,1],conf.type = "log-log")
#pdf("F:/tp53/TCGA10188/survivalPlot_TP53 mutation or not.pdf")
plot(kmsurvival, lty = c('solid', 'solid'), col=c("#12B826","#1DADB8"),lwd=2,
     xlab='survival time in days',ylab='survival probabilities')
legend('bottomleft', cex=0.8,text.width=0.6,c('NoTP53Mu','TP53Mu'), lty=c('solid','solid'),lwd=2,
       col=c("#12B826","#1DADB8"))
text(2000,0.9,cex=0.8,paste(paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""),dif$n[1],dif$n[2],sep="\n"))
#dev.off()

##TP53 not mutation and mutation race 的生存曲线 
#f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t") 
a=subset(f,f$vital_status.y=="alive" & f$race %in% c("asian","blackorafricanamerican","white"))
b=subset(f,f$vital_status.y=="dead" & f$race %in% c("asian","blackorafricanamerican","white"))
alive=cbind(a$race,a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
dead=cbind(b$race,b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
e=data.frame(rbind(alive,dead))
dif <- survdiff(Surv(e[,4],e[,3])~e[,2]+e[,1])#求生存时间
#p_value <- pchisq(dif$chisq,1,lower.tail=F)  
kmsurvival<-survfit(Surv(e[,4],e[,3])~e[,2]+e[,1],conf.type = "log-log")
#pdf("F:/tp53/TCGA10188/survivalPlot_TP53_race_6.pdf")
plot(kmsurvival, col=c("red","yellow","blue","green","orange","purple"),lty = c('solid', 'solid','solid','dashed', 'dashed', 'dashed'), 
     xlab='survival time in days',ylab='survival probabilities')
legend('bottomleft',cex=0.8,text.width=0.2,c("no asian","no blackorafricanamerican","no white","asian","blackorafricanamerican","white"),
       col=c("red","yellow","blue","green","orange","purple"),lty = c('solid','solid','solid', 'dashed', 'dashed', 'dashed'))
text(2200,1,cex=0.4,paste(paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""),dif$n[1],dif$n[2],dif$n[3],dif$n[4],dif$n[5],dif$n[6],sep="\n"))
#dev.off()

##TP53 not mutation and mutation race 的生存曲线 分开种族c("asian","blackorafricanamerican","white")
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t") 
a=subset(f,f$vital_status.y=="alive" & f$race %in% c("white"))
b=subset(f,f$vital_status.y=="dead" & f$race %in% c("white"))
alive=cbind(as.character(a$race),a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
dead=cbind(as.character(b$race),b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,4]=as.numeric(as.character(e[,4]))
e[,3]=as.numeric(as.character(e[,3]))
e[,2]=as.factor(e[,2])
e[,1]=as.factor(e[,1])
dif <- survdiff(Surv(e[,4],e[,3])~e[,2])#求生存时间
kmsurvival<-survfit(Surv(e[,4],e[,3])~e[,2],conf.type = "log-log")
pdf("F:/tp53/TCGA10188/survivalPlot_TP53_race_6_white.pdf")
plot(kmsurvival, col=c("mediumpurple2","mediumpurple2"),lty = c('solid','dashed'),lwd=2, 
     xlab='survival time in days',ylab='survival probabilities',main="white")
legend('bottomleft',cex=0.8,c("no white","white"),lwd=2,
       col=c("mediumpurple2","mediumpurple2"),lty = c('solid','dashed'))
text(2200,1,cex=1,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
dev.off()




#f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t") 
a=subset(f,f$vital_status.y=="alive" & f$race %in% c("asian","blackorafricanamerican","white"))
b=subset(f,f$vital_status.y=="dead" & f$race %in% c("asian","blackorafricanamerican","white"))
alive=cbind(as.character(a$race),a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
dead=cbind(as.character(b$race),b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[grep(0,e[,2]),1]="NoMutation"
e[,4]=as.numeric(as.character(e[,4]))
e[,3]=as.numeric(as.character(e[,3]))
e[,1]=as.factor(e[,1])
dif <- survdiff(Surv(e[,4],e[,3])~e[,1])#求生存时间
#p_value <- pchisq(dif$chisq,1,lower.tail=F)  
kmsurvival<-survfit(Surv(e[,4],e[,3])~e[,1],conf.type = "log-log")
#pdf("F:/tp53/TCGA10188/survivalPlot_TP53_race_4_1.pdf")
plot(kmsurvival, col=c("orange3","gray20","seagreen3","mediumpurple2"),lty = c( 'dashed', 'dashed','solid', 'dashed'),lwd=2, 
     xlab='survival time in days',ylab='survival probabilities',main="All cancer")
legend('bottomleft',cex=0.8,c("NoMutation","asian","blackorafricanamerican","white"),lwd=2,
       col=c("seagreen3","orange3","gray20","mediumpurple2"),lty = c('solid', 'dashed', 'dashed', 'dashed'))
#text(2200,1,cex=0.5,paste(paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""),dif$n[1],dif$n[2],dif$n[3],dif$n[4],sep="\n"))
text(2200,1,cex=1,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
mtext("all cancer", side = 3, outer = TRUE)
dev.off()

 
###cancer_survivalplot 单个cancer
library("survival")
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
#cancer=unique(f[,2])
for (cancer in unique(f[,2])){
  pdf(paste("F:/tp53/TCGA10188/cancer_survivalplot/",cancer,"_survivalplot.pdf",sep=""))
  opar <- par(no.readonly=T) # 保存画图环境
  par(mfrow=c(1,3))
  par(oma=c(2,0,3,0))
  a=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive=cbind(a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
  dead=cbind(b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
  e=rbind(alive,dead)
  dif <- survdiff(Surv(e[,3],e[,2])~e[,1])#求生存时间
  kmsurvival<-survfit(Surv(e[,3],e[,2])~e[,1],conf.type = "log-log")
  plot(kmsurvival, lty = c('solid', 'dashed'), col=c('black','red'),
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft', cex=0.6,text.width=0.4,c('no TP53 mutation','TP53 mutation'), lty=c('solid','dashed'),
         col=c('black','red'))
  text(2200,1,cex=0.8,paste(paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""),dif$n[1],dif$n[2],sep="\n"))
  
  a1=subset(f,f$vital_status.y=="alive" & f$race %in% c("asian","blackorafricanamerican","white")& f$project_id==cancer)
  b1=subset(f,f$vital_status.y=="dead" & f$race %in% c("asian","blackorafricanamerican","white")& f$project_id==cancer)
  alive1=cbind(a1$race,a1$TP53_mutation,a1$vital_status.y,a1$days_to_last_follow_up)
  dead1=cbind(b1$race,b1$TP53_mutation,b1$vital_status.y,b1$days_to_death.y)	
  e1=rbind(alive1,dead1)
  dif1 <- survdiff(Surv(e1[,4],e1[,3])~e1[,2]+e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(e1[,4],e1[,3])~e1[,2]+e1[,1],conf.type = "log-log")
  plot(kmsurvival1, col=c("red","yellow","blue","green","orange","purple"),lty = c('solid', 'solid','solid','dashed', 'dashed', 'dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no asian","no blackorafricanamerican","no white","asian","blackorafricanamerican","white"),
         col=c("red","yellow","blue","green","orange","purple"),lty = c('solid','solid','solid', 'dashed', 'dashed', 'dashed'))
  text(2200,1,cex=0.4,paste(paste("p=",pchisq(dif1$chisq,1,lower.tail=F),sep=""),dif1$n[1],dif1$n[2],dif1$n[3],dif1$n[4],dif1$n[5],dif1$n[6],sep="\n"))
  
  a2=subset(f,f$vital_status.y=="alive" & f$race %in% c("asian","blackorafricanamerican","white") & f$project_id==cancer)
  b2=subset(f,f$vital_status.y=="dead" & f$race %in% c("asian","blackorafricanamerican","white") & f$project_id==cancer)
  alive2=cbind(a2$race,a2$TP53_mutation,a2$vital_status.y,a2$days_to_last_follow_up)
  dead2=cbind(b2$race,b2$TP53_mutation,b2$vital_status.y,b2$days_to_death.y)	
  e2=rbind(alive2,dead2)
  e2[grep(0,e2[,2]),1]=5
  dif2 <- survdiff(Surv(e2[,4],e2[,3])~e2[,2]+e2[,1])#求生存时间
  kmsurvival2<-survfit(Surv(e2[,4],e2[,3])~e2[,2]+e2[,1],conf.type = "log-log")
  plot(kmsurvival2, col=c("red","blue","green","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no mutation","asian","blackorafricanamerican","white"),
         col=c("red","blue","green","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed'))
  text(2200,1,cex=0.5,paste(paste("p=",pchisq(dif2$chisq,1,lower.tail=F),sep=""),dif2$n[1],dif2$n[2],dif2$n[3],dif2$n[4],sep="\n"))
  mtext(cancer, side = 3, outer = TRUE)
  dev.off()
}


###cancer_survivalplot 单个cancer在不同年龄段间
library("survival")
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
#cancer=unique(f[,2])
for (cancer in unique(f[,2])){
  pdf(paste("F:/tp53/TCGA10188/cancer_survivalplot_age/",cancer,"_survivalplot.pdf",sep=""))
  opar <- par(no.readonly=T) # 保存画图环境
  par(mfrow=c(1,3))
  par(oma=c(2,0,3,0))
  a=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive=cbind(a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
  dead=cbind(b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
  e=as.data.frame(na.omit(rbind(alive,dead)))
  dif <- survdiff(Surv(e[,3],e[,2])~e[,1])#求生存时间
  kmsurvival<-survfit(Surv(e[,3],e[,2])~e[,1],conf.type = "log-log")
  plot(kmsurvival, lty = c('solid', 'dashed'), col=c('black','red'),
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft', cex=0.6,text.width=0.4,c('no TP53 mutation','TP53 mutation'), lty=c('solid','dashed'),
         col=c('black','red'))
  text(2200,1,cex=0.8,paste(paste("p=",round(pchisq(dif$chisq,1,lower.tail=F),5),sep=""),paste(dif$n[1],dif$n[2],sep="_"),sep="\n"))
  
  a1=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b1=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive1=cbind(a1$age_at_diagnosis,a1$TP53_mutation,a1$vital_status.y,a1$days_to_last_follow_up)
  dead1=cbind(b1$age_at_diagnosis,b1$TP53_mutation,b1$vital_status.y,b1$days_to_death.y)	
  e1=as.data.frame(na.omit(rbind(alive1,dead1)))
  for(i in 1:nrow(e1)){
    if(as.numeric(e1[i,1]/365)< 35){
      e1[i,5]=1
    }else if(as.numeric(e1[i,1]/365)>=35 & as.numeric(e1[i,1]/365)<50){
      e1[i,5]=2
    }else if(as.numeric(e1[i,1]/365)>=50 & as.numeric(e1[i,1]/365)<70){
      e1[i,5]=3
    }else if(as.numeric(e1[i,1]/365)>=70){
      e1[i,5]=4
    }
  }
  dif1 <- survdiff(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5])#求生存时间
  kmsurvival1<-survfit(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5],conf.type = "log-log")
  plot(kmsurvival1, col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no<35","no35~50","no50~70","no>70","<35","35~50","50~70",">70"),
         col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'))
  text(2200,1,cex=0.4,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),5),sep=""),
                            paste(dif1$n[1],dif1$n[2],dif1$n[3],dif1$n[4],sep="_"),
                            paste(dif1$n[5],dif1$n[6],dif1$n[7],dif1$n[8],sep="_"),sep="\n"))
  
  a2=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b2=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive2=cbind(a2$age_at_diagnosis,a2$TP53_mutation,a2$vital_status.y,a2$days_to_last_follow_up)
  dead2=cbind(b2$age_at_diagnosis,b2$TP53_mutation,b2$vital_status.y,b2$days_to_death.y)	
  e2=as.data.frame(na.omit(rbind(alive2,dead2)))
  for(i in 1:nrow(e2)){
    if(as.numeric(e2[i,1]/365)< 35){
      e2[i,5]=1
    }else if(as.numeric(e2[i,1]/365)>=35 & as.numeric(e2[i,1]/365)<50){
      e2[i,5]=2
    }else if(as.numeric(e2[i,1]/365)>=50 & as.numeric(e2[i,1]/365)<70){
      e2[i,5]=3
    }else if(as.numeric(e2[i,1]/365)>=70){
      e2[i,5]=4
    }
  }
  e2[grep(0,e2[,2]),5]=5
  dif2 <- survdiff(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5])#求生存时间
  kmsurvival2<-survfit(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5],conf.type = "log-log")
  plot(kmsurvival2, col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no mutation","<35","35~50","50~70",">70"),
         col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'))
  text(2200,1,cex=0.5,paste(paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),5),sep=""),paste(dif2$n[1],dif2$n[2],dif2$n[3],dif2$n[4],dif2$n[5],sep="_"),sep="\n"))
  mtext(cancer, side = 3, outer = TRUE)
  dev.off()
}


###cancer_survivalplot all cancer在不同年龄段间
library("survival")
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
pdf(paste("F:/tp53/TCGA10188/cancer_survivalplot_age/all_survivalplot.pdf",sep=""))
opar <- par(no.readonly=T) # 保存画图环境
par(mfrow=c(1,3))
par(oma=c(2,0,3,0))
a=subset(f,f$vital_status.y=="alive")
b=subset(f,f$vital_status.y=="dead")
alive=cbind(a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
dead=cbind(b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
e=as.data.frame(na.omit(rbind(alive,dead)))
dif <- survdiff(Surv(e[,3],e[,2])~e[,1])#求生存时间
kmsurvival<-survfit(Surv(e[,3],e[,2])~e[,1],conf.type = "log-log")
plot(kmsurvival, lty = c('solid', 'dashed'), col=c('black','red'),
     xlab='survival time in days',ylab='survival probabilities')
legend('bottomleft', cex=0.6,text.width=0.4,c('no TP53 mutation','TP53 mutation'), lty=c('solid','dashed'),col=c('black','red'))
text(2200,1,cex=0.8,paste(paste("p=",round(pchisq(dif$chisq,1,lower.tail=F),5),sep=""),paste(dif$n[1],dif$n[2],sep="_"),sep="\n"))
a1=subset(f,f$vital_status.y=="alive" )
b1=subset(f,f$vital_status.y=="dead")
alive1=cbind(a1$age_at_diagnosis,a1$TP53_mutation,a1$vital_status.y,a1$days_to_last_follow_up)
dead1=cbind(b1$age_at_diagnosis,b1$TP53_mutation,b1$vital_status.y,b1$days_to_death.y)	
e1=as.data.frame(na.omit(rbind(alive1,dead1)))
for(i in 1:nrow(e1)){
  if(as.numeric(e1[i,1]/365)< 35){
    e1[i,5]=1
  }else if(as.numeric(e1[i,1]/365)>=35 & as.numeric(e1[i,1]/365)<50){
    e1[i,5]=2
  }else if(as.numeric(e1[i,1]/365)>=50 & as.numeric(e1[i,1]/365)<70){
    e1[i,5]=3
  }else if(as.numeric(e1[i,1]/365)>=70){
    e1[i,5]=4
  }
}
dif1 <- survdiff(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5])#求生存时间
kmsurvival1<-survfit(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5],conf.type = "log-log")
plot(kmsurvival1, col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'), 
      xlab='survival time in days',ylab='survival probabilities')
legend('bottomleft',cex=0.6,text.width=0.4,c("no<35","no35~50","no50~70","no>70","<35","35~50","50~70",">70"),
        col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'))
text(2200,1,cex=0.4,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),5),sep=""),
                          paste(dif1$n[1],dif1$n[2],dif1$n[3],dif1$n[4],sep="_"),
                          paste(dif1$n[5],dif1$n[6],dif1$n[7],dif1$n[8],sep="_"),sep="\n"))
a2=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
b2=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
alive2=cbind(a2$age_at_diagnosis,a2$TP53_mutation,a2$vital_status.y,a2$days_to_last_follow_up)
dead2=cbind(b2$age_at_diagnosis,b2$TP53_mutation,b2$vital_status.y,b2$days_to_death.y)	
e2=as.data.frame(na.omit(rbind(alive2,dead2)))
for(i in 1:nrow(e2)){
  if(as.numeric(e2[i,1]/365)< 35){
    e2[i,5]=1
  }else if(as.numeric(e2[i,1]/365)>=35 & as.numeric(e2[i,1]/365)<50){
    e2[i,5]=2
  }else if(as.numeric(e2[i,1]/365)>=50 & as.numeric(e2[i,1]/365)<70){
    e2[i,5]=3
  }else if(as.numeric(e2[i,1]/365)>=70){
    e2[i,5]=4
  }
}
e2[grep(0,e2[,2]),5]=5
dif2 <- survdiff(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5])#求生存时间
kmsurvival2<-survfit(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5],conf.type = "log-log")
plot(kmsurvival2, col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'), 
      xlab='survival time in days',ylab='survival probabilities')
legend('bottomleft',cex=0.6,text.width=0.4,c("no mutation","<35","35~50","50~70",">70"),
        col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'))
text(2200,1,cex=0.5,paste(paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),5),sep=""),paste(dif2$n[1],dif2$n[2],dif2$n[3],dif2$n[4],dif2$n[5],sep="_"),sep="\n"))
mtext("all_cancer", side = 3, outer = TRUE)
dev.off()



###cancer_survivalplot 在不同年龄段间每个cancer突变、没突变样本的生存
library("survival")
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
#cancer=unique(f[,2])
for (cancer in unique(f[,2])){
  pdf(paste("F:/tp53/TCGA10188/cancer_survivalplot_age/",cancer,"_survivalplot.pdf",sep=""))
  opar <- par(no.readonly=T) # 保存画图环境
  par(mfrow=c(1,3))
  par(oma=c(2,0,3,0))
  a=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive=cbind(a$TP53_mutation,a$vital_status.y,a$days_to_last_follow_up)
  dead=cbind(b$TP53_mutation,b$vital_status.y,b$days_to_death.y)	
  e=as.data.frame(na.omit(rbind(alive,dead)))
  dif <- survdiff(Surv(e[,3],e[,2])~e[,1])#求生存时间
  kmsurvival<-survfit(Surv(e[,3],e[,2])~e[,1],conf.type = "log-log")
  plot(kmsurvival, lty = c('solid', 'dashed'), col=c('black','red'),
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft', cex=0.6,text.width=0.4,c('no TP53 mutation','TP53 mutation'), lty=c('solid','dashed'),
         col=c('black','red'))
  text(2200,1,cex=0.8,paste(paste("p=",round(pchisq(dif$chisq,1,lower.tail=F),5),sep=""),paste(dif$n[1],dif$n[2],sep="_"),sep="\n"))
  
  a1=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b1=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive1=cbind(a1$age_at_diagnosis,a1$TP53_mutation,a1$vital_status.y,a1$days_to_last_follow_up)
  dead1=cbind(b1$age_at_diagnosis,b1$TP53_mutation,b1$vital_status.y,b1$days_to_death.y)	
  e1=as.data.frame(na.omit(rbind(alive1,dead1)))
  for(i in 1:nrow(e1)){
    if(as.numeric(e1[i,1]/365)< 35){
      e1[i,5]=1
    }else if(as.numeric(e1[i,1]/365)>=35 & as.numeric(e1[i,1]/365)<50){
      e1[i,5]=2
    }else if(as.numeric(e1[i,1]/365)>=50 & as.numeric(e1[i,1]/365)<70){
      e1[i,5]=3
    }else if(as.numeric(e1[i,1]/365)>=70){
      e1[i,5]=4
    }
  }
  dif1 <- survdiff(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5])#求生存时间
  kmsurvival1<-survfit(Surv(e1[,4],e1[,3])~e1[,2]+e1[,5],conf.type = "log-log")
  plot(kmsurvival1, col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no<35","no35~50","no50~70","no>70","<35","35~50","50~70",">70"),
         col=c("red","yellow","blue","black","pink","green","orange","purple"),lty = c('solid','solid','solid','solid','dashed','dashed','dashed','dashed'))
  text(2200,1,cex=0.4,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),5),sep=""),
                            paste(dif1$n[1],dif1$n[2],dif1$n[3],dif1$n[4],sep="_"),
                            paste(dif1$n[5],dif1$n[6],dif1$n[7],dif1$n[8],sep="_"),sep="\n"))
  
  a2=subset(f,f$vital_status.y=="alive" & f$project_id==cancer)
  b2=subset(f,f$vital_status.y=="dead" & f$project_id==cancer)
  alive2=cbind(a2$age_at_diagnosis,a2$TP53_mutation,a2$vital_status.y,a2$days_to_last_follow_up)
  dead2=cbind(b2$age_at_diagnosis,b2$TP53_mutation,b2$vital_status.y,b2$days_to_death.y)	
  e2=as.data.frame(na.omit(rbind(alive2,dead2)))
  for(i in 1:nrow(e2)){
    if(as.numeric(e2[i,1]/365)< 35){
      e2[i,5]=1
    }else if(as.numeric(e2[i,1]/365)>=35 & as.numeric(e2[i,1]/365)<50){
      e2[i,5]=2
    }else if(as.numeric(e2[i,1]/365)>=50 & as.numeric(e2[i,1]/365)<70){
      e2[i,5]=3
    }else if(as.numeric(e2[i,1]/365)>=70){
      e2[i,5]=4
    }
  }
  e2[grep(0,e2[,2]),5]=5
  dif2 <- survdiff(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5])#求生存时间
  kmsurvival2<-survfit(Surv(e2[,4],e2[,3])~e2[,2]+e2[,5],conf.type = "log-log")
  plot(kmsurvival2, col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'), 
       xlab='survival time in days',ylab='survival probabilities')
  legend('bottomleft',cex=0.6,text.width=0.4,c("no mutation","<35","35~50","50~70",">70"),
         col=c("black","pink","green","orange","purple"),lty = c('solid', 'dashed', 'dashed', 'dashed', 'dashed'))
  text(2200,1,cex=0.5,paste(paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),5),sep=""),paste(dif2$n[1],dif2$n[2],dif2$n[3],dif2$n[4],dif2$n[5],sep="_"),sep="\n"))
  mtext(cancer, side = 3, outer = TRUE)
  dev.off()
}



##tp53在不同年龄段间的突变率
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
e=as.data.frame(na.omit(cbind(f$age_at_diagnosis,f$TP53_mutation)))
e$V1=e$V1/365
e2=subset(e, e$V2==1)
b=table(cut(e2$V1,15)) #分为15组
#matrix(b)
pdf(paste("F:/tp53/TCGA10188/","age_at_diagnosis_tp53_1.pdf",sep=""))
f2<-na.omit(as.data.frame(f$age_at_diagnosis))
f2=f2/365
d=table(cut(f2$`f$age_at_diagnosis`,15))
#lines(d,col="red",type="b")
plot(d,col="red",type="b",cex=0.1,xlab='age_at_diagnosis in years',ylab='TP53 mutation frequency',main="age_at_diagnosis_tp53")
text(1:15,matrix(d)[,1]+10,round(matrix(d)[,1]/10188,3)) #分为15组
text(1:15,matrix(d)[,1]+50,matrix(d)[,1])
#plot(b,col="blue",type="b",cex=0.2,xlab='age_at_diagnosis in years',ylab='TP53 mutation frequency')
lines(b,col="blue",type="b",cex=0.1)
text(1:15,matrix(b)[,1]+10,round(matrix(b)[,1]/3956,3)) #分为15组
text(1:15,matrix(b)[,1]-40,matrix(b)[,1]) 
#mtext("age_at_diagnosis_tp53", side = 3, outer = FALSE)
legend('topleft',cex=0.6,text.width=0.4,c("all","tp53"),col=c("red","blue"),lty=2)
dev.off()
#ks.test(matrix(b)[,1],matrix(d)[,1])  #不显著 同一分布

##tp53在不同年龄段间的突变率2
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
e=as.data.frame(na.omit(cbind(f$age_at_diagnosis,f$TP53_mutation)))
e$V1=e$V1/365
e2=subset(e, e$V2==1)
b=table(cut(e2$V1,15)) #分为15组
#pdf(paste("F:/tp53/TCGA10188/","age_at_diagnosis_tp53_2.pdf",sep=""))
pdf(paste("F:/tp53/TCGA10188/","age_at_diagnosis_tp53_22.pdf",sep=""))
e3=subset(e, e$V2==0)
d=table(cut(e3$V1,15))
plot(d/sum(d),col="#12B826",type="b",cex=0.1,lwd=2,xlab='age_at_diagnosis in years',ylab='TP53 mutation frequency',main="age_at_diagnosis_tp53")
text(1:15,(matrix(d)[,1]+10)/sum(d),round(matrix(d)[,1]/10188,3)) #分为15组
text(1:15,(matrix(d)[,1]+50)/sum(d),matrix(d)[,1])
lines(b/sum(b),col="#1DADB8",type="b",lwd=2,cex=0.1)
text(1:15,(matrix(b)[,1]+10)/sum(b),round(matrix(b)[,1]/3956,3)) #分为15组
text(1:15,(matrix(b)[,1]-40)/sum(b),matrix(b)[,1]) 
mtext(paste("p=",signif(ks.test(matrix(b)[,1],matrix(d)[,1])$p.value,3),sep=""), side = 3, outer = FALSE)
legend('topleft',cex=0.6,text.width=0.4,c("NoMu","TP53Mu"),col=c("#12B826","#1DADB8"),lty=2)
dev.off()
#ks.test(matrix(b)[,1],matrix(d)[,1])$p.value  #不显著 同一分布


##tp53在不同年龄段间的突变率3 boxplot\Allcancer\age_at_diagnosis正数
library(ggplot2)
f<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t")
f$age_at_diagnosis<-as.numeric(gsub("NULL",NA,f$age_at_diagnosis))
e=as.data.frame(na.omit(cbind(f$age_at_diagnosis,f$TP53_mutation)))
e$V1=e$V1/365
e[,2]=as.factor(e[,2])
e1=subset(e, e$V2==0)
e2=subset(e, e$V2==1)
p<-ggplot(data=e, aes(x = e[,2],y = e[,1]))+geom_boxplot(aes(fill=V2,color=V2),outlier.shape = ".",lwd=0.6,alpha=0.5)+
  scale_fill_manual(values=c("#12B826","#1DADB8"))+scale_color_manual(values=c("#12B826","#1DADB8"))
p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Different Samples") + ylab("Age in years") + ggtitle("Age Distribution in Different Samples")  
p <- p+ theme(legend.position="top")+
  annotate("text",x=1.5, y=1, colour="black",size=4,label=paste("p=",signif(t.test(e1[,1],e2[,1])$p.value,3),sep=""))
ggsave(file=paste("F:/tp53/TCGA10188/","age_at_diagnosis_tp53_3.pdf",sep=""),width = 10,height=9)
#outlier.colour = NA,outlier.shape = ".",outlier.size=0



 

###十九、CBioportal_TP53(提取tp53突变的样本对应的protein change突变信息)
path="F:/tp53/CBioportal_TP53/"
#dir.create(paste(path,"Analysis",sep="")) 
path1="F:/tp53/CBioportal_TP53/Analysis/"
cancer="coadread"
dir.create(paste(path1,cancer,sep=""))
path2=paste(path1,cancer,"/",sep="")

f<-read.csv(paste(path,"Result2_3956.txt",sep=""),header=T,sep="\t")
fname<-list.files(path,pattern = paste("^",cancer,sep="")) 
f1<-read.csv(paste(path,fname[1],"/",cancer,"/tcga/data_mutations_extended.txt",sep=""),header=T,sep="\t",quote="") ##quote为空，文件里不识别的字符不读入文件，可正常读入文件。
f2<-read.csv(paste(path,fname[2],"/",cancer,"_tcga_pub/data_mutations_extended.csv",sep=""),header=T,sep=",")
#f3<-read.csv(paste(path,fname[3],"/",cancer,"_tcga_pub2013/data_mutations_extended.csv",sep=""),header=T,sep=",") 

#sample_id=subset(f,as.character(f[,2])==toupper(cancer))[,1]
a=subset(f,as.character(f[,2])==toupper("coad"))[,1]
b=subset(f,as.character(f[,2])==toupper("read"))[,1]
sample_id=c(as.character(a),as.character(b))
for( i in sample_id){
  dir.create(paste(path2,i,sep=""))
  c1=subset(f1,substring(f1[,17],1,12)==i)
  c2=subset(f2,substring(f2[,17],1,12)==i)
  #c3=subset(f3,substring(f3[,17],1,12)==i)
  write.table(c1,file =paste(path2,i,"/",fname[1],"_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(c2,file =paste(path2,i,"/",fname[2],"_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  #write.table(c3,file =paste(path2,i,"/",fname[3],"_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(c1,file =paste(path2,i,"/","total_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(c2,file =paste(path2,i,"/","total_mutation.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  #write.table(c3,file =paste(path2,i,"/","total_mutation.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}

 
path="F:/tp53/CBioportal_TP53/" 
path1="F:/tp53/CBioportal_TP53/Analysis/"
f<-read.csv(paste(path,"Result2_3956.txt",sep=""),header=T,sep="\t")
#a=subset(f,as.character(f[,2])==toupper("coad"))
#b=subset(f,as.character(f[,2])==toupper("read"))
#f=rbind(a,b)
##tolower(f[i,2])="coadread"
for (i in 1:nrow(f)){
  fname<-list.files(paste(path1,tolower(f[i,2]),"/",f[i,1],"/",sep=""),pattern = paste("^",tolower(f[i,2]),sep="")) 
  b={}
  for(j in fname){
    f1<-read.csv(paste(path1,tolower(f[i,2]),"/",f[i,1],"/",j,sep=""),header=T,sep="\t")
    b=append(b,paste(gsub("_mutation.txt","",j),nrow(f1),sep=","))
  }
  temp=merge(f[i,],t(as.data.frame(b)),all.x=T)
  write.table(temp,file =paste(path,"Result2_3956_gene1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}


source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)
# 下载数据
# Gene HTSeq - Counts (GRCh37.p13)
cancer <- c("BLCA","BRCA","CESC","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
for (k in cancer){
query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
data.category = "Transcriptome Profiling",  #-········· 数据种类
data.type = "Gene Expression Quantification",  #······· 数据类型
workflow.type = "HTSeq - Counts")  #··················· 工作流类型
GDCdownload(query,method = "client")  #·························································· GDC下载数据
}
setwd("D:/ProgramFilesRRR")
query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf",
                           legacy = TRUE)
GDCdownload(query.maf.hg19)

f1<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/gdc_download_20171109_070146/bbad60f0-3b1a-46e2-93d6-0c2bc41bcf23/1.txt",header=T,sep="\t")
f2<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/gdc_download_20171109_070229/a3d33307-c701-416c-b893-c254530c7e59/2.txt",header=T,sep="\t")
f3<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/gdc_download_20171109_070237/2c89ba82-312a-4d81-a02d-dc15c6a4293c/3.txt",header=T,sep="\t")
f4<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/gdc_download_20171109_070244/aead88d6-f87d-4c4a-8a03-5184de8d8e9d/4.txt",header=T,sep="\t")
#head(f1[,1])
a1=f1[grep("^TCGA-C5-A1M6.*",f1$Tumor_Sample_Barcode),]
a2=f2[grep("^TCGA-C5-A1M6.*",f2$Tumor_Sample_Barcode),]
a3=f3[grep("^TCGA-C5-A1M6.*",f3$Tumor_Sample_Barcode),]
a4=f4[grep("^TCGA-C5-A1M6.*",f4$Tumor_Sample_Barcode),]
length(unique(union(union(union(as.character(a1[,1]),as.character(a2[,1])),as.character(a3[,1])),as.character(a4[,1]))))
#194
length(unique(union(as.character(a1[,1]),as.character(a2[,1]))))#192
length(unique(union(as.character(a1[,1]),as.character(a3[,1]))))#170
length(unique(union(as.character(a1[,1]),as.character(a4[,1]))))#176
length(unique(union(as.character(a2[,1]),as.character(a3[,1]))))#185
length(unique(union(as.character(a2[,1]),as.character(a4[,1]))))#188
length(unique(union(as.character(a3[,1]),as.character(a4[,1]))))#164
b=unique(union(union(union(as.character(a1[,1]),as.character(a2[,1])),as.character(a3[,1])),as.character(a4[,1])))
#194
b1=read.table("F:/7.txt",sep="\t")#206
length(intersect(b,as.character(b1[,1])))


###从TCGA下载27个cancer的4种pipeline的突变信息,取4种pipeline的gene并集去冗余为all_gene
#install.packages("devtools")
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")
cancer="UCS"
setwd("F:/tp53/CBioportal_TP53/TCGAmutationdata")
f1 <- GDCquery_Maf(cancer, pipelines = "muse")
f2 <- GDCquery_Maf(cancer, pipelines = "mutect")
f3 <- GDCquery_Maf(cancer, pipelines = "somaticsniper")
f4 <- GDCquery_Maf(cancer, pipelines = "varscan")

path="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
dir.create(paste(path,cancer,sep="")) 
path1=paste(path,cancer,"/",sep="")
f<-read.csv("F:/tp53/CBioportal_TP53/Result2_3956_gene2.txt",header=F,sep="\t")
sample_id=f[grep(paste("^",cancer,sep=""),f[,2]),1]
for(i in sample_id){
  dir.create(paste(path1,i,sep=""))
  a1=f1[grep(paste("^",i,sep=""),f1$Tumor_Sample_Barcode),]
  a2=f2[grep(paste("^",i,sep=""),f2$Tumor_Sample_Barcode),]
  a3=f3[grep(paste("^",i,sep=""),f3$Tumor_Sample_Barcode),]
  a4=f4[grep(paste("^",i,sep=""),f4$Tumor_Sample_Barcode),]
  write.table(a1,file =paste(path1,i,"/","muse.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(a2,file =paste(path1,i,"/","mutect.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(a3,file =paste(path1,i,"/","somaticsniper.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(a4,file =paste(path1,i,"/","varscan.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  b=unique(union(union(union(a1$Hugo_Symbol,a2$Hugo_Symbol),a3$Hugo_Symbol),a4$Hugo_Symbol))
  write.table(b,file =paste(path1,i,"/","all_gene.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}

i="TCGA-C5-A1M6"
dir.create(paste(path1,i,sep=""))
a1=f1[grep(paste("^",i,sep=""),f1$Tumor_Sample_Barcode),]
a2=f2[grep(paste("^",i,sep=""),f2$Tumor_Sample_Barcode),]
a3=f3[grep(paste("^",i,sep=""),f3$Tumor_Sample_Barcode),]
a4=f4[grep(paste("^",i,sep=""),f4$Tumor_Sample_Barcode),]
write.table(a1,file =paste(path1,i,"/","muse.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
write.table(a2,file =paste(path1,i,"/","mutect.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
write.table(a3,file =paste(path1,i,"/","somaticsniper.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
write.table(a4,file =paste(path1,i,"/","varscan.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
b=unique(union(union(union(a1$Hugo_Symbol,a2$Hugo_Symbol),a3$Hugo_Symbol),a4$Hugo_Symbol))
write.table(b,file =paste(path1,i,"/","all_gene.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
filename<-list.files(PATH)
cancer=na.omit(gsub("GDCdata",NA,filename))
for( i in cancer){
  path1=paste(PATH,i,"/",sep="")
  sample_id<-list.files(path1)
  for(j in sample_id){
    path=paste(path1,j,"/",sep="")
    fname<-list.files(path)
    f1<-read.csv(paste(path,fname[2],sep=""),header=T,sep="\t") 
    f2<-read.csv(paste(path,fname[3],sep=""),header=T,sep="\t")
    f3<-read.csv(paste(path,fname[4],sep=""),header=T,sep="\t") 
    f4<-read.csv(paste(path,fname[5],sep=""),header=T,sep="\t") 
    data1<-f1[,-grep("t_depth",colnames(f1)):-grep("n_depth",colnames(f1))]
    data1<-data1[,-grep("src_vcf_id",colnames(data1))]  #删去列t_depth\t_ref_count\t_alt_count\n_depth\src_vcf_id
    data2<-f2[,-grep("t_depth",colnames(f2)):-grep("n_depth",colnames(f2))]
    data2<-data2[,-grep("src_vcf_id",colnames(data2))]  
    data3<-f3[,-grep("t_depth",colnames(f3)):-grep("n_depth",colnames(f3))]
    data3<-data3[,-grep("src_vcf_id",colnames(data3))]  
    data4<-f4[,-grep("t_depth",colnames(f4)):-grep("n_depth",colnames(f4))]
    data4<-data4[,-grep("src_vcf_id",colnames(data4))]  
    a=rbind(data1,data2,data3,data4)
    write.table(unique(a),file =paste(path,"all_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  }
}


path1=paste(PATH,"PAAD","/",sep="")
path=paste(path1,"TCGA-IB-7647","/",sep="")
fname<-list.files(path)
f1<-read.csv(paste(path,fname[2],sep=""),header=T,sep="\t") 
f2<-read.csv(paste(path,fname[3],sep=""),header=T,sep="\t")
f3<-read.csv(paste(path,fname[4],sep=""),header=T,sep="\t") 
f4<-read.csv(paste(path,fname[5],sep=""),header=T,sep="\t") 
data1<-f1[,-grep("t_depth",colnames(f1)):-grep("n_depth",colnames(f1))]
data1<-data1[,-grep("src_vcf_id",colnames(data1))]  #删去列t_depth\t_ref_count\t_alt_count\n_depth\src_vcf_id
data2<-f2[,-grep("t_depth",colnames(f2)):-grep("n_depth",colnames(f2))]
data2<-data2[,-grep("src_vcf_id",colnames(data2))]  
data3<-f3[,-grep("t_depth",colnames(f3)):-grep("n_depth",colnames(f3))]
data3<-data3[,-grep("src_vcf_id",colnames(data3))]  
data4<-f4[,-grep("t_depth",colnames(f4)):-grep("n_depth",colnames(f4))]
data4<-data4[,-grep("src_vcf_id",colnames(data4))]  
a=rbind(data1,data2,data3,data4)
#dim(unique(a))
write.table(unique(a),file =paste(path,"all_mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"ALL_samplesWithMutationsData.csv",sep=""),header=T,sep=",") 
for( i in 1:nrow(f)){
  path1=paste(PATH,f[i,5],"/",sep="")
  if(as.character(f[i,1]) %in% list.files(path1)){
    f1<-read.csv(paste(path1,as.character(f[i,1]),"/","all_mutation.txt",sep=""),header=T,sep="\t") 
    f1<-cbind(f1[,1:10],f1[,35:38])
    f1<-unique(f1)
    f[i,4]=nrow(f1)
    if(length(grep("^TP53$",as.character(f1[,1])))!=0){
      f[i,3]=length(grep("^TP53$",as.character(f1[,1])))
      f[i,2]=1
    }else{
      f[i,3]=0
      f[i,2]=0
    }
  }
}
write.table(f,file =paste(PATH,"ALL_samplesWithMutationsData.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"ALL_samplesWithMutationsData.txt",sep=""),header=T,sep="\t")
cancer="CESC"
f1<-subset(f,as.character(f[,5])==cancer)
f1<-subset(f1,as.character(f1[,4])!="NA")
#table(f1[,3])
b=paste(as.data.frame(table(f1[,3]))[,2],round(as.data.frame(table(f1[,3]))[,2]/sum(as.data.frame(table(f1[,3]))[,2]),2),sep="/")
pdf(paste(PATH,"折线图/",cancer,"_tp53.pdf",sep=""))
plot(as.matrix(table(f1[,3])),type="l",pch=20,col= "darkblue",cex = 1,las=1,main =paste(cancer,"(",dim(f1)[1],")",sep=""),xlab="comutation number in tp53",ylab="number in patient")
text(as.data.frame(table(f1[,3]))[,1],as.data.frame(table(f1[,3]))[,2],b)
dev.off()


library(TCGAbiolinks)
packageVersion("TCGAbiolinks")
d=c("BLCA",list.files("F:/tp53/CBioportal_TP53/TCGAmutationdata")[9:13],list.files("F:/tp53/CBioportal_TP53/TCGAmutationdata")[15:27],list.files("F:/tp53/CBioportal_TP53/TCGAmutationdata")[29:36])
for(cancer in d){
  #cancer="BLCA"
  setwd("F:/tp53/CBioportal_TP53/TCGAmutationdata")
  f1 <- GDCquery_Maf(cancer, pipelines = "muse")
  f2 <- GDCquery_Maf(cancer, pipelines = "mutect")
  f3 <- GDCquery_Maf(cancer, pipelines = "somaticsniper")
  f4 <- GDCquery_Maf(cancer, pipelines = "varscan")
  f<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/SampleWithNoMutation.txt",header=F,sep="\t")
  sample_id=f[grep(paste("^",cancer,sep=""),as.character(f[,5])),1]
  for(i in sample_id){
    a1=f1[grep(paste("^",i,sep=""),f1$Tumor_Sample_Barcode),]
    #a1=a1[,-grep("t_depth",colnames(a1)):-grep("n_depth",colnames(a1))]
    #a1=a1[,-grep("src_vcf_id",colnames(a1))]
    #a1=a1[,-grep("all_effects",colnames(a1))]
    a1=cbind(a1[,1:10],a1[,35:38])
    a2=f2[grep(paste("^",i,sep=""),f2$Tumor_Sample_Barcode),]
    #a2=a2[,-grep("t_depth",colnames(a2)):-grep("n_depth",colnames(a2))]
    #a2=a2[,-grep("src_vcf_id",colnames(a2))]
    #a2=a2[,-grep("all_effects",colnames(a2))]
    a2=cbind(a2[,1:10],a2[,35:38])
    a3=f3[grep(paste("^",i,sep=""),f3$Tumor_Sample_Barcode),]
    #a3=a3[,-grep("t_depth",colnames(a3)):-grep("n_depth",colnames(a3))]
    #a3=a3[,-grep("src_vcf_id",colnames(a3))]
    #a3=a3[,-grep("all_effects",colnames(a3))]
    a3=cbind(a3[,1:10],a3[,35:38])
    a4=f4[grep(paste("^",i,sep=""),f4$Tumor_Sample_Barcode),]
    #a4=a4[,-grep("t_depth",colnames(a4)):-grep("n_depth",colnames(a4))]
    #a4=a4[,-grep("src_vcf_id",colnames(a4))]
    #a4=a4[,-grep("all_effects",colnames(a4))]
    a4=cbind(a4[,1:10],a4[,35:38])
    b=rbind(a1,a2,a3,a4)
    f[grep(paste("^",i,sep=""),as.character(f[,1])),4]=nrow(unique(b))
    write.table(f[grep(paste("^",i,sep=""),as.character(f[,1])),],file ="F:/tp53/CBioportal_TP53/TCGAmutationdata/SampleWithNoMutation1.txt",row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    
  }
}


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
a=f[4:3894,]
b=rbind(f[1:3,],f[3895:nrow(f),])
dataset <- data.frame(value = c(a[,4],b[,4]), group = factor(rep(c("2Tp53 Mutation", "1No Tp53 Mutation"),times = c(nrow(a), nrow(b)))))
pdf(paste(PATH,"boxplot/","Tp53 Mutation or No.pdf",sep=""))
boxplot( value ~ group, dataset, border= c("green","red"),main="The Distribution of Mutation Sample",ylab="Number of Mutation Sample",las=1, font.lab=2)
dev.off()
e1=subset(a,a[,3]==1)
e2=subset(a,a[,3]>1)
dataset1 <- data.frame(value = c(e1[,4],e2[,4]), group = factor(rep(c("3Tp53 Single Mutation", "4Tp53 CoMutation"), times = c(nrow(e1), nrow(e2)))))
pdf(paste(PATH,"boxplot/","Tp53 CoMutation or No.pdf",sep=""))
boxplot( value ~ group, dataset1, border= c("green","red"),main="The Distribution of Mutation Sample",ylab="Number of Mutation Sample",las=1, font.lab=2)
dev.off()
dataset2<-rbind(dataset,dataset1) 
#pdf(paste(PATH,"boxplot/","Tp53integrate.pdf",sep=""))
pdf(paste(PATH,"boxplot/","Tp53integrate1.pdf",sep=""))
boxplot( value ~ group, dataset2,log = "y", pch=".",border= c("green","red","purple","orange"),main="The Distribution of Mutation Sample",ylab="Number of Mutation Sample",las=1)
text(1.5,10000,labels=t.test(matrix(a[,4]),matrix(b[,4]))$p.value,cex=1,pos=3)
text(3.5,10000,labels=t.test(matrix(e1[,4]),matrix(e2[,4]))$p.value,cex=1,pos=3)
dev.off()
#ks.test(matrix(a[,4]),matrix(b[,4]))$p.value #显著 不同一分布 p-value < 2.2e-16
#ks.test(matrix(e1[,4]),matrix(e2[,4]))$p.value  #显著 不同一分布p-value = 2.577e-08


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
a=f[4:3894,]
b=rbind(f[1:3,],f[3895:nrow(f),])
dataset <- data.frame(value = c(a[,4],b[,4]), group = factor(rep(c("TP53Mu", "NoMu"),times = c(nrow(a), nrow(b)))))
pdf(paste(PATH,"boxplot/","Tp53 Mutation or No.pdf",sep=""))
boxplot( value ~ group, dataset, border= c("green","red"),main="The Distribution of Mutation Sample",ylab="Number of Mutation Sample",las=1, font.lab=2)
dev.off()
e1=subset(a,a[,3]==1)
e2=subset(a,a[,3]>1)
dataset1 <- data.frame(value = c(e1[,4],e2[,4]), group = factor(rep(c("SingleMu", "CoMu"), times = c(nrow(e1), nrow(e2)))))
pdf(paste(PATH,"boxplot/","Tp53 CoMutation or No.pdf",sep=""))
boxplot( value ~ group, dataset1, border= c("green","red"),main="The Distribution of Mutation Sample",ylab="Number of Mutation Sample",las=1, font.lab=2)
dev.off()
dataset2<-rbind(dataset,dataset1) 
dataset2[,2]=factor(dataset2[,2],c("NoMu","TP53Mu","SingleMu", "CoMu"))
#pdf(paste(PATH,"boxplot/","Tp53integrate.pdf",sep=""))
#pdf(paste(PATH,"boxplot/","Tp53integrate1.pdf",sep=""))
#pdf(paste(PATH,"boxplot/","Tp53integrate3.pdf",sep=""))
pdf(paste(PATH,"boxplot/","Tp53integrate33.pdf",sep=""))
boxplot( value ~ group, dataset2,log = "y", pch=".",lwd=2,border= c("#12B826","#1DADB8","#5e227f","#d22780"),main="The number of Mutation Genes",ylab="The number of Mutation Genes",las=2)
text(1.5,10000,labels=signif(t.test(matrix(a[,4]),matrix(b[,4]))$p.value,3),cex=1,pos=3)
text(3.5,10000,labels=signif(t.test(matrix(e1[,4]),matrix(e2[,4]))$p.value,3),cex=1,pos=3)
dev.off()
#ks.test(matrix(a[,4]),matrix(b[,4]))$p.value #显著 不同一分布 p-value < 2.2e-16
#ks.test(matrix(e1[,4]),matrix(e2[,4]))$p.value  #显著 不同一分布p-value = 2.577e-08

##做没突变、单突变、多突变生存差异 3class
library("survival")
f1<-read.csv(file ="F:/tp53/TCGA10188/Result2.txt",header=T,sep="\t") 
a1=subset(f1,f1$vital_status.y=="alive")
b1=subset(f1,f1$vital_status.y=="dead")
alive=cbind(as.character(a1$submitter_id),a1$vital_status.y,a1$days_to_last_follow_up)
dead=cbind(as.character(b1$submitter_id),b1$vital_status.y,b1$days_to_death.y)	
ff=data.frame(rbind(alive,dead))
f<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/3894_5886=9780.txt",header=T,sep="\t")
a=f[4:3894,]
b=rbind(f[1:3,],f[3895:nrow(f),])
dataset <- data.frame(Sample=c(as.character(b[,1])), group = factor(rep(c( "NoMu"),times = c(nrow(b)))))
data=merge(dataset,ff,by.x="Sample",by.y="X1",all.x=T)
e1=subset(a,a[,3]==1)
e2=subset(a,a[,3]>1)
dataset1 <- data.frame(Sample=c(as.character(e1[,1]),as.character(e2[,1])), group = factor(rep(c("SingleMu", "CoMu"), times = c(nrow(e1), nrow(e2)))))
data1=merge(dataset1,ff,by.x="Sample",by.y="X1",all.x=T)
data2<-rbind(data,data1) 
data2[,2]=factor(data2[,2],c("NoMu","SingleMu","CoMu"))
data2[,4]=as.numeric(as.character(data2[,4]))
data2[,3]=as.numeric(as.character(data2[,3]))
dif <- survdiff(Surv(data2[,4],data2[,3])~data2[,2])#求生存时间
kmsurvival<-survfit(Surv(data2[,4],data2[,3])~data2[,2],conf.type = "log-log")
pdf("F:/tp53/TCGA10188/survivalPlot_TP53_NoMu_Si_Co.pdf")
plot(kmsurvival, col=c("#12B826","#5e227f","#d22780"),lty = c('solid','solid','solid'),lwd=2, 
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft',cex=0.8,c("NoMu","SingleMu","CoMu"),lwd=2,
       col=c("#12B826","#5e227f","#d22780"),lty = c('solid','solid','solid'))
text(2200,1,cex=1,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
dev.off()


#data1[,4]=as.numeric(as.character(data1[,4]))
#data1[,3]=as.numeric(as.character(data1[,3]))
#dif1 <- survdiff(Surv(data1[,4],data1[,3])~data1[,2])
#kmsurvival1<-survfit(Surv(data1[,4],data1[,3])~data1[,2],conf.type = "log-log")



PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
e={}
e=cbind(e,as.character(unique(f[,5])),0)
for(i in 1:nrow(e)){
  b=subset(f,as.character(f[,5])==as.character(e[i,1]))
  e[i,2]=as.numeric(sum(b[,4]))
}
row.names(e)=e[,1]
pdf(paste(PATH,"boxplot/","mutation number of cancer type1.pdf",sep=""),width = 9,height=9)
barplot(t(e[,2]),main="mutation number of cancer type",cex.names=0.5,col="red",width=0.1,las=2,
        xlab="cancer type",ylab="mutation number")
for (i in 1:nrow(e)){
  text(i/36+i*0.1-0.15,as.numeric(e[i,2])-1800,e[i,2],cex=0.7)
}
dev.off()


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
b1=factor(f[,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b2=factor(f[1:3894,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b3=factor(f[3895:nrow(f),5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
dataset<- data.frame(value = f[,4], group=b1)
dataset1<- data.frame(value = f[1:3894,4], group=b2 )
dataset2<- data.frame(value = f[3895:nrow(f),4], group =b3)
pdf(paste(PATH,"boxplot/","Mutation Number of Different Cancers3.pdf",sep=""),width = 9,height=9)
x1<-boxplot(value ~ group,dataset,log = "y",notch=T,border= "red",pch=".",boxwex = 0.2, at = 1:27-0.5,main="Mutation Number of Different Cancers",xlab="Cancer Type",ylab="Mutation Number",las=2,font.lab=2)
boxplot(value ~ group,dataset2,log = "y",notch=T,border= "green",pch=".",boxwex = 0.2, at = 1:27-0.2,add= TRUE,las=2,axes=FALSE)
boxplot(value ~ group,dataset1,log = "y",notch=T,border= "orange",pch=".",boxwex = 0.2, at = 1:27+ 0.1,add= TRUE,las=2,axes=FALSE)
legend('topleft',cex=0.6,text.width=0.4,c("All","NO TP53 Mutation","TP53 Mutation"),fill = c("red","green","orange"))
mn.t <- tapply(dataset$value,dataset$group,median)
abline(h=sum(mn.t)/dim(mn.t),col = "black") #128.0556
dev.off()
#[order(dataset$group),]


x<-matrix(c(492,2,36,1),nrow=2,byrow = T)
d=chisq.test(x)
d1=fisher.test(x)
paste(d$p.value,"_",d1$p.value,sep="")


library(ggplot2)
library(ggpubr)
PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
b1=factor(f[,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b2=factor(f[1:3894,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b3=factor(f[3895:nrow(f),5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
dataset0<- data.frame(value = f[,4], group=b1)
dataset1<- data.frame(value = f[1:3894,4], group=b2 )
dataset2<- data.frame(value = f[3895:nrow(f),4], group =b3)
a0=cbind(class="All",dataset0)
a1=cbind(class="TP53 Mutation",dataset1)
a2=cbind(class="No TP53 Mutation",dataset2)
data=rbind(a0,a2,a1)
mn.t <- tapply(log(dataset0$value),dataset0$group,median)
mycolour<-scale_fill_manual(values=c("pink","lightgreen","yellow"))
p<-ggplot(data=data, aes(x=group,y=log(value)))+geom_boxplot(aes(fill=class),outlier.colour = NA)+mycolour
p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Cancer Type") + ylab("Mutation Number") + ggtitle("Mutation Number of Different Cancers") #添加横纵坐标，添加title
p <- p + guides(fill=guide_legend(title="Type"))+ theme(legend.position="top")+geom_hline(yintercept = sum(mn.t)/dim(mn.t), col="black")
p <- p + stat_compare_means(aes(group=data$class),label = "p.signif", label.y = 10)
p <- p+ stat_compare_means(aes(group=data$class),label = "p.format", label.y = rep(c(9.5,9.3,9),9),size=2)
ggsave(file=paste(PATH,"boxplot/","Mutation Number of Different Cancers4.pdf",sep=""),width = 10,height=9)



PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894-3=3891.txt",sep=""),header=T,sep="\t")
for( i in 1:nrow(f)){
  path1=paste(PATH,f[i,5],"/",sep="")
  if(as.character(f[i,1]) %in% list.files(path1)){
    f1<-read.csv(paste(path1,as.character(f[i,1]),"/","all_mutation.txt",sep=""),header=T,sep="\t") 
    f1<-cbind(f1[,1:10],f1[,35:38])
    f1<-unique(f1)
    if((length(grep("^TP53$",as.character(f1[,1])))!=0)){
      id=grep("^TP53$",as.character(f1[,1]))
      for(k in id){
        j=append(f[i,],f1[k,])
        write.table(j,file =paste(PATH,"3891.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
    }
  }
}


PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3891.txt",sep=""),header=F,sep="\t")
a<-as.matrix(table(f[,28]))  #all_type
pdf(paste(PATH,"柱状图/","3891_Mutation Type.pdf",sep=""))
x<-barplot(t(a[,1]/sum(a[,1])),main="The Distribution of 3891 Mutation Sample",cex.names=0.6,col="red",beside=TRUE, space=c(0,1),width=0.1,las=2,ylim=c(0,0.65),
        xlab="TP53 Mutation Type",ylab="Frequency")
lbls<-paste(" ",a)  ##factor de zuobiao
text(x,t(a[,1]/sum(a[,1])),labels=lbls,cex=0.6,pos=2)
text(x,t(a[,1]/sum(a[,1])),labels=round(t(a[,1]/sum(a[,1])),3),cex=0.5,pos=3)
dev.off()


b=as.matrix(table(f[,28],f[,3]))  #all_type_mutatin number
for(i in 1:nrow(b)){
  b[i,2]=sum(b[i,2:5])
}
pdf(paste(PATH,"柱状图/","3891_Mutation Type_Mutation number.pdf",sep=""))
x1<-barplot(t(b[,1:2]/sum(b[,1:2])),main="The Distribution of 3891 Mutation Sample",cex.names=0.6,col=c("red","green"),beside=TRUE, space=c(0,1),width=0.1,las=2,ylim=c(0,0.55),
        xlab="TP53 Mutation Type_Mutation number",ylab="Frequency")
text(x1,t(b[,1:2]/sum(b[,1:2])),labels=t(b[,1:2]),cex=0.6,pos=2)
text(x1,t(b[,1:2]/sum(b[,1:2])),labels=round(t(b[,1:2]/sum(b[,1:2])),3),cex=0.5,pos=3)
legend('topleft',cex=0.6,text.width=0.4,c("TP53 Mutation=1","TP53 Mutation>=2"),fill = c("red","green"))
dev.off()


f[,5]=factor(f[,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
d=as.matrix(table(f[,5],f[,28]))  #all_cancer_type
pdf(paste(PATH,"柱状图/","3891_Cancer_Mutation Type.pdf",sep=""))
#opar <- par(no.readonly=T) # 保存画图环境
#par(mfrow=c(3,9))
#par(oma=c(2,0,3,0))
#layout(matrix(1:27,3,9))
par(mar=c(1,1,1,1))
split.screen(c(9, 3)) 
for(j in 1:nrow(d)){
  #subplot(ceiling(j/9),j%%9)
  screen(j)
  x2<-barplot(t(d[j,1:ncol(d)]/sum(d[j,1:ncol(d)])),main=paste(as.character(rownames(d)[j]),sum(d[j,1:ncol(d)]),sep="/"),cex.main=0.2,cex.names=0.2,col=rainbow(13),beside=TRUE, space=c(0,0.5),width=0.01,las=2,axes=FALSE,ylim=c(0,1.2))
  text(x2,t(d[j,1:ncol(d)]/sum(d[j,1:ncol(d)])),labels=paste(t(d[j,1:ncol(d)]),round(t(d[j,1:ncol(d)]/sum(d[j,1:ncol(d)])),3),sep="/"),cex=0.2,pos=3)
  #text(x2,t(d[j,1:ncol(d)]/sum(d[,1:ncol(d)])),labels=round(t(d[j,1:ncol(d)]/sum(d[,1:ncol(d)])),3),cex=0.2,pos=3)
}
close.screen(all = TRUE) 
dev.off()

 
x2<-barplot(t(d[,1:ncol(d)]/sum(d[,1:ncol(d)])),main="The Distribution of 3891 Mutation Sample",cex.names=0.6,col=rainbow(13),beside=TRUE, space=c(0,1),width=0.1,las=2,ylim=c(0,0.1),
            xlab="TP53 Mutation Type_Cancer",ylab="Frequency")
text(x2,t(d[,1:ncol(d)]/sum(d[,1:ncol(d)])),labels=t(d[,1:ncol(d)]),cex=0.6,pos=2)
text(x2,t(d[,1:ncol(d)]/sum(d[,1:ncol(d)])),labels=round(t(d[,1:ncol(d)]/sum(d[,1:ncol(d)])),3),cex=0.5,pos=3)
legend('topleft',cex=0.6,text.width=0.4,c(unique(f[,28])),fill = rainbow(13))


library("survival")
PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
pdf(paste(PATH,"柱状图/","3894-3=3891_survival_the.number.of.TP53.mutation.pdf",sep=""))
f<-read.csv(paste(PATH,"3894-3=3891.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
#f<-read.csv(file ="F:/tp53/TCGA10188/Result1.txt",header=T,sep="\t",stringsAsFactors = F)
f$days_to_death.y<-as.numeric(gsub("NULL",NA,f$days_to_death.y))
f$vital_status.y<-gsub("alive",0,f$vital_status.y)
f$vital_status.y<-as.numeric(gsub("dead",1,f$vital_status.y))
dif <- survdiff(Surv(na.omit(f)$days_to_death.y,na.omit(f)$vital_status.y)~na.omit(f)$the.number.of.TP53.mutation)#求生存时间
#p_value <- pchisq(dif$chisq,1,lower.tail=F) #p= 0.3645801
kmsurvival<-survfit(Surv(na.omit(f)$days_to_death.y,na.omit(f)$vital_status.y)~na.omit(f)$the.number.of.TP53.mutation, conf.type = "log-log")
#summary(kmsurvival)
plot(kmsurvival, col=c('black','red',"blue"),lty="solid",
     xlab='survival time in days',ylab='survival probabilities')
legend('topright', c('the.number.of.TP53.mutation=1','the.number.of.TP53.mutation=2',"the.number.of.TP53.mutation=3"),lty="solid",col=c('black','red',"blue"))
text(2200,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
dev.off()




###从TCGA下载27个cancer的CNV信息\甲基化\Gene Expression数据
library(TCGAbiolinks)
cancer="BLCA"
setwd("F:/tp53/Methylation")
f1 <- GDCquery_Maf(cancer, pipelines = "muse")

query.maf.hg19 <- GDCquery(project = "TCGA-UCS", 
                           data.category = "Copy Number Variation", 
                           data.type = "Copy Number Segment",
                           platform =c("Affymetrix SNP Array 6.0"),
                           access = "open", 
                           legacy = FALSE)
GDCdownload(query.maf.hg19)
maf <- GDCprepare(query.maf.hg19)


query <- GDCquery(project ="TCGA-BLCA",
                  data.category = "DNA Methylation",
                  legacy = FALSE,
                  platform = c("Illumina Human Methylation 450","Illumina Human Methylation 27")
)
GDCdownload(query)
 

###TCGA hub的CNV数据处理
PATH="F:/tp53/CNV/"
f<-read.csv(paste(PATH,"ALL_samplesWithCNVData.txt",sep=""),header=T,sep="\t")
f[,2]=NA
for(i in 1:nrow(f)){
  cancer<-list.files("F:/tp53/CNV/TCGA hub/",pattern =paste(as.character(f[i,3]),'.*._genes',sep="")) 
  f1<-read.csv(paste("F:/tp53/CNV/TCGA hub/",cancer,"/","Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",sep=""),header=T,sep="\t")
  if(as.character(f[i,1]) %in% gsub("\\.","-",substring(colnames(f1),1,12))){
    if(length(grep("^TP53$",as.character(f1[,1])))!=0){
      f[i,2]=f1[grep("^TP53$",as.character(f1[,1])),grep(as.character(f[i,1]),gsub("\\.","-",substring(colnames(f1),1,12)))]
    }else{
      f[i,2]=5
    }
  }
}
write.table(f,file =paste(PATH,"ALL_samplesWithCNVData1.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


f1<-read.csv("F:/tp53/CBioportal_TP53/TCGAmutationdata/ALL_samplesWithMutationsData.txt",header=T,sep="\t")
a1=subset(f1,f1[,2]==1)
f2<-read.csv("F:/tp53/CNV/ALL_samplesWithCNVData1.txt",header=T,sep="\t")
a2=subset(f2,f2[,2]==1|f2[,2]==-1|f2[,2]==2|f2[,2]==-2)
d=merge(a2,a1,all=T) 
#d=union(as.character(a[,1]),as.character(f1[,1]))
write.table(d,"F:/1.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


###从TCGA下载27个cancer的Gene Expression count数据
#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
#"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
library(TCGAbiolinks)
setwd("F:/tp53/Gene Expression")
cancer <- "UCS"
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "Transcriptome Profiling",  #-········· 数据种类
                    data.type = "Gene Expression Quantification",  #······· 数据类型
                    workflow.type = "HTSeq - Counts")  #··················· 工作流类型HTSeq - FPKM
  GDCdownload(query)  #·························································· GDC下载数据
}

#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
#-- 将基因count整合成data.frame数据
library(R.utils)  #--解压.gz文件
library(TCGAbiolinks)
cancer <- "UCS"
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "Transcriptome Profiling",  #·········· 数据种类
                    data.type = "Gene Expression Quantification",  #······· 数据类型
                    workflow.type = "HTSeq - Counts")  #··················· 工作流类型
  setwd(paste("F:/tp53/Gene Expression/GDCdata/TCGA-",k,"/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification",sep = ""))
  gunzip(paste(query$results[[1]]$file_id[1],"/",query$results[[1]]$file_name[1],sep = ""))
  genecount <- read.table(paste(query$results[[1]]$file_id[1],"/",gsub(".gz","",query$results[[1]]$file_name[1]),sep = ""),sep = "\t",header = F)
  colnames(genecount) <- c("gene_ID",query$results[[1]]$cases[1])
  for (i in 2:length(query$results[[1]]$file_id)){
    gunzip(paste(query$results[[1]]$file_id[i],"/",query$results[[1]]$file_name[i],sep = ""))
    genecount_1 <- read.table(paste(query$results[[1]]$file_id[i],"/",gsub(".gz","",query$results[[1]]$file_name[i]),sep = ""),sep = "\t",header = F)
    colnames(genecount_1) <- c("gene_ID",query$results[[1]]$cases[i])
    genecount <- merge(genecount,genecount_1,all = T)
  }
  write.table(genecount,file = paste("F:/tp53/Gene Expression/Gene_count/",k,"_Gene_count.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
}


#处理基因count,文件为F:/tp53/Gene Expression/ENSG_Name2.txt
f<-read.table("F:/tp53/Gene Expression/gencode.v22.annotation.gtf/gencode.v22.annotation1.gtf",header=F,sep="\t")
write.table(unique(f[,9]),"F:/tp53/Gene Expression/ENSG_Name.txt",row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)

f1<-read.table("F:/tp53/Gene Expression/ENSG_Name.txt",header=F,sep="\t")
i=1
m=grep("gene_id",strsplit(as.character(f1[i,1]),";")[[1]])
n=grep("gene_name",strsplit(as.character(f1[i,1]),";")[[1]])
cat(gsub("gene_id ","",strsplit(as.character(f1[i,1]),";")[[1]][m]),"\t",file ="F:/tp53/Gene Expression/ENSG_Name1.txt",append=TRUE)
cat(gsub(" gene_name ","",strsplit(as.character(f1[i,1]),";")[[1]][n]),"\n",file ="F:/tp53/Gene Expression/ENSG_Name1.txt",append=TRUE)
f2<-read.table("F:/tp53/Gene Expression/ENSG_Name1.txt",header=F,sep="\t")   
for(i in 2:nrow(f1)){
  m=grep("gene_id",strsplit(as.character(f1[i,1]),";")[[1]])
  n=grep("gene_name",strsplit(as.character(f1[i,1]),";")[[1]])
  if(! gsub("gene_id ","",strsplit(as.character(f1[i,1]),";")[[1]][m]) %in% f2[,1]){
    cat(gsub("gene_id ","",strsplit(as.character(f1[i,1]),";")[[1]][m]),"\t",file ="F:/tp53/Gene Expression/ENSG_Name1.txt",append=TRUE)
    cat(gsub(" gene_name ","",strsplit(as.character(f1[i,1]),";")[[1]][n]),"\n",file ="F:/tp53/Gene Expression/ENSG_Name1.txt",append=TRUE)
  }
}

f2<-read.table("F:/tp53/Gene Expression/ENSG_Name1.txt",header=F,sep="\t")   
write.table(unique(f2),"F:/tp53/Gene Expression/ENSG_Name2.txt",row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)

#处理基因count,根据文件ENSG_Name2.txt ##for没有merge快
path="F:/tp53/Gene Expression/"
f<-read.table(paste(path,"ENSG_Name2.txt",sep=""),header=F,sep="\t")
filename<-list.files(paste(path,"Gene_count/",sep=""))
for(i in filename){
  f1<-read.table(paste(path,"Gene_count/",i,sep=""),header=T,sep="\t")
  colnames(f1)<-substring(colnames(f1),1,15)
  cat(colnames(f1),"\t",file=paste(path,"Gene_count1/",i,sep=""),append=TRUE)
  cat(file=paste(path,"Gene_count1/",i,sep=""),sep = "\n",append=TRUE)
  for(j in 1:nrow(f1)){
    if(as.character(f1[j,1]) %in% gsub(" ","",as.character(f[,1]))){
      k<-as.character(f[grep(as.character(f1[j,1]),gsub(" ","",as.character(f[,1]))),2])
      write.table(append(k,f1[j,2:ncol(f1)]),file=paste(path,"Gene_count1/",i,sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}

path="F:/tp53/Gene Expression/"
f<-read.table(paste(path,"ENSG_Name2.txt",sep=""),header=F,sep="\t")
f[,1]<-gsub(" ","",as.character(f[,1]))
f[,2]<-gsub(" ","",as.character(f[,2]))
filename<-list.files(paste(path,"Gene_FPKM/",sep=""))
for(i in filename[10:27]){
  f1<-read.table(paste(path,"Gene_FPKM/",i,"/Merge_matrix.txt",sep=""),header=T,sep="\t")
  colnames(f1)<-substring(colnames(f1),1,15)
  a=merge(f,f1[1:nrow(f1),],by.x="V1",by.y="Keys")
  colnames(a)[2]="gene_Name"
  write.table(a[,-1],file=paste(path,"Gene_FPKM1/",i,"_gene_fpkm.txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}


##找出p53上下文基因的fpkm值
path="D:/tp53/Gene Expression/"
f<-read.table(paste(path,"p53_context2.txt",sep=""),header=T,sep="\t")
filename<-list.files(paste(path,"Gene_FPKM1/",sep=""))
for(i in filename){
  f1<-read.table(paste(path,"Gene_FPKM1/",i,sep=""),header=T,sep="\t")
  a=merge(f,f1)
  write.table(a,file=paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm/",gsub("_gene","_p53_context",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}


##对每个样本（突变、WT没突变）的p53上下文基因的fpkm计算foldchange值=p53/normal_mean、WT/normal_mean
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f<-read.table("D:/tp53/Gene Expression/3894_5886=9780.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm/",sep=""))
for(i in filename){
  f1<-read.table(paste(path,"p53_context_fpkm/",i,sep=""),header=T,sep="\t")
  a1<-f1[,c(1,grep(".*.11$",colnames(f1)))]
  mutation<-subset(f,as.character(f[,2])==1 & as.character(f[,5])==gsub("_p53_context_fpkm.txt","",i))[,1]
  nomutation<-subset(f,is.na(as.character(f[,2]))& as.character(f[,5])==gsub("_p53_context_fpkm.txt","",i))[,1]
  colnames(f1)<-gsub("\\.","-",substring(colnames(f1),1,12))
  a2<-f1[,c("gene_Name",intersect(as.character(mutation),colnames(f1)))]
  a3<-f1[,c("gene_Name",intersect(as.character(nomutation),colnames(f1)))]
  write.table(a1,file=paste(path,"p53_context_normal_fpkm/",gsub("_fpkm","_normal_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a2,file=paste(path,"p53_context_mutation_fpkm/",gsub("_fpkm","_mutation_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a3,file=paste(path,"p53_context_nomutation_fpkm/",gsub("_fpkm","_nomutation_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}


#foldchange值=p53/normal_mean、WT/normal_mean取log2
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
for(cancer in  c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                 "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")){
  f1<-read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3<-read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(colnames(f2),"\t",file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
  cat(paste(colnames(f3)[-1],"_no",sep=""),"\n",file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    n=log2(f3[i,2:ncol(f3)]/k[i])
    write.table(append(m,n),file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}

path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
cancer<-"SKCM" ##normal=1
f1<-read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
f2<-read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
f3<-read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
cat(colnames(f2),"\t",file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
cat(paste(colnames(f3)[-1],"_no",sep=""),"\n",file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
for(i in 1:nrow(f2)){
  m=cbind(f2[i,1],log2(f2[i,2:ncol(f2)]/f1[i,2]))
  n=log2(f3[i,2:ncol(f3)]/f1[i,2])
  write.table(append(m,n),file=paste(path,"p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
}


#foldchange值=p53/WT_mean取log2
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC","LAML","LGG","OV","TGCT","UCS","SKCM")){
  f1<-read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(paste(colnames(f2),"_p53/WT_mean",sep=""),"\n",file=paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    write.table(m,file=paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}


##对每个cancer的p53上下文基因的foldchange值=p53/normal_mean、WT/normal_mean做热图heatmap.2和pheatmap
#install.packages("gplots")
library(gplots)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t")
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1<-na.omit(f1)[-1,]
  row.names<-f1$gene_Name
  x<-f1[,2:ncol(f1)]
  y<-as.matrix(x)
  pdf(file=paste(path,"p53_context_fpkm_foldchange_heatmap/",gsub(".txt",".pdf",i),sep=""),width=8, height=8)
  heatmap.2(y,Rowv=FALSE,Colv=FALSE,col=colorsChoice(3),breaks=c(min(y),-1,1,max(y)),key=T,keysize=1.5,margins=c(5,10),labRow=row.names,cexRow=0.3,cexCol=0.3,trace="none")
  dev.off()
}

#install.packages("pheatmap")
library(pheatmap)
library(gplots)
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("F:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1<-na.omit(f1)[-1,]
  rownames(f1)<-f1$gene_Name
  f2=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),]
  x<-f2[,2:ncol(f2)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation", "NoMutation"), c(ncol(f1)-length(grep("_no",colnames(f1)))-1, length(grep("_no",colnames(f1)))))),row.names=colnames(f1)[2:ncol(f1)])  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/1heatmap/",gsub(".txt",".pdf",i),sep=""),width=8,height=8)  
}


##对每个cancer的p53上下文基因的cor\foldchange值=p53/normal_mean、WT/normal_mean、p53/WT_mean做热图pheatmap
library(pheatmap)
library(gplots)
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("F:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  #f1<-na.omit(f1)[-1,]
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),]
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  #f2<-na.omit(f2)[-1,]
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]

  e1<-read.table(paste(path,"p53_context_normal_fpkm/",gsub("_fc.txt","",i),"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  e2<-read.table(paste(path,"p53_context_mutation_fpkm/",gsub("_fc.txt","",i),"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  e3<-read.table(paste(path,"p53_context_nomutation_fpkm/",gsub("_fc.txt","",i),"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  e4<-cbind(e2,e3)
  a1=e1[grep("^TP53$",e1[,1]),2:ncol(e1)]
  a2=e2[grep("^TP53$",e2[,1]),2:ncol(e2)]
  a3=e3[grep("^TP53$",e3[,1]),2:ncol(e3)]
  a4=e4[grep("^TP53$",e4[,1]),2:ncol(e4)]
  f3=data.frame()
  for(j in 1:nrow(e1)){
      f3[j,1]=2*cor(as.numeric(a1),as.numeric(e1[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e1[,1]),2:ncol(e1)]))
      f3[j,2]=2*cor(as.numeric(a2),as.numeric(e2[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e2[,1]),2:ncol(e2)]))
      f3[j,3]=2*cor(as.numeric(a3),as.numeric(e3[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e3[,1]),2:ncol(e3)]))
      f3[j,4]=2*cor(as.numeric(a4),as.numeric(e4[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e4[,1]),2:ncol(e4)]))
  }
  rownames(f3)=e1[,1]
  f3[sapply(f3,is.na)]<-0
  p=data.frame(rep(c(data.frame(f3[,1]),data.frame(f3[,2]),data.frame(f3[,3]),data.frame(f3[,4])),each=50))
  colnames(p)=rep(c("cor_normal","cor_Mutation", "cor_NoMutation","cor_cancer"),each=50)
  p[,201]=e1[,1]
  colnames(p)[201]="gene_Name"
  rownames(p)=p[,201]
  p<-na.omit(p)
  f33=p[intersect(as.character(f[,1]),as.character(p[,201])),]

  d=merge(f11,f22,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
  rownames(d)<-d[,1]
  dd=merge(f33,d,by="gene_Name")
  rownames(dd)<-dd[,1]
  d1=dd[intersect(as.character(f[,1]),as.character(dd[,1])),]
  d1[sapply(d1,is.na)]<-0
  x<-d1[,2:ncol(d1)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep(c("cor_normal","cor_Mutation", "cor_NoMutation","cor_cancer","Mutation", "NoMutation","Mutation/NoMutation"), c(50,50,50,50,ncol(f1)-length(grep("_no",colnames(f1)))-1, length(grep("_no",colnames(f1))),ncol(f2)-1))),row.names=colnames(dd)[2:ncol(dd)])  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/2heatmap/",gsub("_fc.txt","_cor0.5_p53_WT.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
}

 
##对每个cancer22\Allcancer，根据TP53基因的foldchange值=p53/normal_mean、WT/normal_mean
##把样本分为p53WT、p53突变上调、p53突变下调、p53突变不调四种类型，分别做生存分析
##去GBM6\KIRP9\PCPG14\PRAD15\SKCM18\THCA20=22-6=16
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
f$days_to_death.y<-as.numeric(gsub("NULL",NA,f$days_to_death.y))
f$vital_status.y<-gsub("alive",0,f$vital_status.y)
f$vital_status.y<-as.numeric(gsub("dead",1,f$vital_status.y))
rownames(f)=f[,1]
f=na.omit(f)
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
h1={}
h2={}
h3={}
for(i in filename[c(1:5,7:8,10:13,16:17,19,21:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1<-f1[-1,]
  f1[,2]<-as.numeric(f1[,2])
  a=subset(f1,f1[,1]=="TP53")
  WT=gsub("_no","",gsub("\\.","-",colnames(a)[grep("_no",as.character(colnames(a)))]))
  a1=data.frame(a[,-grep("_no",as.character(colnames(a)))])
  p53_N=gsub("\\.","-",colnames(a1)[which(a1[1,]> -1 & a1[1,]< 1)]) 
  p53_up=gsub("\\.","-",colnames(a1)[which(a1[1,]>= 1)][-1])
  p53_down=gsub("\\.","-",colnames(a1)[which(a1[1,]<= -1)])
  
  pdf(paste(path,"p53_context_fpkm_foldchange_survival/",gsub("_fc.txt","",i),"_survivalplot.pdf",sep=""))
  opar <- par(no.readonly=T) # 保存画图环境
  par(mfrow=c(1,3))
  par(oma=c(2,0,3,0))
  ff1=rbind(cbind(class=1,f[WT,]),cbind(class=5,f[p53_N,]),cbind(class=5,f[p53_up,]),cbind(class=5,f[p53_down,]))
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h1=rbind(h1,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  plot(kmsurvival1, lty = 'solid', col=c('black','orange'),
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53'), lty='solid',
         col=c('black','orange'))
  text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""),paste(dif1$n[1],dif1$n[2],sep="_"),sep="\n"))
  
  ff2=rbind(cbind(class=2,f[p53_N,]),cbind(class=3,f[p53_up,]),cbind(class=4,f[p53_down,]))
  b12=subset(ff2,ff2$vital_status.y==0)
  b22=subset(ff2,ff2$vital_status.y==1)
  alive2=cbind(b12$class,b12$vital_status.y,b12$days_to_last_follow_up)
  dead2=cbind(b22$class,b22$vital_status.y,b22$days_to_death.y)	
  e2=rbind(alive2,dead2)
  e2=na.omit(e2)
  h2=rbind(h2,e2)
  dif2 <- survdiff(Surv(as.numeric(e2[,3]),as.numeric(e2[,2]))~e2[,1])#求生存时间
  kmsurvival2<-survfit(Surv(as.numeric(e2[,3]),as.numeric(e2[,2]))~e2[,1],conf.type = "log-log")
  plot(kmsurvival2, lty = 'solid', col=c('red','green','purple'),
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c('p53_N','p53_up','p53_down'), lty='solid',
         col=c('red','green','purple'))
  text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),4),sep=""),paste(dif2$n[1],dif2$n[2],dif2$n[3],sep="_"),sep="\n"))
  
  ff3=rbind(cbind(class=1,f[WT,]),cbind(class=2,f[p53_N,]),cbind(class=3,f[p53_up,]),cbind(class=4,f[p53_down,]))
  b13=subset(ff3,ff3$vital_status.y==0)
  b23=subset(ff3,ff3$vital_status.y==1)
  alive3=cbind(b13$class,b13$vital_status.y,b13$days_to_last_follow_up)
  dead3=cbind(b23$class,b23$vital_status.y,b23$days_to_death.y)	
  e3=rbind(alive3,dead3)
  e3=na.omit(e3)
  h3=rbind(h3,e3)
  dif3 <- survdiff(Surv(as.numeric(e3[,3]),as.numeric(e3[,2]))~e3[,1])#求生存时间
  kmsurvival3<-survfit(Surv(as.numeric(e3[,3]),as.numeric(e3[,2]))~e3[,1],conf.type = "log-log")
  plot(kmsurvival3, lty = 'solid', col=c('black','red','green','purple'),
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53_N','p53_up','p53_down'), lty='solid',
         col=c('black','red','green','purple'))
  text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif3$chisq,1,lower.tail=F),4),sep=""),paste(dif3$n[1],dif3$n[2],dif3$n[3],dif3$n[4],sep="_"),sep="\n"))
  dev.off()
}
pdf(paste(path,"p53_context_fpkm_foldchange_survival/Allcancer_survivalplot.pdf",sep=""))
opar <- par(no.readonly=T) # 保存画图环境
par(mfrow=c(1,3))
par(oma=c(2,0,3,0))
dif4 <- survdiff(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,1])#求生存时间
kmsurvival4<-survfit(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,1],conf.type = "log-log")
plot(kmsurvival4, lty = 'solid', col=c('black','orange'),
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53'), lty='solid',
       col=c('black','orange'))
text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif4$chisq,1,lower.tail=F),4),sep=""),paste(dif4$n[1],dif4$n[2],sep="_"),sep="\n"))
dif5 <- survdiff(Surv(as.numeric(h2[,3]),as.numeric(h2[,2]))~h2[,1])#求生存时间
kmsurvival5<-survfit(Surv(as.numeric(h2[,3]),as.numeric(h2[,2]))~h2[,1],conf.type = "log-log")
plot(kmsurvival5, lty = 'solid', col=c('red','green','purple'),
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c('p53_N','p53_up','p53_down'), lty='solid',
       col=c('red','green','purple'))
text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif5$chisq,1,lower.tail=F),4),sep=""),paste(dif5$n[1],dif5$n[2],dif5$n[3],sep="_"),sep="\n"))
dif6 <- survdiff(Surv(as.numeric(h3[,3]),as.numeric(h3[,2]))~h3[,1])#求生存时间
kmsurvival6<-survfit(Surv(as.numeric(h3[,3]),as.numeric(h3[,2]))~h3[,1],conf.type = "log-log")
plot(kmsurvival6, lty = 'solid', col=c('black','red','green','purple'),
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53_N','p53_up','p53_down'), lty='solid',
       col=c('black','red','green','purple'))
text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif6$chisq,1,lower.tail=F),4),sep=""),paste(dif6$n[1],dif6$n[2],dif6$n[3],dif6$n[4],sep="_"),sep="\n"))
dev.off()




 
i=filename[4] 
f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
f1[,2]<-as.numeric(f1[,2])
f1<-f1[-1,]
a=subset(f1,f1[,1]=="TP53")
WT=gsub("_no","",gsub("\\.","-",colnames(a)[grep("_no",as.character(colnames(a)))]))
a1=data.frame(a[,-grep("_no",as.character(colnames(a)))])
p53_N=gsub("\\.","-",colnames(a1)[which(a1[1,]> -1 & a1[1,]< 1)]) 
p53_up=gsub("\\.","-",colnames(a1)[which(a1[1,]>= 1)][-1])
p53_down=gsub("\\.","-",colnames(a1)[which(a1[1,]<= -1)])
pdf(paste(path,"p53_context_fpkm_foldchange_survival/",gsub("_fc.txt","",i),"_survivalplot.pdf",sep=""))
opar <- par(no.readonly=T) # 保存画图环境
par(mfrow=c(1,3))
par(oma=c(2,0,3,0))
ff1=rbind(cbind(class=1,f[WT,]),cbind(class=5,f[p53_N,]),cbind(class=5,f[p53_up,]),cbind(class=5,f[p53_down,]))
b11=subset(ff1,ff1$vital_status.y==0)
b21=subset(ff1,ff1$vital_status.y==1)
alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
plot(kmsurvival1, lty = 'solid', col=c('black','orange'),
     xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53'), lty='solid',
       col=c('black','orange'))
text(2200,1,cex=0.8,paste(paste("p=",pchisq(dif1$chisq,1,lower.tail=F),sep=""),paste(dif1$n[1],dif1$n[2],sep="_"),sep="\n"))
ff2=rbind(cbind(class=2,f[p53_N,]),cbind(class=3,f[p53_up,]),cbind(class=4,f[p53_down,]))
b12=subset(ff2,ff2$vital_status.y==0)
b22=subset(ff2,ff2$vital_status.y==1)
alive2=cbind(b12$class,b12$vital_status.y,b12$days_to_last_follow_up)
dead2=cbind(b22$class,b22$vital_status.y,b22$days_to_death.y)	
e2=rbind(alive2,dead2)
dif2 <- survdiff(Surv(as.numeric(e2[,3]),as.numeric(e2[,2]))~e2[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(e2[,3]),as.numeric(e2[,2]))~e2[,1],conf.type = "log-log")
plot(kmsurvival2, lty = 'solid', col=c('red','green','purple'),
     xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
legend('bottomleft', cex=0.6,text.width=0.4,c('p53_N','p53_up','p53_down'), lty='solid',
       col=c('red','green','purple'))
text(2200,1,cex=0.8,paste(paste("p=",pchisq(dif2$chisq,1,lower.tail=F),sep=""),paste(dif2$n[1],dif2$n[2],dif2$n[3],sep="_"),sep="\n"))
ff3=rbind(cbind(class=1,f[WT,]),cbind(class=2,f[p53_N,]),cbind(class=3,f[p53_up,]),cbind(class=4,f[p53_down,]))
b13=subset(ff3,ff3$vital_status.y==0)
b23=subset(ff3,ff3$vital_status.y==1)
alive3=cbind(b13$class,b13$vital_status.y,b13$days_to_last_follow_up)
dead3=cbind(b23$class,b23$vital_status.y,b23$days_to_death.y)	
e3=rbind(alive3,dead3)
dif3 <- survdiff(Surv(as.numeric(e3[,3]),as.numeric(e3[,2]))~e3[,1])#求生存时间
kmsurvival3<-survfit(Surv(as.numeric(e3[,3]),as.numeric(e3[,2]))~e3[,1],conf.type = "log-log")
plot(kmsurvival3, lty = 'solid', col=c('black','red','green','purple'),
     xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53_N','p53_up','p53_down'), lty='solid',
       col=c('black','red','green','purple'))
text(2200,1,cex=0.8,paste(paste("p=",pchisq(dif3$chisq,1,lower.tail=F),sep=""),paste(dif3$n[1],dif3$n[2],dif3$n[3],dif3$n[4],sep="_"),sep="\n"))
dev.off()


##对每个cancer27\Allcancer，根据TP53基因的foldchange_p53%WTmean值=p53/WTmean 
##把样本分为p53突变上调、p53突变下调、p53突变不调3种类型，分别做生存分析
##去KIRC8\KIRP9\PCPG17\PRAD18\SKCM21\TGCT23\THCA24=27-7=20cancer 
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
f$days_to_death.y<-as.numeric(gsub("NULL",NA,f$days_to_death.y))
f$vital_status.y<-gsub("alive",0,f$vital_status.y)
f$vital_status.y<-as.numeric(gsub("dead",1,f$vital_status.y))
rownames(f)=f[,1]
f=na.omit(f)
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
h1={}
for(i in filename[c(1:7,10:16,19:20,22,25:27)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]] 
  f1<-f1[-1,]
  f1[,2]<-as.numeric(f1[,2])
  a=subset(f1,f1[,1]=="TP53")
  a1=data.frame(a[,2:ncol(a)])
  p53_N=gsub("\\.","-",gsub("_p53.WT_mean","",colnames(a1)[which(a1[1,]> -1 & a1[1,]< 1)]))
  p53_up=gsub("\\.","-",gsub("_p53.WT_mean","",colnames(a1)[which(a1[1,]>= 1)]))
  p53_down=gsub("\\.","-",gsub("_p53.WT_mean","",colnames(a1)[which(a1[1,]<= -1)]))
  
  pdf(paste(path,"p53_context_fpkm_foldchange_p53_WTmean_survival/",gsub("_fc.txt","",i),"_survivalplot.pdf",sep=""))
  ff1=rbind(cbind(class=1,f[p53_N,]),cbind(class=2,f[p53_up,]),cbind(class=3,f[p53_down,]))
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h1=rbind(h1,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  plot(kmsurvival1, lty = 'solid', col=c('black','red',"green"),
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_fc.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c('p53_N','p53_up',"p53_down"), lty='solid',
         col=c('black','red',"green"))
  text(1000,1,cex=0.8,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""),paste(dif1$n[1],dif1$n[2],dif1$n[3],sep="_"),sep="\n"))
  dev.off()
}
pdf(paste(path,"p53_context_fpkm_foldchange_p53_WTmean_survival/Allcancer_survivalplot.pdf",sep=""))
dif4 <- survdiff(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,1])#求生存时间
kmsurvival4<-survfit(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,1],conf.type = "log-log")
plot(kmsurvival4, lty = 'solid', col=c('black','red',"green"),
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c('p53_N','p53_up',"p53_down"), lty='solid',
       col=c('black','red',"green"))
text(1000,1,cex=0.8,paste(paste("p=",round(pchisq(dif4$chisq,1,lower.tail=F),4),sep=""),paste(dif4$n[1],dif4$n[2],dif4$n[3],sep="_"),sep="\n"))
dev.off()




#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
#"LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
#"LAML","LGG","OV","TGCT","UCS"no normal sample
#"SKCM",1 normal sample
##对每个cancer，求TP53和其他gene的cor相关系数（在normal、mutation、nomutation、cancer样本中）
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
cancer="SKCM"
f1<-read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
f2<-read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
f3<-read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
f4<-cbind(f2,f3)
a1=f1[grep("^TP53$",f1[,1]),2:ncol(f1)]
a2=f2[grep("^TP53$",f2[,1]),2:ncol(f2)]
a3=f3[grep("^TP53$",f3[,1]),2:ncol(f3)]
a4=f4[grep("^TP53$",f4[,1]),2:ncol(f4)]
for(i in 1:nrow(f1)){
  #as.character(f1[i,1])
  b1=cor(as.numeric(a1),as.numeric(f1[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f1[,1]),2:ncol(f1)]))
  b2=cor(as.numeric(a2),as.numeric(f2[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f2[,1]),2:ncol(f2)]))
  b3=cor(as.numeric(a3),as.numeric(f3[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f3[,1]),2:ncol(f3)]))
  b4=cor(as.numeric(a4),as.numeric(f4[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f4[,1]),2:ncol(f4)]))
  pdf(file=paste(path,"p53_context_fpkm_foldchange_heatmap/TP53_Gene_cor/",cancer,"/",paste("TP53_",as.character(f1[i,1]),sep=""),".pdf",sep=""),width=8, height=8)
  plot(1:4,c(b1,b2,b3,b4),pch=20,las=1,xlab="class",ylab="cor_value",col=c("black","red","blue","green"),main=paste("TP53_",as.character(f1[i,1]),sep=""),ylim=c(-1,1))
  text(1:4,c(b1,b2,b3,b4),cex=0.6,c(b1,b2,b3,b4))
  legend('bottomright', c("normal","mutation","nomutation","cancer"),pch=20,col=c("black","red","blue","green"))
  #plot(1:3,c(b2,b3,b4),pch=20,las=1,xlab="class",ylab="cor_value",col=c("red","blue","green"),main=paste("TP53_",as.character(f1[i,1]),sep=""),ylim=c(-1,1))
  #text(1:3,c(b2,b3,b4),cex=0.6,c(b2,b3,b4))
  #legend('bottomright', c("mutation","nomutation","cancer"),pch=20,col=c("red","blue","green"))
  dev.off()
}

for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")){
  f1<-read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3<-read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  f4<-cbind(f2,f3)
  a1=f1[grep("^TP53$",f1[,1]),2:ncol(f1)]
  a2=f2[grep("^TP53$",f2[,1]),2:ncol(f2)]
  a3=f3[grep("^TP53$",f3[,1]),2:ncol(f3)]
  a4=f4[grep("^TP53$",f4[,1]),2:ncol(f4)]
  for(i in 1:nrow(f1)){
    #as.character(f1[i,1])
    b1=cor(as.numeric(a1),as.numeric(f1[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f1[,1]),2:ncol(f1)]))
    b2=cor(as.numeric(a2),as.numeric(f2[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f2[,1]),2:ncol(f2)]))
    b3=cor(as.numeric(a3),as.numeric(f3[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f3[,1]),2:ncol(f3)]))
    b4=cor(as.numeric(a4),as.numeric(f4[grep(paste("^",as.character(f1[i,1]),"$",sep=""),f4[,1]),2:ncol(f4)]))
    pdf(file=paste(path,"p53_context_fpkm_foldchange_heatmap/TP53_Gene_cor/",cancer,"/",paste("TP53_",as.character(f1[i,1]),sep=""),".pdf",sep=""),width=8, height=8)
    plot(1:4,c(b1,b2,b3,b4),pch=20,las=1,xlab="class",ylab="cor_value",col=c("black","red","blue","green"),main=paste("TP53_",as.character(f1[i,1]),sep=""),ylim=c(-1,1))
    text(1:4,c(b1,b2,b3,b4),cex=0.6,c(b1,b2,b3,b4))
    legend('bottomright', c("normal","mutation","nomutation","cancer"),pch=20,col=c("black","red","blue","green"))
    #plot(1:3,c(b2,b3,b4),pch=20,las=1,xlab="class",ylab="cor_value",col=c("red","blue","green"),main=paste("TP53_",as.character(f1[i,1]),sep=""),ylim=c(-1,1))
    #text(1:3,c(b2,b3,b4),cex=0.6,c(b2,b3,b4))
    #legend('bottomright', c("mutation","nomutation","cancer"),pch=20,col=c("red","blue","green"))
    dev.off()
  }
}


##对每个cancer，用pheatmap画TP53和其他gene的cor相关系数的热图（在normal、mutation、nomutation、cancer样本中）
library(pheatmap)
library(gplots)
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("F:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1<-na.omit(f1)[-1,]
  rownames(f1)<-f1$gene_Name
  f2=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),]
  x<-f2[,2:ncol(f2)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  e1<-read.table(paste(path,"p53_context_normal_fpkm/",gsub("_fc.txt","",i),"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  e2<-read.table(paste(path,"p53_context_mutation_fpkm/",gsub("_fc.txt","",i),"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  e3<-read.table(paste(path,"p53_context_nomutation_fpkm/",gsub("_fc.txt","",i),"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  e4<-cbind(e2,e3)
  a1=e1[grep("^TP53$",e1[,1]),2:ncol(e1)]
  a2=e2[grep("^TP53$",e2[,1]),2:ncol(e2)]
  a3=e3[grep("^TP53$",e3[,1]),2:ncol(e3)]
  a4=e4[grep("^TP53$",e4[,1]),2:ncol(e4)]
  d=data.frame()
  for(j in 1:nrow(e1)){
    d[j,1]=cor(as.numeric(a1),as.numeric(e1[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e1[,1]),2:ncol(e1)]))
    d[j,2]=cor(as.numeric(a2),as.numeric(e2[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e2[,1]),2:ncol(e2)]))
    d[j,3]=cor(as.numeric(a3),as.numeric(e3[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e3[,1]),2:ncol(e3)]))
    d[j,4]=cor(as.numeric(a4),as.numeric(e4[grep(paste("^",as.character(e1[j,1]),"$",sep=""),e4[,1]),2:ncol(e4)]))
  }
  row.names(d)=e1[,1]
  x1=d[rownames(y),]
  y1<-as.matrix(x1)
  colnames(y1)=c("normal","Mutation", "NoMutation","cancer")
  pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(-1,-0.5,0.5,1),breaks=c(-1,-0.5,0.5,1),annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/heatmap_cor_0.5/",gsub("_fc.txt","_cor.pdf",i),sep=""),border_color=NA,width=8,height=8)  
}


##对每个cancer\all cancer，对p53上下文基因，计算在mutation、nomutation、cancer样本中上下调基因的个数
path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
p=data.frame()
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1<-na.omit(f1)[-1,]
  rownames(f1)<-f1$gene_Name
  f1=f1[,-1]
  e1=data.frame(f1[,-grep("_no",as.character(colnames(f1)))])
  e2=f1[,grep("_no",as.character(colnames(f1)))]
  n=data.frame() 
  for(j in 1: nrow(f1)){
    cancerup=length(which(f1[j,]>= 1))
    cancerdown=length(which(f1[j,]<= -1))
    mutationup=length(which(e1[j,]>= 1))
    mutationdown=length(which(e1[j,]<= -1))
    nomutationup=length(which(e2[j,]>= 1))
    nomutationdown=length(which(e2[j,]<= -1))
    m=cbind(rownames(f1)[j],cancerup,cancerdown,mutationup,mutationdown,nomutationup,nomutationdown)
    n=rbind(n,m)
  }
  write.table(n,file=paste(path,"p53_context_fpkm_foldchange_updown/",gsub("_fc","_fc_updown",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  p=rbind(p,n)
}
g=data.frame()
for(ii in unique(p[,1])){
  b=subset(p,p[,1]==ii)
  b <- apply(b[,2:ncol(b)],2,as.numeric)
  if(ncol(data.frame(b))==1){
    g=rbind(g,data.frame(t(data.frame(append(ii, b)))))
  }else{
    g=rbind(g,data.frame(t(data.frame(append(ii, apply(b,2,sum))))))
  }
}
write.table(g,file=paste(path,"p53_context_fpkm_foldchange_updown/","Allcancer_fc_updown.txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)


##对每个cancer，计算TP53在mutation、nomutation、mutation/nomutation样本中total、差异表达、上、下调样本的个数
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
m=paste("TP53Cancer",'Mutotal','Mudiff','Muup','Mudown','Nototal','Nodiff','Noup','Nodown','MNtotal','MNdiff','MNup','MNdown',sep="\t")
cat(m,file=paste(path,"Allcancer_fc_updown.txt",sep=""),sep="\n",append=TRUE)
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename[15:22]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  f11=f1[grep("^TP53$",f1[,1]),2:ncol(f1)]
  a=f11[,-grep("_no",as.character(colnames(f11)))]
  b=f11[,grep("_no",as.character(colnames(f11)))]
  a1=length(which(a[1,]>= 1))
  a2=length(which(a[1,]<= -1))
  b1=length(which(b[1,]>= 1))
  b2=length(which(b[1,]<= -1))
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  f22=f2[grep("^TP53$",f2[,1]),2:ncol(f2)]
  c1=length(which(f22[1,]>= 1))
  c2=length(which(f22[1,]<= -1))
  m1=paste(gsub("_fc.txt","",i),ncol(a),a1+a2,a1,a2,ncol(b),b1+b2,b1,b2,ncol(f22),c1+c2,c1,c2,sep="\t")
  write.table(m1,file=paste(path,"Allcancer_fc_updown.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
}

  
#对在mutation、mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename[15:22]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  f1=f1[,-grep("_no",as.character(colnames(f1)))]
  a=colnames(sort(f1[grep("^TP53$",f1[,1]),2:ncol(f1)]))
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a)]
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),c("gene_Name_p53/WT_mean",b)]
  
  d=merge(f11,f22,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
  rownames(d)<-d[,1]
  d1=d[intersect(as.character(f[,1]),as.character(d[,1])),]
  d1[sapply(d1,is.na)]<-0
  x<-d1[,2:ncol(d1)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation","Mutation/NoMutation"), c(ncol(f1)-1,ncol(f2)-1))),row.names=colnames(d)[2:ncol(d)])  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/3heatmap/",gsub("_fc.txt","_M_MN.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
}



#对在mutation、nomutation、mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
##除去PCPG14
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice1<- colorRampPalette(c("green","white","red"))
colorsChoice2<- colorRampPalette(c("green","white","pink","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))  
h1=data.frame(gene_Name=as.character(f[,1]))
h2=data.frame(gene_Name=as.character(f[,1]))
h3=data.frame(as.character(f[,1]))
colnames(h3)="gene_Name_p53/WT_mean"
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=colnames(sort(f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]))
  a2=colnames(sort(f1["TP53",grep("_no",as.character(colnames(f1)))]))
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a1)]
  f12=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a2)]
  h1=merge(h1,f11)
  h2=merge(h2,f12)
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),c("gene_Name_p53/WT_mean",b)]
  h3=merge(h3,f22)
  
  dd=merge(f11,f12,by.x="gene_Name",by.y="gene_Name")
  d=merge(dd,f22,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
  rownames(d)<-d[,1]
  d1=d[intersect(as.character(f[,1]),as.character(d[,1])),]
  d1[sapply(d1,is.na)]<-0
  x<-d1[,2:ncol(d1)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation","NoMutation","Mutation/NoMutation"), c(length(a1),length(a2),ncol(f2)-1))),row.names=colnames(d)[2:ncol(d)])  
  #pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice2(4),legend_breaks=c(min(y),-1,0,1,max(y)),breaks=c(min(y),-1,0,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_fpkm_foldchange_heatmap/4heatmap1/",gsub("_fc.txt","_M_No_MN.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice1(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/4heatmap/",gsub("_fc.txt","_M_No_MN.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  #pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(20),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_fpkm_foldchange_heatmap/5heatmap/",gsub("_fc.txt","_M_No_MN.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  #write.table(y,file =paste(path,"p53_context_fpkm_foldchange_heatmap/data/",gsub("_fc.txt",".txt",i),sep=""), row.names = T,col.names=T, quote = F,sep = "\t",append=TRUE)
} 
rownames(h1)<-h1[,1]
e1=colnames(sort(h1[grep("^TP53$",as.character(h1[,1])),2:ncol(h1)]))
h11=h1[intersect(as.character(f[,1]),as.character(h1[,1])),c("gene_Name",e1)]
rownames(h2)<-h2[,1]
e2=colnames(sort(h2[grep("^TP53$",as.character(h2[,1])),2:ncol(h2)]))
h22=h2[intersect(as.character(f[,1]),as.character(h2[,1])),c("gene_Name",e2)]
rownames(h3)<-h3[,1]
e3=colnames(sort(h3[grep("^TP53$",as.character(h3[,1])),2:ncol(h3)]))
h33=h3[intersect(as.character(f[,1]),as.character(h3[,1])),c("gene_Name_p53/WT_mean",e3)]
dd1=merge(h11,h22,by.x="gene_Name",by.y="gene_Name")
h=merge(dd1,h33,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
rownames(h)<-h[,1]
m=h[intersect(as.character(f[,1]),as.character(h[,1])),]
m[sapply(m,is.na)]<-0
x1<-m[,2:ncol(m)]
y1<-as.matrix(x1)
row_anno1 = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
col_anno1 = data.frame(SampleClass = factor(rep(c("Mutation","NoMutation","Mutation/NoMutation"), c(length(e1),length(e2),length(e3)))),row.names=colnames(h)[2:ncol(h)])  
pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice1(3),legend_breaks=c(min(y1),-1,1,max(y1)),breaks=c(min(y1),-1,1,max(y1)),annotation_col = col_anno1, annotation_row = row_anno1,fontsize_row=5,fontsize_col=3,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/4heatmap/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
#pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice2(4),legend_breaks=c(min(y1),-1,0,1,max(y1)),breaks=c(min(y1),-1,0,1,max(y1)),annotation_col = col_anno1, annotation_row = row_anno1,fontsize_row=5,fontsize_col=3,
         #file=paste(path,"p53_context_fpkm_foldchange_heatmap/4heatmap1/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  



#ceiling()是上取
#floor() 是下取

#对在mutation、nomutation、mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
##除去PCPG14   列样本进行聚类cluster_col = TRUE  #6heatmap
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))  
h1=data.frame(gene_Name=as.character(f[,1]))
h2=data.frame(gene_Name=as.character(f[,1]))
h3=data.frame(as.character(f[,1]))
colnames(h3)="gene_Name_p53/WT_mean"
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=colnames(sort(f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]))
  a2=colnames(sort(f1["TP53",grep("_no",as.character(colnames(f1)))]))
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a1)]
  f12=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a2)]
  h1=merge(h1,f11)
  h2=merge(h2,f12)
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),c("gene_Name_p53/WT_mean",b)]
  h3=merge(h3,f22)
  
  dd=merge(f11,f12,by.x="gene_Name",by.y="gene_Name")
  d=merge(dd,f22,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
  rownames(d)<-d[,1]
  d1=d[intersect(as.character(f[,1]),as.character(d[,1])),]
  d1[sapply(d1,is.na)]<-0
  x<-d1[,2:ncol(d1)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)), annotation_row = row_anno,fontsize_row=5,fontsize_col=1.5,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
} 
rownames(h1)<-h1[,1]
e1=colnames(sort(h1[grep("^TP53$",as.character(h1[,1])),2:ncol(h1)]))
h11=h1[intersect(as.character(f[,1]),as.character(h1[,1])),c("gene_Name",e1)]
rownames(h2)<-h2[,1]
e2=colnames(sort(h2[grep("^TP53$",as.character(h2[,1])),2:ncol(h2)]))
h22=h2[intersect(as.character(f[,1]),as.character(h2[,1])),c("gene_Name",e2)]
rownames(h3)<-h3[,1]
e3=colnames(sort(h3[grep("^TP53$",as.character(h3[,1])),2:ncol(h3)]))
h33=h3[intersect(as.character(f[,1]),as.character(h3[,1])),c("gene_Name_p53/WT_mean",e3)]
dd1=merge(h11,h22,by.x="gene_Name",by.y="gene_Name")
h=merge(dd1,h33,by.x="gene_Name",by.y="gene_Name_p53/WT_mean")
rownames(h)<-h[,1]
m=h[intersect(as.character(f[,1]),as.character(h[,1])),]
m[sapply(m,is.na)]<-0
x1<-m[,2:ncol(m)]
y1<-as.matrix(x1)
row_anno1 = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y1),-1,1,max(y1)),breaks=c(min(y1),-1,1,max(y1)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  


#对在mutation、nomutation、mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
##除去PCPG14  列样本分开进行聚类cluster_col = TRUE  在mutation、nomutation、mutation/nomutation样本
#6heatmap1"green","white","red"三色 #6heatmap2 "green","white","pink","red"四色分为上调、正常上调、正常下调、下调
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","white","red"))
#colorsChoice<- colorRampPalette(c("green","white","pink","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))  
h1=data.frame(gene_Name=as.character(f[,1]))
h2=data.frame(gene_Name=as.character(f[,1]))
h3=data.frame(as.character(f[,1]))
colnames(h3)="gene_Name_p53/WT_mean"
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=colnames(sort(f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]))
  a2=colnames(sort(f1["TP53",grep("_no",as.character(colnames(f1)))]))
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a1)]
  f12=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name",a2)]
  h1=merge(h1,f11)
  h2=merge(h2,f12)
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),c("gene_Name_p53/WT_mean",b)]
  h3=merge(h3,f22)
  
  x1<-f11[,2:ncol(f11)]
  x2<-f12[,2:ncol(f12)]
  x3<-f22[,2:ncol(f22)]
  y1<-as.matrix(x1)
  y2<-as.matrix(x2)
  y3<-as.matrix(x3)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  #pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y1),-1,0,1,max(y1)),breaks=c(min(y1),-1,0,1,max(y1)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/Mu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  #pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y2),-1,0,1,max(y2)),breaks=c(min(y2),-1,0,1,max(y2)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  #pheatmap(y3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y3),-1,0,1,max(y3)),breaks=c(min(y3),-1,0,1,max(y3)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/Mu_NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y1),-1,1,max(y1)),breaks=c(min(y1),-1,1,max(y1)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/Mu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y2),-1,1,max(y2)),breaks=c(min(y2),-1,1,max(y2)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  pheatmap(y3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y3),-1,1,max(y3)),breaks=c(min(y3),-1,1,max(y3)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/Mu_NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
} 
rownames(h1)<-h1[,1]
rownames(h2)<-h2[,1]
rownames(h3)<-h3[,1]
m1=h1[intersect(as.character(f[,1]),as.character(h1[,1])),2:ncol(h1)]
m2=h2[intersect(as.character(f[,1]),as.character(h2[,1])),2:ncol(h2)]
m3=h3[intersect(as.character(f[,1]),as.character(h3[,1])),2:ncol(h3)]
n1<-as.matrix(m1)
n2<-as.matrix(m2)
n3<-as.matrix(m3)
row_anno1 = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
#pheatmap(n1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(n1),-1,0,1,max(n1)),breaks=c(min(n1),-1,0,1,max(n1)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/Mu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
#pheatmap(n2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(n2),-1,0,1,max(n2)),breaks=c(min(n2),-1,0,1,max(n2)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/NoMu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
#pheatmap(n3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(n3),-1,0,1,max(n3)),breaks=c(min(n3),-1,0,1,max(n3)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         #file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap2/Mu_NoMu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
pheatmap(n1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(n1),-1,1,max(n1)),breaks=c(min(n1),-1,1,max(n1)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/Mu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
pheatmap(n2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(n2),-1,1,max(n2)),breaks=c(min(n2),-1,1,max(n2)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/NoMu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  
pheatmap(n3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(n3),-1,1,max(n3)),breaks=c(min(n3),-1,1,max(n3)),annotation_row = row_anno1,fontsize_row=5,fontsize_col=1.5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap1/Mu_NoMu/Allcancer.pdf",sep=""),width=8,height=8,main="Allcancer")  


#对在mutation、nomutation、mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
##除去PCPG14  列样本分开进行聚类cluster_col = TRUE  在mutation、nomutation、mutation/nomutation样本
#6heatmap1"green","white","red"三色 #6heatmap2 "green","white","pink","red"四色分为上调、正常上调、正常下调、下调
#6heatmap3 四色\对每个cancer行为经过筛选的基因 p值小于0.05\0.04\0.03\0.02\0.01
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","lightgreen","pink","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))  
for(i in filename[c(1:13,15:22)]){
  ff<-read.csv(paste(path,"p53_context_fpkm_foldchange_heatmap/gene__Mu_No_ks_0.01/",gsub("_fc.txt",".txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=colnames(sort(f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]))
  a2=colnames(sort(f1["TP53",grep("_no",as.character(colnames(f1)))]))
  f11=f1[intersect(as.character(f[,1]),as.character(ff[,1])),c("gene_Name",a1)]
  f12=f1[intersect(as.character(f[,1]),as.character(ff[,1])),c("gene_Name",a2)]
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(ff[,1])),c("gene_Name_p53/WT_mean",b)]
  
  x1<-f11[,2:ncol(f11)]
  x2<-f12[,2:ncol(f12)]
  x3<-f22[,2:ncol(f22)]
  y1<-as.matrix(x1)
  y2<-as.matrix(x2)
  y3<-as.matrix(x3)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y1),-1,0,1,max(y1)),breaks=c(min(y1),-1,0,1,max(y1)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap3_0.01/Mu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y2),-1,0,1,max(y2)),breaks=c(min(y2),-1,0,1,max(y2)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap3_0.01/NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
  pheatmap(y3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(4),legend_breaks=c(min(y3),-1,0,1,max(y3)),breaks=c(min(y3),-1,0,1,max(y3)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/6heatmap3_0.01/Mu_NoMu/",gsub("_fc.txt",".pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
} 
#clustering_method 表示聚类方法，值可以是hclust的任何一种，如"ward.D","single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"。

 

#"LAML","LGG","OV","TGCT","UCS"no normal sample
#"TGCT"1mutation/nomutation样本delete
#对这5-1=4cancer在mutation/nomutation样本中样本按 TP53表达值为normal、上下调的顺序排序pheatmap
#列样本进行聚类cluster_col = TRUE 
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in c("LAML","LGG","OV","UCS")){
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  b=colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)]))
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),c("gene_Name_p53/WT_mean",b)]
  
  f22[sapply(f22,is.na)]<-0
  x<-f22[,2:ncol(f22)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep("Mutation/NoMutation",ncol(f22)-1)),row.names=colnames(f22)[2:ncol(f22)])  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/4Nonormal_heatmap/",i,"_MN.pdf",sep=""),width=8,height=8,main=i)  
  pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)), annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_fpkm_foldchange_heatmap/4Nonormal_heatmap/",i,"_MN_col.pdf",sep=""),width=8,height=8,main=i)  
}



#对在mutation、nomutation、mutation/nomutation样本中按gene中主要是上调还是下调（负），标出相应的比例 
##除去PCPG14
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
h=data.frame(gene_Name=as.character(f[,1]))
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  h1=data.frame()
  for(j in 1:nrow(f1)){
    h1[j,1]=as.character(f1[j,1])
    a1=f1[j,-grep("_no",as.character(colnames(f1)))][-1]
    if(length(a1[a1[1,]>=1]) > length(a1[a1[1,]<= -1])){
      h1[j,2]=length(a1[a1[1,]>=1])/ncol(a1)
    }else{
      h1[j,2]=-length(a1[a1[1,]<= -1])/ncol(a1)
    }
    a2=f1[j,grep("_no",as.character(colnames(f1)))]
    if(length(a2[a2[1,]>=1]) > length(a2[a2[1,]<= -1])){
      h1[j,3]=length(a2[a2[1,]>=1])/ncol(a2)
    }else{
      h1[j,3]=-length(a2[a2[1,]<= -1])/ncol(a2)
    }
  }
  colnames(h1)[1]="gene_Name"
  colnames(h1)[2]=gsub("_fc.txt","_Mu",i)
  colnames(h1)[3]=gsub("_fc.txt","_No",i)
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  h2=data.frame()
  for(m in 1:nrow(f2)){
    h2[m,1]=as.character(f2[m,1])
    a3=f2[m,2:ncol(f2)] 
    if(length(a3[a3[1,]>=1]) > length(a3[a3[1,]<= -1])){
      h2[m,2]=length(a3[a3[1,]>=1])/ncol(a3)
    }else{
      h2[m,2]=-length(a3[a3[1,]<= -1])/ncol(a3)
    }
  }
  colnames(h2)[1]="gene_Name"
  colnames(h2)[2]=gsub("_fc.txt","_MuNo",i)

  dd=merge(h1,h2,by.x="gene_Name",by.y="gene_Name")
  h=merge(h,dd)
} 
rownames(h)<-h[,1]
d1=h[intersect(as.character(f[,1]),as.character(h[,1])),]
write.table(d1,file =paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
 


#对在mutation、nomutation、mutation/nomutation样本中按gene中主要是上调还是下调（负），标出相应的比例 
##除去PCPG14 #列样本进行聚类cluster_col = TRUE 
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
rownames(f2)<-f2[,1]
x=f2[intersect(as.character(f[,1]),as.character(f2[,1])),2:ncol(f2)]
y<-as.matrix(x)
row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
col_anno = data.frame(CancerClass = factor(rep("Cancer",ncol(f2)-1)),row.names=colnames(f2)[2:ncol(f2)])  
pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(2),legend_breaks=c(min(y),0,max(y)),breaks=c(min(y),0,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate.pdf",sep=""),width=8,height=8,main="gene_regulate")  
pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(2),legend_breaks=c(min(y),0,max(y)),breaks=c(min(y),0,max(y)), annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate1.pdf",sep=""),width=8,height=8,main="gene_regulate")  

#对在mutation、nomutation、mutation/nomutation样本中按gene中主要是上调还是下调（负），标出相应的比例 
##除去PCPG14 #列样本进行聚类cluster_col = TRUE   
##在gene_regulate中列分别提取cancer的Mu和MuNo，得到gene_regulate_Mu、gene_regulate_MuNo
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("green","red"))     
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
rownames(f2)<-f2[,1]
a1=f2[intersect(as.character(f[,1]),as.character(f2[,1])),grep("_Mu$",as.character(colnames(f2)))]
a2=f2[intersect(as.character(f[,1]),as.character(f2[,1])),grep("_MuNo",as.character(colnames(f2)))]
a3=f2[intersect(as.character(f[,1]),as.character(f2[,1])),grep("_No$",as.character(colnames(f2)))]
y1<-as.matrix(a1)
y2<-as.matrix(a2)
y3<-as.matrix(a3)
row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(2),legend_breaks=c(min(y1),0,max(y1)),breaks=c(min(y1),0,max(y1)),annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_Mu.pdf",sep=""),width=8,height=8,main="gene_regulate")  
pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(2),legend_breaks=c(min(y2),0,max(y2)),breaks=c(min(y2),0,max(y2)),annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_MuNo.pdf",sep=""),width=8,height=8,main="gene_regulate")  
pheatmap(y3,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(2),legend_breaks=c(min(y3),0,max(y3)),breaks=c(min(y3),0,max(y3)),annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_No.pdf",sep=""),width=8,height=8,main="gene_regulate")  
a=cbind(a1,a2,a3)
#a[,sort(colnames(a))]
y<-as.matrix(a[,sort(colnames(a))])
pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(2),legend_breaks=c(min(y),0,max(y)),breaks=c(min(y),0,max(y)),annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_Mu_No_MuNo.pdf",sep=""),width=8,height=8,main="gene_regulate")  




#对在mutation、nomutation、mutation/nomutation样本中按gene中主要是上调还是下调（负），标出相应的比例 
##除去PCPG14 #列样本进行聚类cluster_col = TRUE 
##在gene_regulate中列分别提取cancer的Mu和MuNo，把行标为一个功能得到gene_regulate_Mu1、gene_regulate_MuNo1
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
d=merge(f,f2,by="gene_Name")
d1=d[,c(2,grep("_Mu$",as.character(colnames(d))))]
d2=d[,c(2,grep("_MuNo",as.character(colnames(d))))]
d3=d[,c(2,grep("_No$",as.character(colnames(d))))]
d1[,1]=as.character(d1[,1])
d2[,1]=as.character(d2[,1])
d3[,1]=as.character(d3[,1])

d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(d1[,1])){
  m=d1[d1$Function==i,]
  for(j in 2:ncol(m)){
    if(length(which(m[,j]>=0)) >= length(which(m[,j]<0))){
      d11[i,j]=length(which(m[,j]>=0))/length((m[,j]))
    }else if(length(which(m[,j]>=0)) < length(which(m[,j]<0))){
      d11[i,j]=-length(which(m[,j]<0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(d1)

d22=data.frame(Function=unique(f[,2]))
rownames(d22)<-d22[,1]
for(i in unique(d2[,1])){
  m=d2[d2$Function==i,]
  for(j in 2:ncol(m)){
    if(length(which(m[,j]>=0)) >= length(which(m[,j]<0))){
      d22[i,j]=length(which(m[,j]>=0))/length((m[,j]))
    }else if(length(which(m[,j]>=0)) < length(which(m[,j]<0))){
      d22[i,j]=-length(which(m[,j]<0))/length((m[,j]))
    }
  }
}
colnames(d22)=colnames(d2)

d33=data.frame(Function=unique(f[,2]))
rownames(d33)<-d33[,1]
for(i in unique(d3[,1])){
  m=d3[d3$Function==i,]
  for(j in 2:ncol(m)){
    if(length(which(m[,j]>=0)) >= length(which(m[,j]<0))){
      d33[i,j]=length(which(m[,j]>=0))/length((m[,j]))
    }else if(length(which(m[,j]>=0)) < length(which(m[,j]<0))){
      d33[i,j]=-length(which(m[,j]<0))/length((m[,j]))
    }
  }
}
colnames(d33)=colnames(d3)

y1<-as.matrix(d11[,2:ncol(d11)])
y2<-as.matrix(d22[,2:ncol(d22)])
y3<-as.matrix(d33[,2:ncol(d33)])
pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,fontsize_row=7,fontsize_col=7,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_Mu1.pdf",sep=""),width=8,height=8,main="gene_regulate")  
pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,fontsize_row=7,fontsize_col=7,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_MuNo1.pdf",sep=""),width=8,height=8,main="gene_regulate")  
pheatmap(y3,cluster_row = FALSE,cluster_col = TRUE,fontsize_row=7,fontsize_col=7,
         file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_No1.pdf",sep=""),width=8,height=8,main="gene_regulate")  


library(ggplot2) 
da1=data.frame()
for(i in 1:nrow(d11)){
  da=data.frame()
  for(j in 2:ncol(d11)){
    da[j-1,1]=rownames(d11)[i]
    da[j-1,2]=colnames(d11)[j]
    da[j-1,3]=d11[i,j]
  }
  da1=rbind(da1,da)
}

da2=data.frame()
for(i in 1:nrow(d22)){
  da=data.frame()
  for(j in 2:ncol(d22)){
    da[j-1,1]=rownames(d22)[i]
    da[j-1,2]=colnames(d22)[j]
    da[j-1,3]=d22[i,j]
  }
  da2=rbind(da2,da)
}

da3=data.frame()
for(i in 1:nrow(d33)){
  da=data.frame()
  for(j in 2:ncol(d33)){
    da[j-1,1]=rownames(d33)[i]
    da[j-1,2]=colnames(d33)[j]
    da[j-1,3]=d33[i,j]
  }
  da3=rbind(da3,da)
}

da1[da1[,3]>=0,4]="up"
da1[da1[,3]<0,4]="down"
da1[,5]=abs(da1[,3])
da1[,1]=factor(da1[,1],unique(f[,2]))
#da1[,2]=gsub("_Mu","",da1[,2])
da1[,2]=gsub("_Mu","2",da1[,2])
colnames(da1)=c("Function","Cancer","Exp","State","Expression")
p <- ggplot(da1,aes(x=factor(1),y=da1[,5],fill=State))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da1[,1]~da1[,2])+
  scale_fill_manual(values=c("green","red"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_Mu2.pdf",sep=""),width = 10,height=9)
p <- ggplot(da1,aes(x=factor(1),y=da1[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da1[,1]~da1[,2])+
  scale_fill_gradient(low="#12B826",high = "#FF331C")  #"#12B826","#FF331C"
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_Mu3.pdf",sep=""),width = 10,height=9)
 
da2[da2[,3]>=0,4]="up"
da2[da2[,3]<0,4]="down"
da2[,5]=abs(da2[,3])
da2[,1]=factor(da2[,1],unique(f[,2]))
#da2[,2]=gsub("_MuNo","",da2[,2])
da2[,2]=gsub("_MuNo","3",da2[,2])
colnames(da2)=c("Function","Cancer","Exp","State","Expression")
p <- ggplot(da2,aes(x=factor(1),y=da2[,5],fill=State))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da2[,1]~da2[,2])+
  scale_fill_manual(values=c("green","red"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_MuNo2.pdf",sep=""),width = 10,height=9)
p <- ggplot(da2,aes(x=factor(1),y=da2[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da2[,1]~da2[,2])+
  scale_fill_gradient(low="#12B826",high = "#FF331C")
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_MuNo3.pdf",sep=""),width = 10,height=9)

da3[da3[,3]>=0,4]="up"
da3[da3[,3]<0,4]="down"
da3[,5]=abs(da3[,3])
da3[,1]=factor(da3[,1],unique(f[,2]))
#da3[,2]=gsub("_No","",da3[,2])
da3[,2]=gsub("_No","1",da3[,2])
colnames(da3)=c("Function","Cancer","Exp","State","Expression")
p <- ggplot(da3,aes(x=factor(1),y=da3[,5],fill=State))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da3[,1]~da3[,2])+
  scale_fill_manual(values=c("green","red"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_No2.pdf",sep=""),width = 10,height=9)
p <- ggplot(da3,aes(x=factor(1),y=da3[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da3[,1]~da3[,2])+
  scale_fill_gradient(low="#12B826",high = "#FF331C")
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_No3.pdf",sep=""),width = 10,height=9)
##整合Mu为2，No为1，MuNo为3
m=rbind(da1,da2,da3)
m1=m[grep("^BLCA.*|^BRCA.*|^CESC.*|^COAD.*|^ESCA.*|^GBM.*|^HNSC.*",m[,2]),]
m2=m[grep("^KIRC.*|^KIRP.*|^LIHC.*|^LUAD.*|^LUSC.*|^PAAD.*|^PRAD.*",m[,2]),]
m3=m[grep("^READ.*|^SARC.*|^SKCM.*|^STAD.*|^THCA.*|^THYM.*|^UCEC.*",m[,2]),]
p <- ggplot(m1,aes(x=factor(1),y=m1[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(m1[,1]~m1[,2])+
  scale_fill_gradient(low="green",high = "red")
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_all1.pdf",sep=""),width = 10,height=9)
p <- ggplot(m2,aes(x=factor(1),y=m2[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(m2[,1]~m2[,2])+
  scale_fill_gradient(low="green",high = "red")
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_all2.pdf",sep=""),width = 10,height=9)
p <- ggplot(m3,aes(x=factor(1),y=m3[,5],fill=Exp))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(m3[,1]~m3[,2])+
  scale_fill_gradient(low="green",high = "red")
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("Function") + ylab("Cancer") + ggtitle("Cancer_Function")
ggsave(file=paste(path,"p53_context_fpkm_foldchange_heatmap/gene_regulate_all3.pdf",sep=""),width = 10,height=9)




#对gene在mutation、nomutation样本中的foldchange值做ks检验，p值小于0.05
##除去PCPG14
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(gene_Name=as.character(f[,1]))
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  h1=data.frame()
  for(j in 1:nrow(f1)){
    h1[j,1]=as.character(f1[j,1])
    a1=f1[as.character(f1[j,1]),-grep("_no",as.character(colnames(f1)))][-1]
    a2=f1[as.character(f1[j,1]),grep("_no",as.character(colnames(f1)))]
    if(ks.test(as.numeric(a1),as.numeric(a2))$p.value<=0.05){
      h1[j,2]=ks.test(as.numeric(a1),as.numeric(a2))$p.value
    }else{
      h1[j,2]=NA
    }
  }
  colnames(h1)[1]="gene_Name"
  colnames(h1)[2]=gsub("_fc.txt","_Mu_No_ks",i)
  h=merge(h,h1)
}
rownames(h)<-h[,1]
d1=h[intersect(as.character(f[,1]),as.character(h[,1])),]
write.table(d1,file =paste(path,"p53_context_fpkm_foldchange_heatmap/gene__Mu_No_ks.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


#对gene在mutation、nomutation样本中的foldchange值做ks检验，p值小于0.05找出对应的每个cancer的基因
##除去PCPG14  p值小于0.05\0.04\0.03\0.02\0.01
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  for(j in 1:nrow(f1)){
    a1=f1[as.character(f1[j,1]),-grep("_no",as.character(colnames(f1)))][-1]
    a2=f1[as.character(f1[j,1]),grep("_no",as.character(colnames(f1)))]
    if(ks.test(as.numeric(a1),as.numeric(a2))$p.value<=0.05){
      write.table(as.character(f1[j,1]),file =paste(path,"p53_context_fpkm_foldchange_heatmap/gene__Mu_No_ks_0.05/",gsub("_fc.txt",".txt",i),sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
  }
}
 


#######################################3
##对每个癌症22的每个样本在每个基因中有3个表达值，分别为Mu/normal,Mu/nomu,nomu/normal的均值（已定），
##根据3个值的表达，可以把样本标为红红红1，红红绿2，绿红绿3，绿绿绿-1，绿绿红-2，红绿红-3，这6种模式+0为其他。
##14PCPG1Mu 删去=21cancer  fc=0.5\2(log2=-1\1)  fc=2/3\1.5(log2=-0.5849625\0.5849625)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(cancer in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",cancer,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=f1[intersect(as.character(f[,1]),as.character(f1[,1])),-grep("_no",as.character(colnames(f1)))][,-1]
  a2=f1[intersect(as.character(f[,1]),as.character(f1[,1])),grep("_no",as.character(colnames(f1)))]
  a3=data.frame(NoMu=apply(a2,1,mean))
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",cancer,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]  
  a4=f2[intersect(as.character(f[,1]),as.character(f2[,1])),-1]  
  
  da=data.frame(matrix(NA,nrow(a1),ncol(a1)))
  rownames(da)=rownames(a1)  
  colnames(da)=colnames(a1)  
  for(i in 1:nrow(da)){
    for(j in 1:ncol(da)){
      #da[i,j]
      #rownames(da)[i]
      #colnames(da)[j]
      #b1=a1[rownames(da)[i],colnames(da)[j]]
      #b2=a4[rownames(da)[i],paste(colnames(da)[j],"_p53/WT_mean",sep="")]
      #b3=a3[rownames(da)[i],1] #a1[rownames(da)[i],colnames(da)[j]]>=1
      if(a1[rownames(da)[i],colnames(da)[j]]>=0.5849625){
        b1=1
      }else if(a1[rownames(da)[i],colnames(da)[j]]<=-0.5849625){
        b1=-1
      }else if((a1[rownames(da)[i],colnames(da)[j]]>-0.5849625) & (a1[rownames(da)[i],colnames(da)[j]]<0.5849625)){
        b1=0
      }
      if(a4[rownames(da)[i],paste(colnames(da)[j],"_p53/WT_mean",sep="")]>=0.5849625){
        b2=1
      }else if(a4[rownames(da)[i],paste(colnames(da)[j],"_p53/WT_mean",sep="")]<=-0.5849625){
        b2=-1
      }else if((a4[rownames(da)[i],paste(colnames(da)[j],"_p53/WT_mean",sep="")]>-0.5849625) & (a4[rownames(da)[i],paste(colnames(da)[j],"_p53/WT_mean",sep="")]<0.5849625)){
        b2=0
      }
      if(a3[rownames(da)[i],1]>=0.5849625){
        b3=1
      }else if(a3[rownames(da)[i],1]<=-0.5849625){
        b3=-1
      }else if((a3[rownames(da)[i],1]>-0.5849625) & (a3[rownames(da)[i],1]<0.5849625)){
        b3=0
      }
      #b=c(b1,b2,b3)
      #print(b)
      if(TRUE %in% ((b1==1) & (b2==1) & (b3==1))){
        da[i,j]=1
      }else if(TRUE %in% ((b1==1) & (b2==1) & (b3==-1))){
        da[i,j]=5
      }else if(TRUE %in% ((b1==-1) & (b2==1) & (b3==-1))){
        da[i,j]=10
      }else if(TRUE %in% ((b1==-1) & (b2==-1) & (b3==-1))){
        da[i,j]=-1
      }else if(TRUE %in% ((b1==-1) & (b2==-1) & (b3==1))){
        da[i,j]=-5
      }else if(TRUE %in% ((b1==1) & (b2==-1) & (b3==1))){
        da[i,j]=-10
      }else{
        da[i,j]=0 
      }
      
    }
  }  
  da1=cbind(name=rownames(da),da)
  write.table(da1,file =paste(path,"p53_context_discreteMode_fc1.5/1510/",gsub("_fc","",cancer),sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
}
########################################

##根据p53_context_discreteMode对每个癌症做热图  
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  da=f1[intersect(as.character(f[,1]),as.character(f1[,1])),-1]
  da1=na.omit(da)
  y<-as.matrix(da1)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  #pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="euclidean",clustering_method = "ward.D",
           #col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           #file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",gsub(".txt",".pdf",cancer),sep=""),width=8,height=8,main=gsub(".txt","",cancer))  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])
y1<-as.matrix(h1)
pheatmap(y1,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="euclidean",clustering_method = "ward.D",
         col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
         file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/Allcancer1.pdf",sep=""),width=8,height=8,main="Allcancer")  

for(number in 2:60){
  dir.create(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",number,sep="")) 
  fu<-function(x) sum(x!=0) ##删去0比较多的基因337-202-120\112\101
  a1=data.frame(apply(h1,1,fu))
  a2=rownames(a1)[which(a1[,1]>=number)] ##10\15\20
  h2=na.omit(h[a2,-1])
  y2<-as.matrix(h2)
  pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="manhattan",clustering_method = "ward.D",
           col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,fontsize_row=3,show_colnames=FALSE,
           file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",number,"/Allcancer",number,"_manhattan_ward.D.pdf",sep=""),width=8,height=8,main="Pancancer")  
  
  pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="euclidean",clustering_method = "ward.D",
           col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",number,"/Allcancer",number,".pdf",sep=""),width=8,height=8,main="Allcancer")  
  for(nu in 2:10){
    d = dist(t(h2), method = "euclidean")
    hcward = hclust(d, method="ward.D")
    n= cutree(hcward,k=nu) 
    f=cbind(V1=rownames(as.data.frame(n)),as.data.frame(n))
    e2=na.omit(merge(f,e1,by.x="V1",by.y="X1",all.x=T))
    dif1 <- survdiff(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2])#求生存时间
    kmsurvival1<-survfit(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2],conf.type = "log-log")
    if(round(pchisq(dif1$chisq,1,lower.tail=F),4)<=0.05){
      pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/1510/",number,"/Allcancer_class",nu,"_survival.pdf",sep=""))  
      plot(kmsurvival1, lty = 'solid', col=c("red","orange","blue","purple","burlywood4","deeppink","gold","green","salmon4","black"),lwd=1.1,
           xlab='survival time in days',ylab='survival probabilities',main=paste("Allcancer_class",nu,sep=""))
      legend('bottomleft', cex=0.6,text.width=0.4,gsub(".*=","",as.character(data.frame(dif1$n)[,1])), lty='solid',
             col=c("red","orange","blue","purple","burlywood4","deeppink","gold","green","salmon4","black"))
      text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
      for(j in 1:length(as.numeric(unique(e2[,2]))) ){
        text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
      }
      dev.off() 
      write.table(as.data.frame(n),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/1510/",number,"/Allcancer_class",nu,".txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)
      
    }
  }
}


number=36
fu<-function(x) sum(x!=0) ##删去0比较多的基因337-202-120\112\101
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] ##10\15\20
h2=na.omit(h[a2,-1])
#c("manhattan","euclidean","minkowski","canberra")
#c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
##对所有癌症hclust，2,3,4,看每类的生存
d = dist(t(h2), method = "manhattan")
hcward = hclust(d, method="complete")
n= cutree(hcward,k=3)  
write.table(as.data.frame(n),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/36/Allcancer36_class3_manhattan_complete.txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)

pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/1/Allcancer_class4.pdf",sep=""))
plot(hcward,cex = 0.2,hang=-1)   #hang是表明谱系图中各类所在的位置 当hang取负值时，谱系图中的类从底部画起  生成谱系图
re<-rect.hclust(hcward,k=4,border="red") 
dev.off()


a= read.table("D:/go_all.txt",header=T,sep="\t")
rownames(a)=a[,1]
d = dist(t(a[,2:3120]), method = "manhattan")
hc= hclust(d, method="complete")
plot(hc,cex = 0.2,hang=-1)

n= cutree(hc,k=3)
write.table(as.data.frame(n),
            file =paste(path,"p53/Allcancer36_class3_manhattan_complete.txt",sep=""), 
            row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)

##对Allcancer26_class3_manhattan_ward.D 
##Allcancer2_class4_manhattan_ward.D 
##Allcancer2_class4_minkowski_ward.D2 
library("survival")
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
b11=subset(ff,ff$vital_status.y==0)
b21=subset(ff,ff$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(data.frame(e1))

path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",sep=""),header=F,sep="\t")
e2=na.omit(merge(f,e1,by.x="V1",by.y="X1",all.x=T))
e2[,2]=factor(e2[,2])
e2[,3]=as.numeric(as.character(e2[,3]))
e2[,4]=as.numeric(as.character(e2[,4]))
pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_survival.pdf",sep=""))  
dif1 <- survdiff(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2])#求生存时间
kmsurvival1<-survfit(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2],conf.type = "log-log")
plot(kmsurvival1, lty = 'solid', col=c("deeppink3","goldenrod3","palegreen4","royalblue4"),lwd=1.2,
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer2_class4_manhattan_ward.D")
legend('bottomleft', cex=0.6,text.width=0.6,gsub(".*=","",as.character(data.frame(dif1$n)[,1])), lty='solid',lwd=1.2,
       col=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
for(j in 1:length(as.numeric(unique(e2[,2]))) ){
  text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
}
for(jj in 1:length(as.numeric(unique(e2[,2]))) ){
  text(2000,0.9-(jj-1)*0.03,cex=0.8,round(mean(e2[as.numeric(as.character(e2[,2]))==jj,4]),1))
}
dev.off() 


##对Allcancer26_class3_manhattan_ward.D
##Allcancer2_class4_manhattan_ward.D
##Allcancer2_class4_minkowski_ward.D2
##对所有癌症hclust，2,3,4,看每类的生存时间的差异boxplot
library(ggplot2)
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
b11=subset(ff,ff$vital_status.y==0)
b21=subset(ff,ff$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(data.frame(e1))

path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",sep=""),header=F,sep="\t")
e2=na.omit(merge(f1,e1,by.x="V1",by.y="X1",all.x=T))
e2[,2]=factor(e2[,2])
e2[,4]=as.numeric(as.character(e2[,4]))
#write.table(e2[,1:2],
            #file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

a1=e2[as.numeric(as.character(e2[,2]))==1,4]
a2=e2[as.numeric(as.character(e2[,2]))==2,4]
a3=e2[as.numeric(as.character(e2[,2]))==3,4]
a4=e2[as.numeric(as.character(e2[,2]))==4,4]

p<-ggplot(data=e2, aes(x = e2[,2],y = e2[,4]))+geom_boxplot(aes(fill=X2),outlier.colour = NA)+
  scale_fill_manual(values=c("pink","lightgreen"))+geom_jitter() 
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("class_state") + ylab("time") + ggtitle("Allcancer_class") 
p <- p+annotate("text",x=2.5, y=7500,colour="red",size=3,label=paste(round(mean(e2[as.numeric(as.character(e2[,2]))==1,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==2,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==3,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==4,4]),1),sep="_"))+
  annotate("text",x=1.5, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a2)$p.value,3))+
  annotate("text",x=2, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a3)$p.value,3))+
  annotate("text",x=2.5, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a4)$p.value,3))+
  annotate("text",x=2.5, y=max(e2[,4])-500,colour="red",size=3,label=round(t.test(a2,a3)$p.value,3))+
  annotate("text",x=3, y=max(e2[,4])-500,colour="red",size=3,label=round(t.test(a2,a4)$p.value,3))+
  annotate("text",x=3.5, y=max(e2[,4])-1000,colour="red",size=3,label=round(t.test(a3,a4)$p.value,3))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_boxplot.pdf",sep=""),width = 10,height=9)

p<-ggplot(data=e2, aes(x = e2[,2],y = e2[,4]))+geom_violin(aes(fill=V2),alpha=0.6)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("class_state") + ylab("time") + ggtitle("Allcancer_class") 
p <- p+annotate("text",x=2.5, y=7500,colour="red",size=3,label=paste(round(mean(e2[as.numeric(as.character(e2[,2]))==1,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==2,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==3,4]),1),round(mean(e2[as.numeric(as.character(e2[,2]))==4,4]),1),sep="_"))+
  annotate("text",x=1.5, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a2)$p.value,3))+
  annotate("text",x=2, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a3)$p.value,3))+
  annotate("text",x=2.5, y=max(e2[,4]),colour="red",size=3,label=round(t.test(a1,a4)$p.value,3))+
  annotate("text",x=2.5, y=max(e2[,4])-500,colour="red",size=3,label=round(t.test(a2,a3)$p.value,3))+
  annotate("text",x=3, y=max(e2[,4])-500,colour="red",size=3,label=round(t.test(a2,a4)$p.value,3))+
  annotate("text",x=3.5, y=max(e2[,4])-1000,colour="red",size=3,label=round(t.test(a3,a4)$p.value,3))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_minkowski_ward.D2_violin.pdf",sep=""),width = 10,height=9)



############对Allcancer分类及cancer 
library(ggplot2)
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
b11=subset(ff,ff$vital_status.y==0)
b21=subset(ff,ff$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$project_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$project_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(data.frame(e1))

path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",sep=""),header=F,sep="\t")
e2=na.omit(merge(f1,e1,by.x="V1",by.y="X1",all.x=T))
e2[,2]=factor(e2[,2])
e2[,5]=as.numeric(as.character(e2[,5]))
write.table(e2[,1:3],
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",sep=""),header=F,sep="\t")
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff[,1]=gsub("-",".",ff[,1])
e1=ff[,1:2]
e2=na.omit(merge(f1,e1,by.x="V1",by.y="submitter_id",all.x=T))
write.table(e2,
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)





#######################################
#position:dodge并列；stack堆叠默认；fill方式和stack类似，但Y轴不再是计数，而是以百分比显示；identity方式是不做任何改变直接显示出来，所以需要设置透明度才能看得清楚。
#geom_bar函数中加了一个参数stat="identity"，表示不要像之前一样去查数，而是就使用数据Freq本身作为频数
##对Allcancer26_class3_manhattan_ward.D
##Allcancer2_class4_manhattan_ward.D
##Allcancer2_class4_minkowski_ward.D2
#install.packages("zoo")
library("zoo")
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
a=ff[,1:2]
a[,1]=gsub("-",".",a[,1])
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_minkowski_ward.D2_1075.txt",sep=""),header=F,sep="\t")
a1=merge(f1,a,by.x="V1",by.y="submitter_id",all.x=T)
a2=data.frame(table(a1[,2:3]))
a3=data.frame(table(a1[,3]))
a4=merge(a2,a3,by.x="project_id",by.y="Var1",all.x=T)
a4[,5]=a4[,3]/a4[,4]
a4[,2]=factor(a4[,2],c(4,3,2,1))
p<-ggplot(data=a4, aes(x = a4[,1],y = a4[,5],fill=V2))+geom_bar(stat="identity",alpha=0.5)+
  scale_fill_manual(values=c("royalblue4","palegreen4","goldenrod3","deeppink3"))
  #scale_fill_manual(values=c("royalblue4","palegreen4","goldenrod3","deeppink3"))
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Cancer") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of different classes in Cancers")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_minkowski_ward.D2_barplot.pdf",sep=""),width = 10,height=9)

b=rollapply(a4[,5],4,cumsum)
for(i in c(seq(1,nrow(a4),4))){
  a4[i,6]=b[i,1]
  a4[i+1,6]=b[i,2]
  a4[i+2,6]=b[i,3]
  a4[i+3,6]=b[i,4]
}
p<-ggplot(data=a4, aes(x = a4[,1],y = a4[,5],fill=V2))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Cancer") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of different classes in Cancers")+
  annotate("text",x=a4[,1],y=a4[,6], colour="black",size=2.5,label=paste(round(a4[,5]*100,1),"%",sep=""))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_barplot1.pdf",sep=""),width = 10,height=9)


##对Allcancer26_class3_manhattan_ward.D_1075对每类计算轮廓系数silhouette。
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075
library(ConsensusClusterPlus)
library(fpc) # sil
library(cluster)
library(factoextra)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
h2=h1[,as.character(f2[,1])]
number=26
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h2,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h3=na.omit(h2[a2,])
nu=3
m1="manhattan" 
m2="ward.D"
d = dist(t(h3), method =m1)
hcward = hclust(d, method=m2)
sil <- silhouette(f2[,2], d)
pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_silhouette.pdf",sep=""))  
plot(sil,main ="Allcancer26_class3_manhattan_ward.D_silhouette",col=c(terrain.colors(10)[c(1,5,7,9,3)],topo.colors(10)[c(2,4,6,8,10)])[1:nu])
abline(v = mean(sil[,3]),lty = 2,col = "gray30")
dev.off() 



##对Allcancer26_class3_manhattan_ward.D_1075对每类加表型堆叠柱状图
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075
library("zoo")
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
a=ff[,1:2]
a[,1]=gsub("-",".",a[,1])
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
a1=merge(f1,a,by.x="V1",by.y="submitter_id",all.x=T)
da=data.frame()
for(i in unique(a1[,3])){
  a2=a1[a1$project_id==i,]
  f2= read.csv(paste("D:/tp53/GDC_phenotype/",i,"_Mutation.txt",sep=""),header=T,sep="\t")
  f3=f2[grep(".*-01A$",as.character(f2[,1])),]
  f3[,1]=gsub("-",".",substring(as.character(f3[,1]),0,12)) #"neoplasm_histologic_grade", 
  f4=f3[,c("submitter_id.samples","tumor_stage.diagnoses","gender.demographic","race.demographic","age_at_initial_pathologic_diagnosis","days_to_new_tumor_event_after_initial_treatment","initial_weight.samples")]
  a3=unique(merge(a2,f4,by.x="V1",by.y="submitter_id.samples",all.x=T))
  da=rbind(da,a3)
}
##离散型 "tumor_stage.diagnoses","gender.demographic","race.demographic"
for(j in 4:6){
  da1=da[,c(2,j)]
  d1=data.frame(table(da1))
  d11=d1[nchar(as.character(d1[,2]))!=0,]
  d11[,2]=factor(d11[,2],unique(sort(as.character(d11[,2]))))
  a4=data.frame(table(da1[!is.na(as.character(da1[,2])),1]))
  a5=merge(d11,a4,by.x="V2",by.y="Var1",all.x=T)
  a5[,5]=a5[,3]/a5[,4]
  a6=na.omit(a5)
  a7=a6[order(as.character(a6[,1]),as.character(a6[,2]),decreasing=F),]
  p<-ggplot(data=a7, aes(x = a7[,1],y = a7[,5],fill=a7[,2]))+geom_bar(stat="identity",alpha=0.5)
  p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle(colnames(d11)[2])
  ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/",colnames(d11)[2],".pdf",sep=""),width = 10,height=9)
  b=rollapply(a7[,5],length(unique(a7[,2])),cumsum)
  for(i in c(seq(1,nrow(a7),length(unique(a7[,2]))))){
    for(j in 1:length(unique(a7[,2]))){
      a7[i+j-1,6]=b[i,j]
    }
  }
  p<-ggplot(data=a7, aes(x = a7[,1],y = a7[,5],fill=a7[,2]))+geom_bar(stat="identity",alpha=0.5)
  p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle(colnames(d11)[2])+
    annotate("text",x=a7[,1],y=1-a7[,6], colour="black",size=2.5,label=paste(round(a7[,5]*100,1),"%",sep=""))
  ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/",colnames(d11)[2],"1.pdf",sep=""),width = 10,height=9)
}
##连续型"age_at_initial_pathologic_diagnosis","days_to_new_tumor_event_after_initial_treatment","initial_weight.samples"
for(j in 7:ncol(da)){
  da1=na.omit(da[,c(2,j)])
  #c1=da1[da1[,1]==1,2]
  #c2=da1[da1[,1]==2,2]
  #c3=da1[da1[,1]==3,2]
  #c4=da1[da1[,1]==4,2]
  #signif(t.test(c3,c4)$p.value,2)
  p<-ggplot(data=da1, aes(x = da1[,1],y = da1[,2]))+geom_violin(aes(fill=factor(V2)),alpha=0.6)+geom_boxplot(aes(fill=factor(V2)),alpha=0.2,width=0.1,outlier.colour = NA)+
    scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
  p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Class") + ylab(colnames(da1)[2]) + ggtitle(colnames(da1)[2])
  p <- p+annotate("text",x=2.5, y=max(da1[,2])+1,colour="black",size=3,label=paste(round(mean(da1[da1[,1]==1,2]),1),round(mean(da1[da1[,1]==2,2]),1),round(mean(da1[da1[,1]==3,2]),1),round(mean(da1[da1[,1]==4,2]),1),sep="_"))
  ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/",colnames(da1)[2],".pdf",sep=""),width = 10,height=9)
}


##离散型da$tumor_stage.diagnoses
j=4
daa=da[,c(2,j)]
daa[,3]=as.character(daa[,2])
daa[which(as.character(daa$tumor_stage.diagnoses)=="not reported"),3]=NA
daa[,3]=gsub("[a-c]$","",as.character(daa[,3]))
da1=daa[,c(1,3)]

d1=data.frame(table(da1))
d11=d1[nchar(as.character(d1[,2]))!=0,]
d11[,2]=factor(d11[,2],unique(sort(as.character(d11[,2]))))
a4=data.frame(table(da1[!is.na(as.character(da1[,2])),1]))
a5=merge(d11,a4,by.x="V2",by.y="Var1",all.x=T)
a5[,5]=a5[,3]/a5[,4]
a6=na.omit(a5)
a7=a6[order(as.character(a6[,1]),as.character(a6[,2]),decreasing=F),]
p<-ggplot(data=a7, aes(x = a7[,1],y = a7[,5],fill=a7[,2]))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle(colnames(d11)[2])
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D2/",colnames(d11)[2],"_stage.pdf",sep=""),width = 10,height=9)
b=rollapply(a7[,5],length(unique(a7[,2])),cumsum)
for(i in c(seq(1,nrow(a7),length(unique(a7[,2]))))){
  for(j in 1:length(unique(a7[,2]))){
    a7[i+j-1,6]=b[i,j]
  }
}
p<-ggplot(data=a7, aes(x = a7[,1],y = a7[,5],fill=a7[,2]))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle(colnames(d11)[2])+
  annotate("text",x=a7[,1],y=1-a7[,6], colour="black",size=2.5,label=paste(round(a7[,5]*100,1),"%",sep=""))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D2/",colnames(d11)[2],"_stage1.pdf",sep=""),width = 10,height=9)



##对Allcancer26_class3_manhattan_ward.D_1075对每类及表型加HR
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075
library("survival")
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
b11=subset(ff,ff$vital_status.y==0)
b21=subset(ff,ff$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(data.frame(e1))

a=ff[,1:2]
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
a1=merge(f1,a,by.x="V1",by.y="submitter_id",all.x=T)
e2=merge(a1,e1,by.x="V1",by.y="X1",all.x=T)
da=data.frame()
for(i in unique(e2[,3])){
  a2=e2[e2$project_id==i,]
  f2= read.csv(paste("D:/tp53/GDC_phenotype/",i,"_Mutation.txt",sep=""),header=T,sep="\t")
  f3=f2[grep(".*-01A$",as.character(f2[,1])),]
  f3[,1]=gsub("-",".",substring(as.character(f3[,1]),0,12)) #"neoplasm_histologic_grade", 
  f4=f3[,c("submitter_id.samples","tumor_stage.diagnoses","gender.demographic","age_at_initial_pathologic_diagnosis","days_to_new_tumor_event_after_initial_treatment","initial_weight.samples")]
  a3=unique(merge(a2,f4,by.x="V1",by.y="submitter_id.samples",all.x=T))
  da=rbind(da,a3)
}
da[,4]=as.numeric(as.character(da[,4]))
da[,5]=as.numeric(as.character(da[,5]))
write.table(paste("class","number","HR","lower95","high95","p",sep="\t"),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
for(i in sort(unique(da$V2))){ #Class
  e3=da[,2:5]
  e3[,5]=0
  e3[e3$V2==i,5]=i
  cox1=coxph(Surv(e3[,4],e3[,3])~e3[,5])
  write.table(paste(i,nrow(e3[e3$V2==i,]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
              file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
for(i in as.character(sort(unique(da$gender.demographic)))){ #male VS female
  e3=na.omit(da[,c("X2","X3","gender.demographic")])
  e3[,1]=as.numeric(as.character(e3[,1]))
  e3[,2]=as.numeric(as.character(e3[,2]))
  e3[,3]=as.character(e3[,3])
  cox1=coxph(Surv(e3[,2],e3[,1])~e3[,3])
  write.table(paste(i,nrow(e3[e3[,3]==i,]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
              file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
e4=na.omit(da[,c("X2","X3","age_at_initial_pathologic_diagnosis")])
cox1=coxph(Surv(e4[,2],e4[,1])~e4[,3])
write.table(paste(colnames(e4)[3],nrow(e4),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
e5=na.omit(da[,c("X2","X3","days_to_new_tumor_event_after_initial_treatment")])
cox2=coxph(Surv(e5[,2],e5[,1])~e5[,3])
write.table(paste(colnames(e5)[3],nrow(e5),summary(cox2)$coefficients[,2],summary(cox2)$conf.int[,3],summary(cox2)$conf.int[,4],summary(cox2)$coefficients[,5],sep="\t"),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
e6=na.omit(da[,c("X2","X3","initial_weight.samples")])
cox3=coxph(Surv(e6[,2],e6[,1])~e6[,3])
write.table(paste(colnames(e6)[3],nrow(e6),summary(cox3)$coefficients[,2],summary(cox3)$conf.int[,3],summary(cox3)$conf.int[,4],summary(cox3)$coefficients[,5],sep="\t"),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
e7=na.omit(da[,c("X2","X3","tumor_stage.diagnoses")])
e7[,4]=gsub("[a-c]$","",as.character(e7[,3]))
e7[which(as.character(e7$tumor_stage.diagnoses)=="not reported"),4]=NA
e8=na.omit(e7[,c(1,2,4)])
cox4=coxph(Surv(e8[,2],e8[,1])~e8[,3])
write.table(cbind(colnames(e7)[3],rownames(data.frame(summary(cox4)$coefficients)),data.frame(HR=summary(cox4)$coefficients[,2]),data.frame(lower95=summary(cox4)$conf.int[,3]),data.frame(high95=summary(cox4)$conf.int[,4]),data.frame(p=summary(cox4)$coefficients[,5])),
            file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



##对Allcancer26_class3_manhattan_ward.D_1075对每类HR加boxplot
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/class1.txt",sep=""),header=T,sep="\t")
f[,1]=as.character(f[,1])
f1=f[f$class!="female",]
f1[,7]=1:nrow(f1)
pdf(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/classHR.pdf",sep=""),width = 10,height=9)
plot(f1[,c(3,7)],pch=15,col= "black",cex = 2,las=1,xlim=c(min(f1[,4])-0.2,max(f1[,5])+0.6),ylim=c(min(f1[,7])-0.2,max(f1[,7])+0.2),
     yaxt="n",main ="2/D2",xlab="Hazard Ratio(95% CI)",ylab="Class")
abline(v=1,lwd=1,col="grey",lty=2)#虚线
text(f1[,5]+0.2,f1[,7],cex=0.8,round(f1[,6],4))
text(f1[,4]-0.25,f1[,7],cex=0.8,paste(substring(f1[,1],0,4),"=",gsub(".*stage","",as.character(f1[,2])),sep=""))
text(f1[,5]+0.5,f1[,7],cex=0.8,paste(round(f1[,3],2),"(",round(f1[,4],2),",",round(f1[,5],2),")",sep=""))
segments(f1[,4],f1[,7],f1[,5],f1[,7],col="black",lwd=2)
dev.off() 
#abline(v=3,lwd=4,col="blue")#添加一条垂直直线x=3，线宽为4，颜色蓝色
#abline(h=3,lwd=4,col="blue")#添加一条水平直线y=3，线宽为4，颜色蓝色
#abline(lm(y~x), lwd=4, col="red")#添加一条一元线性回归拟合直线，线宽为4，颜色为红色

 


##对Allcancer26_class3_manhattan_ward.D_1075对每类加免疫亚型246boxplot
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075 
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
a=ff[,1:2]
a[,1]=gsub("-",".",a[,1])
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
a1=merge(f,a,by.x="V1",by.y="submitter_id",all.x=T)
a2=data.frame()
for(i in unique(a1$project_id)){
  if(TRUE %in% (paste("AGP-",tolower(i),".txt",sep="") %in% list.files("D:/tp53/AGPall/"))){
    f1= read.csv(paste("D:/tp53/AGPall/AGP-",tolower(i),".txt",sep=""),header=T,sep="\t")
    f2=f1[grep(".*01A$",as.character(f1$sampleid)),c("sampleid","purity","Euploidy")]
    f2[,1]=substring(as.character(f2[,1]),0,12)
    a2=rbind(a2,f2)
  }
}
a3=merge(a1,a2,by.x="V1",by.y="sampleid",all.x=T)
a3[,2]=factor(a3[,2])
a3=a3[!is.na(a3$Euploidy),]
p<-ggplot(data=a3, aes(x = a3$V2,y = a3$purity,fill=a3$V2))+geom_boxplot(alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+facet_grid(a3$Euploidy ~ .)
p <- p+theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("purity") + ggtitle("Different Immune Cell Types")  
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/immun246.pdf",sep=""),width = 10,height=9)
 
for(j in c(2,4,6)){
  a4=a3[a3$Euploidy==j,]
  b1=a4[as.numeric(as.character(a4$V2))==1,4]
  b2=a4[as.numeric(as.character(a4$V2))==2,4]
  b3=a4[as.numeric(as.character(a4$V2))==3,4]
  b4=a4[as.numeric(as.character(a4$V2))==4,4]
  p<-ggplot(data=a4, aes(x = a4$V2,y = a4$purity,fill=a4$V2))+geom_boxplot(alpha=0.4,outlier.colour = NA)+
    scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
  p <- p+theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Class") + ylab("purity") + ggtitle(paste("Immune Cell Types=",j,sep=""))+
    annotate("text",x=1.5, y=max(a4[,4]),colour="red",size=3,label=round(t.test(b1,b2)$p.value,3))+
    annotate("text",x=2, y=max(a4[,4]),colour="red",size=3,label=round(t.test(b1,b3)$p.value,3))+
    annotate("text",x=2.5, y=max(a4[,4]),colour="red",size=3,label=round(t.test(b1,b4)$p.value,3))+
    annotate("text",x=2.5, y=max(a4[,4])-0.1,colour="red",size=3,label=round(t.test(b2,b3)$p.value,3))+
    annotate("text",x=3, y=max(a4[,4])-0.1,colour="red",size=3,label=round(t.test(b2,b4)$p.value,3))+
    annotate("text",x=3.5, y=max(a4[,4])-0.2,colour="red",size=3,label=round(t.test(b3,b4)$p.value,3))
  ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/immun",j,".pdf",sep=""),width = 10,height=9)
}           


##对Allcancer26_class3_manhattan_ward.D_1075对每类做6种免疫细胞类型的purity纯度的boxplot整合为一图 
##Allcancer2_class4_manhattan_ward.D_1075
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1= read.csv("D:/tp53/immuneEstimation.txt",header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-",".",gsub("-01","",a1[,1]))
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
a=merge(f,a1,by.x="V1",by.y="barcode",all.x=T)
da=data.frame()
for(i in 3:ncol(a)){
  a2=a[,c(2,i)]
  colnames(a2)[2]="purity"
  a3=cbind(immuno=colnames(a)[i],a2)
  da=rbind(da,a3)
}
da1=na.omit(da)
da1[,2]=factor(da1[,2])
m=data.frame()
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==1),3]
  g2=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==2),3]
  g3=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==3),3]
  m1=cbind(i,paste(round(t.test(g1,g2)$p.value,4),round(t.test(g1,g3)$p.value,4),round(t.test(g2,g3)$p.value,4),sep="_"))
  m=rbind(m,m1)
}
da2=merge(da1,m,by.x="immuno",by.y="i",all.x=T)
p<-ggplot(data=da2, aes(x = da2[,1],y = da2[,3]))+geom_boxplot(aes(fill=da2[,2]),alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+coord_flip()
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Immuo")+ylab("purity")+ggtitle("Different Classes")
p <- p+annotate("text",x=da2[,1], y=1.2, colour="red",size=3,label=as.character(da2[,4]))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/6.pdf",sep=""),width = 10,height=9)  



##对Allcancer26_class3_manhattan_ward.D_1075对每类做其他免疫特征ImmunoSignature的分析
##Allcancer2_class4_manhattan_ward.D_1075
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f[,1]=as.character(f[,1])
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/ImmunoSignature.txt",sep=""),header=T,sep="\t")
f1[,1]=gsub("-",".",as.character(f1[,1]))
f2=f1[,!(colnames(f1) %in% c("OS","OS.Time","PFI","PFI.Time"))]
a=merge(f,f2,by.x="V1",by.y="TCGA.Participant.Barcode",all.x=T)
da=data.frame()
for(i in 6:ncol(a)){
  a2=a[,c(2,i)]
  colnames(a2)[2]="purity"
  a3=cbind(immuno=colnames(a)[i],a2)
  da=rbind(da,a3)
}
da1=na.omit(da)
da1[,2]=factor(da1[,2])
m=data.frame()
for(i in unique(as.character(da1[,1]))){
  g1=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==1),3]
  g2=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==2),3]
  g3=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==3),3]
  m1=cbind(i,paste(round(t.test(g1,g2)$p.value,4),round(t.test(g1,g3)$p.value,4),round(t.test(g2,g3)$p.value,4),sep="_"))
  m=rbind(m,m1)
}
da2=merge(da1,m,by.x="immuno",by.y="i",all.x=T)
p<-ggplot(data=da2, aes(x = da2[,1],y = da2[,3]))+geom_boxplot(aes(fill=da2[,2]),alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+coord_flip()
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Immuo")+ylab("purity")+ggtitle("Different Classes")
p <- p+annotate("text",x=da2[,1], y=5700, colour="red",size=2.5,label=as.character(da2[,4]))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/ImmunoSignature.pdf",sep=""),width = 10,height=9)  


##Allcancer2_class4_manhattan_ward.D_1075对每类做6种免疫细胞类型的purity纯度+其他免疫特征ImmunoSignature的整合为一图 
##Immunointegrate_p4是亚型一和非亚型一的p值，##Immunointegrate_puritymedian肿瘤纯度的中值也是亚型一和非亚型一
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1= read.csv("D:/tp53/immuneEstimation.txt",header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-",".",gsub("-01","",a1[,1]))
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f[,1]=as.character(f[,1])
a=merge(f,a1,by.x="V1",by.y="barcode",all.x=T)
da=data.frame()
for(i in 3:ncol(a)){
  a2=a[,c(2,i)]
  colnames(a2)[2]="purity"
  a3=cbind(immuno=colnames(a)[i],a2)
  da=rbind(da,a3)
}
da1=na.omit(da)
da1[,2]=factor(da1[,2])
m=data.frame()
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==1),3]
  g2=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==2),3]
  g3=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==3),3]
  g4=da1[(as.character(da1[,1])==i)&(as.character(da1[,2])==4),3]
  m1=cbind(i,signif(t.test(g1,g2)$p.value,2),signif(t.test(g1,g3)$p.value,2),signif(t.test(g1,g4)$p.value,2),signif(t.test(g2,g3)$p.value,2),signif(t.test(g2,g4)$p.value,2),signif(t.test(g3,g4)$p.value,2))
  m=rbind(m,m1)
  #write.table(m1,file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  m2=cbind(i,signif(t.test(g1,c(g2,g3,g4))$p.value,2),signif(t.test(g2,c(g1,g3,g4))$p.value,2),signif(t.test(g3,c(g1,g2,g4))$p.value,2),signif(t.test(g4,c(g1,g2,g3))$p.value,2))
  #write.table(m2,file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_p4.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  if(median(g1)>=median(c(g2,g3,g4))){
    h1=1
  }else{
    h1=-1
  }
  if(median(g2)>=median(c(g1,g3,g4))){
    h2=1
  }else{
    h2=-1
  }
  if(median(g3)>=median(c(g1,g2,g4))){
    h3=1
  }else{
    h3=-1
  }
  if(median(g4)>=median(c(g1,g2,g3))){
    h4=1
  }else{
    h4=-1
  }
  #write.table(paste(i,h1,h2,h3,h4,sep="\t"),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedian.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
da2=merge(da1,m,by.x="immuno",by.y="i",all.x=T)

f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/ImmunoSignature.txt",sep=""),header=T,sep="\t")
f2[,1]=gsub("-",".",as.character(f2[,1]))
f3=f2[,!(colnames(f2) %in% c("OS","OS.Time","PFI","PFI.Time"))]
a3=merge(f,f3,by.x="V1",by.y="TCGA.Participant.Barcode",all.x=T)
da3=data.frame()
for(i in 6:ncol(a3)){
  a4=a3[,c(2,i)]
  colnames(a4)[2]="purity"
  a5=cbind(immuno=colnames(a3)[i],a4)
  da3=rbind(da3,a5)
}
da4=na.omit(da3)
da4[,2]=factor(da4[,2])
mm=data.frame()
for(i in unique(as.character(da4[,1]))){
  g1=da4[(as.character(da4[,1])==i)&(as.character(da4[,2])==1),3]
  g2=da4[(as.character(da4[,1])==i)&(as.character(da4[,2])==2),3]
  g3=da4[(as.character(da4[,1])==i)&(as.character(da4[,2])==3),3]
  g4=da4[(as.character(da4[,1])==i)&(as.character(da4[,2])==4),3]
  m1=cbind(i,signif(t.test(g1,g2)$p.value,2),signif(t.test(g1,g3)$p.value,2),signif(t.test(g1,g4)$p.value,2),signif(t.test(g2,g3)$p.value,2),signif(t.test(g2,g4)$p.value,2),signif(t.test(g3,g4)$p.value,2))
  mm=rbind(mm,m1)
  #write.table(m1,file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  m2=cbind(i,signif(t.test(g1,c(g2,g3,g4))$p.value,2),signif(t.test(g2,c(g1,g3,g4))$p.value,2),signif(t.test(g3,c(g1,g2,g4))$p.value,2),signif(t.test(g4,c(g1,g2,g3))$p.value,2))
  #write.table(m2,file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_p4.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  if(median(g1)>=median(c(g2,g3,g4))){
    h1=1
  }else{
    h1=-1
  }
  if(median(g2)>=median(c(g1,g3,g4))){
    h2=1
  }else{
    h2=-1
  }
  if(median(g3)>=median(c(g1,g2,g4))){
    h3=1
  }else{
    h3=-1
  }
  if(median(g4)>=median(c(g1,g2,g3))){
    h4=1
  }else{
    h4=-1
  }
  #write.table(paste(i,h1,h2,h3,h4,sep="\t"),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedian.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
da5=merge(da4,mm,by.x="immuno",by.y="i",all.x=T)
da6=rbind(da2,da5)
p<-ggplot(data=da6, aes(x = da6[,1],y = da6[,3],color=da6[,2]))+geom_boxplot(aes(fill=da6[,2]),lwd=0.1,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+
  scale_color_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+coord_flip()
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Immuo")+ylab("purity")+ggtitle("Different Classes")+ylim(-2500,2500)
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate.pdf",sep=""),width = 10,height=9)  



#write.table(paste(i,j,da7,sep="\t"),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedian.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


##根据肿瘤纯度和p值画图
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_p4.txt",sep=""),header=F,sep="\t")
for(i in 1:nrow(f1)){
  #write.table(paste(f1[i,1],1:4,f1[i,2:5],sep="\t"),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_pvalue.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedian.txt",sep=""),header=F,sep="\t")
for(i in 1:nrow(f)){
  #write.table(paste(f[i,1],1:4,f[i,2:5],sep="\t"),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedianup.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}


f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_puritymedianup.txt",sep=""),header=F,sep="\t")
colnames(f)[3]="purity"
f[,1]=as.character(f[,1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_pvalue.txt",sep=""),header=F,sep="\t")
colnames(f2)[3]="pvalue"
f2[,1]=as.character(f2[,1])
da1=merge(f2,f,by.x=c("V1","V2"),by.y=c("V1","V2"))
da1[which(da1[,3]>=0.05),3]=1
da1[,3]=abs(log10(da1[,3]))
ma=max(da1[,3])
mi=min(da1[,3])
da1[,3]=(da1[,3]-mi)/(ma-mi)
da1[,3]=da1[,3]+0.1
p <- ggplot(da1,aes(x=factor(1),y=da1[,3],fill=purity))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da1[,1]~da1[,2])+
  scale_fill_gradient(low="#12B826",high = "#FF331C")  #"#12B826","#FF331C"
p <- p+theme(legend.position="right",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=6,colour="black"))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank())
p <- p + xlab("ImmunoSignature") + ylab("Subtype") + ggtitle("Immuno Analysis")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Immunointegrate_heatmap2.pdf",sep=""),width = 4,height=12)









##对Allcancer26_class3_manhattan_ward.D_1075对每类做突变负荷的分析，样本对应的突变基因的数目 
##Allcancer2_class4_manhattan_ward.D_1075
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f[,1]=as.character(f[,1])
f1=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/3891.txt",header=T,sep="\t")
f1[,1]=gsub("-",".",as.character(f1[,1]))
f2=f1[,c("submitter_id","total.number.of.mutation.gene.of.the.sample")]
a=na.omit(unique(merge(f,f2,by.x="V1",by.y="submitter_id",all.x=T)))
a[,2]=factor(a[,2])
g1=a[as.numeric(as.character(a[,2]))==1,3]
g2=a[as.numeric(as.character(a[,2]))==2,3]
g3=a[as.numeric(as.character(a[,2]))==3,3]
g4=a[as.numeric(as.character(a[,2]))==4,3]
#signif(t.test(g3,g4)$p.value,2)
p<-ggplot(data=a, aes(x = a[,2],y = a[,3]))+geom_boxplot(aes(fill=a[,2]),alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4")) 
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Class")+ylab("Mutation load")+ggtitle("Different Classes")
p <- p+annotate("text",x=1.5, y=max(a[,3]), colour="red",size=2.5,label=round(t.test(g1,g2)$p.value,4))+
  annotate("text",x=2, y=max(a[,3]), colour="red",size=2.5,label=round(t.test(g1,g3)$p.value,4))+
  annotate("text",x=2.5, y=max(a[,3])-500, colour="red",size=2.5,label=round(t.test(g2,g3)$p.value,4))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Mutationload.pdf",sep=""),width = 10,height=9)  

p<-ggplot(data=a, aes(x = a[,2],y = a[,3]))+geom_violin(aes(fill=a[,2]),alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4")) 
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Class")+ylab("Mutation load")+ggtitle("Different Classes")
p <- p+annotate("text",x=1.5, y=max(a[,3]), colour="red",size=2.5,label=round(t.test(g1,g2)$p.value,4))+
  annotate("text",x=2, y=max(a[,3]), colour="red",size=2.5,label=round(t.test(g1,g3)$p.value,4))+
  annotate("text",x=2.5, y=max(a[,3])-500, colour="red",size=2.5,label=round(t.test(g2,g3)$p.value,4))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Mutationload1.pdf",sep=""),width = 10,height=9)  

p<-ggplot(data=a, aes(x = a[,2],y = a[,3]))+geom_violin(aes(fill=a[,2]),alpha=0.4,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4")) 
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Class")+ylab("Mutation load")+ggtitle("Different Classes")+ylim(0,1000)
p <- p+annotate("text",x=1, y=median(g1)+20, colour="black",size=2.5,label=median(g1))+
  annotate("text",x=2, y=median(g2)+20, colour="black",size=2.5,label=median(g2))+
  annotate("text",x=3, y=median(g3)+20, colour="black",size=2.5,label=median(g3))+
  annotate("text",x=4, y=median(g4)+20, colour="black",size=2.5,label=median(g4))+
  annotate("text",x=1, y=median(g1), colour="black",size=2.5,labe l="____")+
  annotate("text",x=2, y=median(g2), colour="black",size=2.5,label="____")+
  annotate("text",x=3, y=median(g3), colour="black",size=2.5,label="____")+
  annotate("text",x=4, y=median(g4), colour="black",size=2.5,label="____")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/immuno/Mutationload2.pdf",sep=""),width = 10,height=9)  


##对Allcancer26_class3_manhattan_ward.D_1075对每类做突变类型的堆叠柱状图，样本对应的突变类型5UTR、3UTR、A>T
##Allcancer2_class4_manhattan_ward.D_1075
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f[,1]=as.character(f[,1])
f1=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/3891.txt",header=T,sep="\t")
f1[,1]=gsub("-",".",as.character(f1[,1]))
f2=data.frame(na.omit(cbind(as.character(f1[,1]),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f1$HGVSc),4)))))
colnames(f2)=c("submitter_id","Mutype")
f2[,1]=as.character(f2[,1])
f2[,2]=gsub("\\+","",as.character(f2[,2]))
f3=data.frame(na.omit(cbind(as.character(f1[,1]),as.character(f1$Variant_Classification))))
colnames(f3)=c("submitter_id","Variant_Classification")
f3[,1]=as.character(f3[,1])
f3[,2]=as.character(f3[,2])
f4=data.frame(na.omit(cbind(as.character(f1[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f1$HGVSp_Short),4))))))
colnames(f4)=c("submitter_id","Muposition")
f4[,1]=as.character(f4[,1])
f4[,2]=as.numeric(as.character(f4[,2]))

a1=na.omit(unique(merge(f,f2,by.x="V1",by.y="submitter_id",all.x=T)))  #A>T a1=a1[nchar(a1[,3])==3,]
a1[,2]=factor(a1[,2])
a2=data.frame(table(a1[,2:3]))
a3=data.frame(table(a1[,2]))
a4=merge(a2,a3,by.x="V2",by.y="Var1",all.x=T)
a4[,5]=a4[,3]/a4[,4]
p<-ggplot(data=a4, aes(x = a4[,1],y = a4[,5],fill=Mutype))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of Mutype in different classes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Mutype.pdf",sep=""),width = 10,height=9)  

b1=na.omit(unique(merge(f,f3,by.x="V1",by.y="submitter_id",all.x=T)))  #Missense_Mutation
b1[,2]=factor(b1[,2])
b1[,3]=gsub("^3UTR.*","3UTR",as.character(b1[,3]))
b1[,3]=gsub("^5UTR.*","5UTR",as.character(b1[,3]))
b2=data.frame(table(b1[,2:3]))
b3=data.frame(table(b1[,2]))
b4=merge(b2,b3,by.x="V2",by.y="Var1",all.x=T)
b4[,5]=b4[,3]/b4[,4]
p<-ggplot(data=b4, aes(x = b4[,1],y = b4[,5],fill=Variant_Classification))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of Variant_Classification in different classes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Variant_Classification.pdf",sep=""),width = 10,height=9)  

d1=na.omit(unique(merge(f,f4,by.x="V1",by.y="submitter_id",all.x=T)))  #Muposition
d1[,2]=factor(d1[,2])
d2=data.frame(table(d1[,2:3]))
d3=data.frame(table(d1[,2]))
d4=merge(d2,d3,by.x="V2",by.y="Var1",all.x=T)
d4[,5]=d4[,3]/d4[,4]
p<-ggplot(data=d4, aes(x = d4[,1],y = d4[,5],fill=Muposition))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of Muposition in different classes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Muposition.pdf",sep=""),width = 10,height=9)  

d1=na.omit(unique(merge(f,f4,by.x="V1",by.y="submitter_id",all.x=T)))  #Muposition hotspot
d1[,4]="NO"
d1[(d1[,3]<=40),4]="TAD1"
d1[(d1[,3]>40 & d1[,3]<=60),4]="TAD2"
d1[(d1[,3]>60 & d1[,3]<=95),4]="PRD"
d1[(d1[,3]>= 100 & d1[,3]<=300),4]="DNAbinding"
d1[(d1[,3]>=325 & d1[,3]<=356),4]="Tet"
d1[(d1[,3]>356 & d1[,3]<=393),4]="Basic"
d1[,2]=factor(d1[,2])
d2=data.frame(table(d1[,c(2,4)]))
d3=data.frame(table(d1[,2]))
d4=merge(d2,d3,by.x="V2",by.y="Var1",all.x=T)
d4[,5]=d4[,3]/d4[,4]
d4[,2]=factor(d4[,2],c("NO","TAD1","TAD2","PRD","DNAbinding","Tet","Basic"))
p<-ggplot(data=d4, aes(x = d4[,1],y = d4[,5],fill=V4))+geom_bar(stat="identity",alpha=0.5)
p <- p+ theme(legend.position="top")+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p + xlab("Class") + ylab("Sample Frequency") + ggtitle("The Sample Frequency of Muhotspot in different classes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Muhotspot.pdf",sep=""),width = 10,height=9)  

d1=na.omit(unique(merge(f,f4,by.x="V1",by.y="submitter_id",all.x=T)))  #Muposition hotspot1 slide
d2=data.frame(table(d1[,2:3]))
d2[,2]=as.numeric(as.character(d2[,2]))
d2[,4]=paste(d2[,1],d2[,2],sep="_")
df=data.frame(cbind(class=rep(c(1,2,3),each=400),Var1=rep(1:400,3)))
#df=data.frame(cbind(class=rep(c(1,2,3,4),each=400),Var1=rep(1:400,4)))
df[,3]=paste(df[,1],df[,2],sep="_")
d3=merge(df,d2,by.x="V3",by.y="V4",all.x=T)
d3[,2]=factor(d3[,2])
p<-ggplot(d3, aes(x = d3[,3], y = d3[,6], fill = d3[,2]))+geom_point(aes(shape=d3[,2]),colour="red") + geom_line()+ facet_grid(d3[,2] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+geom_rect(xmin=1,xmax=40,ymin=-5,ymax=-1,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-5,ymax=-1,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-5,ymax=-1,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-5,ymax=-1,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-5,ymax=-1,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-5,ymax=-1,fill="green",alpha=0.5)
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Muhotspot1.pdf",sep=""),width = 10,height=9)  

 
d1=na.omit(unique(merge(f,f4,by.x="V1",by.y="submitter_id",all.x=T)))  #Muposition hotspot2  
d2=data.frame(table(d1[,2:3]))
d2[,2]=as.numeric(as.character(d2[,2]))
d2[,4]=paste(d2[,1],d2[,2],sep="_")
df=data.frame(cbind(class=rep(c(1,2,3),each=400),Var1=rep(1:400,3)))
#df=data.frame(cbind(class=rep(c(1,2,3,4),each=400),Var1=rep(1:400,4)))
df[,3]=paste(df[,1],df[,2],sep="_")
d3=merge(df,d2,by.x="V3",by.y="V4",all.x=T)
d3[,2]=factor(d3[,2])
p<-ggplot(d3, aes(x = d3[,3], y = d3[,6],fill=d3[,2]))+geom_bar(width=0.1,alpha=0.8,stat="identity", position="dodge")+ facet_grid(d3[,2] ~ .)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+geom_rect(xmin=1,xmax=40,ymin=-5,ymax=-1,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-5,ymax=-1,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-5,ymax=-1,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-5,ymax=-1,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-5,ymax=-1,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-5,ymax=-1,fill="green",alpha=0.5)
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Muhotspot2.pdf",sep=""),width = 10,height=9)  

d1=na.omit(unique(merge(f,f4,by.x="V1",by.y="submitter_id",all.x=T)))  #Muposition hotspot3  
d2=data.frame(table(d1[,2:3]))
d2[,2]=as.numeric(as.character(d2[,2]))
d2[,4]=paste(d2[,1],d2[,2],sep="_")
df=data.frame(cbind(class=rep(c(1,2,3),each=400),Var1=rep(1:400,3)))
#df=data.frame(cbind(class=rep(c(1,2,3,4),each=400),Var1=rep(1:400,4)))
df[,3]=paste(df[,1],df[,2],sep="_")
d3=merge(df,d2,by.x="V3",by.y="V4",all.x=T)
d3[,2]=factor(d3[,2])
d4=data.frame(table(d1[,2]))
d5=merge(d3,d4,by.x="class",by.y="Var1",all.x=T)
d5[,8]=(d5[,6]/d5[,7])*100
p<-ggplot(d5, aes(x = d5[,3], y = d5[,8],fill=d5[,1]))+geom_bar(width=0.1,alpha=0.8,stat="identity", position="dodge")+ facet_grid(d5[,1] ~ .)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+geom_rect(xmin=1,xmax=40,ymin=-5,ymax=-0.2,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-5,ymax=-0.2,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-5,ymax=-0.2,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-5,ymax=-0.2,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-5,ymax=-0.2,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-5,ymax=-0.2,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-5,ymax=-0.2,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-5,ymax=-0.2,fill="green",alpha=0.5)
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/immuno/Muhotspot3.pdf",sep=""),width = 10,height=9)  



          

##Allcancer26_class3_manhattan_ward.D_1075做热图 筛基因筛样本
##Allcancer2_class4_manhattan_ward.D_1075
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
f3=f2[order(f2$V2,decreasing=F),]
h3=h2[,f3[,1]]
y1<-as.matrix(h3)
row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
col_anno = data.frame(Class=as.character(f2[,2]), row.names=as.character(f2[,1]))
#ann_colors = list(Time = c("white", "firebrick"), #连续数值型分组可设置成渐变
                  #CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                  #GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
ann_colors = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"))
pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,annotation_col = col_anno,fontsize_row=5,fontsize_col=3,
         annotation_colors = ann_colors,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075.pdf",sep=""),width=8,height=8,main="Allcancer")  

#write.table(rownames(h3),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_gene.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
#write.table(h3,file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_manhattan_ward.D_1075.txt",sep=""), row.names = T,col.names=T, quote = F,sep = "\t",append=TRUE)


##Allcancer26_class3_manhattan_ward.D_1075做热图 筛基因筛样本 把基因整合为功能，得到样本对功能的不同模式。筛功能
##Allcancer2_class4_manhattan_ward.D_1075 把样本分类并标出相应的癌症
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])
f3=f22[order(f22$V2,decreasing=F),]
h3=h2[,f3[,1]]
h3$name=rownames(h3)
h4=merge(f,h3,by.x="gene_Name",by.y="name",all.y=T)
d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(h4[,2])){
  m=h4[h4$Function==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(h4)[-1]

fu1<-function(x) sum(x==0)  
a3=data.frame(apply(d11[,2:ncol(d11)],1,fu1))
a4=rownames(a3)[which(a3[,1]!=1075)] 
d12=na.omit(d11[a4,-1])

y1<-as.matrix(d11[,2:ncol(d11)])
#row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
col_anno = data.frame(Class=as.character(f22[,2]), row.names=as.character(f22[,1]))
ann_colors = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"))
pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,annotation_col = col_anno,fontsize_row=6,fontsize_col=3,
         annotation_colors = ann_colors,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function.pdf",sep=""),width=8,height=8,main="Allcancer")  

y2<-as.matrix(d12[,2:ncol(d12)]) 
col_anno = data.frame(Class=as.character(f22[,2]), row.names=as.character(f22[,1]))
ann_colors = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"))
pheatmap(y2,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,annotation_col = col_anno,fontsize_row=6,fontsize_col=3,
         annotation_colors = ann_colors,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function1.pdf",sep=""),width=8,height=8,main="Allcancer")  

y2<-as.matrix(d12[,2:ncol(d12)]) 
col_anno2 = data.frame(Class=as.character(f22[,2]),Cancer=as.character(f22[,3]),row.names=as.character(f22[,1]))
ann_colors2 = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"),
                   Cancer = c(GBM ="deeppink3",
                              COAD ="goldenrod3",READ ="goldenrod3",
                              KIRC ="royalblue4",KIRP ="royalblue4",
                              LUSC ="palegreen4",LUAD ="palegreen4",
                              UCEC ="black",
                              ESCA ="blueviolet",BLCA ="brown3",
                              PAAD ="chocolate2",LIHC ="coral4",SARC ="cornflowerblue",
                              THYM ="bisque4",SKCM ="chartreuse3",BRCA ="cadetblue3",HNSC ="azure",
                              STAD ="darkslateblue",CESC ="bisque2",PRAD ="lightcoral"))
pheatmap(y2,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,annotation_col = col_anno2,fontsize_row=6,fontsize_col=3,
         annotation_colors = ann_colors2,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function2.pdf",sep=""),width=8,height=8,main="Allcancer")  

 
 


##Allcancer26_class3_manhattan_ward.D_1075做热图 筛基因筛样本 把基因整合为功能，得到样本对功能的不同模式。筛功能
##Allcancer2_class4_manhattan_ward.D_1075 把每一亚型下的所有样本整合成一列，比例max，最终结果为4亚型x功能或3亚型x功能的热图。
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
f3=f2[order(f2$V2,decreasing=F),]
h3=h2[,f3[,1]]
h3$name=rownames(h3)
h4=merge(f,h3,by.x="gene_Name",by.y="name",all.y=T)
d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(h4[,2])){
  m=h4[h4$Function==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(h4)[-1]

fu1<-function(x) sum(x==0)  
a3=data.frame(apply(d11[,2:ncol(d11)],1,fu1))
a4=rownames(a3)[which(a3[,1]!=1075)] 
d12=na.omit(d11[a4,-1])
d13=data.frame(t(d12))
d13$Sample=rownames(d13)
d14=merge(f2,d13,by.x="V1",by.y="Sample",all.x=T)
d15=data.frame(Class=unique(f2[,2]))
rownames(d15)<-d15[,1]
for(i in unique(d14[,2])){
  m=d14[d14$V2==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d15[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d15[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d15[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d15)=colnames(d14)[-1]
d16=d15[order(d15$V2,decreasing=F),]
y1<-as.matrix(t(d16[,2:ncol(d16)]))
pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,cellwidth = 30,fontsize_row=6,fontsize_col=8,
         file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function3.pdf",sep=""),width=6,height=6,main="Allcancer")  



##Allcancer26_class3_manhattan_ward.D_1075做热图 筛基因筛样本 把基因整合为功能，得到样本对功能的不同模式。筛功能
##Allcancer2_class4_manhattan_ward.D_1075 把每一亚型下的所有样本整合成一列，mean，最终结果为4亚型x功能或3亚型x功能的热图。
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
f3=f2[order(f2$V2,decreasing=F),]
h3=h2[,f3[,1]]
h3$name=rownames(h3)
h4=merge(f,h3,by.x="gene_Name",by.y="name",all.y=T)
d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(h4[,2])){
  m=h4[h4$Function==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(h4)[-1]

fu1<-function(x) sum(x==0)  
a3=data.frame(apply(d11[,2:ncol(d11)],1,fu1))
a4=rownames(a3)[which(a3[,1]!=1075)] 
d12=na.omit(d11[a4,-1])
d13=data.frame(t(d12))
d13$Sample=rownames(d13)
d14=merge(f2,d13,by.x="V1",by.y="Sample",all.x=T)
d15=data.frame(Class=unique(f2[,2]))
rownames(d15)<-d15[,1]
for(i in unique(d14[,2])){
  m=d14[d14$V2==i,]
  for(j in 3:ncol(m)){
    d15[i,j-1]=mean(m[,j])
  }
}
colnames(d15)=colnames(d14)[-1]
d16=d15[order(d15$V2,decreasing=F),]
y1<-as.matrix(t(d16[,2:ncol(d16)]))
pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,cellwidth = 30,fontsize_row=6,fontsize_col=8,
         file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function3mean.pdf",sep=""),width=6,height=6,main="Allcancer")  

pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,cellheight=20,cellwidth = 30,fontsize_row=6,fontsize_col=8,
         col=colorsChoice(9),legend_breaks=c(-0.3,-0.2,-0.1,-0.01,-0.005,0,0.01,0.04,0.07,0.1),breaks=c(-0.3,-0.2,-0.1,-0.01,-0.005,0,0.01,0.04,0.07,0.1),
         border_color=NA,display_numbers = TRUE,number_color="black",fontsize_number=6,number_format='%.2f',
         file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function3mean_2.pdf",sep=""),width=6,height=6,main="Pancancer")  



##Allcancer26_class3_manhattan_ward.D_1075做热图 筛基因筛样本 把基因整合为功能，得到样本对功能的不同模式。筛功能
##Allcancer2_class4_manhattan_ward.D_1075 把样本分类并标出相应的癌症,进行单个癌症的分析
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=26
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])
f3=f22[order(f22$V2,decreasing=F),]
h3=h2[,f3[,1]]
h3$name=rownames(h3)
h4=merge(f,h3,by.x="gene_Name",by.y="name",all.y=T)
d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(h4[,2])){
  m=h4[h4$Function==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(h4)[-1]

fu1<-function(x) sum(x==0)  
a3=data.frame(apply(d11[,2:ncol(d11)],1,fu1))
a4=rownames(a3)[which(a3[,1]!=1075)] 
d12=na.omit(d11[a4,-1])
y2<-as.matrix(d12[,2:ncol(d12)]) 
col_anno2 = data.frame(Class=as.character(f22[,2]),Cancer=as.character(f22[,3]),row.names=as.character(f22[,1]))
ann_colors2 = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"),
                   Cancer = c(GBM ="azure",
                              COAD ="deeppink3",READ ="goldenrod3",
                              KIRC ="azure",KIRP ="azure",
                              LUSC ="azure",LUAD ="azure",
                              UCEC ="azure",
                              ESCA ="azure",BLCA ="azure",
                              PAAD ="azure",LIHC ="azure",SARC ="azure",
                              THYM ="azure",SKCM ="azure",BRCA ="azure",HNSC ="azure",
                              STAD ="azure",CESC ="azure",PRAD ="azure"))
pheatmap(y2,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,annotation_col = col_anno2,fontsize_row=6,fontsize_col=3,
         annotation_colors = ann_colors2,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/D/Allcancer1075_function_COADREAD.pdf",sep=""),width=8,height=8,main="Allcancer")  



##Allcancer26_class3_manhattan_ward.D_1075 进行单个癌症的分析##SKCM\PAAD\COADREAD\GBM  ##ESCA\HNSC  ##LUSCLUAD\BRCA
##Allcancer2_class4_manhattan_ward.D_1075  
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])
f3=f22[which(f22$project_id=="GBM"),]  ##SKCM\PAAD\COADREAD\GBM  ##ESCA\HNSC  ##LUSCLUAD\BRCA
f3=f3[order(f3$V2,decreasing=F),]
h3=h2[,f3[,1]]
h3$name=rownames(h3)
h4=merge(f,h3,by.x="gene_Name",by.y="name",all.y=T)
d11=data.frame(Function=unique(f[,2]))
rownames(d11)<-d11[,1]
for(i in unique(h4[,2])){
  m=h4[h4$Function==i,]
  for(j in 3:ncol(m)){
    m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j-1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
    }
  }
}
colnames(d11)=colnames(h4)[-1]

fu1<-function(x) sum(x==0)  
a3=data.frame(apply(d11[,2:ncol(d11)],1,fu1))
a4=rownames(a3)[which(a3[,1]!=1075)] 
d12=na.omit(d11[a4,-1])
y2<-as.matrix(d12[,2:ncol(d12)]) 
col_anno2 = data.frame(Class=as.character(f3[,2]),Cancer=as.character(f3[,3]),row.names=as.character(f3[,1]))
ann_colors2 = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"),
                   Cancer = c(GBM ="deeppink3"))
pheatmap(y2,cluster_row = FALSE,cluster_col = FALSE,color =colorsChoice(256),cellheight=20,annotation_col = col_anno2,fontsize_row=6,fontsize_col=3,
         annotation_colors = ann_colors2,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/Allcancer1075_function_GBM1.pdf",sep=""),width=8,height=8,main="Allcancer")  




##Allcancer2_class4_manhattan_ward.D_1075  进行单个癌症的分析##SKCM\PAAD\COAD\GBM  ##ESCA\HNSC  ##LUSC\LUAD\BRCA
##根据样本对gene，最终形成cancer对gene、cancer对功能的热图 
##"SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA"9个 比例max
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])
##比例max
for(cancer in c("SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA")){
  f3=f22[which(f22$project_id==cancer),]  ##SKCM\PAAD\COAD\GBM  ##ESCA\HNSC  ##LUSC\LUAD\BRCA
  f3=f3[order(f3$V2,decreasing=F),]
  h3=data.frame(t(h2[,f3[,1]]))
  h3$name=rownames(h3)   ##样本对gene
  h4=merge(f3,h3,by.x="V1",by.y="name",all.x=T)
  d1=data.frame(Class=unique(h4[,2]))
  rownames(d1)<-d1[,1]
  for(i in unique(h4[,2])){
    m=h4[h4$V2==i,]
    for(j in 4:ncol(m)){   ##比例max
      m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
      if(which(m1[1,]==max(m1)) ==3){
        d1[i,j-2]=0
      }else if(which(m1[1,]==max(m1)) ==2){
        d1[i,j-2]=-length(which(m[,j]<0))/length((m[,j]))
      }else if(which(m1[1,]==max(m1)) ==1){
        d1[i,j-2]=length(which(m[,j]>0))/length((m[,j]))
      }
    }
  }
  colnames(d1)=colnames(h4)[c(-1,-3)]
  write.table(t(d1),file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/max/",cancer,"_Gene.txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)  
  d2=data.frame(t(d1)) 
  d2$name=rownames(d2)
  d3=d2[-1,]
  d4=merge(f,d3,by.x="gene_Name",by.y="name",all.y=T)   
  d11=data.frame(Function=unique(d4[,2]))
  rownames(d11)<-d11[,1]
  for(i in unique(d4[,2])){
    m=d4[d4$Function==i,]
    for(j in 3:ncol(m)){  ##比例max
      m[,j]=as.numeric(as.character(m[,j]))
      m1=data.frame(up=length(which(m[,j]>0)),down=length(which(m[,j]<0)),normal=length(which(m[,j]==0)))
      if(which(m1[1,]==max(m1)) ==3){
        d11[i,j-1]=0
      }else if(which(m1[1,]==max(m1)) ==2){
        d11[i,j-1]=-length(which(m[,j]<0))/length((m[,j]))
      }else if(which(m1[1,]==max(m1)) ==1){
        d11[i,j-1]=length(which(m[,j]>0))/length((m[,j]))
      }
    }
  }
  colnames(d11)=colnames(d4)[-1]
  write.table(d11,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/max/",cancer,"_Function.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)  
}

##mean
for(cancer in c("SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA")){
  f3=f22[which(f22$project_id==cancer),]   
  f3=f3[order(f3$V2,decreasing=F),]
  h3=data.frame(t(h2[,f3[,1]]))
  h3$name=rownames(h3)   ##样本对gene
  h4=merge(f3,h3,by.x="V1",by.y="name",all.x=T)
  d1=data.frame(Class=unique(h4[,2]))
  rownames(d1)<-d1[,1]
  for(i in unique(h4[,2])){
    m=h4[h4$V2==i,]
    for(j in 4:ncol(m)){  ##mean
      d1[i,j-2]=mean(as.numeric(as.character(m[,j]))) 
    }
  }
  colnames(d1)=colnames(h4)[c(-1,-3)]
  write.table(t(d1),file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/mean/",cancer,"_Gene.txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)  
  d2=data.frame(t(d1)) 
  d2$name=rownames(d2)
  d3=d2[-1,]
  d4=merge(f,d3,by.x="gene_Name",by.y="name",all.y=T)   
  d11=data.frame(Function=unique(d4[,2]))
  rownames(d11)<-d11[,1]
  for(i in unique(d4[,2])){
    m=d4[d4$Function==i,]
    for(j in 3:ncol(m)){  ##mean
      d11[i,j-1]=mean(as.numeric(as.character(m[,j])))
    }
  }
  colnames(d11)=colnames(d4)[-1]
  write.table(d11,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/mean/",cancer,"_Function.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)  
}


################################
##网络图构建单个癌症
##对"SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA"9个，
##点为每个基因，大小为在各亚型中的值，颜色为对应的功能，
##正方形为上调，圆形为下调，边为2个基因间的关系，<->为促进，
##->为抑制，-为NA或修改（不确定的）。
#install.packages("igraph")
library(igraph) 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in c("SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA")){
  f1<-read.csv(paste(path,"mean/",i,"_Gene.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$V2
  for(k in 2:ncol(f1)){
    f1$expressmean=f1[,k]
    f11=f1[intersect(as.character(f[,1]),as.character(f1$V2)),c("V2","expressmean")]
    m=merge(f11,f,by.x="V2",by.y="gene_Name")
    fu_colour=data.frame(Function=unique(m[,3]),
                         colour=c("orangered1","yellow3","snow3","palegoldenrod","maroon1","dodgerblue4","plum",
                                  "mediumvioletred","royalblue3","yellow","steelblue3","slateblue2","deepskyblue",
                                  "lightseagreen","rosybrown2","sienna1","tomato3","purple4","red","mediumpurple1"))
    m1=m[,c("V2","Function","expressmean")]
    #m1[which(m1[,3]>=0),4]="square"
    #m1[which(m1[,3]<0),4]="circle"
    m1[,4]="circle"
    colnames(m1)[4]="shape"
    #m1[which(m1[,3]>=0),5]="red"
    #m1[which(m1[,3]<0),5]="green"
    #colnames(m1)[5]="framecolor"
    m1=merge(m1,fu_colour)
    m1[,6]=as.character(m1$colour)
    m1[which(m1[,3]>0),6]="red4"
    m1[which(m1[,3]<0),6]="darkgreen"
    colnames(m1)[6]="framecolor"
    m1=m1[,c("V2","Function","expressmean","shape","colour","framecolor")]
    f2= read.csv(paste("F:/tp53/3generegulate/new_relationship1/",i,".txt",sep=""),header=T,sep="\t")
    f2[,2]=as.character((f2[,2]))
    f2[,3]=as.character((f2[,3]))
    f2[,4]=as.character((f2[,4]))
    n1=unique(f2)
    for(j in 1:nrow(n1)){
      if(n1[j,4]=="-"){
        a=n1[j,2] 
        n1[j,2]=n1[j,3]
        n1[j,3]=a
        n1[j,4]=paste(n1[j,4],0,sep="")
        n1[j,"count"]=0.5
        n1[j,"linetype"]="->"
        n1[j,"colour"]="blue"
      }
      else if(n1[j,4]=="+"){
        n1[j,"count"]=0.5
        n1[j,"linetype"]="<->"
        n1[j,"colour"]="red"
      }
      else{
        n1[j,"count"]=0.5
        n1[j,"linetype"]="-"
        n1[j,"colour"]="gray"
      }
    }
    n11=n1[,c("gene1","gene2","function.","count","linetype","colour")]
    data1=n11[(n11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))) & (n11[,2] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))),]
    data2=m1[m1[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1])),]
    
    g1 <- graph_from_data_frame(data1, directed=TRUE, vertices=data2)
    a=as.data.frame(degree(g1))
    a1=rownames(a)[which(a$`degree(g1)`!=0)]
    data3=data1[(data1[,1] %in% a1) & (data1[,2] %in% a1),]
    data4=data2[data2[,1] %in% a1,]
    g2 <- graph_from_data_frame(data3, directed=TRUE, vertices=data4)
    
    pdf(paste(path,"mean_graph0/",i,"_",colnames(f1)[k],".pdf",sep=""))  
    plot(g2, layout =layout.kamada.kawai ,vertex.label=V(g2)$V2,vertex.label.cex=0.3,
         vertex.size=3*(abs(V(g2)$expressmean)+0.5),vertex.label.dist=-0.4,vertex.label.color="black",vertex.color=V(g2)$colour,
         vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,vertex.frame.color=V(g2)$framecolor,
         edge.width=E(g2)$count/10,edge.color=E(g2)$colour,edge.arrow.size=0.05,main=paste(i,"_",colnames(f1)[k],sep=""))
    legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
    dev.off()
  }
}
#layout=layout.circle\layout.fruchterman.reingold\layout.kamada.kawai 
tkplot(g2, layout =layout.circle ,vertex.label=V(g2)$V2,vertex.label.cex=0.3,
       vertex.size=3*(abs(V(g2)$expressmean)+0.5),vertex.label.dist=-0.4,vertex.label.color="black",vertex.color=V(g2)$colour,
       vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,vertex.frame.color=V(g2)$framecolor,
       edge.width=E(g2)$count/10,edge.color=E(g2)$colour,edge.arrow.size=0.05,main=paste(i,"_",colnames(f1)[k],sep=""))


##############################
##单个癌症circlize图
library(circlize)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/D/SingleCancer/"
for(i in c("SKCM","PAAD","COAD","GBM","ESCA","HNSC","LUSC","LUAD","BRCA")){
  f3<-read.csv(paste(path,"mean/",i,"_Function.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  f1<-read.csv(paste(path,"mean/",i,"_Gene.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$V2
  for(k in 2:ncol(f1)){
    f1$expressmean=f1[,k]
    f11=f1[intersect(as.character(f[,1]),as.character(f1$V2)),c("V2","expressmean")]
    m=merge(f11,f,by.x="V2",by.y="gene_Name")
    fu_colour=data.frame(Function=unique(m[,3]),
                         colour=c("orangered1","yellow3","snow3","palegoldenrod","maroon1","dodgerblue4","plum",
                                  "mediumvioletred","royalblue3","yellow","steelblue3","slateblue2","deepskyblue",
                                  "lightseagreen","rosybrown2","sienna1","tomato3","purple4","red","mediumpurple1"))
    a2=merge(fu_colour,f3[,c(1,k)])
    rownames(a2)=a2$Function
    a3=a2[fu_colour[,1],]
    a3[a3[,3]>0,4]="red4"
    a3[a3[,3]==0,4]="seashell3"
    a3[a3[,3]<0,4]="darkgreen"
    a4=cbind(a3,start=1:nrow(a3),end=2:(nrow(a3)+1))
    a4$Function=as.character(a4$Function)
    a4$colour=as.character(a4$colour)
    
    circos.initializeWithIdeogram(plotType = NULL) #去原始圈
    circos.par("gap.degree"=0.01) #各板块间隔
    circos.genomicInitialize(a4[,c(1,5,6)])  
    circos.genomicTrackPlotRegion(a4[,c(1,5,6)],ylim=c(-1,1),track.height=0.1,bg.col=a4$colour,bg.border = NA)
    circos.genomicTrackPlotRegion(a4[,c(1,5,6)],ylim=c(-1,1),track.height=0.1,bg.col="white",bg.border ="gray") 
    circos.genomicPoints(a4[,c(5,6)], 0, cex = 1.5*(abs(a4[,3])+0.5), pch = 16, col =a4[,4])
    circos.clear() 
    #colnames(f1)[k]
  }
} 



########################################################
##Allcancer2_class4_manhattan_ward.D_1075  进行PanCancer癌症的分析
##根据样本对gene，最终形成PanCancerr对gene、PanCancer对功能的热图 mean
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])

number=2
fu<-function(x) sum(x!=0)  
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] 
h2=na.omit(h[a2,-1])
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])
##mean
f3=f22   
f3=f3[order(f3$V2,decreasing=F),]
h3=data.frame(t(h2[,f3[,1]]))
h3$name=rownames(h3)   ##样本对gene
h4=merge(f3,h3,by.x="V1",by.y="name",all.x=T)
d1=data.frame(Class=unique(h4[,2]))
rownames(d1)<-d1[,1]
for(i in unique(h4[,2])){
  m=h4[h4$V2==i,]
  for(j in 4:ncol(m)){  ##mean
    d1[i,j-2]=mean(as.numeric(as.character(m[,j]))) 
  }
}
colnames(d1)=colnames(h4)[c(-1,-3)]
write.table(t(d1),file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/mean/All_Gene.txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)  
d2=data.frame(t(d1)) 
d2$name=rownames(d2)
d3=d2[-1,]
d4=merge(f,d3,by.x="gene_Name",by.y="name",all.y=T)   
d11=data.frame(Function=unique(d4[,2]))
rownames(d11)<-d11[,1]
for(i in unique(d4[,2])){
  m=d4[d4$Function==i,]
  for(j in 3:ncol(m)){  ##mean
    d11[i,j-1]=mean(as.numeric(as.character(m[,j])))
  }
}
colnames(d11)=colnames(d4)[-1]
write.table(d11,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/mean/All_Function.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)  

##Gene均值，功能max 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f= read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/mean/All_Gene.txt",sep=""),header=T,sep="\t",row.names=1)
f1= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f1[,1]=as.character(f1[,1])
f1[,2]=as.character(f1[,2])
d11=data.frame(Function=unique(f1[,2]))
rownames(d11)<-d11[,1]
for(i in unique(f1[,2])){
  f2=f1[f1[,2]==i,1]
  for(j in 1:4){
    m1=data.frame(up=length(which(f[f2,j]>0)),down=length(which(f[f2,j]<0)),normal=length(which(f[f2,j]==0)))
    if(which(m1[1,]==max(m1)) ==3){
      d11[i,j+1]=0
    }else if(which(m1[1,]==max(m1)) ==2){
      d11[i,j+1]=-length(which(f[f2,j]<0))/length((f[f2,j]))
    }else if(which(m1[1,]==max(m1)) ==1){
      d11[i,j+1]=length(which(f[f2,j]>0))/length((f[f2,j]))
    }

  }
}
colnames(d11)=paste("Function",colnames(f),sep="\t")
write.table(d11,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/mean/All_Function_max.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)  





################################
##网络图构建PanCancer，
##点为每个基因，大小为在各亚型中的值，颜色为对应的功能，
##正方形为上调，圆形为下调，边为2个基因间的关系，<->为促进，
##->为抑制，-为NA或修改（不确定的）。
#install.packages("igraph")
library(igraph) 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in c("All")){
  f1<-read.csv(paste(path,"mean/",i,"_Gene.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$V2
  for(k in 2:ncol(f1)){
    f1$expressmean=f1[,k]
    f11=f1[intersect(as.character(f[,1]),as.character(f1$V2)),c("V2","expressmean")]
    m=merge(f11,f,by.x="V2",by.y="gene_Name")
    fu_colour=data.frame(Function=unique(m[,3]),
                         colour=c("orangered1","yellow3","snow3","palegoldenrod","maroon1","dodgerblue4","plum",
                                  "mediumvioletred","royalblue3","yellow","steelblue3","slateblue2","deepskyblue",
                                  "lightseagreen","rosybrown2","sienna1","tomato3","purple4","red","mediumpurple1"))
    m1=m[,c("V2","Function","expressmean")]
    #m1[which(m1[,3]>=0),4]="square"
    #m1[which(m1[,3]<0),4]="circle"
    m1[,4]="circle"
    colnames(m1)[4]="shape"
    #m1[which(m1[,3]>=0),5]="red"
    #m1[which(m1[,3]<0),5]="green"
    #colnames(m1)[5]="framecolor"
    m1=merge(m1,fu_colour)
    m1[,6]=as.character(m1$colour)
    m1[which(m1[,3]>0),6]="red4"
    m1[which(m1[,3]<0),6]="darkgreen"
    colnames(m1)[6]="framecolor"
    m1=m1[,c("V2","Function","expressmean","shape","colour","framecolor")]
    f2= read.csv(paste("D:/tp53/regulate.txt",sep=""),header=T,sep="\t")
    f2[,2]=gsub("^p53$","TP53",as.character((f2[,2])))
    f2[,3]=gsub("^p53$","TP53",as.character((f2[,3])))
    f2[,4]=as.character((f2[,4]))
    n1=unique(f2)
    for(j in 1:nrow(n1)){
      if(n1[j,4]=="-"){
        a=n1[j,2] 
        n1[j,2]=n1[j,3]
        n1[j,3]=a
        n1[j,4]=paste(n1[j,4],0,sep="")
        n1[j,"count"]=0.5
        n1[j,"linetype"]="->"
        n1[j,"colour"]="gray"
      }
      else if(n1[j,4]=="+"){
        n1[j,"count"]=0.5
        n1[j,"linetype"]="<->"
        n1[j,"colour"]="gray"
      }
      else{
        n1[j,"count"]=0.5
        n1[j,"linetype"]="-"
        n1[j,"colour"]="gray"
      }
    }
    n11=n1[,c("gene1","gene2","function.","count","linetype","colour")]
    data1=n11[(n11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))) & (n11[,2] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))),]
    data2=m1[m1[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1])),]
    
    g1 <- graph_from_data_frame(data1, directed=TRUE, vertices=data2)
    a=as.data.frame(degree(g1))
    a1=rownames(a)[which(a$`degree(g1)`!=0)]
    data3=data1[(data1[,1] %in% a1) & (data1[,2] %in% a1),]
    data4=data2[data2[,1] %in% a1,]
    g2 <- graph_from_data_frame(data3, directed=TRUE, vertices=data4)
    #layout=layout.circle\layout.fruchterman.reingold\layout.kamada.kawai 
    pdf(paste(path,"mean_graph/",i,"_",colnames(f1)[k],"_layout.kamada.kawai",".pdf",sep=""))  
    plot(g2, layout =layout.kamada.kawai ,vertex.label=V(g2)$V2,vertex.label.cex=0.3,
         vertex.size=3*(abs(V(g2)$expressmean)+0.5),vertex.label.dist=-0.4,vertex.label.color="black",vertex.color=V(g2)$colour,
         vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,vertex.frame.color=V(g2)$framecolor,
         edge.width=E(g2)$count/100,edge.color=E(g2)$colour,edge.arrow.size=0.05,main=paste(i,"_",colnames(f1)[k],sep=""))
    legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
    dev.off()
  }
}

tkplot(g2, layout =layout.circle ,vertex.label=V(g2)$V2,vertex.label.cex=0.3,
       vertex.size=3*(abs(V(g2)$expressmean)+0.5),vertex.label.dist=-0.4,vertex.label.color="black",vertex.color=V(g2)$colour,
       vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,vertex.frame.color=V(g2)$framecolor,
       edge.width=E(g2)$count/10,edge.color=E(g2)$colour,edge.arrow.size=0.05,main=paste(i,"_",colnames(f1)[k],sep=""))

#########################
##网络图构建PanCancer，删去不是上下调的基因
##点为每个基因，大小为在各亚型中的值，颜色为对应的功能，
##正方形为上调，圆形为下调，边为2个基因间的关系，<->为促进，
##->为抑制，-为NA或修改（不确定的）。
#install.packages("igraph")
library(igraph) 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in c("All")){
  f1<-read.csv(paste(path,"mean/",i,"_Gene.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$V2
  for(k in 2:ncol(f1)){
    f1$expressmean=f1[,k]
    f11=f1[intersect(as.character(f[,1]),as.character(f1$V2)),c("V2","expressmean")]
    m=merge(f11,f,by.x="V2",by.y="gene_Name")
    fu_colour=data.frame(Function=unique(m[,3]),
                         colour=c("orangered1","yellow3","snow3","palegoldenrod","maroon1","dodgerblue4","plum",
                                  "mediumvioletred","royalblue3","yellow","steelblue3","slateblue2","deepskyblue",
                                  "lightseagreen","rosybrown2","sienna1","tomato3","purple4","red","mediumpurple1"))
    m1=m[,c("V2","Function","expressmean")]
    m1[,4]="circle"
    colnames(m1)[4]="shape"
    #m1[which(m1[,3]>=0),5]="red"
    #m1[which(m1[,3]<0),5]="green"
    #colnames(m1)[5]="framecolor"
    m1=merge(m1,fu_colour)
    m1[,6]=as.character(m1$colour)
    m1[which(m1[,3]>0),6]="red4"
    m1[which(m1[,3]<0),6]="darkgreen"
    colnames(m1)[6]="framecolor"
    m1=m1[,c("V2","Function","expressmean","shape","colour","framecolor")]
    m11=m1[which(m1[,3]!=0 | as.character(m1[,2])=="TP53"),]
    f2= read.csv(paste("D:/tp53/regulate.txt",sep=""),header=T,sep="\t")
    f2[,2]=gsub("^p53$","TP53",as.character((f2[,2])))
    f2[,3]=gsub("^p53$","TP53",as.character((f2[,3])))
    f2[,4]=as.character((f2[,4]))
    n1=unique(f2)
    for(j in 1:nrow(n1)){
      if(n1[j,4]=="-"){
        a=n1[j,2] 
        n1[j,2]=n1[j,3]
        n1[j,3]=a
        n1[j,4]=paste(n1[j,4],0,sep="")
        n1[j,"count"]=0.5
        n1[j,"linetype"]="->"
        n1[j,"colour"]="gray"
      }
      else if(n1[j,4]=="+"){
        n1[j,"count"]=0.5
        n1[j,"linetype"]="<->"
        n1[j,"colour"]="gray"
      }
      else{
        n1[j,"count"]=0.5
        n1[j,"linetype"]="-"
        n1[j,"colour"]="gray"
      }
    }
    n11=n1[,c("gene1","gene2","function.","count","linetype","colour")]
    data1=n11[(n11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m11[,1]))) & (n11[,2] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m11[,1]))),]
    data2=m11[m11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m11[,1])),]
    
    g1 <- graph_from_data_frame(data1, directed=TRUE, vertices=data2)
    a=as.data.frame(degree(g1))
    a1=rownames(a)[which(a$`degree(g1)`!=0)]
    data3=data1[(data1[,1] %in% a1) & (data1[,2] %in% a1),]
    data4=data2[data2[,1] %in% a1,]
    g2 <- graph_from_data_frame(data3, directed=TRUE, vertices=data4)
    #layout=layout.circle\layout.fruchterman.reingold\layout.kamada.kawai 
    pdf(paste(path,"mean_graph0/",i,"_",colnames(f1)[k],"_layout.kamada.kawai",".pdf",sep=""))  
    plot(g2, layout =layout.kamada.kawai ,vertex.label=V(g2)$V2,vertex.label.cex=0.3,
         vertex.size=3*(abs(V(g2)$expressmean)+0.5),vertex.label.dist=-0.4,vertex.label.color="black",vertex.color=V(g2)$colour,
         vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,vertex.frame.color=V(g2)$framecolor,
         edge.width=E(g2)$count/100,edge.color=E(g2)$colour,edge.arrow.size=0.05,main=paste(i,"_",colnames(f1)[k],sep=""))
    legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
    dev.off()
  }
}


##############################
##PanCancer的circlize图  ##D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/
##123/2/D/PanCancer/mean_graph1/
##All_C1.pdf所有功能
##C1.pdf删去TP53功能
##S1.pdf删去TP53功能且删去正常所占比例，只考虑上下调所占最大比例
library(circlize)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/D/PanCancer/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
for(i in c("All")){
  f3<-read.csv(paste(path,"mean/",i,"_Function.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  f1<-read.csv(paste(path,"mean/",i,"_Gene.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$V2
  for(k in 2:ncol(f1)){
    f1$expressmean=f1[,k]
    f11=f1[intersect(as.character(f[,1]),as.character(f1$V2)),c("V2","expressmean")]
    m=merge(f11,f,by.x="V2",by.y="gene_Name")
    fu_colour=data.frame(Function=unique(m[,3]),
                         colour=c("orangered1","yellow3","snow3","palegoldenrod","maroon1","dodgerblue4","plum",
                                  "mediumvioletred","royalblue3","yellow","steelblue3","slateblue2","deepskyblue",
                                  "lightseagreen","rosybrown2","sienna1","tomato3","purple4","red","mediumpurple1"))
    a2=merge(fu_colour,f3[,c(1,k)])
    rownames(a2)=a2$Function
    a3=a2[fu_colour[,1],]
    a3[a3[,3]>0,4]="red4"
    a3[a3[,3]==0,4]="seashell3"
    a3[a3[,3]<0,4]="darkgreen"
    a5=a3[-grep("TP53",a3[,1]),]  #删去TP53功能
    a4=cbind(a5,start=1:nrow(a5),end=2:(nrow(a5)+1))
    a4$Function=as.character(a4$Function)
    a4$colour=as.character(a4$colour)
    
    #colnames(f1)[k]
    circos.initializeWithIdeogram(plotType = NULL) #去原始圈
    circos.par("gap.degree"=0.01) #各板块间隔
    circos.genomicInitialize(a4[,c(1,5,6)])  
    circos.genomicTrackPlotRegion(a4[,c(1,5,6)],ylim=c(-1,1),track.height=0.1,bg.col=a4$colour,bg.border = NA)
    circos.genomicTrackPlotRegion(a4[,c(1,5,6)],ylim=c(-1,1),track.height=0.1,bg.col="white",bg.border ="gray") 
    circos.genomicPoints(a4[,c(5,6)], 0, cex = 1.5*(abs(a4[,3])+0.5), pch = 16, col =a4[,4])
    circos.clear() 
  }
} 


















 
                        

###################################################
##对Allcancer26_class3_manhattan_ward.D_1075对每类加圆circlize加癌症亚型
##Allcancer2_class4_manhattan_ward.D_1075
##Allcancer2_class4_minkowski_ward.D2_1075
#install.packages("circlize")
library(circlize)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f1=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/26/Allcancer26_class3_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
f1[,1]=as.character(f1[,1])
f<-read.csv("D:/tp53/phenotype_Analysis/subtype/subtype.txt",header=T,sep="\t")
f[,1]=substring(as.character(f[,1]),1,15)
a=f[-grep("-11",f[,1]),c("pan.samplesID","Subtype_Selected")]
a[,1]=substring(gsub("-",".",a[,1]),0,12)
a[,2]=as.character(a[,2])

a2=f1
for(i in unique(gsub("\\..*","",as.character(a$Subtype_Selected)))[c(-2,-6)]){#18ge
  a1=a[grep(i,as.character(a$Subtype_Selected)),]
  colnames(a1)[2]=i
  b=data.frame(cbind(unique(a1[,2]),rainbow(length(unique(a1[,2])))))
  colnames(b)[1]=i
  colnames(b)[2]=paste(i,"_colour",sep="")
  b[,1]=as.character(b[,1])
  b[,2]=as.character(b[,2])
  b1=merge(a1,b)
  a2=unique(merge(a2,b1,by.x="V1",by.y="pan.samplesID",all.x=T))
}  

fu<-function(x) sum(is.na(x))  
m=data.frame(apply(a2,2,fu))
m1=rownames(m)[which(m[,1]!=1076)] 
a3=a2[order(a2$V2,decreasing=F),m1]
a4=cbind(V=1:nrow(a3),start=1:nrow(a3),end=2:(nrow(a3)+1),a3)
 
#par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
#circos.par(track.height = 0.1)
circos.initializeWithIdeogram(plotType = NULL) #去原始圈
circos.par("gap.degree"=0.01) #各板块间隔
circos.genomicInitialize(a4[,1:3]) #原始圈250\655\812\950
#circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=c(rep("pink",data.frame(table(a4$V2))[1,2]),rep("yellow",data.frame(table(a4$V2))[2,2]),rep("lightgreen",data.frame(table(a4$V2))[3,2]),rep("lightblue",data.frame(table(a4$V2))[4,2])),bg.border=NA,track.height=0.05)
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=c(rep("pink",data.frame(table(a4$V2))[1,2]),rep("yellow",data.frame(table(a4$V2))[2,2]),rep("lightgreen",data.frame(table(a4$V2))[3,2])),bg.border=NA,track.height=0.05)
circos.text(250, 0.5, "C1", sector.index = "1000", track.index = 2) 
circos.text(655, 0.5, "C2", sector.index = "1000", track.index = 2) 
circos.text(950, 0.5, "C3", sector.index = "1000", track.index = 2) 
circos.text(1050, 0.5, "C4", sector.index = "1000", track.index = 2) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$BLCA_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "BLCA", sector.index = "1000", track.index = 3) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$BRCA_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "BRCA", sector.index = "1000", track.index = 4) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$GBM_LGG_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "GBM_LGG", sector.index = "1000", track.index = 5) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$HNSC_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "HNSC", sector.index = "1000", track.index = 6) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$KIRC_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "KIRC", sector.index = "1000", track.index = 7) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$KIRP_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "KIRP", sector.index = "1000", track.index = 8) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$LIHC_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "LIHC", sector.index = "1000", track.index = 9) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$LUAD_colour,bg.border=NA,track.height=0.05)
circos.text(260, 0.5, "LUAD", sector.index = "500", track.index = 10) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$LUSC_colour,bg.border=NA,track.height=0.05)
circos.text(260, 0.5, "LUSC", sector.index = "1000", track.index = 11) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$PRAD_colour,bg.border=NA,track.height=0.05)
circos.text(220, 0.5, "PRAD", sector.index = "1000", track.index = 12) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$SKCM_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "SKCM", sector.index = "1000", track.index = 13) 
circos.genomicTrackPlotRegion(ylim=c(0,1),bg.col=a4$UCEC_colour,bg.border=NA,track.height=0.05)
circos.text(910, 0.5, "UCEC", sector.index = "1000", track.index = 14) 
circos.clear()


#e2=data.frame(class=sort(unique(f1[,2])),start=c(1,1+nrow(f1[f1$V2==1,]),1+nrow(f1[f1$V2==1,])+nrow(f1[f1$V2==2,]),1+nrow(f1[f1$V2==1,])+nrow(f1[f1$V2==2,])+nrow(f1[f1$V2==3,])),end=c(nrow(f1[f1$V2==1,]),nrow(f1[f1$V2==1,])+nrow(f1[f1$V2==2,]),nrow(f1[f1$V2==1,])+nrow(f1[f1$V2==2,])+nrow(f1[f1$V2==3,]),nrow(f1)))
circos.genomicInitialize(e2) #原始圈
circos.initializeWithIdeogram(plotType = NULL) #去原始圈
circos.par("gap.degree"=0.01) #各板块间隔、起始角度
circos.par("start.degree"=30)
circos.clear()


library(circlize)
# 简单创建一个数据集
set.seed(999)
n <- 1000
a <- data.frame(factors = sample(letters[1:8], n, replace = TRUE), x = rnorm(n), y = runif(n))   
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.6)
circos.par(track.height = 0.1)
circos.initialize(factors = a$factors, x = a$x)
#初始化，factors来控制track数目，初始化里只有x， 没有y。这一步相当于ggplot()
circos.trackPlotRegion(factors = a$factors, y = a$y,
                       panel.fun = function(x, y) {
                         circos.axis()})
col <- rep(c("#FF0000", "#00FF00"), 4) #自定义一下颜色
# 这里先解释一下，一个track有好几个cell，具体数目由factors决定的，向本数据集中factors有八个，
#因此绘制一个track，其包含八个cell。含有前缀circos.track的函数会在所有的cel里添加基本元素，
#而只有前缀circos.的函数可以在特定的track、cell里添加基本元素。具体看下演示。
circos.trackPoints(a$factors, a$x, a$y, col = col, pch = 16, cex = 0.5) #所有的cell里都绘制点图
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1) #在track 1中的标记为a的cell里添加
circos.text(1, 0.5, "right", sector.index = "a")
bg.col <- rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factors, a$x, bg.col = bg.col, col = NA)
circos.trackPlotRegion(factors = a$factors, x = a$x, y = a$y,
                       panel.fun = function(x, y) {
                         grey = c("#FFFFFF", "#CCCCCC", "#999999")
                         sector.index = get.cell.meta.data("sector.index") #这个是第三个track，因为我们刚刚创建，这里这一步不用也可。
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), mean(ylim), sector.index)
                         circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
                         circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)})
# update第2个track中标记为d的sector
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim <- get.cell.meta.data("xlim")
ylim <- get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.trackPlotRegion(factors = a$factors, y = a$y)
circos.trackLines(a$factors[1:100], a$x[1:100], a$y[1:100], type = "h")
circos.link("a", 0, "b", 0, h = 0.3) #point to point
circos.link("c", c(-0.5, 0.5), "d", c(-0.5, 0.5), col = "red", border = NA, h = 0.2) #intreval to interval
circos.link("e", 0, "g", c(-1, 1), col = "green", border = "black", lwd = 2, lty = 2) #point to interval
circos.clear()
















 
#########################################
for(clu in c("hc","pam","km","kmdist")){
  for(dis in c("pearson","spearman","euclidean","binary","maximum","canberra","minkowski","manhattan")){
    tryCatch({
      results <- ConsensusClusterPlus(cancerimm,maxK = 10,reps = 100,pItem = 0.8,pFeature = 1,title = paste0(clu,"-",dis),
                                      clusterAlg = clu,distance = dis,seed = 1262118388.71279,
                                      plot = "pdf",writeTable = TRUE)
      icl <- calcICL(results,title = paste0(clu,"-",dis),plot = "pdf",writeTable = TRUE)
    },
    error = function(e){cat("ERROR :",clu," ",dis," ",conditionMessage(e),"\n")}) 
  }
}

######################################
#对123模式Allcancer基因卡从2:60，聚类从2:10，采取不同聚类方法，不同距离，对每类做生存加时间均值。
library("survival")
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
b11=subset(ff,ff$vital_status.y==0)
b21=subset(ff,ff$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(data.frame(e1))

path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])
for(number in 2:60){
  fu<-function(x) sum(x!=0)  
  a1=data.frame(apply(h1,1,fu))
  a2=rownames(a1)[which(a1[,1]>=number)] 
  h2=na.omit(h[a2,-1])
  for(nu in 2:10){
    for(m1 in c("manhattan","euclidean","minkowski","canberra")){
       tryCatch({
         for(m2 in c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")){
           d = dist(t(h2), method =m1)
           hcward = hclust(d, method=m2)
           n= cutree(hcward,k=nu) 
           f2=cbind(V1=rownames(as.data.frame(n)),as.data.frame(n))
           e2=na.omit(merge(f2,e1,by.x="V1",by.y="X1",all.x=T))
           e2[,2]=factor(e2[,2])
           e2[,3]=as.numeric(as.character(e2[,3]))
           e2[,4]=as.numeric(as.character(e2[,4]))
           if( TRUE %in% ((length(as.numeric(unique(e2[,2])))>1) & (length(as.numeric(unique(e2[,3])))==2))){
             dif1 <- survdiff(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2])
             kmsurvival1<-survfit(Surv(as.numeric(e2[,4]),as.numeric(e2[,3]))~e2[,2],conf.type = "log-log")
             if(round(pchisq(dif1$chisq,1,lower.tail=F),4)<=0.05){
               pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123_new/survival/",number,"_Allcancer_class",nu,"_",m1,"_",m2,"_survival.pdf",sep=""))  
               plot(kmsurvival1, lty = 'solid', col=c("red","orange","blue","purple","burlywood4","deeppink","gold","green","salmon4","black"),lwd=1.1,
                    xlab='survival time in days',ylab='survival probabilities',main=paste(number,"Allcancer_class",nu,m1,m2,sep="_"))
               legend('bottomleft', cex=0.6,text.width=0.4,gsub(".*=","",as.character(data.frame(dif1$n)[,1])), lty='solid',
                      col=c("red","orange","blue","purple","burlywood4","deeppink","gold","green","salmon4","black"))
               text(600,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
               for(j in 1:length(as.numeric(unique(e2[,2]))) ){
                 text(600,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
               }
               for(jj in 1:length(as.numeric(unique(e2[,2]))) ){
                 text(2000,0.9-(jj-1)*0.03,cex=0.8,round(mean(e2[as.numeric(as.character(e2[,2]))==jj,4]),1))
               }
               dev.off() 
               write.table(as.data.frame(n),file =paste(path,"p53_context_discreteMode_fc1.5_heatmap/123_new/data/",number,"_Allcancer_class",nu,"_",m1,"_",m2,".txt",sep=""), row.names = T,col.names=F, quote = F,sep = "\t",append=TRUE)
             }
             
           }
         }

       },
      error = function(e){cat("ERROR :",m1,m2,"  ",conditionMessage(e),"\n")})
    }
  }
}



 


######################################
#对123模式Allcancer基因卡从2:60，聚类从2:10，采取不同聚类方法，不同距离，对每类计算轮廓系数silhouette。
#install.packages("fpc")
#install.packages("factoextra")
library(ConsensusClusterPlus)
library(fpc) # sil
library(cluster)
library(factoextra)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
filename<-list.files(paste(path,"p53_context_discreteMode_fc1.5/123/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
h=data.frame(name=as.character(f[,1]))
for(cancer in filename){
  f1<-read.csv(paste(path,"p53_context_discreteMode_fc1.5/123/",cancer,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1[,1]  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])
for(number in 2:60){
  fu<-function(x) sum(x!=0)  
  a1=data.frame(apply(h1,1,fu))
  a2=rownames(a1)[which(a1[,1]>=number)] 
  h2=na.omit(h[a2,-1])
  for(nu in 2:10){
    for(m1 in c("manhattan","euclidean","minkowski","canberra")){
      tryCatch({
        for(m2 in c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")){
          d = dist(t(h2), method =m1)
          hcward = hclust(d, method=m2)
          n= cutree(hcward,k=nu) 
          class=as.data.frame(n)
          sil <- silhouette(class[,1], d)
          pdf(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123_new/silhouette/",number,"_Allcancer_class",nu,"_",m1,"_",m2,".pdf",sep=""))  
          plot(sil,main =paste(number,"_Allcancer_class",nu,"_",m1,"_",m2,"_Silhouette",sep=""),col=c(terrain.colors(10)[c(1,3,5,7,9)],topo.colors(10)[c(2,4,6,8,10)])[1:nu])
          abline(v = mean(sil[,3]),lty = 2,col = "gray30")
          dev.off() 
        }
      },
      error = function(e){cat("ERROR :",m1,m2,"  ",conditionMessage(e),"\n")})
    }
  }
}

 










##对每个cancer\all cancer上下调基因的个数
##对每个cancer某个基因的上下调的个数、分页画图
##将F盘tp53/Gene Expression移动到D盘tp53/Gene Expression
#path="F:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_updown/",sep=""))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_updown/",filename[2],sep=""),header=T,sep="\t",stringsAsFactors = F)
for(i in filename[3:length(filename)]){
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_updown/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  f1=merge(f1,f2,all=T,by="V1")
}
f1[sapply(f1,is.na)]<-0
rownames(f1)=f1[,1]
ff=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
x1=ff[,2:ncol(ff)]
y1<-as.matrix(x1)
row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
colnames(y1)=paste(rep(c(gsub("_fc_updown.txt","",filename[2:length(filename)])),each=6),rep(1:6,6),sep="")
#col_anno = data.frame(CancerClass = factor(rep(c(gsub("_fc_updown.txt","",filename[2:length(filename)])),each=6)),row.names=colnames(ff)[2:ncol(ff)])
#pheatmap(y1,cluster_row = FALSE,cluster_col = FALSE,annotation_row = row_anno,fontsize_row=5,fontsize_col=5,
         #display_numbers = TRUE,fontsize = 1,border_color=NA)  
#for(i in 1:nrow(y1)){
  #pdf(file=paste(path,"p53_context_fpkm_foldchange_updown_plot/",rownames(y1)[i],".pdf",sep=""),width=10,height=8)
  #symbols(1:ncol(y1),y1[i,],circle = rep(c(sqrt(1.5/pi),sqrt(1.5/pi),sqrt(1/pi),sqrt(1/pi),sqrt(0.5/pi),sqrt(0.5/pi)),ncol(y1)/6),bg=rep(c("red","green"),ncol(y1)/2),inches = FALSE,xlab="cancer",ylab="number_updown",main=rownames(y1)[i])
  #text(1:ncol(y1),y1[i,],y1[i,],cex=0.6)
  #text(seq(1,ncol(y1),6),-10,c(gsub("_fc_updown.txt","",filename[2:length(filename)])),cex=0.5)
  #legend('topright', c("up","down"),pch=20,col=c("red","green"),cex=0.5)
  #legend('top', c("CancerSample","MuSample","NoMuSample"),pch=c(19,16,20),col="red",cex=0.5)
  #dev.off()
#}

b=data.frame()
for(i in 1:nrow(y1)){
  b=rbind(b,cbind(rep(rownames(y1)[i],ncol(y1)),colnames(y1),rep(c("above","down"),3),rep(c("CancerSample","MuSample","NoMuSample"),each=2),y1[i,]))
} 
colnames(b)=c("gene","cancer","class1","class2","number_updown")
b$cancer=gsub("1","",b$cancer)
for(j in as.character(unique(f[,2])[-1])){
  d=rbind(subset(f,f[,2]=="TP53"),subset(f,f[,2]==j))
  d1=subset(b,as.character(b[,1]) %in% as.character(d[,1]))
  Expression_values=factor(d1$class1)
  Sample_Type=factor(d1$class2)
  p <- ggplot(d1,aes(x = d1$cancer,y = as.numeric(as.character(d1$number_updown))))
  p + facet_grid(gene ~ .) +
    geom_point(aes(colour = Expression_values,shape=Sample_Type),size =0.7) +
    scale_color_manual(values=c("red","green")) +
    theme(axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 3),
          strip.text = element_text(size = 5),
          panel.background = NULL)+
    geom_vline(xintercept = seq(1,nrow(d1),6)-0.6, col="grey",lwd=0.01)+ 
    xlab("The Type of Cancer") + ylab("The Number of Disregulation Genes") + ggtitle(gsub("/","&",as.character(unique(d[,2])))[2])+
    theme(plot.title=element_text(hjust = 0.5),legend.position="top")+
    scale_x_discrete(breaks = d1$cancer[seq(1,nrow(d1),6)]) 
  ggsave(file=paste(path,"p53_context_fpkm_foldchange_updown_plot/2/",gsub("/","&",as.character(unique(d[,2])))[1],"_",gsub("/","&",as.character(unique(d[,2])))[2],".pdf",sep=""),width = 10,height=9)
}
 





###TCGAbiolinks GDCdata的Gene_count数据处理
PATH="F:/tp53/Gene Expression/"
f<-read.csv(paste(PATH,"ALL_samplesWithGeneCountData.txt",sep=""),header=T,sep="\t")
f[,2]=NA
a=data.frame()
cancer<-c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
for (k in cancer){
  b<-subset(f,as.character(f[,3])==k)
  if(paste(as.character(k),'_Gene_count.txt',sep="") %in% list.files("F:/tp53/Gene Expression/Gene_count1/")){
    f1<-read.csv(paste("F:/tp53/Gene Expression/Gene_count1/",paste(as.character(k),'_Gene_count.txt',sep=""),sep=""),header=T,sep="\t")
    for(i in 1:nrow(b)){
      if(as.character(b[i,1]) %in% gsub("\\.","-",substring(colnames(f1),1,12))){
        if(length(grep("^TP53$",as.character(f1[,1])))!=0){
          b[i,2]=1
        }else{
          b[i,2]=0
        }
      }
    }
  }
  a=rbind(a,b)
}
write.table(a,file =paste(PATH,"ALL_samplesWithGeneCountData1.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


###从TCGA下载27个cancer的miRNA数据
#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
#"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
library(TCGAbiolinks)
setwd("F:/tp53/miRNA")
cancer <- "UCS"
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "Transcriptome Profiling",  #-········· 数据种类
                    data.type = "miRNA Expression Quantification") #······· 数据类型)   
  GDCdownload(query)  #·························································· GDC下载数据
}


#-- 将基因miRNA_count整合成data.frame数据
library(TCGAbiolinks)
cancer <- c("BLCA","BRCA","CESC","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "Transcriptome Profiling",  #·········· 数据种类
                    data.type = "miRNA Expression Quantification")  #······ 数据类型
  setwd(paste("F:/tp53/miRNA/GDCdata/TCGA-",k,"/harmonized/Transcriptome_Profiling/miRNA_Expression_Quantification",sep = ""))
  genecount <- read.table(paste(query$results[[1]]$file_id[1],"/",query$results[[1]]$file_name[1],sep = ""),sep = "\t",header = T)
  colnames(genecount)[2] <- query$results[[1]]$cases[1]
  genecount <- genecount[,1:2]
  for (i in 2:length(query$results[[1]]$file_id)){
    genecount_1 <- read.table(paste(query$results[[1]]$file_id[i],"/",query$results[[1]]$file_name[i],sep = ""),sep = "\t",header = T)
    colnames(genecount_1)[2] <- query$results[[1]]$cases[i]
    genecount_1 <- genecount_1[,1:2]
    genecount <- merge(genecount,genecount_1,all = T)
  }
  write.table(genecount,file = paste("F:/tp53/miRNA/miRNA_count/",k,"_miRNA_count.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
}


#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
#-- 将基因miRNA_fpkm整合成data.frame数据
library(TCGAbiolinks)
cancer <- "BRCA"
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "Transcriptome Profiling",  #·········· 数据种类
                    data.type = "miRNA Expression Quantification")  #······ 数据类型
  setwd(paste("F:/tp53/miRNA/GDCdata/TCGA-",k,"/harmonized/Transcriptome_Profiling/miRNA_Expression_Quantification",sep = ""))
  genecount <- read.table(paste(query$results[[1]]$file_id[1],"/",query$results[[1]]$file_name[1],sep = ""),sep = "\t",header = T)
  genecount <-genecount[,-2]
  colnames(genecount)[2] <- query$results[[1]]$cases[1]
  genecount <- genecount[,1:2]
  for (i in 2:length(query$results[[1]]$file_id)){
    genecount_1 <- read.table(paste(query$results[[1]]$file_id[i],"/",query$results[[1]]$file_name[i],sep = ""),sep = "\t",header = T)
    genecount_1 <-genecount_1[,-2]
    colnames(genecount_1)[2] <- query$results[[1]]$cases[i]
    genecount_1 <- genecount_1[,1:2]
    genecount <- merge(genecount,genecount_1,all = T)
  }
  write.table(genecount,file = paste("F:/tp53/miRNA/miRNA_fpkm/",k,"_miRNA_fpkm.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = T)
}



###从TCGA下载27个cancer的DNA Methylation数据
#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
#"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
library(TCGAbiolinks)
setwd("F:/tp53/Methylation/450")
cancer <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
            "LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
            "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
for (k in cancer){
  query <- GDCquery(project = paste("TCGA-",k,sep = ""),     #············· 癌症类型(project id)
                    data.category = "DNA Methylation",  #-········· 数据种类
                    platform = "Illumina Human Methylation 450") 
  GDCdownload(query)  #·························································· GDC下载数据
}
























##protein change对COSMIC位置，得到patient对应的表型、tp53
##选出所有的protein change
filenames<-list.files("F:/tp53/TCGA/",pattern = '*_tp53.csv') 
da={}
for(i in filenames){
  f<-read.csv(paste("F:/tp53/TCGA/",i,sep=""),header=T)
  da=append(da,as.character(f[,3]))
}
#da[!duplicated(da)] 
write.table(unique(da), file = "F:/tp53/protein change.txt", row.names = F,col.names=F, quote = F,sep = "", append = TRUE)

##protein change对旧的CosmicMutantExportCensus.tsv位置
f1<-read.table("F:/tp53/protein change.txt",sep="\t",header=F) 
f2<-read.csv("F:/tp53/CosmicMutantExportCensus.tsv",sep="\t",header=T)  
f3<-unique(cbind(as.character(f2[,17]),as.character(f2[,19])))
#strsplit(as.character(f2[,15]), "p.")
#substr(as.character(f2[,15]), 3,3:1000)
#gsub("^..|\\?","",f2[,15])#字符串分割
i<-intersect(paste("p.",as.character(f1[,1]),sep=""),f3[,2])
row.names(f1)=paste("p.",as.character(f1[,1]),sep="")
row.names(f3)=f3[,2]
f1[i,2]=f3[i,1]
write.table(f1, file="F:/tp53/protein change_COSM1.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

##protein change对新的CosmicMutantExport.tsv位置
f1<-read.table("F:/tp53/protein change.txt",sep="\t",header=F) 
f2<-read.csv("F:/tp53/CosmicMutantExport.tsv",sep="\t",header=T)  
#17Mutation.ID COSM18232 
#19Mutation.AA p.N363N
#24Mutation.genome.position 10:26541626-26541626
f3<-unique(cbind(as.character(f2[,17]),as.character(f2[,19])))
i<-intersect(paste("p.",as.character(f1[,1]),sep=""),f3[,2])
row.names(f1)=paste("p.",as.character(f1[,1]),sep="")
row.names(f3)=f3[,2]
f1[i,2]=f3[i,1]
write.table(f1, file="F:/tp53/protein change_COSM2.txt",row.names = F,col.names=F,quote = F,sep = "\t",append=TRUE)

##网页提取信息
f<-read.table("F:/tp53/protein change_COSM1.txt",sep="\t",header=F)
for (i in 1:dim(f)[1]){
  if(!is.na(gsub("^....","",as.character(f[i,2])))){
    URL = paste("http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",gsub("^....","",as.character(f[i,2])),sep="")
    web <- readLines(URL,encoding="UTF-8")
    a<-strsplit(gsub("^..","",strsplit(strsplit(web[176], ">")[[1]][3],"<")[[1]][1]),",")[[1]][1]#"Substitution - Missense"
    b<-strsplit(gsub("^..","",strsplit(strsplit(web[176], ">")[[1]][3],"<")[[1]][1]),",")[[1]][2]#" position "
    c<-strsplit(strsplit(web[176], ">")[[1]][4],"<")[[1]][1]#"194"
    d<-strsplit(strsplit(web[177], ">")[[1]][2],"<")[[1]][1]#"L"
    e<-strsplit(strsplit(web[177], ">")[[1]][6],"<")[[1]][1]#"F"
    f[i,3]=a
    f[i,4]=paste(b,c,sep="")
    f[i,5]=paste(d,"_",e,sep="")
  }
}
write.table(f, file="F:/tp53/protein change_COSM1_Mut.txt", row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

#############################################

#   红All、蓝Mu、绿NoMu
#  "#FF331C","#1DADB8","#12B826"
#"#12B826"没突变,"#1DADB8"突变,"#5e227f"单,"#d22780"多

#"indianred3","ivory4"


#正常样本，疾病样本，突变样本，不突变样本
#"mediumblue","violetred3","wheat3","seagreen3"
#"paleturquoise","pink","grey","springgreen1"

#4yax
#"deeppink3","goldenrod3","palegreen4","royalblue4"
