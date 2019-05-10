####################################
######对Allcancer2_class4_manhattan_ward.D_1075_1进行多组学数据处理1075
##1、mu对突变数据进行处理整合成Allcancer对应的全基因组的突变信息
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/TCGA-",i,".txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f3$name=rownames(f3)
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("name",f4)])))
  f1=f5
}
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  f3 <- read.table(paste0("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/TCGA-",i,".txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f3$name=rownames(f3)
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("name",f4)])))
  f1=merge(f1,f5,all=T)
}
ff=unique(f1[-grep("^TP53$",as.character(f1[,1])),])
ff[sapply(ff,is.na)]<-0
write.table(ff,file ="F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)

##2、mu对突变数据进行处理、用超几何检验
rm(list = ls())  #删去上面的程序
class <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu.txt",header=T,sep="\t")
rownames(mutation)=as.character(mutation[,1])
mutation1 <- mutation[,intersect(colnames(mutation),class$V1)]
class <- class[class$V1 %in% colnames(mutation1),]
classlist <- split(class$V1,class$V2)
test <- matrix(0,nrow = nrow(mutation1),ncol = 4*2,dimnames = list(rownames(mutation1),c(paste("Subtype",c(1:4)),paste("Sample",c(1:4)))))
for (subtype in names(classlist)){
  test[,paste("Subtype",subtype)] <- apply(mutation1,1,function(x){phyper(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"]),
                                                                          length(classlist[[subtype]]),
                                                                          length(x)-length(classlist[[subtype]]),
                                                                          as.numeric(table(as.character(x))["1"]),lower.tail = FALSE)})
  test[,paste("Sample",subtype)] <- apply(mutation1,1,function(x){paste(ifelse(is.na(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                               0,as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                        length(classlist[[subtype]]),
                                                                        length(x)-length(classlist[[subtype]]),
                                                                        ifelse(is.na(as.numeric(table(as.character(x))["1"])),
                                                                               0,as.numeric(table(as.character(x))["1"])),
                                                                        sep = ";")})
}
write.table(test,file = "F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper.txt",col.names = T,row.names = T,quote = F,sep = "\t")

##3、mu突变数据筛选突变率>= 0.1 & mutation$phyper_pvalue <= 0.05\0.01 样本数>=10
rm(list = ls())
library(reshape2)
library(stringr)
mutationselect <- data.frame()
#mutationall <- data.frame()
mutation <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper.txt",sep = "\t",header = T)
mutation$Gene <- rownames(mutation)
mutation$Cancer <- "All"
mutation$Method <- "manhattan"
mutation$Group <- as.character(4)
mutation <- melt(mutation,variable.name = "Subgroup",value.name = "phyper_pvalue")
mutation$Sample <- apply(mutation,1,function(x){as.character(x[grep(paste0(str_sub(x[grep("Subgroup",colnames(mutation))],-1,-1),"$"),colnames(mutation))])})
mutation <- mutation[,-grep("Sample\\.",colnames(mutation))]
mutation <- na.omit(mutation)
mutation$ratio <- unlist(lapply(str_split(mutation$Sample,";"),function(x){ifelse(as.numeric(x[2]) >= 10 & as.numeric(x[3]) >= 10,
                                                                                  as.numeric(x[1])/as.numeric(x[2]),
                                                                                  NA)}))
mutation <- na.omit(mutation)
#mutationselect <- rbind(mutationselect,mutation[mutation$ratio >= 0.1 & mutation$phyper_pvalue <= 0.05,])
#mutationall <- rbind(mutationall,mutation) 
a=rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.05,])
b=rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.01,])
d=rbind(mutationselect,mutation)
write.table(a,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper_select_0.05.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(b,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper_select_0.01.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(d,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper_select_All_1075.txt",col.names = T,row.names = F,quote = F,sep = "\t")


##4、mu画图
## 引用分组时每组样本个数大于10，秩和检验<0.01基因、频率为这个基因在此亚型中突变的频率
rm(list = ls())
library(ggplot2)
library(grid)
# 突变谱及亚型分类
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
class=data.frame(manhattan.4=f[,"V2"])
rownames(class)=as.character(f[,"V1"])
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu.txt",header=T,sep="\t",row.names = 1)
# 删选后的基因
gene <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper_select_0.01.txt",header = T,stringsAsFactors = F)
gene <- gene[unlist(lapply(lapply(lapply(strsplit(as.character(gene$Sample),";"),as.numeric),function(i){i >= 10}),all)),]
gene <- gene[!duplicated(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  genelist1 <- split(genelist[[i]],genelist[[i]]$Subgroup)
  gene1 <- lapply(genelist1,function(x){x$Gene})
  sample <- lapply(genelist1,function(x){rownames(class)[class[,gsub("-",".",i)] %in% gsub("Subtype.","",unique(x$Subgroup))]})
  for(j in names(genelist1)){
    genelist1[[j]]$Freq <- apply(mutation[gene1[[j]],sample[[j]]],1,sum)/length(mutation[gene1[[j]],sample[[j]]])
  }
  genelist[[i]] <- unsplit(genelist1,genelist[[i]]$Subgroup)
  
  p1 <- ggplot(data = genelist[[i]],aes(x = Gene,y = Freq)) +
    geom_bar(fill = "steelblue",stat = "identity") +
    geom_bar(color = "steelblue",fill = "transparent",stat = "identity",position = "fill",size = 1) +
    coord_flip()  +
    scale_y_continuous(labels = scales::percent) +
    #geom_text(aes(label = paste0(signif(Freq,2)*100,"%")), position = position_fill(vjust = 0.5)) +
    facet_grid(~Subgroup) + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=rel(1.1),face="bold"),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  
  p2 <- ggplot(data = genelist[[i]],aes(x = Freq)) +
    geom_density(aes(fill = Subgroup),alpha = 0.8)+
    facet_grid(~Subgroup) + 
    scale_fill_manual(values = c("#3A94BB","#30B89D","#A7C46A","#F9B042","#E05B3E","#9D5F74","#ECC9C9","#E88080","#E8D380","#D87A80")[1:length(unique(genelist[[i]]$Subgroup))]) +
    scale_x_continuous(limits = c(0,1)) +
    labs(title = i) +
    theme(legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())
  
  ########新建画图页面###########
  #grid.newpage()  ##新建页面
  pdf("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_figure_0.01.pdf",width = 15,height = 15)
  pushViewport(viewport(layout = grid.layout(9,1))) ####将页面分成2*2矩阵
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(p1, vp = vplayout(2:9,1))   ###将（1,1)和(1,2)的位置画图c
  print(p2, vp = vplayout(1,1))   ###将(2,1)的位置画图b
  dev.off() ##画下一幅图，记得关闭窗口
}
#gene <- unsplit(genelist,paste(gene$Method,gene$Group,sep = "."))
#write.table(gene,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/All_mu_p.txt",col.names = T,row.names = F,sep = "\t")





#######Allcancer2_class4_manhattan_ward.D_1进行多组学数据处理3009
##1、mu对突变数据进行处理整合成Allcancer对应的全基因组的突变信息
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/TCGA-",i,".txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f3$name=rownames(f3)
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("name",f4)])))
  f1=f5
}
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  f3 <- read.table(paste0("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/TCGA-",i,".txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f3$name=rownames(f3)
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("name",f4)])))
  f1=merge(f1,f5,all = T)
}
ff=unique(f1[-grep("^TP53$",as.character(f1[,1])),])
ff[sapply(ff,is.na)]<-0
write.table(ff,file ="F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)

##2、mu对突变数据进行处理、用超几何检验
rm(list = ls())  #删去上面的程序
class <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu.txt",header=T,sep="\t")
rownames(mutation)=as.character(mutation[,1])
mutation1 <- mutation[,intersect(colnames(mutation),class$V1)]
class <- class[class$V1 %in% colnames(mutation1),]
classlist <- split(class$V1,class$V2)
test <- matrix(0,nrow = nrow(mutation1),ncol = 4*2,dimnames = list(rownames(mutation1),c(paste("Subtype",c(1:4)),paste("Sample",c(1:4)))))
for (subtype in names(classlist)){
  test[,paste("Subtype",subtype)] <- apply(mutation1,1,function(x){phyper(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"]),
                                                                          length(classlist[[subtype]]),
                                                                          length(x)-length(classlist[[subtype]]),
                                                                          as.numeric(table(as.character(x))["1"]),lower.tail = FALSE)})
  test[,paste("Sample",subtype)] <- apply(mutation1,1,function(x){paste(ifelse(is.na(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                               0,as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                        length(classlist[[subtype]]),
                                                                        length(x)-length(classlist[[subtype]]),
                                                                        ifelse(is.na(as.numeric(table(as.character(x))["1"])),
                                                                               0,as.numeric(table(as.character(x))["1"])),
                                                                        sep = ";")})
}
write.table(test,file = "F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper.txt",col.names = T,row.names = T,quote = F,sep = "\t")

##3、mu突变数据筛选突变率>= 0.1 & mutation$phyper_pvalue <= 0.05\0.01  样本数>=10
rm(list = ls())
library(reshape2)
library(stringr)
mutationselect <- data.frame()
#mutationall <- data.frame()
mutation <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper.txt",sep = "\t",header = T)
mutation$Gene <- rownames(mutation)
mutation$Cancer <- "All"
mutation$Method <- "manhattan"
mutation$Group <- as.character(4)
mutation <- melt(mutation,variable.name = "Subgroup",value.name = "phyper_pvalue")
mutation$Sample <- apply(mutation,1,function(x){as.character(x[grep(paste0(str_sub(x[grep("Subgroup",colnames(mutation))],-1,-1),"$"),colnames(mutation))])})
mutation <- mutation[,-grep("Sample\\.",colnames(mutation))]
mutation <- na.omit(mutation)
mutation$ratio <- unlist(lapply(str_split(mutation$Sample,";"),function(x){ifelse(as.numeric(x[2]) >= 10 & as.numeric(x[3]) >= 10,
                                                                                  as.numeric(x[1])/as.numeric(x[2]),
                                                                                  NA)}))
mutation <- na.omit(mutation)
#mutationselect <- rbind(mutationselect,mutation[mutation$ratio >= 0.1 & mutation$phyper_pvalue <= 0.05,])
#mutationselect <- rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.05,])
#mutationselect <- rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.01,])
a=rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.05,])
b=rbind(mutationselect,mutation[mutation$phyper_pvalue <= 0.01,])
d=rbind(mutationselect,mutation)
write.table(a,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_select_0.05.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(b,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_select_0.01.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(d,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_select_All_3009.txt",col.names = T,row.names = F,quote = F,sep = "\t")


##4、mu画图
## 引用分组时每组样本个数大于10，秩和检验<0.05\0.01基因、频率为这个基因在此亚型中突变的频率
rm(list = ls())
library(ggplot2)
library(grid)
# 突变谱及亚型分类
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
class=data.frame(manhattan.4=f[,"V2"])
rownames(class)=as.character(f[,"V1"])
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu.txt",header=T,sep="\t",row.names = 1)
# 删选后的基因
gene <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_select_0.05.txt",header = T,stringsAsFactors = F)
gene <- gene[unlist(lapply(lapply(lapply(strsplit(as.character(gene$Sample),";"),as.numeric),function(i){i >= 10}),all)),]
gene <- gene[!duplicated(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  genelist1 <- split(genelist[[i]],genelist[[i]]$Subgroup)
  gene1 <- lapply(genelist1,function(x){x$Gene})
  sample <- lapply(genelist1,function(x){rownames(class)[class[,gsub("-",".",i)] %in% gsub("Subtype.","",unique(x$Subgroup))]})
  for(j in names(genelist1)){
    genelist1[[j]]$Freq <- apply(mutation[gene1[[j]],sample[[j]]],1,sum)/length(mutation[gene1[[j]],sample[[j]]])
  }
  genelist[[i]] <- unsplit(genelist1,genelist[[i]]$Subgroup)
  
  p1 <- ggplot(data = genelist[[i]],aes(x = Gene,y = Freq)) +
    geom_bar(fill = "steelblue",stat = "identity") +
    geom_bar(color = "steelblue",fill = "transparent",stat = "identity",position = "fill",size = 1) +
    coord_flip()  +
    scale_y_continuous(labels = scales::percent) +
    #geom_text(aes(label = paste0(signif(Freq,2)*100,"%")), position = position_fill(vjust = 0.5)) +
    facet_grid(~Subgroup) + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=rel(1.1),face="bold"),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  
  p2 <- ggplot(data = genelist[[i]],aes(x = Freq)) +
    geom_density(aes(fill = Subgroup),alpha = 0.8)+
    facet_grid(~Subgroup) + 
    scale_fill_manual(values = c("#3A94BB","#30B89D","#A7C46A","#F9B042","#E05B3E","#9D5F74","#ECC9C9","#E88080","#E8D380","#D87A80")[1:length(unique(genelist[[i]]$Subgroup))]) +
    scale_x_continuous(limits = c(0,1)) +
    labs(title = i) +
    theme(legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())
  
  ########新建画图页面###########
  #grid.newpage()  ##新建页面
  pdf("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_figure_0.05.pdf",width = 15,height = 15)
  pushViewport(viewport(layout = grid.layout(9,1))) ####将页面分成2*2矩阵
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(p1, vp = vplayout(2:9,1))   ###将（1,1)和(1,2)的位置画图c
  print(p2, vp = vplayout(1,1))   ###将(2,1)的位置画图b
  dev.off() ##画下一幅图，记得关闭窗口
}
#gene <- unsplit(genelist,paste(gene$Method,gene$Group,sep = "."))
#write.table(gene,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_p1.txt",col.names = T,row.names = F,sep = "\t")


##突变1075和3009整合为一张图
library(reshape2)
library(stringr)
library(ggplot2)
library(grid)
f1<- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/1075/All_mu_phyper_1075.txt",sep = "\t",header = T,stringsAsFactors = F)
colnames(f1)=paste0(colnames(f1),"_1075")
f2<- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_3009.txt",sep = "\t",header = T,stringsAsFactors = F)
colnames(f2)=paste0(colnames(f2),"_3009")
mutation=cbind(f1,f2)
mutationselect <- data.frame()
mutation$Gene <- rownames(mutation)
mutation$Cancer <- "All"
mutation$Method <- "manhattan"
mutation$Group <- as.character(4)
mutation <- melt(mutation,variable.name = "Subgroup",value.name = "phyper_pvalue")
mutation$Sample <- apply(mutation,1,function(x){as.character(x[grep(paste0(str_sub(as.character(x[grep("Subgroup",colnames(mutation))]),-6,-1),"$"),colnames(mutation))])})
mutation <- mutation[,-grep("Sample\\.",colnames(mutation))]
mutation <- na.omit(mutation)
mutation$ratio <- unlist(lapply(str_split(mutation$Sample,";"),function(x){ifelse(as.numeric(x[2]) >= 10 & as.numeric(x[3]) >= 10,
                                                                                  as.numeric(x[1])/as.numeric(x[2]),
                                                                                  NA)}))
mutation <- na.omit(mutation)
d=rbind(mutationselect,mutation)
write.table(d,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_integrate.txt",col.names = T,row.names = F,quote = F,sep = "\t")
rownames(d)=paste(d$Gene,d$Subgroup,sep="_")
ff=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/mutation_significant.txt",header = T,stringsAsFactors = F)
d1=na.omit(d[paste(ff$Gene,ff$Subgroup,sep="_"),])
write.table(d1,"F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_integrate_mutation_significant.txt",col.names = T,row.names = F,quote = F,sep = "\t")



f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
class=data.frame(manhattan.4=f[,"V2"])
rownames(class)=as.character(f[,"V1"])
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_3009.txt",header=T,sep="\t",row.names = 1)
gene <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_integrate_mutation_significant.txt",header = T,stringsAsFactors = F)
#gene <- gene[unlist(lapply(lapply(lapply(strsplit(as.character(gene$Sample),";"),as.numeric),function(i){i >= 10}),all)),]
gene <- gene[!duplicated(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  genelist1 <- split(genelist[[i]],genelist[[i]]$Subgroup)
  gene1 <- lapply(genelist1,function(x){x$Gene})
  sample <- lapply(genelist1,function(x){rownames(class)[class[,gsub("-",".",i)] %in% gsub("_.*","",gsub("Subtype.","",unique(x$Subgroup)))]})
  for(j in names(genelist1)){
    genelist1[[j]]$Freq <- apply(mutation[gene1[[j]],sample[[j]]],1,sum)/length(mutation[gene1[[j]],sample[[j]]])
  }
  genelist[[i]] <- unsplit(genelist1,genelist[[i]]$Subgroup)
  
  p1 <- ggplot(data = genelist[[i]],aes(x = Gene,y = Freq)) +
    geom_bar(fill = "steelblue",stat = "identity") +
    geom_bar(color = "steelblue",fill = "transparent",stat = "identity",position = "fill",size = 1) +
    coord_flip()  +
    scale_y_continuous(labels = scales::percent) +
    #geom_text(aes(label = paste0(signif(Freq,2)*100,"%")), position = position_fill(vjust = 0.5)) +
    facet_grid(~Subgroup) + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=rel(1.1),face="bold"),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  
  p2 <- ggplot(data = genelist[[i]],aes(x = Freq)) +
    geom_density(aes(fill = Subgroup),alpha = 0.8)+
    facet_grid(~Subgroup) + 
    scale_fill_manual(values = c("deeppink3","deeppink3","goldenrod3","goldenrod3","palegreen4","palegreen4","royalblue4","royalblue4")[1:length(unique(genelist[[i]]$Subgroup))]) +
    scale_x_continuous(limits = c(0,1)) +
    labs(title = i) +
    theme(legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())
  
  ########新建画图页面###########
  #grid.newpage()  ##新建页面
  pdf("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/mutation_significant.pdf",width = 15,height = 15)
  pushViewport(viewport(layout = grid.layout(9,1))) ####将页面分成2*2矩阵
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(p1, vp = vplayout(2:9,1))   ###将（1,1)和(1,2)的位置画图c
  print(p2, vp = vplayout(1,1))   ###将(2,1)的位置画图b
  dev.off() ##画下一幅图，记得关闭窗口
}


#############只取1075
library(ggplot2)
library(grid)
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
class=data.frame(manhattan.4=f[,"V2"])
rownames(class)=as.character(f[,"V1"])
mutation=read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_3009.txt",header=T,sep="\t",row.names = 1)
gene <- read.table("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/All_mu_phyper_integrate_mutation_significant.txt",header = T,stringsAsFactors = F)
gene=gene[grep(".*_1075$",as.character(gene$Subgroup)),]
gene <- gene[!duplicated(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  genelist1 <- split(genelist[[i]],genelist[[i]]$Subgroup)
  gene1 <- lapply(genelist1,function(x){x$Gene})
  sample <- lapply(genelist1,function(x){rownames(class)[class[,gsub("-",".",i)] %in% gsub("_.*","",gsub("Subtype.","",unique(x$Subgroup)))]})
  for(j in names(genelist1)){
    genelist1[[j]]$Freq <- apply(mutation[gene1[[j]],sample[[j]]],1,sum)/length(mutation[gene1[[j]],sample[[j]]])
  }
  genelist[[i]] <- unsplit(genelist1,genelist[[i]]$Subgroup)
  
  p1 <- ggplot(data = genelist[[i]],aes(x = Gene,y = Freq)) +
    geom_bar(fill = "steelblue",stat = "identity") +
    geom_bar(color = "steelblue",fill = "transparent",stat = "identity",position = "fill",size = 1) +
    coord_flip()  +
    scale_y_continuous(labels = scales::percent) +
    #geom_text(aes(label = paste0(signif(Freq,2)*100,"%")), position = position_fill(vjust = 0.5)) +
    facet_grid(~Subgroup) + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=rel(1.1),face="bold"),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  
  p2 <- ggplot(data = genelist[[i]],aes(x = Freq)) +
    geom_density(aes(fill = Subgroup),alpha = 0.8)+
    facet_grid(~Subgroup) + 
    scale_fill_manual(values = c("deeppink3","goldenrod3","palegreen4","royalblue4")[1:length(unique(genelist[[i]]$Subgroup))]) +
    scale_x_continuous(limits = c(0,1)) +
    labs(title = i) +
    theme(legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())
  
  ########新建画图页面###########
  #grid.newpage()  ##新建页面
  pdf("F:/tp53/CBioportal_TP53/TCGAmutationdata/Mutation/3009/mutation_significant1.pdf",width = 15,height = 15)
  pushViewport(viewport(layout = grid.layout(9,1))) ####将页面分成2*2矩阵
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(p1, vp = vplayout(2:9,1))   ###将（1,1)和(1,2)的位置画图c
  print(p2, vp = vplayout(1,1))   ###将(2,1)的位置画图b
  dev.off() ##画下一幅图，记得关闭窗口
}





################################################
##Allcancer2_class4_manhattan_ward.D_1\3009
##Allcancer2_class4_manhattan_ward.D_1075_1\1075
##1、cnv对拷贝数变异数据进行处理整合成Allcancer的全基因组的cnv信息 
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("F:/tp53/CNV/TCGA hub/",i,"_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("Gene.Symbol",intersect(f4,colnames(f3)))])))
  f1=f5
}
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  f3 <- read.table(paste0("F:/tp53/CNV/TCGA hub/",i,"_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("Gene.Symbol",intersect(f4,colnames(f3)))])))
  #f1=merge(f1,f5,by.x="Gene.Symbol",by.y="Gene.Symbol")
  f1=merge(f1,f5,all=T)
}
ff=unique(f1[-grep("^TP53$",as.character(f1[,1])),])
ff[sapply(ff,is.na)]<-0
write.table(ff,file ="F:/tp53/CNV/1075/All_cnv.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)


##2、cnv用超几何检验
rm(list = ls())  #删去上面的程序
class <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
ff=read.table("F:/tp53/CNV/3009/All_cnv.txt",header=T,sep="\t")
rownames(ff)=as.character(ff[,1])
ff1 <- ff[,intersect(colnames(ff),class$V1)]
ff1[ff1 != 0] <- 1
class <- class[class$V1 %in% colnames(ff1),]
classlist <- split(class$V1,class$V2)
test <- matrix(0,nrow = nrow(ff1),ncol = 4*2,dimnames = list(rownames(ff1),c(paste("Subtype",c(1:4)),paste("Sample",c(1:4)))))
for (subtype in names(classlist)){
  test[,paste("Subtype",subtype)] <- apply(ff1,1,function(x){phyper(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"]),
                                                                          length(classlist[[subtype]]),
                                                                          length(x)-length(classlist[[subtype]]),
                                                                          as.numeric(table(as.character(x))["1"]),lower.tail = FALSE)})
  test[,paste("Sample",subtype)] <- apply(ff1,1,function(x){paste(ifelse(is.na(as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                               0,as.numeric(table(as.character(x[classlist[[subtype]]]))["1"])),
                                                                        length(classlist[[subtype]]),
                                                                        length(x)-length(classlist[[subtype]]),
                                                                        ifelse(is.na(as.numeric(table(as.character(x))["1"])),
                                                                               0,as.numeric(table(as.character(x))["1"])),
                                                                        sep = ";")})
}
write.table(test,file = "F:/tp53/CNV/1075/All_cnv_phyper.txt",col.names = T,row.names = T,quote = F,sep = "\t")

##3、cnv数据筛选
rm(list = ls())
library(reshape2)
library(stringr)
mutationselect <- data.frame()
mutationselect30 <- data.frame()
#mutationselect40 <- data.frame()
#mutationselect50 <- data.frame()
mutationall <- data.frame() 
mutation <- read.table("F:/tp53/CNV/3009/All_cnv_phyper.txt",sep = "\t",header = T,stringsAsFactors = F)
mutation$Gene <- rownames(mutation)
mutation$Cancer <- "All"
mutation$Method <- "manhattan"
mutation$Group <- as.character(4)
mutation <- melt(mutation,variable.name = "Subgroup",value.name = "phyper_pvalue")
mutation$Sample <- apply(mutation,1,function(x){x[grep(str_sub(x[grep("Subgroup",colnames(mutation))],-1,-1),colnames(mutation))]})
mutation <- mutation[,-grep("Sample\\.",colnames(mutation))]
mutation <- na.omit(mutation)
mutation$ratio <- unlist(lapply(str_split(mutation$Sample,";"),function(x){ifelse(as.numeric(x[2]) >= 10 & as.numeric(x[3]) >= 10,
                                                                                  as.numeric(x[1])/as.numeric(x[2]),
                                                                                  NA)}))
mutation <- na.omit(mutation)
mutationselect30 <- rbind(mutationselect,mutation[mutation$ratio >= 0.3 & mutation$phyper_pvalue <= 0.05,])
#mutationselect40 <- rbind(mutationselect,mutation[mutation$ratio >= 0.4 & mutation$phyper_pvalue <= 0.05,])
#mutationselect50 <- rbind(mutationselect,mutation[mutation$ratio >= 0.5 & mutation$phyper_pvalue <= 0.05,])
mutationall <- rbind(mutationall,mutation)
write.table(mutationall,file = "F:/tp53/CNV/3009/All_cnv_phyper_selectall.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(mutationselect30,file ="F:/tp53/CNV/3009/All_cnv_phyper_select30.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(mutationselect40,file ="F:/tp53/CNV/3009/All_cnv_phyper_select40.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(mutationselect50,file ="F:/tp53/CNV/3009/All_cnv_phyper_select50.txt",col.names = T,row.names = F,quote = F,sep = "\t")

##4、cnv画图
## 引用分组时ratio>0.3，秩和检验<0.01基因
## 频率为这个基因在此亚型中突变的频率
rm(list = ls())
library(ggplot2)
library(grid)
library(reshape2) 
color <- data.frame(row.names = c("Censored2","Censored1","nomutate","Amplification1","Amplification2","Censored","Amplification"),color = c("#45315D","#665C84","#FaF9F9","#FF9280","#FF2400","#45315D","#FF2400"),stringsAsFactors = F)
# 突变谱及亚型分类
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
class=data.frame(manhattan.4=f[,"V2"])
rownames(class)=as.character(f[,"V1"])
mutation=read.table("F:/tp53/CNV/3009/All_cnv.txt",header=T,sep="\t",row.names = 1)
# 删选后的基因
gene <- read.table("F:/tp53/CNV/3009/All_cnv_phyper_select30.txt",sep = "\t",header = T,stringsAsFactors = F)
gene <- gene[!duplicated(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = ".")) 
for (i in names(genelist)){
  genelist1 <- split(genelist[[i]],genelist[[i]]$Subgroup)
  gene1 <- lapply(genelist1,function(x){x$Gene})
  sample <- lapply(genelist1,function(x){rownames(class)[class[,gsub("-",".",i)] %in% gsub("Subtype.","",unique(x$Subgroup))]})
  for(j in names(genelist1)){
    genelist1[[j]]$Censored2 <- apply(mutation[gene1[[j]],intersect(colnames(mutation),sample[[j]])],1,function(x){length(which((x == -2)))})/length(intersect(colnames(mutation),sample[[j]]))
    genelist1[[j]]$Censored1 <- apply(mutation[gene1[[j]],intersect(colnames(mutation),sample[[j]])],1,function(x){length(which((x == -1)))})/length(intersect(colnames(mutation),sample[[j]]))
    genelist1[[j]]$nomutate <- apply(mutation[gene1[[j]],intersect(colnames(mutation),sample[[j]])],1,function(x){length(which((x == 0)))})/length(intersect(colnames(mutation),sample[[j]]))
    genelist1[[j]]$Amplification1 <- apply(mutation[gene1[[j]],intersect(colnames(mutation),sample[[j]])],1,function(x){length(which((x == 1)))})/length(intersect(colnames(mutation),sample[[j]]))
    genelist1[[j]]$Amplification2 <- apply(mutation[gene1[[j]],intersect(colnames(mutation),sample[[j]])],1,function(x){length(which((x == 2)))})/length(intersect(colnames(mutation),sample[[j]]))
  }
  genelist[[i]] <- unsplit(genelist1,genelist[[i]]$Subgroup)
  genelist[[i]] <- melt(genelist[[i]],id.vars = c("Gene","Cancer","Method","Group","Subgroup","phyper_pvalue","Sample","ratio"),value.name = "mutation_ratio",variable.name = "mutation_method")
  genelist[[i]]$mutation_method1 <- gsub("[0-9]","",genelist[[i]]$mutation_method)
  genelist[[i]]$mutation_method1[genelist[[i]]$mutation_method1 == "nomutate"] <- NA
  genelist[[i]]$mutation_method1 <- factor(genelist[[i]]$mutation_method1,levels = unique(genelist[[i]]$mutation_method1))
 
  p1 <- ggplot(data = genelist[[i]],aes(x = Gene,y = mutation_ratio)) +
    geom_bar(color = "grey30",fill = "transparent",stat = "identity",position = "fill") +
    geom_bar(aes(fill = mutation_method),stat = "identity") +
    scale_fill_manual(values = color[levels(genelist[[i]]$mutation_method),]) +
    coord_flip()  +
    scale_y_continuous(labels = scales::percent) +
    #geom_text(aes(label = paste0(signif(Freq,2)*100,"%")), position = position_fill(vjust = 0.5)) +
    facet_grid(~Subgroup,scales = "free_x") + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=0.1,face="bold"),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  
  p2 <- ggplot(data = genelist[[i]],aes(x = mutation_ratio)) +
    geom_density(aes(color = mutation_method1))+
    facet_wrap(~Subgroup,scales = "free_y",nrow = 1) + 
    scale_color_manual(values = color[levels(genelist[[i]]$mutation_method1),]) +
    scale_x_continuous(limits = c(0,1)) +
    labs(title = i) +
    theme(legend.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())
  ########新建画图页面###########
  #grid.newpage()  ##新建页面
  pdf("F:/tp53/CNV/3009/All_cnv_figure.pdf",width = 15,height = 15)
  pushViewport(viewport(layout = grid.layout(9,1))) ####将页面分成2*2矩阵
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(p1, vp = vplayout(2:9,1))   ###将（1,1)和(1,2)的位置画图c
  print(p2, vp = vplayout(1,1))   ###将(2,1)的位置画图b
  dev.off() ##画下一幅图，记得关闭窗口
}    
    

###########################################################
##Allcancer2_class4_manhattan_ward.D_1\3009
##Allcancer2_class4_manhattan_ward.D_1075_1\1075
##1、对mRNA-FPKM进行处理整合成Allcancer对应的全基因组的mRNA信息
##文件太大，分段，数据调用思宇的mRNA的F:\tp53\mRNA\data下的pancancer-1.txt共20个分段
rm(list = ls())
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("D:/tp53/Gene Expression/Gene_FPKM1/",i,"_gene_fpkm.txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("gene_Name",f4)])))
  f1=f5
}
f11=f1[1:5000,]##文件太大分段
#f12=f1[15001:30000,]
#f13=f1[30001:45000,]
#f14=f1[45001:nrow(f1),]
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  f3 <- read.table(paste0("D:/tp53/Gene Expression/Gene_FPKM1/",i,"_gene_fpkm.txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("gene_Name",f4)])))
  f11=merge(f11,f5,by.x="gene_Name",by.y="gene_Name")
}
#write.table(unique(na.omit(f11)),file ="F:/tp53/mRNA/All_mRNA_part1.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)

##对每个分段进行检验，然后rbind为一个文件
##2、mRNA\wilcoxon test && fold change、wilcoxon p value < 0.05;fold change > 1 & < -1
library(reshape2)
library(stringr)
for(i in 1:20){
  gene <- read.table(paste("F:/tp53/mRNA/data/pancancer-",i,".txt",sep=""),row.names = 1,sep = "\t",header = T)
  colnames(gene)=gsub("-",".",substr(colnames(gene),1,12))
  gene <- gene[rowSums(gene) > 0,]
  f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
  classall=data.frame(manhattan.4=f[,"V2"])
  rownames(classall)=as.character(f[,"V1"])
  testall <- data.frame() 
  for (method in colnames(classall)){
    tryCatch({
      gene1 <- gene[,intersect(colnames(gene),rownames(classall))]
      class <- data.frame(gene = colnames(gene1),method = classall[colnames(gene1),method])
      classlist <- split(class$gene,class$method)
      test <- matrix(0,nrow = nrow(gene1),ncol = as.numeric(gsub("^[a-z\\.]+","",method))*2,dimnames = list(rownames(gene1),c(paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method))),paste("FC",c(1:gsub("^[a-z\\.]+","",method))))))
      for (subtype in names(classlist)){
        test[,paste("Wilcoxon",subtype)] <- apply(gene1,1,function(x){wilcox.test(x[factor(classlist[[subtype]],levels = colnames(gene1))],x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])$p.value})
        test[,paste("FC",subtype)] <- apply(gene1,1,function(x){mean(x[factor(classlist[[subtype]],levels = colnames(gene1))])/mean(x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])})
      }
      test <- data.frame(Gene = rownames(test),test)
      testwilcoxon <- melt(test[,c("Gene",paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "Wilcoxon",value.name = "Wilpvalue")
      testwilcoxon$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testwilcoxon$Wilcoxon)))
      testfc <- melt(test[,c("Gene",paste("FC",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "FC",value.name = "FCp")
      testfc$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testfc$FC)))
      test <- merge(testwilcoxon[,c(1,3,4)],testfc[,c(1,3,4)],by = c("Gene","Subgroup"),all = T)
      test$Cancer <- "All"
      test$Method <- gsub("\\.","-",gsub("[^a-z]+$","",method))
      test$Group <- gsub("^[a-z\\.]+","",method)
      test$Sample <- str_c(unlist(lapply(classlist,length)),collapse = ";")
      # 条件
      test <- test[!is.nan(test$Wilpvalue) & !is.nan(test$FCp) & test$Wilpvalue < 0.05,]
      test$FCp[is.infinite(test$FCp)] <- 10
      test <- test[test$FCp > 2 | test$FCp < 0.5,]
      testall <- rbind(testall,test)
    },
    error = function(e){cat("ERROR :",method,"  ",conditionMessage(e),"\n")})
  }
  write.table(testall,file = paste("F:/tp53/mRNA/3009/pancancer_",i,"_Wil_FC.txt",sep=""),col.names = T,row.names = F,quote = F,sep = "\t")
}

#合并文件
ff=data.frame()
for(i in 1:20){
  f <- read.table(paste("F:/tp53/mRNA/3009/pancancer_",i,"_Wil_FC.txt",sep=""),sep = "\t",header = T)
  ff=rbind(ff,f)
}
ff[,1]=as.character(ff[,1])
write.table(ff,file ="F:/tp53/mRNA/3009/All_mRNA_Wil_FC0.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#f1=read.table("D:/tp53/Gene Expression/ENSG_Name2.txt",sep = "\t",header = F)
#f1[,1]=gsub(" ","",as.character(f1[,1]))
#f1[,2]=as.character(f1[,2])
#ff1=merge(ff,f1,by.x="Gene",by.y="V1",all.x=T)
#colnames(ff1)[1]="V1"
#colnames(ff1)[9]="Gene"
#ff2=ff1[,c("Gene","Subgroup","Wilpvalue","FCp","Cancer","Method","Group","Sample")]
#write.table(ff2,file ="F:/tp53/mRNA/1075/All_mRNA_Wil_FC.txt",col.names = T,row.names = F,quote = F,sep = "\t")


##3、mRNA画图 ratio>0.1
rm(list = ls())
library(ggplot2)
library(reshape2)
gene1 <- read.table("F:/tp53/mRNA/1075/All_mRNA_Wil_FC0.txt",sep = "\t",header = T,stringsAsFactors = F)
#gene1 <-gene1[(gene1$Wilpvalue<0.00001)&(gene1$FCp > 6 | gene1$FCp < (1/6)),]
gene1 <-gene1[(gene1$Wilpvalue<0.001)&(gene1$FCp > 3 | gene1$FCp < (1/3)),]  #DAVID
f1=read.table("D:/tp53/Gene Expression/ENSG_Name2.txt",sep = "\t",header = F)
f1[,1]=gsub(" ","",as.character(f1[,1]))
f1[,2]=as.character(f1[,2])
ff1=merge(gene1,f1,by.x="Gene",by.y="V1",all.x=T)
colnames(ff1)[1]="V1"
colnames(ff1)[9]="Gene"
ff2=ff1[,c("Gene","Subgroup","Wilpvalue","FCp","Cancer","Method","Group","Sample")]
ge=data.frame(table(ff2[,1]))
gene <- ff2[(ff2[,1] %in% as.character(ge[which(ge[,2]<=4),1])),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  a <- genelist[[i]]
  a$FCp <- log2(a$FCp)
  a$FCp[a$FCp == -Inf] <- -max(abs(a$FCp[a$FCp != -Inf]))
  dis <- dcast(a[,c(1,2,4)],Gene ~ Subgroup,value.var = "FCp")
  dis[is.na(dis)] <- 0
  hc <- hclust(dist(dis[,-1]),method = "average")
  dis <- dis[hc$order,]
  geneclust <- unique(dis$Gene)
  dis <- melt(dis,id.vars = "Gene",variable.name = "Subgroup",value.name = "log2(FC)")
  a <- merge(dis,a,all.x = T)
  a$Wilpvalue[is.na(a$Wilpvalue)] <- 0
  a$Gene <- factor(a$Gene,levels = geneclust)
  p <- ggplot(a,aes(Gene,-log10(Wilpvalue)))
  p + #geom_tile(aes(fill = FCp)) +
    #geom_bar(color = "grey30",fill = "transparent",stat = "identity",position = "fill") +
    geom_bar(aes(fill = `log2(FC)`),stat = "identity") +
    scale_fill_gradient2(low = "#05004E",high = "#FF0000",mid = "#F9F9F9") +
    coord_flip()  +
    facet_grid(~Subgroup,scales = "free_x") + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=10,face="bold"),
      axis.text.y = element_text(size = 1),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  #ggsave("F:/tp53/mRNA/1075/All_mRNA_figure_0.00001_6_1.pdf",width = 10,height = 10)
}
write.table(a,file ="F:/tp53/mRNA/1075/mRNA.txt",col.names = T,row.names = F,quote = F,sep = "\t")




###########################################################
##Allcancer2_class4_manhattan_ward.D_1\3009
##Allcancer2_class4_manhattan_ward.D_1075_1\1075
##1、对miRNA-RPKM进行处理整合成Allcancer对应的全基因组的miRNA信息
rm(list = ls())
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("F:/tp53/miRNA/miRNA_fpkm/",i,"_miRNA_fpkm.txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("miRNA_ID",intersect(f4,colnames(f3)))])))
  f1=f5
}
colnames(f1)="miRNA_ID"
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  f3 <- read.table(paste0("F:/tp53/miRNA/miRNA_fpkm/",i,"_miRNA_fpkm.txt"),sep = "\t",header = T)
  colnames(f3)=gsub("-",".",substr(colnames(f3),1,12))
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("miRNA_ID",intersect(f4,colnames(f3)))])))
  f1=merge(f1,f5,all=T)
}
f1[sapply(f1,is.na)]<-0
write.table(f1,file ="F:/tp53/miRNA/3009/All_miRNA.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)

#2、miRNA数据处理
# wilcoxon test && fold change、wilcoxon p value < 0.01;fold change > 1 & < -1
library(reshape2)
library(stringr)
gene <- read.table("F:/tp53/miRNA/3009/All_miRNA.txt",row.names = 1,sep = "\t",header = T)
gene <- gene[rowSums(gene) > 0,]
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
classall=data.frame(manhattan.4=f[,"V2"])
rownames(classall)=as.character(f[,"V1"])
testall <- data.frame()
for (method in colnames(classall)){
  tryCatch({
    gene1 <- gene[,intersect(colnames(gene),rownames(classall))]
    class <- data.frame(gene = colnames(gene1),method = classall[colnames(gene1),method])
    classlist <- split(class$gene,class$method)
    test <- matrix(NA,nrow = nrow(gene1),ncol = as.numeric(gsub("^[a-z\\.]+","",method))*2,dimnames = list(rownames(gene1),c(paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method))),paste("FC",c(1:gsub("^[a-z\\.]+","",method))))))
    for (subtype in names(classlist)){
      test[,paste("Wilcoxon",subtype)] <- apply(gene1,1,function(x){wilcox.test(x[factor(classlist[[subtype]],levels = colnames(gene1))],x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])$p.value})
      test[,paste("FC",subtype)] <- apply(gene1,1,function(x){mean(x[factor(classlist[[subtype]],levels = colnames(gene1))])/mean(x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])})
    }
    test <- data.frame(Gene = rownames(test),test)
    testwilcoxon <- melt(test[,c("Gene",paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "Wilcoxon",value.name = "Wilpvalue")
    testwilcoxon$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testwilcoxon$Wilcoxon)))
    testfc <- melt(test[,c("Gene",paste("FC",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "FC",value.name = "FCp")
    testfc$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testfc$FC)))
    test <- merge(testwilcoxon[,c(1,3,4)],testfc[,c(1,3,4)],by = c("Gene","Subgroup"),all = T)
    test$Cancer <- "All"
    test$Method <- gsub("\\.","-",gsub("[^a-z]+$","",method))
    test$Group <- gsub("^[a-z\\.]+","",method)
    test$Sample <- str_c(unlist(lapply(classlist,length)),collapse = ";")
    # 条件
    test <- test[!is.nan(test$Wilpvalue) & !is.nan(test$FCp) & test$Wilpvalue < 0.01,]
    test$FCp[is.infinite(test$FCp)] <- 10
    test <- test[test$FCp > 2 | test$FCp < 0.5,]
    testall <- rbind(testall,test)
    
  },
  error = function(e){cat("ERROR :",method,"  ",conditionMessage(e),"\n")})
}
write.table(testall,file ="F:/tp53/miRNA/3009/All_miRNA_Wil_FC.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#3、miRNA画图
library(ggplot2)
library(reshape2)
gene <- read.table("F:/tp53/miRNA/1075/All_miRNA_Wil_FC.txt",sep = "\t",header = T,stringsAsFactors = F)
gene <-gene[(gene$Wilpvalue<0.00001)&(gene$FCp > 6 | gene$FCp < (1/6)),]
gene <- gene[complete.cases(gene),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  a <- genelist[[i]]
  a$FCp <- log2(a$FCp)
  a$FCp[a$FCp == -Inf] <- -max(abs(a$FCp[a$FCp != -Inf]))
  dis <- dcast(a[,c(1,2,4)],Gene ~ Subgroup,value.var = "FCp")
  dis[is.na(dis)] <- 0
  hc <- hclust(dist(dis[,-1]),method = "average")
  dis <- dis[hc$order,]
  geneclust <- unique(dis$Gene)
  dis <- melt(dis,id.vars = "Gene",variable.name = "Subgroup",value.name = "log2(FC)")
  a <- merge(dis,a,all.x = T)
  a$Wilpvalue[is.na(a$Wilpvalue)] <- 0
  a$Gene <- factor(a$Gene,levels = geneclust)
  p <- ggplot(a,aes(Gene,-log10(Wilpvalue)))
  p + #geom_tile(aes(fill = FCp)) +
    #geom_bar(color = "grey30",fill = "transparent",stat = "identity",position = "fill") +
    geom_bar(aes(fill = `log2(FC)`),stat = "identity") +
    scale_fill_gradient2(low = "#05004E",high = "#FF0000",mid = "#F9F9F9") +
    coord_flip()  +
    facet_grid(~Subgroup,scales = "free_x") + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=10,face="bold"),
      axis.text.y = element_text(size = 8),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  #ggsave("F:/tp53/miRNA/1075/All_miRNA_figure_0.00001_6.pdf",width = 10,height = 10)
}
write.table(a,file ="F:/tp53/miRNA/1075/38_miRNA.txt",col.names = T,row.names = F,quote = F,sep = "\t")




##############
f <- read.table("F:/tp53/miRNA/1075/38_miRNA.txt",sep = "\t",header = T,stringsAsFactors = F)
f1 <- read.table("F:/tp53/mRNA/starBase-miRNA_mRNA.txt",sep = "\t",header = T,stringsAsFactors = F)
f1[,1]=gsub("miR","mir",f1[,1])
for( i in 1:4){
  a=f[which((f[,2]==paste("Subtype",i,sep=" ")) & (f[,3]!=0)),1]
  for(j in 1:length(a)){
    a1=f1[grep(a[j],f1[,1]),2]
    write.table(a1,file =paste("F:/tp53/miRNA/1075/C",i,".txt",sep=""),col.names = F,row.names = F,quote = F,sep = "\t",append=T)
  }
}






###########################################################
##Allcancer2_class4_manhattan_ward.D_1\3009
##Allcancer2_class4_manhattan_ward.D_1075_1\1075
##1、对lncRNA-RPkM进行处理整合成Allcancer对应的全基因组的lncRNA信息
rm(list = ls())
library(stringr)
f2<-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
for(i in unique(as.character(f2[,3]))[1]){
  f3 <- read.table(paste0("F:/tp53/lncRNA/TCGA-",i,"-rnaexpr/TCGA-",i,"-rnaexpr.tsv"),sep = "\t",header = T)
  colnames(f3)=str_sub(colnames(f3),-18,-1)
  colnames(f3)[1]="Gene_ID"
  f4=as.character(f2[as.character(f2[,3])==i,1])
  f5=unique(na.omit(data.frame(f3[,c("Gene_ID",intersect(paste0("Tumor.",f4),colnames(f3)))])))
  f1=f5
}
for(i in unique(as.character(f2[,3]))[2:length(unique(as.character(f2[,3])))]){
  if(dir.exists(paste0("F:/tp53/lncRNA/TCGA-",i,"-rnaexpr/"))){
    f3 <- read.table(paste0("F:/tp53/lncRNA/TCGA-",i,"-rnaexpr/TCGA-",i,"-rnaexpr.tsv"),sep = "\t",header = T)
    colnames(f3)=str_sub(colnames(f3),-18,-1)
    colnames(f3)[1]="Gene_ID"
    f4=as.character(f2[as.character(f2[,3])==i,1])
    f5=unique(na.omit(data.frame(f3[,c("Gene_ID",intersect(paste0("Tumor.",f4),colnames(f3)))])))
    if(ncol(f5)==1){
      colnames(f5)[1]="Gene_ID"
    }
    f1=merge(f1,f5,all=T)
  }
}
f1[sapply(f1,is.na)]<-0
colnames(f1)=gsub("Tumor.","",colnames(f1))
write.table(f1,file ="F:/tp53/lncRNA/3009/All_lncRNA.txt", row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)



#2、lncRNA数据处理 1075
# wilcoxon test && fold change、wilcoxon p value < 0.01;fold change > 1 & < -1
rm(list = ls())
library(reshape2)
library(stringr)
gene <-read.table("F:/tp53/lncRNA/3009/All_lncRNA.txt",row.names = 1,sep = "\t",header = T)
gene <- gene[rowSums(gene) > 0,]
f <-read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1.txt",header=F,sep="\t")
classall=data.frame(manhattan.4=f[,"V2"])
rownames(classall)=as.character(f[,"V1"])
testall <- data.frame()
for (method in colnames(classall)){
  tryCatch({
    gene1 <- gene[,intersect(colnames(gene),rownames(classall))]
    class <- data.frame(gene = colnames(gene1),method = classall[colnames(gene1),method])
    classlist <- split(class$gene,class$method)
    test <- matrix(0,nrow = nrow(gene1),ncol = as.numeric(gsub("^[a-z\\.]+","",method))*2,dimnames = list(rownames(gene1),c(paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method))),paste("FC",c(1:gsub("^[a-z\\.]+","",method))))))
    for (subtype in names(classlist)){
      test[,paste("Wilcoxon",subtype)] <- apply(gene1,1,function(x){wilcox.test(x[factor(classlist[[subtype]],levels = colnames(gene1))],x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])$p.value})
      test[,paste("FC",subtype)] <- apply(gene1,1,function(x){mean(x[factor(classlist[[subtype]],levels = colnames(gene1))])/mean(x[factor(setdiff(colnames(gene1),classlist[[subtype]]),levels = colnames(gene1))])})
    }
    test <- data.frame(Gene = rownames(test),test)
    testwilcoxon <- melt(test[,c("Gene",paste("Wilcoxon",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "Wilcoxon",value.name = "Wilpvalue")
    testwilcoxon$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testwilcoxon$Wilcoxon)))
    testfc <- melt(test[,c("Gene",paste("FC",c(1:gsub("^[a-z\\.]+","",method)),sep = "."))],variable.name = "FC",value.name = "FCp")
    testfc$Subgroup <- paste("Subtype",gsub("^[A-Za-z\\.]+","",as.character(testfc$FC)))
    test <- merge(testwilcoxon[,c(1,3,4)],testfc[,c(1,3,4)],by = c("Gene","Subgroup"),all = T)
    test$Cancer <- "All"
    test$Method <- gsub("\\.","-",gsub("[^a-z]+$","",method))
    test$Group <- gsub("^[a-z\\.]+","",method)
    test$Sample <- str_c(unlist(lapply(classlist,length)),collapse = ";")
    # 条件
    test <- test[!is.nan(test$Wilpvalue) & !is.nan(test$FCp) & test$Wilpvalue < 0.01,]
    test$FCp[is.infinite(test$FCp)] <- 10
    test <- test[test$FCp > 2 | test$FCp < 0.5,]
    testall <- rbind(testall,test)
  },
  error = function(e){cat("ERROR :",method,"  ",conditionMessage(e),"\n")})
}
write.table(testall,file ="F:/tp53/lncRNA/3009/All_lncRNA_Wil_FC.txt",col.names = T,row.names = F,quote = F,sep = "\t")



#3、lncRNA画图
rm(list = ls())
library(ggplot2)
library(reshape2)
gene1 <- read.table("F:/tp53/lncRNA/1075/All_lncRNA_Wil_FC.txt",sep = "\t",header = T,stringsAsFactors = F)
gene1 <-gene1[(gene1$Wilpvalue<0.00001)&(gene1$FCp > 6 | gene1$FCp < (1/6)),]
f1=read.table("D:/tp53/Gene Expression/ENSG_Name2.txt",sep = "\t",header = F)
f1[,1]=gsub(" ","",as.character(f1[,1]))
f1[,2]=as.character(f1[,2])
ff1=merge(gene1,f1,by.x="Gene",by.y="V1",all.x=T)
colnames(ff1)[1]="V1"
colnames(ff1)[9]="Gene"
ff2=na.omit(ff1[,c("Gene","Subgroup","Wilpvalue","FCp","Cancer","Method","Group","Sample")])
ge=data.frame(table(ff2[,1]))
gene <- ff2[(ff2[,1] %in% as.character(ge[which(ge[,2]<=4),1])),]
genelist <- split(gene,paste(gene$Method,gene$Group,sep = "."))
for (i in names(genelist)){
  a <- genelist[[i]]
  a$FCp <- log2(a$FCp)
  a$FCp[a$FCp == -Inf] <- -max(abs(a$FCp[a$FCp != -Inf]))
  dis <- dcast(a[,c(1,2,4)],Gene ~ Subgroup,value.var = "FCp")
  dis[is.na(dis)] <- 0
  hc <- hclust(dist(dis[,-1]),method = "average")
  dis <- dis[hc$order,]
  geneclust <- unique(dis$Gene)
  dis <- melt(dis,id.vars = "Gene",variable.name = "Subgroup",value.name = "log2(FC)")
  a <- merge(dis,a,all.x = T)
  a$Wilpvalue[is.na(a$Wilpvalue)] <- 0
  a$Gene <- factor(a$Gene,levels = geneclust)
  p <- ggplot(a,aes(Gene,-log10(Wilpvalue)))
  p + #geom_tile(aes(fill = FCp)) +
    #geom_bar(color = "grey30",fill = "transparent",stat = "identity",position = "fill") +
    geom_bar(aes(fill = `log2(FC)`),stat = "identity") +
    scale_fill_gradient2(low = "#05004E",high = "#FF0000",mid = "#F9F9F9") +
    coord_flip()  +
    facet_grid(~Subgroup,scales = "free_x") + 
    #labs(title = i) +
    theme(#legend.position = "top",
      legend.background=element_blank(),
      axis.text.x=element_text(size=10,face="bold"),
      axis.text.y = element_text(size = 1),
      axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank())
  #ggsave("F:/tp53/lncRNA/1075/All_lncRNA_figure_0.00001_6_1.pdf",width = 10,height = 10)
}
write.table(a,file ="F:/tp53/lncRNA/1075/447_lncRNA.txt",col.names = T,row.names = F,quote = F,sep = "\t")










#################################
install.packages("Rtsne")
library(Rtsne)
d1075<-read.table("E:/图表/t-SNE data/Allcancer2_manhattan_ward.D_1075.txt",header=T)
d1075<-t(d1075)#转置
d1075<-normalize_input(d1075)#标准化
feature<-read.table("E:/图表/t-SNE data/feature1075.txt",header=F)#feature
set.seed(100)
feature$V2<-as.factor(feature$V2)
colors = col=c("deeppink3","goldenrod3","palegreen4","royalblue4")
names(colors)<-unique(feature$V2)
#excuting
tsne <- Rtsne(d1075,dim=2,perplexity=30,theta=0.0,max_iter=500)
plot(tsne$Y,pch=16,cex=0.5,col=colors[feature$V2],main="tsne-perplexity=30",asp=0)
##text(tsne$Y,labels=feature$V2,col=colors[feature$V2])
savePlot("E:/图表/t-SNE data/1075_1",type="pdf",device=dev.cur(),restoreConsole=TRUE)

##
library(Rtsne)
d3009<-read.table("E:/图表/t-SNE data/Allcancer2_manhattan_ward.D_3009.txt",header=T)
d3009<-t(d3009)#转置
d3009<-normalize_input(d3009)#标准化
feature<-read.table("E:/图表/t-SNE data/3009feature.txt",header=F)#feature
set.seed(100)
feature$V2<-as.factor(feature$V2)
colors =c("deeppink3","goldenrod3","palegreen4","royalblue4")
names(colors)<-unique(feature$V2)
#excuting
tsne <- Rtsne(d3009,dim=2,perplexity=30,theta=0.0,max_iter=500)
plot(tsne$Y,pch=16,cex=0.5,col=colors[feature$V2],main="tsne-perplexity=30",asp=0)
##text(tsne$Y,labels=feature$V2,col=colors[feature$V2])
savePlot("E:/图表/t-SNE data/3009_1",type="pdf",device=dev.cur(),restoreConsole=TRUE)

