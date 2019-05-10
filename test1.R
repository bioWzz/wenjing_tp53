#' 使用R基本绘图函数绘制y轴不连续的柱形图
#' 绘制y轴不连续的柱形图，具有误差线添加功能。断点位置通过btm和top参数设置，如果不设置，函数可自动计算合适的断点位置。
#' @title gap.barplot function
#' @param df 长格式的data.frame，即数据框中每一列为一组绘图数据。
#' @param y.cols 用做柱形图y值的数据列（序号或名称），一列为一组。
#' @param sd.cols 与y值列顺序对应的误差值的数据列（序号或名称）。
#' @param btm 低位断点。如果btm和top均不设置，程序将自动计算和设置断点位置。
#' @param top 高位断点。
#' @param min.range 自动计算断点的阈值：最大值与最小值的最小比值
#' @param max.fold 自动计算断点时最大值与下方数据最大值的最大倍数比
#' @param ratio 断裂后上部与下部y轴长度的比例。
#' @param gap.width y轴断裂位置的相对物理宽度（非坐标轴实际刻度）
#' @param brk.type 断点类型，可设为normal或zigzag
#' @param brk.bg 断点处的背景颜色
#' @param brk.srt 断点标记线旋转角度
#' @param brk.size 断点标记线的大小（长度）
#' @param brk.col 断点标记线的颜色
#' @param brk.lwd 断点标记线的线宽
#' @param cex.error 误差线相对长度，默认为1
#' @param ... 其他传递给R基本绘图函数barplot的参数
#' @return 返回barplot的原始返回值，即柱形图的x坐标
#' @examples
#' datax <- na.omit(airquality)[,1:4]
#' cols <- cm.colors(ncol(datax))
#' layout(matrix(1:6, ncol=2))
#' set.seed(0)
#' for (ndx in 1:6){
#'     dt <- datax[sample(rownames(datax), 10), ]
#'     par(mar=c(0.5,2,0.5,0.5))
#'     brkt <- sample(c('normal', 'zigzag'), 1)
#'     gap.barplot(dt, col=cols, brk.type=brkt, max.fold=5, ratio=2)
#' }

####一、自己定义的函数gap.barplot，先运行一遍
gap.barplot <- function(df, y.cols = 1:ncol(df), sd.cols = NULL, btm = NULL,
                        top = NULL, min.range = 10, max.fold = 5, ratio = 1, gap.width = 1, brk.type = "normal",
                        brk.bg = "white", brk.srt = 135, brk.size = 1, brk.col = "black", brk.lwd = 1,
                        cex.error = 1, ...) {
  if (missing(df))
    stop("No data provided.")
  if (is.numeric(y.cols))
    ycol <- y.cols else ycol <- colnames(df) == y.cols
    if (!is.null(sd.cols))
      if (is.numeric(sd.cols))
        scol <- sd.cols else scol <- colnames(df) == sd.cols
        ## Arrange data
        opts <- options()
        options(warn = -1)
        y <- t(df[, ycol])
        colnames(y) <- NULL
        if (missing(sd.cols))
          sdx <- 0 else sdx <- t(df[, scol])
        sdu <- y + sdx
        sdd <- y - sdx
        ylim <- c(0, max(sdu) * 1.05)
        ## 如果没有设置btm或top，自动计算
        if (is.null(btm) | is.null(top)) {
          autox <- .auto.breaks(dt = sdu, min.range = min.range, max.fold = max.fold)
          if (autox$flag) {
            btm <- autox$btm
            top <- autox$top
          } else {
            xx <- barplot(y, beside = TRUE, ylim = ylim, ...)
            if (!missing(sd.cols))
              errorbar(xx, y, sdu - y, horiz = FALSE, cex = cex.error)
            box()
            return(invisible(xx))
          }
        }
        ## Set up virtual y limits
        halflen <- btm - ylim[1]
        xlen <- halflen * 0.1 * gap.width
        v_tps1 <- btm + xlen  # virtual top positions
        v_tps2 <- v_tps1 + halflen * ratio
        v_ylim <- c(ylim[1], v_tps2)
        r_tps1 <- top  # real top positions
        r_tps2 <- ylim[2]
        ## Rescale data
        lmx <- summary(lm(c(v_tps1, v_tps2) ~ c(r_tps1, r_tps2)))
        lmx <- lmx$coefficients
        sel1 <- y > top
        sel2 <- y >= btm & y <= top
        y[sel1] <- y[sel1] * lmx[2] + lmx[1]
        y[sel2] <- btm + xlen/2
        sel1 <- sdd > top
        sel2 <- sdd >= btm & sdd <= top
        sdd[sel1] <- sdd[sel1] * lmx[2] + lmx[1]
        sdd[sel2] <- btm + xlen/2
        sel1 <- sdu > top
        sel2 <- sdu >= btm & sdu <= top
        sdu[sel1] <- sdu[sel1] * lmx[2] + lmx[1]
        sdu[sel2] <- btm + xlen/2
        ## bar plot
        xx <- barplot(y, beside = TRUE, ylim = v_ylim, axes = FALSE, names.arg = NULL,
                      ...)
        ## error bars
        if (!missing(sd.cols))
          errorbar(xx, y, sdu - y, horiz = FALSE, cex = cex.error)
        ## Real ticks and labels
        brks1 <- pretty(seq(0, btm, length = 10), n = 4)
        brks1 <- brks1[brks1 >= 0 & brks1 < btm]
        brks2 <- pretty(seq(top, r_tps2, length = 10), n = 4)
        brks2 <- brks2[brks2 > top & brks2 <= r_tps2]
        labx <- c(brks1, brks2)
        ## Virtual ticks
        brks <- c(brks1, brks2 * lmx[2] + lmx[1])
        axis(2, at = brks, labels = labx)
        box()
        ## break marks
        pos <- par("usr")
        xyratio <- (pos[2] - pos[1])/(pos[4] - pos[3])
        xlen <- (pos[2] - pos[1])/50 * brk.size
        px1 <- pos[1] - xlen
        px2 <- pos[1] + xlen
        px3 <- pos[2] - xlen
        px4 <- pos[2] + xlen
        py1 <- btm
        py2 <- v_tps1
        rect(px1, py1, px4, py2, col = brk.bg, xpd = TRUE, border = brk.bg)
        x1 <- c(px1, px1, px3, px3)
        x2 <- c(px2, px2, px4, px4)
        y1 <- c(py1, py2, py1, py2)
        y2 <- c(py1, py2, py1, py2)
        px <- .xy.adjust(x1, x2, y1, y2, xlen, xyratio, angle = brk.srt * pi/90)
        if (brk.type == "zigzag") {
          x1 <- c(x1, px1, px3)
          x2 <- c(x2, px2, px4)
          if (brk.srt > 90) {
            y1 <- c(y1, py2, py2)
            y2 <- c(y2, py1, py1)
          } else {
            y1 <- c(y1, py1, py1)
            y2 <- c(y2, py2, py2)
          }
        }
        if (brk.type == "zigzag") {
          px$x1 <- c(pos[1], px2, px1, pos[2], px4, px3)
          px$x2 <- c(px2, px1, pos[1], px4, px3, pos[2])
          mm <- (v_tps1 - btm)/3
          px$y1 <- rep(c(v_tps1, v_tps1 - mm, v_tps1 - 2 * mm), 2)
          px$y2 <- rep(c(v_tps1 - mm, v_tps1 - 2 * mm, btm), 2)
        }
        par(xpd = TRUE)
        segments(px$x1, px$y1, px$x2, px$y2, lty = 1, col = brk.col, lwd = brk.lwd)
        options(opts)
        par(xpd = FALSE)
        invisible(xx)
}
## 绘制误差线的函数
errorbar <- function(x, y, sd.lwr, sd.upr, horiz = FALSE, cex = 1, ...) {
  if (missing(sd.lwr) & missing(sd.upr))
    return(NULL)
  if (missing(sd.upr))
    sd.upr <- sd.lwr
  if (missing(sd.lwr))
    sd.lwr <- sd.upr
  if (!horiz) {
    arrows(x, y, y1 = y - sd.lwr, length = 0.1 * cex, angle = 90, ...)
    arrows(x, y, y1 = y + sd.upr, length = 0.1 * cex, angle = 90, ...)
  } else {
    arrows(y, x, x1 = y - sd.lwr, length = 0.1 * cex, angle = 90, ...)
    arrows(y, x, x1 = y + sd.upr, length = 0.1 * cex, angle = 90, ...)
  }
}
.xy.adjust <- function(x1, x2, y1, y2, xlen, xyratio, angle) {
  xx1 <- x1 - xlen * cos(angle)
  yy1 <- y1 + xlen * sin(angle)/xyratio
  xx2 <- x2 + xlen * cos(angle)
  yy2 <- y2 - xlen * sin(angle)/xyratio
  return(list(x1 = xx1, x2 = xx2, y1 = yy1, y2 = yy2))
}
## 自动计算断点位置的函数
.auto.breaks <- function(dt, min.range, max.fold) {
  datax <- sort(as.vector(dt))
  flags <- FALSE
  btm <- top <- NULL
  if (max(datax)/min(datax) < min.range)
    return(list(flag = flags, btm = btm, top = top))
  m <- max(datax)
  btm <- datax[2]
  i <- 3
  while (m/datax[i] > max.fold) {
    btm <- datax[i]
    flags <- TRUE
    i <- i + 1
  }
  if (flags) {
    btm <- btm + 0.05 * btm
    x <- 2
    top <- datax[i] * (x - 1)/x
    while (top < btm) {
      x <- x + 1
      top <- datax[i] * (x - 1)/x
      if (x > 100) {
        flags <- FALSE
        break
      }
    }
  }
  return(list(flag = flags, btm = btm, top = top))
}


####一、对27cancer的所有样本、tp53突变样本、不突变样本做柱状断点图gap，标显著性
#install.packages('plotrix')
library(plotrix)
PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
b1=factor(f[,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b2=factor(f[1:3894,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b3=factor(f[3895:nrow(f),5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
dataset0<- data.frame(value = f[,4], group=b1)
dataset1<- data.frame(value = f[1:3894,4], group=b2 )
dataset2<- data.frame(value = f[3895:nrow(f),4], group =b3)
dataset1[sapply(dataset1,is.na)]<-0
dataset2[sapply(dataset2,is.na)]<-0
a0=data.frame(tapply(dataset0$value,dataset0$group,median))
colnames(a0)="All"
a1=data.frame(tapply(dataset1$value,dataset1$group,median))
colnames(a1)="TP53 Mutation"
a2=data.frame(tapply(dataset2$value,dataset2$group,median))
colnames(a2)="No TP53 Mutation"
a3=data.frame(tapply(dataset0$value,dataset0$group,sd))
colnames(a3)="All"
a4=data.frame(tapply(dataset1$value,dataset1$group,sd))
colnames(a4)="TP53 Mutation"
a5=data.frame(tapply(dataset2$value,dataset2$group,sd))
colnames(a5)="No TP53 Mutation"
data=cbind(a0,a1,a2,a3/100,a4/100,a5/100)
data[sapply(data,is.na)]<-0
#brkt <- sample(c("normal", "zigzag"), 1)
pdf(file=paste(PATH,"boxplot/","Mutation Number of Different Cancers6.pdf",sep=""),width = 10,height=9)
gap.barplot(data, y.cols = 1:3, sd.cols = 4:6, col =c("#FF331C","#1DADB8","#12B826"), brk.type = "zigzag",gap.width=0.2,
            brk.size = 0.2, brk.lwd = 0.1, max.fold = 5, ratio = 2, cex.error = 0.2,border = NA,
            xlab="Cancer Type", ylab="Mutation Number",main="Mutation Number of Different Cancers" )
legend('topleft',cex=0.6,text.width=0.4,c("All","TP53 Mutation","No TP53 Mutation"),fill = c("#FF331C","#1DADB8","#12B826"))
axis(1,labels=rownames(data),at=seq(2,27*4,4),las=3)  
##i=25 TGCT、i=26 PCPG yes只有一个值，不能做t.test
#b4=c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA')
#for( i in 1:length(b4)){
  #if(i!=25&i!=26){
    #yes=dataset1[dataset1$group==as.character(b4[i]),1]
    #no=dataset2[dataset2$group==as.character(b4[i]),1]
    #round(t.test(yes,no)$p.value,3)
    #text(4*i-2,max(data[i,1:3]),round(t.test(yes,no)$p.value,4),cex=0.5,pos=1)
  #}
#}
signa=c("ns","ns","*","ns","ns","**","***","*","****","*","****","ns","**","*","****","ns","***","****","*","**","**","ns","ns","ns","ns","ns","ns")
#signa=c("","","*","","","**","***","*","****","*","****","","**","*","****","","***","****","*","**","**","","","","","","")
for( i in c(3,6,7,8,9,10,11,13,14,15,17,18,19,20,21)){
  #text(4*i-2,max(data[i,1:3]),signa[i],cex=0.8,pos=1)
  text(4*i-1.5,data[i,3],signa[i],cex=0.8,pos=3)
  #segments(4*i-2,data[i,3],4*i,data[i,3],col="black",lwd=1)
}
dev.off()


####一、对27cancer的tp53突变样本、不突变样本做柱状断点图gap，标显著性
####不要All，只有TP53 Mutation、No TP53 Mutation
library(plotrix)
PATH="F:/tp53/CBioportal_TP53/TCGAmutationdata/"
f<-read.csv(paste(PATH,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
b2=factor(f[1:3894,5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
b3=factor(f[3895:nrow(f),5],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP','TGCT','PCPG','THCA'))
dataset1<- data.frame(value = f[1:3894,4], group=b2 )
dataset2<- data.frame(value = f[3895:nrow(f),4], group =b3)
dataset1[sapply(dataset1,is.na)]<-0
dataset2[sapply(dataset2,is.na)]<-0
a1=data.frame(tapply(dataset1$value,dataset1$group,median))
colnames(a1)="TP53 Mutation"
a2=data.frame(tapply(dataset2$value,dataset2$group,median))
colnames(a2)="No TP53 Mutation"
a4=data.frame(tapply(dataset1$value,dataset1$group,sd))
colnames(a4)="TP53 Mutation"
a5=data.frame(tapply(dataset2$value,dataset2$group,sd))
colnames(a5)="No TP53 Mutation"
data=cbind(a1,a2,a4/100,a5/100)
data[sapply(data,is.na)]<-0
pdf(file=paste(PATH,"boxplot/","Mutation Number of Different Cancers7.pdf",sep=""),width = 10,height=9)
gap.barplot(data, y.cols = 1:2, sd.cols = 3:4, col =c("#1DADB8","#12B826"), brk.type = "zigzag",gap.width=0.2,
            brk.size = 0.2, brk.lwd = 0.1, max.fold = 5, ratio = 2, cex.error = 0.2,border = NA,
            xlab="Cancer Type", ylab="Mutation Number",main="Mutation Number of Different Cancers" )
legend('topleft',cex=0.6,text.width=0.4,c("TP53 Mutation","No TP53 Mutation"),fill = c("#1DADB8","#12B826"))
axis(1,labels=rownames(data),at=seq(2,27*3,3),las=3)  
signa=c("ns","ns","*","ns","ns","**","***","*","****","*","****","ns","**","*","****","ns","***","****","*","**","**","ns","ns","ns","ns","ns","ns")
for( i in c(3,6,7,8,9,10,11,13,14,15,17,18,19,20,21)){
  text(3*i-1,data[i,3],signa[i],cex=0.8,pos=3)
}
dev.off()


###########################################################################
#install.packages("reshape")
#install.packages("sciplot")
library(reshape)
library(sciplot)##画mean+sd
bargraph.CI(group, log(value), group = class, data = data, col = c("black", "red","green"), err.col = c("black", "red","green"), ci.fun = function(x) c(mean(x)-sd(x), mean(x)+sd(x)), xlab = "Dose", ylab = "Activity", ylim = c(0, 10), lwd = 2 ) #作图
legend("topright", legend = c("A", "B","C"), bty = "n", horiz = T, fill = c("black", "red","green"))#加标注
library(plotrix)##画断点gap
gap.barplot(data,gap=c(3,5),col=c("black", "red","green"),main='barplot with gap')
axis.break(2,3,breakcol='snow',style='gap')
axis.break(2,3*(1+0.02),breakcol='black',style='slash')
axis.break(4,3*(1+0.02),breakcol='black',style='slash')
axis(2,at=3)
###########################################################################


####二、对每个cancer,样本对应的hotspot的上下调基因的个数，bar上调在x上面，下调在x下面
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
for(i in filename[21:22]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  a1=cbind(as.character(a[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(a$HGVSp_Short),4))))
  a1=na.omit(a1)
  b1=data.frame(unique(a1[,2])) #up
  b2=data.frame(unique(a1[,2])) #down
  b3=data.frame(unique(a1[,2])) #normal
  b1[,2]=0
  b2[,2]=0
  b3[,2]=0
  for(j in 1:nrow(b1)){
    #gsub("-",".",a1[a1[,2]==j,1])
    if(TRUE %in% (gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1))){
      if(length(which(gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1)==TRUE))!=1){
        s1=0
        s2=0
        s3=0
        for(m in 1:length(gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1))){
          if(length(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]])!=0){
            s1=s1+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]>=1))
            s2=s2+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]<=-1))
            s3=s3+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]<1 & f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]>-1))
          }
        }
        b1[j,2]=s1
        b2[j,2]=s2
        b3[j,2]=s3
      }else{
        if(length(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])])!=0){
          b1[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]>=1))
          b2[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]<=-1))
          b3[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]<1 & f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]>-1))
        }
      }
    }
  }
  b1[,1]=as.numeric(as.character(b1[,1]))
  b2[,1]=as.numeric(as.character(b2[,1]))
  b3[,1]=as.numeric(as.character(b3[,1]))
  b1$pos=T
  b2[,2]=-b2[,2]
  b2$pos=F
  data=rbind(b1,b2)
  Expression_Value=factor(data[,3])
  p<-ggplot(data, aes(x=data[,1], y=data[,2], fill=Expression_Value)) +geom_bar(stat="identity", position="identity")+
    scale_fill_manual(values=c("green","red"))#position="identity"是为了关闭负值直方图没有定义的警告
  p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Position") + ylab("Gene Change Number") + ggtitle(gsub("_fc.txt","",i)) #添加横纵坐标，添加title
  p <- p+ theme(legend.position="top")+geom_hline(yintercept = 0, col="grey")+
    geom_text(aes(label=data[,2]),hjust=0.5, vjust=-0.5,size=2)
  ggsave(file=paste(path,"Position_GeneNumber/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
  #p1<-ggplot(b3, aes(x=b3[,1], y=b3[,2])) +geom_bar(stat="identity", position="identity")+geom_text(aes(label=b3[,2]),hjust=0.5, vjust=-0.5,size=2)+
    #theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))+ 
    #xlab("Position") + ylab("Gene Change Number") + ggtitle(gsub("_fc.txt","",i))
}


####二、对每个cancer,样本对应的hotspot的上下调基因的个数，bar上调、下调、正常都在x坐标上面
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  a1=cbind(as.character(a[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(a$HGVSp_Short),4))))
  a1=na.omit(a1)
  b1=data.frame(unique(a1[,2])) #up
  b2=data.frame(unique(a1[,2])) #down
  b3=data.frame(unique(a1[,2])) #normal
  b1[,2]=0
  b2[,2]=0
  b3[,2]=0
  for(j in 1:nrow(b1)){
    #gsub("-",".",a1[a1[,2]==j,1])
    if(TRUE %in% (gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1))){
      if(length(which(gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1)==TRUE))!=1){
        s1=0
        s2=0
        s3=0
        for(m in 1:length(gsub("-",".",a1[a1[,2]==b1[j,1],1]) %in% colnames(f1))){
          if(length(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]])!=0){
            s1=s1+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]>=1))
            s2=s2+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]<=-1))
            s3=s3+length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]<1 & f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])[m]]>-1))
          }
        }
        b1[j,2]=s1
        b2[j,2]=s2
        b3[j,2]=s3
      }else{
        if(length(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])])!=0){
          b1[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]>=1))
          b2[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]<=-1))
          b3[j,2]=length(which(f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]<1 & f1[,colnames(f1)==gsub("-",".",a1[a1[,2]==b1[j,1],1])]>-1))
        }
      }
    }
  }
  b1[,1]=as.numeric(as.character(b1[,1]))
  b2[,1]=as.numeric(as.character(b2[,1]))
  b3[,1]=as.numeric(as.character(b3[,1]))
  da=rbind(cbind(class="up",b1),cbind(class="down",b2),cbind(class="normal",b3))
  Expression_Value=factor(da[,1])    # 条形图函数：position设置条形图类型为簇状
  p<-ggplot(da, aes(x = da[,2], y = da[,3], fill = Expression_Value)) +geom_bar(position = "dodge", stat = "identity")+scale_fill_manual(values=c("red","green","grey"))
  p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Position") + ylab("Gene Change Number") + ggtitle(gsub("_fc.txt","",i)) #添加横纵坐标，添加title
  p <- p+ theme(legend.position="top")+geom_text(aes(label=da[,3]),hjust=0.5, vjust=-0.5,size=2)+
    annotate("text",x=da[,2],y=-10, colour="black",size=1.5,label=da[,2])
  ggsave(file=paste(path,"Position_GeneNumber/1/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}


####二、对每个cancer,hotspot对应的样本中tp53基因上下调的样本的个数，bar上调、下调、正常都在x坐标上面
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
ff= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t") 
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  colnames(f1)<-gsub("_no","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  a1=cbind(as.character(a[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(a$HGVSp_Short),4))))
  a1=na.omit(a1)
  a1=a1[a1[,2] %in% names( table(a1[,2])[table(a1[,2])>=6]),]#选择突变样本数>=6的位点
  b=data.frame(rep(unique(a1[,2]),each=3))#up#down#normal
  b[,2]=rep(c("up","down","normal"),length(unique(a1[,2])))
  for(n in as.character(unique(ff[,2]))[2:length(as.character(unique(ff[,2])))]){
    e=c("TP53",as.character(ff[as.character(ff[,2])==n,1]))#每个功能对应的基因
    da=data.frame()
    for(m in 1:length(e)){
      for(j in 1:nrow(b)){
        if(TRUE %in% (gsub("-",".",a1[a1[,2]==as.character(b[j,1]),1]) %in% colnames(f1))){
          if(e[m] %in% f1[,1] & gsub("-",".",a1[a1[,2]==as.character(b[j,1]),1]) %in% colnames(f1) ){
            b1=f1[f1[,1]==e[m],intersect(gsub("-",".",a1[a1[,2]==as.character(b[j,1]),1]),colnames(f1))]
            b[,3]=e[m]
            b[b[,1]==as.character(b[j,1]) & b[,2]=="up",4]=length(which(b1[1,]>=1))
            b[b[,1]==as.character(b[j,1]) & b[,2]=="down",4]=length(which(b1[1,]<=-1))
            b[b[,1]==as.character(b[j,1]) & b[,2]=="normal",4]=length(which(b1[1,]<1 & b1[1,]>-1))
            #b=na.omit(b)
          }
  
        }
      }
      #b=na.omit(b)
      da=rbind(da,b)
    }
    da[,1]=as.numeric(as.character(da[,1]))
    if(nrow(da)!=0){
      Expression_Value=factor(da[,2])    # 条形图函数：position设置条形图类型为簇状
      p<-ggplot(da, aes(x = da[,1], y = da[,4], fill = Expression_Value)) +geom_bar(position = "dodge", stat = "identity")+scale_fill_manual(values=c("green","grey","red"))
      p <- p+ facet_grid(V3 ~ .)
      p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
      p <- p+xlab("Position")+ylab("Gene Change Number")+ggtitle(paste(gsub("_fc.txt","",i),n,sep="_"))+xlim(0,400) #添加横纵坐标，添加title
      ggsave(file=paste(path,"Position_GeneNumber/2/",paste(gsub("_fc.txt","",i),gsub("/","&",n),sep="_"),".pdf",sep=""),width = 10,height=9)
      
    }
    
  }
}







####三、27cancer的phenotype分析
path="D:/tp53/"
f= read.table(paste(path,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
#cancer="CESC"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  a=f[as.character(f[,5])==cancer,c(1,3)]
  a[is.na(as.character(a[,2])),2]<-0
  e1<-as.character(a[a[,2]==0,1])   #NoMutation
  e2<-as.character(a[a[,2]!=0,1])   #Mutation
  e3<-as.character(a[a[,2]==1,1])   #SingleMutation
  e4<-as.character(a[a[,2]>1,1])    #CoMutation
  f1= read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
  #f1$TCGA=substring(as.character(f1[,1]),1,12) 
  #f1[substring(as.character(f1[,1]),1,12) %in% e1,]
  #f1[substring(as.character(f1[,1]),1,12) %in% e2,]
  #f1[substring(as.character(f1[,1]),1,12) %in% e3,]
  #f1[substring(as.character(f1[,1]),1,12) %in% e4,]
  write.table(f1[substring(as.character(f1[,1]),1,12) %in% e1,],file =paste(path,"GDC_phenotype/",cancer,"_NoMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(f1[substring(as.character(f1[,1]),1,12) %in% e2,],file =paste(path,"GDC_phenotype/",cancer,"_Mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(f1[substring(as.character(f1[,1]),1,12) %in% e3,],file =paste(path,"GDC_phenotype/",cancer,"_SingleMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(f1[substring(as.character(f1[,1]),1,12) %in% e4,],file =paste(path,"GDC_phenotype/",cancer,"_CoMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  f2= read.csv(paste(path,"phenotype/TCGA-",cancer,".survival.tsv",sep=""),header=T,sep="\t")
  write.table(unique(f2[substring(as.character(f2[,1]),1,12) %in% e1,2:ncol(f2)]),file =paste(path,"GDC_survival/",cancer,"_NoMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(unique(f2[substring(as.character(f2[,1]),1,12) %in% e2,2:ncol(f2)]),file =paste(path,"GDC_survival/",cancer,"_Mutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(unique(f2[substring(as.character(f2[,1]),1,12) %in% e3,2:ncol(f2)]),file =paste(path,"GDC_survival/",cancer,"_SingleMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(unique(f2[substring(as.character(f2[,1]),1,12) %in% e4,2:ncol(f2)]),file =paste(path,"GDC_survival/",cancer,"_CoMutation.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
}




#######################################################
####27cancerGDC_phenotype要去冗余、GDC_survival不需要
path="D:/tp53/"
#cancer= "COAD","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PCPG","READ","SKCM","TGCT","UCEC",
#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
#"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
###GDC_phenotype
for(cancer in c("TGCT") ){
  f1= read.csv(paste(path,"GDC_phenotype/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
  f2= read.csv(paste(path,"GDC_phenotype/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
  f3= read.csv(paste(path,"GDC_phenotype/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
  f4= read.csv(paste(path,"GDC_phenotype/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
  write.table(paste("Variables","NoMutation","Mutation","SingleMutation","CoMutation","No_Mutation_pValue","Single_CoMutation_pValue",sep="\t"),file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  for (i in colnames(f1)[2:ncol(f1)]){
    a1=unique(cbind(substring(as.character(f1[,1]),1,12),gsub("^$|i",NA,as.character(f1[,colnames(f1)==i]))))[,2]
    a2=unique(cbind(substring(as.character(f2[,1]),1,12),gsub("^$|i",NA,as.character(f2[,colnames(f2)==i]))))[,2]
    a3=unique(cbind(substring(as.character(f3[,1]),1,12),gsub("^$|i",NA,as.character(f3[,colnames(f3)==i]))))[,2]
    a4=unique(cbind(substring(as.character(f4[,1]),1,12),gsub("^$|i",NA,as.character(f4[,colnames(f4)==i]))))[,2]
    a11=a1[!is.na(a1)]
    a21=a2[!is.na(a2)]
    a31=a3[!is.na(a3)]
    a41=a4[!is.na(a4)]
    if(i %in% c("additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","anatomic_treatment_site",
                "birth_control_pill_history_usage_category","histological_type","lost_follow_up","measure_of_response","neoplasm_histologic_grade",
                "person_neoplasm_cancer_status","primary_therapy_outcome_success","new_neoplasm_event_occurrence_anatomic_site","new_neoplasm_event_type",                             
                "new_tumor_event_additional_surgery_procedure","new_tumor_event_after_initial_treatment","radiation_therapy","radiation_treatment_ongoing",                         
                "route_of_administration","targeted_molecular_therapy","ethnicity.demographic",                              
                "gender.demographic","race.demographic","tumor_stage.diagnoses","vital_status.diagnoses","bcr_id.tissue_source_site","sample_type.samples",                                 
                "state.samples","oct_embedded.samples") ){
      ##离散型
      if(length(table(a11))>1 & length(table(a21))>1){
        p1=fisher.test(rbind(table(a11),table(a21)))$p.value  #NoMutation_Mutation
      }else{p1=-5}
      if(length(table(a31))>1 & length(table(a41))>1){
        p2=fisher.test(rbind(table(a31),table(a41)))$p.value  #SingleMutation_CoMutation
      }else{p2=-5}
      g1=paste(i,"","","","",round(p1,4),round(p2,4),sep="\t")
      g2=paste(names(table(a21)),table(a11),table(a21),table(a31),table(a41),"","",sep="\t")
      write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      write.table(g2,file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
    if(i=="clinical_stage"){
      ##离散型
      b=data.frame(V1=union(union(union(names(table(a11)),names(table(a21))),names(table(a31))),names(table(a41))))
      rownames(b)=b[,1]
      b[names(table(a11)),2]=data.frame(table(a11))[,2]
      b[names(table(a21)),3]=data.frame(table(a21))[,2]
      b[names(table(a31)),4]=data.frame(table(a31))[,2]
      b[names(table(a41)),5]=data.frame(table(a41))[,2]
      #p1=fisher.test(b[,2],b[,3])$p.value  #NoMutation_Mutation
      #p2=fisher.test(b[,4],b[,5])$p.value  #SingleMutation_CoMutation
      #g1=paste(i,"","","","",round(p1,4),round(p2,4),sep="\t")
      g1=paste(i,"","","","",-2,-2,sep="\t")
      write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      write.table(b,file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
    if(i %in% c("age_at_initial_pathologic_diagnosis",   
                "days_to_drug_therapy_end","days_to_drug_therapy_start",
                "days_to_new_tumor_event_after_initial_treatment","days_to_radiation_therapy_end","days_to_radiation_therapy_start",
                "radiation_dosage","age_at_diagnosis.diagnoses","days_to_birth.diagnoses","days_to_death.diagnoses",                             
                "days_to_last_follow_up.diagnoses","days_to_collection.samples","initial_weight.samples")){
      ##连续型
      if(length(as.numeric(a11))>=2 & length(as.numeric(a21))>=2){
        p1=t.test(as.numeric(a11),as.numeric(a21))$p.value  #NoMutation_Mutation
      }else{p1=-5}
      if(length(as.numeric(a31))>=2 & length(as.numeric(a41))>=2){
        p2=t.test(as.numeric(a31),as.numeric(a41))$p.value  #SingleMutation_CoMutation
      }else{p2=-5}
      g1=paste(i,paste(round(mean(as.numeric(a11)),4),round(sd(as.numeric(a11)),4),sep="+"),paste(round(mean(as.numeric(a21)),4),round(sd(as.numeric(a21)),4),sep="+"),paste(round(mean(as.numeric(a31)),4),round(sd(as.numeric(a31)),4),sep="+"),paste(round(mean(as.numeric(a41)),4),round(sd(as.numeric(a41)),4),sep="+"),round(p1,4),round(p2,4),sep="\t")
      write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      d<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段
      NO1<-(max(na.omit(as.numeric(a11)))-min(na.omit(as.numeric(a11))))%/%3+min(na.omit(as.numeric(a11))) #min~NO1
      NO2<-(max(na.omit(as.numeric(a11)))-min(na.omit(as.numeric(a11))))%/%3*2+min(na.omit(as.numeric(a11)))#NO1~NO2#NO2~max
      row.names(d)=c(paste(min(na.omit(as.numeric(a11))),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(na.omit(as.numeric(a11)))+1))
      #NO1<-(max(na.omit(as.numeric(a21)))-min(na.omit(as.numeric(a21))))%/%3+min(na.omit(as.numeric(a21))) #min~NO1
      #NO2<-(max(na.omit(as.numeric(a21)))-min(na.omit(as.numeric(a21))))%/%3*2+min(na.omit(as.numeric(a21)))#NO1~NO2#NO2~max
      #row.names(d)=c(paste(min(na.omit(as.numeric(a21))),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(na.omit(as.numeric(a21)))+1))
      if(length(a11)!=0){
        d1=as.data.frame(as.numeric(a11))
        d1[,2]=1
        da1={}
        for(j in row.names(d)){
          data=d1[as.numeric(d1[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d1[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
          da1=append(da1,sum(data[,2]))
        }
        names(da1)=row.names(d)
        da1<-as.data.frame(da1)
      }else{
        da1<-as.data.frame(-1,-1,-1)
      }
      if(length(a21)!=0){
        d2=as.data.frame(as.numeric(a21))
        d2[,2]=1
        da2={}
        for(j in row.names(d)){
          data=d2[as.numeric(d2[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d2[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
          da2=append(da2,sum(data[,2]))
        }
        names(da2)=row.names(d)
        da2<-as.data.frame(da2)
      }else{
        da2<-as.data.frame(-1,-1,-1)
      }
      if(length(a31)!=0){
        d3=as.data.frame(as.numeric(a31))
        d3[,2]=1
        da3={}
        for(j in row.names(d)){
          data=d3[as.numeric(d3[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d3[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
          da3=append(da3,sum(data[,2]))
        }
        names(da3)=row.names(d)
        da3<-as.data.frame(da3)
      }else{
        da3<-as.data.frame(-1,-1,-1)
      }
      if(length(a41)!=0){
        d4=as.data.frame(as.numeric(a41))
        d4[,2]=1
        da4={}
        for(j in row.names(d)){
          data=d4[as.numeric(d4[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d4[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
          da4=append(da4,sum(data[,2]))
        }
        names(da4)=row.names(d)
        da4<-as.data.frame(da4)
      }else{
        da4<-as.data.frame(-1,-1,-1)
      }
      da=cbind(da1,da2,da3,da4)
      write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
  }
}



######################################################
###GDC_survival
path="D:/tp53/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC",
                  "KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD",
                  "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  f5= read.csv(paste(path,"GDC_survival/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
  f6= read.csv(paste(path,"GDC_survival/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
  f7= read.csv(paste(path,"GDC_survival/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
  f8= read.csv(paste(path,"GDC_survival/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
  write.table(paste("Variables","NoMutation","Mutation","SingleMutation","CoMutation","No_Mutation_pValue","Single_CoMutation_pValue",sep="\t"),file =paste(path,"phenotype_Analysis/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  for (i in colnames(f5)[c(1,3:ncol(f5))]){
    a5=as.character(f5[,colnames(f5)==i])
    a6=as.character(f6[,colnames(f6)==i])
    a7=as.character(f7[,colnames(f7)==i])
    a8=as.character(f8[,colnames(f8)==i])
    if(i %in% c("X_EVENT","X_OS_IND")){
      ##离散型
      if(length(table(a5))>1 & length(table(a6))>1){
        p1=fisher.test(rbind(table(a5),table(a6)))$p.value  #NoMutation_Mutation
      }else{p1=-5}
      if(length(table(a7))>1 & length(table(a8))>1){
        p2=fisher.test(rbind(table(a7),table(a8)))$p.value  #SingleMutation_CoMutation
      }else{p2=-5}
      g1=paste(i,"","","","",round(p1,4),round(p2,4),sep="\t")
      g2=paste(names(table(a5)),table(a5),table(a6),table(a7),table(a8),"","",sep="\t")
      write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      write.table(g2,file =paste(path,"phenotype_Analysis/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
    if(i %in% c("X_TIME_TO_EVENT","X_OS")){
      ##连续型
      if(length(as.numeric(a5))>=2 & length(as.numeric(a6))>=2){
        p1=t.test(as.numeric(a5),as.numeric(a6))$p.value  #NoMutation_Mutation
      }else{p1=-5}
      if(length(as.numeric(a7))>=2 & length(as.numeric(a8))>=2){
        p2=t.test(as.numeric(a7),as.numeric(a8))$p.value  #SingleMutation_CoMutation
      }else{p2=-5}
      g1=paste(i,paste(round(mean(as.numeric(a5)),4),round(sd(as.numeric(a5)),4),sep="+"),paste(round(mean(as.numeric(a6)),4),round(sd(as.numeric(a6)),4),sep="+"),paste(round(mean(as.numeric(a7)),4),round(sd(as.numeric(a7)),4),sep="+"),paste(round(mean(as.numeric(a8)),4),round(sd(as.numeric(a8)),4),sep="+"),round(p1,4),round(p2,4),sep="\t")
      write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      d<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段
      NO1<-(max(as.numeric(a5))-min(as.numeric(a5)))%/%3+min(as.numeric(a5)) #min~NO1
      NO2<-(max(as.numeric(a5))-min(as.numeric(a5)))%/%3*2+min(as.numeric(a5))#NO1~NO2#NO2~max
      row.names(d)=c(paste(min(as.numeric(a5)),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(as.numeric(a5))))
      d1=as.data.frame(as.numeric(a5))
      d1[,2]=1
      da1={}
      for(j in row.names(d)){
        data=d1[as.numeric(d1[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d1[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
        da1=append(da1,sum(data[,2]))
      }
      names(da1)=row.names(d)
      da1<-as.data.frame(da1)
      d2=as.data.frame(as.numeric(a6))
      d2[,2]=1
      da2={}
      for(j in row.names(d)){
        data=d2[as.numeric(d2[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d2[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
        da2=append(da2,sum(data[,2]))
      }
      names(da2)=row.names(d)
      da2<-as.data.frame(da2)
      d3=as.data.frame(as.numeric(a7))
      d3[,2]=1
      da3={}
      for(j in row.names(d)){
        data=d3[as.numeric(d3[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d3[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
        da3=append(da3,sum(data[,2]))
      }
      names(da3)=row.names(d)
      da3<-as.data.frame(da3)
      if(length(a8)!=0){
        d4=as.data.frame(as.numeric(a8))
        d4[,2]=1
        da4={}
        for(j in row.names(d)){
          data=d4[as.numeric(d4[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d4[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
          da4=append(da4,sum(data[,2]))
        }
        names(da4)=row.names(d)
        da4<-as.data.frame(da4)
      }else{
        da4<-as.data.frame(-1,-1,-1)
      }
      da=cbind(da1,da2,da3,da4)
      write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
    }
  }
}


##########################################################

#purrr::map_chr(b, ~class(.x))         
#purrr::map_lgl(a, ~is.factor(.x))
#purrr::map_lgl(b , ~is.double(.x))
#purrr::map_lgl(a , ~is.integer(.x))
#purrr::map_lgl(a , ~is.character(.x))
##离散型
c("additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure","anatomic_treatment_site",
"birth_control_pill_history_usage_category","clinical_stage","histological_type","lost_follow_up","measure_of_response","neoplasm_histologic_grade",
"person_neoplasm_cancer_status","primary_therapy_outcome_success","new_neoplasm_event_occurrence_anatomic_site","new_neoplasm_event_type",                             
"new_tumor_event_additional_surgery_procedure","new_tumor_event_after_initial_treatment","radiation_therapy","radiation_treatment_ongoing",                         
"radiation_type","route_of_administration","system_version","targeted_molecular_therapy","vial_number","ethnicity.demographic",                              
"gender.demographic","race.demographic","tumor_stage.diagnoses","vital_status.diagnoses","bcr_id.tissue_source_site","sample_type.samples",                                 
"state.samples")
##连续型
c("age_at_initial_pathologic_diagnosis","day_of_dcc_upload","day_of_form_completion","days_to_additional_surgery_metastatic_procedure",     
"days_to_drug_therapy_end","days_to_drug_therapy_start","days_to_initial_pathologic_diagnosis","days_to_new_tumor_event_additional_surgery_procedure",
"days_to_new_tumor_event_after_initial_treatment","days_to_radiation_therapy_end","days_to_radiation_therapy_start","month_of_form_completion", 
"number_cycles","numfractions","radiation_dosage","age_at_diagnosis.diagnoses","days_to_birth.diagnoses","days_to_death.diagnoses",                             
"days_to_last_follow_up.diagnoses","days_to_collection.samples","initial_weight.samples","oct_embedded.samples")   
                                                           

         
                         
                             
  


####四、对20个cancer的突变、不突变样本做6种免疫细胞类型的purity纯度的boxplot
path="D:/tp53/"
#cancer="CESC"
f2= read.csv(paste(path,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
for(cancer in c("BLCA","BRCA","CESC","COAD","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV",
                "PRAD","READ","SKCM","STAD","THCA","UCEC","UCS")){
  f1= read.csv(paste(path,"AGPall/AGP-",tolower(cancer),".txt",sep=""),header=T,sep="\t")
  a=f2[as.character(f2[,5])==cancer,]
  Mu=as.character(na.omit(a[a[,2]==1,1]))
  NoMu=as.character(a[is.na(a[,2]),1])
  a1=data.frame(cbind(gsub("\\.","-",substring(as.character(f1[,1]),1,12)),f1[,5],f1[,2]))
  d1=a1[as.character(a1[,1]) %in% Mu,2:3]
  d2=a1[as.character(a1[,1]) %in% NoMu,2:3]
  d1[,1]=as.numeric(as.character(d1[,1]))
  d1[,2]=as.numeric(as.character(d1[,2]))
  d2[,1]=as.numeric(as.character(d2[,1]))
  d2[,2]=as.numeric(as.character(d2[,2]))
  d1=d1[order(d1$X2, -d1$X3),]
  d2=d2[order(d2$X2, -d2$X3),]
  pdf(paste(path,"AGPall_boxplot/",cancer,".pdf",sep=""),width = 9,height=9)
  x1<-boxplot(d1[,2] ~ d1[,1],d1,border= "red",pch=".",boxwex = 0.2, at = unique(d1[,1])-0.2,main=cancer,xlab="Immune Cell Types",ylab="purity",las=1,font.lab=2,xlim=c(0,8),ylim=c(0,1))
  boxplot(d2[,2] ~ d2[,1],d2,border= "green",pch=".",boxwex = 0.2, at = unique(d2[,1])+0.2,add= TRUE,las=1,axes=TRUE)
  legend('topleft',cex=0.6,text.width=0.4,c("TP53 Mutation","No TP53 Mutation"),fill = c("red","green"))
  for(i in 1:nrow(data.frame(table(d1[,1])))){
    text(as.numeric(as.character(data.frame(table(d1[,1]))[i,1]))-0.3,0.9,labels=data.frame(table(d1[,1]))[i,2],cex=0.6,pos=2)
    text(as.numeric(as.character(data.frame(table(d2[,1]))[i,1]))+0.3,0.9,labels=data.frame(table(d2[,1]))[i,2],cex=0.6,pos=2)
    j=as.numeric(as.character(data.frame(table(d1[,1]))[i,1]))
    if(length(d1[d1$X2==j,2])>=2 & length(d2[d2$X2==j,2])>=2){
      text(as.numeric(as.character(data.frame(table(d2[,1]))[i,1])),0.95,labels=t.test(d1[d1$X2==j,2],d2[d2$X2==j,2])$p.value,cex=0.6,pos=2)
    }
  }
  dev.off()
}
 


####四、对20个cancer的突变、不突变样本做2\4\6种免疫细胞类型的purity纯度的boxplot整合为一图
library(ggplot2)
path="D:/tp53/"
f2= read.csv(paste(path,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
g=data.frame()
for(cancer in c("BLCA","BRCA","CESC","COAD","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV",
                "PRAD","READ","SKCM","STAD","THCA","UCEC","UCS")){
  f1= read.csv(paste(path,"AGPall/AGP-",tolower(cancer),".txt",sep=""),header=T,sep="\t")
  a=f2[as.character(f2[,5])==cancer,]
  Mu=as.character(na.omit(a[a[,2]==1,1]))
  NoMu=as.character(a[is.na(a[,2]),1])
  a1=data.frame(cbind(gsub("\\.","-",substring(as.character(f1[,1]),1,12)),f1[,5],f1[,2]))
  d1=a1[as.character(a1[,1]) %in% Mu,2:3]
  d2=a1[as.character(a1[,1]) %in% NoMu,2:3]
  d1[,1]=as.numeric(as.character(d1[,1]))
  d1[,2]=as.numeric(as.character(d1[,2]))
  d2[,1]=as.numeric(as.character(d2[,1]))
  d2[,2]=as.numeric(as.character(d2[,2]))
  d1=cbind(cancer,state="Mu",d1)
  d2=cbind(cancer,state="NoMu",d2)
  g=rbind(g,d1,d2)
}
for(i in c(2,4,6)){
  g1=g[g[,3]==i,c(1,2,4)]
  g2=data.frame()
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])>=2){
      g2=rbind(g2,cbind(j,t.test(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3],g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","green"))+coord_flip()
  p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p + xlab("Cancer") + ylab("purity") + ggtitle(paste("Immune Cell Types=",i,sep="")) #添加横纵坐标，添加title
  p <- p+ theme(legend.position="top")+
    annotate("text",x=g2[,1], y=1, colour="black",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/1/",i,".pdf",sep=""),width = 10,height=9)
}
#fisher.test(d1[d1$X2==j,2],d2[d2$X2==j,2])$p.value
#wilcox.test(d1[d1$X2==j,2],d2[d2$X2==j,2])
#cor.test(d1[d1$X2==j,2],d2[d2$X2==j,2], method = "spearman") 
#t.test(d1[d1$X2==j,2],d2[d2$X2==j,2])$p.value



####四、对27个cancer的突变、不突变样本、单突变样本、多突变样本做6种免疫细胞类型的purity纯度的boxplot整合为一图
library(ggplot2)
path="D:/tp53/"
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
f2= read.csv(paste(path,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
g=data.frame()
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
                "LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  a=f2[as.character(f2[,5])==cancer,]
  Mu=as.character(na.omit(a[a[,2]==1,1]))
  NoMu=as.character(a[is.na(a[,2]),1])
  SigleMu=as.character(na.omit(a[a[,3]==1,1]))
  CoMu=as.character(na.omit(a[a[,3]>1,1]))
  a1=f1[grep(".*-01$",f1[,1]),]
  a1[,1]=gsub("-01","",a1[,1])
  d1=a1[as.character(a1[,1]) %in% Mu,2:7]
  d2=a1[as.character(a1[,1]) %in% NoMu,2:7]
  d3=a1[as.character(a1[,1]) %in% SigleMu,2:7]
  d4=a1[as.character(a1[,1]) %in% CoMu,2:7]
  d11=data.frame()
  d21=data.frame()
  d31=data.frame()
  d41=data.frame()
  for(m in 1:ncol(d1)){
    d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
  }
  for(m in 1:ncol(d2)){
    d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
  }
  for(m in 1:ncol(d3)){
    d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
  }
  for(m in 1:ncol(d4)){
    d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
  }
  if(nrow(d11)!=0){
  d12=cbind(cancer,state="Mu",d11)
  }
  if(nrow(d21)!=0){
  d22=cbind(cancer,state="NoMu",d21)
  }
  if(nrow(d31)!=0){
  d32=cbind(cancer,state="SigleMu",d31)
  }
  if(nrow(d41)!=0){
    d42=cbind(cancer,state="CoMu",d41)
  }
  if(nrow(d12)!=0){
    g=rbind(g,d12)
  }
  if(nrow(d22)!=0){
    g=rbind(g,d22)
  }
  if(nrow(d32)!=0){
    g=rbind(g,d32)
  }
  if(nrow(d42)!=0){
    g=rbind(g,d42)
  }
  #g=rbind(g,d12,d22,d32,d42)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #Mu_NoMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #SigleMu_CoMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="SigleMu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="CoMu",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="SigleMu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="CoMu",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","grey","blue","green"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=g2[,1], y=0.7, colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))+
    annotate("text",x=g3[,1], y=0.9, colour="blue",size=2,label=paste("p=",as.numeric(as.character(g3[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/2/",i,".pdf",sep=""),width = 10,height=9)
}
#lbls<-paste(" ",a)


####四、对27个cancer的不突变样本、突变样本、tp53表达上调、下调、正常的样本做6种免疫细胞类型的purity纯度的boxplot
path="D:/tp53/"
f= read.csv(paste(path,"3894_5886=9780+1.txt",sep=""),header=T,sep="\t") ##添加标签tp53表达上调下调正常
f[,5]=NA
for(i in 1:nrow(f)){
  if(paste(as.character(f[i,6]),"_fc.txt",sep="") %in% list.files("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/")){
    f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",as.character(f[i,6]),"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
    colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
    colnames(f1)=gsub("_no","",colnames(f1))
    f1[,2]<-as.numeric(f1[,2])
    f1[sapply(f1,is.infinite)]<-NA
    f1[sapply(f1,is.na)]<-0
    f1<-f1[-1,]
    a=f1[as.character(f1[,1])=="TP53",gsub("-",".",as.character(f[i,1]))]
    if(length(a)!=0){
      if(as.numeric(a)>=1){
        f[i,5]="up"
      }
      if(as.numeric(a)<=-1){
        f[i,5]="down"
      }
      if(as.numeric(a)>-1 & as.numeric(a)<1){
        f[i,5]="normal"
      }
    }
  } 
}
write.table(f,file =paste(path,"3894_5886=9780+1_1.txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)

library(ggplot2)
path="D:/tp53/"
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
f2= read.csv(paste(path,"3894_5886=9780+1_1.txt",sep=""),header=T,sep="\t")
g=data.frame()
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD",
                "LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  a=f2[as.character(f2[,6])==cancer,]
  Mu=as.character(na.omit(a[a[,2]==1,1]))
  NoMu=as.character(a[is.na(a[,2]),1])
  UpMu=as.character(na.omit(a[a[,2]==1 & as.character(a[,5])=="up",1]))
  DownMu=as.character(na.omit(a[a[,2]==1 & as.character(a[,5])=="down",1]))
  NormalMu=as.character(na.omit(a[a[,2]==1 & as.character(a[,5])=="normal",1]))
  a1=f1[grep(".*-01$",f1[,1]),]
  a1[,1]=gsub("-01","",a1[,1])
  d1=a1[as.character(a1[,1]) %in% Mu,2:7]
  d2=a1[as.character(a1[,1]) %in% NoMu,2:7]
  d3=a1[as.character(a1[,1]) %in% UpMu,2:7]
  d4=a1[as.character(a1[,1]) %in% DownMu,2:7]
  d5=a1[as.character(a1[,1]) %in% NormalMu,2:7]
  d11=data.frame()
  d21=data.frame()
  d31=data.frame()
  d41=data.frame()
  d51=data.frame()
  for(m in 1:ncol(d1)){
    d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
  }
  for(m in 1:ncol(d2)){
    d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
  }
  for(m in 1:ncol(d3)){
    d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
  }
  for(m in 1:ncol(d4)){
    d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
  }
  for(m in 1:ncol(d5)){
    d51=rbind(d51,cbind(rep(colnames(d5)[m],nrow(d5)),d5[,m]))
  }
  if(nrow(d11)!=0){
    d12=cbind(cancer,state="Mu",d11)
  }
  if(nrow(d21)!=0){
    d22=cbind(cancer,state="NoMu",d21)
  }
  if(nrow(d31)!=0){
    d32=cbind(cancer,state="UpMu",d31)
  }
  if(nrow(d41)!=0){
    d42=cbind(cancer,state="DownMu",d41)
  }
  if(nrow(d51)!=0){
    d52=cbind(cancer,state="NormalMu",d51)
  }
  if(nrow(d12)!=0){
    g=rbind(g,d12)
  }
  if(nrow(d22)!=0){
    g=rbind(g,d22)
  }
  if(nrow(d32)!=0){
    g=rbind(g,d32)
  }
  if(nrow(d42)!=0){
    g=rbind(g,d42)
  }
  if(nrow(d52)!=0){
    g=rbind(g,d52)
  }
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #Mu_NoMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NoMu",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #Mu_UpMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="UpMu",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="UpMu",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g4=data.frame()  #Mu_DownMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="DownMu",3])>=2){
      g4=rbind(g4,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="DownMu",3])))$p.value))
    }else{
      g4=rbind(g4,cbind(j,-5))
    }
  }
  g5=data.frame()  #Mu_NormalMu$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NormalMu",3])>=2){
      g5=rbind(g5,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="Mu",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="NormalMu",3])))$p.value))
    }else{
      g5=rbind(g5,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","grey","blue","green","purple"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=g2[,1], y=0.7, colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))+
    annotate("text",x=g3[,1], y=0.85, colour="blue",size=2,label=paste("p=",as.numeric(as.character(g3[,2])),sep=""))+
    annotate("text",x=g4[,1], y=0.98, colour="green",size=2,label=paste("p=",as.numeric(as.character(g4[,2])),sep=""))+
    annotate("text",x=g5[,1], y=1.1, colour="purple",size=2,label=paste("p=",as.numeric(as.character(g5[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/3/",i,".pdf",sep=""),width = 10,height=9)
}



####四、对BRCA、HNSC的不同亚型的样本做6种免疫细胞类型的purity纯度的boxplot
library(ggplot2)
path="D:/tp53/"
cancer="BRCA"  ##4亚型
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$breast_carcinoma_estrogen_receptor_status),as.character(f$breast_carcinoma_progesterone_receptor_status),as.character(f$lab_proc_her2_neu_immunohistochemistry_receptor_status))))
a[,1]=as.character(a[,1])
a[,2]=as.character(a[,2])
a[,3]=as.character(a[,3])
a[,4]=as.character(a[,4])
luminalA=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Negative",1] ##和或
luminalB=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Positive",1] ##和或
her2=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Positive",1]
basal=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Negative",1]
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-01","",a1[,1])
d1=a1[as.character(a1[,1]) %in% luminalA,2:7]
d2=a1[as.character(a1[,1]) %in% luminalB,2:7]
d3=a1[as.character(a1[,1]) %in% her2,2:7]
d4=a1[as.character(a1[,1]) %in% basal,2:7]
d11=data.frame()
d21=data.frame()
d31=data.frame()
d41=data.frame()
for(m in 1:ncol(d1)){
  d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
}
for(m in 1:ncol(d2)){
  d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
}
for(m in 1:ncol(d3)){
  d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
}
for(m in 1:ncol(d4)){
  d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
}
if(nrow(d11)!=0){
  d12=cbind(cancer,state="luminalA",d11)
}
if(nrow(d21)!=0){
  d22=cbind(cancer,state="luminalB",d21)
}
if(nrow(d31)!=0){
  d32=cbind(cancer,state="her2",d31)
}
if(nrow(d41)!=0){
  d42=cbind(cancer,state="basal",d41)
}
g=data.frame()
if(nrow(d12)!=0){
  g=rbind(g,d12)
}
if(nrow(d22)!=0){
  g=rbind(g,d22)
}
if(nrow(d32)!=0){
  g=rbind(g,d32)
}
if(nrow(d42)!=0){
  g=rbind(g,d42)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #basal_luminalA$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #basal_luminalB$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g4=data.frame()  #basal_her2$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])>=2){
      g4=rbind(g4,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])))$p.value))
    }else{
      g4=rbind(g4,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","blue","green","grey"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=g2[,1], y=0.7, colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))+
    annotate("text",x=g3[,1], y=0.85, colour="blue",size=2,label=paste("p=",as.numeric(as.character(g3[,2])),sep=""))+
    annotate("text",x=g4[,1], y=0.98, colour="green",size=2,label=paste("p=",as.numeric(as.character(g4[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/BRCA/1/",i,".pdf",sep=""),width = 10,height=9)
}


library(ggplot2)
path="D:/tp53/"
cancer="BRCA"  ##4亚型+normal+mu+nomu
m1=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
m2=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
m3=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
mutation=gsub("\\.","-",colnames(m1)[-1])
nomutation=gsub("\\.","-",colnames(m2)[-1])
normal=gsub("\\.","-",substring(as.character(colnames(m3)[-1]),1,12))
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$breast_carcinoma_estrogen_receptor_status),as.character(f$breast_carcinoma_progesterone_receptor_status),as.character(f$lab_proc_her2_neu_immunohistochemistry_receptor_status))))
a[,1]=as.character(a[,1])
a[,2]=as.character(a[,2])
a[,3]=as.character(a[,3])
a[,4]=as.character(a[,4])
luminalA=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Negative",1] ##和或
luminalB=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Positive",1] ##和或
her2=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Positive",1]
basal=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Negative",1]
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-01","",a1[,1])
d1=a1[as.character(a1[,1]) %in% luminalA,2:7]
d2=a1[as.character(a1[,1]) %in% luminalB,2:7]
d3=a1[as.character(a1[,1]) %in% her2,2:7]
d4=a1[as.character(a1[,1]) %in% basal,2:7]
d5=a1[as.character(a1[,1]) %in% mutation,2:7]
d6=a1[as.character(a1[,1]) %in% nomutation,2:7]
d7=a1[as.character(a1[,1]) %in% normal,2:7]
d11=data.frame()
d21=data.frame()
d31=data.frame()
d41=data.frame()
d51=data.frame()
d61=data.frame()
d71=data.frame()
for(m in 1:ncol(d1)){
  d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
}
for(m in 1:ncol(d2)){
  d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
}
for(m in 1:ncol(d3)){
  d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
}
for(m in 1:ncol(d4)){
  d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
}
for(m in 1:ncol(d5)){
  d51=rbind(d51,cbind(rep(colnames(d5)[m],nrow(d5)),d5[,m]))
}
for(m in 1:ncol(d6)){
  d61=rbind(d61,cbind(rep(colnames(d6)[m],nrow(d6)),d6[,m]))
}
for(m in 1:ncol(d7)){
  d71=rbind(d71,cbind(rep(colnames(d7)[m],nrow(d7)),d7[,m]))
}
if(nrow(d11)!=0){
  d12=cbind(cancer,state="luminalA",d11)
}
if(nrow(d21)!=0){
  d22=cbind(cancer,state="luminalB",d21)
}
if(nrow(d31)!=0){
  d32=cbind(cancer,state="her2",d31)
}
if(nrow(d41)!=0){
  d42=cbind(cancer,state="basal",d41)
}
if(nrow(d51)!=0){
  d52=cbind(cancer,state="mutation",d51)
}
if(nrow(d61)!=0){
  d62=cbind(cancer,state="nomutation",d61)
}
if(nrow(d71)!=0){
  d72=cbind(cancer,state="normal",d71)
}
g=data.frame()
if(nrow(d12)!=0){
  g=rbind(g,d12)
}
if(nrow(d22)!=0){
  g=rbind(g,d22)
}
if(nrow(d32)!=0){
  g=rbind(g,d32)
}
if(nrow(d42)!=0){
  g=rbind(g,d42)
}
if(nrow(d52)!=0){
  g=rbind(g,d52)
}
if(nrow(d62)!=0){
  g=rbind(g,d62)
}
if(nrow(d72)!=0){
  g=rbind(g,d72)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #basal_luminalA$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #basal_luminalB$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g4=data.frame()  #basal_her2$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])>=2){
      g4=rbind(g4,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])))$p.value))
    }else{
      g4=rbind(g4,cbind(j,-5))
    }
  }
  g5=data.frame()  #luminalA_luminalB$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])>=2){
      g5=rbind(g5,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])))$p.value))
    }else{
      g5=rbind(g5,cbind(j,-5))
    }
  }
  g6=data.frame()  #luminalA_her2$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])>=2){
      g6=rbind(g6,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])))$p.value))
    }else{
      g6=rbind(g6,cbind(j,-5))
    }
  }
  g7=data.frame()  #luminalB_her2$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])>=2){
      g7=rbind(g7,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2",3])))$p.value))
    }else{
      g7=rbind(g7,cbind(j,-5))
    }
  }
  g8=data.frame()  #mutation_normal$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])>=2){
      g8=rbind(g8,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])))$p.value))
    }else{
      g8=rbind(g8,cbind(j,-5))
    }
  }
  g9=data.frame()  #nomutation_normal$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])>=2){
      g9=rbind(g9,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])))$p.value))
    }else{
      g9=rbind(g9,cbind(j,-5))
    }
  }
  g10=data.frame()  #mutation_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])>=2){
      g10=rbind(g10,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])))$p.value))
    }else{
      g10=rbind(g10,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","blue","green","grey","pink","lightgreen","orange"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=-Inf, y=-Inf,hjust =-2,vjust =-5,colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))+
    annotate("text",x=-Inf, y=-Inf,hjust = -3.5,vjust = -5,colour="blue",size=2,label=paste("p=",as.numeric(as.character(g3[,2])),sep=""))+
    annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5,colour="green",size=2,label=paste("p=",as.numeric(as.character(g4[,2])),sep=""))+
    annotate("text",x=g5[,1], y=0.7, colour="blue",size=2,label=paste("p=",as.numeric(as.character(g5[,2])),sep=""))+
    annotate("text",x=g6[,1], y=0.85, colour="black",size=2,label=paste("p=",as.numeric(as.character(g6[,2])),sep=""))+
    annotate("text",x=g7[,1], y=1, colour="blue",size=2,label=paste("p=",as.numeric(as.character(g7[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust=4,vjust=5,colour="pink",size=2,label=paste("p=",as.numeric(as.character(g8[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust=6,vjust=5,colour="lightgreen",size=2,label=paste("p=",as.numeric(as.character(g9[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust=8,vjust=5,colour="black",size=2,label=paste("p=",as.numeric(as.character(g10[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/BRCA/2/",i,".pdf",sep=""),width = 10,height=9)
}



library(ggplot2)
path="D:/tp53/"
cancer="BRCA"  ##4亚型+p53\WT
m1=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
m2=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
m3=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
mutation=gsub("\\.","-",colnames(m1)[-1])
nomutation=gsub("\\.","-",colnames(m2)[-1])
#normal=gsub("\\.","-",substring(as.character(colnames(m3)[-1]),1,12))
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$breast_carcinoma_estrogen_receptor_status),as.character(f$breast_carcinoma_progesterone_receptor_status),as.character(f$lab_proc_her2_neu_immunohistochemistry_receptor_status))))
a[,1]=as.character(a[,1])
a[,2]=as.character(a[,2])
a[,3]=as.character(a[,3])
a[,4]=as.character(a[,4])
luminalA=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Negative",1] ##和或
luminalB=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Positive",1] ##和或
her2=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Positive",1]
basal=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Negative",1]
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-01","",a1[,1])
d1=a1[as.character(a1[,1]) %in% intersect(luminalA,mutation),2:7]
d2=a1[as.character(a1[,1]) %in% intersect(luminalB,mutation),2:7]
d3=a1[as.character(a1[,1]) %in% intersect(her2,mutation),2:7]
d4=a1[as.character(a1[,1]) %in% intersect(basal,mutation),2:7]
d5=a1[as.character(a1[,1]) %in% intersect(luminalA,nomutation),2:7]
d6=a1[as.character(a1[,1]) %in% intersect(luminalB,nomutation),2:7]
d7=a1[as.character(a1[,1]) %in% intersect(her2,nomutation),2:7]
d8=a1[as.character(a1[,1]) %in% intersect(basal,nomutation),2:7]
d11=data.frame()
d21=data.frame()
d31=data.frame()
d41=data.frame()
d51=data.frame()
d61=data.frame()
d71=data.frame()
d81=data.frame()
for(m in 1:ncol(d1)){
  d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
}
for(m in 1:ncol(d2)){
  d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
}
for(m in 1:ncol(d3)){
  d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
}
for(m in 1:ncol(d4)){
  d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
}
for(m in 1:ncol(d5)){
  d51=rbind(d51,cbind(rep(colnames(d5)[m],nrow(d5)),d5[,m]))
}
for(m in 1:ncol(d6)){
  d61=rbind(d61,cbind(rep(colnames(d6)[m],nrow(d6)),d6[,m]))
}
for(m in 1:ncol(d7)){
  d71=rbind(d71,cbind(rep(colnames(d7)[m],nrow(d7)),d7[,m]))
}
for(m in 1:ncol(d8)){
  d81=rbind(d81,cbind(rep(colnames(d8)[m],nrow(d8)),d8[,m]))
}
if(nrow(d11)!=0){
  d12=cbind(cancer,state="luminalA_mutation",d11)
}
if(nrow(d21)!=0){
  d22=cbind(cancer,state="luminalB_mutation",d21)
}
if(nrow(d31)!=0){
  d32=cbind(cancer,state="her2_mutation",d31)
}
if(nrow(d41)!=0){
  d42=cbind(cancer,state="basal_mutation",d41)
}
if(nrow(d51)!=0){
  d52=cbind(cancer,state="luminalA_nomutation",d51)
}
if(nrow(d61)!=0){
  d62=cbind(cancer,state="luminalB_nomutation",d61)
}
if(nrow(d71)!=0){
  d72=cbind(cancer,state="her2_nomutation",d71)
}
if(nrow(d81)!=0){
  d82=cbind(cancer,state="basal_nomutation",d81)
}
g=data.frame()
if(nrow(d12)!=0){
  g=rbind(g,d12)
}
if(nrow(d52)!=0){
  g=rbind(g,d52)
}
if(nrow(d22)!=0){
  g=rbind(g,d22)
}
if(nrow(d62)!=0){
  g=rbind(g,d62)
}
if(nrow(d32)!=0){
  g=rbind(g,d32)
}
if(nrow(d72)!=0){
  g=rbind(g,d72)
}
if(nrow(d42)!=0){
  g=rbind(g,d42)
}
if(nrow(d82)!=0){
  g=rbind(g,d82)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #luminalA_mutation\luminalA_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA_mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA_nomutation",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA_mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalA_nomutation",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #luminalB_mutation\luminalB_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB_mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB_nomutation",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB_mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="luminalB_nomutation",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g4=data.frame()  #her2_mutation\her2_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2_mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2_nomutation",3])>=2){
      g4=rbind(g4,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2_mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="her2_nomutation",3])))$p.value))
    }else{
      g4=rbind(g4,cbind(j,-5))
    }
  }
  g5=data.frame()  #basal_mutation\basal_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_nomutation",3])>=2){
      g5=rbind(g5,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_nomutation",3])))$p.value))
    }else{
      g5=rbind(g5,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","pink","blue","lightblue","green","lightgreen","black","grey"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=-Inf, y=-Inf,hjust =-2,vjust =-5,colour="red",size=2,label=paste(as.numeric(as.character(g2[,2])),nrow(g1[g1[,2]=="luminalA_mutation",]),nrow(g1[g1[,2]=="luminalA_nomutation",]),sep="_"))+
    annotate("text",x=-Inf, y=-Inf,hjust = -3.5,vjust = -5,colour="blue",size=2,label=paste(as.numeric(as.character(g3[,2])),nrow(g1[g1[,2]=="luminalB_mutation",]),nrow(g1[g1[,2]=="luminalB_nomutation",]),sep="_"))+
    annotate("text",x=Inf, y=Inf,hjust=4,vjust=5,colour="green",size=2,label=paste(as.numeric(as.character(g4[,2])),nrow(g1[g1[,2]=="her2_mutation",]),nrow(g1[g1[,2]=="her2_nomutation",]),sep="_"))+
    annotate("text",x=Inf, y=Inf,hjust=6,vjust=5,colour="black",size=2,label=paste(as.numeric(as.character(g5[,2])),nrow(g1[g1[,2]=="basal_mutation",]),nrow(g1[g1[,2]=="basal_nomutation",]),sep="_"))
    ggsave(file=paste(path,"AGPall_boxplot/BRCA/3/",i,".pdf",sep=""),width = 10,height=9)
}
 
x1=as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_mutation",3]))
x2=as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="basal_nomutation",3]))
da1=cbind(1:length(x1),x1,1)
da2=cbind(1:length(x2),x2,2)
da=rbind(da1,da2)
ggplot(data.frame(da),aes(x=da[,2],fill=factor(da[,3])))+geom_density()
t.test(x1,x2,alternative ="greater")
t.test(x1,x2,alternative ="two.sided")
t.test(x1,x2,alternative ="less")
ks.test(x1,x2)

write.table(hpvneg,file =paste(path,"hpvneg.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

library(ggplot2)
path="D:/tp53/"
cancer="HNSC" ##2亚型hpvpos\hpvneg
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a1=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$hpv_status_by_p16_testing))))
a1[,1]=as.character(a1[,1])
a1[,2]=as.character(a1[,2])
hpvpos=a1[a1[,2]=="Positive",1]
hpvneg=a1[a1[,2]=="Negative",1]
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-01","",a1[,1])
d1=a1[as.character(a1[,1]) %in% hpvpos,2:7]
d2=a1[as.character(a1[,1]) %in% hpvneg,2:7]
d11=data.frame()
d21=data.frame()
for(m in 1:ncol(d1)){
  d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
}
for(m in 1:ncol(d2)){
  d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
}
if(nrow(d11)!=0){
  d12=cbind(cancer,state="hpvpos",d11)
}
if(nrow(d21)!=0){
  d22=cbind(cancer,state="hpvneg",d21)
}
g=data.frame()
if(nrow(d12)!=0){
  g=rbind(g,d12)
}
if(nrow(d22)!=0){
  g=rbind(g,d22)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #hpvpos_hpvneg$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvpos",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvneg",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvpos",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvneg",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","blue"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=g2[,1], y=0.7, colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/HNSC/1/",i,".pdf",sep=""),width = 10,height=9)
}



library(ggplot2)
path="D:/tp53/"
cancer="HNSC" ##2亚型+normal+mu+nomu
m1=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
m2=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
m3=read.csv(paste(path,"Gene Expression/Gene_FPKM1_Analysis1/p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
mutation=gsub("\\.","-",colnames(m1)[-1])
nomutation=gsub("\\.","-",colnames(m2)[-1])
normal=gsub("\\.","-",substring(as.character(colnames(m3)[-1]),1,12))
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a1=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$hpv_status_by_p16_testing))))
a1[,1]=as.character(a1[,1])
a1[,2]=as.character(a1[,2])
hpvpos=a1[a1[,2]=="Positive",1]
hpvneg=a1[a1[,2]=="Negative",1]
f1= read.csv(paste(path,"immuneEstimation.txt",sep=""),header=T,sep="\t")
f1[,1]=as.character(f1[,1])
a1=f1[grep(".*-01$",f1[,1]),]
a1[,1]=gsub("-01","",a1[,1])
d1=a1[as.character(a1[,1]) %in% hpvpos,2:7]
d2=a1[as.character(a1[,1]) %in% hpvneg,2:7]
d3=a1[as.character(a1[,1]) %in% normal,2:7]
d4=a1[as.character(a1[,1]) %in% mutation,2:7]
d5=a1[as.character(a1[,1]) %in% nomutation,2:7]
d11=data.frame()
d21=data.frame()
d31=data.frame()
d41=data.frame()
d51=data.frame()
for(m in 1:ncol(d1)){
  d11=rbind(d11,cbind(rep(colnames(d1)[m],nrow(d1)),d1[,m]))
}
for(m in 1:ncol(d2)){
  d21=rbind(d21,cbind(rep(colnames(d2)[m],nrow(d2)),d2[,m]))
}
for(m in 1:ncol(d3)){
  d31=rbind(d31,cbind(rep(colnames(d3)[m],nrow(d3)),d3[,m]))
}
for(m in 1:ncol(d4)){
  d41=rbind(d41,cbind(rep(colnames(d4)[m],nrow(d4)),d4[,m]))
}
for(m in 1:ncol(d5)){
  d51=rbind(d51,cbind(rep(colnames(d5)[m],nrow(d5)),d5[,m]))
}
if(nrow(d11)!=0){
  d12=cbind(cancer,state="hpvpos",d11)
}
if(nrow(d21)!=0){
  d22=cbind(cancer,state="hpvneg",d21)
}
if(nrow(d31)!=0){
  d32=cbind(cancer,state="normal",d31)
}
if(nrow(d41)!=0){
  d42=cbind(cancer,state="mutation",d41)
}
if(nrow(d51)!=0){
  d52=cbind(cancer,state="nomutation",d51)
}
g=data.frame()
if(nrow(d12)!=0){
  g=rbind(g,d12)
}
if(nrow(d22)!=0){
  g=rbind(g,d22)
}
if(nrow(d32)!=0){
  g=rbind(g,d32)
}
if(nrow(d42)!=0){
  g=rbind(g,d42)
}
if(nrow(d52)!=0){
  g=rbind(g,d52)
}
for(i in c("B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")){
  g1=g[as.character(g[,3])==i,c(1,2,4)]
  g2=data.frame()  #hpvpos_hpvneg$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvpos",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvneg",3])>=2){
      g2=rbind(g2,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvpos",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="hpvneg",3])))$p.value))
    }else{
      g2=rbind(g2,cbind(j,-5))
    }
  }
  g3=data.frame()  #normal_mutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])>=2){
      g3=rbind(g3,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])))$p.value))
    }else{
      g3=rbind(g3,cbind(j,-5))
    }
  }
  g4=data.frame()  #normal_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])>=2){
      g4=rbind(g4,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="normal",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])))$p.value))
    }else{
      g4=rbind(g4,cbind(j,-5))
    }
  }
  g5=data.frame()  #mutation_nomutation$p.value
  for(j in as.character(unique(g1[,1]))){
    if(length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])>=2 & length(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])>=2){
      g5=rbind(g5,cbind(j,t.test(as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="mutation",3])),as.numeric(as.character(g1[as.character(g1[,1])==j & as.character(g1[,2])=="nomutation",3])))$p.value))
    }else{
      g5=rbind(g5,cbind(j,-5))
    }
  }
  g1[,1]=as.character(g1[,1])
  g1[,3]=as.numeric(as.character(g1[,3]))
  p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
    scale_fill_manual(values=c("red","blue","grey","pink","lightgreen"))+coord_flip()
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Cancer")+ylab("purity")+ggtitle(i)+ylim(0,1.1) 
  p <- p+annotate("text",x=-Inf, y=-Inf,hjust = -3.5,vjust = -5,colour="red",size=2,label=paste("p=",as.numeric(as.character(g2[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5,colour="pink",size=2,label=paste("p=",as.numeric(as.character(g3[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 4,vjust = 5, colour="lightgreen",size=2,label=paste("p=",as.numeric(as.character(g4[,2])),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 6,vjust = 5, colour="black",size=2,label=paste("p=",as.numeric(as.character(g5[,2])),sep=""))
  ggsave(file=paste(path,"AGPall_boxplot/HNSC/2/",i,".pdf",sep=""),width = 10,height=9)
}





####五、对于基因上下文的功能的调控机制
library(ggplot2)
path="D:/tp53/"
cancer="CESC"
f1= read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
f1[,2]<-as.numeric(f1[,2])
f1[sapply(f1,is.infinite)]<-NA
f1[sapply(f1,is.na)]<-0
f1<-f1[-1,]

f= read.csv(paste(path,"regulate.txt",sep=""),header=T,sep="\t")
for( i in as.character(unique(f[,1]))){ #i="Core" 
  a=f[as.character(f[,1])==i,]
  for(j in 1:nrow(a)){ # j=8
    x=f1[f1[,1]==as.character(a[j,2]),2:ncol(f1)]
    y=f1[f1[,1]==as.character(a[j,3]),2:ncol(f1)]
    if(nrow(y)!=0 & nrow(y)!=0){
      data=data.frame(cbind(t(x),t(y)))
      ggplot(data,aes(x = data[,1],y =data[,2]))+geom_point(size =0.2,aes(colour = "red"))
      p<-ggplot(data=g1, aes(x = g1[,1],y = g1[,3]))+geom_boxplot(aes(fill=state),outlier.colour = NA)+
        scale_fill_manual(values=c("red","green"))+coord_flip()
      p <- p+theme(plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
      p <- p + xlab("Cancer") + ylab("purity") + ggtitle(paste("Immune Cell Types=",i,sep="")) #添加横纵坐标，添加title
      p <- p+ theme(legend.position="top")
      ggsave(file=paste(path,"AGPall_boxplot/1/",i,".pdf",sep=""),width = 10,height=9)
      
    }
    
  }
}








####六、对19个cancer的lncRNA的数据处理,结果去了LGG9、OV13没有normal样本的17cancer
path="F:/tp53/"
f= read.csv(paste(path,"lncRNA/Geneid_symbol_type.txt",sep=""),header=T,sep="\t")
f[,1]=as.character(f[,1])
f[,2]=as.character(f[,2])
filename<-list.files(paste(path,"lncRNA/",sep=""),pattern="^TCGA-")
for(i in filename[14:19]){
  cancer=gsub("TCGA-|-rnaexpr","",i)
  f1= read.csv(paste(path,"lncRNA/",i,"/",i,".tsv",sep=""),header=T,sep="\t")
  f1[,1]=gsub("\\..*","",as.character(f1[,1]))
  a=merge(f,f1,by.x="Gene.stable.ID",by.y="Gene_ID")
  e0=cbind(gene_Name=a[,2],a[,grep(paste(cancer,".Normal.",sep=""),colnames(a))])
  write.table(e0,file =paste(path,"lncRNA/normal/",cancer,".txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  f2<-read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-c(strsplit(as.character(f2[1,1])," ")[[1]],strsplit(as.character(f2[1,2])," ")[[1]])
  nomutation=paste(cancer,".Tumor.",gsub("_no","",colnames(f2)[grep("_no",colnames(f2))]),sep="")
  mutation=paste(cancer,".Tumor.",colnames(f2)[-grep("_no",colnames(f2))][-1],sep="")
  e1=cbind(gene_Name=a[,2],a[,intersect(colnames(a),mutation)])
  e2=cbind(gene_Name=a[,2],a[,intersect(colnames(a),nomutation)])
  write.table(e1,file =paste(path,"lncRNA/mutation/",cancer,".txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
  write.table(e2,file =paste(path,"lncRNA/nomutation/",cancer,".txt",sep=""), row.names = F,col.names=T, quote = F,sep = "\t",append=TRUE)
}

#lncRNA foldchange值=p53/normal_mean、WT/normal_mean取log2
#"COAD","GBM","READ","SKCM",没有normal样本，17-4=13cancer
path="F:/tp53/lncRNA/"
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PRAD","STAD","THCA","UCEC")){
  f1<-read.table(paste(path,"normal/",cancer,".txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"mutation/",cancer,".txt",sep=""),header=T,sep="\t")
  f3<-read.table(paste(path,"nomutation/",cancer,".txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(gsub(".Tumor.","",gsub(cancer,"",colnames(f2))),"\t",file=paste(path,"foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
  cat(paste(gsub(".Tumor.","",gsub(cancer,"",colnames(f3)))[-1],"_no",sep=""),"\n",file=paste(path,"foldchange/",cancer,"_fc.txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    n=log2(f3[i,2:ncol(f3)]/k[i])
    write.table(append(m,n),file=paste(path,"foldchange/",cancer,"_fc.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}


##对13cancer的lncRNA的foldchange值=p53/normal_mean、WT/normal_mean取log2扩充，colnames和GeneExpression的foldchange值一致
path="F:/tp53/lncRNA/"
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PRAD","STAD","THCA","UCEC")){
  f1<-read.csv(paste(path,"foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  f1[,1]<-gsub("LINC-","lincRNA-",f1[,1])
  f2<-read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-c(strsplit(as.character(f2[1,1])," ")[[1]],strsplit(as.character(f2[1,2])," ")[[1]])
  d=data.frame(matrix(0,nrow(f1),ncol(f2)))
  colnames(d)=colnames(f2)
  d$gene_Name=f1$gene_Name
  d[,colnames(f1)]=f1[,colnames(f1)]
  write.table(d,file=paste(path,"foldchange1/",cancer,"_fc.txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

#lncRNA的foldchange值=p53/WT_mean取log2
path="F:/tp53/lncRNA/"
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PRAD","STAD","THCA","UCEC")){
  f1<-read.table(paste(path,"nomutation/",cancer,".txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"mutation/",cancer,".txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(paste(gsub(".Tumor.","",gsub(cancer,"",colnames(f2))),"_p53/WT_mean",sep=""),"\n",file=paste(path,"foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    write.table(m,file=paste(path,"foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}


##对13cancer的lncRNA的foldchange值=p53/WT_mean取log2扩充，colnames和GeneExpression的foldchange值一致
path="F:/tp53/lncRNA/"
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PRAD","STAD","THCA","UCEC")){
  f1<-read.csv(paste(path,"foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  f1[,1]<-gsub("LINC-","lincRNA-",f1[,1])
  f2<-read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2<-f2[-1,]
  d=data.frame(matrix(0,nrow(f1),ncol(f2)))
  colnames(d)=colnames(f2)
  d[,1]=f1[,1]
  d[,colnames(f1)]=f1[,colnames(f1)]
  write.table(d,file=paste(path,"foldchange_p53%WTmean1/",cancer,"_fc.txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

#############################没做
##对13cancer的lncRNA的cor\foldchange值=p53/normal_mean、WT/normal_mean、p53/WT_mean做热图pheatmap
library(pheatmap)
library(gplots)
path="F:/tp53/lncRNA/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"foldchange1/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"foldchange1/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$gene_Name
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),]
  
  f2<-read.csv(paste(path,"foldchange_p53%WTmean1/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  
  e1<-read.table(paste(path,"normal/",gsub("_fc","",i),sep=""),header=T,sep="\t")
  e2<-read.table(paste(path,"mutation/",gsub("_fc","",i),sep=""),header=T,sep="\t")
  e3<-read.table(paste(path,"nomutation/",gsub("_fc","",i),sep=""),header=T,sep="\t")
  e1[,1]<-gsub("LINC-","lincRNA-",e1[,1])
  e2[,1]<-gsub("LINC-","lincRNA-",e2[,1])
  e3[,1]<-gsub("LINC-","lincRNA-",e3[,1])
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
           file=paste(path,"foldchange_heatmap/2heatmap/",gsub("_fc.txt","_cor0.5_p53_WT.pdf",i),sep=""),width=8,height=8,main=gsub("_fc.txt","",i))  
}
#########################################################



####六、对27cancer的miRNA的数据处理
path="F:/tp53/miRNA/"
f<-read.table("D:/tp53/Gene Expression/3894_5886=9780.txt",header=T,sep="\t")
filename<-list.files(paste(path,"miRNA_fpkm/",sep=""))
for(i in filename){
  f1<-read.table(paste(path,"miRNA_fpkm/",i,sep=""),header=T,sep="\t")
  colnames(f1)=substring(as.character(colnames(f1)),1,15)
  a1<-f1[,c(1,grep(".*.11$",colnames(f1)))]
  mutation<-subset(f,as.character(f[,2])==1 & as.character(f[,5])==gsub("_miRNA_fpkm.txt","",i))[,1]
  nomutation<-subset(f,is.na(as.character(f[,2]))& as.character(f[,5])==gsub("_miRNA_fpkm.txt","",i))[,1]
  colnames(f1)<-gsub("\\.","-",substring(colnames(f1),1,12))
  a2<-f1[,c("miRNA_ID",intersect(as.character(mutation),colnames(f1)))]
  a3<-f1[,c("miRNA_ID",intersect(as.character(nomutation),colnames(f1)))]
  write.table(a1,file=paste(path,"normal/",gsub("_miRNA_fpkm","",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a2,file=paste(path,"mutation/",gsub("_miRNA_fpkm","",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a3,file=paste(path,"nomutation/",gsub("_miRNA_fpkm","",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

#27cancer foldchange值=p53/normal_mean、WT/normal_mean取log2==20cancer
#6个no normal:"LAML","LGG","OV","SARC","TGCT","UCS",
#1个no mutation："GBM",
path="F:/tp53/miRNA/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","HNSC","KIRC","KIRP",
                "LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","STAD",
                "THCA","THYM","UCEC","SKCM")){
  f1<-read.table(paste(path,"normal/",cancer,".txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"mutation/",cancer,".txt",sep=""),header=T,sep="\t")
  f3<-read.table(paste(path,"nomutation/",cancer,".txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(colnames(f2),"\t",file=paste(path,"foldchange/",cancer,".txt",sep=""),append=TRUE)
  cat(paste(colnames(f3)[-1],"_no",sep=""),"\n",file=paste(path,"foldchange/",cancer,".txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    n=log2(f3[i,2:ncol(f3)]/k[i])
    write.table(append(m,n),file=paste(path,"foldchange/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}

#27cancer foldchange值=p53/WT_mean取log2==26cancer
#1个no mutation："GBM",
path="F:/tp53/miRNA/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC","LAML","LGG","OV","TGCT","UCS","SKCM")){
  f1<-read.table(paste(path,"nomutation/",cancer,".txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"mutation/",cancer,".txt",sep=""),header=T,sep="\t")
  k=apply(f1[,2:ncol(f1)],1,mean)
  cat(paste(colnames(f2),"_p53/WT_mean",sep=""),"\n",file=paste(path,"foldchange_p53%WTmean/",cancer,".txt",sep=""),append=TRUE)
  for(i in 1:length(k)){
    m=cbind(as.character(f2[i,1]),log2(f2[i,2:ncol(f2)]/k[i]))
    write.table(m,file=paste(path,"foldchange_p53%WTmean/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}

##对20cancer的miRNA的foldchange值=p53/normal_mean、WT/normal_mean取log2扩充，colnames和GeneExpression的foldchange值一致
path="F:/tp53/miRNA/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","HNSC","KIRC","KIRP",
                "LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","READ","STAD",
                "THCA","THYM","UCEC","SKCM")){
  f1<-read.csv(paste(path,"foldchange/",cancer,".txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  f1[,1]<-gsub("mir","miR",gsub("hsa-","",f1[,1]))
  f2<-read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-c(strsplit(as.character(f2[1,1])," ")[[1]],strsplit(as.character(f2[1,2])," ")[[1]])
  d=data.frame(matrix(0,nrow(f1),ncol(f2)))
  colnames(d)=gsub("gene_Name","miRNA_ID",colnames(f2))
  d[,1]=f1[,1]
  d[,intersect(colnames(f1)[-1],colnames(f2)[-1])]=f1[,intersect(colnames(f1)[-1],colnames(f2)[-1])]
  write.table(d,file=paste(path,"foldchange1/",cancer,".txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

##对26cancer的miRNA的foldchange值=p53/WT_mean取log2扩充，colnames和GeneExpression的foldchange值一致
path="F:/tp53/miRNA/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC","LAML","LGG","OV","TGCT","UCS","SKCM")){
  f1<-read.csv(paste(path,"foldchange_p53%WTmean/",cancer,".txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  f1[,1]<-gsub("mir","miR",gsub("hsa-","",f1[,1]))
  f2<-read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2<-f2[-1,]
  d=data.frame(matrix(0,nrow(f1),ncol(f2)))
  colnames(d)=colnames(f2)
  d[,1]=f1[,1]
  d[,intersect(colnames(f1)[-1],colnames(f2)[-1])]=f1[,intersect(colnames(f1)[-1],colnames(f2)[-1])]
  write.table(d,file=paste(path,"foldchange_p53%WTmean1/",cancer,".txt",sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

###############################没做
##对20cancer的miRNA的cor\foldchange值=p53/normal_mean、WT/normal_mean、p53/WT_mean做热图pheatmap
library(pheatmap)
library(gplots)
path="F:/tp53/miRNA/"
colorsChoice<- colorRampPalette(c("green","white","red"))
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"foldchange1/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"foldchange1/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f1)<-f1$miRNA_ID
  f11=f1[intersect(as.character(f[,1]),as.character(f1$miRNA_ID)),]
  
  f2<-read.csv(paste(path,"foldchange_p53%WTmean1/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  
  e1<-read.table(paste(path,"normal/",i,sep=""),header=T,sep="\t")
  e2<-read.table(paste(path,"mutation/",i,sep=""),header=T,sep="\t")
  e3<-read.table(paste(path,"nomutation/",i,sep=""),header=T,sep="\t")
  e1[,1]<-gsub("mir","miR",gsub("hsa-","",e1[,1]))
  e2[,1]<-gsub("mir","miR",gsub("hsa-","",e2[,1]))
  e3[,1]<-gsub("mir","miR",gsub("hsa-","",e3[,1]))
  e4<-cbind(e2,e3)
  h1<-read.table(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_normal_fpkm/",gsub(".txt","_p53_context_normal_fpkm.txt",i),sep=""),header=T,sep="\t")
  h2<-read.table(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_mutation_fpkm/",gsub(".txt","_p53_context_mutation_fpkm.txt",i),sep=""),header=T,sep="\t")
  h3<-read.table(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_nomutation_fpkm/",gsub(".txt","_p53_context_nomutation_fpkm.txt",i),sep=""),header=T,sep="\t")
  h4<-cbind(h2,h3)
  a1=h1[grep("^TP53$",h1[,1]),2:ncol(h1)]
  a2=h2[grep("^TP53$",h2[,1]),2:ncol(h2)]
  a3=h3[grep("^TP53$",h3[,1]),2:ncol(h3)]
  a4=h4[grep("^TP53$",h4[,1]),2:ncol(h4)]
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
  colnames(p)[201]="miRNA_ID"
  rownames(p)=p[,201]
  p<-na.omit(p)
  f33=p[intersect(as.character(f[,1]),as.character(p[,201])),]
  
  d=merge(f11,f22,by.x="miRNA_ID",by.y="gene_Name_p53.WT_mean")
  rownames(d)<-d[,1]
  dd=merge(f33,d,by="miRNA_ID")
  rownames(dd)<-dd[,1]
  d1=dd[intersect(as.character(f[,1]),as.character(dd[,1])),]
  d1[sapply(d1,is.na)]<-0
  x<-d1[,2:ncol(d1)]
  y<-as.matrix(x)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  col_anno = data.frame(SampleClass = factor(rep(c("cor_normal","cor_Mutation", "cor_NoMutation","cor_cancer","Mutation", "NoMutation","Mutation/NoMutation"), c(50,50,50,50,ncol(f1)-length(grep("_no",colnames(f1)))-1, length(grep("_no",colnames(f1))),ncol(f2)-1))),row.names=colnames(dd)[2:ncol(dd)])  
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-1,1,max(y)),breaks=c(min(y),-1,1,max(y)),annotation_col = col_anno, annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"foldchange_heatmap/2heatmap/",gsub(".txt","_cor0.5_p53_WT.pdf",i),sep=""),width=8,height=8,main=gsub(".txt","",i))  
}
#########################################################






####七、对BRCA、HNSC分亚型、SKCM等其他25cancer分"Primary","Metastasis"亚型，统计对应的突变数、tp53突变数，突变率等
path="D:/tp53/"
ff=read.csv(paste(path,"3894_5886=9780.txt",sep=""),header=T,sep="\t")
rownames(ff)=as.character(ff[,1])
cancer="BRCA"
f=read.csv(paste(path,"phenotype/TCGA-",cancer,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a=data.frame(unique(cbind(substring(as.character(f$submitter_id.samples),1,12),as.character(f$breast_carcinoma_estrogen_receptor_status),as.character(f$breast_carcinoma_progesterone_receptor_status),as.character(f$lab_proc_her2_neu_immunohistochemistry_receptor_status))))
a[,1]=as.character(a[,1])
a[,2]=as.character(a[,2])
a[,3]=as.character(a[,3])
a[,4]=as.character(a[,4])
luminalA=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Negative",1] ##和或
luminalB=a[!(a[,2]=="Negative"& a[,3]=="Negative")&a[,4]=="Positive",1] ##和或
her2=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Positive",1]
basal=a[a[,2]=="Negative"& a[,3]=="Negative"&a[,4]=="Negative",1]
g=t(rbind(class=c("luminalA","luminalB","her2","basal"),all=c(length(luminalA),length(luminalB),length(her2),length(basal))))
g=data.frame(g)
g[,1]=as.character(g[,1])
g[,2]=as.numeric(as.character(g[,2]))
g$Mutation[1]=length(intersect(as.character(ff[,1]),luminalA))
g$Mutation[2]=length(intersect(as.character(ff[,1]),luminalB))
g$Mutation[3]=length(intersect(as.character(ff[,1]),her2))
g$Mutation[4]=length(intersect(as.character(ff[,1]),basal))
g$TP53Mutation[1]=length(na.omit(ff[intersect(as.character(ff[,1]),luminalA),2]))
g$TP53Mutation[2]=length(na.omit(ff[intersect(as.character(ff[,1]),luminalB),2]))
g$TP53Mutation[3]=length(na.omit(ff[intersect(as.character(ff[,1]),her2),2]))
g$TP53Mutation[4]=length(na.omit(ff[intersect(as.character(ff[,1]),basal),2]))

cancer1="HNSC"
f1=read.csv(paste(path,"phenotype/TCGA-",cancer1,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
a1=data.frame(unique(cbind(substring(as.character(f1$submitter_id.samples),1,12),as.character(f1$hpv_status_by_p16_testing))))
a1[,1]=as.character(a1[,1])
a1[,2]=as.character(a1[,2])
hpvpos=a1[a1[,2]=="Positive",1]
hpvneg=a1[a1[,2]=="Negative",1]
g1=t(rbind(class=c("hpvpos","hpvneg"),all=c(length(hpvpos),length(hpvneg))))
g1=data.frame(g1)
g1[,1]=as.character(g1[,1])
g1[,2]=as.numeric(as.character(g1[,2]))
g1$Mutation[1]=length(intersect(as.character(ff[,1]),hpvpos))
g1$Mutation[2]=length(intersect(as.character(ff[,1]),hpvneg))
g1$TP53Mutation[1]=length(na.omit(ff[intersect(as.character(ff[,1]),hpvpos),2]))
g1$TP53Mutation[2]=length(na.omit(ff[intersect(as.character(ff[,1]),hpvneg),2]))

#cancer2="SKCM"
for(cancer2 in c("BLCA","CESC","COAD","ESCA","GBM","KIRC","KIRP","LAML",
"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  f2=read.csv(paste(path,"phenotype/TCGA-",cancer2,".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
  a2=data.frame(cbind(as.character(f2$submitter_id.samples),as.character(f2$sample_type.samples)))
  a2[,1]=as.character(a2[,1])
  a2[,2]=as.character(a2[,2])
  Primary=substring(as.character(a2[a2[,2]=="Primary Tumor",1]),1,12)
  Metastasis=substring(as.character(a2[a2[,2]=="Metastatic",1]),1,12)
  g2=t(rbind(class=c("Primary","Metastasis"),all=c(length(Primary),length(Metastasis))))
  g2=data.frame(g2)
  g2[,1]=as.character(g2[,1])
  g2[,2]=as.numeric(as.character(g2[,2]))
  g2$Mutation[1]=length(intersect(as.character(ff[,1]),Primary))
  g2$Mutation[2]=length(intersect(as.character(ff[,1]),Metastasis))
  g2$TP53Mutation[1]=length(na.omit(ff[intersect(as.character(ff[,1]),Primary),2]))
  g2$TP53Mutation[2]=length(na.omit(ff[intersect(as.character(ff[,1]),Metastasis),2]))
  write.table(cbind(cancer2,g2),file =paste(path,"Primary_Metastasis.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}

x<-matrix(c(492,2,8,0),nrow=2,byrow = T)
d1=fisher.test(x)



####八、对每一cancer，把gene，miRNA，lincRNA的foldchange（Mu/normal,NoMu/normal）整合在一起，行为gene，列为样本。
#然后确定两两gene关系在样本中的关系：如在sample1中gene1上调，gene2下调。
#取p53上下文里的gene，关系对也是p53上下文里的调控关系。
#D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/有22cancer=13+7+"GBM","SARC".
#F:/tp53/miRNA/foldchange1/有20cancer=13+"COAD","ESCA","PAAD","PCPG","READ","SKCM","THYM".
#F:/tp53/lncRNA/foldchange1/有13cancer
###13cancer，mRNA、miRNA、lncRNA都有。
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP",
                "LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  f2= read.table(paste("F:/tp53/miRNA/foldchange1/",cancer,".txt",sep=""),header=T,sep="\t")
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  colnames(f22)<-gsub("miRNA_ID","gene_Name",colnames(f22))
  f3= read.table(paste("F:/tp53/lncRNA/foldchange1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  m=rbind(f11,f22,f33)
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}


cancer="BRCA"
ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
ff[,2]=gsub("^p53$","TP53",ff[,2])
ff[,3]=gsub("^p53$","TP53",ff[,3])
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
f1[,2]<-as.numeric(f1[,2])
f1[sapply(f1,is.infinite)]<-NA
f1[sapply(f1,is.na)]<-0
f1<-f1[-1,]
rownames(f1)<-f1[,1]
f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
f2= read.table(paste("F:/tp53/miRNA/foldchange1/",cancer,".txt",sep=""),header=T,sep="\t")
rownames(f2)<-f2[,1]
f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
colnames(f22)<-gsub("miRNA_ID","gene_Name",colnames(f22))
f3= read.table(paste("F:/tp53/lncRNA/foldchange1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
m=rbind(f11,f22,f33)
cat(colnames(ff),file=paste("F:/tp53/3generegulate/relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
cat("\t",file=paste("F:/tp53/3generegulate/relationship/",cancer,".txt",sep=""),append=TRUE)
cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
for(i in 1:nrow(ff)){
  if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
    for(j in 2:ncol(m)){
      if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
        ff[i,j+3]=1
      }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
        ff[i,j+3]=2
      }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
        ff[i,j+3]=3
      }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
        ff[i,j+3]=4
      }else{
        ff[i,j+3]=0
      }
    }
    write.table(ff[i,],file=paste("F:/tp53/3generegulate/relationship/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}
  
###7cancer，有mRNA、miRNA。
#F:/tp53/miRNA/foldchange1/有20cancer=13+"COAD","ESCA","PAAD","PCPG","READ","SKCM","THYM".
for(cancer in c("COAD","ESCA","PAAD","PCPG","READ","SKCM","THYM")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  f2= read.table(paste("F:/tp53/miRNA/foldchange1/",cancer,".txt",sep=""),header=T,sep="\t")
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  colnames(f22)<-gsub("miRNA_ID","gene_Name",colnames(f22))
  #f3= read.table(paste("F:/tp53/lncRNA/foldchange1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  #f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  #m=rbind(f11,f22,f33)
  m=rbind(f11,f22)
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}


###2cancer只有mRNA。
#D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/有22cancer=13+7+"GBM","SARC".
for(cancer in c("GBM","SARC")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  #f2= read.table(paste("F:/tp53/miRNA/foldchange1/",cancer,".txt",sep=""),header=T,sep="\t")
  #rownames(f2)<-f2[,1]
  #f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  #colnames(f22)<-gsub("miRNA_ID","gene_Name",colnames(f22))
  #f3= read.table(paste("F:/tp53/lncRNA/foldchange1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  #f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  #m=rbind(f11,f22,f33)
  m=f11
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}



##、对每一cancer13，把gene，miRNA，lincRNA的foldchange_p53%WTmean
#（Mu/NoMu）整合在一起，行为gene，列为样本。
###13cancer，mRNA、miRNA、lncRNA都有。
for(cancer in c("BLCA","BRCA","CESC","HNSC","KIRC","KIRP",
                "LIHC","LUAD","LUSC","PRAD","STAD","THCA","UCEC")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  colnames(f11)<-gsub("_p53/WT_mean","_p53.WT_mean",colnames(f11))
  f2= read.table(paste("F:/tp53/miRNA/foldchange_p53%WTmean1/",cancer,".txt",sep=""),header=T,sep="\t")
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  f3= read.table(paste("F:/tp53/lncRNA/foldchange_p53%WTmean1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  m=rbind(f11,f22,f33)
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}

###7cancer，有mRNA、miRNA。
for(cancer in c("COAD","ESCA","PAAD","PCPG","READ","SKCM","THYM")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  colnames(f11)<-gsub("_p53/WT_mean","_p53.WT_mean",colnames(f11))
  f2= read.table(paste("F:/tp53/miRNA/foldchange_p53%WTmean1/",cancer,".txt",sep=""),header=T,sep="\t")
  rownames(f2)<-f2[,1]
  f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  #f3= read.table(paste("F:/tp53/lncRNA/foldchange_p53%WTmean1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  #f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  m=rbind(f11,f22)
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}

###2cancer只有mRNA。
#D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/p53_context_fpkm_foldchange/有22cancer=13+7+"GBM","SARC".
for(cancer in c("GBM","SARC")){
  ff= read.csv("D:/tp53/regulate.txt",header=T,sep="\t")
  ff[,2]=gsub("^p53$","TP53",ff[,2])
  ff[,3]=gsub("^p53$","TP53",ff[,3])
  f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
  f1=read.csv(paste("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_fpkm_foldchange_p53%WTmean/",cancer,"_fc.txt",sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1[,1]
  f11=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  colnames(f11)<-gsub("_p53/WT_mean","_p53.WT_mean",colnames(f11))
  #f2= read.table(paste("F:/tp53/miRNA/foldchange_p53%WTmean1/",cancer,".txt",sep=""),header=T,sep="\t")
  #rownames(f2)<-f2[,1]
  #f22=f2[intersect(as.character(f[,1]),as.character(f2[,1])),]
  #f3= read.table(paste("F:/tp53/lncRNA/foldchange_p53%WTmean1/",cancer,"_fc.txt",sep=""),header=T,sep="\t")
  #f33=f3[grep(intersect(as.character(f[,1]),as.character(f3[,1])),f3[,1]),]
  m=f11
  cat(colnames(ff),file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  cat("\t",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),append=TRUE)
  cat(colnames(m)[2:ncol(m)],"\n",file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),sep="\t",append=TRUE)
  for(i in 1:nrow(ff)){
    if(as.character(ff[i,2]) %in% rownames(m) & as.character(ff[i,3]) %in% rownames(m)){
      for(j in 2:ncol(m)){
        if(m[as.character(ff[i,2]),j] >1 & m[as.character(ff[i,3]),j] >1){
          ff[i,j+3]=1
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=2
        }else if(m[as.character(ff[i,2]),j] > 1 & m[as.character(ff[i,3]),j] < -1){
          ff[i,j+3]=3
        }else if(m[as.character(ff[i,2]),j] < -1 & m[as.character(ff[i,3]),j] > 1){
          ff[i,j+3]=4
        }else{
          ff[i,j+3]=0
        }
      }
      write.table(ff[i,],file=paste("F:/tp53/3generegulate/new_relationship1/",cancer,".txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}



##对F:/tp53/3generegulate/relationship/下的每一文件进行分析
##heatmap
library(pheatmap)
library(gplots)
colorsChoice<- colorRampPalette(c("white","red","green"))
path="F:/tp53/3generegulate/"
filename<-list.files(paste(path,"relationship/",sep=""))
for(i in filename){
  f= read.csv(paste(path,"relationship/",i,sep=""),header=T,sep="\t")
  f=unique(f)
  m=as.matrix(f[,5:ncol(f)])
  row_anno = data.frame(GeneFunction=as.character(f[,1]), row.names=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep="."))
  rownames(m)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
  pheatmap(m,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,
           file=paste(path,"heatmap/",gsub(".txt",".pdf",i),sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  
}

f= read.csv(paste(path,"relationship/",i,sep=""),header=T,sep="\t")
f=unique(f)
m=as.matrix(f[,5:ncol(f)])
row_anno = data.frame(GeneFunction=as.character(f[,1]), row.names=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep="."))
rownames(m)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
pheatmap(m,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,
         file=paste(path,"heatmap/",gsub(".txt",".pdf",i),sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  


##对F:/tp53/3generegulate/relationship/下gene-gene的关系
#对在mutation、nomutation样本中样本
#按 TP53表达值为normal、上下调的顺序排序pheatmap
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("white","red","green"))
filename<-list.files("F:/tp53/3generegulate/new_relationship/")
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",gsub(".txt","_fc.txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1[,-grep("_no",as.character(colnames(f1)))]
  b=f1["TP53",grep("_no",as.character(colnames(f1)))]
  a1=colnames(sort(a[grep("^TP53$",a[,1]),2:ncol(a)]))
  b1=colnames(sort(b))
  
  f= read.csv(paste("F:/tp53/3generegulate/new_relationship/",i,sep=""),header=T,sep="\t")
  f=unique(f)
  f=f[,-grep("^X$",colnames(f))]
  m=f[,5:ncol(f)]
  m1=as.matrix(m[,append(a1,b1)])
  row_anno = data.frame(GeneFunction=as.character(f[,1]), row.names=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep="."))
  rownames(m1)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation", "NoMutation"), c(ncol(f)-length(grep("_no",colnames(f)))-4, length(grep("_no",colnames(f)))))),row.names=colnames(f)[5:ncol(f)])  
  pheatmap(m1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,annotation_col = col_anno, 
           file=paste("F:/tp53/3generegulate/new_heatmap_sort/",gsub(".txt",".pdf",i),sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  
}


##对F:/tp53/3generegulate/relationship/下gene-gene的关系
#对在mutation、nomutation、Mu/nomu样本中样本
#按 TP53表达值为normal、上下调的顺序排序pheatmap,gene关系按字母排序结果为xx1.pdf
##22cancer：PCPG14只有一个Mu样本 delete 共21cancer
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
colorsChoice<- colorRampPalette(c("white","red","green"))
filename<-list.files("F:/tp53/3generegulate/new_relationship/")
for(i in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",gsub(".txt","_fc.txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1[,-grep("_no",as.character(colnames(f1)))]
  b=f1["TP53",grep("_no",as.character(colnames(f1)))]
  a1=colnames(sort(a[grep("^TP53$",a[,1]),2:ncol(a)]))
  b1=colnames(sort(b))
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",gsub(".txt","_fc.txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  c1=gsub("_p53/WT_mean","_p53.WT_mean",colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)])))
  
  f= read.csv(paste("F:/tp53/3generegulate/new_relationship/",i,sep=""),header=T,sep="\t")
  f=unique(f)
  f=f[,-grep("^X$",colnames(f))]
  f=f[order(f$function.,f$gene1,f$gene2),]
  m1=f[,5:ncol(f)]
  rownames(m1)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
  
  ff= read.csv(paste("F:/tp53/3generegulate/new_relationship1/",i,sep=""),header=T,sep="\t")
  ff=unique(ff)
  ff=ff[,-grep("^X$",colnames(ff))]
  ff=ff[order(ff$function.,ff$gene1,ff$gene2),]
  m2=ff[,5:ncol(ff)]
  rownames(m2)=paste(substring(as.character(ff[,1]),0,4),as.character(ff[,2]),as.character(ff[,3]),sep=".")
  
  m3=as.matrix(m1[,append(a1,b1)])
  m4=as.matrix(m2[,c1])
  mm=cbind(m3,m4)
  n=cbind(f1["TP53",a1],f1["TP53",b1],f2["TP53",gsub("_p53.WT_mean","_p53/WT_mean",c1)])
  n[1,which(n[1,]>=1)]=2
  n[1,which(n[1,]<= -1)]=4
  n[1,which(n[1,]<1 & n[1,]> -1)]=0
  colnames(n)=gsub("_p53/WT_mean","_p53.WT_mean",colnames(n))
  row_anno = data.frame(GeneFunction=c("TP53",as.character(f[,1])), row.names=c("TP53",paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")))
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation", "NoMutation","Mu/NoMu"), c(ncol(f)-length(grep("_no",colnames(f)))-4, length(grep("_no",colnames(f))),length(c1)))),row.names=colnames(mm)[1:ncol(mm)])  
  pheatmap(rbind(n,mm),cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,annotation_col = col_anno, 
           file=paste("F:/tp53/3generegulate/new_heatmap_sort2/1/",gsub(".txt","1.pdf",i),sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  
}


#############################################nono
##对F:/tp53/3generegulate/relationship/下gene-gene的关系
#对在mutation、nomutation、Mu/nomu样本中样本
#按 TP53表达值为normal、上下调的样本中突变位点的顺序排序pheatmap,gene关系按字母排序结果为2/xx.pdf
##22cancer：PCPG只有一个Mu样本 delete 共21cancer
library(pheatmap)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
fi= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
fi[,1]=gsub("-",".",fi[,1])
fi$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(fi$HGVSp_Short),4)))

cancer=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD",
         "READ","SARC","SKCM","STAD","THCA","THYM","UCEC")

a=subset(fi,as.character(fi$project_id) %in% cancer)
position=data.frame(a$submitter_id,a$project_id,a$HGVSp_Short)

 

colorsChoice<- colorRampPalette(c("white","red","green"))
filename<-list.files("F:/tp53/3generegulate/relationship/")
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",gsub(".txt","_fc.txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",as.character(colnames(f1)))]
  b=f1["TP53",grep("_no",as.character(colnames(f1)))]
  
  a1=colnames(sort(a[,2:ncol(a)]))
  
  
  colnames(a[1,which(a[1,]<= -1)]) 
  colnames(a[1,which(a[1,]<1 & a[1,]> -1)]) 
  colnames(a[1,which(a[1,]>=1)])[-1]
  
  
  b1=colnames(sort(b))
   
  
  
  
  f2<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",gsub(".txt","_fc.txt",i),sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f2)<-strsplit(as.character(f2[1,1])," ")[[1]]
  f2[sapply(f2,is.infinite)]<-NA
  f2[sapply(f2,is.na)]<-0
  f2<-f2[-1,]
  rownames(f2)<-f2[,1]
  c1=gsub("_p53/WT_mean","_p53.WT_mean",colnames(sort(f2[grep("^TP53$",f2[,1]),2:ncol(f2)])))
  
  f= read.csv(paste("F:/tp53/3generegulate/relationship/",i,sep=""),header=T,sep="\t")
  f=unique(f)
  f=f[,-grep("^X$",colnames(f))]
  f=f[order(f$function.,f$gene1,f$gene2),]
  m1=f[,5:ncol(f)]
  rownames(m1)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
  
  ff= read.csv(paste("F:/tp53/3generegulate/relationship1/",i,sep=""),header=T,sep="\t")
  ff=unique(ff)
  ff=ff[,-grep("^X$",colnames(ff))]
  ff=ff[order(ff$function.,ff$gene1,ff$gene2),]
  m2=ff[,5:ncol(ff)]
  rownames(m2)=paste(substring(as.character(ff[,1]),0,4),as.character(ff[,2]),as.character(ff[,3]),sep=".")
  
  m3=as.matrix(m1[,append(a1,b1)])
  m4=as.matrix(m2[,c1])
  mm=cbind(m3,m4)
  n=cbind(f1["TP53",a1],f1["TP53",b1],f2["TP53",gsub("_p53.WT_mean","_p53/WT_mean",c1)])
  n[1,which(n[1,]>=1)]=2
  n[1,which(n[1,]<= -1)]=4
  n[1,which(n[1,]<1 & n[1,]> -1)]=0
  colnames(n)=gsub("_p53/WT_mean","_p53.WT_mean",colnames(n))
  row_anno = data.frame(GeneFunction=c("TP53",as.character(f[,1])), row.names=c("TP53",paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")))
  col_anno = data.frame(SampleClass = factor(rep(c("Mutation", "NoMutation","Mu/NoMu"), c(ncol(f)-length(grep("_no",colnames(f)))-4, length(grep("_no",colnames(f))),length(c1)))),row.names=colnames(mm)[1:ncol(mm)])  
  pheatmap(rbind(n,mm),cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,annotation_col = col_anno, 
           file=paste("F:/tp53/3generegulate/heatmap_sort2/2/",gsub(".txt",".pdf",i),sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  
}

#####################################

##一致聚类
library(ConsensusClusterPlus)
path="F:/tp53/3generegulate/"
filename<-list.files(paste(path,"relationship/",sep=""))
for(i in filename[12:13]){
  f= read.csv(paste(path,"relationship/",i,sep=""),header=T,sep="\t")
  f=unique(f)
  m=as.matrix(f[,5:ncol(f)])
  n=scale(m[,1:(grep(".*_no$",colnames(m))[1]-1)]) ##选取突变样本
  #path1=paste(path,"ConsensusClus/",gsub(".txt","",i),"/",sep="")
  setwd("F:/tp53/3generegulate/ConsensusClus/")
  results = ConsensusClusterPlus(n,maxK=10,reps=100,pItem=0.8,pFeature=1,title=gsub(".txt","",i),clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
}


setwd("F:/tp53/3generegulate/ConsensusClus/")
results = ConsensusClusterPlus(n,maxK=10,reps=100,pItem=0.8,pFeature=1,title="1",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png",writeTable = TRUE)

f= read.csv(paste(path,"relationship/",i,sep=""),header=T,sep="\t")
f=unique(f)
m=as.matrix(f[,5:ncol(f)])
n=scale(m[,1:(grep(".*_no$",colnames(m))[1]-1)]) 
f1= read.csv(paste(path,"ConsensusClus/1/1.k=2.consensusClass.csv",sep=""),header=F,sep=",")
n1=n[,as.character(f1[,1])]
row_anno = data.frame(GeneFunction=as.character(f[,1]), row.names=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep="."))
rownames(n1)=paste(substring(as.character(f[,1]),0,4),as.character(f[,2]),as.character(f[,3]),sep=".")
library(pheatmap)
colorsChoice<- colorRampPalette(c("white","red","green"))
pheatmap(n1,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(0,0.5,2.5,4.5),breaks=c(0,0.5,2.5,4.5),fontsize_row=4,fontsize_col=4,annotation_row = row_anno,
         file=paste(path,"ConsensusClus/1/1.pdf",sep=""),main=gsub(".txt","",i),border_color=NA,width=8,height=8)  



for(clu in c("hc","pam","km","kmdist")){
  for(dis in c("pearson","spearman","euclidean","binary","maximum","minkowski","manhattan")){
    results <- ConsensusClusterPlus(m,maxK = 10,reps = 100,pItem = 0.8,pFeature = 1,title = paste0(clu,"-",dis),
                                    clusterAlg = clu,distance = dis,seed = 1262118388.71279,
                                    plot = "pdf",writeTable = TRUE)
    icl <- calcICL(results,title = paste0(clu,"-",dis),plot = "pdf",writeTable = TRUE)
  }
}


##层次聚类
#n=cale(m[,1:(grep(".*_no$",colnames(m))[1]-1)]) 
n=scale(m[,1:201]) #对数据做中心化或者标准化处理
d=dist(t(n))             #计算距离默认欧式距离euclidean
# 聚类方法："ward.D", "ward.D2", "single", "complete", 
#"average" (= UPGMA), "mcquitty" (= WPGMA), 
#"median" (= WPGMC) or "centroid" (= UPGMC)重心法.
hc1<-hclust(d,"ward.D2")  
pdf("F:/tp53/3generegulate/hclust.pdf")
opar<-par(mfrow=c(3,3), mar=c(5.2,4,1,0))  #生成3行3列  图像距离边界的距离
plot(hc1,cex = 0.2,hang=-1)   #hang是表明谱系图中各类所在的位置 当hang取负值时，谱系图中的类从底部画起  生成谱系图
re3<-rect.hclust(hc1,k=2,border="red")  #将分类结果分成2类 用红色矩形笔迹标记
plot(hc1,cex = 0.2,hang=-1)
re4<-rect.hclust(hc1,k=3,border="red")
plot(hc1,cex = 0.2,hang=-1)
re5<-rect.hclust(hc1,k=4,border="red")
plot(hc1,cex = 0.2,hang=-1)
re6<-rect.hclust(hc1,k=5,border="red")
plot(hc1,cex = 0.2,hang=-1)
re7<-rect.hclust(hc1,k=6,border="red")
plot(hc1,cex = 0.2,hang=-1)
re8<-rect.hclust(hc1,k=7,border="red")
plot(hc1,cex = 0.2,hang=-1)
re9<-rect.hclust(hc1,k=8,border="red")
plot(hc1,cex = 0.2,hang=-1)
re10<-rect.hclust(hc1,k=9,border="red")
plot(hc1,cex = 0.2,hang=-1)
re11<-rect.hclust(hc1,k=10,border="red")
#在R软件中 与确定类的个数有关的函数是rect.hclust()函数 它的本质是由给定的个数或给定的阈值来确定聚类的情况
#tree是由hclust生成的结构 k是类的个数 border是矩形框的颜色
par(opar)  #在活动设备中返回所有图形参数和他们的值
dev.off()

#install.packages('ape')
library(ape)
# vector of colors
mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
# cutting dendrogram in 5 clusters
clus5 = cutree(hc1, 5)
# plot
op = par(bg = "#E8DDCB")
# Size reflects miles per gallon
plot(as.phylo(hc1), type = "fan", tip.color = mypal[clus5], label.offset = 1, 
     cex = 0.4, col = "red")
plot(as.phylo(hc1), tip.color = mypal[clus5],cex = 0.4, label.offset = 1)

# install.packages('sparcl')
library(sparcl)
# colors the leaves of a dendrogram
y = cutree(hc1, 3)
ColorDendrogram(hc1, y = y, labels = names(y), main = "My Simulated Data", 
                branchlength = 10)



####九、对Allcancer_fc_updown.txt文件处理
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.csv(paste(path,"Allcancer_fc_updown.txt",sep=""),header=T,sep="\t")
write.table(paste("cancer","Mu_NoPvalue",sep="\t"),file =paste(path,"Allcancer_fisher.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
for(i in 1:nrow(f)){
  x<-matrix(c(f[i,2],f[i,3],f[i,9],f[i,10]),nrow=2,byrow = T)
  d=fisher.test(x)$p.value
  write.table(paste(f[i,1],d,sep="\t"),file =paste(path,"Allcancer_fisher.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}

###对每个cancer，计算tp53突变总数、上、下调、正常样本数
###以及对应的在总的、上、下调、正常中单突变、双突变的样本数目
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
ff=unique(f[,1:5])
m=paste("TP53Cancer",'total','up','down','normal',sep="\t")
cat(m,file=paste(path,"Allcancer_fc_updown_double.txt",sep=""),sep="\n",append=TRUE)
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",as.character(colnames(f1)))]
  if(length(colnames(a[1,which(a[1,]>=1)]))!=0){
    up=colnames(a[1,which(a[1,]>=1)])
  }
  if(length(colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))!=0){
    normal=colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])
  }
  if(length(colnames(a[1,which(a[1,]<= -1)]))!=0){
    down=colnames(a[1,which(a[1,]<= -1)])
  }
  m1=paste(gsub("_fc.txt","_Mu",i),ncol(a)-1,length(up)-1,length(down),length(normal),sep="\t")
  write.table(m1,file=paste(path,"Allcancer_fc_updown_double.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  
  singletotal=nrow(ff[(ff[,1] %in% colnames(a)[-1]) & (ff[,3]==1),])
  doubletotal=nrow(ff[(ff[,1] %in% colnames(a)[-1]) & (ff[,3]>1),])
  singleup=nrow(ff[(ff[,1] %in% up[-1]) & (ff[,3]==1),])
  doubleup=nrow(ff[(ff[,1] %in% up[-1]) & (ff[,3]>1),])
  singledown=nrow(ff[(ff[,1] %in% down) & (ff[,3]==1),])
  doubledown=nrow(ff[(ff[,1] %in% down) & (ff[,3]>1),])
  singlenormal=nrow(ff[(ff[,1] %in% normal) & (ff[,3]==1),])
  doublenormal=nrow(ff[(ff[,1] %in% normal) & (ff[,3]>1),])
  
  m2=paste(gsub("_fc.txt","_singleMu",i),singletotal,singleup,singledown,singlenormal,sep="\t")
  write.table(m2,file=paste(path,"Allcancer_fc_updown_double.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  m3=paste(gsub("_fc.txt","_doubleMu",i),doubletotal,doubleup,doubledown,doublenormal,sep="\t")
  write.table(m3,file=paste(path,"Allcancer_fc_updown_double.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
}

##计算22Allcancer_fc_updown_double.txt的Mu_singleMu、Mu_doubleMu的fisher.test(x2)$p.value
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Allcancer_fc_updown_double.txt",sep=""),header=T,sep="\t")
m=paste("TP53Cancer",'total_up','total_down','total_normal',sep="\t")
cat(m,file=paste(path,"Allcancer_fc_updown_double1.txt",sep=""),sep="\n",append=TRUE)
for(i in seq(1,nrow(f),3)){
  x1<-matrix(c(f[i,2],f[i,3],f[i+1,2],f[i+1,3]),nrow=2,byrow = T)
  x2<-matrix(c(f[i,2],f[i,4],f[i+1,2],f[i+1,4]),nrow=2,byrow = T)
  x3<-matrix(c(f[i,2],f[i,5],f[i+1,2],f[i+1,5]),nrow=2,byrow = T)
  Mu_singleMu=paste(fisher.test(x1)$p.value,fisher.test(x2)$p.value,fisher.test(x3)$p.value,sep="\t")
  write.table(paste(gsub("_Mu","_Mu_singleMu",f[i,1]),Mu_singleMu,sep="\t"),file=paste(path,"Allcancer_fc_updown_double1.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  
  y1<-matrix(c(f[i,2],f[i,3],f[i+2,2],f[i+2,3]),nrow=2,byrow = T)
  y2<-matrix(c(f[i,2],f[i,4],f[i+2,2],f[i+2,4]),nrow=2,byrow = T)
  y3<-matrix(c(f[i,2],f[i,5],f[i+2,2],f[i+2,5]),nrow=2,byrow = T)
  Mu_doubleMu=paste(fisher.test(y1)$p.value,fisher.test(y2)$p.value,fisher.test(y3)$p.value,sep="\t")
  write.table(paste(gsub("_Mu","_Mu_doubleMu",f[i,1]),Mu_doubleMu,sep="\t"),file=paste(path,"Allcancer_fc_updown_double1.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
}


##计算22Allcancer_fc_updown_double.txt的singleMu_doubleMu的fisher.test(x2)$p.value
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Allcancer_fc_updown_double.txt",sep=""),header=T,sep="\t")
m=paste("TP53Cancer",'total_up','total_down','total_normal',sep="\t")
cat(m,file=paste(path,"Allcancer_fc_updown_double2.txt",sep=""),sep="\n",append=TRUE)
for(i in seq(2,nrow(f),3)){
  x1<-matrix(c(f[i,2],f[i,3],f[i+1,2],f[i+1,3]),nrow=2,byrow = T)
  x2<-matrix(c(f[i,2],f[i,4],f[i+1,2],f[i+1,4]),nrow=2,byrow = T)
  x3<-matrix(c(f[i,2],f[i,5],f[i+1,2],f[i+1,5]),nrow=2,byrow = T)
  Mu_singleMu=paste(fisher.test(x1)$p.value,fisher.test(x2)$p.value,fisher.test(x3)$p.value,sep="\t")
  write.table(paste(gsub("_singleMu","_singleMu_doubleMu",f[i,1]),Mu_singleMu,sep="\t"),file=paste(path,"Allcancer_fc_updown_double2.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  
}

##对BRCA、LUAD、BLCA、LIHC、SKCM、PRAD6cancer计算total_up_down,diff_up_down的p值。
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Allcancer_fc_updown.txt",sep=""),header=T,sep="\t")
m=paste("cancer", "Mu_NoMu_total_up", "Mu_NoMu_total_down", "Mu_NoMu_diff_up"," Mu_NoMu_diff_down",sep="\t")
cat(m,file=paste(path,"Allcancer_fc_updown1.txt",sep=""),sep="\n",append=TRUE)
for(i in c(3,4,7,10,14,15)){
  x1<-matrix(c(f[i,2],f[i,5],f[i,9],f[i,12]),nrow=2,byrow = T)
  x2<-matrix(c(f[i,2],f[i,7],f[i,9],f[i,14]),nrow=2,byrow = T)
  x3<-matrix(c(f[i,3],f[i,5],f[i,10],f[i,12]),nrow=2,byrow = T)
  x4<-matrix(c(f[i,3],f[i,7],f[i,10],f[i,14]),nrow=2,byrow = T)
  p=paste(fisher.test(x1)$p.value,fisher.test(x2)$p.value,fisher.test(x3)$p.value,fisher.test(x4)$p.value,sep="\t")
  write.table(paste(f[i,1],p,sep="\t"),file=paste(path,"Allcancer_fc_updown1.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
}



####十、对每个cancer22，allcancer,根据TP53表达值，把样本分为上调下调和不变样本，
#根据每个样本对应的position，找出在对应position上样本改变的数目，画barplot。
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0){
    up=gsub("_no","",colnames(a[1,which(a[1,]>=1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0){
    normal=gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0){
    down=gsub("_no","",colnames(a[1,which(a[1,]<= -1)])) ##GBM的down为空
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
  }
  
  da[,2]=as.numeric(as.character(da[,2]))
  Expression_Value=factor(da[,1])     
  p<-ggplot(da, aes(x = da[,2], y = da[,3], fill = Expression_Value)) +geom_bar(width=0.5,position = "dodge", stat = "identity")+scale_fill_manual(values=c("red","grey","green"))
  p <- p+ facet_grid(da[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da[,2],y=-0.1*max(da[,3]), colour="black",size=2,label=da[,2],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/1/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
dd[,2]=as.numeric(as.character(dd[,2]))
Expression_Value=factor(dd[,1])     
p<-ggplot(dd, aes(x = dd[,2], y = dd[,3], fill = Expression_Value)) +geom_bar(width=0.5,position = "dodge", stat = "identity")+scale_fill_manual(values=c("red","grey","green","black"))
p <- p+ facet_grid(dd[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd[,2],y=-0.1*max(dd[,3]), colour="black",size=2,label=dd[,2],angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/1/Allcancer.pdf",sep=""),width = 10,height=9)

 


###滑动窗口，求和，求均值
slidesum<- function(df, num) {
  if (missing(df) & missing(num))
    return(NULL)
  data=data.frame()
  for(j in as.character(unique(df[,1]))){
    df1=subset(df,as.character(df[,1])==j)
    for(m in 1:(nrow(df1)-num)){
      data=rbind(data,cbind(as.character(df1[m,1]),as.numeric(as.character(df1[m,2])),sum(df1[m:(m+num),3])))
    }
  }
  return(data)
}


slidemean<- function(df, num) {
  if (missing(df) & missing(num))
    return(NULL)
  data=data.frame()
  for(j in as.character(unique(df[,1]))){
    df1=subset(df,as.character(df[,1])==j)
    for(m in 1:(nrow(df1)-num)){
      data=rbind(data,cbind(as.character(df1[m,1]),as.numeric(as.character(df1[m,2])),mean(df1[m:(m+num),3])))
    }
  }
  return(data)
}

##slidesum  用的foldchange值 22cancer-6:CESC\KIRC\KIRP\PCPG\SARC\THYM\"GBM_fc.txt","PRAD_fc.txt","READ_fc.txt",
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
for(i in c("BLCA_fc.txt","BRCA_fc.txt","COAD_fc.txt","ESCA_fc.txt",
           "HNSC_fc.txt","LIHC_fc.txt",
           "LUAD_fc.txt","LUSC_fc.txt","PAAD_fc.txt",
           "SKCM_fc.txt","STAD_fc.txt","THCA_fc.txt",
           "UCEC_fc.txt")){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0){
    up=gsub("_no","",colnames(a[1,which(a[1,]>=1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0){
    normal=gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0){
    down=gsub("_no","",colnames(a[1,which(a[1,]<= -1)])) ##GBM的down为空
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  b4={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b1)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b2)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b3)
  }
  da=rbind(da,cbind(class="all",as.data.frame(table(b4[,2]))))
  da1=slidesum(da,3)
  da1[,2]=as.numeric(as.character(da1[,2]))
  da1[,3]=as.numeric(as.character(da1[,3]))
  Expression_Value=factor(da1[,1])  
  p<-ggplot(da1, aes(x = da1[,2], y = da1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
  p <- p+ facet_grid(da1[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da1[,2],y=0, colour="black",size=2,label=da1[,2],angle=-90)
  p <- p+geom_rect(xmin=1,xmax=40,ymin=-3,ymax=-0.5,fill="red",alpha=0.2)+
    geom_rect(xmin=40,xmax=60,ymin=-3,ymax=-0.5,fill="violet",alpha=0.5)+
    geom_rect(xmin=60,xmax=95,ymin=-3,ymax=-0.5,fill="yellow",alpha=0.6)+
    geom_rect(xmin=95,xmax=100,ymin=-3,ymax=-0.5,fill="grey",alpha=0.8)+
    geom_rect(xmin=100,xmax=300,ymin=-3,ymax=-0.5,fill="blue",alpha=0.5)+
    geom_rect(xmin=300,xmax=325,ymin=-3,ymax=-0.5,fill="grey",alpha=0.8)+
    geom_rect(xmin=325,xmax=356,ymin=-3,ymax=-0.5,fill="orange",alpha=0.5)+
    geom_rect(xmin=356,xmax=393,ymin=-3,ymax=-0.5,fill="green",alpha=0.5)
  ggsave(file=paste(path,"Position_SampleNumber/slidesum3f/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
dd1=slidesum(dd,3)
dd1[,2]=as.numeric(as.character(dd1[,2]))
dd1[,3]=as.numeric(as.character(dd1[,3]))
Expression_Value=factor(dd1[,1])  
p<-ggplot(dd1, aes(x = dd1[,2], y = dd1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
p <- p+ facet_grid(dd1[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p+ xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd1[,2],y=-0.1, colour="black",size=2,label=dd1[,2],angle=-90)
p <- p+geom_rect(xmin=1,xmax=40,ymin=-10,ymax=-7,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-10,ymax=-7,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-10,ymax=-7,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-10,ymax=-7,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-10,ymax=-7,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-10,ymax=-7,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-10,ymax=-7,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-10,ymax=-7,fill="green",alpha=0.5)
ggsave(file=paste(path,"Position_SampleNumber/slidesum3f/Allcancer.pdf",sep=""),width = 10,height=9)

 
##slidesum  用的foldchange值 22cancer-6:CESC\KIRC\KIRP\PCPG\SARC\THYM\"GBM_fc.txt","PRAD_fc.txt","READ_fc.txt",
##扩充0~400，补0
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
df=data.frame(cbind(class=rep(c("up","normal","down","all"),each=400),Var1=rep(1:400,4),Freq=rep(0,1600)))
df[,1]=factor(df[,1],c("up","normal","down","all"))
df[,3]=as.numeric(as.character((df[,3])))
for(i in c("BLCA_fc.txt","BRCA_fc.txt","COAD_fc.txt","ESCA_fc.txt",
           "HNSC_fc.txt","LIHC_fc.txt",
           "LUAD_fc.txt","LUSC_fc.txt","PAAD_fc.txt",
           "SKCM_fc.txt","STAD_fc.txt","THCA_fc.txt",
           "UCEC_fc.txt")){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0){
    up=gsub("_no","",colnames(a[1,which(a[1,]>=1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0){
    normal=gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0){
    down=gsub("_no","",colnames(a[1,which(a[1,]<= -1)])) ##GBM的down为空
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  b4={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b1)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b2)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b3)
  }
  da=rbind(da,cbind(class="all",as.data.frame(table(b4[,2]))))
  daa=merge(df,da,by=c("class","Var1"),all.x=T)[,-3] 
  daa[,2]=as.numeric(as.character((daa[,2])))
  daa[sapply(daa,is.na)]<-0
  daa=daa[order(daa[,1],daa[,2]),]
  da1=slidesum(daa,3)
  da1[,2]=as.numeric(as.character(da1[,2]))
  da1[,3]=as.numeric(as.character(da1[,3]))
  Expression_Value=factor(da1[,1])  
  p<-ggplot(da1, aes(x = da1[,2], y = da1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
  p <- p+ facet_grid(da1[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da1[,2],y=-0.5, colour="black",size=2,label=da1[,2],angle=-90)
  p <- p+geom_rect(xmin=1,xmax=40,ymin=-5,ymax=-1,fill="red",alpha=0.2)+
    geom_rect(xmin=40,xmax=60,ymin=-5,ymax=-1,fill="violet",alpha=0.5)+
    geom_rect(xmin=60,xmax=95,ymin=-5,ymax=-1,fill="yellow",alpha=0.6)+
    geom_rect(xmin=95,xmax=100,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
    geom_rect(xmin=100,xmax=300,ymin=-5,ymax=-1,fill="blue",alpha=0.5)+
    geom_rect(xmin=300,xmax=325,ymin=-5,ymax=-1,fill="grey",alpha=0.8)+
    geom_rect(xmin=325,xmax=356,ymin=-5,ymax=-1,fill="orange",alpha=0.5)+
    geom_rect(xmin=356,xmax=393,ymin=-5,ymax=-1,fill="green",alpha=0.5)
  ggsave(file=paste(path,"Position_SampleNumber/slidesum3f0/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
ddd=merge(df,dd,by=c("class","Var1"),all.x=T)[,-3] 
ddd[,2]=as.numeric(as.character((ddd[,2])))
ddd[sapply(ddd,is.na)]<-0
ddd=ddd[order(ddd[,1],ddd[,2]),]
dd1=slidesum(ddd,3)
dd1[,2]=as.numeric(as.character(dd1[,2]))
dd1[,3]=as.numeric(as.character(dd1[,3]))
Expression_Value=factor(dd1[,1])  
p<-ggplot(dd1, aes(x = dd1[,2], y = dd1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
p <- p+ facet_grid(dd1[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p+ xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd1[,2],y=-3, colour="black",size=2,label=dd1[,2],angle=-90)
p <- p+geom_rect(xmin=1,xmax=40,ymin=-10,ymax=-7,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-10,ymax=-7,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-10,ymax=-7,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-10,ymax=-7,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-10,ymax=-7,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-10,ymax=-7,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-10,ymax=-7,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-10,ymax=-7,fill="green",alpha=0.5)
ggsave(file=paste(path,"Position_SampleNumber/slidesum3f0/Allcancer.pdf",sep=""),width = 10,height=9)








##slidesum  用的_foldchange_p53%WTmean值 
##27cancer-12:CESC\"COAD_fc.txt","KIRC_fc.txt","GBM_fc.txt","PAAD_fc.txt","STAD_fc.txt","LAML_fc.txt",
##"UCEC_fc.txt","PCPG_fc.txt","TGCT_fc.txt","THYM_fc.txt","UCS_fc.txt"
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
for(i in c("BLCA_fc.txt","BRCA_fc.txt","ESCA_fc.txt",
           "HNSC_fc.txt","LIHC_fc.txt","PRAD_fc.txt","READ_fc.txt",
           "LUAD_fc.txt","LUSC_fc.txt","KIRP_fc.txt","OV_fc.txt",
           "SKCM_fc.txt","THCA_fc.txt","LGG_fc.txt",
           "SARC_fc.txt")){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  colnames(f1)<-gsub("_p53/WT_mean","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(colnames(a[1,which(a[1,]>=1)]))!=0){
    up=colnames(a[1,which(a[1,]>=1)])
  }
  if(length(colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))!=0){
    normal=colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])
  }
  if(length(colnames(a[1,which(a[1,]<= -1)]))!=0){
    down=colnames(a[1,which(a[1,]<= -1)])
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(colnames(a[1,which(a[1,]>=1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(colnames(a[1,which(a[1,]<= -1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  b4={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b1)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b2)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b3)
  }
  da=rbind(da,cbind(class="all",as.data.frame(table(b4[,2]))))
  da1=slidesum(da,3)
  da1[,2]=as.numeric(as.character(da1[,2]))
  da1[,3]=as.numeric(as.character(da1[,3]))
  Expression_Value=factor(da1[,1])  
  p<-ggplot(da1, aes(x = da1[,2], y = da1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
  p <- p+ facet_grid(da1[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da1[,2],y=0, colour="black",size=2,label=da1[,2],angle=-90)
  p <- p+geom_rect(xmin=1,xmax=40,ymin=-3,ymax=-0.5,fill="red",alpha=0.2)+
    geom_rect(xmin=40,xmax=60,ymin=-3,ymax=-0.5,fill="violet",alpha=0.5)+
    geom_rect(xmin=60,xmax=95,ymin=-3,ymax=-0.5,fill="yellow",alpha=0.6)+
    geom_rect(xmin=95,xmax=100,ymin=-3,ymax=-0.5,fill="grey",alpha=0.8)+
    geom_rect(xmin=100,xmax=300,ymin=-3,ymax=-0.5,fill="blue",alpha=0.5)+
    geom_rect(xmin=300,xmax=325,ymin=-3,ymax=-0.5,fill="grey",alpha=0.8)+
    geom_rect(xmin=325,xmax=356,ymin=-3,ymax=-0.5,fill="orange",alpha=0.5)+
    geom_rect(xmin=356,xmax=393,ymin=-3,ymax=-0.5,fill="green",alpha=0.5)
  ggsave(file=paste(path,"Position_SampleNumber/slidesum3p53_WTmean/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
dd1=slidesum(dd,3)
dd1[,2]=as.numeric(as.character(dd1[,2]))
dd1[,3]=as.numeric(as.character(dd1[,3]))
Expression_Value=factor(dd1[,1])  
p<-ggplot(dd1, aes(x = dd1[,2], y = dd1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
p <- p+ facet_grid(dd1[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p+ xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd1[,2],y=-0.1, colour="black",size=2,label=dd1[,2],angle=-90)
p <- p+geom_rect(xmin=1,xmax=40,ymin=-25,ymax=-5,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-25,ymax=-5,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-25,ymax=-5,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-25,ymax=-5,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-25,ymax=-5,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-25,ymax=-5,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-25,ymax=-5,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-25,ymax=-5,fill="green",alpha=0.5)
ggsave(file=paste(path,"Position_SampleNumber/slidesum3p53_WTmean/Allcancer.pdf",sep=""),width = 10,height=9)


 

##slidesum  用的_foldchange_p53%WTmean值、扩充补0 
##27cancer-12:CESC\"COAD_fc.txt","KIRC_fc.txt","GBM_fc.txt","PAAD_fc.txt","STAD_fc.txt","LAML_fc.txt",
##"UCEC_fc.txt","PCPG_fc.txt","TGCT_fc.txt","THYM_fc.txt","UCS_fc.txt"
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
df=data.frame(cbind(class=rep(c("up","normal","down","all"),each=400),Var1=rep(1:400,4),Freq=rep(0,1600)))
df[,1]=factor(df[,1],c("up","normal","down","all"))
df[,3]=as.numeric(as.character((df[,3])))
for(i in c("BLCA_fc.txt","BRCA_fc.txt","ESCA_fc.txt",
           "HNSC_fc.txt","LIHC_fc.txt","PRAD_fc.txt","READ_fc.txt",
           "LUAD_fc.txt","LUSC_fc.txt","KIRP_fc.txt","OV_fc.txt",
           "SKCM_fc.txt","THCA_fc.txt","LGG_fc.txt",
           "SARC_fc.txt")){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  colnames(f1)<-gsub("_p53/WT_mean","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(colnames(a[1,which(a[1,]>=1)]))!=0){
    up=colnames(a[1,which(a[1,]>=1)])
  }
  if(length(colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))!=0){
    normal=colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])
  }
  if(length(colnames(a[1,which(a[1,]<= -1)]))!=0){
    down=colnames(a[1,which(a[1,]<= -1)])
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(colnames(a[1,which(a[1,]>=1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(colnames(a[1,which(a[1,]<= -1)]))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  b4={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b1)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b2)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
    b4=rbind(b4,b3)
  }
  da=rbind(da,cbind(class="all",as.data.frame(table(b4[,2]))))
  daa=merge(df,da,by=c("class","Var1"),all.x=T)[,-3] 
  daa[,2]=as.numeric(as.character((daa[,2])))
  daa[sapply(daa,is.na)]<-0
  daa=daa[order(daa[,1],daa[,2]),]
  da1=slidesum(daa,3)
  da1[,2]=as.numeric(as.character(da1[,2]))
  da1[,3]=as.numeric(as.character(da1[,3]))
  Expression_Value=factor(da1[,1])  
  p<-ggplot(da1, aes(x = da1[,2], y = da1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
  p <- p+ facet_grid(da1[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da1[,2],y=-0.5, colour="black",size=2,label=da1[,2],angle=-90)
  p <- p+geom_rect(xmin=1,xmax=40,ymin=-3,ymax=-1,fill="red",alpha=0.2)+
    geom_rect(xmin=40,xmax=60,ymin=-3,ymax=-1,fill="violet",alpha=0.5)+
    geom_rect(xmin=60,xmax=95,ymin=-3,ymax=-1,fill="yellow",alpha=0.6)+
    geom_rect(xmin=95,xmax=100,ymin=-3,ymax=-1,fill="grey",alpha=0.8)+
    geom_rect(xmin=100,xmax=300,ymin=-3,ymax=-1,fill="blue",alpha=0.5)+
    geom_rect(xmin=300,xmax=325,ymin=-3,ymax=-1,fill="grey",alpha=0.8)+
    geom_rect(xmin=325,xmax=356,ymin=-3,ymax=-1,fill="orange",alpha=0.5)+
    geom_rect(xmin=356,xmax=393,ymin=-3,ymax=-1,fill="green",alpha=0.5)
  ggsave(file=paste(path,"Position_SampleNumber/slidesum3p53_WTmean0/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
ddd=merge(df,dd,by=c("class","Var1"),all.x=T)[,-3] 
ddd[,2]=as.numeric(as.character((ddd[,2])))
ddd[sapply(ddd,is.na)]<-0
ddd=ddd[order(ddd[,1],ddd[,2]),]
dd1=slidesum(ddd,3)
dd1[,2]=as.numeric(as.character(dd1[,2]))
dd1[,3]=as.numeric(as.character(dd1[,3]))
Expression_Value=factor(dd1[,1])  
p<-ggplot(dd1, aes(x = dd1[,2], y = dd1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="grey") + geom_line()  
p <- p+ facet_grid(dd1[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p+ xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd1[,2],y=-1, colour="black",size=2,label=dd1[,2],angle=-90)
p <- p+geom_rect(xmin=1,xmax=40,ymin=-20,ymax=-5,fill="red",alpha=0.2)+
  geom_rect(xmin=40,xmax=60,ymin=-20,ymax=-5,fill="violet",alpha=0.5)+
  geom_rect(xmin=60,xmax=95,ymin=-20,ymax=-5,fill="yellow",alpha=0.6)+
  geom_rect(xmin=95,xmax=100,ymin=-20,ymax=-5,fill="grey",alpha=0.8)+
  geom_rect(xmin=100,xmax=300,ymin=-20,ymax=-5,fill="blue",alpha=0.5)+
  geom_rect(xmin=300,xmax=325,ymin=-20,ymax=-5,fill="grey",alpha=0.8)+
  geom_rect(xmin=325,xmax=356,ymin=-20,ymax=-5,fill="orange",alpha=0.5)+
  geom_rect(xmin=356,xmax=393,ymin=-20,ymax=-5,fill="green",alpha=0.5)
ggsave(file=paste(path,"Position_SampleNumber/slidesum3p53_WTmean0/Allcancer1.pdf",sep=""),width = 10,height=9)





##slidemean 22cancer-6:CESC\KIRC\KIRP\PCPG\SARC\THYM\
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
d1={}
d2={}
d3={}
for(i in c("BLCA_fc.txt","BRCA_fc.txt","COAD_fc.txt","ESCA_fc.txt",
           "GBM_fc.txt","HNSC_fc.txt","LIHC_fc.txt",
           "LUAD_fc.txt","LUSC_fc.txt","PAAD_fc.txt","PRAD_fc.txt",
           "READ_fc.txt","SKCM_fc.txt","STAD_fc.txt","THCA_fc.txt",
           "UCEC_fc.txt")){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0){
    up=gsub("_no","",colnames(a[1,which(a[1,]>=1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0){
    normal=gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0){
    down=gsub("_no","",colnames(a[1,which(a[1,]<= -1)])) ##GBM的down为空
  }
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  if(length(gsub("_no","",colnames(a[1,which(a[1,]>=1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% up)){
    b1=na.omit(b[as.character(b[,1]) %in% up,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% normal)){
    b2=na.omit(b[as.character(b[,1]) %in% normal,c("submitter_id","HGVSp_Short")])
  }
  if(length(gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))!=0 & TRUE %in% (as.character(b[,1]) %in% down)){
    b3=na.omit(b[as.character(b[,1]) %in% down,c("submitter_id","HGVSp_Short")])
  }
  
  if(nrow(b1)!=0){
    d1=rbind(d1,b1)
  }
  if(nrow(b2)!=0){
    d2=rbind(d2,b2)
  }
  if(nrow(b3)!=0){
    d3=rbind(d3,b3)
  }
  
  da={}
  if(nrow(b1)!=0){
    d=cbind(class="up",as.data.frame(table(b1[,2])))
    da=rbind(da,d)
  }
  if(nrow(b2)!=0){
    d=cbind(class="normal",as.data.frame(table(b2[,2])))
    da=rbind(da,d)
  }
  if(nrow(b3)!=0){
    d=cbind(class="down",as.data.frame(table(b3[,2])))
    da=rbind(da,d)
  }
  da1=slidemean(da,3)
  da1[,2]=as.numeric(as.character(da1[,2]))
  da1[,3]=as.numeric(as.character(da1[,3]))
  Expression_Value=factor(da1[,1])  
  p<-ggplot(da1, aes(x = da1[,2], y = da1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
  p <- p+ facet_grid(da1[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle(gsub("_fc.txt","",i))+xlim(0,400) 
  p <- p+annotate("text",x=da1[,2],y=-0.1, colour="black",size=2,label=da1[,2],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/slidemean/",gsub("_fc.txt","",i),".pdf",sep=""),width = 10,height=9)
}
d4=rbind(d1,d2,d3)
dd=rbind(cbind(class="up",as.data.frame(table(d1[,2]))),cbind(class="normal",as.data.frame(table(d2[,2]))),cbind(class="down",as.data.frame(table(d3[,2]))),cbind(class="all",as.data.frame(table(d4[,2]))))
dd1=slidemean(dd,3)
dd1[,2]=as.numeric(as.character(dd1[,2]))
dd1[,3]=as.numeric(as.character(dd1[,3]))
Expression_Value=factor(dd1[,1])  
p<-ggplot(dd1, aes(x = dd1[,2], y = dd1[,3], fill = Expression_Value))+geom_point(aes(shape=Expression_Value),colour="red") + geom_line()  
p <- p+ facet_grid(dd1[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Allcancer")+xlim(0,400) 
p <- p+annotate("text",x=dd1[,2],y=-0.1, colour="black",size=2,label=dd1[,2],angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/slidemean/Allcancer.pdf",sep=""),width = 10,height=9)






####十一、对每个cancer27，看tp53突变样本，以及在各个domain(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)的样本的生存分析
###PCPG17\PRAD18\TGCT23\THCA24\UCS27只有一组数据删去 1:16 19:22 25:26
###27-5=22cancer Mu+6domain=7class
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
ff<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
df=unique(na.omit(f[,c("submitter_id","project_id","HGVSp_Short")]))
df[,2]=as.character(df[,2])
filename<-list.files(paste(path,"p53_context_mutation_fpkm/",sep=""))
h={}
for(i in filename[c(1:16,19:22,25:26)]){
  f1<-read.csv(paste(path,"p53_context_mutation_fpkm/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  Mu=colnames(f1)[-1]
  TAD1=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]<=40),1]
  TAD2=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>40 & df[,3]<=60),1]
  PRD=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>=325 & df[,3]<=356),1]
  Basic=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>356 & df[,3]<=393),1]
  
  pdf(paste(path,"p53Mu_domain_survival/p53Mu_domain/",gsub("_p53_context_mutation_fpkm.txt",".pdf",i),sep=""))
  ff1={}
  if(nrow(ff[Mu,])!=0){
    d=cbind(class=1,ff[Mu,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD1,])!=0){
    d=cbind(class=2,ff[TAD1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD2,])!=0){
    d=cbind(class=3,ff[TAD2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD,])!=0){
    d=cbind(class=4,ff[PRD,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding,])!=0){
    d=cbind(class=5,ff[DNAbinding,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet,])!=0){
    d=cbind(class=6,ff[Tet,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic,])!=0){
    d=cbind(class=7,ff[Basic,])
    ff1=rbind(ff1,d)
  }
  ff1=na.omit(ff1)
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h=rbind(h,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  co=c("black","red","orange","blue","purple","yellow","green")
  plot(kmsurvival1, lty = 'solid', col=co[as.numeric(unique(e1[,1]))],
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_p53_context_mutation_fpkm.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c("tp53Mu","TAD1","TAD2","PRD","DNAbinding","Tet","Basic"), lty='solid',
         col=c("black","red","orange","blue","purple","yellow","green"))
  text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
  for(j in as.numeric(unique(e1[,1]))){
    text(500,0.9-(j-1)*0.03,cex=0.8,j)
  }
  for(j in 1:length(as.numeric(unique(e1[,1]))) ){
    text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
  }
  dev.off()
}
pdf(paste(path,"p53Mu_domain_survival/p53Mu_domain/Allcancer.pdf",sep=""))  
dif2 <- survdiff(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1],conf.type = "log-log")
co2=c("black","red","orange","blue","purple","yellow","green")
plot(kmsurvival2, lty = 'solid', col=co2[as.numeric(unique(h[,1]))],
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c("tp53Mu","TAD1","TAD2","PRD","DNAbinding","Tet","Basic"), lty='solid',
       col=c("black","red","orange","blue","purple","yellow","green"))
text(700,1,cex=0.8,paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),4),sep=""))
for(j in as.numeric(unique(h[,1]))){
  text(500,0.9-(j-1)*0.03,cex=0.8,j)
}
for(j in 1:length(as.numeric(unique(h[,1]))) ){
  text(800,0.9-(j-1)*0.03,cex=0.8,dif2$n[j])
}
dev.off()


####十一、对每个cancer27，看tp53突变样本，以及在各个domain(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)的样本的生存分析
###KIRC8\KIRP9\LAML10\PCPG17\PRAD18\READ19\SKCM21\TGCT23\THCA24\THYM25\UCS27只有一组数据删去 1:7 11:16 20,22 26
###27-11=16cancer 只有6domain=6class
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
ff<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
df=unique(na.omit(f[,c("submitter_id","project_id","HGVSp_Short")]))
df[,2]=as.character(df[,2])
filename<-list.files(paste(path,"p53_context_mutation_fpkm/",sep=""))
h={}
for(i in filename[c(1:7,11:16,20,22,26)]){
  f1<-read.csv(paste(path,"p53_context_mutation_fpkm/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  Mu=colnames(f1)[-1]
  TAD1=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]<=40),1]
  TAD2=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>40 & df[,3]<=60),1]
  PRD=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>=325 & df[,3]<=356),1]
  Basic=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>356 & df[,3]<=393),1]
  
  pdf(paste(path,"p53Mu_domain_survival/domain/",gsub("_p53_context_mutation_fpkm.txt",".pdf",i),sep=""))
  ff1={}
  if(nrow(ff[TAD1,])!=0){
    d=cbind(class=1,ff[TAD1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD2,])!=0){
    d=cbind(class=2,ff[TAD2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD,])!=0){
    d=cbind(class=3,ff[PRD,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding,])!=0){
    d=cbind(class=4,ff[DNAbinding,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet,])!=0){
    d=cbind(class=5,ff[Tet,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic,])!=0){
    d=cbind(class=6,ff[Basic,])
    ff1=rbind(ff1,d)
  }
  ff1=na.omit(ff1)
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h=rbind(h,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  co=c("red","orange","blue","purple","yellow","green")
  plot(kmsurvival1, lty = 'solid', col=co[as.numeric(unique(e1[,1]))],
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_p53_context_mutation_fpkm.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1","TAD2","PRD","DNAbinding","Tet","Basic"), lty='solid',
         col=c("red","orange","blue","purple","yellow","green"))
  text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
  for(j in as.numeric(unique(e1[,1]))){
    text(500,0.9-(j-1)*0.03,cex=0.8,j)
  }
  for(j in 1:length(as.numeric(unique(e1[,1]))) ){
    text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
  }
  dev.off()
}
pdf(paste(path,"p53Mu_domain_survival/domain/Allcancer.pdf",sep=""))  
dif2 <- survdiff(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1],conf.type = "log-log")
co2=c("red","orange","blue","purple","yellow","green")
plot(kmsurvival2, lty = 'solid', col=co2[as.numeric(unique(h[,1]))],
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1","TAD2","PRD","DNAbinding","Tet","Basic"), lty='solid',
       col=c("red","orange","blue","purple","yellow","green"))
text(700,1,cex=0.8,paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),4),sep=""))
for(j in as.numeric(unique(h[,1]))){
  text(500,0.9-(j-1)*0.03,cex=0.8,j)
}
for(j in 1:length(as.numeric(unique(h[,1]))) ){
  text(800,0.9-(j-1)*0.03,cex=0.8,dif2$n[j])
}
dev.off()


####十一、对每个cancer27，看tp53突变样本，以及在各个domain(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)的样本的生存分析
###KIRC8\KIRP9\LAML10\PCPG17\PRAD18\READ19\SKCM21\TGCT23\THCA24\THYM25\UCS27只有一组数据删去 1:7 11:16 20,22 26
###27-11=16cancer 只有6domain, 在domain和no_domain=2class,出6图
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
ff<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
df=unique(na.omit(f[,c("submitter_id","project_id","HGVSp_Short")]))
df[,2]=as.character(df[,2])
filename<-list.files(paste(path,"p53_context_mutation_fpkm/",sep=""))
h={}
for(i in filename[c(1:7,11:16,20,22,26)]){
  f1<-read.csv(paste(path,"p53_context_mutation_fpkm/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  Mu=colnames(f1)[-1]
  TAD1=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]<=40),1]
  TAD2=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>40 & df[,3]<=60),1]
  PRD=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>=325 & df[,3]<=356),1]
  Basic=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>356 & df[,3]<=393),1]
  ff1={}
  if(nrow(ff[TAD1,])!=0){
    d=cbind(class=1,ff[TAD1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD2,])!=0){
    d=cbind(class=2,ff[TAD2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD,])!=0){
    d=cbind(class=3,ff[PRD,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding,])!=0){
    d=cbind(class=4,ff[DNAbinding,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet,])!=0){
    d=cbind(class=5,ff[Tet,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic,])!=0){
    d=cbind(class=6,ff[Basic,])
    ff1=rbind(ff1,d)
  }
  ff1=na.omit(ff1)
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=data.frame(rbind(alive1,dead1))
  e1=na.omit(e1)
  h=rbind(h,e1)
}
pdf(paste(path,"p53Mu_domain_survival/domain/Allcancer1.pdf",sep=""))  
dif2 <- survdiff(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1],conf.type = "log-log")
co2=c("red","orange","blue","purple","yellow4","darkgreen")
plot(kmsurvival2, lty = 'solid', col=co2[as.numeric(unique(h[,1]))],lwd=2,
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,c("TAD1","TAD2","PRD","DNAbinding","Tet","Basic"), lty='solid',lwd=2,
       col=c("red","orange","blue","purple","yellow4","darkgreen"))
text(2700,1,cex=0.8,paste("p=",signif(pchisq(dif2$chisq,1,lower.tail=F),3),sep=""))
for(j in as.numeric(unique(h[,1]))){
  text(2500,0.9-(j-1)*0.03,cex=0.8,j)
}
for(j in 1:length(as.numeric(unique(h[,1]))) ){
  text(2800,0.9-(j-1)*0.03,cex=0.8,dif2$n[j])
}
dev.off()

##Allcancer在domain和no_domain=2class,出6图(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)
h1=data.frame(h) 
h1[h1[,1]==6,4]="Basic"
h1[h1[,1]!=6,4]="No_Basic"
h1[,2]=as.numeric(as.character(h1[,2]))
h1[,3]=as.numeric(as.character(h1[,3]))
pdf(paste(path,"p53Mu_domain_survival/domain/Allcancer_Basic.pdf",sep=""))  
dif1 <- survdiff(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,4])#求生存时间
kmsurvival1<-survfit(Surv(as.numeric(h1[,3]),as.numeric(h1[,2]))~h1[,4],conf.type = "log-log")
#c("red","orange","blue","purple","yellow4","darkgreen")
plot(kmsurvival1, lty = c('dashed','solid'),lwd=2, col="darkgreen",
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,c("Basic","No_Basic"),lwd=2,lty=c('dashed','solid'),col="darkgreen")
text(2700,1,cex=0.8,paste("p=",signif(pchisq(dif1$chisq,1,lower.tail=F),3),sep=""))
for(j in 1:length(as.numeric(unique(h1[,4]))) ){
  text(2800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
}
dev.off()










####十一、对每个cancer27，看tp53突变样本，以及在各个domain(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)的样本的生存分析
###KIRC8\KIRP9\LAML10\PCPG17\PRAD18\READ19\SKCM21\TGCT23\THCA24\THYM25\UCS27只有一组数据删去 1:7 11:16 20,22 26
###27-11=16cancer 只考虑前5个domain=5class
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
ff<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
df=unique(na.omit(f[,c("submitter_id","project_id","HGVSp_Short")]))
df[,2]=as.character(df[,2])
filename<-list.files(paste(path,"p53_context_mutation_fpkm/",sep=""))
h={}
for(i in filename[c(1:7,11:16,20,22,26)]){
  f1<-read.csv(paste(path,"p53_context_mutation_fpkm/",i,sep=""),header=T,sep="\t",stringsAsFactors = F)
  Mu=colnames(f1)[-1]
  TAD1=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]<=40),1]
  TAD2=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>40 & df[,3]<=60),1]
  PRD=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>=325 & df[,3]<=356),1]
  #Basic=df[(df[,2]==gsub("_p53_context_mutation_fpkm.txt","",i)) & (df[,1] %in% Mu) & (df[,3]>356 & df[,3]<=393),1]
  
  pdf(paste(path,"p53Mu_domain_survival/domain/1/",gsub("_p53_context_mutation_fpkm.txt",".pdf",i),sep=""))
  ff1={}
  if(nrow(ff[TAD1,])!=0){
    d=cbind(class=1,ff[TAD1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD2,])!=0){
    d=cbind(class=2,ff[TAD2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD,])!=0){
    d=cbind(class=3,ff[PRD,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding,])!=0){
    d=cbind(class=4,ff[DNAbinding,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet,])!=0){
    d=cbind(class=5,ff[Tet,])
    ff1=rbind(ff1,d)
  }
  ff1=na.omit(ff1)
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h=rbind(h,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  co=c("red","orange","blue","purple","yellow")
  plot(kmsurvival1, lty = 'solid', col=co[as.numeric(unique(e1[,1]))],
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_p53_context_mutation_fpkm.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1","TAD2","PRD","DNAbinding","Tet"), lty='solid',
         col=c("red","orange","blue","purple","yellow"))
  text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
  for(j in as.numeric(unique(e1[,1]))){
    text(500,0.9-(j-1)*0.03,cex=0.8,j)
  }
  for(j in 1:length(as.numeric(unique(e1[,1]))) ){
    text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
  }
  dev.off()
}
pdf(paste(path,"p53Mu_domain_survival/domain/1/Allcancer.pdf",sep=""))  
dif2 <- survdiff(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1],conf.type = "log-log")
co2=c("red","orange","blue","purple","yellow")
plot(kmsurvival2, lty = 'solid', col=co2[as.numeric(unique(h[,1]))],
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1","TAD2","PRD","DNAbinding","Tet"), lty='solid',
       col=c("red","orange","blue","purple","yellow"))
text(700,1,cex=0.8,paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),4),sep=""))
for(j in as.numeric(unique(h[,1]))){
  text(500,0.9-(j-1)*0.03,cex=0.8,j)
}
for(j in 1:length(as.numeric(unique(h[,1]))) ){
  text(800,0.9-(j-1)*0.03,cex=0.8,dif2$n[j])
}
dev.off()



####十一、对每个cancer22\allcancer，把各个domain(6个TAD1、TAD2、PRD、DNAbinding、Tet、Basic)的样本
####分为上调、下调、正常的样本，6x3=18class做生存分析22cancer
##只有1group删去KIRC8\KIRP9\PCPG14\PRAD15\READ16\THCA20\THYM21=22-7=15cancer
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
ff<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
ff[,1]=gsub("-",".",ff[,1])
rownames(ff)=ff[,1]
ff=na.omit(ff)
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f[,1]=gsub("-",".",f[,1])
f$HGVSp_Short=as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4)))
df=unique(na.omit(f[,c("submitter_id","project_id","HGVSp_Short")]))
df[,2]=as.character(df[,2])
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
h={}
for(i in filename[c(1:7,10:13,17:19,22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]
  a1=colnames(a[1,which(a[1,]>=1)])
  a2=colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])
  a3=colnames(a[1,which(a[1,]<= -1)])
  
  TAD11=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]<=40),1]
  TAD21=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]>40 & df[,3]<=60),1]
  PRD1=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding1=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet1=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]>=325 & df[,3]<=356),1]
  Basic1=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a1) & (df[,3]>356 & df[,3]<=393),1]
  
  TAD12=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]<=40),1]
  TAD22=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]>40 & df[,3]<=60),1]
  PRD2=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding2=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet2=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]>=325 & df[,3]<=356),1]
  Basic2=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a2) & (df[,3]>356 & df[,3]<=393),1]
  
  TAD13=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]<=40),1]
  TAD23=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]>40 & df[,3]<=60),1]
  PRD3=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]>60 & df[,3]<=95),1]
  DNAbinding3=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]>= 100 & df[,3]<=300),1]
  Tet3=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]>=325 & df[,3]<=356),1]
  Basic3=df[(df[,2]==gsub("_fc.txt","",i)) & (df[,1] %in% a3) & (df[,3]>356 & df[,3]<=393),1]
  
  pdf(paste(path,"p53Mu_domain_survival/domain_updown/",gsub("_fc.txt",".pdf",i),sep=""))
  ff1={}
  if(nrow(ff[TAD11,])!=0){
    d=cbind(class=1,ff[TAD11,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD12,])!=0){
    d=cbind(class=2,ff[TAD12,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD13,])!=0){
    d=cbind(class=3,ff[TAD13,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD21,])!=0){
    d=cbind(class=4,ff[TAD21,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD22,])!=0){
    d=cbind(class=5,ff[TAD22,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[TAD23,])!=0){
    d=cbind(class=6,ff[TAD23,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD1,])!=0){
    d=cbind(class=7,ff[PRD1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD2,])!=0){
    d=cbind(class=8,ff[PRD2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[PRD3,])!=0){
    d=cbind(class=9,ff[PRD3,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding1,])!=0){
    d=cbind(class=10,ff[DNAbinding1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding2,])!=0){
    d=cbind(class=11,ff[DNAbinding2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[DNAbinding3,])!=0){
    d=cbind(class=12,ff[DNAbinding3,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet1,])!=0){
    d=cbind(class=13,ff[Tet1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet2,])!=0){
    d=cbind(class=14,ff[Tet2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Tet3,])!=0){
    d=cbind(class=15,ff[Tet3,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic1,])!=0){
    d=cbind(class=16,ff[Basic1,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic2,])!=0){
    d=cbind(class=17,ff[Basic2,])
    ff1=rbind(ff1,d)
  }
  if(nrow(ff[Basic3,])!=0){
    d=cbind(class=18,ff[Basic3,])
    ff1=rbind(ff1,d)
  }
  ff1=na.omit(ff1)
  b11=subset(ff1,ff1$vital_status.y==0)
  b21=subset(ff1,ff1$vital_status.y==1)
  alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
  dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
  e1=rbind(alive1,dead1)
  e1=na.omit(e1)
  h=rbind(h,e1)
  dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
  kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
  co=c("red4","red2","red","orange4","orange3","orange","maroon4", "maroon1","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue2","royalblue1")
  plot(kmsurvival1, lty = 'solid', col=co[as.numeric(unique(e1[,1]))],
       xlab='survival time in days',ylab='survival probabilities',main=gsub("_p53_context_mutation_fpkm.txt","",i))
  legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1_up","TAD1_normal","TAD1_down","TAD21","TAD22","TAD23","PRD1","PRD2","PRD3","DNAbinding1","DNAbinding2","DNAbinding3","Tet1","Tet2","Tet3","Basic1","Basic2","Basic3"), lty='solid',
         col=c("red4","red2","red","orange4","orange3","orange","maroon4", "maroon1","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue2","royalblue1"))
  text(700,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
  for(j in as.numeric(unique(e1[,1]))){
    text(500,0.9-(j-1)*0.03,cex=0.8,j)
  }
  for(j in 1:length(as.numeric(unique(e1[,1]))) ){
    text(800,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
  }
  dev.off()
}
pdf(paste(path,"p53Mu_domain_survival/domain_updown/Allcancer.pdf",sep=""))  
dif2 <- survdiff(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1])#求生存时间
kmsurvival2<-survfit(Surv(as.numeric(h[,3]),as.numeric(h[,2]))~h[,1],conf.type = "log-log")
co2=c("red4","red2","red","orange4","orange3","orange","maroon4", "maroon1","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue2","royalblue1")
plot(kmsurvival2, lty = 'solid', col=co2[as.numeric(unique(h[,1]))],
     xlab='survival time in days',ylab='survival probabilities',main="Allcancer")
legend('bottomleft', cex=0.6,text.width=0.4,c("TAD1_up","TAD1_normal","TAD1_down","TAD21","TAD22","TAD23","PRD1","PRD2","PRD3","DNAbinding1","DNAbinding2","DNAbinding3","Tet1","Tet2","Tet3","Basic1","Basic2","Basic3"), lty='solid',
       col=c("red4","red2","red","orange4","orange3","orange","maroon4", "maroon1","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue2","royalblue1"))
text(3000,1,cex=0.8,paste("p=",round(pchisq(dif2$chisq,1,lower.tail=F),4),sep=""))
for(j in as.numeric(unique(h[,1]))){
  text(2000,0.9-(j-1)*0.03,cex=0.8,j)
}
for(j in 1:length(as.numeric(unique(h[,1]))) ){
  text(3000,0.9-(j-1)*0.03,cex=0.8,dif2$n[j])
}
dev.off()


####对每个cancer27，把样本分为正常样本，疾病样本，突变样本，不突变样本，根据TP53的表达值做boxplot
##没有normal："LAML","LGG","OV","TGCT","UCS",27-7=20cancer
##1normal："SKCM",
##1mutation："PCPG",
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
h={}
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")){
  f1= read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  f2= read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3= read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  a1=f1[as.character(f1[,1])=="TP53",2:ncol(f1)]
  a2=f2[as.character(f2[,1])=="TP53",2:ncol(f2)]
  a3=f3[as.character(f3[,1])=="TP53",2:ncol(f3)] 
  a4=cbind(a2,a3)
  da=data.frame(rbind(cbind(class="normal",t(a1)),cbind(class="cancer",t(a4)),cbind(class="mutation",t(a2)),cbind(class="nomutation",t(a3))))
  da[,1]=factor(da[,1],c("normal","cancer","mutation","nomutation"))
  da[,2]=as.numeric(as.character(da[,2]))
  h=rbind(h,da)
  p<-ggplot(da, aes(x = da[,1],y = da[,2]))+geom_boxplot(aes(fill=class),outlier.colour = NA)+
    scale_fill_manual(values=c("grey","red","lightblue","lightgreen"))+geom_jitter()  
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Sample Type")+ylab("fpkm")+ggtitle(cancer) 
  p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="red",size=2.5,label=paste("p=",t.test(a1,a4)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(a1,a2)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(a1,a3)$p.value,sep=""))+
    annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(a2,a3)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust =5,vjust = 8,colour="lightblue",size=2.5,label=paste("p=",t.test(a4,a2)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 8, colour="lightgreen",size=2.5,label=paste("p=",t.test(a4,a3)$p.value,sep=""))
  ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_boxplot/",cancer,".pdf",sep=""),width = 10,height=9)
}
h[,1]=factor(h[,1],c("normal","cancer","mutation","nomutation"))
h[,2]=as.numeric(as.character(h[,2]))
p<-ggplot(h, aes(x = h[,1],y = h[,2]))+geom_boxplot(aes(fill=class),outlier.colour = NA)+
  scale_fill_manual(values=c("grey","red","lightblue","lightgreen"))+geom_jitter()  
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Sample Type")+ylab("fpkm")+ggtitle("Allcancer") 
p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="red",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="cancer",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
  annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="mutation",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 8,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 8, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_boxplot/Allcancer.pdf",sep=""),width = 10,height=9)



####对每个cancer27，把样本分为正常样本，疾病样本，突变样本，不突变样本，根据TP53的表达值做boxplot整合为一张图
##没有normal："LAML","LGG","OV","TGCT","UCS",27-7=20cancer
##1normal："SKCM",
##1mutation："PCPG",
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
h={}
h1=data.frame(Name=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                     "LUSC","PAAD","PRAD","READ","SARC","STAD","THCA","THYM","UCEC"))
rownames(h1)=h1$Name
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")){
  f1= read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  f2= read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3= read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  a1=f1[as.character(f1[,1])=="TP53",2:ncol(f1)]
  a2=f2[as.character(f2[,1])=="TP53",2:ncol(f2)]
  a3=f3[as.character(f3[,1])=="TP53",2:ncol(f3)] 
  a4=cbind(a2,a3)
  da=data.frame(rbind(cbind(class="normal",t(a1)),cbind(class="cancer",t(a4)),cbind(class="mutation",t(a2)),cbind(class="nomutation",t(a3))))
  da[,1]=factor(da[,1],c("normal","cancer","mutation","nomutation"))
  da[,2]=as.numeric(as.character(da[,2]))
  da[,3]=cancer
  h=rbind(h,da)
  h1[cancer,2]=median(as.numeric(a1[1,]))
}
h[,1]=factor(h[,1],c("normal","cancer","mutation","nomutation"))
h[,2]=as.numeric(as.character(h[,2]))
m=h
m[,3]="All"
n=rbind(m,h)
n[,3]=factor(n[,3],c("All",as.character(h1[order(h1[,2]),1])))
for(i in 1:length(c("All",as.character(h1[order(h1[,2]),1])))){
  n[which(as.character(n[,3])==c("All",as.character(h1[order(h1[,2]),1]))[i]),4]=i
}
n[which(as.character(n[,1])=="normal"),5]=-0.4
n[which(as.character(n[,1])=="cancer"),5]=-0.1
n[which(as.character(n[,1])=="mutation"),5]=0.1
n[which(as.character(n[,1])=="nomutation"),5]=0.4
p<-ggplot(n, aes(x = n[,3],y = n[,2],color = class))+geom_boxplot(aes(fill=class),alpha=0.8,lwd=0.1,width=0.8,position = position_dodge(1),outlier.colour = NA,outlier.shape = ".",outlier.size=0)+
  scale_fill_manual(values=c("mediumblue","violetred3","#1DADB8","#12B826"))+
  scale_color_manual(values=c("paleturquoise","pink","grey","springgreen1"))+
  geom_jitter(aes(V4 + V5,X132,colour = class,fill = class),
              position = position_jitter(width = 0.1,height = 0),
              alpha=0.2,size=0.001)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Cancer Type")+ylab("fpkm")+ggtitle("Allcancer") 
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/cancer_fpkm2.pdf",sep=""),width = 10,height=9)

n[,3]=factor(n[,3],c("All","LUAD","STAD","BLCA","PRAD","COAD","UCEC","READ","LIHC","SARC","ESCA","GBM","LUSC","CESC","BRCA",
                     "HNSC","PAAD","KIRC","KIRP","THCA","THYM"))
n1=n[which(as.character(n[,1])=="normal" | as.character(n[,1])=="cancer"),]
n2=n[which(as.character(n[,1])=="mutation" | as.character(n[,1])=="nomutation"),]
p<-ggplot(n1, aes(x = n1[,3],y = n1[,2],color = class))+geom_boxplot(aes(fill=class),lwd=0.1,width=0.8,position = position_dodge(1),outlier.colour = NA,outlier.shape = ".",outlier.size=0)+
  scale_fill_manual(values=c("violetred3","mediumblue"))+
  scale_color_manual(values=c("grey","grey"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Cancer Type")+ylab("fpkm")+ggtitle("Allcancer") 
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/cancer_fpkm2_1.pdf",sep=""),width = 10,height=9)
p<-ggplot(n2, aes(x = n2[,3],y = n2[,2],color = class))+geom_boxplot(aes(fill=class),lwd=0.1,width=0.8,position = position_dodge(1),outlier.colour = NA,outlier.shape = ".",outlier.size=0)+
  scale_fill_manual(values=c("#e87c25","#27727b"))+
  scale_color_manual(values=c("grey","grey"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Cancer Type")+ylab("fpkm")+ggtitle("Allcancer") 
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/cancer_fpkm2_2.pdf",sep=""),width = 10,height=9)


 





 

#为了看数据分布是否有偏，我们还可以增加均值与中值进行比较，主要用stat_summary把均值以菱形相展示。
#ggplot(birthwt, aes(x=factor(race), y=bwt)) + geom_boxplot() +
  #stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white")



####对每个cancer27，把样本分为正常样本，疾病样本，突变样本，不突变样本，根据TP53的表达值做density
##没有normal："LAML","LGG","OV","TGCT","UCS",27-7=20cancer
##1normal："SKCM",
##1mutation："PCPG",
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
h={}
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
                "LUSC","PAAD","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")){
  f1= read.table(paste(path,"p53_context_normal_fpkm/",cancer,"_p53_context_normal_fpkm.txt",sep=""),header=T,sep="\t")
  f2= read.table(paste(path,"p53_context_mutation_fpkm/",cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3= read.table(paste(path,"p53_context_nomutation_fpkm/",cancer,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  a1=f1[as.character(f1[,1])=="TP53",2:ncol(f1)]
  a2=f2[as.character(f2[,1])=="TP53",2:ncol(f2)]
  a3=f3[as.character(f3[,1])=="TP53",2:ncol(f3)] 
  a4=cbind(a2,a3)
  da=data.frame(rbind(cbind(class="normal",t(a1)),cbind(class="cancer",t(a4)),cbind(class="mutation",t(a2)),cbind(class="nomutation",t(a3))))
  da[,1]=factor(da[,1],c("normal","cancer","mutation","nomutation"))
  da[,2]=as.numeric(as.character(da[,2]))
  h=rbind(h,da)
  #p<-ggplot(da,aes(x = da[,2]))+geom_density(aes(fill=class),alpha=0.4)+scale_fill_manual(values=c("grey","red","lightblue","lightgreen"))
  #p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  #p <- p+xlab("fpkm")+ylab("density")+ggtitle(cancer) 
  #p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="red",size=2.5,label=paste("p=",t.test(a1,a4)$p.value,sep=""))+
  #  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(a1,a2)$p.value,sep=""))+
  #  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(a1,a3)$p.value,sep=""))+
  #  annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(a2,a3)$p.value,sep=""))+
  #  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 8,colour="lightblue",size=2.5,label=paste("p=",t.test(a4,a2)$p.value,sep=""))+
  #  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 8, colour="lightgreen",size=2.5,label=paste("p=",t.test(a4,a3)$p.value,sep=""))
  #ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_density/",cancer,".pdf",sep=""),width = 10,height=9)
  p<-ggplot(da,aes(x = da[,2],colour=class))+geom_density(aes(fill=class),alpha=0.4)+facet_grid(da[,1] ~ .) 
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("fpkm")+ylab("density")+ggtitle(cancer) 
  ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_density/1/",cancer,".pdf",sep=""),width = 10,height=9)
}
h[,1]=factor(h[,1],c("normal","cancer","mutation","nomutation"))
h[,2]=as.numeric(as.character(h[,2]))
#p<-ggplot(h,aes(x = h[,2]))+geom_density(aes(fill=class),alpha=0.4)+scale_fill_manual(values=c("grey","red","lightblue","lightgreen"))  
#p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
#p <- p+xlab("fpkm")+ylab("density")+ggtitle("Allcancer") 
#p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="red",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="cancer",2])$p.value,sep=""))+
#  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
#  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="normal",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
#  annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="mutation",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
#  annotate("text",x=Inf, y=Inf,hjust =5,vjust = 8,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
#  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 8, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))
#ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_density/Allcancer.pdf",sep=""),width = 10,height=9)
p<-ggplot(h,aes(x = h[,2],colour=class))+geom_density(aes(fill=class),alpha=0.4)+facet_grid(h[,1] ~ .) 
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("fpkm")+ylab("density")+ggtitle("Allcancer")  
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/4fpkm_density/1/Allcancer.pdf",sep=""),width = 10,height=9)



####对每个cancer22，把样本分为疾病样本(突变样本，不突变样本)，根据TP53的foldchange值做boxplot
##1mutation："PCPG"14,22-1=21cancer
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
h={}
for(cancer in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",cancer,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=f1[as.character(f1[,1])=="TP53",2:ncol(f1)]
  a2=f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]
  a3=f1["TP53",grep("_no",as.character(colnames(f1)))]
  
  da=data.frame(rbind(cbind(class="cancer",t(a1)),cbind(class="mutation",t(a2)),cbind(class="nomutation",t(a3))))
  da[,1]=factor(da[,1],c("cancer","mutation","nomutation"))
  da[,2]=as.numeric(as.character(da[,2]))
  h=rbind(h,da)
  p<-ggplot(da, aes(x = da[,1],y = da[,2]))+geom_boxplot(aes(fill=class),outlier.colour = NA)+
    scale_fill_manual(values=c("red","lightblue","lightgreen"))+geom_jitter()  
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("Sample Type")+ylab("foldchange")+ggtitle(gsub("_fc.txt","",cancer))
  p <- p+annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(a1,a2)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(a1,a3)$p.value,sep=""))+
    annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(a2,a3)$p.value,sep=""))
  ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/3foldchange_boxplot/",gsub("_fc.txt",".pdf",cancer),sep=""),width = 10,height=9)
}
h[,1]=factor(h[,1],c("cancer","mutation","nomutation"))
h[,2]=as.numeric(as.character(h[,2]))
p<-ggplot(h, aes(x = h[,1],y = h[,2]))+geom_boxplot(aes(fill=class),outlier.colour = NA)+
  scale_fill_manual(values=c("red","lightblue","lightgreen"))  
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("Sample Type")+ylab("foldchange")+ggtitle("Allcancer") 
p <- p+annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
  annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="mutation",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/3foldchange_boxplot/Allcancer.pdf",sep=""),width = 10,height=9)


####对每个cancer22，把样本分为疾病样本(突变样本，不突变样本)，根据TP53的foldchange值做density
##1mutation："PCPG"14,22-1=21cancer 
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
h={}
for(cancer in filename[c(1:13,15:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",cancer,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a1=f1[as.character(f1[,1])=="TP53",2:ncol(f1)]
  a2=f1["TP53",-grep("_no",as.character(colnames(f1)))][-1]
  a3=f1["TP53",grep("_no",as.character(colnames(f1)))]
  
  da=data.frame(rbind(cbind(class="cancer",t(a1)),cbind(class="mutation",t(a2)),cbind(class="nomutation",t(a3))))
  da[,1]=factor(da[,1],c("cancer","mutation","nomutation"))
  da[,2]=as.numeric(as.character(da[,2]))
  h=rbind(h,da)
  p<-ggplot(da,aes(x = da[,2]))+geom_density(aes(fill=class),alpha=0.5)+scale_fill_manual(values=c("red","lightblue","lightgreen"))  
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
  p <- p+xlab("foldchange")+ylab("density")+ggtitle(gsub("_fc.txt","",cancer))
  p <- p+annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(a1,a2)$p.value,sep=""))+
    annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(a1,a3)$p.value,sep=""))+
    annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(a2,a3)$p.value,sep=""))
  ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/3foldchange_density/",gsub("_fc.txt",".pdf",cancer),sep=""),width = 10,height=9)
}
h[,1]=factor(h[,1],c("cancer","mutation","nomutation"))
h[,2]=as.numeric(as.character(h[,2]))
p<-ggplot(h,aes(x = h[,2]))+geom_density(aes(fill=class),alpha=0.5)+scale_fill_manual(values=c("red","lightblue","lightgreen"))  
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=6,colour="black"))
p <- p+xlab("foldchange")+ylab("density")+ggtitle("Allcancer") 
p <- p+annotate("text",x=Inf, y=Inf,hjust =5,vjust = 5,colour="lightblue",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="mutation",2])$p.value,sep=""))+
  annotate("text",x=Inf, y=Inf,hjust = 2,vjust = 5, colour="lightgreen",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="cancer",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))+
  annotate("text",x=-Inf, y=-Inf,hjust = -6,vjust = -5, colour="purple",size=2.5,label=paste("p=",t.test(h[as.character(h[,1])=="mutation",2],h[as.character(h[,1])=="nomutation",2])$p.value,sep=""))
ggsave(file=paste(path,"normal_cancer_mu_nomu_boxplot/3foldchange_density/Allcancer.pdf",sep=""),width = 10,height=9)



####对Allcancer_fc_updown.txt的17cancer，分突变、不突变样本，
####看TP53基因上调、下调、正常的比例，画饼图pie
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Allcancer_fc_updown.txt",sep=""),header=T,sep="\t")
for(i in 1:nrow(f)){
  pdf(paste(path,"Mu_NoMu_pie/",as.character(f[i,1]),".pdf",sep=""))
  opar<-par(mfrow=c(1,3), mar=c(4,4,4,4))   
  x1 <- c(f[i,5]/f[i,2],f[i,7]/f[i,2],(f[i,2]-f[i,5]-f[i,7])/f[i,2])
  x2 <- c(f[i,12]/f[i,9],f[i,14]/f[i,9],(f[i,9]-f[i,12]-f[i,14])/f[i,9])
  x3 <- c(f[i,19]/f[i,16],f[i,21]/f[i,16],(f[i,16]-f[i,19]-f[i,21])/f[i,16])
  label <- c("up","down", "normal")
  pie(x1, labels=paste(round(x1,4)*100, "%",x1*f[i,2],sep = ""), main=paste(as.character(f[i,1]),"_Mu",sep=""), col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
  legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
  pie(x2, labels=paste(round(x2,4)*100, "%",x2*f[i,9],sep = ""), main=paste(as.character(f[i,1]),"_NoMu",sep=""), col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
  legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
  pie(x3, labels=paste(round(x3,4)*100, "%",x3*f[i,16],sep = ""), main=paste(as.character(f[i,1]),"_Mu/NoMu",sep=""), col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
  legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
  par(opar)  #在活动设备中返回所有图形参数和他们的值
  dev.off()
}
pdf(paste(path,"Mu_NoMu_pie/Allcancer.pdf",sep=""))
opar<-par(mfrow=c(1,3), mar=c(4,4,4,4))   
x4 <- c(sum(f[,5])/sum(f[,2]),sum(f[,7])/sum(f[,2]),(sum(f[,2])-sum(f[,5])-sum(f[,7]))/sum(f[,2]))
x5 <- c(sum(f[,12])/sum(f[,9]),sum(f[,14])/sum(f[,9]),(sum(f[,9])-sum(f[,12])-sum(f[,14]))/sum(f[,9]))
x6 <- c(sum(f[,19])/sum(f[,16]),sum(f[,21])/sum(f[,16]),(sum(f[,16])-sum(f[,19])-sum(f[,21]))/sum(f[,16]))
label <- c("up","down", "normal")
pie(x4, labels=paste(round(x4,4)*100, "%",x4*sum(f[,2]),sep = ""), main="Allcancer_Mu", col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
pie(x5, labels=paste(round(x5,4)*100, "%",x5*sum(f[,9]), sep = ""), main="Allcancer_NoMu", col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
pie(x6, labels=paste(round(x6,4)*100, "%",x6*sum(f[,16]),sep = ""), main="Allcancer_Mu/NoMu", col=c("pink","lightgreen","grey"), clockwise = TRUE) #这是按顺时针方向来绘制饼图
legend("topleft",label, cex=0.8, fill=c("pink","lightgreen","grey"))
par(opar)  #在活动设备中返回所有图形参数和他们的值
dev.off()
#library(plotrix)
#pie3D(x1,labels=label,explode=0.1) #绘制3D 饼图



install.packages('easyGgplot2')
library(easyGgplot2)
ggplot2.multiplot(p2,p3,p4,p5,p6,p7,p8,p9,p10,cols=3) 
multiplot(p1, p2, p3, p4, cols=2)
# Multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

 
####十二1、对每个cancer22,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
##CESC3\GBM6\KIRC8\KIRP9\PAAD13\PCPG14\SKCM18\THCA20\THYM21
#library(ggrepel)文本标签重叠很难看清楚。使用ggrepel可以使绘图区域内的文本标签相互分离。geom_text_repe()
#m=dd[which(dd[,2]>100 & dd[,2]< 200),]
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",gsub("_no","",colnames(a[1,which(a[1,]>=1)])))
  a2=gsub("\\.","-",gsub("_no","",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)])))
  a3=gsub("\\.","-",gsub("_no","",colnames(a[1,which(a[1,]<= -1)])))
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  dd=dd[order(dd[,1]),]
  Expression_Value=factor(dd[,1])  
  p <- ggplot(dd, aes(x=dd[,5], y=dd[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    #+facet_grid(dd[,1] ~ .)
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","",i))  #+xlim(0,400) 
  p <- p+annotate("text",x=dd[,5],y=dd[,4], colour="black",size=2,label=dd[,3],angle=-90)+
    annotate("text",x=dd[,5],y=-0.2, colour="black",size=2,label=dd[,2],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/1/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)
}

h1[,1]=as.numeric(as.character(h1[,1]))
h1[,2]=as.character(h1[,2])
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2[,2]=as.character(h2[,2])
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3[,2]=as.character(h3[,2])
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h=h[order(h[,2]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[order(h[,1]),]
Expression_Value1=factor(h[,1])  
p <- ggplot(h, aes(x=h[,5], y=h[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))    #+facet_grid(dd[,1] ~ .)
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer")  #+xlim(0,400) 
p <- p+annotate("text",x=h[,5],y=h[,4], colour="black",size=2,label=h[,3],angle=-90)+
  annotate("text",x=h[,5],y=-0.2, colour="black",size=2,label=h[,2],angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/1/Allcancer.pdf",sep=""),width = 10,height=9)

 
 

####十二2、对每个cancer22,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
##position分段0-100,100-200,200-300,300-400  Allcancer除去\KIRC8\KIRP9\
#library(ggrepel)文本标签重叠很难看清楚。使用ggrepel可以使绘图区域内的文本标签相互分离。geom_text_repe()
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename[c(1:7,10:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  m1=dd[which(dd[,2]<100),]
  m2=dd[which(dd[,2]>=100 & dd[,2]< 200),]
  m3=dd[which(dd[,2]>=200 & dd[,2]< 300),]
  m4=dd[which(dd[,2]>=300),]
  Expression_Value1=factor(m1[,1]) 
  Expression_Value2=factor(m2[,1]) 
  Expression_Value3=factor(m3[,1]) 
  Expression_Value4=factor(m4[,1])
  p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_1",i))   
  p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
    annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/",gsub("_fc.txt","_1.pdf",i),sep=""),width = 10,height=9)

  p <- ggplot(m2, aes(x=m2[,5], y=m2[,4], fill=Expression_Value2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_2",i))   
  p <- p+annotate("text",x=m2[,5],y=m2[,4], colour="black",size=2,label=m2[,4],angle=-90)+
    annotate("text",x=m2[,5],y=-0.2, colour="black",size=2,label=paste(m2[,2],m2[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/",gsub("_fc.txt","_2.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m3, aes(x=m3[,5], y=m3[,4], fill=Expression_Value3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_3",i))   
  p <- p+annotate("text",x=m3[,5],y=m3[,4], colour="black",size=2,label=m3[,4],angle=-90)+
    annotate("text",x=m3[,5],y=-0.2, colour="black",size=2,label=paste(m3[,2],m3[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/",gsub("_fc.txt","_3.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m4, aes(x=m4[,5], y=m4[,4], fill=Expression_Value4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_4",i))   
  p <- p+annotate("text",x=m4[,5],y=m4[,4], colour="black",size=2,label=m4[,4],angle=-90)+
    annotate("text",x=m4[,5],y=-0.2, colour="black",size=2,label=paste(m4[,2],m4[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/",gsub("_fc.txt","_4.pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
n1=h[which(h[,2]<100),]
n2=h[which(h[,2]>=100 & h[,2]< 135),]
n3=h[which(h[,2]>=135 & h[,2]< 160),]
n4=h[which(h[,2]>=160 & h[,2]< 190),]
n5=h[which(h[,2]>=190 & h[,2]< 220),]
n6=h[which(h[,2]>=220 & h[,2]< 250),]
n7=h[which(h[,2]>=250 & h[,2]< 275),]
n8=h[which(h[,2]>=275 & h[,2]< 300),]
n9=h[which(h[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_1")   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.2, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_2")   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.2, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_3")   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.2, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_4")   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.2, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_4.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_5")   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.2, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_5.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_6")   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.2, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_6.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_7")   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.2, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_7.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_8")   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.2, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_8.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_9")   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.2, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/2_1/Allcancer_9.pdf",sep=""),width = 10,height=9)





####十二3、对每个cancer22,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
##position分段0-100,100-200,200-300,300-400  Allcancer除去CESC3\KIRC8\KIRP9\PRAD15\READ16\SARC17  
##把突变数为1的删去（第四列）存为3文件夹下 22-6=16cancer
#library(ggrepel)文本标签重叠很难看清楚。使用ggrepel可以使绘图区域内的文本标签相互分离。geom_text_repe()
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename[c(1:2,4:7,10:14,18:22)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  m1=dd[which(dd[,4]!=1),]
  Expression_Value1=factor(m1[,1]) 
  p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","",i))   
  p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
    annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1[,2]=as.character(h1[,2])
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2[,2]=as.character(h2[,2])
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3[,2]=as.character(h3[,2])
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h[,3]=gsub("\\+","",as.character(h[,3]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[which(h[,4]!=1),]
n1=h[which(h[,2]<100),]
n2=h[which(h[,2]>=100 & h[,2]< 135),]
n3=h[which(h[,2]>=135 & h[,2]< 160),]
n4=h[which(h[,2]>=160 & h[,2]< 190),]
n5=h[which(h[,2]>=190 & h[,2]< 220),]
n6=h[which(h[,2]>=220 & h[,2]< 250),]
n7=h[which(h[,2]>=250 & h[,2]< 275),]
n8=h[which(h[,2]>=275 & h[,2]< 300),]
n9=h[which(h[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_1")   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.2, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_2")   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.2, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_3")   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.2, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_4")   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.2, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_4.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_5")   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.2, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_5.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_6")   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.2, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_6.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_7")   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.2, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_7.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_8")   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.2, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_8.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_9")   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.2, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/3_1/Allcancer_9.pdf",sep=""),width = 10,height=9)




####十二4、对每个cancer27,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
##用的是p53_context_fpkm_foldchange_p53%WTmean
##position分段0-100,100-200,200-300,300-400    
##把突变数为1的删去（第四列）存为4_WTmean文件夹下 
##Allcancer除去KIRC8\KIRP9\LAML10\PCPG17\PRAD18\TGCT23\UCS27    27-7=20cancer 
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename[c(1:7,11:16,19:22,24:26)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  colnames(f1)<-gsub("_p53/WT_mean","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  m1=dd[which(dd[,4]!=1),]
  Expression_Value1=factor(m1[,1]) 
  p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","",i))   
  p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
    annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1[,2]=as.character(h1[,2])
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2[,2]=as.character(h2[,2])
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3[,2]=as.character(h3[,2])
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h[,3]=gsub("\\+","",as.character(h[,3]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[which(h[,4]!=1),]
n1=h[which(h[,2]<100),]
n2=h[which(h[,2]>=100 & h[,2]< 135),]
n3=h[which(h[,2]>=135 & h[,2]< 160),]
n4=h[which(h[,2]>=160 & h[,2]< 190),]
n5=h[which(h[,2]>=190 & h[,2]< 220),]
n6=h[which(h[,2]>=220 & h[,2]< 250),]
n7=h[which(h[,2]>=250 & h[,2]< 275),]
n8=h[which(h[,2]>=275 & h[,2]< 300),]
n9=h[which(h[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_1")   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.2, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_2")   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.2, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_3")   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.2, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_4")   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.2, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_4.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_5")   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.2, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_5.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_6")   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.2, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_6.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_7")   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.2, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_7.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_8")   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.2, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_8.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_9")   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.2, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/4_WTmean1/Allcancer_9.pdf",sep=""),width = 10,height=9)




####十二5、对每个cancer22\allcancer,对每个cancer27\allcancer,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
####fc\
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  Expression_Value=factor(dd[,1])  
  p <- ggplot(dd, aes(x=dd[,5], y=dd[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))     
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","",i))  
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[order(h[,1]),]
Expression_Value1=factor(h[,1])  
p <- ggplot(h, aes(x=h[,5], y=h[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer")  
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc/Allcancer.pdf",sep=""),width = 10,height=9)

####fc_WT
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
for(i in filename[c(1:7,11:16,19:22,24:26)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  colnames(f1)<-gsub("_p53/WT_mean","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normal",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normal",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="down",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  Expression_Value=factor(dd[,1]) 
  p <- ggplot(dd, aes(x=dd[,5], y=dd[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","green"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","",i))   
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1[,2]=as.character(h1[,2])
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2[,2]=as.character(h2[,2])
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normal",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3[,2]=as.character(h3[,2])
h3=h3[order(h3[,1]),]
h33=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33))
h[,2]=as.numeric(as.character(h[,2]))
h[,3]=gsub("\\+","",as.character(h[,3]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
ExpressionValue1=factor(h[,1])
p <- ggplot(h, aes(x=h[,5], y=h[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","green"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer")   
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT/Allcancer.pdf",sep=""),width = 10,height=9)



####十二6、对每个cancer22\allcancer,对每个cancer27\allcancer,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
####在5的基础上把normal分为normalup（0~1）、normaldown（-1~0），标出normaldown\normal的比例
####fc\整个长度为文件夹fc_normalupdown、分段为fc_normalupdown\1\
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
h4={}
for(i in filename){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]>= 0  & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 0)]))
  a4=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normalup",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normalup",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0 ){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="normaldown",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0 ){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="normaldown",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  
  if(length(a4)!=0){
    d4=data.frame(b1[b1[,1] %in% a4,2:3])
    if(ncol(d4)==1 & nrow(d4)!=0){
      rownames(d4)=c("X1","X2")
      h4=rbind(h4,t(d4))
      d44=data.frame(v1="down",
                     v2=as.numeric(as.character(d4[1,1])),
                     v3=as.character(d4[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d4)==2 & nrow(d4)!=0){
      h4=rbind(h4,d4)
      d4[,1]=as.numeric(as.character(d4[,1]))
      d4[,2]=as.character(d4[,2])
      d4=d4[order(d4[,1]),]
      d44=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v4=as.numeric(table(paste(d4[,1],d4[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33,d44))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  Expression_Value=factor(dd[,1])  
  p <- ggplot(dd, aes(x=dd[,5], y=dd[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)

  m1=dd[which(dd[,2]<100),]
  m2=dd[which(dd[,2]>=100 & dd[,2]< 200),]
  m3=dd[which(dd[,2]>=200 & dd[,2]< 300),]
  m4=dd[which(dd[,2]>=300),]
  Expression_Value1=factor(m1[,1]) 
  Expression_Value2=factor(m2[,1]) 
  Expression_Value3=factor(m3[,1]) 
  Expression_Value4=factor(m4[,1])
  p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
    annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/",gsub("_fc.txt","_1.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m2, aes(x=m2[,5], y=m2[,4], fill=Expression_Value2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m2[,5],y=m2[,4], colour="black",size=2,label=m2[,4],angle=-90)+
    annotate("text",x=m2[,5],y=-0.1, colour="black",size=2,label=paste(m2[,2],m2[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/",gsub("_fc.txt","_2.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m3, aes(x=m3[,5], y=m3[,4], fill=Expression_Value3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m3[,5],y=m3[,4], colour="black",size=2,label=m3[,4],angle=-90)+
    annotate("text",x=m3[,5],y=-0.1, colour="black",size=2,label=paste(m3[,2],m3[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/",gsub("_fc.txt","_3.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m4, aes(x=m4[,5], y=m4[,4], fill=Expression_Value4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m4[,5],y=m4[,4], colour="black",size=2,label=m4[,4],angle=-90)+
    annotate("text",x=m4[,5],y=-0.1, colour="black",size=2,label=paste(m4[,2],m4[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/",gsub("_fc.txt","_4.pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normalup",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3=h3[order(h3[,1]),]
h33=data.frame(v1="normaldown",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h4[,1]=as.numeric(as.character(h4[,1]))
h4=h4[order(h4[,1]),]
h44=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v4=as.numeric(table(paste(h4[,1],h4[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33,h44))
h[,2]=as.numeric(as.character(h[,2]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[order(h[,1]),]
Expression_Value=factor(h[,1])  
p <- ggplot(h, aes(x=h[,5], y=h[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/Allcancer.pdf",sep=""),width = 10,height=9)

n1=h[which(h[,2]<100),]
n2=h[which(h[,2]>=100 & h[,2]< 135),]
n3=h[which(h[,2]>=135 & h[,2]< 160),]
n4=h[which(h[,2]>=160 & h[,2]< 190),]
n5=h[which(h[,2]>=190 & h[,2]< 220),]
n6=h[which(h[,2]>=220 & h[,2]< 250),]
n7=h[which(h[,2]>=250 & h[,2]< 275),]
n8=h[which(h[,2]>=275 & h[,2]< 300),]
n9=h[which(h[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.1, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.1, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.1, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.1, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_4.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.1, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_5.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.1, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_6.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.1, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_7.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.1, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_8.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.1, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/1/Allcancer_9.pdf",sep=""),width = 10,height=9)




####十二6、对每个cancer22\allcancer,对每个cancer27\allcancer,样本对应的hotspot的上下调正常的突变个数+具体突变类型bar
####在5的基础上把normal分为normalup（0~1）、normaldown（-1~0），标出normaldown\normal的比例
####fc\整个长度为文件夹fc_normalupdown、分段为fc_normalupdown\2\  #GBM6、SARC17不分段，PRAD15没有重复突变位点不画图
####ESCA5\SARC17没有ddmix,GBM6PRAD15
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
h4={}
for(i in filename[18:22]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-c(strsplit(as.character(f1[1,1])," ")[[1]],strsplit(as.character(f1[1,2])," ")[[1]])
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",-grep("_no",colnames(f1))][-1]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]>= 0  & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 0)]))
  a4=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normalup",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normalup",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0 ){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="normaldown",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0 ){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="normaldown",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  
  if(length(a4)!=0){
    d4=data.frame(b1[b1[,1] %in% a4,2:3])
    if(ncol(d4)==1 & nrow(d4)!=0){
      rownames(d4)=c("X1","X2")
      h4=rbind(h4,t(d4))
      d44=data.frame(v1="down",
                     v2=as.numeric(as.character(d4[1,1])),
                     v3=as.character(d4[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d4)==2 & nrow(d4)!=0){
      h4=rbind(h4,d4)
      d4[,1]=as.numeric(as.character(d4[,1]))
      d4[,2]=as.character(d4[,2])
      d4=d4[order(d4[,1]),]
      d44=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v4=as.numeric(table(paste(d4[,1],d4[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33,d44))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  
  dd[,6]=paste(dd[,2],dd[,3],sep="_")  #dd分为ddup、dddown、ddmix，dd分为dd1（4（567）、8）、dd2，
  dd1=subset(dd,dd[,6] %in% unique(dd[,6][duplicated(dd[,6])]))    #重复突变位点
  dd2=subset(dd,!(dd[,6] %in% unique(dd[,6][duplicated(dd[,6])]))) #没有重复突变位点分为up、normalup、down、normaldown
  dd3=data.frame(table(dd1[,6]))  #重复突变位点次数
  dd4=subset(dd1,dd1[,6] %in% as.character(dd3[dd3[,2]==2,1]))  #重复突变位点次数为2
  dd5=data.frame() #重复突变位点次数为2且up
  dd6=data.frame() #重复突变位点次数为2且down
  dd7=data.frame() #重复突变位点次数为2且mix
  for(r in seq(1,nrow(dd4),2)){
    if((dd4[r,6]==dd4[r+1,6]) & (as.character(dd4[r,1])=="up") & (as.character(dd4[r+1,1])=="normalup")){
      dd5=rbind(dd5,dd4[r,])
      dd5=rbind(dd5,dd4[r+1,])
    }
    else if((dd4[r,6]==dd4[r+1,6]) & (as.character(dd4[r,1])=="normaldown") & (as.character(dd4[r+1,1])=="down")){
      dd6=rbind(dd6,dd4[r,])
      dd6=rbind(dd6,dd4[r+1,])
    }
    else{
      dd7=rbind(dd7,dd4[r,])
      dd7=rbind(dd7,dd4[r+1,])
    }
  }
  dd8=subset(dd1,dd1[,6] %in% as.character(dd3[dd3[,2]!=2,1]))  #重复突变位点次数不为2且mix
  #ddup为up、normalup（dd2）或upnormalup（dd5）
  ddup=rbind(subset(dd2,as.character(dd2[,1])=="up"|as.character(dd2[,1])=="normalup"),dd5)
  #dddown为down、normaldown（dd2）或downnormaldown（dd6）
  dddown=rbind(subset(dd2,as.character(dd2[,1])=="down"|as.character(dd2[,1])=="normaldown"),dd6)
  #ddmix为混合模式 （dd7+dd8）
  ddmix= rbind(dd7,dd8)
  Expression_Value1=factor(ddup[,1]) 
  Expression_Value2=factor(dddown[,1])
  Expression_Value3=factor(ddmix[,1])
  p <- ggplot(ddup, aes(x=ddup[,5], y=ddup[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey"))     
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_up",i)) +
    annotate("text",x=ddup[,5],y=-0.1, colour="black",size=2,label=ddup[,6],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_up.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(dddown, aes(x=dddown[,5], y=dddown[,4], fill=Expression_Value2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("purple","lightgreen"))     
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_down",i)) +
    annotate("text",x=dddown[,5],y=-0.1, colour="black",size=2,label=dddown[,6],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_down.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(ddmix, aes(x=ddmix[,5], y=ddmix[,4], fill=Expression_Value3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_mix",i)) +
    annotate("text",x=ddmix[,5],y=-0.1, colour="black",size=2,label=ddmix[,6],angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_mix.pdf",i),sep=""),width = 10,height=9)
  
  m1=ddup[which(ddup[,2]<260),]
  m2=ddup[which(ddup[,2]>=260),]
  m3=dddown[which(dddown[,2]<260),]
  m4=dddown[which(dddown[,2]>=260),]
  m5=ddmix[which(ddmix[,2]<260),]
  m6=ddmix[which(ddmix[,2]>=260),]
  Expression_Value11=factor(m1[,1]) 
  Expression_Value12=factor(m2[,1]) 
  Expression_Value21=factor(m3[,1]) 
  Expression_Value22=factor(m4[,1]) 
  Expression_Value31=factor(m5[,1]) 
  Expression_Value32=factor(m6[,1]) 
  p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value11)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_up_1",i))  
  p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
    annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_up_1.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m2, aes(x=m2[,5], y=m2[,4], fill=Expression_Value12)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_up_2",i))  
  p <- p+annotate("text",x=m2[,5],y=m2[,4], colour="black",size=2,label=m2[,4],angle=-90)+
    annotate("text",x=m2[,5],y=-0.1, colour="black",size=2,label=paste(m2[,2],m2[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_up_2.pdf",i),sep=""),width = 10,height=9)

  p <- ggplot(m3, aes(x=m3[,5], y=m3[,4], fill=Expression_Value21)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_down_1",i))  
  p <- p+annotate("text",x=m3[,5],y=m3[,4], colour="black",size=2,label=m3[,4],angle=-90)+
    annotate("text",x=m3[,5],y=-0.1, colour="black",size=2,label=paste(m3[,2],m3[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_down_1.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m4, aes(x=m4[,5], y=m4[,4], fill=Expression_Value22)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_down_2",i))  
  p <- p+annotate("text",x=m4[,5],y=m4[,4], colour="black",size=2,label=m4[,4],angle=-90)+
    annotate("text",x=m4[,5],y=-0.1, colour="black",size=2,label=paste(m4[,2],m4[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_down_2.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m5, aes(x=m5[,5], y=m5[,4], fill=Expression_Value31)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_mix_1",i))  
  p <- p+annotate("text",x=m5[,5],y=m5[,4], colour="black",size=2,label=m5[,4],angle=-90)+
    annotate("text",x=m5[,5],y=-0.1, colour="black",size=2,label=paste(m5[,2],m5[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_mix_1.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m6, aes(x=m6[,5], y=m6[,4], fill=Expression_Value32)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(gsub("_fc.txt","_mix_2",i))  
  p <- p+annotate("text",x=m6[,5],y=m6[,4], colour="black",size=2,label=m6[,4],angle=-90)+
    annotate("text",x=m6[,5],y=-0.1, colour="black",size=2,label=paste(m6[,2],m6[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/",gsub("_fc.txt","_mix_2.pdf",i),sep=""),width = 10,height=9)
  
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normalup",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3=h3[order(h3[,1]),]
h33=data.frame(v1="normaldown",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h4[,1]=as.numeric(as.character(h4[,1]))
h4=h4[order(h4[,1]),]
h44=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v4=as.numeric(table(paste(h4[,1],h4[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33,h44))
h[,2]=as.numeric(as.character(h[,2]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
h=h[order(h[,1]),]
h[,6]=paste(h[,2],h[,3],sep="_")  #h分为hhup、hhdown、hhmix，h分为hh1（4（567）、8）=hh1、hh2=hh2，
hh1=subset(h,h[,6] %in% unique(h[,6][duplicated(h[,6])]))    #重复突变位点
hh2=subset(h,!(h[,6] %in% unique(h[,6][duplicated(h[,6])]))) #没有重复突变位点分为up、normalup、down、normaldown
hh3=data.frame(table(hh1[,6]))  #重复突变位点次数
hh4=subset(hh1,hh1[,6] %in% as.character(hh3[hh3[,2]==2,1]))  #重复突变位点次数为2
hh5=data.frame() #重复突变位点次数为2且up
hh6=data.frame() #重复突变位点次数为2且down
hh7=data.frame() #重复突变位点次数为2且mix
for(r in seq(1,nrow(hh4),2)){
  if((hh4[r,6]==hh4[r+1,6]) & (as.character(hh4[r,1])=="up") & (as.character(hh4[r+1,1])=="normalup")){
    hh5=rbind(hh5,hh4[r,])
    hh5=rbind(hh5,hh4[r+1,])
  }
  else if((hh4[r,6]==hh4[r+1,6]) & (as.character(hh4[r,1])=="normaldown") & (as.character(hh4[r+1,1])=="down")){
    hh6=rbind(hh6,hh4[r,])
    hh6=rbind(hh6,hh4[r+1,])
  }
  else{
    hh7=rbind(hh7,hh4[r,])
    hh7=rbind(hh7,hh4[r+1,])
  }
}
hh8=subset(hh1,hh1[,6] %in% as.character(hh3[hh3[,2]!=2,1]))  #重复突变位点次数不为2且mix
hhup=rbind(subset(hh2,as.character(hh2[,1])=="up"|as.character(hh2[,1])=="normalup"),hh5)
hhdown=rbind(subset(hh2,as.character(hh2[,1])=="down"|as.character(hh2[,1])=="normaldown"),hh6)
hhmix= rbind(hh7,hh8)
ExpressValue1=factor(hhup[,1]) 
ExpressValue2=factor(hhdown[,1])
ExpressValue3=factor(hhmix[,1])

p <- ggplot(hhup, aes(x=hhup[,5], y=hhup[,4], fill=ExpressValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_up") 
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_up.pdf",sep=""),width = 10,height=9)

p <- ggplot(hhdown, aes(x=hhdown[,5], y=hhdown[,4], fill=ExpressValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_down") 
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_down.pdf",sep=""),width = 10,height=9)

p <- ggplot(hhmix, aes(x=hhmix[,5], y=hhmix[,4], fill=ExpressValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_mix") 
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_mix.pdf",sep=""),width = 10,height=9)

n1=hhup[which(hhup[,2]<150),]
n2=hhup[which(hhup[,2]>=150 & hhup[,2]< 300),]
n3=hhup[which(hhup[,2]>=300),]
n4=hhdown[which(hhdown[,2]<150),]
n5=hhdown[which(hhdown[,2]>=150 & hhdown[,2]< 300),]
n6=hhdown[which(hhdown[,2]>=300),]
n7=hhmix[which(hhmix[,2]<150),]
n8=hhmix[which(hhmix[,2]>=150 & hhmix[,2]< 300),]
n9=hhmix[which(hhmix[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_up_1")   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.1, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_up_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_up_2")   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.1, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_up_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_up_3")   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.1, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_up_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_down_1")   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.1, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_down_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_down_2")   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.1, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_down_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_down_3")   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.1, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_down_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_mix_1")   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.1, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_mix_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_mix_2")   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.1, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_mix_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle("Allcancer_mix_3")   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.1, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_normalupdown/2/Allcancer_mix_3.pdf",sep=""),width = 10,height=9)



####fc_WT
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
h1={}
h2={}
h3={}
h4={}
for(i in filename[c(1:7,11:16,19:22,24:26)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  colnames(f1)<-gsub("_p53/WT_mean","",colnames(f1))
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  a=f1["TP53",2:ncol(f1)]
  a1=gsub("\\.","-",colnames(a[1,which(a[1,]>=1)]))
  a2=gsub("\\.","-",colnames(a[1,which(a[1,]>= 0 & a[1,]< 1)]))
  a3=gsub("\\.","-",colnames(a[1,which(a[1,]> -1 & a[1,]< 0)]))
  a4=gsub("\\.","-",colnames(a[1,which(a[1,]<= -1)]))
  b=subset(f,as.character(f[,5])==gsub("_fc.txt","",i))
  b1=na.omit(cbind(as.character(b[,1]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(b$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(b$HGVSc),4))))
  b1[,3]=gsub("\\+","",b1[,3])
  if(length(a1)!=0){
    d1=data.frame(b1[b1[,1] %in% a1,2:3])
    if(ncol(d1)==1 & nrow(d1)!=0){
      rownames(d1)=c("X1","X2")
      h1=rbind(h1,t(d1))
      d11=data.frame(v1="up",
                     v2=as.numeric(as.character(d1[1,1])),
                     v3=as.character(d1[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d1)==2 & nrow(d1)!=0){
      h1=rbind(h1,d1)
      d1[,1]=as.numeric(as.character(d1[,1]))
      d1[,2]=as.character(d1[,2])
      d1=d1[order(d1[,1]),]
      d11=data.frame(v1="up",
                     v2=gsub("_.*","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d1[,1],d1[,2],sep="_")))),
                     v4=as.numeric(table(paste(d1[,1],d1[,2],sep="_"))))
    }
  }
  if(length(a2)!=0){
    d2=data.frame(b1[b1[,1] %in% a2,2:3])
    if(ncol(d2)==1 & nrow(d2)!=0 ){
      rownames(d2)=c("X1","X2")
      h2=rbind(h2,t(d2))
      d22=data.frame(v1="normalup",
                     v2=as.numeric(as.character(d2[1,1])),
                     v3=as.character(d2[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d2)==2 & nrow(d2)!=0 ){
      h2=rbind(h2,d2)
      d2[,1]=as.numeric(as.character(d2[,1]))
      d2[,2]=as.character(d2[,2])
      d2=d2[order(d2[,1]),]
      d22=data.frame(v1="normalup",
                     v2=gsub("_.*","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d2[,1],d2[,2],sep="_")))),
                     v4=as.numeric(table(paste(d2[,1],d2[,2],sep="_"))))
    }
  }
  if(length(a3)!=0){
    d3=data.frame(b1[b1[,1] %in% a3,2:3])
    if(ncol(d3)==1 & nrow(d3)!=0){
      rownames(d3)=c("X1","X2")
      h3=rbind(h3,t(d3))
      d33=data.frame(v1="normaldown",
                     v2=as.numeric(as.character(d3[1,1])),
                     v3=as.character(d3[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d3)==2 & nrow(d3)!=0){
      h3=rbind(h3,d3)
      d3[,1]=as.numeric(as.character(d3[,1]))
      d3[,2]=as.character(d3[,2])
      d3=d3[order(d3[,1]),]
      d33=data.frame(v1="normaldown",
                     v2=gsub("_.*","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d3[,1],d3[,2],sep="_")))),
                     v4=as.numeric(table(paste(d3[,1],d3[,2],sep="_"))))
    }
  }
  
  if(length(a4)!=0){
    d4=data.frame(b1[b1[,1] %in% a4,2:3])
    if(ncol(d4)==1 & nrow(d4)!=0){
      rownames(d4)=c("X1","X2")
      h4=rbind(h4,t(d4))
      d44=data.frame(v1="down",
                     v2=as.numeric(as.character(d4[1,1])),
                     v3=as.character(d4[2,1]),
                     v4=as.numeric(1))
    }
    if(ncol(d4)==2 & nrow(d4)!=0){
      h4=rbind(h4,d4)
      d4[,1]=as.numeric(as.character(d4[,1]))
      d4[,2]=as.character(d4[,2])
      d4=d4[order(d4[,1]),]
      d44=data.frame(v1="down",
                     v2=gsub("_.*","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v3=gsub(".*_","",names(table(paste(d4[,1],d4[,2],sep="_")))),
                     v4=as.numeric(table(paste(d4[,1],d4[,2],sep="_"))))
    }
  }
  dd=na.omit(rbind(d11,d22,d33,d44))
  dd[,2]=as.numeric(as.character(dd[,2]))
  dd=dd[order(dd[,2],dd[,3]),]
  dd[,5]=dd[,2]
  for(j in 1:(nrow(dd)-1)){
    if(dd[j,2]==dd[j+1,2]){
      dd[j+1,5]=dd[j,5]+1
    }
    if((dd[j,2]!=dd[j+1,2]) & (dd[j+1,5]<=dd[j,5])){
      dd[j+1,5]=dd[j,5]+1
    }
  }
  Expression_Value=factor(dd[,1]) 
  p <- ggplot(dd, aes(x=dd[,5], y=dd[,4], fill=Expression_Value)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/",gsub("_fc.txt",".pdf",i),sep=""),width = 10,height=9)

  m1=dd[which(dd[,2]<100),]
  m2=dd[which(dd[,2]>=100 & dd[,2]< 200),]
  m3=dd[which(dd[,2]>=200 & dd[,2]< 300),]
  m4=dd[which(dd[,2]>=300),]
  Expression_Value1=factor(m1[,1]) 
  Expression_Value2=factor(m2[,1]) 
  Expression_Value3=factor(m3[,1]) 
  Expression_Value4=factor(m4[,1])
  if(nrow(m1)!=0){
    p <- ggplot(m1, aes(x=m1[,5], y=m1[,4], fill=Expression_Value1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
      scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
    p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
    p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
    p <- p+annotate("text",x=m1[,5],y=m1[,4], colour="black",size=2,label=m1[,4],angle=-90)+
      annotate("text",x=m1[,5],y=-0.1, colour="black",size=2,label=paste(m1[,2],m1[,3],sep="_"),angle=-90)
    ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/",gsub("_fc.txt","_1.pdf",i),sep=""),width = 10,height=9)
  }
 
  p <- ggplot(m2, aes(x=m2[,5], y=m2[,4], fill=Expression_Value2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m2[,5],y=m2[,4], colour="black",size=2,label=m2[,4],angle=-90)+
    annotate("text",x=m2[,5],y=-0.1, colour="black",size=2,label=paste(m2[,2],m2[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/",gsub("_fc.txt","_2.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m3, aes(x=m3[,5], y=m3[,4], fill=Expression_Value3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m3[,5],y=m3[,4], colour="black",size=2,label=m3[,4],angle=-90)+
    annotate("text",x=m3[,5],y=-0.1, colour="black",size=2,label=paste(m3[,2],m3[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/",gsub("_fc.txt","_3.pdf",i),sep=""),width = 10,height=9)
  
  p <- ggplot(m4, aes(x=m4[,5], y=m4[,4], fill=Expression_Value4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
    scale_fill_manual(values=c("red","grey","purple","lightgreen"))    
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
  p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste(gsub("_fc.txt","",i),sum(d33[,4]),(sum(d22[,4])+ sum(d33[,4])),(sum(d33[,4])/(sum(d22[,4])+ sum(d33[,4]))),sep="_" ))  
  p <- p+annotate("text",x=m4[,5],y=m4[,4], colour="black",size=2,label=m4[,4],angle=-90)+
    annotate("text",x=m4[,5],y=-0.1, colour="black",size=2,label=paste(m4[,2],m4[,3],sep="_"),angle=-90)
  ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/",gsub("_fc.txt","_4.pdf",i),sep=""),width = 10,height=9)
}
h1[,1]=as.numeric(as.character(h1[,1]))
h1[,2]=as.character(h1[,2])
h1=h1[order(h1[,1]),]
h11=data.frame(v1="up",
               v2=gsub("_.*","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h1[,1],h1[,2],sep="_")))),
               v4=as.numeric(table(paste(h1[,1],h1[,2],sep="_"))))
h2[,1]=as.numeric(as.character(h2[,1]))
h2[,2]=as.character(h2[,2])
h2=h2[order(h2[,1]),]
h22=data.frame(v1="normalup",
               v2=gsub("_.*","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h2[,1],h2[,2],sep="_")))),
               v4=as.numeric(table(paste(h2[,1],h2[,2],sep="_"))))
h3[,1]=as.numeric(as.character(h3[,1]))
h3[,2]=as.character(h3[,2])
h3=h3[order(h3[,1]),]
h33=data.frame(v1="normaldown",
               v2=gsub("_.*","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h3[,1],h3[,2],sep="_")))),
               v4=as.numeric(table(paste(h3[,1],h3[,2],sep="_"))))
h4[,1]=as.numeric(as.character(h4[,1]))
h4=h4[order(h4[,1]),]
h44=data.frame(v1="down",
               v2=gsub("_.*","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v3=gsub(".*_","",names(table(paste(h4[,1],h4[,2],sep="_")))),
               v4=as.numeric(table(paste(h4[,1],h4[,2],sep="_"))))
h=na.omit(rbind(h11,h22,h33,h44))
h[,2]=as.numeric(as.character(h[,2]))
h[,3]=gsub("\\+","",as.character(h[,3]))
h=h[order(h[,2],h[,3]),]
h[,5]=h[,2]
for(m in 1:(nrow(h)-1)){
  if(h[m,2]==h[m+1,2]){
    h[m+1,5]=h[m,5]+1
  }
  if((h[m,2]!=h[m+1,2]) & (h[m+1,5]<=h[m,5])){
    h[m+1,5]=h[m,5]+1
  }
}
ExpressionValue=factor(h[,1])
p <- ggplot(h, aes(x=h[,5], y=h[,4], fill=ExpressionValue)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/Allcancer.pdf",sep=""),width = 10,height=9)

n1=h[which(h[,2]<100),]
n2=h[which(h[,2]>=100 & h[,2]< 135),]
n3=h[which(h[,2]>=135 & h[,2]< 160),]
n4=h[which(h[,2]>=160 & h[,2]< 190),]
n5=h[which(h[,2]>=190 & h[,2]< 220),]
n6=h[which(h[,2]>=220 & h[,2]< 250),]
n7=h[which(h[,2]>=250 & h[,2]< 275),]
n8=h[which(h[,2]>=275 & h[,2]< 300),]
n9=h[which(h[,2]>=300),]
ExpressionValue1=factor(n1[,1]) 
ExpressionValue2=factor(n2[,1]) 
ExpressionValue3=factor(n3[,1]) 
ExpressionValue4=factor(n4[,1]) 
ExpressionValue5=factor(n5[,1]) 
ExpressionValue6=factor(n6[,1]) 
ExpressionValue7=factor(n7[,1])
ExpressionValue8=factor(n8[,1])
ExpressionValue9=factor(n9[,1])
p <- ggplot(n1, aes(x=n1[,5], y=n1[,4], fill=ExpressionValue1)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n1[,5],y=n1[,4], colour="black",size=2,label=n1[,4],angle=-90)+
  annotate("text",x=n1[,5],y=-0.1, colour="black",size=2,label=paste(n1[,2],n1[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_1.pdf",sep=""),width = 10,height=9)
p <- ggplot(n2, aes(x=n2[,5], y=n2[,4], fill=ExpressionValue2)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n2[,5],y=n2[,4], colour="black",size=2,label=n2[,4],angle=-90)+
  annotate("text",x=n2[,5],y=-0.1, colour="black",size=2,label=paste(n2[,2],n2[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_2.pdf",sep=""),width = 10,height=9)
p <- ggplot(n3, aes(x=n3[,5], y=n3[,4], fill=ExpressionValue3)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n3[,5],y=n3[,4], colour="black",size=2,label=n3[,4],angle=-90)+
  annotate("text",x=n3[,5],y=-0.1, colour="black",size=2,label=paste(n3[,2],n3[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_3.pdf",sep=""),width = 10,height=9)
p <- ggplot(n4, aes(x=n4[,5], y=n4[,4], fill=ExpressionValue4)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n4[,5],y=n4[,4], colour="black",size=2,label=n4[,4],angle=-90)+
  annotate("text",x=n4[,5],y=-0.1, colour="black",size=2,label=paste(n4[,2],n4[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_4.pdf",sep=""),width = 10,height=9)
p <- ggplot(n5, aes(x=n5[,5], y=n5[,4], fill=ExpressionValue5)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n5[,5],y=n5[,4], colour="black",size=2,label=n5[,4],angle=-90)+
  annotate("text",x=n5[,5],y=-0.1, colour="black",size=2,label=paste(n5[,2],n5[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_5.pdf",sep=""),width = 10,height=9)
p <- ggplot(n6, aes(x=n6[,5], y=n6[,4], fill=ExpressionValue6)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n6[,5],y=n6[,4], colour="black",size=2,label=n6[,4],angle=-90)+
  annotate("text",x=n6[,5],y=-0.1, colour="black",size=2,label=paste(n6[,2],n6[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_6.pdf",sep=""),width = 10,height=9)
p <- ggplot(n7, aes(x=n7[,5], y=n7[,4], fill=ExpressionValue7)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n7[,5],y=n7[,4], colour="black",size=2,label=n7[,4],angle=-90)+
  annotate("text",x=n7[,5],y=-0.1, colour="black",size=2,label=paste(n7[,2],n7[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_7.pdf",sep=""),width = 10,height=9)
p <- ggplot(n8, aes(x=n8[,5], y=n8[,4], fill=ExpressionValue8)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n8[,5],y=n8[,4], colour="black",size=2,label=n8[,4],angle=-90)+
  annotate("text",x=n8[,5],y=-0.1, colour="black",size=2,label=paste(n8[,2],n8[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_8.pdf",sep=""),width = 10,height=9)
p <- ggplot(n9, aes(x=n9[,5], y=n9[,4], fill=ExpressionValue9)) +geom_bar(width=0.1,stat="identity", position="dodge")+
  scale_fill_manual(values=c("red","grey","purple","lightgreen"))     
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Mutation Change Number") + ggtitle(paste("Allcancer",sum(h33[,4]),(sum(h22[,4])+ sum(h33[,4])),(sum(h33[,4])/(sum(h22[,4])+ sum(h33[,4]))),sep="_" ))   
p <- p+annotate("text",x=n9[,5],y=n9[,4], colour="black",size=2,label=n9[,4],angle=-90)+
  annotate("text",x=n9[,5],y=-0.1, colour="black",size=2,label=paste(n9[,2],n9[,3],sep="_"),angle=-90)
ggsave(file=paste(path,"Position_SampleNumber/MutationChange/fc_WT_normalupdown/1/Allcancer_9.pdf",sep=""),width = 10,height=9)




#########################################################
##网络图构建，构建突变、没突变网络。p53_igraph/
##对每个cancer22，点为每个基因，大小为表达值是取在突变、没突变样本中的fc的均值，颜色为对应的功能，
##正方形为上调，圆形为下调，边为2个基因间的关系，<->为促进，权重为同上调1同下调2的个数，
##->为抑制，权重为异调3、4的个数，-为NA或修改（不确定的），权重=1。
#PCPG14只有一个突变样本删去
#install.packages("igraph")
library(igraph) 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
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
  f1$Muexpressmean=apply(f1[,2:(min(grep("_no",colnames(f1)))-1)],1,mean)
  f1$NoMuexpressmean=apply(f1[,min(grep("_no",colnames(f1))):max(grep("_no",colnames(f1)))],1,mean)
  
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name","Muexpressmean","NoMuexpressmean")]
  m=merge(f11,f)
  fu_colour=data.frame(Function=unique(m[,4]),colour=rainbow(length(unique(m[,4]))))
  m1=m[,c("gene_Name","Function","Muexpressmean")]
  m1[which(m1[,3]>=0),4]="square"
  m1[which(m1[,3]<0),4]="circle"
  colnames(m1)[4]="shape"
  m1=merge(m1,fu_colour)
  m1=m1[,c("gene_Name","Function","Muexpressmean","shape","colour")]
  
  m2=m[,c("gene_Name","Function","NoMuexpressmean")]
  m2[which(m2[,3]>=0),4]="square"
  m2[which(m2[,3]<0),4]="circle"
  colnames(m2)[4]="shape"
  m2=merge(m2,fu_colour)
  m2=m2[,c("gene_Name","Function","NoMuexpressmean","shape","colour")]
  
  f2= read.csv(paste("F:/tp53/3generegulate/new_relationship/",gsub("_fc.txt",".txt",i),sep=""),header=T,sep="\t")
  f2[,2]=as.character((f2[,2]))
  f2[,3]=as.character((f2[,3]))
  f2[,4]=as.character((f2[,4]))
  f2=unique(f2)
  n1=f2[,-grep("_no",colnames(f2))]
  n2=f2[,c(1:4,grep("_no",colnames(f2)))]
  
  for(j in 1:nrow(n1)){
    if(n1[j,4]=="-"){
      a=n1[j,2] 
      n1[j,2]=n1[j,3]
      n1[j,3]=a
      n1[j,4]=paste(n1[j,4],0,sep="")
      n1[j,"Mucount"]=length(n1[j,][,which(n1[j,]==3|n1[j,]==4)])
      n1[j,"linetype"]="->"
    }
    else if(n1[j,4]=="+"){
      n1[j,"Mucount"]=length(n1[j,][,which(n1[j,]==1|n1[j,]==2)])
      n1[j,"linetype"]="<->"
    }
    else{
      n1[j,"Mucount"]=1
      n1[j,"linetype"]="-"
    }
  }
  n11=n1[,c("gene1","gene2","function.","Mucount","linetype")]
  n11=merge(n11,fu_colour,by.x="function.",by.y="Function")
  n11=n11[,c("gene1","gene2","function.","Mucount","linetype","colour")]
  
  for(j in 1:nrow(n2)){
    if(n2[j,4]=="-"){
      a=n2[j,2] 
      n2[j,2]=n2[j,3]
      n2[j,3]=a
      n2[j,4]=paste(n2[j,4],0,sep="")
      n2[j,"NoMucount"]=length(n2[j,][,which(n2[j,]==3|n2[j,]==4)])
      n2[j,"linetype"]="->"
    }
    else if(n2[j,4]=="+"){
      n2[j,"NoMucount"]=length(n2[j,][,which(n2[j,]==1|n2[j,]==2)])
      n2[j,"linetype"]="<->"
    }
    else{
      n2[j,"Mucount"]=1
      n2[j,"linetype"]="-"
    }
  }
  n22=n2[,c("gene1","gene2","function.","NoMucount","linetype")]
  n22=merge(n22,fu_colour,by.x="function.",by.y="Function")
  n22=n22[,c("gene1","gene2","function.","NoMucount","linetype","colour")]
  
  data1=n11[(n11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))) & (n11[,2] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))),]
  data2=m1[m1[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1])),]
  g1 <- graph_from_data_frame(data1, directed=TRUE, vertices=data2)
  pdf(paste(path,"p53_igraph/",gsub("_fc.txt","",i),"_Mu.pdf",sep=""))  
  plot(g1, layout = layout.fruchterman.reingold,vertex.label=V(g1)$gene_Name,vertex.label.cex=0.3,
       vertex.size=abs(V(g1)$Muexpressmean),vertex.label.dist=-0.5,vertex.label.color=V(g1)$colour,vertex.color=V(g1)$colour,
       vertex.shape=V(g1)$shape,edge.arrow.mode = E(g1)$linetype,
       edge.width=E(g1)$Mucount/100,edge.color=E(g1)$colour,edge.arrow.size=0.1,main=paste(gsub("_fc.txt","",i),"_Mu",sep=""))
  legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
  dev.off()
  
  data3=n22[(n22[,1] %in% intersect(unique(c(n22[,1],n22[,2])),unique(m2[,1]))) & (n22[,2] %in% intersect(unique(c(n22[,1],n22[,2])),unique(m2[,1]))),]
  data4=m2[m2[,1] %in% intersect(unique(c(n22[,1],n22[,2])),unique(m2[,1])),]
  g2 <- graph_from_data_frame(data3, directed=TRUE, vertices=data4)
  pdf(paste(path,"p53_igraph/",gsub("_fc.txt","",i),"_NoMu.pdf",sep=""))  
  plot(g2, layout = layout.fruchterman.reingold,vertex.label=V(g2)$gene_Name,vertex.label.cex=0.3,
       vertex.size=abs(V(g2)$NoMuexpressmean),vertex.label.dist=-0.5,vertex.label.color=V(g2)$colour,vertex.color=V(g2)$colour,
       vertex.shape=V(g2)$shape,edge.arrow.mode = E(g2)$linetype,
       edge.width=E(g2)$NoMucount/100,edge.color=E(g2)$colour,edge.arrow.size=0.1,main=paste(gsub("_fc.txt","",i),"_NoMu",sep=""))
  legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
  dev.off()
}


##网络图构建，构建突变/没突变网络p53%WTmean。p53_igraph/
##对每个cancer22，点为每个基因，大小为表达值是取在突变/没突变样本中的fc的均值，颜色为对应的功能，
##正方形为上调，圆形为下调，边为2个基因间的关系，<->为促进，权重为同上调1同下调2的个数，
##->为抑制，权重为异调3、4的个数，-为NA或修改（不确定的），权重=1。
#PCPG17只有一个突变样本删去 #vertex.frame.color
#install.packages("igraph")
library(igraph) 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",sep=""))
for(i in filename[c(1:9,12:14,16,18:22,24:26)]){
  f1<-read.csv(paste(path,"p53_context_fpkm_foldchange_p53%WTmean/",i,sep=""),header=F,sep="\t",stringsAsFactors = F)
  colnames(f1)<-strsplit(as.character(f1[1,1])," ")[[1]]
  f1[,2]<-as.numeric(f1[,2])
  f1[sapply(f1,is.infinite)]<-NA
  f1[sapply(f1,is.na)]<-0
  f1<-f1[-1,]
  rownames(f1)<-f1$gene_Name
  f1$Mu_NoMuexpressmean=apply(f1[,2:ncol(f1)],1,mean)
  
  f11=f1[intersect(as.character(f[,1]),as.character(f1$gene_Name)),c("gene_Name_p53/WT_mean","Mu_NoMuexpressmean")]
  m=merge(f11,f,by.x="gene_Name_p53/WT_mean",by.y="gene_Name")
  fu_colour=data.frame(Function=unique(m[,3]),colour=rainbow(length(unique(m[,3]))))
  m1=m[,c("gene_Name_p53/WT_mean","Function","Mu_NoMuexpressmean")]
  m1[which(m1[,3]>=0),4]="square"
  m1[which(m1[,3]<0),4]="circle"
  colnames(m1)[4]="shape"
  m1=merge(m1,fu_colour)
  m1=m1[,c("gene_Name_p53/WT_mean","Function","Mu_NoMuexpressmean","shape","colour")]
  
  f2= read.csv(paste("F:/tp53/3generegulate/new_relationship1/",gsub("_fc.txt",".txt",i),sep=""),header=T,sep="\t")
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
      n1[j,"Mu_NoMucount"]=length(n1[j,][,which(n1[j,]==3|n1[j,]==4)])
      n1[j,"linetype"]="->"
    }
    else if(n1[j,4]=="+"){
      n1[j,"Mu_NoMucount"]=length(n1[j,][,which(n1[j,]==1|n1[j,]==2)])
      n1[j,"linetype"]="<->"
    }
    else{
      n1[j,"Mu_NoMucount"]=1
      n1[j,"linetype"]="-"
    }
  }
  n11=n1[,c("gene1","gene2","function.","Mu_NoMucount","linetype")]
  n11=merge(n11,fu_colour,by.x="function.",by.y="Function")
  n11=n11[,c("gene1","gene2","function.","Mu_NoMucount","linetype","colour")]
  
  data1=n11[(n11[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))) & (n11[,2] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1]))),]
  data2=m1[m1[,1] %in% intersect(unique(c(n11[,1],n11[,2])),unique(m1[,1])),]
  g1 <- graph_from_data_frame(data1, directed=TRUE, vertices=data2)
  pdf(paste(path,"p53_igraph1/",gsub("_fc.txt","",i),"_Mu_NoMu1.pdf",sep=""))  
  plot(g1, layout =layout.circle ,vertex.label=V(g1)$gene_Name,vertex.label.cex=0.3,
       vertex.size=abs(V(g1)$Mu_NoMuexpressmean),vertex.label.dist=-0.5,vertex.label.color=V(g1)$colour,vertex.color=V(g1)$colour,
       vertex.shape=V(g1)$shape,edge.arrow.mode = E(g1)$linetype,vertex.frame.color="gray",
       edge.width=E(g1)$Mu_NoMucount/100,edge.color=E(g1)$colour,edge.arrow.size=0.1,main=paste(gsub("_fc.txt","",i),"_Mu_NoMu",sep=""))
  legend('topright',cex=0.3,text.width=0.2,as.character(fu_colour[,1]),fill = as.character(fu_colour[,2]))
  dev.off()
}

#vertex.label.degree=pi/2,
#layout=layout.circle\layout.fruchterman.reingold\layout.kamada.kawai 
#edge.arrow.mode

# all vertex shapes, minus "raster", that might not be available
shapes <- setdiff(shapes(), "")
g1 <- make_ring(length(shapes))
set.seed(42)
plot(g1, vertex.shape=shapes, vertex.label=shapes, vertex.label.dist=1,
     vertex.size=15, vertex.size2=15,
     vertex.pie=lapply(shapes, function(x) if (x=="pie") 2:6 else 0),
     vertex.pie.color=list(heat.colors(5)))
####
tkplot(g1, layout = layout.circle,vertex.label=V(g1)$gene_Name,vertex.label.cex=0.3,
       vertex.size=abs(V(g1)$Muexpressmean),vertex.label.dist=-0.5,vertex.label.color=V(g1)$colour,vertex.color=V(g1)$colour,
       vertex.shape=V(g1)$shape,
       edge.width=E(g1)$Mucount/100,edge.color=E(g1)$colour,edge.arrow.size=0.1,main=paste(gsub("_fc.txt","",i),"_Mu",sep=""))

id=tkplot(g, layout = layout.fruchterman.reingold,vertex.label=V(g)$gene_Name,vertex.label.cex=0.5,
          vertex.size=V(g)$TCGA.FD.A3B8,vertex.label.dist=-0.5,vertex.label.color=V(g)$colour,vertex.color=V(g)$colour,
          edge.width=E(g)$TCGA.FD.A3B8,edge.color="gray",edge.arrow.size=0.12)
coords <- tkplot.getcoords(id)

 

install.packages("networkD3")
#install.packages("network") 
library("networkD3")
simpleNetwork(n[,1:2])

 # Load data
data(MisLinks)
data(MisNodes)
# Plot
forceNetwork(Links = MisLinks, Nodes = MisNodes,            
             Source = "source", Target = "target",            
             Value = "value", NodeID = "name",            
             Group = "group", opacity = 0.8)
html <-forceNetwork(Links = data1,#线性质数据框  
             Nodes = data2,#节点性质数据框  
             Source = "gene1",#连线的源变量  
             Target = "gene2",#连线的目标变量  
             Value = "TCGA.FD.A3B8",#连线的粗细值  
             NodeID = "gene_Name",#节点名称  
             Group = "Function",#节点的分组  
             Nodesize = "TCGA.FD.A3B8", opacity = 0.8, zoom = TRUE)#节点大小，节点数据框中  
saveNetwork(html,"D:/tp53/networkD3.html",selfcontained=TRUE)

pdf("D:/tp53/1.pdf")
forceNetwork(Links = n,#线性质数据框  
             Nodes = m,#节点性质数据框  
             Source = "gene1",#连线的源变量  
             Target = "gene2",#连线的目标变量  
             Value = "TCGA.FD.A3B8",#连线的粗细值  
             NodeID = "gene_Name",#节点名称  
             Group = "Function",#节点的分组  
             Nodesize = "TCGA.FD.A3B8" ,#节点大小，节点数据框中  
             ###美化部分  
             fontFamily="宋体",#字体设置如"华文行楷" 等  
             fontSize = 20, #节点文本标签的数字字体大小（以像素为单位）。  
             linkColour="black",#连线颜色,black,red,blue,    
             #colourScale ,linkWidth,#节点颜色,red，蓝色blue,cyan,yellow等  
             charge = -100,#数值表示节点排斥强度（负值）或吸引力（正值）    
             opacity = 0.9,  
             legend=T,#显示节点分组的颜色标签  
             arrows=T,#是否带方向  
             bounded=F,#是否启用限制图像的边框  
             #opacityNoHover=1.0,#当鼠标悬停在其上时，节点标签文本的不透明度比例的数值  
             zoom = T)#允许放缩，双击放大  
dev.off()
 
 
####################################################
##突变热点对生存的影响,分Allcancer和cancer在multibase中
library(openxlsx) ##在R3.4.2里安装 打开xlsx文件
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f2=data.frame(cbind(class=paste(f1[,3],f1[,4],sep="_"),as.character(f1$project_id),f1$vital_status,f1$days_to_last_follow_up,f1$days_to_death.y))
f2[,1]=gsub("\\+","",as.character(f2[,1]))
f2[,2]=as.character(f2[,2])
f2[,3]=as.numeric(as.character(f2[,3]))
f2[,4]=as.numeric(as.character(f2[,4]))
f2[,5]=as.numeric(as.character(f2[,5]))
colnames(f2)=c("class","project_id","vital_status","days_to_last_follow_up","days_to_death.y")
b1=subset(f2,f2$vital_status==1)
b2=subset(f2,f2$vital_status==2)
alive=cbind(b1$class,b1$project_id,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$class,b2$project_id,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=gsub("\\+","",as.character(e[,1]))
e[,2]=as.character(e[,2])
e[,3]=as.numeric(as.character(e[,3]))
e[,4]=as.numeric(as.character(e[,4])) 

cancer23=c("Allcancer","BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
           "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
for(i in 1:length(cancer23)){
  ff=read.xlsx(paste(path,"Position_SampleNumber/Muhotspot_survival/1.xlsx",sep=""),sheet=i) 
  if(length(na.omit(ff[,1]))!=0){
    dir.create(paste(path,"Position_SampleNumber/Muhotspot_survival/",cancer23[i],sep="")) 
    if(i==1){
      for(j in na.omit(ff[,1])){
        e[grep(paste("^",j,sep=""),e[,1]),5]=j
        e[-grep(paste("^",j,sep=""),e[,1]),5]=paste(j,"_no",sep="")
        pdf(paste(path,"Position_SampleNumber/Muhotspot_survival/",cancer23[i],"/",gsub(">","",j),"_survival.pdf",sep="")) 
        dif <- survdiff(Surv(as.numeric(e[,4]),as.numeric(e[,3]))~e[,5])#求生存时间
        kmsurvival<-survfit(Surv(as.numeric(e[,4]),as.numeric(e[,3]))~e[,5],conf.type = "log-log")
        plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(j,"_multibase",sep=""))
        legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
               lty='solid',col=c("red","green"))
        text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
        for(j in 1:length(unique(e[,5])) ){
          text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
        }
        dev.off()
      }
    }
    if(i!=1){
      for(j in na.omit(ff[,1])){
        e1=e[e[,2]==cancer23[i],]
        if(nrow(e1)!=0){
          e1[grep(paste("^",j,sep=""),e1[,1]),5]=j
          e1[-grep(paste("^",j,sep=""),e1[,1]),5]=paste(j,"_no",sep="")
        }
        if(nrow(e1)!=0 & length(unique(e1[,5]))!=1){
          pdf(paste(path,"Position_SampleNumber/Muhotspot_survival/",cancer23[i],"/",gsub(">","",j),"_survival.pdf",sep="")) 
          dif <- survdiff(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5])#求生存时间
          kmsurvival<-survfit(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5],conf.type = "log-log")
          plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(j,"_multibase",sep=""))
          legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
                 lty='solid',col=c("red","green"))
          text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
          for(j in 1:length(unique(e1[,5])) ){
            text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
          }
          dev.off()
        }
      }
    }
  }
}


##p值显著在文件夹significant下
for(i in 1:length(cancer23)){
  ff=read.xlsx(paste(path,"Position_SampleNumber/Muhotspot_survival/1.xlsx",sep=""),sheet=i)
  if(i==1){
    for(j in na.omit(ff[,1])){
      e[grep(paste("^",j,sep=""),e[,1]),5]=j
      e[-grep(paste("^",j,sep=""),e[,1]),5]=paste(j,"_no",sep="")
      dif <- survdiff(Surv(as.numeric(e[,4]),as.numeric(e[,3]))~e[,5])#求生存时间
      kmsurvival<-survfit(Surv(as.numeric(e[,4]),as.numeric(e[,3]))~e[,5],conf.type = "log-log")
      if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
        #dir.create(paste(path,"Position_SampleNumber/Muhotspot_survival/significant/",cancer23[i],sep="")) 
        pdf(paste(path,"Position_SampleNumber/Muhotspot_survival/significant/",cancer23[i],"/",gsub(">","",j),"_survival.pdf",sep="")) 
        plot(kmsurvival, lty = 'solid',lwd=2,col=c("indianred3","ivory4"),xlab='Survival time(days)',ylab='Survival(%)',main=paste(j,"_multibase",sep=""))
        legend('bottomleft',cex=0.5,lwd=2,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
               lty='solid',col=c("indianred3","ivory4"))
        text(600,1,cex=0.8,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
        for(j in 1:length(unique(e[,5])) ){
          text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
        }
        dev.off()
      }
    }
  }
  if(i!=1 & length(na.omit(ff[,1]))!=0){
    for(j in na.omit(ff[,1])){
      e1=e[e[,2]==cancer23[i],]
      if(nrow(e1)!=0){
        e1[grep(paste("^",j,sep=""),e1[,1]),5]=j
        e1[-grep(paste("^",j,sep=""),e1[,1]),5]=paste(j,"_no",sep="")
      }
      dif <- survdiff(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5])#求生存时间
      kmsurvival<-survfit(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5],conf.type = "log-log")
      if(nrow(e1)!=0 & length(unique(e1[,5]))!=1 & pchisq(dif$chisq,1,lower.tail=F)<0.05){
        #dir.create(paste(path,"Position_SampleNumber/Muhotspot_survival/significant/",cancer23[i],sep="")) 
        pdf(paste(path,"Position_SampleNumber/Muhotspot_survival/significant/",cancer23[i],"/",gsub(">","",j),"_survival.pdf",sep="")) 
        plot(kmsurvival, lty = 'solid',lwd=2,col=c("indianred3","ivory4"),xlab='Survival time(days)',ylab='Survival(%)',main=paste(j,"_multibase",sep=""))
        legend('bottomleft',cex=0.5,lwd=2,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
               lty='solid',col=c("indianred3","ivory4"))
        text(600,1,cex=0.8,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
        for(j in 1:length(unique(e1[,5])) ){
          text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
        }
        dev.off()
      }
    }
  }
}

#第一种方法：读取csv
a<-read.csv("exercise1.csv",header = T)
#第二种方法：RODBC包
#安装载入RODBC包，如果已安装，请跳过第一句语句
install.packages(RODBC)
library(RODBC)
ab<-odbcConnectExcel2007("exercise1.xls")#连接excel，32位系统使用odbcConnectExcel函数
sqlTables(ab)
a<-sqlFetch(ab,"Sheet1$")
odbcClose(ab)#关闭句柄，此句是必须。
#第三种方法：openxlsx
install.packages(openxlsx)
library(openxlsx)
a<-read.xlsx("exercise1.xlsx",sheet=1)#文件名+sheet的序号，简单粗暴

 
####################################################
##突变类型对生存的影响，不考虑dup和del类型的突变，只有单个碱基突变，如A>T，c>G,整体为Mutype_survival
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,2]=gsub("\\+","",as.character(f1[,2]))
f1[,3]=as.numeric(f1[,3])
f1[,4]=as.numeric(f1[,4])
f1[,5]=as.numeric(f1[,5])
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=gsub("\\+","",as.character(e[,1]))
e[,2]=as.numeric(e[,2])
e[,3]=as.numeric(e[,3])
e1=e[nchar(e[,1])==3,]

pdf(paste(path,"Position_SampleNumber/Mutype_survival/Mutype_survival.pdf",sep="")) 
dif <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
kmsurvival<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
co=c("red","orange4","maroon4","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue1")
plot(kmsurvival, lty = 'solid', col=co,xlab='survival time in days',ylab='survival probabilities',main="Mutype_survival")
legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
       lty='solid',col=co)
text(1000,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
for(j in 1:length(unique(e1[,1])) ){
  text(300,0.9-(j-1)*0.03,cex=0.8,dif$n[j])
}
dev.off()

##1、看具体突变类型对生存的影响，只有单个碱基突变，如A>T，和非A>T样本的生存分析。
##AT_survival为AT和非AT在singlebase
##AT_survival1为AT和非AT在multibase(包括dup和del)
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,2]=gsub("\\+","",as.character(f1[,2]))
f1[,3]=as.numeric(f1[,3])
f1[,4]=as.numeric(f1[,4])
f1[,5]=as.numeric(f1[,5])

b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=gsub("\\+","",as.character(e[,1]))
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.numeric(as.character(e[,3]))
e1=e[nchar(e[,1])==3,]
##singlebase
for(i in unique(e1[,1])){
  e1[e1[,1]==i,4]=i
  e1[e1[,1]!=i,4]=paste(i,"_no",sep="")
  pdf(paste(path,"Position_SampleNumber/Mutype_survival/",gsub(">","",i),"_survival.pdf",sep="")) 
  dif <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,4],conf.type = "log-log")
  plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,"_singlebase",sep=""))
  legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
         lty='solid',col=c("red","green"))
  text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
  for(j in 1:length(unique(e1[,4])) ){
    text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
  }
  dev.off()
}
##multibase
for(i in unique(e1[,1])){
  e[e[,1]==i,4]=i
  e[e[,1]!=i,4]=paste(i,"_no",sep="")
  pdf(paste(path,"Position_SampleNumber/Mutype_survival/",gsub(">","",i),"_survival1.pdf",sep="")) 
  dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
  plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,"_multibase",sep=""))
  legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
         lty='solid',col=c("red","green"))
  text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
  for(j in 1:length(unique(e[,4])) ){
    text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
  }
  dev.off()
}

##2、看具体突变类型对生存的影响，只有dup和del，如A>T，和非A>T样本的生存分析。
##dup_survival1为dup和非dup在multibase(包括dup和del)
##del_survival1为dup和非del在multibase(包括dup和del)
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,2]=gsub("\\+","",as.character(f1[,2]))
f1[,3]=as.numeric(f1[,3])
f1[,4]=as.numeric(f1[,4])
f1[,5]=as.numeric(f1[,5])

b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=gsub("\\+","",as.character(e[,1]))
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.numeric(as.character(e[,3]))
e1=e[nchar(e[,1])>3,]
for(i in unique(e1[,1])){
  e[e[,1]==i,4]=i
  e[e[,1]!=i,4]=paste(i,"_no",sep="")
  pdf(paste(path,"Position_SampleNumber/Mutype_survival/multibase/",gsub(">","",i),"_survival1.pdf",sep="")) 
  dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
  plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,"_multibase",sep=""))
  legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
         lty='solid',col=c("red","green"))
  text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
  for(j in 1:length(unique(e[,4])) ){
    text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
  }
  dev.off()
}
##p值显著在文件夹multibase_significant下
for(i in unique(e1[,1])){
  e[e[,1]==i,4]=i
  e[e[,1]!=i,4]=paste(i,"_no",sep="")
  dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
  if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
    pdf(paste(path,"Position_SampleNumber/Mutype_survival/multibase_significant/",gsub(">","",i),"_survival1.pdf",sep="")) 
    plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,"_multibase",sep=""))
    legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
           lty='solid',col=c("red","green"))
    text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
    for(j in 1:length(unique(e[,4])) ){
      text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
    }
    dev.off()
  }
}

##分为dup、del、ins和no四类做生存分析
e[grep("^dup",e[,1]),4]="dup"
e[grep("^del",e[,1]),4]="del"
e[grep("^ins",e[,1]),4]="ins"
e[nchar(e[,1])==3,4]="no"
dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
pdf(paste(path,"Position_SampleNumber/Mutype_survival/4_survival1.pdf",sep="")) 
plot(kmsurvival, lty = 'solid', col=c("red","green","purple","blue"),xlab='survival time in days',ylab='survival probabilities',main="4_multibase")
legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
       lty='solid',col=c("red","green","purple","blue"))
text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
for(j in 1:length(unique(e[,4])) ){
  text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
}
dev.off()

##3、看具体突变类型Variant_Classification对生存的影响，如Missense_Mutation，和非Missense_Mutation样本的生存分析
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f$Variant_Classification),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","Variant_Classification","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,2]=gsub("^3UTR.*","3UTR",as.character(f1[,2]))
f1[,2]=gsub("^5UTR.*","5UTR",as.character(f1[,2]))
f1[,3]=as.numeric(f1[,3])
f1[,4]=as.numeric(f1[,4])
f1[,5]=as.numeric(f1[,5])
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(b1$Variant_Classification,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$Variant_Classification,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.numeric(as.character(e[,3]))
for(i in unique(e[,1])){
  e[e[,1]==i,4]=i
  e[e[,1]!=i,4]=paste(i,"_no",sep="")
  dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
  if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
    pdf(paste(path,"Position_SampleNumber/Mutype_survival/Variant_Classification/",strtrim(gsub("_","",i),10),".pdf",sep="")) 
    plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=i)
    legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
           lty='solid',col=c("red","green"))
    text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
    for(j in 1:length(unique(e[,4])) ){
      text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
    }
    dev.off()
  }
}










#####################################################
##突变区域细分，根据每个滑窗的宽度和步长，对domain分段，找出对应分段的样本和非样本的logp值<0.05，并画生存图。
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,3]=as.numeric(as.character(f1[,3]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
f1[,7]=as.numeric(as.character(f1[,7]))
f2=f1[,c("submitter_id","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")]
f2=na.omit(f2) 
b1=subset(f2,f2$vital_status==1)
b2=subset(f2,f2$vital_status==2)
alive=cbind(b1$Muposition,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(b2$Muposition,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)

wide=10
step=3
dir.create(paste(path,"Position_SampleNumber/Mudomain_survival/",paste(wide,step,sep="_"),sep="")) 
for(i in seq(0,max(e[,1]),step)){
  e[(e[,1])>=i & (e[,1])<(i+wide),4]=1
  e[!((e[,1])>=i & (e[,1])<(i+wide)),4]=2
  if(length(unique(e[,4]))!=1){
    dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
    kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
    if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
      pdf(paste(path,"Position_SampleNumber/Mudomain_survival/",paste(wide,step,sep="_"),"/",paste(i,(i+wide),sep="_"),".pdf",sep="")) 
      plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,(i+wide),"base",sep="_"))
      legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
             lty='solid',col=c("red","green"))
      text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
      for(j in 1:length(unique(e[,4])) ){
        text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
      }
      dev.off()
    }
  }
}

wide=1
step=0
dir.create(paste(path,"Position_SampleNumber/Mudomain_survival/",paste(wide,step,sep="_"),sep="")) 
for(i in 0:max(e[,1])){
  e[(e[,1])>=i & (e[,1])<(i+wide),4]=1
  e[!((e[,1])>=i & (e[,1])<(i+wide)),4]=2
  if(length(unique(e[,4]))!=1 & nrow(e[which(e[,4]==1),])>=3 ){
    dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
    kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
    if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
      pdf(paste(path,"Position_SampleNumber/Mudomain_survival/",paste(wide,step,sep="_"),"/",paste(i,(i+wide),sep="_"),".pdf",sep="")) 
      plot(kmsurvival, lty = 'solid', col=c("red","green"),xlab='survival time in days',ylab='survival probabilities',main=paste(i,(i+wide),"base",sep="_"))
      legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
             lty='solid',col=c("red","green"))
      text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
      for(j in 1:length(unique(e[,4])) ){
        text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
      }
      dev.off()
    }
  }
}


##根据10_3、10_2、5_3、3_1确定最终分段及1_0 确定重要点，文件为overlap，根据所得最终分段做生存分析
library(ggplot2)
wide=1
step=0
a1=data.frame()
a2=data.frame()
for(i in 0:max(e[,1])){
  e[(e[,1])>=i & (e[,1])<(i+wide),4]=1
  e[!((e[,1])>=i & (e[,1])<(i+wide)),4]=2
  if(length(unique(e[,4]))!=1 & nrow(e[which(e[,4]==1),])>=3 ){
    dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
    kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
    if(pchisq(dif$chisq,1,lower.tail=F)<0.05){
      a1=rbind(a1,i)
      a2=rbind(a2,dif$n[1])
    }
  }
}
a=cbind(a1,a2)
b1=data.frame(cbind(c(12,42,123,174,306,324,357),c(25,64,133,184,316,340,367)))
b2=data.frame(cbind(c(12,42,124,166,260,272,298,322,356),c(26,64,134,184,270,282,316,340,366)))
b3=data.frame(cbind(c(45,123,147,171,261,303,327,360),c(59,128,152,176,266,311,335,365)))
b4=data.frame(cbind(c(46,66,124,142,148,164,173,251,259,296,305,323,329,360),c(57,69,129,145,151,167,176,254,265,301,310,327,334,365)))
p<-ggplot(a, aes(x = a[,1], y = a[,2]))+geom_bar(width=0.1,stat="identity", position="dodge",colour="pink")  
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black"),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Position") + ylab("Sample Change Number") + ggtitle("Position_Number")  
p <- p+annotate("text",x=a[,1],y=a[,2], colour="black",size=2,label=a[,2],angle=0)+
  annotate("text",x=a[,1],y=-0.2, colour="black",size=1.8,label=paste(a[,1],a[,2],sep="_"),angle=-90)+
  annotate("text",x=0,y=-2, colour="black",size=1.8,label="10_3",angle=0)+
  annotate("text",x=0,y=-3, colour="black",size=1.8,label="10_2",angle=0)+
  annotate("text",x=0,y=-4, colour="black",size=1.8,label="5_3",angle=0)+
  annotate("text",x=0,y=-5, colour="black",size=1.8,label="3_1",angle=0)+
  annotate("text",x=unique(c(b1[,1],b1[,2],b2[,1],b2[,2],b3[,1],b3[,2],b4[,1],b4[,2])),y=-7, colour="black",size=1.8,label=unique(c(b1[,1],b1[,2],b2[,1],b2[,2],b3[,1],b3[,2],b4[,1],b4[,2])),angle=-90)
for(j in 1:nrow(b1)){
  p <- p+geom_rect(xmin=b1[j,1],xmax=b1[j,2],ymin=-3,ymax=-2,fill="plum2",alpha=0.2)
}
for(j in 1:nrow(b2)){
  p <- p+geom_rect(xmin=b2[j,1],xmax=b2[j,2],ymin=-4,ymax=-3,fill="plum2",alpha=0.2)
}
for(j in 1:nrow(b3)){
  p <- p+geom_rect(xmin=b3[j,1],xmax=b3[j,2],ymin=-5,ymax=-4,fill="plum2",alpha=0.2)
}
for(j in 1:nrow(b4)){
  p <- p+geom_rect(xmin=b4[j,1],xmax=b4[j,2],ymin=-6,ymax=-5,fill="plum2",alpha=0.2)
}
ggsave(file=paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/overlap.pdf",sep=""),width = 10,height=9)


##对于最终分段之间做生存分析
m=data.frame(cbind(c(12,42,66,123,142,164,171,251,259,296,322,360),c(26,64,69,134,152,167,176,254,266,310,334,365)))
dir.create(paste(path,"Position_SampleNumber/Mudomain_survival/result",sep="")) 
for(i in 1:nrow(m)){
  e[(e[,1])>=m[i,1] & (e[,1])<m[i,2],4]=paste(m[i,1],m[i,2],sep="_")
}
dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
co=c("red","orange4","maroon4","maroon","limegreen","olivedrab1","olivedrab","purple4","purple","plum2","royalblue4","royalblue1")
pdf(paste(path,"Position_SampleNumber/Mudomain_survival/result/part.pdf",sep="")) 
plot(kmsurvival, lty = 'solid', col=co,xlab='survival time in days',ylab='survival probabilities',main="part")
legend('bottomleft',cex=0.5,text.width=0.4,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
       lty='solid',col=co)
text(600,1,cex=0.8,paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""))
for(j in 1:length(unique(e[,4])) ){
  text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
}
dev.off()

 
##对于每一个最终分段，取在这分段和不在这分段里的样本做生存分析 
m=data.frame(cbind(c(12,42,66,123,142,148,164,171,251,259,296,322),c(26,64,69,134,145,151,167,176,254,266,310,334)))
for(i in 1:nrow(m)){
  e[(e[,1])>=m[i,1] & (e[,1])<m[i,2],4]=paste(m[i,1],m[i,2],sep="_")
  e[!((e[,1])>=m[i,1] & (e[,1])<m[i,2]),4]=paste(m[i,1],m[i,2],"no",sep="_")
  dif <- survdiff(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e[,3]),as.numeric(e[,2]))~e[,4],conf.type = "log-log")
  co=c("indianred3","ivory4")
  pdf(paste(path,"Position_SampleNumber/Mudomain_survival/result/",paste(m[i,1],m[i,2],sep="_"),".pdf",sep="")) 
  plot(kmsurvival, lty = 'solid',lwd=2,col=co,xlab='survival time in days',ylab='survival probabilities',main=paste(m[i,1],m[i,2],sep="_"))
  legend('bottomleft',cex=0.6,lwd=2,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
         lty='solid',col=co)
  text(600,1,cex=1,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
  for(j in 1:length(unique(e[,4])) ){
    text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
  }
  dev.off()
}


##对于每一个最终分段，取在这分段和不在这分段里的样本做表型分析,先统计在每段的cancer数目分布 
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,3]=as.numeric(as.character(f1[,3]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
f1[,7]=as.numeric(as.character(f1[,7]))
f2=f1[,c("project_id","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")]
f2=na.omit(f2) 
b1=subset(f2,f2$vital_status==1)
b2=subset(f2,f2$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
m=data.frame(cbind(c(12,42,66,123,142,148,164,171,251,259,296,322),c(26,64,69,134,145,151,167,176,254,266,310,334)))
for(i in 1:nrow(m)){
  a=cbind(paste(m[i,1],m[i,2],sep="_"),data.frame(table(e[(e[,2])>=m[i,1] & (e[,2])<m[i,2],1])))
  write.table(a,file=paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
}

ff= read.table(paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer.txt",sep=""),header=F,sep="\t")
ExpressionValue=factor(ff[,2])
p<-ggplot(ff, aes(x = ff[,1], y = ff[,3], fill=ExpressionValue))+geom_bar(width=0.6,stat="identity", position="dodge")  
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("domain")+ylab("number")+ggtitle("cancer")
p <- p+annotate("text",x=ff[,1],y=ff[,3], colour="black",size=1.8,label=paste(ff[,2],ff[,3],sep=""),angle=0)
ggsave(file=paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer.pdf",sep=""),width = 10,height=9)

 
##对于每一个最终分段，对所有cancer发生突变的数目，在分段的数目，每个cancer突变数目、每个cancer在分段的数目做fisher检验。 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))))))
colnames(f1)=c("project_id","Muposition")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
ff= read.table(paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer.txt",sep=""),header=F,sep="\t")
for(i in 1:nrow(ff)){
  a=nrow(f1[f1[,2]>=as.numeric(strsplit(as.character(ff[i,1]),"_")[[1]][1]) & f1[,2]<as.numeric(strsplit(as.character(ff[i,1]),"_")[[1]][2]),])
  b=nrow(f1[f1[,1]==as.character(ff[i,2]),])
  x<-matrix(c(nrow(f1),a,b,ff[i,3]),nrow=2,byrow = T)
  d=paste(ff[i,1],ff[i,2],fisher.test(x)$p.value,sep="\t")
  #write.table(d,file=paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer_p.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  if(fisher.test(x)$p.value<0.05){
    write.table(d,file=paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer_p_significant.txt",sep=""),row.names = F,col.names =F, quote = F,sep = "\t",append=TRUE)
  }
}


##对于每一分段，对显著的cancer和其他cancer做生存分析
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=as.numeric(as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1=na.omit(f1) 
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$vital_status,b2$days_to_death.y)	
e=rbind(alive,dead)
e=data.frame(e)
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.numeric(as.character(e[,3]))
e[,4]=as.numeric(as.character(e[,4]))
ff= read.table(paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/cancer_p_significant.txt",sep=""),header=F,sep="\t")
for(i in 1:nrow(ff)){
  e1=e[e[,2]>=as.numeric(strsplit(as.character(ff[i,1]),"_")[[1]][1]) & e[,2]<as.numeric(strsplit(as.character(ff[i,1]),"_")[[1]][2]),]
  e1[e1[,1]==as.character(ff[i,2]),5]=paste(as.character(ff[i,1]),as.character(ff[i,2]),sep="_")
  e1[e1[,1]!=as.character(ff[i,2]),5]=paste(as.character(ff[i,1]),as.character(ff[i,2]),"no",sep="_")
  dif <- survdiff(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5])#求生存时间
  kmsurvival<-survfit(Surv(as.numeric(e1[,4]),as.numeric(e1[,3]))~e1[,5],conf.type = "log-log")
  co=c("indianred3","ivory4")
  if(length(unique(e1[,3]))!=1 ){
    pdf(paste(path,"Position_SampleNumber/Mudomain_survival/overlapdomain/",paste(as.character(ff[i,1]),as.character(ff[i,2]),sep="_"),".pdf",sep="")) 
    plot(kmsurvival, lty = 'solid',lwd=2, col=co,xlab='survival time in days',ylab='survival probabilities',main=paste(as.character(ff[i,1]),as.character(ff[i,2]),sep="_"))
    legend('bottomleft',cex=0.5,lwd=2,gsub(".*=","",as.character(data.frame(dif$n)[,1])), 
           lty='solid',col=co)
    text(600,1,cex=0.8,paste("p=",signif(pchisq(dif$chisq,1,lower.tail=F),3),sep=""))
    for(j in 1:length(unique(e[,4])) ){
      text(300,0.6-(j-1)*0.05,cex=0.8,dif$n[j])
    }
    dev.off()
  }
}
  



#####对于每一最终分段，对表型做分析
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))))))
colnames(f1)=c("submitter_id","project_id","Muposition")
f1[,1]=as.character(f1[,1])
f1[,2]=as.character(f1[,2])
f1[,3]=as.numeric(as.character(f1[,3]))
f1=na.omit(f1) 
f1[(f1[,3]>=12 & f1[,3]<=26),4]="D1"
f1[(f1[,3]>=42 & f1[,3]<=64),4]="D2"
f1[(f1[,3]>=66 & f1[,3]<=69),4]="D3"
f1[(f1[,3]>=123 & f1[,3]<=134),4]="D4"
f1[(f1[,3]>=142 & f1[,3]<=145),4]="D5"
f1[(f1[,3]>=148 & f1[,3]<=151),4]="D6"
f1[(f1[,3]>=164 & f1[,3]<=167),4]="D7"
f1[(f1[,3]>=171 & f1[,3]<=176),4]="D8"
f1[(f1[,3]>=251 & f1[,3]<=254),4]="D9"
f1[(f1[,3]>=259 & f1[,3]<=266),4]="D10"
f1[(f1[,3]>=296 & f1[,3]<=310),4]="D11"
f1[(f1[,3]>=322 & f1[,3]<=334),4]="D12"
f2=na.omit(f1)
f2[,5:9]="NA"
colnames(f2)[5:9]=c("gender.demographic","initial_weight.samples","days_to_new_tumor_event_after_initial_treatment",
                    "neoplasm_histologic_grade","route_of_administration")
for(i in 1:nrow(f2)){
  ff= read.csv(paste("D:/tp53/phenotype/TCGA-",f2[i,2],".GDC_phenotype.tsv",sep=""),header=T,sep="\t")
  if(TRUE %in% (colnames(f2)[5] %in% colnames(ff)) ){
    ff1=ff[,c("submitter_id.samples",colnames(f2)[5])]
    ff1=ff1[grep("-01A$",as.character(ff1[,1])),]
    ff1[,1]=substring(as.character(ff1[,1]),1,12)
    ff1=unique(ff1)
    rownames(ff1)=ff1[,1]
    f2[i,5]=as.character(ff1[as.character(f2[i,1]),colnames(f2)[5]])
  }
  if(TRUE %in% (colnames(f2)[6] %in% colnames(ff)) ){
    ff1=ff[,c("submitter_id.samples",colnames(f2)[6])]
    ff1=ff1[grep("-01A$",as.character(ff1[,1])),]
    ff1[,1]=substring(as.character(ff1[,1]),1,12)
    ff1=unique(ff1)
    rownames(ff1)=ff1[,1]
    f2[i,6]=as.character(ff1[as.character(f2[i,1]),colnames(f2)[6]])
  }
  if(TRUE %in% (colnames(f2)[7] %in% colnames(ff)) ){
    ff1=ff[,c("submitter_id.samples",colnames(f2)[7])]
    ff1=ff1[grep("-01A$",as.character(ff1[,1])),]
    ff1[,1]=substring(as.character(ff1[,1]),1,12)
    ff1=unique(ff1)
    rownames(ff1)=ff1[,1]
    f2[i,7]=as.character(ff1[as.character(f2[i,1]),colnames(f2)[7]])
  }
  
  if(TRUE %in% (colnames(f2)[8] %in% colnames(ff)) ){
    ff1=ff[,c("submitter_id.samples",colnames(f2)[8])]
    ff1=ff1[grep("-01A$",as.character(ff1[,1])),]
    ff1[,1]=substring(as.character(ff1[,1]),1,12)
    ff1=unique(ff1)
    rownames(ff1)=ff1[,1]
    f2[i,8]=as.character(ff1[as.character(f2[i,1]),colnames(f2)[8]])
  }
  
  if(TRUE %in% (colnames(f2)[9] %in% colnames(ff)) ){
    ff1=ff[,c("submitter_id.samples",colnames(f2)[9])]
    ff1=ff1[grep("-01A$",as.character(ff1[,1])),]
    ff1[,1]=substring(as.character(ff1[,1]),1,12)
    ff1=unique(ff1)
    rownames(ff1)=ff1[,1]
    f2[i,9]=as.character(ff1[as.character(f2[i,1]),colnames(f2)[9]])
  }
} 
write.table(f2,file=paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno.txt",sep=""),row.names = F,col.names =T, quote = F,sep = "\t",append=TRUE)


####对5种表型进行分析 Table2
write.table(paste("Variables","Domain","in_domain","out_domain","pValue",sep="\t"),
            file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

##离散型c("gender.demographic","neoplasm_histologic_grade","route_of_administration")
#"gender.demographic"
f5=na.omit(f2[,c("V4","gender.demographic")])
f5=f5[(f5[,2]!="NA"),]
for(j in unique(f5[,1])){
  a1=f5[(f5[,1]==j),2]
  a2=f5[(f5[,1]!=j),2]
  if(length(table(a1))>1 & length(table(a2))>1){
    p1=fisher.test(rbind(table(a1),table(a2)))$p.value  
  }else{p1=-5}
  g1=paste("gender.demographic",j,"","",signif(p1,3),sep="\t")
  g2=paste("",names(table(a2)),table(a1),table(a2),"",sep="\t")
  write.table(g1,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  write.table(g2,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
#"neoplasm_histologic_grade"
f6=na.omit(f2[,c("V4","neoplasm_histologic_grade")])
f6=f6[(f6[,2]!="NA"),]
for(j in unique(f6[,1])){
  a1=f6[(f6[,1]==j),2]
  a2=f6[(f6[,1]!=j),2]
  if(length(table(a1))>1 & length(table(a2))>1){
    p1=fisher.test(rbind(table(a1),table(a2)))$p.value  
  }else{p1=-5}
  g1=paste("neoplasm_histologic_grade",j,"","",signif(p1,3),sep="\t")
  g2=paste("",names(table(a2)),table(a1),table(a2),"",sep="\t")
  write.table(g1,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  write.table(g2,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
#"route_of_administration"
f7=na.omit(f2[,c("V4","route_of_administration")])
f7=f7[(f7[,2]!="NA"),]
for(j in unique(f7[,1])){
  a1=f7[(f7[,1]==j),2]
  a2=f7[(f7[,1]!=j),2]
  if(length(table(a1))>1 & length(table(a2))>1){
    p1=fisher.test(rbind(table(a1),table(a2)))$p.value  
  }else{p1=-5}
  g1=paste("route_of_administration",j,"","",signif(p1,3),sep="\t")
  g2=paste("",names(table(a2)),table(a1),table(a2),"",sep="\t")
  write.table(g1,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  write.table(g2,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
##连续型c("initial_weight.samples","days_to_new_tumor_event_after_initial_treatment")
#"initial_weight.samples"
f3=na.omit(f2[,c("V4","initial_weight.samples")])
f3=f3[(f3[,2]!="NA"),]
for(j in unique(f3[,1])){
  a1=as.numeric(f3[(f3[,1]==j),2])
  a2=as.numeric(f3[(f3[,1]!=j),2])
  if(length(a1)>=2 & length(a2)>=2){
    p1=t.test(a1,a2)$p.value  
  }else{p1=-5}
  g1=paste("initial_weight.samples",j,paste(round(mean(a1),2),round(sd(a1),2),sep="+"),paste(round(mean(a2),2),round(sd(a2),2),sep="+"),signif(p1,3),sep="\t")
  write.table(g1,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
#"days_to_new_tumor_event_after_initial_treatment"
f4=na.omit(f2[,c("V4","days_to_new_tumor_event_after_initial_treatment")])
f4=f4[(f4[,2]!="NA"),]
for(j in unique(f4[,1])){
  a1=as.numeric(f4[(f4[,1]==j),2])
  a2=as.numeric(f4[(f4[,1]!=j),2])
  if(length(a1)>=2 & length(a2)>=2){
    p1=t.test(a1,a2)$p.value  
  }else{p1=-5}
  g1=paste("days_to_new_tumor_event_after_initial_treatment",j,paste(round(mean(a1),2),round(sd(a1),2),sep="+"),paste(round(mean(a2),2),round(sd(a2),2),sep="+"),signif(p1,3),sep="\t")
  write.table(g1,file =paste(path,"Position_SampleNumber/Mudomain_survival/result/12domain_pheno_p.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}




###########################################################
##Allcancer对12domain,Missense_Mutation,3’UTR,Splice_site,G>C每类加HR  
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),as.character(f$Variant_Classification),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Mutype","Variant_Classification","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=gsub("\\+","",as.character(f1[,2]))
f1[,3]=as.character(f1[,3])
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
f1[,7]=as.numeric(as.character(f1[,7]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$Variant_Classification,b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$Variant_Classification,b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.character(e[,3])
e[,4]=as.character(e[,4])
e[,5]=as.numeric(as.character(e[,5]))
e[,6]=as.numeric(as.character(e[,6]))
m=data.frame(cbind(c(12,42,66,123,142,164,171,251,259,296,322,360),c(26,64,69,134,152,167,176,254,266,310,334,365)))
e[,7]="-10_0"
for(i in 1:nrow(m)){
  e[(e[,2])>=m[i,1] & (e[,2])<m[i,2],7]=paste(m[i,1],m[i,2],sep="_")
}
write.table(paste("class","in_number","out_number","HR","lower95","high95","p",sep="\t"),
            file =paste(path,"Position_SampleNumber/MuclassHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
for(i in as.character(sort(unique(e[,7])))){ ##12domain in VS out
  e[which(e[,7]==i),8]=i
  e[which(e[,7]!=i),8]=paste(i,"no",sep="_")
  e3=na.omit(e[,c("X5","X6","V8")])
  e3[,1]=as.numeric(as.character(e3[,1]))
  e3[,2]=as.numeric(as.character(e3[,2]))
  e3[,3]=as.character(e3[,3])
  cox1=coxph(Surv(e3[,2],e3[,1])~e3[,3])
  write.table(paste(i,nrow(e3[e3[,3]==i,]),nrow(e3[e3[,3]!=i,]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
              file =paste(path,"Position_SampleNumber/MuclassHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
for(i in c("Missense_Mutation","Splice_Site")){ ##Variant_Classification\3’UTR NO
  e[which(e[,3]==i),8]=i
  e[which(e[,3]!=i),8]=paste(i,"no",sep="_")
  e3=na.omit(e[,c("X5","X6","V8")])
  e3[,1]=as.numeric(as.character(e3[,1]))
  e3[,2]=as.numeric(as.character(e3[,2]))
  e3[,3]=as.character(e3[,3])
  cox1=coxph(Surv(e3[,2],e3[,1])~e3[,3])
  write.table(paste(i,nrow(e3[e3[,3]==i,]),nrow(e3[e3[,3]!=i,]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
              file =paste(path,"Position_SampleNumber/MuclassHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}
for(i in c("G>C")){ ##Mutype
  e[which(e[,4]==i),8]=i
  e[which(e[,4]!=i),8]=paste(i,"no",sep="_")
  e3=na.omit(e[,c("X5","X6","V8")])
  e3[,1]=as.numeric(as.character(e3[,1]))
  e3[,2]=as.numeric(as.character(e3[,2]))
  e3[,3]=as.character(e3[,3])
  cox1=coxph(Surv(e3[,2],e3[,1])~e3[,3])
  write.table(paste(i,nrow(e3[e3[,3]==i,]),nrow(e3[e3[,3]!=i,]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
              file =paste(path,"Position_SampleNumber/MuclassHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}


##Allcancer对每类显著的HR加boxplot 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Position_SampleNumber/"
f=read.table(paste(path,"MuclassHR_signif.txt",sep=""),header=T,sep="\t")
f[,1]=as.character(f[,1])
f[,8]=1:nrow(f)
pdf(file=paste(path,"MuclassHR_signif.pdf",sep=""),width = 10,height=9)
plot(f[,c(4,8)],pch=15,col= "black",cex = 2,las=1,xlim=c(min(f[,5])-2,max(f[,6])+5),ylim=c(min(f[,8])-0.2,max(f[,8])+0.2),
     yaxt="n",main ="MuClass",xlab="Hazard Ratio(95% CI)",ylab="Class")
abline(v=1,lwd=1,col="grey",lty=2)#虚线
text(f[,6]+1,f[,8],cex=0.8,signif(f[,7],3))
text(f[,5]-1.8,f[,8],cex=0.8,paste(f[,1],"=",as.character(f[,2]),sep=""))
text(f[,6]+3.5,f[,8],cex=0.8,paste(round(f[,4],2),"(",round(f[,5],2),",",round(f[,6],2),")",sep=""))
segments(f[,5],f[,8],f[,6],f[,8],col="black",lwd=2)
dev.off() 



##HNSC对domain296_310加HR  
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=as.numeric(as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.numeric(as.character(e[,3]))
e[,4]=as.numeric(as.character(e[,4]))
e1=e[(e[,2])>=296 & (e[,2])<310,]
e1[,5]="296_310_HNSC_no"
e1[e1[,1]=="HNSC",5]="296_310_HNSC"
write.table(paste("class","in_number","out_number","HR","lower95","high95","p",sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
cox1=coxph(Surv(e1[,4],e1[,3])~e1[,5])
write.table(paste("296_310_HNSC",nrow(e1[e1[,5]=="296_310_HNSC",]),nrow(e1[e1[,5]!="296_310_HNSC",]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)

##HNSC对hotspot298_G>T加HR  
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=gsub("\\+","",as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.character(e[,3])
e[,4]=as.numeric(as.character(e[,4]))
e[,5]=as.numeric(as.character(e[,5]))
e[,6]=paste(e[,2],e[,3],sep="_")
e1=e[e[,1]=="HNSC",]
e1[,7]="HNSC_298_G>T_no"
e1[e1[,6]=="298_G>T",7]="HNSC_298_G>T"
cox1=coxph(Surv(e1[,5],e1[,4])~e1[,7])
write.table(paste("HNSC_298_G>T",nrow(e1[e1[,7]=="HNSC_298_G>T",]),nrow(e1[e1[,7]!="HNSC_298_G>T",]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)


##ESCA对hotspot248_C>T\248_G>A加HR  
#LUSC_220_A>G\LUSC_234_A>G
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=gsub("\\+","",as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.character(e[,3])
e[,4]=as.numeric(as.character(e[,4]))
e[,5]=as.numeric(as.character(e[,5]))
e[,6]=paste(e[,2],e[,3],sep="_")
e1=e[e[,1]=="LUSC",]
e1[,7]="LUSC_220_A>G_no"
e1[e1[,6]=="220_A>G",7]="LUSC_220_A>G"
cox1=coxph(Surv(e1[,5],e1[,4])~e1[,7])
write.table(paste("LUSC_220_A>G",nrow(e1[e1[,7]=="LUSC_220_A>G",]),nrow(e1[e1[,7]!="LUSC_220_A>G",]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
e1[,7]="LUSC_234_A>G_no"
e1[e1[,6]=="234_A>G",7]="LUSC_234_A>G"
cox2=coxph(Surv(e1[,5],e1[,4])~e1[,7])
write.table(paste("LUSC_234_A>G",nrow(e1[e1[,7]=="LUSC_234_A>G",]),nrow(e1[e1[,7]!="LUSC_234_A>G",]),summary(cox2)$coefficients[,2],summary(cox2)$conf.int[,3],summary(cox2)$conf.int[,4],summary(cox2)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



##LUAD对hotspot158_G>T加HR
#BRCA_220_A>G
#STAD_261_G>A
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=gsub("\\+","",as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.character(e[,3])
e[,4]=as.numeric(as.character(e[,4]))
e[,5]=as.numeric(as.character(e[,5]))
e[,6]=paste(e[,2],e[,3],sep="_")
e1=e[e[,1]=="STAD",]
e1[,7]="STAD_261_G>A_no"
e1[e1[,6]=="261_G>A",7]="STAD_261_G>A"
cox1=coxph(Surv(e1[,5],e1[,4])~e1[,7])
write.table(paste("STAD_261_G>A",nrow(e1[e1[,7]=="STAD_261_G>A",]),nrow(e1[e1[,7]!="STAD_261_G>A",]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)



##Allcancer对hotspot加HR 113_T>G\126_G>A\132_G>T\158_G>T\164_A>G\237_G>T\248_G>C\261_G>A\298_G>T\307_A>G
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),gsub("[0-9]|_|+|+|-|-]{1}.*","",substring(as.character(f$HGVSc),4)),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("project_id","Muposition","Mutype","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.numeric(as.character(f1[,2]))
f1[,3]=gsub("\\+","",as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
f1[,6]=as.numeric(as.character(f1[,6]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$project_id),b1$Muposition,b1$Mutype,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$project_id),b2$Muposition,b2$Mutype,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.numeric(as.character(e[,2]))
e[,3]=as.character(e[,3])
e[,4]=as.numeric(as.character(e[,4]))
e[,5]=as.numeric(as.character(e[,5]))
e[,6]=paste(e[,2],e[,3],sep="_")
e1=e
e1[,7]="All_307_A>G_no"
e1[e1[,6]=="307_A>G",7]="All_307_A>G"
cox1=coxph(Surv(e1[,5],e1[,4])~e1[,7])
write.table(paste("All_307_A>G",nrow(e1[e1[,7]=="All_307_A>G",]),nrow(e1[e1[,7]!="All_307_A>G",]),summary(cox1)$coefficients[,2],summary(cox1)$conf.int[,3],summary(cox1)$conf.int[,4],summary(cox1)$coefficients[,5],sep="\t"),
            file =paste(path,"Position_SampleNumber/MuhotspotHR.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)




##对每类显著的HR加boxplot 
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Position_SampleNumber/"
f=read.table(paste(path,"MuhotspotHR_signif.txt",sep=""),header=T,sep="\t")
f[,1]=as.character(f[,1])
f[,8]=1:nrow(f)
pdf(file=paste(path,"MuhotspotHR_signif.pdf",sep=""),width = 10,height=9)
plot(f[,c(4,8)],pch=15,col= "black",cex = 2,las=1,xlim=c(min(f[,5])-2,max(f[,6])+5),ylim=c(min(f[,8])-0.2,max(f[,8])+0.2),
     yaxt="n",main ="Muhotspot",xlab="Hazard Ratio(95% CI)",ylab="Class")
abline(v=1,lwd=1,col="grey",lty=2)#虚线
text(f[,6]+1,f[,8],cex=0.8,signif(f[,7],3))
text(f[,5]-1.8,f[,8],cex=0.8,paste(f[,1],"=",as.character(f[,2]),sep=""))
text(f[,6]+3.5,f[,8],cex=0.8,paste(round(f[,4],2),"(",round(f[,5],2),",",round(f[,6],2),")",sep=""))
segments(f[,5],f[,8],f[,6],f[,8],col="black",lwd=2)
dev.off() 


##HNSC对domain296_310的样本号
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"3891.txt",sep=""),header=T,sep="\t")
f1=data.frame(na.omit(cbind(as.character(f[,1]),as.character(f[,5]),as.numeric(gsub("[a-z|A-Z|\\*|_]{1}.*","",substring(as.character(f$HGVSp_Short),4))),f$vital_status.y,f$days_to_last_follow_up,f$days_to_death.y)))
colnames(f1)=c("submitter_id","project_id","Muposition","vital_status","days_to_last_follow_up","days_to_death.y")
f1[,1]=as.character(f1[,1])
f1[,2]=as.character(f1[,2])
f1[,6]=as.numeric(as.character(f1[,6]))
f1[,3]=as.numeric(as.character(f1[,3]))
f1[,4]=as.numeric(as.character(f1[,4]))
f1[,5]=as.numeric(as.character(f1[,5]))
b1=subset(f1,f1$vital_status==1)
b2=subset(f1,f1$vital_status==2)
alive=cbind(as.character(b1$submitter_id),as.character(b1$project_id),b1$Muposition,b1$vital_status,b1$days_to_last_follow_up)
dead=cbind(as.character(b2$submitter_id),as.character(b2$project_id),b2$Muposition,b2$vital_status,b2$days_to_death.y)	
e=data.frame(rbind(alive,dead))
e[,1]=as.character(e[,1])
e[,2]=as.character(e[,2])
e[,5]=as.numeric(as.character(e[,5]))
e[,3]=as.numeric(as.character(e[,3]))
e[,4]=as.numeric(as.character(e[,4]))
e1=e[((e[,3])>=296) & ((e[,3])<310),]
write.table(e1,file =paste(path,"79sample.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)











#############################################
##cancer22的年龄分布,对每个cancer分为突变1、不突变样本0，用ks检验年龄是否一致
library(ggplot2)
cancer22=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
           "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Result2.txt",sep=""),header=T,sep="\t")
f1=f[as.character(f$project_id ) %in% cancer22,c("project_id","TP53_mutation","days_to_birth")]
f2=data.frame(cbind(class=paste(f1[,1],f1[,2],sep=""),age=-as.numeric(as.character(f1[,3]))/365))
f2[,2]=as.numeric(as.character(f2[,2]))
f2[,1]=as.character(f2[,1])
f2=na.omit(f2)
for(i in 1:length(cancer22)){
  f3=f2[f2[,1] %in% c(paste(cancer22[i],0,sep=""),paste(cancer22[i],1,sep="")),]
  d=ks.test(f2[f2[,1] %in% paste(cancer22[i],0,sep=""),2],f2[f2[,1] %in% paste(cancer22[i],1,sep=""),2])
  #p<-ggplot(f3,aes(x = f3[,2],colour=class))+geom_density(aes(fill=class),alpha=0.1)+
    #scale_fill_manual(values=c("#12B826","#1DADB8"))+scale_color_manual(values=c("#12B826","#1DADB8"))
  p<-ggplot(f3,aes(x = f3[,2],color=class))+geom_density(lwd=1.5)+
    scale_color_manual(values=c("#12B826","#1DADB8"))
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
  p <- p+xlab("age")+ylab("density")+ggtitle("cancer")
  p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="black",size=4,label=paste("p=",signif(d$p.value,3),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust =12,vjust = 5,colour="black",size=4,label=paste(length(f2[f2[,1] %in% paste(cancer22[i],0,sep=""),2]),length(f2[f2[,1] %in% paste(cancer22[i],1,sep=""),2]),sep="_"))
  ggsave(file=paste(path,"Position_SampleNumber/age_density/",cancer22[i],".pdf",sep=""),width = 10,height=9)
}

##cancer7的ks的p值显著
library(ggplot2)
cancer7=c("BRCA","CESC","ESCA","GBM","HNSC","LUAD","UCEC")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Result2.txt",sep=""),header=T,sep="\t")
f1=f[as.character(f$project_id ) %in% cancer7,c("project_id","TP53_mutation","days_to_birth")]
f2=data.frame(cbind(as.character(f1[,1]),class=paste(f1[,1],f1[,2],sep=""),age=-as.numeric(as.character(f1[,3]))/365))
f2[,3]=as.numeric(as.character(f2[,3]))
f2[,1]=as.character(f2[,1])
f2[,2]=as.character(f2[,2])
f2=na.omit(f2)
p<-ggplot(f2,aes(x = f2[,3],colour=class))+geom_density(aes(fill=class),alpha=0.6)+facet_grid(f2[,1]+f2[,2] ~ .)+
  scale_fill_manual(values=rep(c("#929591","#650021"),7))+scale_color_manual(values=rep(c("#929591","#650021"),7))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("age")+ylab("density")+ggtitle("cancer")
ggsave(file=paste(path,"Position_SampleNumber/age_density/significant1.pdf",sep=""),width = 10,height=9)


f3=f2[grep("1",f2[,2]),]
f4=f2[grep("0",f2[,2]),]
p<-ggplot(f3,aes(x = f3[,3],colour=class))+geom_density(aes(fill=class),alpha=0.6)+facet_grid(f3[,1]+f3[,2] ~ .)+
  scale_fill_manual(values=rep("#650021",7))+scale_color_manual(values=rep("#650021",7))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("age")+ylab("density")+ggtitle("cancer")
ggsave(file=paste(path,"Position_SampleNumber/age_density/significant1_mu.pdf",sep=""),width = 10,height=9)
p<-ggplot(f4,aes(x = f4[,3],colour=class))+geom_density(aes(fill=class),alpha=0.6)+facet_grid(f4[,1]+f4[,2] ~ .)+
  scale_fill_manual(values=rep("#929591",7))+scale_color_manual(values=rep("#929591",7))+scale_y_reverse()
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("age")+ylab("density")+ggtitle("cancer")
ggsave(file=paste(path,"Position_SampleNumber/age_density/significant1_nomu.pdf",sep=""),width = 10,height=9)


#############################################
##cancer22的年龄分布,对每个cancer+Allcancer分为突变1、不突变样本0，用t检验days_to_birth负数=age_at_diagnosis正数
library(ggplot2)
cancer22=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
           "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Result2.txt",sep=""),header=T,sep="\t")
f1=f[as.character(f$project_id ) %in% cancer22,c("project_id","TP53_mutation","days_to_birth")]
f2=data.frame(cbind(class=paste(f1[,1],f1[,2],sep=""),age=-as.numeric(as.character(f1[,3]))/365))
f2[,2]=as.numeric(as.character(f2[,2]))
f2[,1]=as.character(f2[,1])
f2=na.omit(f2)
for(i in 1:length(cancer22)){
  f3=f2[f2[,1] %in% c(paste(cancer22[i],0,sep=""),paste(cancer22[i],1,sep="")),]
  d=t.test(f2[f2[,1] %in% paste(cancer22[i],0,sep=""),2],f2[f2[,1] %in% paste(cancer22[i],1,sep=""),2])
  p<-ggplot(f3,aes(x = f3[,1],y = f3[,2]))+geom_boxplot(aes(fill=class,color=class),alpha=0.5,lwd=0.6,outlier.colour = NA)+
    scale_fill_manual(values=c("#12B826","#1DADB8"))+scale_color_manual(values=c("#12B826","#1DADB8"))
  p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
  p <- p+xlab("Different Samples")+ylab("Age in years")+ggtitle("Age Distribution in Different Samples")
  p <- p+annotate("text",x=Inf, y=Inf,hjust =7,vjust = 5,colour="black",size=4,label=paste("p=",signif(d$p.value,3),sep=""))+
    annotate("text",x=Inf, y=Inf,hjust =12,vjust = 5,colour="black",size=2.5,label=paste(length(f2[f2[,1] %in% paste(cancer22[i],0,sep=""),2]),length(f2[f2[,1] %in% paste(cancer22[i],1,sep=""),2]),sep="_"))
  ggsave(file=paste(path,"Position_SampleNumber/age_density_boxplot/",cancer22[i],".pdf",sep=""),width = 10,height=9)
}

##Allcancer\days_to_birth
library(ggplot2)
cancer22=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
  "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Result2.txt",sep=""),header=T,sep="\t")
f1=f[as.character(f$project_id ) %in% cancer22,c("project_id","TP53_mutation","days_to_birth")]
f1[,1]="Allcancer"
f1[,3]=-as.numeric(as.character(f1[,3]))/365
f1=na.omit(f1)
d=t.test(f1[which(f1[,2]==1),3],f1[which(f1[,2]==0),3])
f1[,2]=as.factor(f1[,2])
p<-ggplot(f1,aes(x = f1[,2],y = f1[,3]))+geom_violin(aes(fill=TP53_mutation,color=TP53_mutation),alpha=0.5)+
  geom_boxplot(aes(color=TP53_mutation),alpha=0.5,lwd=0.6,width=0.3,outlier.colour = NA)+
  scale_fill_manual(values=c("#12B826","#1DADB8"))+scale_color_manual(values=c("#12B826","#1DADB8"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("Different Samples")+ylab("Age in years")+ggtitle("Age Distribution in Different Samples")
p <- p+annotate("text",x=Inf, y=Inf,hjust =3,vjust = 5,colour="black",size=4,label=paste("p=",signif(d$p.value,3),sep=""))
ggsave(file=paste(path,"Position_SampleNumber/age_density_boxplot/Allcancer1.pdf",sep=""),width = 10,height=9)



##显著的有Allcancer\CESC\ESCA\GBM\LUAD\UCEC#cancer22=c("CESC","ESCA","GBM","LUAD","UCEC")
library(ggplot2)
cancer22=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
           "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f= read.table(paste(path,"Result2.txt",sep=""),header=T,sep="\t")
f1=f[as.character(f$project_id ) %in% cancer22,c("project_id","TP53_mutation","days_to_birth")]
f1[,1]="Allcancer"
f1[,3]=-as.numeric(as.character(f1[,3]))/365
f1=na.omit(f1)
f1[,2]=as.factor(f1[,2])  #Allcancer
f2=f[as.character(f$project_id ) %in% c("CESC","ESCA","GBM","LUAD","UCEC") ,c("project_id","TP53_mutation","days_to_birth")]
f2[,3]=-as.numeric(as.character(f2[,3]))/365
f2[,2]=as.factor(f2[,2]) 
f2=na.omit(f2)
f3=rbind(f1,f2)
p<-ggplot(f3,aes(x = f3[,1],y = f3[,3]))+geom_boxplot(aes(fill=TP53_mutation,color=TP53_mutation),width=0.5,alpha=0.5,lwd=0.6,outlier.colour = NA)+
  scale_fill_manual(values=c("#12B826","#1DADB8"))+scale_color_manual(values=c("#12B826","#1DADB8"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("Different Samples in Cancers")+ylab("Age in years")+ggtitle("Signifcant Age Distribution in Different Samples")
ggsave(file=paste(path,"Position_SampleNumber/age_density_boxplot/signifcant.pdf",sep=""),width = 10,height=9)

 



####根据分好的4亚型，看病人的年龄分布days_to_birth
##Allcancer2_class4_manhattan_ward.D_1075 
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f= read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
f1=f[,c("submitter_id","days_to_birth")]
f1[,1]=gsub("-",".",as.character(f1[,1]))
f1[,2]=-as.numeric(as.character(f1[,2]))/365
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
m=merge(f2,f1,by.x="V1",by.y="submitter_id",all.x=T)
m[,2]=as.factor(m[,2])
p<-ggplot(m,aes(x = m[,2],y = m[,3]))+geom_violin(aes(fill=V2,color=V2),alpha=0.5)+
  geom_boxplot(aes(fill=V2,color=V2),width=0.1,alpha=0.5,lwd=0.6,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+scale_color_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("Different Subtypes in Cancers")+ylab("Age in years")+ggtitle("Age Distribution in Different Subtypes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/24Ddays_to_birth.pdf",sep=""),width = 10,height=9)


####根据分好的4亚型，看病人的年龄分布days_to_birth +Allcancer的Nomu 5类
##Allcancer2_class4_manhattan_ward.D_1075 
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f= read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
f1=f[,c("submitter_id","days_to_birth")]
f1[,1]=gsub("-",".",as.character(f1[,1]))
f1[,2]=-as.numeric(as.character(f1[,2]))/365
f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",sep=""),header=F,sep="\t")
m=merge(f2,f1,by.x="V1",by.y="submitter_id",all.x=T)
m[,2]=as.factor(m[,2])
p<-ggplot(m,aes(x = m[,2],y = m[,3]))+geom_violin(aes(fill=V2,color=V2),alpha=0.5)+
  geom_boxplot(aes(fill=V2,color=V2),width=0.1,alpha=0.5,lwd=0.6,outlier.colour = NA)+
  scale_fill_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))+scale_color_manual(values=c("deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("Different Subtypes in Cancers")+ylab("Age in years")+ggtitle("Age Distribution in Different Subtypes")
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/24Ddays_to_birth.pdf",sep=""),width = 10,height=9)


cancer22=c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
           "LUSC","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","THYM","UCEC")
f3=f[as.character(f$project_id ) %in% cancer22,c("TP53_mutation","days_to_birth")]
f3[,2]=-as.numeric(as.character(f3[,2]))/365
f4=f3[which(f3[,1]==0),]
f4[,1]=as.factor(as.character(f4[,1]))
colnames(f4)[1]="V2"
m1=rbind(m[,2:3],f4)
m1[,1]=factor(m1[,1],c(0,1,2,3,4))
d1=t.test(m1[which(as.numeric(as.character(m1[,1]))==0),2],m1[which(as.numeric(as.character(m1[,1]))==1),2])
d2=t.test(m1[which(as.numeric(as.character(m1[,1]))==0),2],m1[which(as.numeric(as.character(m1[,1]))==2),2])
d3=t.test(m1[which(as.numeric(as.character(m1[,1]))==0),2],m1[which(as.numeric(as.character(m1[,1]))==3),2])
d4=t.test(m1[which(as.numeric(as.character(m1[,1]))==0),2],m1[which(as.numeric(as.character(m1[,1]))==4),2])
p<-ggplot(m1,aes(x = m1[,1],y = m1[,2]))+geom_violin(aes(fill=V2,color=V2),alpha=0.5)+
  geom_boxplot(aes(fill=V2,color=V2),width=0.1,alpha=0.5,lwd=0.6,outlier.colour = NA)+
  scale_fill_manual(values=c("#12B826","deeppink3","goldenrod3","palegreen4","royalblue4"))+scale_color_manual(values=c("#12B826","deeppink3","goldenrod3","palegreen4","royalblue4"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=6,colour="black"),strip.text=element_text(size=5,colour="black"))
p <- p+xlab("Different Subtypes in Cancers")+ylab("Age in years")+ggtitle("Age Distribution in Different Subtypes")+
  annotate("text",x=Inf, y=Inf,hjust =10,vjust = 4,colour="black",size=4,label=signif(d1$p.value,3))+
  annotate("text",x=Inf, y=Inf,hjust =8,vjust = 4,colour="black",size=4,label=signif(d2$p.value,3))+
  annotate("text",x=Inf, y=Inf,hjust =6.5,vjust = 4,colour="black",size=4,label=signif(d3$p.value,3))+
  annotate("text",x=Inf, y=Inf,hjust =3,vjust = 4,colour="black",size=4,label=signif(d4$p.value,3))
ggsave(file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/24Ddays_to_birth1.pdf",sep=""),width = 10,height=9)





