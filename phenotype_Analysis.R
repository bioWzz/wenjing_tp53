
####一、27cancerGDC_phenotype要去冗余、GDC_survival不需要##取显著p值
path="D:/tp53/"
#cancer= "COAD","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PCPG","READ","SKCM","TGCT","UCEC",
#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML",
#"LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG",
#"PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
###GDC_phenotype
for(cancer in c("THCA","THYM","UCEC","UCS")){
  f1= read.csv(paste(path,"GDC_phenotype/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
  f2= read.csv(paste(path,"GDC_phenotype/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
  f3= read.csv(paste(path,"GDC_phenotype/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
  f4= read.csv(paste(path,"GDC_phenotype/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
  write.table(paste("Variables","NoMutation","Mutation","SingleMutation","CoMutation","No_Mutation_pValue","Single_CoMutation_pValue",sep="\t"),file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
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
      if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
        write.table(g1,file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
        write.table(g2,file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
      
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
      if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
        write.table(g1,file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
        write.table(b,file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
      
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
      if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
        write.table(g1,file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
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
        write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/11/",cancer,"1.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
        
    }
  }
}


#######################################################
##对Allcancer检查突变、没突变、单突变、多突变的initial_weight.samples的boxplot
path="D:/tp53/"
data1=data.frame()
data2=data.frame()
data3=data.frame()
data4=data.frame()
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC",
                 "LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA",
                "THYM","UCEC","UCS")){
  f1= read.csv(paste(path,"GDC_phenotype/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
  f2= read.csv(paste(path,"GDC_phenotype/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
  f3= read.csv(paste(path,"GDC_phenotype/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
  f4= read.csv(paste(path,"GDC_phenotype/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
  data1=append(data1,f1$initial_weight.samples)
  data2=append(data2,f2$initial_weight.samples)
  data3=append(data3,f3$initial_weight.samples)
  data4=append(data4,f4$initial_weight.samples)
}
m1=na.omit(data.frame(cbind(class="NoMu",t(data.frame(data1)))))
m2=na.omit(data.frame(cbind(class="TP53Mu",t(data.frame(data2)))))
m3=na.omit(data.frame(cbind(class="SingleMu",t(data.frame(data3)))))
m4=na.omit(data.frame(cbind(class="CoMu",t(data.frame(data4)))))
m=rbind(m1,m2,m3,m4)
m[,2]=as.numeric(as.character(m[,2]))
pdf(paste(path,"phenotype_Analysis/","Allcancer_initial_weight.samples.pdf",sep=""))
boxplot( V2~ class,m, pch=".",lwd=2,border= c("#12B826","#1DADB8","#5e227f","#d22780"),main="The initial weight in samples",ylab="weight",xlab="Sample Types",las=1)
text(1.5,3000,labels=signif(t.test(as.numeric(as.character(m1[,2])),as.numeric(as.character(m2[,2])))$p.value,3),cex=1,pos=3)
text(3.5,3000,labels=signif(t.test(as.numeric(as.character(m3[,2])),as.numeric(as.character(m4[,2])))$p.value,3),cex=1,pos=3)
dev.off()








######################################################
###GDC_survival##取显著p值
path="D:/tp53/"
for(cancer in c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC",
                "KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD",
                "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")){
  f5= read.csv(paste(path,"GDC_survival/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
  f6= read.csv(paste(path,"GDC_survival/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
  f7= read.csv(paste(path,"GDC_survival/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
  f8= read.csv(paste(path,"GDC_survival/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
  write.table(paste("Variables","NoMutation","Mutation","SingleMutation","CoMutation","No_Mutation_pValue","Single_CoMutation_pValue",sep="\t"),file =paste(path,"phenotype_Analysis/22/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
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
      if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
        write.table(g1,file =paste(path,"phenotype_Analysis/22/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
        write.table(g2,file =paste(path,"phenotype_Analysis/22/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
      
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
      if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
        write.table(g1,file =paste(path,"phenotype_Analysis/22/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
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
        write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/22/",cancer,"2.txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
      }
      
    }
  }
}


###############################################
##额外的分析、对ESCA的性别划分后再分年龄
path="D:/tp53/"
cancer="ESCA"
f1= read.csv(paste(path,"GDC_phenotype/",cancer,"_NoMutation.txt",sep=""),header=T,sep="\t")
f2= read.csv(paste(path,"GDC_phenotype/",cancer,"_Mutation.txt",sep=""),header=T,sep="\t")
f3= read.csv(paste(path,"GDC_phenotype/",cancer,"_SingleMutation.txt",sep=""),header=T,sep="\t")
f4= read.csv(paste(path,"GDC_phenotype/",cancer,"_CoMutation.txt",sep=""),header=T,sep="\t")
write.table(paste("Variables","NoMutation","Mutation","SingleMutation","CoMutation","No_Mutation_pValue","Single_CoMutation_pValue",sep="\t"),file =paste(path,"phenotype_Analysis/",cancer,".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
i="gender.demographic"
j="age_at_initial_pathologic_diagnosis"
a1=data.frame(unique(cbind(f1[,colnames(f1)==j],gsub("^$|i",NA,as.character(f1[,colnames(f1)==i])))))
a2=data.frame(unique(cbind(f2[,colnames(f2)==j],gsub("^$|i",NA,as.character(f2[,colnames(f2)==i]))))) 
a3=data.frame(unique(cbind(f3[,colnames(f3)==j],gsub("^$|i",NA,as.character(f3[,colnames(f3)==i]))))) 
a4=data.frame(unique(cbind(f4[,colnames(f4)==j],gsub("^$|i",NA,as.character(f4[,colnames(f4)==i]))))) 
b1=as.numeric(as.character(a1[as.character(a1[,2])=="male",1]))
b2=as.numeric(as.character(a2[as.character(a2[,2])=="male",1]))
b3=as.numeric(as.character(a3[as.character(a3[,2])=="male",1]))
b4=as.numeric(as.character(a4[as.character(a4[,2])=="male",1]))
e1=as.numeric(as.character(a1[as.character(a1[,2])=="female",1]))
e2=as.numeric(as.character(a2[as.character(a2[,2])=="female",1]))
e3=as.numeric(as.character(a3[as.character(a3[,2])=="female",1]))
e4=as.numeric(as.character(a4[as.character(a4[,2])=="female",1]))  
 
##连续型male_age
if(length(as.numeric(b1))>=2 & length(as.numeric(b2))>=2){
  p1=t.test(as.numeric(b1),as.numeric(b2))$p.value  #NoMutation_Mutation
}else{p1=-5}
if(length(as.numeric(b3))>=2 & length(as.numeric(b4))>=2){
  p2=t.test(as.numeric(b3),as.numeric(b4))$p.value  #SingleMutation_CoMutation
}else{p2=-5}
g1=paste("male_age",paste(round(mean(as.numeric(b1)),4),round(sd(as.numeric(b1)),4),sep="+"),paste(round(mean(as.numeric(b2)),4),round(sd(as.numeric(b2)),4),sep="+"),paste(round(mean(as.numeric(b3)),4),round(sd(as.numeric(b3)),4),sep="+"),paste(round(mean(as.numeric(b4)),4),round(sd(as.numeric(b4)),4),sep="+"),round(p1,4),round(p2,4),sep="\t")
if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
  write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  d<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段
  NO1<-(max(as.numeric(b1))-min(as.numeric(b1)))%/%3+min(as.numeric(b1)) #min~NO1
  NO2<-(max(as.numeric(b1))-min(as.numeric(b1)))%/%3*2+min(as.numeric(b1))#NO1~NO2#NO2~max
  row.names(d)=c(paste(min(as.numeric(b1)),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(as.numeric(b1))))
  d1=as.data.frame(as.numeric(b1))
  d1[,2]=1
  da1={}
  for(j in row.names(d)){
    data=d1[as.numeric(d1[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d1[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da1=append(da1,sum(data[,2]))
  }
  names(da1)=row.names(d)
  da1<-as.data.frame(da1)
  d2=as.data.frame(as.numeric(b2))
  d2[,2]=1
  da2={}
  for(j in row.names(d)){
    data=d2[as.numeric(d2[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d2[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da2=append(da2,sum(data[,2]))
  }
  names(da2)=row.names(d)
  da2<-as.data.frame(da2)
  d3=as.data.frame(as.numeric(b3))
  d3[,2]=1
  da3={}
  for(j in row.names(d)){
    data=d3[as.numeric(d3[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d3[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da3=append(da3,sum(data[,2]))
  }
  names(da3)=row.names(d)
  da3<-as.data.frame(da3)
  if(length(b4)!=0){
    d4=as.data.frame(as.numeric(b4))
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
  write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/",cancer,".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
}  


##连续型female_age
if(length(as.numeric(e1))>=2 & length(as.numeric(e2))>=2){
  p1=t.test(as.numeric(e1),as.numeric(e2))$p.value  #NoMutation_Mutation
}else{p1=-5}
if(length(as.numeric(e3))>=2 & length(as.numeric(e4))>=2){
  p2=t.test(as.numeric(e3),as.numeric(e4))$p.value  #SingleMutation_CoMutation
}else{p2=-5}
g1=paste("female_age",paste(round(mean(as.numeric(e1)),4),round(sd(as.numeric(e1)),4),sep="+"),paste(round(mean(as.numeric(e2)),4),round(sd(as.numeric(e2)),4),sep="+"),paste(round(mean(as.numeric(e3)),4),round(sd(as.numeric(e3)),4),sep="+"),paste(round(mean(as.numeric(e4)),4),round(sd(as.numeric(e4)),4),sep="+"),round(p1,4),round(p2,4),sep="\t")
if((p1>=0 & p1<0.05) | (p2>=0 & p2 <0.05)){
  write.table(g1,file =paste(path,"phenotype_Analysis/",cancer,".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
  d<-as.data.frame(matrix(numeric(0),nrow =3,ncol=4)) #分为3段
  NO1<-(max(as.numeric(e1))-min(as.numeric(e1)))%/%3+min(as.numeric(e1)) #min~NO1
  NO2<-(max(as.numeric(e1))-min(as.numeric(e1)))%/%3*2+min(as.numeric(e1))#NO1~NO2#NO2~max
  row.names(d)=c(paste(min(as.numeric(e1)),"~",NO1),paste(NO1,"~",NO2),paste(NO2,"~",max(as.numeric(e1))))
  d1=as.data.frame(as.numeric(e1))
  d1[,2]=1
  da1={}
  for(j in row.names(d)){
    data=d1[as.numeric(d1[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d1[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da1=append(da1,sum(data[,2]))
  }
  names(da1)=row.names(d)
  da1<-as.data.frame(da1)
  d2=as.data.frame(as.numeric(e2))
  d2[,2]=1
  da2={}
  for(j in row.names(d)){
    data=d2[as.numeric(d2[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d2[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da2=append(da2,sum(data[,2]))
  }
  names(da2)=row.names(d)
  da2<-as.data.frame(da2)
  d3=as.data.frame(as.numeric(e3))
  d3[,2]=1
  da3={}
  for(j in row.names(d)){
    data=d3[as.numeric(d3[,1])>=as.numeric(strsplit(j,"~")[[1]][1]) & as.numeric(d3[,1])<as.numeric(strsplit(j,"~")[[1]][2]),] 
    da3=append(da3,sum(data[,2]))
  }
  names(da3)=row.names(d)
  da3<-as.data.frame(da3)
  if(length(e4)!=0){
    d4=as.data.frame(as.numeric(e4))
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
  write.table(cbind(rownames(da),da),file =paste(path,"phenotype_Analysis/",cancer,".txt",sep=""), row.names = F,col.names=F, quote = F,sep = "\t",append=TRUE)
} 



########################################################
x<-matrix(c(61,28,77,153),nrow=2,byrow = T)
signif(fisher.test(x)$p.value,3)


###################################################
##补充
##对OV，根据OV_p53_context_mutation_fpkm\OV_p53_context_nomutation_fpkm
##把样本分为p53\WT做生存分析
library("survival")
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/"
f<-read.csv(paste(path,"Result2.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
f$days_to_death.y<-as.numeric(gsub("NULL",NA,f$days_to_death.y))
f$vital_status.y<-gsub("alive",0,f$vital_status.y)
f$vital_status.y<-as.numeric(gsub("dead",1,f$vital_status.y))
rownames(f)=f[,1]
f=na.omit(f)

f1<-read.csv(paste(path,"p53_context_mutation_fpkm/OV_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
f2<-read.csv(paste(path,"p53_context_nomutation_fpkm/OV_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
WT=gsub("\\.","-",colnames(f2))[-1]
p53=gsub("\\.","-",colnames(f1))[-1]

pdf(paste(path,"p53_context_fpkm_foldchange_survival/OV_survivalplot.pdf",sep=""))
ff1=rbind(cbind(class=1,f[WT,]),cbind(class=5,f[p53,]))
b11=subset(ff1,ff1$vital_status.y==0)
b21=subset(ff1,ff1$vital_status.y==1)
alive1=cbind(b11$class,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$class,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=na.omit(e1)
dif1 <- survdiff(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1])#求生存时间
kmsurvival1<-survfit(Surv(as.numeric(e1[,3]),as.numeric(e1[,2]))~e1[,1],conf.type = "log-log")
plot(kmsurvival1, lty = 'solid', col=c('black','orange'),
     xlab='survival time in days',ylab='survival probabilities',main="OV")
legend('bottomleft', cex=0.6,text.width=0.4,c('WT','p53'), lty='solid',
       col=c('black','orange'))
text(1800,1,cex=0.8,paste(paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""),paste(dif1$n[1],dif1$n[2],sep="_"),sep="\n"))
dev.off()
 





####################################
##找出154临床基因的fpkm值
path="D:/tp53/clinical_gene_drug/"
path1="D:/tp53/Gene Expression/"
f<-read.table(paste(path,"154gene.txt",sep=""),header=T,sep="\t")
filename<-list.files(paste(path1,"Gene_FPKM1/",sep=""))
for(i in filename){
  f1<-read.table(paste(path1,"Gene_FPKM1/",i,sep=""),header=T,sep="\t")
  a=f1[as.character(f1[,1]) %in% as.character(f[,1]),]
  write.table(a,file=paste("D:/tp53/clinical_gene_drug/Gene_FPKM1_Analysis/clinic_gene_fpkm/",gsub("_gene","_clinic_gene",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}

##找出154临床基因在突变、没突变、正常样本中的fpkm值 
path="D:/tp53/clinical_gene_drug/Gene_FPKM1_Analysis/"
f<-read.table("D:/tp53/Gene Expression/3894_5886=9780.txt",header=T,sep="\t")
filename<-list.files(paste(path,"clinic_gene_fpkm/",sep=""))
for(i in filename){
  f1<-read.table(paste(path,"clinic_gene_fpkm/",i,sep=""),header=T,sep="\t")
  a1<-f1[,c(1,grep(".*.11$",colnames(f1)))]
  mutation<-subset(f,as.character(f[,2])==1 & as.character(f[,5])==gsub("_clinic_gene_fpkm.txt","",i))[,1]
  nomutation<-subset(f,is.na(as.character(f[,2]))& as.character(f[,5])==gsub("_clinic_gene_fpkm.txt","",i))[,1]
  colnames(f1)<-gsub("\\.","-",substring(colnames(f1),1,12))
  a2<-f1[,c("gene_Name",intersect(as.character(mutation),colnames(f1)))]
  a3<-f1[,c("gene_Name",intersect(as.character(nomutation),colnames(f1)))]
  write.table(a1,file=paste(path,"clinic_gene_normal_fpkm/",gsub("_fpkm","_normal_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a2,file=paste(path,"clinic_gene_mutation_fpkm/",gsub("_fpkm","_mutation_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
  write.table(a3,file=paste(path,"clinic_gene_nomutation_fpkm/",gsub("_fpkm","_nomutation_fpkm",i),sep=""),row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
}




#c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LIHC","LUAD",
#"LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","SKCM","THCA","THYM","UCEC","LAML","LGG","OV","TGCT","UCS")
##对每个cancer，求TP53上下文基因和154临床基因的cor相关系数（在mutation样本中），行列都为基因，画热图
##PCPG、TGCT只有一个样本，删去
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_mutation_fpkm/"
path1="D:/tp53/clinical_gene_drug/Gene_FPKM1_Analysis/clinic_gene_mutation_fpkm/"
f= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
#cancer="SKCM"
for(cancer in c("UCS")){
  f1<-read.table(paste(path,cancer,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  rownames(f1)=f1$gene_Name
  f2<-read.table(paste(path1,cancer,"_clinic_gene_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f3=data.frame(c(gene_Name="ICOSLG",apply(f2[as.character(f2[,1])=="ICOSLG",2:ncol(f2)],2,mean)))  ##重复行取均值
  f4=rbind(f2[as.character(f2[,1])!="ICOSLG",],t(f3))
  rownames(f4)=f4$gene_Name
  da=as.data.frame(matrix(numeric(0),nrow =nrow(f1),ncol=nrow(f4))) 
  rownames(da)=f1$gene_Name
  colnames(da)=f4$gene_Name
  for(i in 1:nrow(da)){
    for(j in 1:ncol(da)){
      da[i,j]=cor(as.numeric(f1[rownames(da)[i],2:ncol(f1)]),as.numeric(f4[colnames(da)[j],2:ncol(f4)]))
      #f1[rownames(da)[i],2:ncol(f1)]
      #f2[as.character(f2[,1])==colnames(da)[j],2:ncol(f2)]
    }
  }
  da1=da[intersect(as.character(f[,1]),rownames(da)),]
  y<-as.matrix(da1)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  pheatmap(y,cluster_row = FALSE,cluster_col = FALSE,col=colorsChoice(3),legend_breaks=c(min(y),-0.5,0.5,max(y)),breaks=c(min(y),-0.5,0.5,max(y)),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste("D:/tp53/clinical_gene_drug/Gene_FPKM1_Analysis/heatmap/",cancer,".pdf",sep=""),width=8,height=8)  
}

 
#############################################
##TCGA molecular subtypes下载，用TCGAbiolinks包里的功能PanCancerAtlas_subtypes()下载每个cancer的样本所对应的亚型。
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")
install.packages("DT")
subtypes <- PanCancerAtlas_subtypes()
#DT::datatable(subtypes,
              #filter = 'top',
              #options = list(scrollX = TRUE, keys = TRUE),
              #rownames = FALSE)  #没用
write.table(subtypes,file="D:/tp53/phenotype_Analysis/subtype/subtype.txt",row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)

#统计每个cancer所对应的亚型的总的样本数和发生TP53突变的样本数。
path="D:/tp53/phenotype_Analysis/subtype/"
f1<-read.csv(paste(path,"subtype.txt",sep=""),header=T,sep="\t")
f1[,1]=substring(as.character(f1[,1]),1,12)
f<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
f2=merge(f,f1,by.x="submitter_id",by.y="pan.samplesID",all.x=T)
write.table(paste("cancer","subtype","AllMu","TP53Mu",sep="\t"),file="D:/tp53/phenotype_Analysis/subtype/subtypecount.txt",row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
for(i in unique(f2[,2])){
  a=f2[f2[,2]==i,]
  for(j in unique(a$Subtype_Selected)){
    b1=nrow(a[a$Subtype_Selected==j,])
    b2=nrow(a[(a$TP53_mutation==1)&(a$Subtype_Selected==j),])
    b3=paste(i,j,b1,b2,sep="\t")
    write.table(b3,file="D:/tp53/phenotype_Analysis/subtype/subtypecount.txt",row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  }
}
 

###对显著的癌症GBM、LGG、BRCA的亚型做生存分析
library("survival")
path="D:/tp53/phenotype_Analysis/subtype/"
ff<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
ff$days_to_death.y<-as.numeric(gsub("NULL",NA,ff$days_to_death.y))
ff$vital_status.y<-gsub("alive",0,ff$vital_status.y)
ff$vital_status.y<-as.numeric(gsub("dead",1,ff$vital_status.y))
rownames(ff)=ff[,1]
ff=na.omit(ff)
#cancer="GBM"
#cancer="LGG"
cancer="BRCA"
ff1=ff[ff$project_id==cancer,c("submitter_id","vital_status.y","days_to_last_follow_up","days_to_death.y")]
b11=subset(ff1,ff1$vital_status.y==0)
b21=subset(ff1,ff1$vital_status.y==1)
alive1=cbind(b11$submitter_id,b11$vital_status.y,b11$days_to_last_follow_up)
dead1=cbind(b21$submitter_id,b21$vital_status.y,b21$days_to_death.y)	
e1=rbind(alive1,dead1)
e1=data.frame(na.omit(e1))

f1<-read.csv(paste(path,"subtype.txt",sep=""),header=T,sep="\t")
a=f1[as.character(f1$cancer.type)==cancer,c("pan.samplesID","Subtype_Selected")]
a[,1]=substring(as.character(a[,1]),0,12)
h=merge(a,e1,by.x="pan.samplesID",by.y="X1")
#unique(h[,2])
#h=h[as.character(h$Subtype_Selected)!="BRCA.Normal",]
#h=h[as.character(h$Subtype_Selected)!="GBM_LGG.NA",]
#h[,5]="other"
h[,5]=NA
#h[as.character(h$Subtype_Selected)=="GBM_LGG.G-CIMP-high",5]="GBM_LGG.G-CIMP-high"
#h[as.character(h$Subtype_Selected)=="GBM_LGG.G-CIMP-low",5]="GBM_LGG.G-CIMP-low"
h=na.omit(h)
pdf(paste(path,cancer,"_subtype_survival2.pdf",sep=""))  
dif1 <- survdiff(Surv(as.numeric(h[,4]),as.numeric(h[,3]))~h[,5])#求生存时间
kmsurvival1<-survfit(Surv(as.numeric(h[,4]),as.numeric(h[,3]))~h[,5],conf.type = "log-log")
plot(kmsurvival1,lty = 'solid', col=c("red","orange","blue","purple"),lwd=1.1,
     xlab='survival time in days',ylab='survival probabilities',main=cancer)
legend('bottomleft', cex=0.6,text.width=0.4,gsub(".*=","",as.character(data.frame(dif1$n)[,1])), lty='solid',
       col=c("red","orange","blue","purple"))
text(100,1,cex=0.8,paste("p=",round(pchisq(dif1$chisq,1,lower.tail=F),4),sep=""))
for(j in 1:length(as.numeric(unique(h[,5]))) ){
  text(100,0.9-(j-1)*0.03,cex=0.8,dif1$n[j])
}
dev.off()



 

###########################################
##与风险相关，对于每一癌症，对TP53上下文基因里的每一个基因，根据gene在癌症突变样本中的表达值，做单cox，得到HR（>1\<1）；
##根据gene在癌症突变样本中的表达值的中值分为2组，高表达和低表达，对2组样本做生存得到logrankp值；
##根据HR>1找到坏预后，看表达为高还是低；根据HR<1找到好预后，看表达为高还是低；
##最终得到矩阵画热饼图（行为基因列为癌症）：第一列为TP53突变状态+TP53上下文基因名字，第二列为HR值（边框红蓝），
##第三列为logrankp值（-logp为圆圈大小,p值越小越显著圆圈越大），
##第四列为gene在癌症突变、没突变样本中的表达值均值的fc（圆另一半标fc>1.5\<0.5为紫色，中间为绿色），
##第五列为高表达或低表达（圆一半标红蓝）
##共27-2=25cancer：PCPG17 1突变样本 TGCT23 2突变样本(删去，没法分高低表达组)
library("survival") #alive0dead1
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/"
f<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t",stringsAsFactors = F)
a=subset(f,f$vital_status.y=="alive")
b=subset(f,f$vital_status.y=="dead")
alive=cbind(a$submitter_id,a$vital_status.y,a$days_to_last_follow_up,"state"=0)
dead=cbind(b$submitter_id,b$vital_status.y,b$days_to_death.y,"state"=1)	
e=data.frame(rbind(alive,dead))
e[,1]=gsub("-",".",as.character(e[,1]))

filename<-gsub("_.*","",list.files(paste(path,"p53_context_mutation_fpkm/",sep="")))
for(i in filename[24:27]){
  f1<-read.table(paste(path,"p53_context_mutation_fpkm/",i,"_p53_context_mutation_fpkm.txt",sep=""),header=T,sep="\t")
  f2<-read.table(paste(path,"p53_context_nomutation_fpkm/",i,"_p53_context_nomutation_fpkm.txt",sep=""),header=T,sep="\t")
  dir.create(paste(path,"riskAnalysis/",i,sep="")) 
  write.table(paste("name","HR","p","fc",sep="\t"),file=paste(path,"riskAnalysis/",i,"/gene_risk.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
  for(j in 1:nrow(f1)){
    h1=colnames(f1)[f1[j,2:ncol(f1)]>=median(as.numeric(f1[j,2:ncol(f1)]))] ##高表达
    h2=colnames(f1)[f1[j,2:ncol(f1)]<median(as.numeric(f1[j,2:ncol(f1)]))][-1] ##低表达
    e1=e[e[,1] %in% c(h1,h2),]
    e1[e1[,1] %in% h1,5]="h1"
    e1[e1[,1] %in% h2,5]="h2"
    e1[,3]=as.numeric(as.character(e1[,3]))
    e1[,4]=as.numeric(as.character(e1[,4]))
    d=data.frame(t(f1[j,2:ncol(f1)]))
    d$name=rownames(d)
    e2=merge(e1,d,by.x="V1",by.y="name")
    d1=(summary(coxph(Surv(e2[,3],e2[,4])~e2[,6]))$coefficients)[,2] ##HR
    
    if(length(unique(e2[,5]))!=1){
      dif <- survdiff(Surv(e2[,3],e2[,4])~e2[,5])#求生存时间
      kmsurvival<-survfit(Surv(e2[,3],e2[,4])~e2[,5],conf.type = "log-log")
      d2=pchisq(dif$chisq,1,lower.tail=F)
      pdf(paste(path,"riskAnalysis/",i,"/",f1[j,1],".pdf",sep=""))
      plot(kmsurvival, lty = c('solid', 'dashed'), col=c('red','green'),main=paste(i,f1[j,1],sep="_"),
           xlab='survival time in days',ylab='survival probabilities')
      legend('bottomleft', cex=0.8,text.width=0.6,c('high','low'), lty=c('solid','dashed'),
             col=c('red','green'))
      text(1500,0.9,cex=0.5,paste(paste("p=",pchisq(dif$chisq,1,lower.tail=F),sep=""),dif$n[1],dif$n[2],sep="\n"))
      dev.off()
      
      d3=mean(as.numeric(f1[j,2:ncol(f1)]))/mean(as.numeric(f2[f2[,1]==f1[j,1],2:ncol(f2)]))
      write.table(paste(f1[j,1],d1,d2,d3,sep="\t"),file=paste(path,"riskAnalysis/",i,"/gene_risk.txt",sep=""),row.names = F,col.names = F, quote = F,sep = "\t",append=TRUE)
    }
  }
}


##共27-2=25cancer：PCPG17 1突变样本 TGCT23 2突变样本(删去，没法分高低表达组)
##根据矩阵画热饼图
#THCA23有1个样本，删去##结果为24cancer
library(ggplot2)
path="D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/riskAnalysis/"
ff= read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
filename<-list.files(path)
da=data.frame()
for(i in filename[c(2:22,24:26)]){
  f<-read.table(paste(path,i,"/gene_risk1.txt",sep=""),header=T,sep="\t")
  f[f[,2]>=1,6]="highHR"
  f[f[,2]<1,6]="lowHR"
  colnames(f)[6]="HR_state"
  f[(f[,4]>=1.5 | f[,4]<=0.5),7]="highfc"
  f[(f[,4]<1.5 & f[,4]>0.5),7]="lowfc"
  colnames(f)[7]="fc_state"
  f[,8]=0.5
  f[f[,3]<=0.05,9]="psign"
  f[f[,3]>0.05,9]="pnosign"
  colnames(f)[9]="Pvalue"
  f=merge(f,ff,by.x="name",by.y="gene_Name")
  f1=f[,c(1,6,5,9,8,10)]
  f2=f[,c(1,6,7,9,8,10)]
  colnames(f1)[3]="state"
  colnames(f2)[3]="state"
  f3=rbind(f1,f2)
  f4=cbind(cancer=i,f3)
  da=rbind(da,f4)
}
da[,2]=factor(da[,2],ff[,1])
da[,4]=factor(da[,4],c("high","low","lowfc","highfc"))
da[,1]=factor(da[,1],c("UCS","OV","LUSC","ESCA","READ",'HNSC','PAAD','COAD','LUAD','STAD','BLCA','LGG','UCEC','SARC',"BRCA",'GBM',"LIHC",'SKCM',"PRAD",'LAML','CESC','THYM','KIRC','KIRP'))
# 当width < 1 时饼图将变成饼环  ## width >= 1 时中心的杂点将消失#direction=-1
p <- ggplot(da,aes(x=factor(1),y=da[,6],fill=state,colour=HR_state))+
  geom_bar(stat="identity",width=1)+coord_polar(theta = "y")+facet_grid(da[,2]+da[,7]~da[,1])+
  scale_color_manual(values=c("red","green"))+
  geom_point(aes(shape=Pvalue),size =0.6) 
p <- p+theme(legend.position="right",legend.box="vertical",legend.title = element_text(size = 5),legend.text = element_text(size = 3,face = "bold"),plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "white"),axis.text=element_text(size=3,colour="black"),strip.text=element_text(size=3,colour="black",angle=45))+
  theme(axis.text.y = element_blank())+theme(axis.text.x = element_blank(),axis.ticks=element_line(colour = "white"),strip.background=element_rect(fill ="white", colour = "grey"),panel.spacing=unit(0,"lines"))
p <- p + xlab("Gene") + ylab("Cancer") + ggtitle("RiskAnalysis")
ggsave(file=paste(path,"1.pdf",sep=""),width = 5,height=25)

 
 
 


###############################################
install.packages("corrplot")
require(corrplot)
#method参数主要包括以下几种："circle", "square", "ellipse", "number", "shade", "color", "pie"
corrplot(corr,method="shade",shade.col=NA,tl.col="black",tl.srt=45,is.corr=FALSE)
corrplot.mixed(cor(mtcars))
#混合方法之上三角为圆形，下三角为方形
corrplot(corr = corr,order="AOE",type="upper",tl.pos="d")
corrplot(corr = corr,add=TRUE, type="lower", method="square",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n")

###############################################
##做p53上下文基因对药物的热图#对所有细胞系16cancer heatmap1
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/clinical_gene_drug/"
f<-read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f[,1]=as.character(f[,1])

ff<-read.table("D:/tp53/Gene Expression/ENSG_Name2.txt",header=F,sep="\t")
ff[,1]=gsub("\\..*$","",as.character(ff[,1]))
ff[,2]=gsub(" ", "",as.character(ff[,2]))
m=data.frame(table(ff[,2]))
m1=ff[ff$V2 %in% m[m[,2]==1,1],]
rownames(m1)<-m1[,2]
m2=na.omit(m1[f[,1],])

f1<-read.csv(paste(path,"drug_xx/Cell_Lines_Details.txt",sep=""),header=T,sep="\t") 
f2<-read.csv(paste(path,"drug_xx/dose_response.txt",sep=""),header=T,sep="\t")
a=unique(f2[,c("COSMIC_ID","CELL_LINE_NAME")])
f3<-read.csv(paste(path,"drug_xx/sanger1018_brainarray_ensemblgene_rma.txt",sep=""),header=T,sep="\t")
d=na.omit(merge(m2,f3,by.x="V1",by.y="ensembl_gene",all.x=T))
colnames(d)=gsub("X","",colnames(d))
rownames(d)=d$V2
d1=data.frame(t(d[,-(1:2)]))
d1$id=rownames(d1)
d2=na.omit(merge(d1,a,by.x="id",by.y="COSMIC_ID",all.x=T))
rownames(d2)=d2$CELL_LINE_NAME
d3=data.frame(t(d2[,-c(1,193)]))
colnames(d3)=gsub("\\.","",colnames(d3))
#21-2-3=16cancer PRAD、ACC不cluster_col 去LIHC\LUSC\UCEC只有2个细胞系
#cancer=c("HNSC","ESCA","LAML","BRCA","COAD/READ","LIHC","STAD","PRAD","OV","UCEC","CESC","BLCA","THCA","SKCM","PAAD","LGG","GBM","LUSC","LUAD","KIRC","ACC")
for(cancer in c("HNSC","ESCA","LAML","BRCA","COAD/READ","STAD","OV","CESC","BLCA","THCA","SKCM","PAAD","LGG","GBM","LUAD","KIRC")){
  b1=gsub("-","",as.character(f1[as.character(f1$Cancer.Type..matching.TCGA.label.)==cancer,"Sample.Name"]))
  b2=f2[as.character(f2$CELL_LINE_NAME) %in% b1,c("DRUG_NAME","CELL_LINE_NAME","AUC")]
  da=data.frame(matrix(NA,length(unique(b2[,1])),nrow(d3))) ##drug_gene cor相关性矩阵
  rownames(da)=unique(b2[,1])
  colnames(da)=rownames(d3)
  #daP=data.frame(matrix(-2,length(unique(b2[,1])),nrow(d3)))
  #rownames(daP)=unique(b2[,1])
  #colnames(daP)=rownames(d3)
  for(i in 1:nrow(da)){
    for(j in 1:ncol(da)){
      n1=b2[b2$DRUG_NAME==rownames(da)[i],3]
      n2=d3[colnames(da)[j],intersect(colnames(d3),as.character(b2[b2$DRUG_NAME==rownames(da)[i],2]))]
      if((length(n1)!=1) & (length(as.numeric(n2))!=1) & (length(n1)=length(as.numeric(n2)))){
        da[i,j]=cor(n1,as.numeric(n2))
        #daP[i,j]=cor.test(n1,as.numeric(n2))$p.value
      }
    }
  }
  da1=t(na.omit(da))
  da1=da1[intersect(as.character(f[,1]),as.character(rownames(da1))),]
  y<-as.matrix(da1)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(-1,-0.5,0.5,1),breaks=c(-1,-0.5,0.5,1),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"drug_xx/heatmap1/",gsub("/","_",cancer),".pdf",sep=""),width=8,height=8,main=cancer)  
}


###把细胞系分为tp53突变和没突变的，然后做p53上下文基因对药物的热图
#heatmap2/1/  9cancer有mu\nomu cor>0.5为红
#heatmap2/2/  9cancer有mu\nomu cor>0.5且p<0.05为红number
library(pheatmap)
colorsChoice<- colorRampPalette(c("green","white","red"))
path="D:/tp53/clinical_gene_drug/"
f<-read.table("D:/tp53/Gene Expression/p53_context1.txt",header=T,sep="\t")
f[,1]=as.character(f[,1])

ff<-read.table("D:/tp53/Gene Expression/ENSG_Name2.txt",header=F,sep="\t")
ff[,1]=gsub("\\..*$","",as.character(ff[,1]))
ff[,2]=gsub(" ", "",as.character(ff[,2]))
m=data.frame(table(ff[,2]))
m1=ff[ff$V2 %in% m[m[,2]==1,1],]
rownames(m1)<-m1[,2]
m2=na.omit(m1[f[,1],])

f1<-read.csv(paste(path,"drug_xx/Cell_Lines_Details.txt",sep=""),header=T,sep="\t") 
f2<-read.csv(paste(path,"drug_xx/dose_response.txt",sep=""),header=T,sep="\t")
a=unique(f2[,c("COSMIC_ID","CELL_LINE_NAME")])
f3<-read.csv(paste(path,"drug_xx/sanger1018_brainarray_ensemblgene_rma.txt",sep=""),header=T,sep="\t")
d=na.omit(merge(m2,f3,by.x="V1",by.y="ensembl_gene",all.x=T))
colnames(d)=gsub("X","",colnames(d))
rownames(d)=d$V2
d1=data.frame(t(d[,-(1:2)]))
d1$id=rownames(d1)
d2=na.omit(merge(d1,a,by.x="id",by.y="COSMIC_ID",all.x=T))
rownames(d2)=d2$CELL_LINE_NAME
d3=data.frame(t(d2[,-c(1,193)]))
colnames(d3)=gsub("\\.","",colnames(d3))

f4=read.csv(paste(path,"drug_xx/TP53_mutations.txt",sep=""),header=T,sep="\t") ##所有发生tp53突变的细胞系
d4=gsub("_.*","",as.character(f4$Tumor.Sample.Barcode))
##"ESCA","LAML","THCA","PAAD","LGG"的nomu不cluster_col、"CESC","LUAD"的mu不cluster_col
#cancer=c("HNSC","BRCA","COAD/READ","STAD","OV","BLCA","SKCM","GBM","KIRC")##9cancer有mu\nomu
#3去BLCA
number=3
for(cancer in c("SKCM","GBM","KIRC")){
  b1=gsub("-","",as.character(f1[as.character(f1$Cancer.Type..matching.TCGA.label.)==cancer,"Sample.Name"]))
  mu=intersect(d4,b1)
  nomu=setdiff(b1,mu)
  b2=f2[as.character(f2$CELL_LINE_NAME) %in% mu,c("DRUG_NAME","CELL_LINE_NAME","AUC")]
  b3=f2[as.character(f2$CELL_LINE_NAME) %in% nomu,c("DRUG_NAME","CELL_LINE_NAME","AUC")]
  da=data.frame(matrix(NA,length(unique(b2[,1])),nrow(d3))) ##drug_gene cor相关性矩阵mu
  rownames(da)=unique(b2[,1])
  colnames(da)=rownames(d3)
  dano=data.frame(matrix(NA,length(unique(b3[,1])),nrow(d3))) ##drug_gene cor相关性矩阵nomu
  rownames(dano)=unique(b3[,1])
  colnames(dano)=rownames(d3)
  for(i in 1:nrow(da)){
    for(j in 1:ncol(da)){
      n1=b2[b2$DRUG_NAME==rownames(da)[i],3]
      n2=d3[colnames(da)[j],intersect(colnames(d3),as.character(b2[b2$DRUG_NAME==rownames(da)[i],2]))]
      ##细胞系个数大于3
      if((length(n1)>=number) & (length(as.numeric(n2))>=number) & (length(n1)=length(as.numeric(n2)))){
        da[i,j]=cor(n1,as.numeric(n2))
      }
    }
  }
  for(ii in 1:nrow(dano)){
    for(jj in 1:ncol(dano)){
      n11=b3[b3$DRUG_NAME==rownames(dano)[ii],3]
      n22=d3[colnames(dano)[jj],intersect(colnames(d3),as.character(b3[b3$DRUG_NAME==rownames(dano)[ii],2]))]
      if((length(n11)>=number) & (length(as.numeric(n22))>=number) & (length(n11)=length(as.numeric(n22)))){
        dano[ii,jj]=cor(n11,as.numeric(n22))
      }
    }
  }
  da1=t(na.omit(da))
  da1=da1[intersect(as.character(f[,1]),as.character(rownames(da1))),]
  y<-as.matrix(da1)
  row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(-1,-0.5,0.5,1),breaks=c(-1,-0.5,0.5,1),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"drug_xx/heatmap2/2/",gsub("/","_",cancer),"_mu_",number,".pdf",sep=""),width=8,height=8,main=paste(cancer,"mu",sep="_"))  
  dano1=t(na.omit(dano))
  dano1=dano1[intersect(as.character(f[,1]),as.character(rownames(dano1))),]
  yno<-as.matrix(dano1)
  pheatmap(yno,cluster_row = FALSE,cluster_col = TRUE,col=colorsChoice(3),legend_breaks=c(-1,-0.5,0.5,1),breaks=c(-1,-0.5,0.5,1),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
           file=paste(path,"drug_xx/heatmap2/2/",gsub("/","_",cancer),"_nomu_",number,".pdf",sep=""),width=8,height=8,main=paste(cancer,"nomu",sep="_"))  
}

 
#######################
f<-read.csv("D:/tp53/GO/goa_human.txt",header=F,sep="\t")
f[,3]=as.character(f[,3])
f[,5]=as.character(f[,5])
f[,9]=as.character(f[,9])
a=unique(f[,c(3,5,9)])
a1=data.frame(table(a[,2]))
a2=unique(a[,2:3])
da=merge(a1,a2,by.x="Var1",by.y="V5") 
for(i in 1:nrow(da)){
  b={}
  for(j in a[a$V5==as.character(da[i,1]),1][1:length(a[a$V5==as.character(da[i,1]),1])]){
    b=paste(b,j,sep=",")
  }
  da[i,4]=b
}
f1<-read.csv("D:/tp53/GO/goose",header=F,sep="\t")
f1[,2]=as.character(f1[,2])
f1[,4]=as.character(f1[,4]) 
f1[,9]=apply(f1[,5:8],1,max)
a3=unique(f1[,c(2,4)])
a4=data.frame(GOid=unique(f1[,4]))
for(i in unique(f1[,4])){
  a4[a4$GOid==i,2]=max(f1[f1$V4==i,9])
}
a4[,1]=as.character(a4[,1])
da1=unique(merge(da,a3,by.x="Var1",by.y="V4",all.x=T))
da2=merge(da1,a4,by.x="Var1",by.y="GOid",all.x=T)
write.table(da2,file="D:/tp53/GO/GOid.txt",row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)

##
f<-read.csv("D:/tp53/GO/GOid.txt",header=T,sep="\t") 
f[,3]=gsub("P","BP",as.character(f[,3]))
f[,3]=gsub("C","CC",as.character(f[,3]))
f[,3]=gsub("F","MF",as.character(f[,3]))
a=na.omit(unique(f[(f[,2]>=5) & (f[,2]<=100) & (f[,6]>=7),]))
write.table(a,file="D:/tp53/GO/GOid1.txt",row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)
##
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(GO.db)
#ancestor <- as.list(GOBPANCESTOR)
#ancestor <- toTable(GOBPANCESTOR)
#ancestor <- GOBPANCESTOR$"GO:0000012"
#ancestor1 <- GOBPPARENTS$"GO:0000012"
BP <- toTable(GOBPANCESTOR)
BP1 <- toTable(GOBPPARENTS)
MF <- toTable(GOMFANCESTOR)
MF1 <- toTable(GOMFPARENTS)
CC <- toTable(GOCCANCESTOR)
CC1 <- toTable(GOCCPARENTS)
f<-read.csv("D:/tp53/GO/GOid1.txt",header=T,sep="\t") 
a1=as.character(f[as.character(f[,3])=="BP",1])
a2=as.character(f[as.character(f[,3])=="MF",1])
a3=as.character(f[as.character(f[,3])=="CC",1])
b1=setdiff(a1,unique(c(BP[,2],BP1[,2])))
b2=setdiff(a2,unique(c(MF[,2],MF1[,2])))
b3=setdiff(a3,unique(c(CC[,2],CC1[,2])))
da=c(b1,b2,b3)
da1=f[as.character(f$Var1) %in% da,]
write.table(da1,file="D:/tp53/GO/GOid2.txt",row.names = F,col.names = T, quote = F,sep = "\t",append=TRUE)



#####另加图
library(ggplot2)
path="D:/tp53/phenotype_Analysis/"
f=data.frame(Project=c("LGG","GBM","READ","COAD","KIRC","KIRP","LUSC","LUAD","UCS","UCEC"),
             PrimarySite=c("Brain","Brain","Colorectal","Colorectal","Kidney","Kidney","Lung","Lung","Uterus","Uterus"),
             Mu=c(509,393,137,399,336,281,492,567,57,530),
             TP53Mu=c(243,123,110,227,9,5,427,298,53,212)) 
f[,5]=f[,4]/f[,3]
f[,1]=factor(f[,1],c("LGG","GBM","READ","COAD","KIRC","KIRP","LUSC","LUAD","UCS","UCEC"))
p<-ggplot(f, aes(x = f[,1], y = f[,5],fill=f$PrimarySite))+geom_bar(width=0.3,alpha=0.7,stat="identity",position="dodge")+
  scale_fill_manual(values=c("coral4","darkgoldenrod4","forestgreen","dodgerblue3","darkorchid3"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black",angle=90),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Cancer") + ylab("TP53 Mutation Ratio") + ggtitle("In Same Primary Site")  
ggsave(file=paste(path,"SamePrimarySite.pdf",sep=""),width = 10,height=9)  



library(openxlsx) 
ff=read.xlsx(paste(path,"subtype/subtypecount2.xlsx",sep=""),sheet=1) 
ff[4:5,2]=gsub("GI","COAD",ff[4:5,2])
ff[6:7,2]=gsub("_LGG","",ff[6:7,2])
ff[8:12,2]=gsub("GBM_","",ff[8:12,2])
ff[,9]=gsub("\\..*","",ff$`#Cancer.Subtype`)
ff[,2]=as.factor(ff[,2])
ff[,9]=as.factor(ff[,9])
p<-ggplot(ff, aes(x = ff[,2], y = ff[,5],fill=ff$V9))+geom_bar(width=0.3,alpha=0.7,stat="identity",position="dodge")+
  scale_fill_manual(values=c("cadetblue4","darkgoldenrod4","coral4","coral","darkorchid3"))
p <- p+theme(legend.position="top",plot.title=element_text(hjust = 0.5),panel.border=element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = NULL,axis.text=element_text(size=8,colour="black",angle=90),strip.text=element_text(size=8,colour="black"))
p <- p + xlab("Cancer Subtype") + ylab("TP53 Mutation Ratio") + ggtitle("In Same Cancer")  
ggsave(file=paste(path,"SameCancer.pdf",sep=""),width = 10,height=9)  



##################################################
##

f<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",header=F,sep="\t") 
f1<-read.csv("D:/tp53/Gene Expression/Gene_FPKM1_Analysis2/p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D_1075.txt",header=F,sep="\t") 
a=as.character(f[f[,2]==4,1])
a1=as.character(f1[f1[,2]==4,1])
length(intersect(a,a1))



##################
##根据p53_context_discreteMode对All癌症做热图+注释
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
  #row_anno = data.frame(GeneFunction=as.character(f[,2]), row.names=as.character(f[,1]))
  #pheatmap(y,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="euclidean",clustering_method = "ward.D",
  #col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,fontsize_row=5,fontsize_col=3,
  #file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",gsub(".txt",".pdf",cancer),sep=""),width=8,height=8,main=gsub(".txt","",cancer))  
  a=f1[intersect(as.character(f[,1]),as.character(f1[,1])),]
  h=merge(h,a,by.x="name",by.y="name",all.x=T)
}
rownames(h)=h[,1]
h1=na.omit(h[intersect(as.character(f[,1]),as.character(h[,1])),-1])
number=2
fu<-function(x) sum(x!=0) ##删去0比较多的基因337-202-120\112\101
a1=data.frame(apply(h1,1,fu))
a2=rownames(a1)[which(a1[,1]>=number)] ##10\15\20
h2=na.omit(h[a2,-1])
#write.table(h2,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",number,"/Allcancer",number,"_manhattan_ward.D_3009.txt",sep=""),row.names = T,col.names = T, quote = F,sep = "\t",append=TRUE)

f2=read.table(paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/2/Allcancer2_class4_manhattan_ward.D.txt",sep=""),header=F,sep="\t")
f2[,2]=paste("C",f2[,2],sep="")
f2[,1]=as.character(f2[,1])
ff=read.table("D:/tp53/Gene Expression/Gene_FPKM1_Analysis1/Result2.txt",header=T,sep="\t")
ff[,1]=gsub("-",".",as.character(ff[,1]))
f22=unique(merge(f2,ff[,c(1,2)],by.x="V1",by.y="submitter_id",all.x=T))
f22[,3]=as.character(f22[,3])

y2<-as.matrix(h2) 
row_anno = data.frame(Function=as.character(f[,2]), row.names=as.character(f[,1]))
col_anno2 = data.frame(Class=as.character(f22[,2]),Cancer=as.character(f22[,3]),row.names=as.character(f22[,1]))
ann_colors2 = list(Class = c(C1 ="deeppink3",C2 ="goldenrod3",C3 ="palegreen4",C4 ="royalblue4"),
                   Cancer = c(GBM ="deeppink3",
                              COAD ="goldenrod3",READ ="goldenrod2",
                              KIRC ="royalblue4",KIRP ="royalblue3",
                              LUSC ="palegreen4",LUAD ="palegreen3",
                              UCEC ="grey",
                              ESCA ="blueviolet",BLCA ="brown3",
                              PAAD ="chocolate2",LIHC ="coral4",SARC ="cornflowerblue",
                              THYM ="bisque4",SKCM ="chartreuse3",BRCA ="cadetblue3",HNSC ="brown1",
                              STAD ="darkslateblue",CESC ="bisque2",PRAD ="lightcoral",THCA="black"))
pheatmap(y2,cluster_row = FALSE,cluster_col = TRUE,clustering_distance_cols="manhattan",clustering_method = "ward.D",
         col=colorsChoice(7),legend_breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),breaks=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5),annotation_row = row_anno,annotation_col = col_anno2,fontsize_row=3,show_colnames=FALSE,
         fontsize=6,annotation_colors = ann_colors2,file=paste(path,"p53_context_discreteMode_fc1.5_heatmap/123/",number,"/Allcancer",number,"_manhattan_ward.D_1.pdf",sep=""),width=8,height=8,main="Pancancer")  

