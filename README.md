# test
step123-test
#===========================library==================================
library(tidyverse)
library(VIM)
library(pheatmap)
library(scales)
library(mice)
library(pcaMethods)

library(missForest)
library(VennDiagram)
library(reshape2)
library(MASS)
library(impute)
#==========================mac data dir==============================
input_dir <- "F:/F16FTSNWKF1186_nongkeyuan_serum_18_mice/F16FTSNWKF1186_nongkeyuan_serum_18_mice/CSV"
meaf <- paste(input_dir,"20160622_mics_18_measurements_neg.csv",sep="/")
listf <- paste(input_dir,"20160622_mice_18_neg_spl.txt",sep="/")
out_dir <- "F:/F16FTSNWKF1186_nongkeyuan_serum_18_mice/result_w"
#==========================commonlly used functions===================
ratio_count <- function(x){
  sum(is.na(x))/length(x)
}
#==========================intata===================

  #========================read in file==================================================
  list <- read.table(listf,header = T,sep="\t")
  list$class <- as.character(list$class)
  original_data <- read_csv(meaf)
  #==========================================merge list and peak data================================================
  list$mark <- "sample"
  list$mark[list$class=="QC"] <- "QC"
  original_data[original_data==0] <- NA
  original_datat <- as.data.frame(t(original_data[,-c(1,2)]))
  original_datat[c(1:5),c(1:10)]
  names(original_datat) <- original_data$mz
  original_datat$sample <- row.names(original_datat)
  original_tl <- merge(list,original_datat,by.x = "sample",by.y = "sample")
  original_tl[c(1:5),c(1:15)]
  #========================================return==========================================
  setClass("Indata",
           slots = list(
             listf = "data.frame",
             peak = "data.frame",
             merged="data.frame",
             filter_by_group="data.frame",
             traditional_filter="data.frame",
             total_group_miss_imp="data.frame",
             impute_low="data.frame",
             current_dataset="data.frame",
             mod_filter="data.frame"
           )
  )
  indatas <- new("Indata", listf = list, peak=original_data, merged=original_tl,current_dataset=original_tl)

  
  
  #==============================percent of missing in sample mz and in total================

    data <- indatas@current_dataset
    df <- data.frame(
      group = c("no_miss","has_miss"),
      value = c("FALSE", "TRUE")
    )
    pdf(paste(out_dir,"miss_percent.pdf",sep="/"),width = 6,height = 3)
    op=par(mfrow=c(1,2))
    #mz miss percent
    mz_miss <- table(sapply(data[,-c(1,2,3,4,5)], function(x){any(is.na(x))})) %>%
      as.data.frame() %>%
      merge(df,by.x="Var1",by.y="value")
    
    pie(mz_miss$Freq,mz_miss$group,main = "mz contain missing", col = rainbow(length(mz_miss$Freq)))
    #pie(mz_miss$Freq,mz_miss$group, col = rainbow(length(mz_miss$Freq)))
    
    #percent of missing data
    miss_compose <- data.frame(
      group = c("no_miss","has_miss"),
      Freq = c(sum(!(is.na(data[,-c(1,2,3,4,5)]))),sum(is.na(data[,-c(1,2,3,4,5)])))
    )
    pie(miss_compose$Freq,miss_compose$group,main = "total missing data percent", col = rainbow(length(miss_compose$Freq)))
    #pie(miss_compose$Freq,miss_compose$group, col = rainbow(length(miss_compose$Freq)))
    par(op)
    dev.off()
    #===================================sample wise missing pattern and distribution==============================

      data <- indatas@current_dataset
      list <- indatas@listf
      sample_name <- data$sample
      data <- as.data.frame(t(data[,-c(1,2,3,4,5)]))
      names(data) <- sample_name
      md_pattern <- md.pattern(data)
      #missing percent in each sample boxplot
      sample_miss_number <- md_pattern[nrow(md_pattern),][-ncol(md_pattern)]/nrow(data)
      sample_miss_number_plot <- data.frame(
        sample=names(sample_miss_number),
        miss_ratio=sample_miss_number) %>%
        merge(list[,c(1,3)],by.x="sample",by.y="sample") %>%
        ggplot(aes(class,miss_ratio))+
        geom_boxplot()+
        geom_jitter(width = 0.2)+
        theme_light()
      pdf(paste(out_dir,"sample_percent.pdf",sep="/"),width = 3,height = 3)
      print(sample_miss_number_plot)
      dev.off()
    

      
      
      #=================================mz wise distribution=========================

        data <- indatas@current_dataset
        #correlation between missing ratio of a mz and the mz intensity
        mean_intensity <- apply(data[,-c(1,2,3,4,5)],2,function(x){mean(as.numeric(x),na.rm = T)})  #head(mean_intensity)
        miss_ratio <- apply(data[,-c(1,2,3,4,5)],2,function(x){sum(is.na(x))/length(x)})  #head(miss_ratio)
        intensity_miss_ratio <- data.frame(mean_intensity,miss_ratio) #head(intensity_miss_ratio)
        pdf(paste(out_dir,"mz_miss_correlation_with_intensity.pdf",sep="/"),width=3,height=2)
        p <- ggplot(subset(intensity_miss_ratio,miss_ratio>0),aes(mean_intensity,miss_ratio))+geom_point()+theme_light()
        print(p)
        dev.off()
        intensity_miss_ratio$mz <- as.numeric(row.names(intensity_miss_ratio))
        pdf(paste(out_dir,"mz_miss_correlation_with_mz_value.pdf",sep="/"),width = 3,height = 2)
        p <- ggplot(subset(intensity_miss_ratio,miss_ratio>0),aes(mz,miss_ratio))+geom_point()+theme_light()
        print(p)
        dev.off()
        #correlation between missing ratio of a mz and the RT value
        rt_intensity <- merge(intensity_miss_ratio,indatas@peak[,c(1,2)],by.x="mz",by.y="mz")
        rt_intensity$RT <- as.numeric(as.character(rt_intensity$RT))
        pdf(paste(out_dir,"mz_miss_correlation_with_rt_value.pdf",sep="/"),width = 3,height = 2)
        p <- ggplot(subset(rt_intensity,miss_ratio>0),aes(RT,miss_ratio))+geom_point()+theme_light()
        print(p)
        dev.off()
        #missing ratio of mz distribution in QC and samples
        miss_ratio <- aggregate(data[,-c(1:5)],by=list(data$mark),ratio_count)
        miss_ratiot <- as_tibble(t(miss_ratio[,-1]))
        names(miss_ratiot) <- as.character(miss_ratio[,1])
        pdf(paste(out_dir,"mz_miss_ratio_distribution.pdf",sep="/"),width = 3,height = 2)
        p <- miss_ratiot %>%
          gather(`QC`,`sample`,key="mark",value="value") %>%
          filter(value>0) %>%
          ggplot(aes(value,fill=mark)) +
          geom_histogram(position = "dodge",bins=10)+
          theme_light()+
          theme(legend.position="top")
        print(p)
        dev.off()
        
        
        #=================================out liear=====================================

          data <- as.matrix(indatas@current_dataset[,-c(1,2,3,4,5)])
          sgroup <- indatas@current_dataset[,3]
          rownames(data) <- indatas@current_dataset$sample
          pc <- pca(data, nPcs=2, method="ppca")
          #slplot(pc, sl=NULL, spch=5)
          plotd <- cbind(sgroup,as.data.frame(pc@scores))
          pdf(paste(out_dir,"outlier_discovery.pdf",sep="/"),width = 3,height = 2)
          p <- ggplot(plotd,aes(PC1,PC2))+
            geom_point(aes(colour=sgroup))+
            stat_ellipse(aes(x=PC1, y=PC2,color=sgroup),type = "norm")+
            theme_light()+
            theme(legend.position="top")
          print(p)
          dev.off()
          
          
          
 #============================correlation analysis===============================       
          


         
            
            
            ###====================miss filter based on group================================

            data <- indatas@current_dataset
            #data <- indata@merged
            #qc_ratio <- 0.5
            #sam_ratio <- 0.5
            miss_ratio_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),ratio_count)   #miss_ratio_by_class[,c(1:5)]
            qc_miss_ratio <- as.data.frame(t(miss_ratio_by_class[miss_ratio_by_class$Group.1=="QC",-1]) >= 0.5)    #head(qc_miss_ratio)
            names(qc_miss_ratio) <- c("qc")
            qc_fliter_mz <- row.names(subset(qc_miss_ratio,qc=="TRUE"))
            
            sam_miss_ratio <- miss_ratio_by_class[!(miss_ratio_by_class$Group.1=="QC"),]  #sam_miss_ratio[,c(1:5)]
            sam_miss_judge <- apply(sam_miss_ratio[,-1],2,function(x){all(x>=0.5)})
            sample_fliter_mz <- names(which(sam_miss_judge=="TRUE"))
            
            mz_fliter <- union(qc_fliter_mz,sample_fliter_mz)   #length(mz_fliter)
            data_flitered <- data[,!(names(data) %in% mz_fliter)]   #data_flitered[c(1:5),c(1:10)]
            indatas@filter_by_group <- data_flitered
            indatas@current_dataset <- data_flitered
 
            
              ##=======================traditional miss fliter==============================

                data <- indatas@current_dataset
                #data <- indatas@merged
                #qc_ratio <- 0.5
                #sam_ratio <- 0.5
                
                
                miss_ratio_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),ratio_count)   #miss_ratio_by_class[,c(1:5)]
                sam_miss_ratio <- miss_ratio_by_class[!(miss_ratio_by_class$Group.1=="QC"),]  #sam_miss_ratio[,c(1:5)]
                sam_miss_judge <- apply(sam_miss_ratio[,-1],2,function(x){(any(x==1) & any(x<0.5))})
                sample_left_mz <- names(which(sam_miss_judge=="TRUE"))       #any mz with one total group missing and another group with less than 50% missing
                
                miss_ratio_by_type <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$mark),ratio_count)   #miss_ratio_by_type[,c(1:10)]
                miss_ratio_by_type2 <- as.data.frame(t(miss_ratio_by_type[,-1]))
                names(miss_ratio_by_type2) <- miss_ratio_by_type$Group.1   
                mz_filter_qc_t <- row.names(subset(miss_ratio_by_type2,QC >= 0.5))   #remove mz based on qc
                mz_filter_sa_t <- row.names(subset(miss_ratio_by_type2,sample >= 0.5))  #remove mz based on sample
                mz_tra_filter<-union(mz_filter_qc_t,mz_filter_sa_t)   #traditional missing filtering
                
                data_flitered<- data[,!(names(data) %in% mz_tra_filter)]
                indatas@traditional_filter <- data_flitered
                indatas@current_dataset <- data_flitered
           
           
              
              ##=======================mod traditional miss fliter============================== 
              
              data <- indatas@current_dataset
              #data <- indatas@merged
              #qc_ratio <- 0.5
              #sam_ratio <- 0.5
              
              
              miss_ratio_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),ratio_count)   #miss_ratio_by_class[,c(1:5)]
              sam_miss_ratio <- miss_ratio_by_class[!(miss_ratio_by_class$Group.1=="QC"),]  #sam_miss_ratio[,c(1:5)]
              sam_miss_judge <- apply(sam_miss_ratio[,-1],2,function(x){(any(x==1) & any(x<0.5))})
              sample_left_mz <- names(which(sam_miss_judge=="TRUE"))       #any mz with one total group missing and another group with less than 50% missing
              
              miss_ratio_by_type <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$mark),ratio_count)   #miss_ratio_by_type[,c(1:10)]
              miss_ratio_by_type2 <- as.data.frame(t(miss_ratio_by_type[,-1]))
              names(miss_ratio_by_type2) <- miss_ratio_by_type$Group.1   
              mz_filter_qc_t <- row.names(subset(miss_ratio_by_type2,QC >= 0.5))   #remove mz based on qc
              mz_filter_sa_t <- row.names(subset(miss_ratio_by_type2,sample >= 0.5))  #remove mz based on sample
              mz_mod_filter<- union(mz_filter_qc_t,setdiff(mz_filter_sa_t,sample_left_mz))  #modified filtering leaves total missing mz
              data_flitered<- data[,!(names(data) %in% mz_mod_filter)]
              indatas@mod_filter <- data_flitered
              indatas@current_dataset <- data_flitered
              
            
              
              
              
              
              
#==============================comparison of the two filtering strategy============================

              
  #===================================== compare_group_filter_traditional_filter==========================================================
            
group_filtered_data <- indatas@filter_by_group   #group_filtered_data[c(1:5),c(1:10)]
traditional_filtered_data <- indatas@traditional_filter     #traditional_filtered_data[c(1:5),c(1:10)]
original <- indatas@merged
group_filter_left_mz <- setdiff(names(group_filtered_data),names(traditional_filtered_data))  #original[c(1:5),c(1:10)]
if(length(group_filter_left_mz)>0){
  group_filter_left <- original[,c(3,which(names(original) %in% group_filter_left_mz))]  #group_filter_left[c(1:5),]
  group_filter_left2 <- t(group_filter_left[,-1])   #group_filter_left2[c(1:5),c(1:5)]
  colnames(group_filter_left2) <- original$sample
  list<-list[order(list$class),]
  sorted_order_matrixplot <- sapply(list$sample,function(x){which(colnames(group_filter_left2)==x)})
  matrixp <- group_filter_left2[,sorted_order_matrixplot]
  colnames(matrixp) <- list$class
  pdf(paste(out_dir,"comparison_tra_with_group.pdf",sep="/"),width = 6,height = 4)
  matrixplot(matrixp)
  dev.off()
  }else{cat("no difference")
        }

              
              
  #==================================compare_mod_filter_traditional_filter==============================================================
mod_filtered_data <- indatas@mod_filter   
traditional_filtered_data <- indatas@traditional_filter     #traditional_filtered_data[c(1:5),c(1:10)]
original <- indatas@merged
mod_filter_left_mz <- setdiff(names(mod_filtered_data),names(traditional_filtered_data))  #original[c(1:5),c(1:10)]
if(length(mod_filter_left_mz)>0){
  mod_filter_left <- original[,c(3,which(names(original) %in% mod_filter_left_mz))]  #group_filter_left[c(1:5),]
  mod_filter_left2 <- t(mod_filter_left[,-1])   #group_filter_left2[c(1:5),c(1:5)]
  colnames(mod_filter_left2) <- original$sample
  list<-list[order(list$class),]
  sorted_order_matrixplot <- sapply(list$sample,function(x){which(colnames(mod_filter_left2)==x)})
  matrixp <- mod_filter_left2[,sorted_order_matrixplot]
  colnames(matrixp) <- list$class
  pdf(paste(out_dir,"comparison_mod_with_tra.pdf",sep="/"),width = 6,height = 4)
  matrixplot(matrixp)
  dev.off()
}else{cat("no difference")
      }

              
              
#==================================compare_mod_filter_filter_bygroup==============================================================
mod_filtered_data <- indatas@mod_filter   
group_filtered_data <- indatas@filter_by_group    
original <- indatas@merged
group_filter_left_mz <- setdiff(names(group_filtered_data),names(traditional_filtered_data))  
if(length(group_filter_left_mz)>0){
  group_filter_left <- original[,c(3,which(names(original) %in% group_filter_left_mz))]  
  group_filter_left2 <- t(group_filter_left[,-1])  
  colnames(group_filter_left2) <- original$sample
  list<-list[order(list$class),]
  sorted_order_matrixplot <- sapply(list$sample,function(x){which(colnames(group_filter_left2)==x)})
  matrixp <- group_filter_left2[,sorted_order_matrixplot]
  colnames(matrixp) <- list$class
  pdf(paste(out_dir,"comparison_mod_with_group.pdf",sep="/"),width = 6,height = 4)
  matrixplot(matrixp)
  dev.off()
}else{cat("no difference")
      }
venn.diagram(list(A=mz_tra_filter,B=mz_mod_filter,C=mz_fliter), fill=c("red","green","blue"), alpha=c(0.5,0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename=paste(out_dir,"VennDiagram.tiff",sep="/"))         

              


#=====================delete_miss=======================================
miss_del_data<- original_data[complete.cases(original_data),] 
miss_del_datat <- as.data.frame(t(miss_del_data[,-c(1,2)]))
names(miss_del_datat) <- miss_del_data$mz
miss_del_datat$sample <- row.names(miss_del_datat)
miss_del_datatl <- merge(list,miss_del_datat,by.x = "sample",by.y = "sample")
write.csv(miss_del_datatl,paste(input_dir,"miss_del_data",sep="/"),row.names = FALSE)
miss_del_data_D <- read.csv(paste(input_dir,"miss_del_data_D.csv",sep="/"))





  #===================misforset===============================
    data <- miss_del_data_D[,-c(1,2,3,4,5)]
    data2 <- missForest(data)$ximp
    data2[data2<0] <- 0.01
    missforest_imp<- cbind(miss_del_data_D[,c(1,2,3,4,5)],data2)

  #==================knn======================================
    data <- miss_del_data_D[,-c(1,2,3,4,5)]
    data2 <- impute.knn(as.matrix(data))
    data3 <- data2$data
    data3[data3<0] <- 0.01
    knn_impt <- cbind(miss_del_data_D[,c(1,2,3,4,5)],data3)

  #========================ppca===============================
    data <- miss_del_data_D[,-c(1,2,3,4,5)]
    row.names(data) <- miss_del_data_D$sample
    data2 <- t(data)
    pc <- pca(data2,nPcs=2,method="ppca")
    imputed <- completeObs(pc)
    imputed[imputed<0] <- 0.01
    ppca_imp <- cbind(miss_del_data_D[,c(1,2,3,4,5)],as.data.frame(t(imputed)))
  
#===================================RMSE=======================================
actual<- miss_del_datatl[,-c(1,2,3,4,5)]
predict<- missforest_imp[,-c(1,2,3,4,5)]
RMSE_missforcet = (mean((predict - actual)^2))^0.5


actual<- miss_del_datatl[,-c(1,2,3,4,5)]
predict<- knn_impt[,-c(1,2,3,4,5)]
RMSE_knn = (mean((predict - actual)^2))^0.5   
    
    
actual<- miss_del_datatl[,-c(1,2,3,4,5)]
predict<- ppca_imp[,-c(1,2,3,4,5)]
RMSE_ppca = (mean((predict - actual)^2))^0.5


#==============================min_knn_impute=========================================================

#==============================min_impute=============================================================
#==============================total group missed handling==================================

data <- miss_del_data_D
all_value <- as.numeric(as.matrix(data[,-c(1,2,3,4,5)]))
low_value <- quantile(all_value,0.01,na.rm = T)
all_low <-all_value[which(all_value<low_value)]
mean_l <- mean(all_low,na.rm = T)
sd_l <- sd(all_low,na.rm = T)
miss_ratio_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),ratio_count) 
total_group_missed_mz <- miss_ratio_by_class %>%  
  melt(id=c("Group.1")) %>%
  filter(value == 1) 

min <- min(as.data.frame(miss_del_data_D[,-c(1,2,3,4,5)]),na.rm=TRUE)  

data2 <- as.matrix(data[,-c(1,2,3,4,5)])
for(i in c(1:nrow(total_group_missed_mz))){
  data2[grep(total_group_missed_mz$Group.1[i],data$class),grep(total_group_missed_mz$variable[i],colnames(data2))] <- min
  
}
data <- cbind(data[,c(1,2,3,4,5)],data2)
total_group_miss_imp <- data

#===================================missing due to low detection limit=============================
data <- total_group_miss_imp
miss_ratio_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),ratio_count) 
mean_intensity_by_class <- aggregate(data[,-c(1,2,3,4,5)],by=list(data$class),function(x){mean(as.numeric(as.character(x)),na.rm = T)})
a <- melt(miss_ratio_by_class,id=c("Group.1"))
names(a) <- c("group","mz","ratio")
a$gmz <- paste(a$group,a$mz,sep="_")
b <- melt(mean_intensity_by_class,id=c("Group.1"))
names(b) <- c("group","mz","intensity")
b$gmz <- paste(b$group,b$mz,sep="_")
ab <- merge(a,b,by.x = "gmz",by.y = "gmz")  
all_value<- as.numeric(as.matrix(data[,-c(1,2,3,4,5)]))
low_value <- quantile(all_value,0.05,na.rm = T)
all_low <- all_value[which(all_value<low_value)]
mean_l <- mean(all_low,na.rm = T)
sd_l <- sd(all_low,na.rm = T)
low_mz <- subset(ab,ratio>=0.5 & intensity<=low_value)  
#----------------
data2 <- as.matrix(data[,-c(1,2,3,4,5)])
for(i in c(1:nrow(low_mz))){
  class_mark <- grep(low_mz$group.x[i],data$class)
  mz_mark <- grep(low_mz$mz.x[i],colnames(data2))
  k<-which(is.na(data2[class_mark,mz_mark]))
  data2[class_mark[k],mz_mark]<- min
}
data <- cbind(data[,c(1,2,3,4,5)],data2) 
impute_low_group <- data
write.csv(impute_low_group,paste(input_dir,"impute_low_group.csv",sep="/"))
#===================================knn_impute===================================
data <- impute_low_group[,-c(1,2,3,4,5)]
data2 <- impute.knn(as.matrix(data))
data3 <- data2$data
data3[data3<0] <- 0.01
min_knn_impute <- cbind(impute_low_group[,c(1,2,3,4,5)],data3)
actual<- miss_del_datatl[,-c(1,2,3,4,5)]
predict<- min_knn_impute[,-c(1,2,3,4,5)]
RMSE_min_knn = (mean((predict - actual)^2))^0.5







#==================================accuracy===================================

pca_lda_classification_accuracy_cal2(min_knn_impute)
pca_lda_classification_accuracy_cal2(knn_impt)






















#===========================================================================


library(ROCR)
library(lda)
library(MASS)


pca_lda_classification_accuracy_cal2 <- function(objection){
  impdata <- as.matrix(objection[,-c(1,2,3,4,5)])
  row.names(impdata) <- objection$sample
  diagnosis <- as.numeric(as.factor(objection$class))
  #pca analysis
  imp.pr <- prcomp(impdata, scale = TRUE, center = TRUE)
  imp.pcs <- imp.pr$x[,1:10]
  imp.pcst <- cbind(imp.pcs, diagnosis)
  #train and test dataset
  N <- nrow(imp.pcst)
  rvec <- runif(N)
  imp.pcst.train <- imp.pcst[rvec < 0.75,]
  imp.pcst.test <- imp.pcst[rvec >= 0.75,]
  nrow(imp.pcst.test)
  #lda analysis
  imp.pcst.train.df <- as.data.frame(imp.pcst.train)
  imp.lda <- lda(diagnosis ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = imp.pcst.train.df)
  imp.pcst.test.df <- as.data.frame(imp.pcst.test)
  imp.lda.predict <- predict(imp.lda, newdata = imp.pcst.test.df)
  pred <- cbind(imp.pcst.test.df,imp.lda.predict$class)
  accuracy <- sum(pred$diagnosis == pred$`imp.lda.predict$class`)/nrow(pred)
  png(paste(out_dir,paste("validation_pca_lda","knn",".png",sep=""),sep="/"))
  plotda <- cbind(as.data.frame(imp.pr$x[, c(1, 2)]),objection$class)
  names(plotda)[3] <- "sgroup"
  ggplot(plotda,aes(PC1,PC2))+
    geom_point(aes(colour=sgroup))+
    stat_ellipse(aes(x=PC1, y=PC2,color=sgroup),type = "norm")
  dev.off()
  return(accuracy)
}



