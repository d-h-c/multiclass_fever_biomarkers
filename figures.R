

# colour pallette ####
lin = 1
titlesize = 1
library(RColorBrewer)
colour_table <- data.frame(row.names =  c(bacteria , "bacterial",viral,"viral",inflammatory,"inflammatory","KD","TB","malaria"))

pallette <- RColorBrewer::brewer.pal(n = 8,name = "Dark2")
viralcol <- pallette[3]
inflamcol <- pallette[7]
bactcol <- pallette[4]


colour_table[c("bacterial","inflammatory","viral"),1] <- c(bactcol,inflamcol,viralcol)
colour_table[bacteria,1] <-  colorRamp2(c(1,10),c(bactcol,"grey"))(2:7)
colour_table[inflammatory,1] <- colorRamp2(c(1,7),c(inflamcol,"grey"))(2:4)
colour_table[viral,1] <-  colorRamp2(c(1,10),c(viralcol,"grey"))(2:7)
colour_table[c("KD","TB","malaria"),1] <-  pallette[c(1,2,5)]

colour_table[c("bacterial","viral","inflammatory"),2] <- 15:17 
colour_table[bacteria,2] <- 1:6
colour_table[viral,2] <- 7:12
colour_table[inflammatory,2] <- c(13,14,8)
colour_table["TB",2] <- 20
colour_table["KD",2] <- 18
colour_table["malaria",2] <- 25
colour_table[,3] <- rep(2:6,12)[1:21]
colour_table[c("bacterial","viral","inflammatory"),3] <- 1
colour_table["malaria",3] <- 6
colour_table[,4] <- rownames(colour_table)

colour_table["bacterial",4] <- "Bacterial"
colour_table["viral",4] <- "Viral"
colour_table["malaria",4] <- "Malaria"
colour_table["TB",4] <- "Tuberculosis"
colour_table["meningococcal",4] <- "N. meningitidis"
colour_table["pneumo",4] <- "S. pneumoniae"
colour_table["adeno",4] <- "Adenovirus"
colour_table["staph",4] <- "S. aureus"
colour_table["gas",4] <- "GAS"
colour_table["gbs",4] <- "GBS"
colour_table["rhino",4] <- "Rhinovirus"
colour_table["enterovirus",4] <- "Enterovirus"
colour_table["ecoli",4] <- "E. coli"
colour_table["flu",4] <- "Influenza"
colour_table["inflammatory",4] <- "Inflammatory"
colnames(colour_table) <- c("col","shape","line","name")



# supplementary tables ####

## microarray dataset table ####

load("result_data_objects/dataset_plot_mat.RData")
ds_tab <- data.frame()
for(dataset in rev(colnames(dataset_plot_mat))[-1]){
  groups <- rownames(dataset_plot_mat)[dataset_plot_mat[,dataset]!=0]
  groups1 <- groups[groups%in% rownames(class_weights)]
  groups2 <- groups[!groups%in% rownames(class_weights)]
  groups2 <- groups2[groups2 != "dataset_total"]
  ds_tab[dataset,1] <- paste(groups1,collapse = " ")
  ds_tab[dataset,2] <- paste(groups2,collapse = " ")  
  numbers <- ""
  for(grp in groups1){
    numbers <- paste( numbers ,"N_",grp,"=", dataset_plot_mat[grp,dataset]," ",sep ="")
  }
  ds_tab[dataset,3] <- numbers
}
ds_tab <- ds_tab[-2,]

ds_tab["GSE63881",4] <- "No controls"
ds_tab["GSE25504",4] <- "Premature infants" # !!!!
ds_tab["GSE65391",4] <- "Only used initial timepoints"
ds_tab["GSE38900",4] <- "excluded due to failed probes"

colnames(ds_tab) <- c("Diseases","Excluded","Patient numbers","note")
write.table(ds_tab[,-2],sep = ",",file = "figures/dataset_table.csv")



## patients numbers ####
# microarray : 
table(c(phenotypes[phenotypes[,"trte"]=="train","group"],phenotypes[phenotypes[,"trte"]=="test","group"]))
#rnaseq :
table(c(rnaseq_validation_pheno[train_rnaseq],rnaseq_validation_pheno[test_rnaseq]))

## cost table ####
disease <- c("HSP",
             "Adenovirus",
             "Enterovirus",
             "HHV6",
             "Rhinovirus",
             "RSV",
             "Influenza",
             "JIA",
             "SLE",
             "Tuberculosis",
             "KD",
             "E. coli",
             "GAS",
             "GBS",
             "Malaria",
             "N. meningitidis",
             "S. pneumoniae",
             "S. aureus",
             "Viral",
             "Inflammatory",
             "Bacterial")

justification <-  c("Typically self-limiting condition, with very few consequences of misdiagnosis at first presentation.",
                    "Usually self-limiting viral infection, no specific treatment, misclassification may result in unnecessary antibiotic use.",
                    "Usually self-limiting viral infection, no specific treatment, misclassification may result in unnecessary antibiotic use.",
                    "Usually self-limiting viral infection, no specific treatment, misclassification may result in unnecessary antibiotic use.",
                    "Usually self-limiting viral infection, no specific treatment, misclassification may result in unnecessary antibiotic use.",
                    "Usually self-limiting viral infection, no specific treatment, misclassification may result in unnecessary antibiotic use.",
                    "Sometimes severe, viral infection. Treatable by early initiation of antiviral therapy.",
                    "Chronic inflammatory condition, predominantly involving joints, delay to diagnosis incurs significant healthcare costs and may result in joint damage.",
                    "Inflammatory condition, multiple organ involvement, delay to diagnosis incurs significant healthcare costs and may result in multiple end-organ organ damage.",
                    "Slowly progressive infection, delay in diagnosis can result in more severe and disseminated disease.",
                    "Inflammatory condition, delay to diagnosis carries high risk of coronary artery aneurysms.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Acute and life-threatening parasitic infection, requiring immediate diagnosis.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis.",
                    "Usually self-limiting viral infections, no specific treatment for most, misclassification may result in unnecessary antibiotic use.",
                    "Mostly slowly progressive but earlier diagnosis may improve outcomes.",
                    "Acute and life-threatening bacterial infection, requiring immediate diagnosis. ")

names(justification) <- disease

micro_weights <- rbind(class_weights, class_weights_toplevel[c("bacterial","viral","inflammatory"),])
micro_weights <- round(micro_weights,digits = 3)
rownames(micro_weights) <- colour_table[rownames(micro_weights) ,"name"]

micro_weights[,"Justification"] <- justification[rownames(micro_weights)]
path_order <- colour_table[c(bacteria,"TB",viral,inflammatory,"KD","malaria","bacterial","viral","inflammatory"),"name"]

rnaseq_weights <- refit_class_weights
rnaseq_weights <- round(rnaseq_weights,digits = 3)
rnaseq_weights <- rnaseq_weights[rownames(rnaseq_weights) != "HC",]
rownames(rnaseq_weights) <- colour_table[rownames(rnaseq_weights) ,"name"]
rnaseq_weights[is.na(rnaseq_weights[,1]),1] <- 0

write.csv(micro_weights[path_order , c("freq","cost","Justification","cost+imbalance")] ,
          file = "figures/class_weights_microarray.csv")
path_order.tmp <- path_order[path_order!="Inflammatory"]
path_order.tmp <- path_order.tmp[rnaseq_weights[path_order.tmp,"freq"]>0]
write.csv(rnaseq_weights[path_order.tmp, c("freq","cost","cost+imbalance")] , 
          file = "figures/class_weights_rnaseq.csv")



# microarray 18C ROC ####
pdf("figures/lowlevel_roc_ma.dec.pdf" ,width = 10,height = 10)

sf <- 1
auc_low <- plot_curve_loc(performance_list$test$full$low.dec$ROCs,
                          more_pr_loc =  performance_list$test$full$low.dec$PR_mod ,
                          type = "ROC",linethickness = F, sf =3 ) #,sf = 50)
title( main =list("161 transcript signature predicting 18 specific disease classes" ,cex =titlesize),line = 0.5,adj = 0)
auc_low
rownames(auc_low) <- auc_low[,1]
confusion.loc <- confusion_microarray_test[[2]]
# decoupling probabilities does not change the max(p) confusion
auc_low[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp","tp","pos")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in colnames(confusion.loc)){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  print(c(class, specificity ,recall, precision ))
  auc_low[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp","tp","pos")] <- c(recall , specificity,precision,tp, tp+fp)
}
auc_low["macro average","sens_at_maxp"]  <- mean(auc_low[,"sens_at_maxp"],na.rm = T)
auc_low["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_low["macro average","spec_at_maxp"]  <- mean(auc_low[,"spec_at_maxp"],na.rm = T)
auc_low["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
auc_low["macro average","prec_at_maxp"]  <- mean(auc_low[,"prec_at_maxp"],na.rm = T)
auc_low["micro average","prec_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fps))
points(auc_low["macro average","spec_at_maxp"] , auc_low["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_low["micro average","spec_at_maxp"] , auc_low["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
dev.off()

reformat_auc_table <- function(auc_low,comp = F){
  auc_low.form <- auc_low[rownames(auc_low)%in%rownames(colour_table),]
  auc_low.form[,1] <- colour_table[rownames(auc_low.form),"name"]
  rownames(auc_low.form) <- auc_low.form[,1]
  auc_low.form <- auc_low.form[path_order[path_order%in% rownames(auc_low.form)],]
  for( i in 2:(ncol(auc_low.form)-1)){
    auc_low.form[,i] <- round(digits = 3 , as.numeric(auc_low.form[,i]))  
  }
  auc_low.form[,"CI"] <- paste(auc_low.form[,"CI.l"],auc_low.form[,"CI.h"], sep ="-")
  auc_low.form[,"bn_CI"] <- paste(auc_low.form[,"bn_CI.l"],auc_low.form[,"bn_CI.h"], sep ="-")
  if(comp){
    return(auc_low.form[,c("disease", "#cases","AUC","CI","sens_at_maxp", "spec_at_maxp","bn_AUC","bn_CI","Composition")])
  }else{
    return(auc_low.form[,c("disease", "#cases","AUC","CI","sens_at_maxp", "spec_at_maxp","bn_AUC","bn_CI")])
  }
}


for(i in names(hier.ma)){
  r <- performance_list$test$full$low.dec$ROCs[[paste0(i,"-binormal")]]
  if(!is.na(r)){
    auc_low[i,"bn_AUC"] <- r$auc
    auc_low[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_low[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
}

write.table(reformat_auc_table(auc_low ) ,row.names = F, file = "figures/lowlevel_roc_ma.dec.tab",sep = "\t")




pdf("figures/training_set_lowlevel_roc_ma.dec.pdf" ,width = 10,height = 10)
sf <- 1
auc_low.train <- plot_curve_loc(performance_list$training$full$low.dec$ROCs,
                                more_pr_loc =  performance_list$training$full$low.dec$PR_mod ,
                                type = "ROC",linethickness = F, sf =3 ) #,sf = 50)
title( main =list(" 161 transcript signature predicting 18 specific disease classes in training data" ,cex =titlesize),line = 0.5,adj = 0)
rownames(auc_low.train) <- auc_low.train[,1]
confusion.loc <- confusion_microarray_train[[2]]

auc_low.train[,c("sens_at_maxp","spec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in colnames(confusion.loc)){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  print(c(class, specificity ,recall, precision ))
  auc_low.train[class,c("sens_at_maxp","spec_at_maxp")] <- c(recall , specificity)
}
auc_low.train["macro average","sens_at_maxp"]  <- mean(auc_low.train[,"sens_at_maxp"],na.rm = T)
auc_low.train["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_low.train["macro average","spec_at_maxp"]  <- mean(auc_low.train[,"spec_at_maxp"],na.rm = T)
auc_low.train["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_low.train["macro average","spec_at_maxp"] , auc_low.train["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_low.train["micro average","spec_at_maxp"] , auc_low.train["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
dev.off()

for(i in names(hier.ma)){
  r <- performance_list$train$full$low.dec$ROCs[[paste0(i,"-binormal")]]
  if(!is.na(r)){
    auc_low.train[i,"bn_AUC"] <- r$auc
    auc_low.train[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_low.train[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
}

write.table(reformat_auc_table(auc_low.train) ,row.names = F, file = "figures/training_lowlevel_roc_ma.dec.tab",sep = "\t")


# microarray 18C ROCs split up ####
pdf("figures/lowlevel_roc_ma_split.dec.pdf" ,width = 8,height = 10)
confusion.loc <- confusion_microarray_test[[2]]
par(mfrow = c(5,4))
lin <- 2
par(mar=c(3.5,3.5,2,2))
for(class in colnames(confusion.loc)){
  plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
  title(main =  colour_table[class,"name"],line =1)
  title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
  axis(1,pos = 0)
  axis(2,las =2,pos = 1)
  points(c(1,0),c(1,1),type = "l")
  points(c(1,0),c(0,1),type = "l",col ="grey")
  par(xpd = T)
  points(performance_list$test$full$low.dec$ROCs[[class]]$specificities,
         performance_list$test$full$low.dec$ROCs[[class]]$sensitivities,
         col = colour_table[class , "col"],
         type = "l",lwd =2)
  points(auc_low[class,"spec_at_maxp"] , auc_low[class,"sens_at_maxp"] , col = colour_table[class , "col"],pch =16 ,cex =2)
  par(xpd = F)
}

plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Macro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")
points(performance_list$test$full$low.dec$PR_mod$specificity.macro,
       performance_list$test$full$low.dec$PR_mod$sensitivity.macro,
       type ="l")
points(auc_low["macro average","spec_at_maxp"] , auc_low["macro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)


plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Micro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")
points(performance_list$test$full$low.dec$PR_mod$specificity.micro,
       performance_list$test$full$low.dec$PR_mod$sensitivity.micro,
       type = "l")
points(auc_low["micro average","spec_at_maxp"] , auc_low["micro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)
dev.off()


pdf("figures/lowlevel_PR_ma_split.dec.pdf" ,width = 8,height = 10)
confusion.loc <- confusion_microarray_test[[2]]

par(mfrow = c(5,4))
lin <- 2
par(mar=c(3.5,3.5,2,2))
for(class in colnames(confusion.loc)){
  
  plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(0,1),axes =F,xlab = "",ylab = "")
  title(main =  colour_table[class,"name"],line =1)
  title(xlab = "sensitivity",ylab = "precision",line =lin)
  axis(1,pos = 0)
  axis(2,las =2,pos = 0)
  points(c(1,0),c(1,1),type = "l")
  par(xpd = T)
  points(performance_list$test$full$low.dec$PR[[class]]$curve[,1],
         performance_list$test$full$low.dec$PR[[class]]$curve[,2],
         col = colour_table[class , "col"],
         type = "l",lwd =2)
  
  points(auc_low[class,"sens_at_maxp"] , auc_low[class,"prec_at_maxp"] ,
         col = colour_table[class , "col"],pch =16 ,cex =2)
  par(xpd = F)
  
  ncase <- sum(test_phenotypes_indicator[,class])
  ncontrol <- nrow(test_phenotypes_indicator)
  abline(h = ncase / ncontrol , col = colour_table[class,"col"],lty = colour_table[class,"line"])
  
}

plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(0,1),axes =F,xlab = "",ylab = "")
title(main = "Macro-average",line =1)
title(xlab = "sensitivity",ylab = "precision",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(performance_list$test$full$low.dec$PR_mod$sensitivity.macro,
       performance_list$test$full$low.dec$PR_mod$precision.macro,
       type ="l")
points(auc_low["macro average","sens_at_maxp"] , auc_low["macro average","prec_at_maxp"] , col = "black",pch =16 ,cex =2)


plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(0,1),axes =F,xlab = "",ylab = "")
title(main = "Micro-average",line =1)
title(xlab = "sensitivity",ylab = "precision",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(performance_list$test$full$low.dec$PR_mod$sensitivity.micro,
       performance_list$test$full$low.dec$PR_mod$precision.micro,
       type = "l")
points(auc_low["micro average","sens_at_maxp"] , auc_low["micro average","prec_at_maxp"] , col = "black",pch =16 ,cex =2)

dev.off()



# PCA in signature space to justify TB and KD separation  ####


pca.sigspace.ma <- prcomp(t(expression_filt[selgenes,training_samples_all ]),center = T,scale. = T)
phenos <- unique(training_phenotypes[,"group"])
phecols <- colour_table[phenos,"col"]
names(phecols)  <-  phenos
plot_pca <- function(xpc,ypc,textlab = T,annot = T,lab_pos ){ 
  par(mar = c(7.5,7.5,3,3))
  groups <- training_phenotypes[training_samples_all,"group"]
  groups[groups%in%bacteria] <- "bacterial"
  groups[groups%in%viral] <- "viral"
  groups[groups%in%inflammatory] <- "inflammatory"
  par(xpd=T)
  plot(pca.sigspace.ma$x[,xpc],pca.sigspace.ma$x[,ypc] , 
       cex =1.5, col = colour_table[groups,"col"] ,axes = F,xlab ="", ylab = "",
       pch=colour_table[training_phenotypes[training_samples_all,"group"],"shape"])
  if(textlab == T)  text(pca.sigspace.ma$x[,xpc],pca.sigspace.ma$x[,ypc] + 0.2, labels = training_phenotypes[training_samples_all,"group"],cex = 0.5)
  groups  <- as.factor(groups)
  df <- cbind.data.frame(groups,pca.sigspace.ma$x[,c(xpc,ypc)])
  bb <- coord.ellipse(df,level.conf = 0.95)
  for(grp in levels(groups)){
    points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =colour_table[grp,"col"],lwd = 2)
    dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,1])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,2])**2
    cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
    if(annot == T){
      text(lab_pos[grp,1],lab_pos[grp,2] , colour_table[grp,"name"],pos = lab_pos[grp,3],col = colour_table[grp,"col"],cex = 3)
      points(c(cp[1],lab_pos[grp,1]) ,c(cp[2],lab_pos[grp,2]) , type ="l",col = colour_table[grp,"col"],lwd = 2) 
    }
  }
  
  axis(side = 1,tick = T,line = 3)
  axis(side = 2,tick=T,lwd =1,lwd.ticks = 1,las = 2,line = 3)
  title(xlab = paste("PC",xpc," (variance explained:", round(digits = 1,summary(pca.sigspace.ma)$importance[2,xpc]*100),"%)"),
        ylab = paste("PC",ypc," (variance explained:", round(digits = 1,summary(pca.sigspace.ma)$importance[2,ypc]*100),"%)"),line = 6,cex.lab =1.5)
  par(xpd=FALSE)
  
}

pdf("figures/signature_space_pca.pdf" ,width = 15,height = 15)
lab_pos <- data.frame(matrix(byrow = T,ncol = 3 ,data = c(-12,4,3,
                                                          -9, 7,3,
                                                          -3, -16,2,
                                                          -7,-12,1,
                                                          -10,-10,1,
                                                          9.5,-5 ,4 )))
rownames(lab_pos) <- c("viral",
                       "inflammatory",
                       "TB",
                       "bacterial",
                       "malaria",
                       "KD")

plot_pca(1,2,textlab = F, lab_pos = lab_pos )

lab_pos["KD",] <- c(6,8,4)
lab_pos["malaria",] <- c(-8,10,2)
lab_pos["TB",] <- c(0,-10,4)
lab_pos["bacterial",] <- c(-7,-10,1)
plot_pca(2,3,textlab = F, lab_pos = lab_pos)

lab_pos["TB",] <- c(0,12,4)
lab_pos["inflammatory",] <- c(-6,-7,2)
lab_pos["bacterial",] <- c(-4,-10,1)
lab_pos["KD",] <- c(6,-10,4)
lab_pos["malaria",] <- c(6,9,3)
lab_pos["viral",] <- c(2,-10,1)
plot_pca(3,4,textlab = F, lab_pos = lab_pos)

lab_pos["KD",] <- c(-6,8,2)
lab_pos["malaria",] <- c(7,7,4)
lab_pos["TB",] <- c(7,-8,1)
lab_pos["bacterial",] <- c(-4,-10,1)
lab_pos["inflammatory",] <- c(-9,-8,1)

plot_pca(4,5,textlab = F, lab_pos = lab_pos)
dev.off()



# microarray 6C ROC #### 

pdf("figures/broadlevel_roc_ma.dec.pdf" ,width = 10,height = 10)
titlesize <- 1
plotclasses.loc <- colnames( confusion_microarray_test_top[[2]]) 
auc_top <- plot_curve_loc(performance_list$test$full$top.dec$ROCs,
                          more_pr_loc = performance_list$test$full$top.dec$PR_mod,
                          plotclasses = plotclasses.loc,type = "ROC",sf = 2,linethickness = F)
title( main =list("161 transcript signature predicting broad disease classes" ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_top) <- auc_top[,1]
confusion.loc <- confusion_microarray_test_top[[2]]

auc_top[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in plotclasses.loc){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  print(c(class, specificity ,recall, precision ))
  auc_top[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity, precision)
}
auc_top["macro average","sens_at_maxp"]  <- mean(auc_top[,"sens_at_maxp"],na.rm = T)
auc_top["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_top["macro average","spec_at_maxp"]  <- mean(auc_top[,"spec_at_maxp"],na.rm = T)
auc_top["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_top["macro average","spec_at_maxp"] , auc_top["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_top["micro average","spec_at_maxp"] , auc_top["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
write.table(reformat_auc_table(auc_top) ,row.names = F, file = "figures/broadlevel_roc_ma.dec.tab",sep = "\t")
dev.off()

for(i in unique(hier.ma)){
  r <- performance_list$test$full$top.dec$ROCs[[paste0(i,"-binormal")]]
  if(!is.na(r)){
    auc_top[i,"bn_AUC"] <- r$auc
    auc_top[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_top[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
}

hier.ma["KD"] <- "KD"
auc_top_mod <- auc_top
for(i in 3:8){
  class <- auc_top_mod[i,1]
  subclass <- names(hier.ma[hier.ma == class])
  if(length(subclass) > 1){
    mat <- auc_low[auc_low[,"disease"] %in% subclass,c(1,3)]
    mat[,1] <- colour_table[mat[,1],"name"]
    mat[,2] <- gsub("^","n=",mat[,2])
    auc_top_mod[i,"Composition"] <- paste(rowpaste(mat), collapse = ", ")
  }else{
    auc_top_mod[i,"Composition"] <- ""
  }
}

write.table(reformat_auc_table(auc_top_mod, comp = T) ,row.names = F, file = "figures/broadlevel_roc_ma.dec_plus_group_composition.tab",sep = "\t")


# training microarray broad
pdf("figures/training_broadlevel_roc_ma.dec.pdf" ,width = 10,height = 10)
titlesize <- 1
plotclasses.loc <- colnames( confusion_microarray_train_top[[2]]) 
auc_top.train <- plot_curve_loc(performance_list$training$full$top.dec$ROCs,
                          more_pr_loc = performance_list$training$full$top.dec$PR_mod,
                          plotclasses = plotclasses.loc,type = "ROC",sf = 2,linethickness = F)
title( main =list("161 transcript signature predicting broad disease classes on training data" ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_top.train) <- auc_top.train[,1]
confusion.loc <- confusion_microarray_train_top[[2]]
auc_top.train[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in plotclasses.loc){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  print(c(class, specificity ,recall, precision ))
  auc_top.train[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity,precision)
}

auc_top.train["macro average","sens_at_maxp"]  <- mean(auc_top.train[,"sens_at_maxp"],na.rm = T)
auc_top.train["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_top.train["macro average","spec_at_maxp"]  <- mean(auc_top.train[,"spec_at_maxp"],na.rm = T)
auc_top.train["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_top.train["macro average","spec_at_maxp"] , auc_top.train["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_top.train["micro average","spec_at_maxp"] , auc_top.train["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
write.table(reformat_auc_table(auc_top.train) ,row.names = F, file = "figures/training_broadlevel_roc_ma.dec.tab",sep = "\t")
dev.off()





## microarray 6C ROCs splitup ####
pdf("figures/broadlevel_roc_ma_split.dec.pdf" ,width = 8,height = 4)
confusion.loc <- confusion_microarray_test_top[[2]]
par(mfrow = c(2,4))
lin <- 2
par(mar=c(3.5,3.5,2,2))
for(class in colnames(confusion.loc)){
  plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
  title(main =  colour_table[class,"name"],line =1)
  title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
  axis(1,pos = 0)
  axis(2,las =2,pos = 1)
  points(c(1,0),c(1,1),type = "l")
  points(c(1,0),c(0,1),type = "l",col ="grey")
  par(xpd = T)
  points(performance_list$test$full$top.dec$ROCs[[class]]$specificities,
         performance_list$test$full$top.dec$ROCs[[class]]$sensitivities,
         col = colour_table[class , "col"],
         type = "l",lwd =2)
  points(auc_top[class,"spec_at_maxp"] , auc_top[class,"sens_at_maxp"] ,
         col = colour_table[class , "col"],pch =16 ,cex =2)
  par(xpd = F)
}

plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Macro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")

points(performance_list$test$full$top.dec$PR_mod$specificity.macro,
       performance_list$test$full$top.dec$PR_mod$sensitivity.macro,
       type ="l")
points(auc_top["macro average","spec_at_maxp"] , auc_top["macro average","sens_at_maxp"] ,
       col = "black",pch =16 ,cex =2)


plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Micro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")
points(performance_list$test$full$top.dec$PR_mod$specificity.micro,
       performance_list$test$full$top.dec$PR_mod$sensitivity.micro,
       type = "l")
points(auc_top["micro average","spec_at_maxp"] , auc_top["micro average","sens_at_maxp"] ,
       col = "black",pch =16 ,cex =2)
dev.off()


## microarray confusions  + level disagreement ####
pdf("figures/microarray_confusion_splitbydisagreement.pdf" ,width = 15,height = 15)
colord <- c(bacteria,viral,inflammatory,"TB","KD","malaria")
confusion_heatmap(confusion_microarray_test[[2]][colord,colord])
title( main =list("all microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(microarray_disagreement$confident.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test[[2]][colord,colord]),square_lab = "number")
title( main =list("confident microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(microarray_disagreement$ambiguous.low[colord,colord],scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test[[2]][colord,colord]) ,square_lab = "number")
title( main =list("ambiguous microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)


colord <- c("bacterial","viral","inflammatory","TB","KD","malaria")
confusion_heatmap(confusion_microarray_test_top[[2]][colord,colord],level = "top" )
title( main =list("all microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(microarray_disagreement$confident.top[colord,colord],level = "top"  ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test_top[[2]][colord,colord]),square_lab = "number")
title( main =list("confident microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(microarray_disagreement$ambiguous.top[colord,colord],level = "top" ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test_top[[2]][colord,colord]),square_lab = "number")
title( main =list("ambiguous microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
dev.off()

# RNA-Seq ROcs direct transfer ####

expr.base <-  expr.allprobe.voom[probe_order, ]
sams.loc <- names(rnaseq_validation_pheno)
sams.ord <- c()
for(i in group_ord){
  sams.ord <- c(sams.ord,sams.loc[rnaseq_validation_pheno[sams.loc] == i])
}

library(pROC)
pdf("figures/rnaseq_validation_rocs_allprobes_top.pdf" ,width = 21,height = 14)
par(mfrow = c(2,3), xpd = F)
for(disease in names(coef(ridge_refit.top))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(ridge_refit.allprobe.top)[[disease]][rownames(expr.loc),1]
  expr.loc <- expr.loc / (2*max(abs(expr.loc), na.rm = T)) + 0.5
  expr.loc <- as.data.frame(expr.loc)
  expr.loc <- expr.loc[probe_order,sams.ord]
  ind <- hier.ma[rnaseq_validation_pheno[colnames(expr.loc)]] == disease
  pretty_plot_area(ytick = c(0,1,1),cols = c("white","lightgray"),
                   xtick = c(1,0,-1),show_x_tick = T,show_y_tick = T,
                   show_x = T,show_y = T,plot_outside_margin = F, margins = rep(6,4))
  abline(1,-1, col = "lightgrey")
  par(xpd = T)
  vals <- colSums(expr.loc, na.rm = T)
  vals <- vals - min(vals)
  vals <- vals / max(vals)
  r <- roc(cases = vals[ind], controls = vals[!ind], ci =T)
  plot(r,add = T, col  = colour_table[disease,"col"])
  plot(ci.sp(r, sensitivities=seq(0, 1, .01)), type="shape", col  = add.alpha(colour_table[disease,"col"],0.6))
  text(0.4,0.1,paste0("AUC=",round(auc(r),digits = 3)," (",paste0(round(digits = 3,ci.auc(r)[c(1,3)]), collapse = "-"),")",
                      "\n n=",sum(ind)), cex = 2,  col  = colour_table[disease,"col"])
  title(main = colour_table[disease,"name"], cex = 2)
}
dev.off()

pdf("figures/rnaseq_validation_rocs_allprobes_low.pdf" ,width = 21,height = 35)
par(mfrow = c(5,3))
for(disease in names(coef(ridge_refit))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(ridge_refit.allprobe.low)[[disease]][rownames(expr.loc),1]
  expr.loc <- expr.loc / (2*max(abs(expr.loc),na.rm = T)) + 0.5
  expr.loc <- as.data.frame(expr.loc)
  expr.loc[probe_order[probe_order%ni%rownames(expr.loc)],] <- NA
  expr.loc <- expr.loc[probe_order,sams.ord]
  ind <- rnaseq_validation_pheno[colnames(expr.loc)] == disease
  if(sum(ind) > 0){
    vals <- colSums(expr.loc, na.rm = T)
    vals <- vals - min(vals)
    vals <- vals / max(vals)
    pretty_plot_area(ytick = c(0,1,1),cols = c("white","lightgray"),
                     xtick = c(1,0,-1),show_x_tick = T,show_y_tick = T,
                     show_x = T,show_y = T,plot_outside_margin = F, margins = rep(6,4))
    abline(1,-1, col = "lightgrey")
    par(xpd = T)
    r <- roc(cases = vals[ind], controls = vals[!ind], ci= T)
    plot(r,add = T, col  = colour_table[disease,"col"])
    plot(ci.sp(r, sensitivities=seq(0, 1, .01)), type="shape", col  = add.alpha(colour_table[disease,"col"],0.6))
    text(0.4,0.1,paste0("AUC=",round(auc(r),digits = 3)," (",paste0(round(digits = 3,ci.auc(r)[c(1,3)]), collapse = "-"),")",
                        "\n n=",sum(ind)), cex = 2,  col  = colour_table[disease,"col"])
    title(main = colour_table[disease,"name"],cex = 2)  
  }
}
dev.off()




# RNA-Seq ROC with new coefficients ####
pdf("figures/rnaseq_rocs.dec.pdf" ,width = 10,height = 10)
plotclasses.loc <- colnames(rnasesq_test_confusion.low[[2]])
auc_low.rnaseq <- plot_curve_loc(rnaseq_validation_performance_mn_low.dec$ROCs ,
                          plotclasses = plotclasses.loc,
                          type = "ROC",sf = 2,linethickness = F,dstype = "rnaseq" ,
                          more_pr_loc = rnaseq_validation_performance_mn_low.dec$PR_mod)
title( main =list("ROC for 145 gene signature predicting in RNA-Seq validation test set"
                  ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_low.rnaseq) <- auc_low.rnaseq[,1]
confusion.loc <- rnasesq_test_confusion.low[[2]]
auc_low.rnaseq[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in colnames(confusion.loc)){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  auc_low.rnaseq[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity,precision)
}

auc_low.rnaseq["macro average","sens_at_maxp"]  <- mean(auc_low.rnaseq[,"sens_at_maxp"],na.rm = T)
auc_low.rnaseq["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_low.rnaseq["macro average","spec_at_maxp"]  <- mean(auc_low.rnaseq[,"spec_at_maxp"],na.rm = T)
auc_low.rnaseq["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_low.rnaseq["macro average","spec_at_maxp"] , auc_low.rnaseq["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_low.rnaseq["micro average","spec_at_maxp"] , auc_low.rnaseq["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)

for(i in 3:15){
  class <- auc_low.rnaseq[i,1]
  r <- rnaseq_validation_performance_mn_low.dec$ROCs[[paste0(class,"-binormal")]]
  if(!is.na(r)){
    auc_low.rnaseq[i,"bn_AUC"] <- r$auc
    auc_low.rnaseq[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_low.rnaseq[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
}

write.table(reformat_auc_table(auc_low.rnaseq) ,row.names = F, file = "figures/lowlevel_roc_rnaseq.dec.tab",sep = "\t")

# broad class
plotclasses.loc <- colnames(rnasesq_test_confusion.top[[2]])
auc_top.rnaseq <- plot_curve_loc(rnaseq_validation_performance_mn_top.dec$ROCs ,
                          plotclasses =  plotclasses.loc,
                          type = "ROC",sf = 2,linethickness = F,dstype = "rnaseq" ,
                          more_pr_loc = rnaseq_validation_performance_mn_top.dec$PR_mod)
title( main =list("ROC for 145 gene signature predicting broad class in RNA-Seq validation test set" ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_top.rnaseq) <- auc_top.rnaseq[,1]
confusion.loc <- rnasesq_test_confusion.top[[2]]
auc_top.rnaseq[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in plotclasses.loc){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch =colour_table[class , "shape"] ,cex =2)
  auc_top.rnaseq[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity,precision)
}
auc_top.rnaseq["macro average","sens_at_maxp"]  <- mean(auc_top.rnaseq[,"sens_at_maxp"],na.rm = T)
auc_top.rnaseq["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_top.rnaseq["macro average","spec_at_maxp"]  <- mean(auc_top.rnaseq[,"spec_at_maxp"],na.rm = T)
auc_top.rnaseq["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_top.rnaseq["macro average","spec_at_maxp"] , auc_top.rnaseq["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_top.rnaseq["micro average","spec_at_maxp"] , auc_top.rnaseq["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
write.table(reformat_auc_table(auc_top.rnaseq) ,row.names = F, file = "figures/broadlevel_roc_rnaseq.dec.tab",sep = "\t")
dev.off()

for(i in 3:8){
  class <- auc_top.rnaseq[i,1]
  r <- rnaseq_validation_performance_mn_top.dec$ROCs[[paste0(class,"-binormal")]]
  if(!is.na(r)){
    auc_top.rnaseq[i,"bn_AUC"] <- r$auc
    auc_top.rnaseq[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_top.rnaseq[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
}

auc_top_mod <- auc_top.rnaseq

for(i in 3:8){
  class <- auc_top_mod[i,1]
  r <- rnaseq_validation_performance_mn_top.dec$ROCs[[paste0(class,"-binormal")]]
  if(!is.na(r)){
    auc_top_mod[i,"bn_AUC"] <- r$auc
    auc_top_mod[i,"bn_CI.l"] <- round( r$ci,digits = 4)[1]
    auc_top_mod[i,"bn_CI.h"] <- round( r$ci,digits = 4)[3]
  }
  
  subclass <- names(hier.rna[hier.rna == class])
  if(length(subclass) > 1){
    mat <- auc_low.rnaseq[auc_low.rnaseq[,"disease"] %in% subclass,c(1,3)]
    mat[,1] <- colour_table[mat[,1],"name"]
    mat[,2] <- gsub("^","n=",mat[,2])
    auc_top_mod[i,"Composition"] <- paste(rowpaste(mat), collapse = ", ")
  }else{
    auc_top_mod[i,"Composition"] <- ""
  }
}

write.table(reformat_auc_table(auc_top_mod,comp = T) ,row.names = F, file = "figures/broadlevel_roc_rnaseq.dec_plus_group_composition.tab",sep = "\t")


## training set ROC ####

pdf("figures/training_rnaseq_rocs.dec.pdf" ,width = 10,height = 10)
plotclasses.loc <- colnames(rnasesq_train_confusion.low[[2]])
auc_low.rnaseq.train <- plot_curve_loc(rnaseq_validation_performance_mn_low.dec.train$ROCs ,
                                 plotclasses = plotclasses.loc,
                                 type = "ROC",sf = 2,linethickness = F,dstype = "rnaseq" ,
                                 more_pr_loc = rnaseq_validation_performance_mn_low.dec.train$PR_mod)
title( main =list("ROC for 145 gene signature predicting in RNA-Seq training set" ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_low.rnaseq.train) <- auc_low.rnaseq.train[,1]
confusion.loc <- rnasesq_train_confusion.low[[2]]
auc_low.rnaseq.train[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in colnames(confusion.loc)){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],
         pch =colour_table[class , "shape"] ,cex =2)
  auc_low.rnaseq.train[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity, precision)
}

auc_low.rnaseq.train["macro average","sens_at_maxp"]  <- mean(auc_low.rnaseq.train[,"sens_at_maxp"],na.rm = T)
auc_low.rnaseq.train["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_low.rnaseq.train["macro average","spec_at_maxp"]  <- mean(auc_low.rnaseq.train[,"spec_at_maxp"],na.rm = T)
auc_low.rnaseq.train["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_low.rnaseq.train["macro average","spec_at_maxp"] , auc_low.rnaseq.train["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_low.rnaseq.train["micro average","spec_at_maxp"] , auc_low.rnaseq.train["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)
write.table(reformat_auc_table(auc_low.rnaseq.train) ,row.names = F, file = "figures/training_lowlevel_roc_rnaseq.dec.tab",sep = "\t")

# broad class
plotclasses.loc <- colnames(rnasesq_train_confusion.top[[2]])
auc_top.rnaseq.train <- plot_curve_loc(rnaseq_validation_performance_mn_top.dec.train$ROCs ,
                                 plotclasses =  plotclasses.loc,
                                 type = "ROC",sf = 2,linethickness = F,dstype = "rnaseq" ,
                                 more_pr_loc = rnaseq_validation_performance_mn_top.dec.train$PR_mod)
title( main =list("ROC for 145 gene signature predicting broad class in RNA-Seq validation training set" ,cex =titlesize),line = 0.2,adj = 0)
rownames(auc_top.rnaseq.train) <- auc_top.rnaseq.train[,1]
confusion.loc <- rnasesq_train_confusion.top[[2]]

auc_top.rnaseq.train[,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- NA
tps <- c()
fps <- c()
tns <- c()
fns <- c()
for(class in plotclasses.loc){
  tp <- confusion.loc[class,class]
  fp <- sum(confusion.loc[,class]) - tp 
  fn <- sum(confusion.loc[class,]) - tp
  tn <- sum(confusion.loc) - tp - fp - fn
  tps <- c(tps,tp)
  fps <- c(fps,fp)
  tns <- c(tns,tn)
  fns <- c(fns,fn)
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  points(specificity , recall , col = colour_table[class , "col"],pch = colour_table[class , "shape"] ,cex =2)
  auc_top.rnaseq.train[class,c("sens_at_maxp","spec_at_maxp","prec_at_maxp")] <- c(recall , specificity, precision)
  
}

auc_top.rnaseq.train["macro average","sens_at_maxp"]  <- mean(auc_top.rnaseq.train[,"sens_at_maxp"],na.rm = T)
auc_top.rnaseq.train["micro average","sens_at_maxp"]  <- mean(tps)/(mean(tps)+mean(fns))
auc_top.rnaseq.train["macro average","spec_at_maxp"]  <- mean(auc_top.rnaseq.train[,"spec_at_maxp"],na.rm = T)
auc_top.rnaseq.train["micro average","spec_at_maxp"]  <- mean(tns)/(mean(tns)+mean(fps))
points(auc_top.rnaseq.train["macro average","spec_at_maxp"] , auc_top.rnaseq.train["macro average","sens_at_maxp"] , col = "black",pch = 1 ,cex =2)
points(auc_top.rnaseq.train["micro average","spec_at_maxp"] , auc_top.rnaseq.train["micro average","sens_at_maxp"] , col = "black",pch = 16 ,cex =2)


write.table(reformat_auc_table(auc_top.rnaseq.train) ,row.names = F, file = "figures/training_broadlevel_roc_rnaseq.dec.tab",sep = "\t")

dev.off()




# RNAseq ROCs splitup ####

## LOWLEVEL
pdf("figures/lowlevel_roc_rnaseq_split.dec.pdf" ,width = 8,height = 8)
confusion.loc <- rnasesq_test_confusion.low[[2]]
par(mfrow = c(4,4)) ## 14 plots
lin <- 2
par(mar=c(3.5,3.5,2,2))
for(class in colnames(confusion.loc)){
  plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
  title(main =  colour_table[class,"name"],line =1)
  title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
  axis(1,pos = 0)
  axis(2,las =2,pos = 1)
  points(c(1,0),c(1,1),type = "l")
  points(c(1,0),c(0,1),type = "l",col ="grey")
  par(xpd = T)
  points(rnaseq_validation_performance_mn_low.dec$ROCs[[class]]$specificities,
         rnaseq_validation_performance_mn_low.dec$ROCs[[class]]$sensitivities,
         col = colour_table[class , "col"],
         type = "l",lwd =2)
  points(auc_low.rnaseq[class,"spec_at_maxp"] , auc_low.rnaseq[class,"sens_at_maxp"] , 
         col = colour_table[class , "col"],pch =16 ,cex =2)
  par(xpd = F)
}
plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Macro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")

points(rnaseq_validation_performance_mn_low.dec$PR_mod$specificity.macro,
       rnaseq_validation_performance_mn_low.dec$PR_mod$sensitivity.macro,
       type ="l")
points(auc_low.rnaseq["macro average","spec_at_maxp"] ,
       auc_low.rnaseq["macro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)

plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Micro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")
points(rnaseq_validation_performance_mn_low.dec$PR_mod$specificity.micro,
       rnaseq_validation_performance_mn_low.dec$PR_mod$sensitivity.micro,
       type = "l")
points(auc_low.rnaseq["micro average","spec_at_maxp"] , auc_low.rnaseq["micro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)
dev.off()

# BROAD level

pdf("figures/broadlevel_roc_rnaseq_split.dec.pdf" ,width = 8,height = 4)
confusion.loc <- rnasesq_test_confusion.top[[2]]
par(mfrow = c(2,4)) ## 12 plots
lin <- 2
par(mar=c(3.5,3.5,2,2))
for(class in colnames(confusion.loc)){
  plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
  title(main =  colour_table[class,"name"],line =1)
  title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
  axis(1,pos = 0)
  axis(2,las =2,pos = 1)
  points(c(1,0),c(1,1),type = "l")
  points(c(1,0),c(0,1),type = "l",col ="grey")
  par(xpd = T)
  points(rnaseq_validation_performance_mn_top.dec$ROCs[[class]]$specificities,
         rnaseq_validation_performance_mn_top.dec$ROCs[[class]]$sensitivities,
         col = colour_table[class , "col"],
         type = "l",lwd =2)
  points(auc_top.rnaseq[class,"spec_at_maxp"] , auc_top.rnaseq[class,"sens_at_maxp"] , col = colour_table[class , "col"],pch =16 ,cex =2)
  par(xpd = F)
}

plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Macro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")

points(rnaseq_validation_performance_mn_top.dec$PR_mod$specificity.macro,
       rnaseq_validation_performance_mn_top.dec$PR_mod$sensitivity.macro,
       type ="l")
points(auc_top.rnaseq["macro average","spec_at_maxp"] , auc_top.rnaseq["macro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)


plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
title(main = "Micro-average",line =1)
title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
axis(1,pos = 0)
axis(2,las =2,pos = 1)
points(c(1,0),c(1,1),type = "l")
points(c(1,0),c(0,1),type = "l",col ="grey")
points(rnaseq_validation_performance_mn_top.dec$PR_mod$specificity.micro,
       rnaseq_validation_performance_mn_top.dec$PR_mod$sensitivity.micro,
       type = "l")
points(auc_top.rnaseq["micro average","spec_at_maxp"] , auc_top.rnaseq["micro average","sens_at_maxp"] , col = "black",pch =16 ,cex =2)
dev.off()

## RNAseq confusions + level disagreement ####
pdf("figures/rnaseq_split_confusions.pdf" ,width = 10,height = 10)
colord <- c(bacteria,viral,inflammatory,"TB","KD","malaria")
colord <- colord[colord%in%colnames(rnasesq_test_confusion.low[[2]])]
confusion_heatmap(rnasesq_test_confusion.low[[2]][colord,colord],square_lab = "both",
                  margins = c(12,12,4,2))
title( main =list("RNA-Seq validation set predictions \nfor 13 disease classes" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$confident.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.low[[2]][colord,colord]),square_lab = "number")
title( main =list("Confident rnaseq predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$ambiguous.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.low[[2]][colord,colord]),square_lab = "number")
title( main =list("Ambiguous rnaseq predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)


colord <- c("bacterial","viral","JIA","TB","KD","malaria")
colord <- colord[colord%in%colnames(rnasesq_test_confusion.top[[2]])]
confusion_heatmap(rnasesq_test_confusion.top[[2]][colord,colord],level = "top_rnaseq",
                  margins = c(12,12,4,2) ,square_lab = "both",textsize = 1.5)
title( main =list("RNA-Seq validation set predictions \nfor 6 broad disease classes" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$confident.top[colord,colord],level = "top_rnaseq"  ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.top[[2]][colord,colord]),square_lab = "number")
title( main =list("Confident rnaseq predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$ambiguous.top[colord,colord],level = "top_rnaseq" ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.top[[2]][colord,colord]),square_lab = "number")
title( main =list("Ambiguous rnaseq predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
dev.off()



# AUC vs number of samples ####

pdf("figures/auc_vs_nsam.all.pdf" ,width = 15,height = 15)
par(mfrow = c(2,2))
auc_mat1 <- auc_low[!is.na(auc_low[,"CI.l"]),]
auc_mat2 <- auc_top[!is.na(auc_top[,"CI.l"]),]
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",
                 ytick = c(0.5,1,0.1),x_lab = "Number of samples",y_lab = "AUC",
                 xtick = c(0, 60,20),xbuffer = 10,ybuffer = 0.05)
for(i in rownames(auc_mat1)){
  points(auc_mat1[i,"#cases"],auc_mat1[i,"AUC"], col = colour_table[i,"col"],pch = 16,cex = 2)
  text(auc_mat1[i,"#cases"],auc_mat1[i,"AUC"],labels  = colour_table[i,"name"],pos = 1)
}
title( main =list("Microarray test set AUCROC and number of cases" ,cex = 1),line = 0.4,adj = 0)

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",
                 ytick = c(0.5,1,0.1),x_lab = "Number of samples",y_lab = "AUC",
                 xtick = c(0, 100,20),xbuffer = 10,ybuffer = 0.05)
for(i in rownames(auc_mat2)){
  points(auc_mat2[i,"#cases"],auc_mat2[i,"AUC"], col = colour_table[i,"col"],pch = 16,cex = 2)
  text(auc_mat2[i,"#cases"],auc_mat2[i,"AUC"],labels  = colour_table[i,"name"],pos = 1)
}
title( main =list("Microarray test set AUCROC and number of cases for broad disease class" ,cex = 1),line = 0.4,adj = 0)

auc_mat1 <- auc_low.rnaseq[!is.na(auc_low.rnaseq[,"CI.l"]),]
auc_mat2 <- auc_top.rnaseq[!is.na(auc_top.rnaseq[,"CI.l"]),]

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",
                 ytick = c(0.5,1,0.1),x_lab = "Number of samples",y_lab = "AUC",
                 xtick = c(0, 60,20),xbuffer = 10,ybuffer = 0.05)

for(i in rownames(auc_mat1)){
  points(auc_mat1[i,"#cases"],auc_mat1[i,"AUC"], col = colour_table[i,"col"],pch = 16,cex = 2)
  text(auc_mat1[i,"#cases"],auc_mat1[i,"AUC"],labels  = colour_table[i,"name"],pos = 1)
}
title( main =list("RNA-Seq validation set AUCROC and number of cases" ,cex = 1),line = 0.4,adj = 0)


pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",
                 ytick = c(0.5,1,0.1),x_lab = "Number of samples",y_lab = "AUC",
                 xtick = c(0, 100,20),xbuffer = 10,ybuffer = 0.05)

for(i in rownames(auc_mat2)){
  points(auc_mat2[i,"#cases"],auc_mat2[i,"AUC"], col = colour_table[i,"col"],pch = 16,cex = 2)
  text(auc_mat2[i,"#cases"],auc_mat2[i,"AUC"],labels  = colour_table[i,"name"],pos = 1)
}
title( main =list("RNA-Seq validation set AUCROC and number of cases for broad disease class" ,cex = 1),line = 0.4,adj = 0)

dev.off()





# PCA + coconut ####
library(factoextra)
ds_cols <- data.frame() 
ds_cols[unique(merged_phenos[,"experiment"]),"col"] <- rainbow(14)
ds_cols[,"names"] <- rownames(ds_cols)
ds_cols["GSE68004","col"] <- "#000000FF"
ds_cols["mega","names"] <- "GSE73464"
ds_cols["GSE38900[[1]]","names"] <- "GSE38900"
ds_cols["GSE38900[[2]]","names"] <- "GSE38900(V3/GPL6884)"

lab_pos <- data.frame()
lab_pos[rownames(ds_cols),1:9] <- NA
for(grp in rownames(ds_cols)){
  sams <- rownames(merged_phenos[merged_phenos[,"experiment"] == grp,])
  lab_pos[grp,1] <- mean(pca.uncorrected$x[sams,1])
  lab_pos[grp,2] <- mean(pca.uncorrected$x[sams,2])
  lab_pos[grp,3] <- 2
  sams <- sams[sams %in% rownames(pca.uncorrected_filtered$x)]
  if(length(sams)>0){
    lab_pos[grp,4] <- mean(pca.uncorrected_filtered$x[sams,1])
    lab_pos[grp,5] <- mean(pca.uncorrected_filtered$x[sams,2])
    lab_pos[grp,6] <- 2
  }
  sams <- sams[sams %in% rownames(pca.coconut$x)]
  if(any(sams %in% rownames(pca.coconut$x))){
    lab_pos[grp,7] <- mean(pca.coconut$x[sams,1])
    lab_pos[grp,8] <- mean(pca.coconut$x[sams,2])
    lab_pos[grp,9] <- 2
  }
}  

lab_pos[,1] <- lab_pos[,1] -40
lab_pos["GSE38900[[2]]",1] <- -500
lab_pos["GSE38900[[2]]",2] <- 60
lab_pos["GSE38900[[2]]",3] <- 4
lab_pos["GSE40396",2] <- 40
lab_pos["GSE39941",2] <- -55
lab_pos["GSE39941",1] <- 140
lab_pos["mega",1] <- 70
lab_pos["mega",2] <- 20
lab_pos["mega",3] <- 4
lab_pos["GSE68004",2] <- -30
lab_pos["GSE29366",2] <- -65
lab_pos["GSE72810",2] <- 30
lab_pos["GSE22098",2] <- -40

lab_pos[,4] <- lab_pos[,4] -20
lab_pos["GSE65391",4] <- -250
lab_pos["GSE65391",6] <- 4
lab_pos["GSE72810",5] <- 105
lab_pos["mega",5] <- -30
lab_pos["GSE42026",4] <- -60
lab_pos["GSE42026",5] <- 20
lab_pos["GSE42026",6] <- 4
lab_pos["GSE34404",5] <- 90
cexval <- 0.7

pdf("figures/pca_pre_and_post_combat.pdf" ,width = 10,height = 5)
par(mfrow = c(2,3))
plotsams <- rownames(merged_phenos)[rownames(merged_phenos) %in% rownames(pca.uncorrected$x)]
plot(pca.uncorrected$x[plotsams,c(1,2)],main ="Raw",las = 2,axes = F,
     col =add.alpha(ds_cols[merged_phenos[plotsams,"experiment"],"col"],0.4),pch =16,cex =cexval)
axis(1,tick=T,lwd =1,lwd.ticks = 1)
axis(2,tick=T,lwd =1,lwd.ticks = 1,las =2)

groups  <- as.factor(merged_phenos[plotsams,"experiment"])
df <- cbind.data.frame(groups,pca.uncorrected$x[,c(1,2)])
bb <- coord.ellipse(df,level.conf = 0.95)
for(grp in levels(groups)){
  points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =ds_cols[grp,1],lwd = 1)
  dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,1])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,2])**2
  cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
  text(lab_pos[grp,1],lab_pos[grp,2] , ds_cols[grp,"names"],pos = lab_pos[grp,3],
       col = ds_cols[grp,"col"],cex = cexval)
  points(c(cp[1],lab_pos[grp,1]) ,c(cp[2],lab_pos[grp,2]) ,
         type ="l",col = ds_cols[grp,"col"],lwd = 1) 
}

plotsams <- plotsams[plotsams%in% rownames(pca.uncorrected_filtered$x)]
plot( pca.uncorrected_filtered$x[plotsams,c(1,2)],
      main ="Raw filtered",las = 2,axes = F,col =add.alpha(ds_cols[merged_phenos[plotsams,"experiment"],"col"],0.4),pch =16,cex =cexval)
axis(1,tick=T,lwd =1,lwd.ticks = 1)
axis(2,tick=T,lwd =1,lwd.ticks = 1,las =2)
groups  <- as.factor(merged_phenos[plotsams,"experiment"])
df <- cbind.data.frame(groups,pca.uncorrected_filtered$x[,c(1,2)])
bb <- coord.ellipse(df,level.conf = 0.95)
for(grp in levels(groups)){
  points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =ds_cols[grp,1],lwd = 1)
  dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,4])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,5])**2
  cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
  text(lab_pos[grp,4],lab_pos[grp,5] , ds_cols[grp,"names"],pos = lab_pos[grp,6],
       col = ds_cols[grp,"col"],cex = cexval)
  points(c(cp[1],lab_pos[grp,4]) ,c(cp[2],lab_pos[grp,5]) ,
         type ="l",col = ds_cols[grp,"col"],lwd = 1) 
}

lab_ord <- rownames(lab_pos)[!is.na(lab_pos[,7])]
lab_ord <- lab_ord[c(2,4,6,7,8,5,9,10,3,1,12,13,11)]
lab_pos[!is.na(lab_pos[,7]),7] <- -180
lab_pos[lab_ord,8] <- seq(200,-200,-26)[1:13]

plot(pca.coconut$x[plotsams,c(1,2)],main ="Corrected",las = 2,axes = F,xlim = c(-300,200),
     col =add.alpha(ds_cols[merged_phenos[plotsams,"experiment"],"col"],0.4),pch =16,cex =cexval)
axis(1,tick=T,lwd =1,lwd.ticks = 1)
axis(2,tick=T,lwd =1,lwd.ticks = 1,las =2)
groups  <- as.factor(merged_phenos[plotsams,"experiment"])
df <- cbind.data.frame(groups,pca.coconut$x[,c(1,2)])
bb <- coord.ellipse(df,level.conf = 0.95)
for(grp in levels(groups)){
  points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =ds_cols[grp,1],lwd = 1)
  dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,7])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,8])**2
  cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
  text(lab_pos[grp,7],lab_pos[grp,8] , ds_cols[grp,"names"],pos = lab_pos[grp,9],
       col = ds_cols[grp,"col"],cex = cexval)
  points(c(cp[1,1],lab_pos[grp,7]) ,c(cp[1,2],lab_pos[grp,8]) ,
         type ="l",col = ds_cols[grp,"col"],lwd = 1) 
}

plot_scree(pca.uncorrected,main ="Raw")
plot_scree(pca.uncorrected_filtered,main ="Raw filtered")
plot_scree(pca.coconut,main = "Corrected",tickint = 1)
dev.off()

# rnaseq  pca ####
cexval <- 1
plotsams <- names(rnaseq_validation_pheno.seen)
pdf("figures/pca_rnaseq.pdf" ,width = 10,height = 10)
par(mfrow = c(2,2))
lab_pos <- data.frame()
groups <- unique(hier.rna)[c(2,1,3:6)]
lab_pos[groups,1] <- seq(-150,1200, 12)[1:6]
lab_pos[groups,2] <- seq(80,600,12)[1:6]
lab_pos[groups,3] <- 2

plot(pca.uncorrected.rnaseq$x[plotsams,c(1,2)],main ="Raw",las = 2,axes = F,xlim = c(-250,260),
     col =add.alpha(colour_table[hier.rna[rnaseq_validation_pheno[plotsams]],"col"],0.4),
     pch =16,cex =cexval,xaxs="r", yaxs="r")
axis(1,tick=T,lwd =1,lwd.ticks = 1)
axis(2,tick=T,lwd =1,lwd.ticks = 1,las =2)
groups  <- as.factor(hier.rna[rnaseq_validation_pheno[plotsams]])
df <- cbind.data.frame(groups,pca.uncorrected.rnaseq$x[plotsams,c(1,2)])
bb <- coord.ellipse(df,level.conf = 0.95)
for(grp in levels(groups)){
  points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =colour_table[grp,"col"],lwd = 1)
  dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,1])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,2])**2
  cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
  text(lab_pos[grp,1],lab_pos[grp,2] , colour_table[grp,"name"],pos = lab_pos[grp,3],
       col = colour_table[grp,"col"],cex = cexval)
  points(c(cp[1,1],lab_pos[grp,1]) ,c(cp[1,2],lab_pos[grp,2]) ,
         type ="l",col = colour_table[grp,"col"],lwd = 1) 
}

lab_pos <- data.frame()
groups <- unique(hier.rna)[c(2,1,3:6)]
lab_pos[groups,1] <- seq(-120,1200, 12)[1:6]
lab_pos[groups,2] <- seq(80,600,12)[1:6]
lab_pos[groups,3] <- 2

plot( pca.corrected.rnaseq$x[plotsams,c(1,2)],main ="Normalised",xlim = c(-200,170),
      las = 2,axes = F,col =add.alpha(colour_table[hier.rna[rnaseq_validation_pheno[plotsams]],"col"],0.4),pch =16,cex =cexval,xaxs="r", yaxs="r")
axis(1,tick=T,lwd =1,lwd.ticks = 1)
axis(2,tick=T,lwd =1,lwd.ticks = 1,las =2)
groups  <- as.factor(hier.rna[rnaseq_validation_pheno[plotsams]])
df <- cbind.data.frame(groups,pca.corrected.rnaseq$x[plotsams,c(1,2)])
bb <- coord.ellipse(df,level.conf = 0.95)
for(grp in levels(groups)){
  points(bb$res[bb$res[,1]==grp,2:3],type ="l",col =colour_table[grp,"col"],lwd = 1)
  dv <- (bb$res[bb$res[,1]==grp,2] - lab_pos[grp,1])**2 + (bb$res[bb$res[,1]==grp,3] - lab_pos[grp,2])**2
  cp <- bb$res[bb$res[,1]==grp,2:3][dv == min(dv),]
  text(lab_pos[grp,1],lab_pos[grp,2] , colour_table[grp,"name"],pos = lab_pos[grp,3],
       col = colour_table[grp,"col"],cex = cexval)
  points(c(cp[1,1],lab_pos[grp,1]) ,c(cp[1,2],lab_pos[grp,2]) ,
         type ="l",col = colour_table[grp,"col"],lwd = 1) 
  
}

plot_scree(pca.uncorrected.rnaseq,main ="Raw")
plot_scree(pca.corrected.rnaseq,main ="Normalised")
par(mfrow = c(1,1))
dev.off()


# relaxed lasso / lasso ridge ####
titlesize <- 1.6
pdf("figures/two_stage_cv.pdf" ,width = 15,height = 10)
par(mar = c(5,5,2,2))
plot_cv_twostage_bridge(relaxed_lasso_fit,type = "contour",yaxl = "phi (LASSO refit)")
title( main =list("Relaxed Lasso cross validation" ,cex =titlesize),line = 0.4,adj = 0)
plot_cv_twostage_bridge(ridge_relaxed_lasso_fit,type = "contour",yaxl = "phi (Ridge refit)")
title( main =list("Lasso-Ridge hybrid cross validation" ,cex =titlesize),line = 0.4,adj = 0)
par(mar = c(4,4,2,2))
plot_ridge_cv(mse_plot_dat.low ,lasso_fit_cv = lasso_fit_cv.low,xvar = "lambda",
                     aux_mse_pltdat =list(mse_plot_dat.low.relaxed_lasso),cexval =1,
                      mode = "LOW",xrange = c(-6,-1.5),colvec = brewer.pal(3,"Dark2"),lw = 2)
title( main =list("Tuning Lambda in Lasso and two stage methods" ,cex =titlesize),line = 0.4,adj = 0)
par(mfrow = c(1,1))
dev.off()

# error change  + cost ####
plot_errchange <- function(m1,m2,pn,seed_n,err_list,lim,title,xlabs =c(1,2),sf=3){
  pretty_plot_area(cols = c("grey90","white"),text_col = "grey30",
                 ytick = c(0,lim,0.001),x_lab = "",y_lab = paste(toupper(pn)," error"),
                 xtick = c(-1,2,1),show_x_tick = F,show_x = F)
  training_groups.loc <- colnames(err_list[[m1]][[pn]])
  pathcols <- rainbow(length(training_groups.loc))
  for(path in training_groups.loc){
    points(c(0,1),c(err_list[[m1]][[pn]][seed_n,path],
                    err_list[[m2]][[pn]][seed_n,path]),
           lwd = sf * class_weights[path,"cost"]/4 ,type = "b",
           lty=colour_table[path,"line"],
           col = colour_table[path,"col"],
           pch =colour_table[path,"shape"])
    text(-0.1 , err_list[[m1]][[pn]][seed_n,path] ,pos = 2 ,labels = colour_table[path,"name"],cex = 0.9)
    text(1.1 , err_list[[m2]][[pn]][seed_n,path] ,pos = 4 ,labels = colour_table[path,"name"],cex = 0.9)
    axis(1,at = c(-0.5,1.5) , labels = xlabs,cex =0.6,las =1)
  }

}


# unweighted error on values:
opt <- c("u","v")
lim <- 11
seed_n <- 1

err_list <- list()
for(meth in c( "M_TT", "M_TF", "M_FT", "M_FF")){
  err_list[[meth]] <- seed_performance.merged[[paste(meth,opt[1],opt[2],sep=".")]]
}
par(mfrow = c(1,2))
par(mar = c(4,4,3,2))
plot_errchange("M_FT","M_TT","fp",seed_n,err_list,lim,xlabs = c("insensitive","cost-sensitive"),
               "FP error change")
plot_errchange("M_FT","M_TT","fn",seed_n,err_list,lim, xlabs = c("insensitive","cost-sensitive"),
               "FN error change")

# plot_delta_error
  
opt <- c("u","v")
m1 <- "M_TT"
m2 <- "M_FT"
training_groups.loc <- colnames(err_list[[m1]]$fp)
ord_paths <- class_weights[,"cost"][order(class_weights[,"cost"])]
delta_errs <- c()
jitterparam <- 1
atvals <- c()
colvals <- c()
for(i in 1:18){
  # fp err
  delta_err <- err_list[[m1]][["fp"]][,names(ord_paths)[i]] - err_list[[m2]][["fp"]][,names(ord_paths)[i]]
  atvals <- c(atvals ,rep(i-0.15,10))
  delta_errs <- c(delta_errs,delta_err)  
  colvals <- c(colvals , rep("red",10))
  # fn err
  delta_err <- err_list[[m1]][["fn"]][,names(ord_paths)[i]] - err_list[[m2]][["fn"]][,names(ord_paths)[i]]
  atvals <- c(atvals ,rep(i+0.15,10))
  delta_errs <- c(delta_errs,delta_err)  
  colvals <- c(colvals , rep("blue",10))
}

pdf("figures/cost_sensitivity_error.pdf" ,width = 8,height = 10)
par(mar = c(2,4,4,2))
par(xpd = F)
plot(0,0,cex =0, axes = F, ylim=c(-5,5),xlim=c(0.5,18), ylab ="cost sensitive - cost insensitive error                     ",xlab = "")
axis(2,las = 2,cex.axis =0.5,lwd = 0,lwd.ticks = 1)
points(jitter(atvals,1.4),delta_errs,col = colvals,cex =0.5)
abline(v = seq(0.5,20.5),col = "lightgray")
abline(h = 0)
par(xpd = T)
text(1:18,rep(3.8,18),paste(colour_table[names(ord_paths),"name"]),
     cex = 0.9,srt = 90 ,adj = 0)
text(1:18,rep(3.3,18),paste(ord_paths),
     cex = 1 )
text(0.4,3.3,"Misclassification Cost", cex =0.6, adj = 1)
title( main =list("Changes in FP and FN error\n due to cost weighting in training",cex =titlesize),line = 1,adj = 0)
dev.off()

# method selection #####
# model size vs mse
plot_ms_v_mse <- function(opt,plotgroups = c("mn_TT", "mn_TF", "mn_FT", "mn_FF", "ova"), colvec = NULL , plotmacro = T){
  if( is.null(colvec)) colvec <- pallette
  ova_err <- seed_performance.merged[[paste("OVA",opt[1],opt[2],sep=".")]]$meansquared / length(seed_performance.merged[[paste("OVA",opt[1],opt[2],sep=".")]]$allerr)
  plot(0,0,cex = 0 ,ylim = c(0,1.4*max(ova_err)),
       xlim = c(0,800),axes = F,xpd=NA,
       xlab = "Signature size",ylab ="Weighted MSE",main = "",las =2)
  axis(1,tick=T,lwd =1,lwd.ticks = 1,pos = 0)
  axis(2,tick=T,lwd =1,lwd.ticks = 2,las =2,pos = 0)
  for(j in 1:length(plotgroups)){
    x.tmp <- c()
    y.tmp <- c()
    for(i in 1:10){
      ms <- mean(seed_modelsize[[i]][,plotgroups[j]])
      M <- toupper(gsub("mn_" ,"M_" ,plotgroups[j]))
      mse <- seed_performance.merged[[paste(M,opt[1],opt[2],sep=".")]]$meansquared[i] / length(seed_performance.merged[[paste(M,opt[1],opt[2],sep=".")]]$allerr)
      points(ms , mse , col = colvec[j],pch = 16)
      
      if(plotmacro){
        mse <- seed_performance.merged[[paste(M,opt[1],opt[2],sep=".")]]$macromwse[i]
        points(ms , mse , col = colvec[j],pch = 1)
      }
      
      x.tmp <- c(x.tmp , ms)
      y.tmp <- c(y.tmp , mse)
    }
  }
}


pallette <- brewer.pal(8,"Dark2")
pdf("figures/method_comparison.pdf" ,width = 10,height = 5)
par(mar = c(5, 4, 4, 8))
par(xpd = T)
opt <- c("cw","v")
plot_ms_v_mse(opt, plotmacro = F)
numsams <- length(pooled_lasso_ridges[[1]][[paste("ridge_lasso",".",opt[1],".",opt[2],sep = "")]]$allerr)
# ridge lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:length(pooled_lasso_ridges)){ 
  points(pooled_lasso_ridges[[seedvalue]]$ms , col = pallette[6],pch =16,
         pooled_lasso_ridges[[seedvalue]][[paste("ridge_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
  
}
# relaxed lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:10){ 
  points(pooled_relaxed_lasso[[seedvalue]]$ms , col =pallette[7],pch =16,
         pooled_relaxed_lasso[[seedvalue]][[paste("relaxed_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
}
# ridge ova
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:10){ 
  points(pooled_ridgeova[[seedvalue]]$ms , col = pallette[8],pch =16,
         pooled_ridgeova[[seedvalue]][[paste("ova_ridge",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
}
legend(450, 1.5, legend=c(colnames(seed_modelsize[[1]]),"Lasso+Ridge","Relaxed Lasso","ova+ridge"),
       col=c(pallette), pch = 16, cex=0.8,bty = "n",pt.cex = 2)
title( main =list("  Method performance comparison",cex =titlesize),line = 0.1,adj = 0)

opt <- c("cw","v")
plot_ms_v_mse(opt,plotgroups = c("mn_TT","ova"),colvec = pallette[c(1,2)])
# ridge lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:length(pooled_lasso_ridges)){ 
  points(pooled_lasso_ridges[[seedvalue]]$ms , col =pallette[3],pch =16,
         pooled_lasso_ridges[[seedvalue]][[paste("ridge_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
}
# relaxed lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:10){ 
  points(pooled_relaxed_lasso[[seedvalue]]$ms , col =pallette[4],pch =16,
         pooled_relaxed_lasso[[seedvalue]][[paste("relaxed_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
}

x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:10){ 
  points(pooled_ridgeova[[seedvalue]]$ms , col = pallette[5],pch =16,
         pooled_ridgeova[[seedvalue]][[paste("ova_ridge",".",opt[1],".",opt[2],sep = "")]]$meansquared/numsams)
}
legend(800, 1.5, legend=c("multinomial Lasso","OVA Lasso","Lasso+Ridge","Relaxed Lasso","OVA LASSO + ridge"),
       col=c(pallette[1:5]), pch = 16, cex=0.8,bty = "n",pt.cex = 2)
title( main =list("  Method performance comparison",cex =titlesize),line = 0.1,adj = 0)


opt <- c("cw","v")
plot_ms_v_mse(opt,plotgroups = c("mn_TT","ova"),colvec = pallette[c(1,2)],plotmacro = F)
# ridge lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:length(pooled_lasso_ridges)){ 
  points(pooled_lasso_ridges[[seedvalue]]$ms , col =pallette[3],pch =16,
         pooled_lasso_ridges[[seedvalue]][[paste("ridge_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared /numsams)
}
# relaxed lasso
x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:length(pooled_relaxed_lasso)){ 
  points(pooled_relaxed_lasso[[seedvalue]]$ms , col =pallette[4],pch =16,
         pooled_relaxed_lasso[[seedvalue]][[paste("relaxed_lasso",".",opt[1],".",opt[2],sep = "")]]$meansquared /numsams)
}

x.tmp <- c()
y.tmp <- c()
for(seedvalue in 1:10){ 
  points(pooled_ridgeova[[seedvalue]]$ms , col = pallette[5],pch =16,
         pooled_ridgeova[[seedvalue]][[paste("ova_ridge",".",opt[1],".",opt[2],sep = "")]]$meansquared /numsams)
}

legend(450, 1.5, legend=c("multinomial Lasso","OVA Lasso","Lasso+Ridge","Relaxed Lasso","OVA LASSO + ridge"),
       col=c(pallette[1:5]), pch = 16, cex=0.8,bty = "n",pt.cex = 2)
title( main =list("  Method performance comparison",cex =titlesize),line = 0.1,adj = 0)
dev.off()
par(xpd = F)




# pairs plots ####

## Microarray ####

hgt <- 17
colord <- c(bacteria,viral,inflammatory,"KD","TB","malaria")
pdf("figures/pairs_plots.low.ma.ci.pdf" ,width = hgt + ((2*hgt)/6),height = hgt)
plot_custom_pairs(inmatrix = prediction_matrices$test$full$low.dec[,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = test_phenotypes_indicator,sf = 0.5,
                  margin_val = 0.75, axcex = 0.7,
                  conf = T)
dev.off()

hgt <- 15
pdf("figures/pairs_plots.top.ma.ci.pdf" ,width = hgt + ((2*hgt)/6),height = hgt)
plot_custom_pairs(inmatrix = prediction_matrices$test$full$top.dec, 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = test_phenotypes_indicator,margin_val = 0.75, conf = T)
dev.off()



## RNA-Seq direct model transfer ####

expr.base <-  expr.allprobe.voom[probe_order, ]

sams.loc <- names(rnaseq_validation_pheno)
sams.ord <- c()
for(i in group_ord){
  sams.ord <- c(sams.ord,sams.loc[rnaseq_validation_pheno[sams.loc] == i])
}
pm <- data.frame()
for(disease in names(coef(ridge_refit))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc[probe_order,] * coef(ridge_refit.allprobe.low)[[disease]][probe_order,1]
  expr.loc[probe_order[probe_order%ni%transferprobes],] <- 0
  expr.loc <- as.data.frame(expr.loc)
  expr.loc <- expr.loc[probe_order,sams.ord]
  pm[sams.ord,disease]  <- colSums(expr.loc)
}
rnaseq_pm.low <- pm
pm <- data.frame()
for(disease in names(coef(ridge_refit.allprobe.top))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc[probe_order,] * coef(ridge_refit.allprobe.top)[[disease]][probe_order,1]
  expr.loc[probe_order[probe_order%ni%transferprobes],] <- 0
  expr.loc <- as.data.frame(expr.loc)
  expr.loc <- expr.loc[probe_order,sams.ord]
  pm[sams.ord,disease]  <- colSums(expr.loc)
}
rnaseq_pm.top <- pm

rnaseq_ind.low <- class.ind(rnaseq_validation_pheno)
rnaseq_ind.top <- class.ind(rnaseq_broad_validation_pheno)
colnames(rnaseq_ind.top)[2] <- "inflammatory"

for(i in 1:ncol(rnaseq_pm.top)){ # transforming to put on 0-1
  cp <- rnaseq_pm.top[,i]
  cp <- cp + abs(min(cp))
  rnaseq_pm.top[,i] <- cp / max(cp)
}

for(i in 1:ncol(rnaseq_pm.low)){ # transforming to put on 0-1
  cp <- rnaseq_pm.low[,i]
  cp <- cp + abs(min(cp))
  rnaseq_pm.low[,i] <- cp / max(cp)
}

rnaseq_pm.top[rnaseq_pm.top == 0] <- 0.00000001 # divide by 0 error
hgt <- 15
colord <- c("bacterial", "viral", "inflammatory", "KD", "TB","malaria")
pdf("figures/pairs_plots.top_rna.norefit.pdf" ,width = hgt + ((3*hgt)/length(colord)),height = hgt)
plot_custom_pairs(inmatrix = rnaseq_pm.top[,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = rnaseq_ind.top,margin_val = 1, conf = T)
dev.off()



## RNA-Seq (new coefficients) ####

hgt <- 17
colord <- c(bacteria,viral,inflammatory,"KD","TB","malaria")
colord <- colord[colord%in%colnames(pred.resp.refit_mn_WCI.low.dec)]
pdf("figures/pairs_plots.low_rna.pdf" ,width = hgt + ((3*hgt)/length(colord)),height = hgt)
plot_custom_pairs(inmatrix = pred.resp.refit_mn_WCI.low.dec[test_rnaseq,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,colord],
                  sf = 0.5,margin_val = 0.4,n_wide_boxplot = 3,scale_to_max = T, regular = F,axcex = 0.7)
dev.off()


hgt <- 15
colord <- c("bacterial","viral","JIA","KD","TB","malaria")
pdf("figures/pairs_plots.top_rna.pdf" ,width = hgt + ((2*hgt)/length(colord)),height = hgt)
plot_custom_pairs(inmatrix = pred.resp.refit_mn_WCI.top.dec[test_rnaseq,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,colord],
                  sf = 1,margin_val = 0.75,scale_to_max = T)
dev.off()



## with confidence intervals on ROCs #### 
hgt <- 15
colord <- c("bacterial","viral","JIA","KD","TB","malaria")
pdf("figures/pairs_plots.top_rna.ci.pdf" ,width = hgt + ((2*hgt)/length(colord)),height = hgt)
plot_custom_pairs(inmatrix = pred.resp.refit_mn_WCI.top.dec[test_rnaseq,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,colord],
                  sf = 1,margin_val = 0.75,scale_to_max = T, conf = T)
dev.off()

hgt <- 17
colord <- c(bacteria,viral,inflammatory,"KD","TB","malaria")
colord <- colord[colord%in%colnames(pred.resp.refit_mn_WCI.low.dec)]
pdf("figures/pairs_plots.low_rna.ci.pdf" ,width = hgt + ((3*hgt)/length(colord)),height = hgt)
plot_custom_pairs(inmatrix = pred.resp.refit_mn_WCI.low.dec[test_rnaseq,colord], 
                  upper_panel = upper.panel_one_thresh,
                  lower_panel = lower.panel,
                  middle_panel = middle.panel,
                  right1_panel = right1.panel,
                  right2_panel = right2.panel,
                  bottom_panel = bottom.panel,
                  indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,colord],
                  sf = 0.5,margin_val = 0.4,n_wide_boxplot = 3,scale_to_max = T,
                  regular = F,axcex = 0.7,conf = T)
dev.off()

# confusion matrices (circle size) ####
pdf("figures/confusions_circle_size.low_ma.pdf" ,width = 18,height = 15)
sf <- 2.5
layout(matrix(c(1,1,1,1,2) ,
              nrow = 1, ncol = 5, byrow = T))
colord <- c(bacteria,"TB",viral,inflammatory,"KD","malaria")

ar <- function(){
  rect(xleft = 0,ybottom = 11,xright = 7,ytop = 18,lwd = 3,col = add.alpha(colour_table["bacterial","col"],0.3),border = F)
  text(3,10.4,"Bacterial",col = add.alpha(colour_table["bacterial","col"],0.3),cex = 5)
  rect(xleft = 7,ybottom = 5,xright = 13,ytop = 11,lwd = 3,col = add.alpha(colour_table["viral","col"],0.3),border = F)
  text(10,4.4,"Viral",col = add.alpha(colour_table["viral","col"],0.3),cex = 5)
  rect(xleft = 13,ybottom = 1,xright = 17,ytop = 5,lwd = 3,col = add.alpha(colour_table["inflammatory","col"],0.3),border = F)
  text(15.3,6.4,"Inflammatory",col = add.alpha(colour_table["inflammatory","col"],0.3),cex = 4)
} 

confusion_heatmap_sizemod(confusion_microarray_test[[2]][colord,colord],square_lab = "number",
                  margins = c(16,17,4,2),sf = 5,addrects = ar,textsize = 2,sf.ax = 2.5)
title(xlab = "Predicted class",ylab = "True class",line = 14,cex.lab = 3.75)
title( main =list("Microarray test set confusion matrix" ,cex = 4),line = 0.5,adj = 0.4)
pretty_plot_area(ytick = c(0,nrow(confusion_microarray_test[[2]]),1),cols = c("white","lightgray"),
                 xtick = c(-0.05,0.37,0.01),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,margins = c(16,0.5,4,2))
text(rep(0 ,length(colord)),seq(18,1,by = -1) - 0.5, ## These are going top to bottom 
     rbind(class_weights,class_weights_toplevel)[colord,"cost"], cex = sf,adj = 0)
text(rep(0.1 ,length(colord)),seq(18,1,by = -1) - 0.5, ## These are going top to bottom 
     rowSums(confusion_microarray_test[[2]][colord,colord]), cex = sf,adj = 0)
text(rep(0.2 ,length(colord)),seq(18,1,by = -1) - 0.5,
     format(round(auc_low[colord,"sens_at_maxp"],digits = 2)), cex = sf,adj = 0)
text(rep(0.3 ,length(colord)),seq(18,1,by = -1) - 0.5,
     format(round(auc_low[colord,"spec_at_maxp"],digits = 2)), cex = sf,adj = 0)

par(xpd = T)
text(c(0,0.1,0.2,0.3)+0.02,rep(0,4), c("Cost","Number of samples","Sensitivity","Specificity"),
     cex  = sf*0.7,srt = 90 , adj = 1)
par(xpd = F)
dev.off()




pdf("figures/confusions_circle_size.top.ma.pdf" ,width = 18,height = 15)
layout(matrix(c(1,1,1,1,2) ,
              nrow = 1, ncol = 5, byrow = T))
colord <- c("bacterial","TB","viral","inflammatory","KD","malaria")
confusion_heatmap_sizemod(confusion_microarray_test_top[[2]][colord,colord],square_lab = "number",
                          margins = c(21,22,4,2),sf = 7, level = "top",textsize = 3,sf.ax = 3)
title(xlab = "Predicted class",ylab = "True class",line = 17,cex.lab = 3.75)
title( main =list("Microarray test set confusion matrix for broad disease category" ,cex =4),
       line = 0.5,adj = 0)
pretty_plot_area(ytick = c(0,length(colord),1),cols = c("white","lightgray"),
                 xtick = c(-0.05,0.45,0.01),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,margins = c(21,0.5,4,2))
text(rep(0 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rbind(class_weights,class_weights_toplevel)[colord,"cost"], cex = 2.5,adj = 0)
text(rep(0.1 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rowSums(confusion_microarray_test_top[[2]][colord,colord]), cex = 2.5,adj = 0)
text(rep(0.2 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_top[colord,"sens_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
text(rep(0.3 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_top[colord,"spec_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
par(xpd = T)
text(c(0,0.1,0.2,0.3)+0.02,rep(0,3), c("Cost","Number of samples","Sensitivity","Specificity" ),
     cex  = 3*0.7,srt = 90 , adj = 1)
par(xpd = F)
dev.off()


pdf("figures/confusions_circle_size.low.rna.pdf" ,width = 18,height = 15)
layout(matrix(c(1,1,1,1,2) ,
              nrow = 1, ncol = 5, byrow = T))
colord <- c(bacteria,"TB",viral,inflammatory,"KD","malaria")
colord <- colord[colord%in%colnames(rnasesq_test_confusion.low[[2]])]
ar <- function(){
  rect(xleft = 0,ybottom = 7,xright = 6,ytop = 13,lwd = 3,col = add.alpha(colour_table["bacterial","col"],0.3),border = F)
  text(2.8,6.4,"Bacterial",col = add.alpha(colour_table["bacterial","col"],0.3),cex = 5)
  rect(xleft = 6,ybottom = 3,xright = 10,ytop = 7,lwd = 3,col = add.alpha(colour_table["viral","col"],0.3),border = F)
  text(7,7.4,"Viral",col = add.alpha(colour_table["viral","col"],0.3),cex = 5)
  rect(xleft = 10,ybottom = 1,xright = 12,ytop = 3,lwd = 3,col = add.alpha(colour_table["inflammatory","col"],0.3),border = F)
  text(8.5,2,"Inflammatory",col = add.alpha(colour_table["inflammatory","col"],0.3),cex = 4)
} 
confusion_heatmap_sizemod(rnasesq_test_confusion.low[[2]][colord,colord],square_lab = "number",
                          margins = c(16,17,4,2),sf = 5,addrects = ar,textsize = 2,sf.ax = 2.5)
title(xlab = "Predicted class",ylab = "True class",line = 14,cex.lab = 3.75)
title( main =list("RNA-Seq validation set confusion matrix" ,cex =4),line = 0.5,adj = 0.4)
pretty_plot_area(ytick = c(0,length(colord),1),cols = c("white","lightgray"),
                 xtick = c(-0.05,0.45,0.01),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,margins = c(16,0.5,4,0.5))
text(rep(0 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rbind(class_weights,class_weights_toplevel)[colord,"cost"], cex = 2.5,adj = 0)
text(rep(0.1 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rowSums(rnasesq_test_confusion.low[[2]][colord,colord]), cex = 2.5,adj = 0)
text(rep(0.2 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_low.rnaseq[colord,"sens_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
text(rep(0.3 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_low.rnaseq[colord,"spec_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
par(xpd = T)
text(c(0,0.1,0.2,0.3)+0.02,rep(0,4), c("Cost","Number of samples","Sensitivity","Specificity" ),
     cex  = 2.1,srt = 90 , adj = 1)
par(xpd = F)
dev.off()


pdf("figures/confusions_circle_size.top.rna.pdf" ,width = 18,height = 15)
layout(matrix(c(1,1,1,1,2) ,
              nrow = 1, ncol = 5, byrow = T))
colord <- c("bacterial","TB","viral","JIA","KD","malaria")
confusion_heatmap_sizemod(rnasesq_test_confusion.top[[2]][colord,colord],square_lab = "number",
                          margins = c(21,22,8,2),sf = 7, level = "top_rna",textsize = 3,sf.ax = 3)
title(xlab = "Predicted class",ylab = "True class",line = 17,cex.lab = 3.75)
title( main =list("RNA-Seq validation set confusion matrix \nfor broad disease category" 
                  ,cex =4),line = 0.5,adj = 0)
pretty_plot_area(ytick = c(0,length(colord),1),cols = c("white","lightgray"),
                 xtick = c(-0.05,0.45,0.01),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,margins = c(21,0.5,8,2))
text(rep(0 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rbind(class_weights,class_weights_toplevel)[colord,"cost"], cex = 2.5,adj = 0)
text(rep(0.1 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
     rowSums(rnasesq_test_confusion.top[[2]][colord,colord]), cex = 2.5,adj = 0)
text(rep(0.2 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_top.rnaseq[colord,"sens_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
text(rep(0.3 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
     format(round(auc_top.rnaseq[colord,"spec_at_maxp"],digits = 2)), cex = 2.5,adj = 0)
par(xpd = T)
text(c(0,0.1,0.2,0.3)+0.02,rep(0,4), c("Cost","Number of samples","Sensitivity","Specificity" ),  
     cex  = 3,srt = 90 , adj = 1)
par(xpd = F)
par(mfrow = c(1,1))
dev.off()



# comparison with previously published signatures ####

es_pallette <- RColorBrewer::brewer.pal(n = 6,name = "Dark2")
es_rocs <- list()
pdf("figures/external_sig_comparison.dec.pdf" ,width = 15,height = 10)
par(mfrow = c(2,3))
external_auc_table <- data.frame(matrix(nrow = 1 , ncol = 5))
colnames(external_auc_table) <- c("comparison",	"AUC",	"#cases",	"CI.l",
                                  "CI.h")
external_sig_tests <- list()
cases <- names(external_perfs$`SWEENEY dn retrained`$TBvALL$vals.cases)
controls <- names(external_perfs$`SWEENEY dn retrained`$TBvALL$vals.controls)
# MCS
r1 <- roc(cases = pred.resp.refit_mn_WCI.top.dec[cases,"TB"] ,
          controls =pred.resp.refit_mn_WCI.top.dec[controls,"TB"],ci =T)
r_sp1 <- ci.sp(r1, sensitivities=seq(0, 1, .01))
external_auc_table[1,] <- c("MCS-TB-vs-ALL",r1$auc,length(cases),r1$ci[c(1,3)])
# SWEENEY
names.tmp <- names(external_perfs)[grep("SWEENEY",names(external_perfs))]
names.tmp <- names.tmp[grep("dn.log retrained",names.tmp)]
r2 <- external_perfs[[names.tmp]]$TBvALL$roc
external_sig_tests[["TBALL_MCS-TB_SWEENEY"]] <-roc.test(r1,r2,method = "bootstrap") 
r_sp2 <- ci.sp(r2, sensitivities=seq(0, 1, .01))
external_auc_table[nrow(external_auc_table)+1,] <- c("SWEENEY-TB-vs-ALL",r2$auc,length(cases),r2$ci[c(1,3)])
# NEJM
names.tmp <- names(external_perfs)[grep("NEJM",names(external_perfs))]
names.tmp <- names.tmp[grep("dn.log retrained",names.tmp)]
r3 <- external_perfs[[names.tmp]]$TBvALL$roc
external_sig_tests[["TBALL_MCS-TB_NEJM"]] <-roc.test(r1,r3,method = "bootstrap") 
r_sp3 <- ci.sp(r3, sensitivities=seq(0, 1, .01))
external_auc_table[nrow(external_auc_table)+1,] <- c("NEJM-TB-vs-ALL",r3$auc,length(cases),r3$ci[c(1,3)])

pretty_plot_area(cols = c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
# MCS
plot(add =T,r_sp1, type  = "shape",col=add.alpha(es_pallette[1],0.3))
plot(add =T,r1,col =es_pallette[1], lty = 1)
es_rocs[["MCS_TB"]] <- r_sp1
#sweeney
plot(add =T,r_sp2, type  = "shape",col=add.alpha(es_pallette[2],0.3))
plot(r2, lty =1,add =T,col =es_pallette[2])
es_rocs[["sweeney_TB"]] <- r_sp2
# NEJM - fit
plot(add =T,r_sp3, type  = "shape",col=add.alpha(es_pallette[3],0.3))
plot(r3, lty =1,add =T,col =es_pallette[3])
es_rocs[["anderson_TB"]] <- r_sp3
title( main =list("Distinguishing Tuberculosis\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)

# KD
names.tmp <- names(external_perfs)[grep("KD",names(external_perfs))]
names.tmp <- names.tmp[grep("dn.log retrained",names.tmp)]
cases <- names(external_perfs$`KD13 dn.voom.ln retrained`$KDvALL$vals.cases)
controls <- names(external_perfs$`KD13 dn.voom.ln retrained`$KDvALL$vals.controls)
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
r <- roc(cases = pred.resp.refit_mn_WCI.top.dec[cases,"KD"] , 
         controls =pred.resp.refit_mn_WCI.top.dec[controls,"KD"],ci = T)
r_sp <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp, type  = "shape",col=add.alpha(es_pallette[1],0.5))
plot(add =T ,r , col =es_pallette[1])
external_auc_table[nrow(external_auc_table)+1,] <- c("MCS-Kawasaki-vs-ALL",r$auc,length(cases),r$ci[c(1,3)])
external_sig_tests[["KDALL_MCS-KD_13GENE"]] <- roc.test(r,external_perfs[[names.tmp]]$KDvALL$roc,method = "bootstrap") 
es_rocs[["MCS_KD"]] <- r_sp
# 13 Gene
r <- external_perfs[[names.tmp]]$KDvALL$roc
r_sp2 <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp2, type  = "shape",col=add.alpha(es_pallette[6],0.5))
plot(r, lty =1,add =T , col =es_pallette[6])
external_auc_table[nrow(external_auc_table)+1,] <- c("KD13-Kawasaki-vs-ALL",r$auc,length(cases),r$ci[c(1,3)])
title( main =list("Distinguishing Kawasaki disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
es_rocs[["wright_KD"]] <- r_sp2

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "",y_lab = "",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05,
                 show_x_tick = F,show_y_tick = F,show_x = F,show_y = F)
legend(x = 1,y = 1,legend = c("Multi-class biomarker panel",
                              "TB (Sweeney et al. 2016)",
                              "TB (Anderson et al. 2014)",
                              "KD (Wright et al. 2018)",
                              "Bacterial-Viral (Herberg et al. 2016)"),
       fill = es_pallette[c(1,2,3,6,4)],cex = 1)

# BV
cases.b <- names(external_perfs$`TWOGENE dn retrained`$BvALL$vals.cases)
cases.v <- names(external_perfs$`TWOGENE dn retrained`$VvALL$vals.cases)
names.tmp <-"TWOGENE dn.log retrained"

# Bacterial vs viral
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",margins = c(8, 8, 4, 1),
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05) 
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")

# as one measure m = gradient
m <- pred.resp.refit_mn_WCI.top.dec[,"bacterial"] /  pred.resp.refit_mn_WCI.top.dec[,"viral"]
names(m) <- rownames(pred.resp.refit_mn_WCI.top.dec)
r.diag <- roc( cases = m[cases.b],  controls = m[cases.v], ci = T)

r_sp2 <- ci.sp(r.diag, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp2, type  = "shape",col=add.alpha(es_pallette[1],0.5))
plot(add =T ,r.diag,col=es_pallette[1],lty =1)
es_rocs[["MCS_BV"]] <- r_sp2
external_auc_table[nrow(external_auc_table)+1,] <- c("MCS-Bacterial-vs-Viral",
                                                     r.diag$auc,length(cases.b),
                                                     r.diag$ci[c(1,3)])
external_sig_tests[["BV_MCS_TWOGENE"]] <- roc.test(r.diag,external_perfs[[names.tmp]]$BvV$roc,
                                                   method = "bootstrap") 

# TWOGENE
r<- external_perfs[[names.tmp]]$BvV$roc
r_sp3 <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp3, type  = "shape",col=add.alpha(es_pallette[4],0.5))
plot(r, lty =1,add =T,col =es_pallette[4])
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENE-Bacterial-vs-Viral",r$auc,length(cases.b),r$ci[c(1,3)])
title( main =list("Distinguishing bacterial from viral illness" ,cex = titlesize),line = 0.2,adj = 0)
es_rocs[["herberg_BV"]] <- r_sp3

# viral vs all
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
controls <- names(external_perfs$`TWOGENE dn retrained`$VvALL$vals.controls)
r <- roc(cases = pred.resp.refit_mn_WCI.top.dec[cases.v,"viral"] ,
         controls =pred.resp.refit_mn_WCI.top.dec[controls,"viral"],ci =T)
r_sp <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp, type  = "shape",col=add.alpha(es_pallette[1],0.5))
plot(add =T ,r,col=es_pallette[1],lty =1)
external_auc_table[nrow(external_auc_table)+1,] <- c("MCS-Viral-vs-ALL",r$auc,length(cases.v),r$ci[c(1,3)])
external_sig_tests[["VALL_MCS-V_TWOGENE"]] <- roc.test(r,external_perfs[[names.tmp]]$VvALL$roc,method = "bootstrap") 
es_rocs[["MCS_VALL"]] <- r_sp

#TWOGENE
r <- external_perfs[[names.tmp]]$VvALL$roc
r_sp2 <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp2, type  = "shape",col=add.alpha(es_pallette[4],0.5))
plot(r, lty =1,add =T,col =es_pallette[4])
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENE-Viral-vs-ALL",r$auc,length(cases.v),r$ci[c(1,3)])
title( main =list("Distinguishing viral illness\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
es_rocs[["herberg_VALL"]] <- r_sp2

# bacterial vs all
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")

controls <- names(external_perfs$`TWOGENE dn retrained`$BvALL$vals.controls)
r <- roc(cases = pred.resp.refit_mn_WCI.top.dec[cases.b,"bacterial"] ,
         controls =pred.resp.refit_mn_WCI.top.dec[controls,"bacterial"],ci =T)
r_sp <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp, type  = "shape",col=add.alpha(es_pallette[1],0.5))
plot(add =T ,r,col=es_pallette[1],lty =1)
external_auc_table[nrow(external_auc_table)+1,] <- c("MCS-Bacterial-vs-ALL",r$auc,length(cases.b),r$ci[c(1,3)])
external_sig_tests[["BALL_MCS-B_TWOGENE"]] <- roc.test(r,external_perfs[[names.tmp]]$BvALL$roc,method = "bootstrap") 
es_rocs[["MCS_BALL"]] <- r_sp

#TWOGENE
r <- external_perfs[[names.tmp]]$BvALL$roc
r_sp2 <- ci.sp(r, sensitivities=seq(0, 1, .01))
plot(add =T,r_sp2, type  = "shape",col=add.alpha(es_pallette[4],0.5))
plot(r, lty =1,add =T,col =es_pallette[4])
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENE-Bacterial-vs-ALL",r$auc,length(cases.b),r$ci[c(1,3)])
title( main =list("Distinguishing bacterial illness\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
es_rocs[["herberg_BALL"]] <- r_sp2

dev.off()


## broken into multiple plots  ####
pdf("figures/external_sig_comparison.dec.pdf" ,width = 15,height = 20)
par(mfrow = c(4,3))
# TB
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Multiclass signature distinguishing TB disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$MCS_TB, type  = "shape",col=add.alpha(es_pallette[1],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Sweeney et al. 2016 signature distinguishing TB disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$sweeney_TB, type  = "shape",col=add.alpha(es_pallette[2],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Anderson et al. 2014 signature distinguishing TB disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$anderson_TB, type  = "shape",col=add.alpha(es_pallette[3],0.5))

# BV
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Multiclass signature distinguishing bacterial \nfrom viral infection" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$MCS_BV, type  = "shape",col=add.alpha(es_pallette[1],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Herberg et al. 2016 signature distinguishing bacterial \nfrom viral infection" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$herberg_BV, type  = "shape",col=add.alpha(es_pallette[4],0.5))
#
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Multiclass signature distinguishing bacterial infection\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$MCS_BALL, type  = "shape",col=add.alpha(es_pallette[1],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Herberg et al. 2016 signature distinguishing bacterial infection\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$herberg_BALL, type  = "shape",col=add.alpha(es_pallette[4],0.5))
#
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Multiclass signature distinguishing viral infection\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$MCS_VALL, type  = "shape",col=add.alpha(es_pallette[1],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Herberg et al. 2016 signature distinguishing viral infection\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$herberg_VALL, type  = "shape",col=add.alpha(es_pallette[4],0.5))

# KD
pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Multiclass signature distinguishing Kawasaki disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$MCS_KD, type  = "shape",col=add.alpha(es_pallette[1],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                 xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
points(c(0,1),c(1,0), lwd = 4, col = "grey90", type = "l")
title( main =list("Wright et al. 2018 signature distinguishing Kawasaki disease\nfrom other febrile illness" ,cex = titlesize),line = 0.2,adj = 0)
plot(add =T,es_rocs$wright_KD, type  = "shape",col=add.alpha(es_pallette[6],0.5))

pretty_plot_area(cols =  c("white","grey90"),text_col = "grey30",margins = c(8, 8, 4, 1),
                 ytick = c(0,1,0.2),x_lab = "",y_lab = "",
                 xtick = c(0,1,0.2),xbuffer = 0.05,ybuffer = 0.05,
                 show_x_tick = F,show_y_tick = F,show_x = F,show_y = F)
legend(x = 0,y = 0.5,legend = c("Multi-class biomarker panel",
                                "TB (Sweeney et al. 2016)",
                                "TB (Anderson et al. 2014)",
                                "KD (Wright et al. 2018)",
                                "Bacterial-Viral (Herberg et al. 2016)"),
       fill = es_pallette[c(1,2,3,6,4)],cex = 1.5)
dev.off()



# add DRS


# SWEENEY DRS
r <- external_perfs[["NEJMTB dn.log givencoefs"]]$TBvALL$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("SWEENEY_givencoef-TB-vs-ALL",r$auc,length(cases),r$ci[c(1,3)])

# NEJM drs
r <- external_perfs[["NEJMTB dn.log givencoefs"]]$TBvALL$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("NEJMgivencoef-TB-vs-ALL",r$auc,length(cases),r$ci[c(1,3)])

# KD13 DRS
r <- external_perfs[["KD13 dn.log givencoefs"]]$KDvALL$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("KD13givencoef-Kawasaki-vs-ALL",r$auc,length(cases),r$ci[c(1,3)])

# twogene
r<- external_perfs[["TWOGENE dn.log givencoefs"]]$BvV$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENEgivencoef-Bacterial-vs-Viral",r$auc,length(cases.v),r$ci[c(1,3)])
r<- external_perfs[["TWOGENE dn.log givencoefs"]]$BvALL$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENEgivencoef-Bacterial-vs-ALL",r$auc,length(cases.v),r$ci[c(1,3)])
r<- external_perfs[["TWOGENE dn.log givencoefs"]]$VvALL$roc
external_auc_table[nrow(external_auc_table)+1,] <- c("TWOGENEgivencoef-Viral-vs-ALL",r$auc,length(cases.v),r$ci[c(1,3)])

write.table(external_auc_table,file = "figures/external_auc_table.dec.tab",sep = "\t",quote = F,row.names = F)

external_sig_table  <- data.frame(matrix(nrow = length(external_sig_tests),ncol = 6))
for(i in 1:length(external_sig_tests)){
  external_sig_table[i,] <- c(names(external_sig_tests)[i],
                              external_sig_tests[[i]]$estimate,
                              external_sig_tests[[i]]$statistic,
                              external_sig_tests[[i]]$p.value,
                              external_sig_tests[[i]]$method)
}
colnames(external_sig_table) <- c("comparison","AUC1","AUC2","Z","p","method")

write.table(external_sig_table,file = "figures/external_sig_comparison.dec.tab",sep = "\t",quote = F,row.names = F)



# parallel plot ####


pdf("figures//parallel_p_ma.pdf")
for(i in c(bacteria,viral,inflammatory,"TB","KD","malaria") ){
  sam_list <- phenotypes[phenotypes[,"trte"] == "test"  & phenotypes[,"group"]==i,1]
  plot_parallel(prediction.matrix = prediction_matrices$test$full,sam_list = sam_list, main = i )
}
dev.off()


pdf("figures/parallel_p_ma_merged.pdf", paper = "a4r",width = 15, height = 16)
par(mfrow = c(1,3), mar = c(0,0,0,0))
for(i in c("bacteria","viral","inflammatory","TB","KD","malaria") ){
  if(i == "bacteria"){
    sam_list <- c()
    for(j in bacteria){
      sam_list <- c(sam_list, "N", phenotypes[phenotypes[,"trte"] == "test"  & phenotypes[,"group"]==j,1] )
    }
    sam_list <- sam_list[2:length(sam_list)]
  }else if(i == "viral"){
    sam_list <- c()
    for(j in viral){
      sam_list <- c(sam_list, "N", phenotypes[phenotypes[,"trte"] == "test"  & phenotypes[,"group"]==j,1] )
    }
    sam_list <- sam_list[2:length(sam_list)]
  }else if(i == "inflammatory"){
    sam_list <- c()
    for(j in inflammatory){
      sam_list <- c(sam_list, "N", phenotypes[phenotypes[,"trte"] == "test"  & phenotypes[,"group"]==j,1] )
    }
    sam_list <- sam_list[2:length(sam_list)]
  }else{
    sam_list <- phenotypes[phenotypes[,"trte"] == "test"  & phenotypes[,"group"]==i,1]  
  }
  plot_parallel(prediction_matrices$test$full,sam_list = sam_list, main = colour_table[i,"name"] ,multiplot = T)
}
dev.off()


# RNASeq

table(rnaseq_validation_pheno[test_rnaseq])

rnaseq_tesst_pred <- list(low = pred.resp.refit_mn_WCI.low,
                          low.dec = pred.resp.refit_mn_WCI.low.dec,
                          top = pred.resp.refit_mn_WCI.top,
                          top.dec = pred.resp.refit_mn_WCI.top.dec)


pdf("figures/parallel_p_rnaseq.pdf")
for(i in c(bacteria,viral,inflammatory,"TB","KD","malaria") ){
  if(i %in% colnames(rnaseq_tesst_pred$low)){
    sam_list <- test_rnaseq[rnaseq_validation_pheno[test_rnaseq] ==i]
    plot_parallel(rnaseq_tesst_pred,sam_list = sam_list, main = i )  
  }
}
dev.off()




pdf("figures/parallel_p_rnaseq_merged.pdf", paper = "a4r",width = 15, height = 16)
par(mfrow = c(1,3), mar = c(0,0,0,0))
for(i in c("bacteria","viral","inflammatory","TB","KD","malaria") ){
  
  if(i == "bacteria"){
    sam_list <- c()
    for(j in bacteria){
      sam_list <- c(sam_list, "N",  test_rnaseq[rnaseq_validation_pheno[test_rnaseq] ==j])
    }
    sam_list <- sam_list[2:length(sam_list)]
  }else if(i == "viral"){
    sam_list <- c()
    for(j in viral){
      sam_list <- c(sam_list, "N", test_rnaseq[rnaseq_validation_pheno[test_rnaseq] ==j])
    }
    sam_list <- sam_list[2:length(sam_list)]
  }else if(i == "inflammatory"){
    sam_list <- c()
    for(j in inflammatory){
      sam_list <- c(sam_list, "N", test_rnaseq[rnaseq_validation_pheno[test_rnaseq] ==j])
    }
    while(sam_list[1] == "N"){
      sam_list <- sam_list[2:length(sam_list)]  
    }
    while(sam_list[length(sam_list)] == "N"){
      sam_list <- sam_list[1:(length(sam_list)-1)]  
    }
    
  }else{
    sam_list <-   test_rnaseq[rnaseq_validation_pheno[test_rnaseq] ==i]
  }
  print(table(rnaseq_validation_pheno[sam_list]))
  if(i %in% c("TB","malaria")){
    plot_parallel(rnaseq_tesst_pred,sam_list = sam_list, main = colour_table[i,"name"] ,multiplot = T,axis_buffer = 0.1)
  }else{
    plot_parallel(rnaseq_tesst_pred,sam_list = sam_list, main = colour_table[i,"name"] ,multiplot = T)  
  }
  
}
dev.off()



# Level disagreement ####

pdf("figures/level_disagreement_confusion.pdf" ,width = 23,height = 15)

par(mfrow = c(2,3))
colord <- c(bacteria,viral,inflammatory,"TB","KD","malaria")
confusion_heatmap(confusion_microarray_test[[2]][colord,colord])
title( main =list("all microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)

confusion_heatmap(microarray_disagreement$confident.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test[[2]][colord,colord]),square_lab = "number")
title( main =list("confident microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)

confusion_heatmap(microarray_disagreement$ambiguous.low[colord,colord],scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test[[2]][colord,colord]) ,square_lab = "number")
title( main =list("ambiguous microarray predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)


colord <- c("bacterial","viral","inflammatory","TB","KD","malaria")
confusion_heatmap(confusion_microarray_test_top[[2]][colord,colord],level = "top" )
title( main =list("all microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)

confusion_heatmap(microarray_disagreement$confident.top[colord,colord],level = "top"  ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test_top[[2]][colord,colord]),square_lab = "number")
title( main =list("confident microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)

confusion_heatmap(microarray_disagreement$ambiguous.top[colord,colord],level = "top" ,scale = "",
                  class_abundance_vector =  rowSums(confusion_microarray_test_top[[2]][colord,colord]),square_lab = "number")
title( main =list("ambiguous microarray predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
##
par(mfrow = c(2,3))
colord <- c(bacteria,viral,inflammatory,"TB","KD","malaria")
colord <- colord[colord%in%colnames(rnasesq_test_confusion.low[[2]])]
confusion_heatmap(rnasesq_test_confusion.low[[2]][colord,colord])
title( main =list("all rnaseq predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$confident.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.low[[2]][colord,colord]),square_lab = "number")
title( main =list("confident rnaseq predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$ambiguous.low[colord,colord] ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.low[[2]][colord,colord]),square_lab = "number")
title( main =list("ambiguous rnaseq predictions lowlevel" ,cex = titlesize),line = 0.2,adj = 0)


colord <- c("bacterial","viral","JIA","TB","KD","malaria")
colord <- colord[colord%in%colnames(rnasesq_test_confusion.top[[2]])]
confusion_heatmap(rnasesq_test_confusion.top[[2]][colord,colord],level = "top_rnaseq" )
title( main =list("all rnaseq predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$confident.top[colord,colord],level = "top_rnaseq"  ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.top[[2]][colord,colord]),square_lab = "number")
title( main =list("confident rnaseq predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
confusion_heatmap(rnaseq_disagreement$ambiguous.top[colord,colord],level = "top_rnaseq" ,scale = "",
                  class_abundance_vector =  rowSums(rnasesq_test_confusion.top[[2]][colord,colord]),square_lab = "number")
title( main =list("ambiguous rnaseq predictions top level" ,cex = titlesize),line = 0.2,adj = 0)
par(mfrow = c(1,1))
dev.off()



# Heatmaps ####

pal <- function(x, colors=c("blue","white","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(0,1, length.out=colsteps)) ] )
}
group_ord <- c(bacteria , viral , inflammatory , "KD" ,"TB","malaria")
sams.loc <- phenotypes[phenotypes[,"trte"] == "test",1]
sams.ord <- c()
for(i in group_ord){
  sams.ord <- c(sams.ord,sams.loc[phenotypes[sams.loc,"group"] == i])
}

## microarray test set ####

# cluster genes prior to plotting
expr.base <- expression_filt[rownames(coef(ridge_fit_final.low)[[1]])[2:162],sams.ord]
d <- dist(scale(expr.base))
clust <- hclust(d)
probe_order <- clust$labels[clust$order]
dend.genes <- as.dendrogram(clust)

# cluster samples
d <- dist(scale(t(expr.base)))
clust <- hclust(d)
sample_order.clust <- clust$labels[clust$order]
dend.samples <- as.dendrogram(clust)

probe_conv.tmp <- probe_conv
probe_conv.tmp[,"probe"] <- rownames(probe_conv)
probe_labels <- probe_conv.tmp[probe_order,c("Symbol")]
probe_labels[is.na(probe_labels)] <- ""
probe_labels[probe_labels == ""] <- probe_order[probe_labels == ""]

# Mulitplying gene expression by the coefficients for each disease 
pdf("figures/microarray_heatmap_multiplycoef.pdf" ,width = 14,height = 14)
for(disease in names(coef(ridge_fit_final.low))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(ridge_fit_final.low)[[disease]][2:162]
  expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
  expr.loc <- expr.loc[probe_order,]
  pretty_plot_area(ytick = c(-10,nrow(expr.loc),10),cols = c("white","lightgray"),
                   xtick = c(0,ncol(expr.loc) + 20,10),show_x_tick = F,show_y_tick = F,
                   show_x = F,show_y = F,plot_outside_margin = T)
  for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = 0 , ybottom = -10, col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
  for(i in 1:nrow(expr.loc)){
    for(j in 1:ncol(expr.loc)){
      rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F, lwd = 0)
    }
  }
  title(main = colour_table[disease,"name"])
  
  tsams <- (1:ncol(expr.loc))[phenotypes[colnames(expr.loc),"group"] == disease]
  rect(xleft = min(tsams), xright = max(tsams), ybottom = 0, ytop = 161,density = 0 , col = "red")
  points(1:ncol(expr.loc),(prediction_matrices$test$full$low[sams.ord,disease]*10) - 20,type = "l", col = "grey")
  ind <-  phenotypes[colnames(expr.loc),"group"] == disease
  vals <- colSums(expr.loc)
  vals <- vals - min(vals)
  vals <- vals / max(vals)
  points(1:ncol(expr.loc),(vals*10) -20,type = "l")
  text( rep(0,161) , (1:161)-0.5 ,labels = probe_labels,tick = F, cex = 0.4 ,adj = 1)
  # coefficient line
  probecoef <- coef(ridge_fit_final.low)[[disease]][probe_order,1]
  probecoef <- probecoef * 20
  points(ncol(expr.loc) + 10 + probecoef,1:nrow(expr.loc),type = "l")
  points(c(ncol(expr.loc) + 10,ncol(expr.loc) + 10 ), c(0,nrow(expr.loc))-0.5, typ= "l", col = "grey")
}
dev.off()


pdf("figures/microarray_heatmap_nocoef.pdf" ,width = 14,height = 14)
# non-coef plot 
layout(matrix(1:2,nrow=2), width = c(10,10),height = c(10,1))
expr.loc <- expr.base
expr.loc <- t(scale(t(expr.loc)))
expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
expr.loc <- expr.loc[probe_order,sample_order.clust]

pretty_plot_area(ytick = c(-10,nrow(expr.loc),10),cols = c("white","lightgray"),
                 xtick = c(0,ncol(expr.loc) + 20,10),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,plot_outside_margin = T)
for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = 0 , ybottom = -10, col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
for(i in 1:nrow(expr.loc)){
  for(j in 1:ncol(expr.loc)){
    rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F, lwd = 0)
  }
}
title(main = "without coefficients")
text( rep(0,161) , (1:161)-0.5 ,labels = probe_labels, cex = 0.4 ,adj = 1)
dev.off()



pdf("figures/microarray_heatmap_dend.pdf" ,width = 14,height = 14)
layout(matrix(1:4,nrow=2, byrow = T), width = c(2,10),height = c(2,10))
plot(1,1, col = "white", axes = F,xlab ="",ylab = "")
par(mar = c(0,0,3,0))
par(xpd = T)
plot(dend.samples , horiz = F, axes = F, leaflab  = "none", frame.plot = F,
     xlim = c(-4,322.5),
     ylim =c(3,40))
par(mar = c(0,3,0,0))
plot(dend.genes , horiz = T, axes = F, leaflab  = "none", frame.plot = F,
     ylim = c(4.5,159.3), xlim = c(75,5))
expr.loc <- expr.base[probe_order , sample_order.clust]
expr.loc <- t(scale(t(expr.loc)))
expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
expr.loc <- expr.loc[probe_order,sample_order.clust]
pretty_plot_area(ytick = c(0,nrow(expr.loc)+4,5),cols = c("white","lightgray"),
                 xtick = c(0,ncol(expr.loc) + 20,10),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,plot_outside_margin = T,margins = c(1,3.5,0,1))
for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = nrow(expr.loc) + 5 ,border = F,lwd = 0,
                                ybottom = nrow(expr.loc), col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
for(i in 1:nrow(expr.loc)){
  for(j in 1:ncol(expr.loc)){
    rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F, lwd = 0)
  }
}
text( rep(-8,161) , (1:161)-0.5 ,labels = probe_labels, cex = 0.35 ,adj = 0.5)
dev.off()


pdf("figures/microarray_heatmap_genedend.pdf" ,width = 20,height = 10)
layout(matrix(1:2, byrow = T, nrow = 1), width = c(2,10))
par(xpd = T)
par(mar = c(0,4,4,0))
plot(dend.genes , horiz = T, axes = F, leaflab  = "none", frame.plot = F,
     ylim = c(3.2,166.6), xlim = c(75,2))
expr.loc <- expr.base[probe_order , sams.ord]
expr.loc <- t(scale(t(expr.loc)))
expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
expr.loc <- expr.loc[probe_order,sams.ord]
pretty_plot_area(ytick = c(0,nrow(expr.loc)+4,5),cols = c("white","lightgray"),
                 xtick = c(0,ncol(expr.loc) + 20,10),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,plot_outside_margin = T,margins = c(1,3.5,6,1))
for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = nrow(expr.loc) + 5 ,border = F,lwd = 0,
                                ybottom = nrow(expr.loc), col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
for(i in 1:nrow(expr.loc)){
  for(j in 1:ncol(expr.loc)){
    rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F, lwd = 0)
  }
}
text( rep(-8,161) , (1:161)-0.5 ,labels = probe_labels, cex = 0.35 ,adj = 0.5)
for(group in group_ord){
  ind <- phenotypes[colnames(expr.loc),"group"] == group
  xrange <- (1:length(ind))[ind]
  xpos <- mean(xrange)
  text(xpos,nrow(expr.loc)+ 6,colour_table[group,"name"], srt = 60,adj = 0)
  rect(xleft = min(xrange) -1 , xright = max(xrange), ytop = nrow(expr.loc) + 5 ,border = T,lwd = 1,
       ybottom = nrow(expr.loc))
}
dev.off()

hier.ma["TB"] <- "TB"
hier.ma["KD"] <-"KD"


pdf("figures/microarray_heatmap_top.pdf" ,width = 14,height = 14)
for(disease in names(coef(ridge_fit_final.top))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(ridge_fit_final.top)[[disease]][2:162]
  expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
  expr.loc <- expr.loc[probe_order,]
  
  pretty_plot_area(ytick = c(-10,nrow(expr.loc),10),cols = c("white","lightgray"),
                   xtick = c(0,ncol(expr.loc)+20,10),show_x_tick = F,show_y_tick = F,
                   show_x = F,show_y = F,plot_outside_margin = T)
  for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = 0 , ybottom = -10, col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
  for(i in 1:nrow(expr.loc)){
    for(j in 1:ncol(expr.loc)){
      
      rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F,lwd = 0)
    }
  }
  title(main = colour_table[disease,"name"])
  # trueclass red box
  tsams <- (1:ncol(expr.loc))[hier.ma[phenotypes[colnames(expr.loc),"group"]] == disease]
  rect(xleft = min(tsams), xright = max(tsams), ybottom = 0, ytop = 161,density = 0 , col = "red")
  points(1:ncol(expr.loc),(prediction_matrices$test$full$top[sams.ord,disease]*10) - 20,type = "l", col = "grey")
  
  # colsums line
  ind <-  hier.ma[phenotypes[colnames(expr.loc),"group"]] == disease
  vals <- colSums(expr.loc)
  vals <- vals - min(vals)
  vals <- vals / max(vals)
  points(1:ncol(expr.loc),(vals*10) -20,type = "l")
  
  # gene labels
  text( rep(0,161) , (1:161)-0.5 ,labels = probe_labels,tick = F, cex = 0.4 ,adj = 1)
  
  # coefficient line
  probecoef <- coef(ridge_fit_final.top)[[disease]][probe_order,1]
  probecoef <- probecoef * 20
  points(ncol(expr.loc) + 10 + probecoef,1:nrow(expr.loc),type = "l")
  points(c(ncol(expr.loc) + 10,ncol(expr.loc) + 10 ), c(0,nrow(expr.loc))-0.5, typ= "l", col = "grey")
}
dev.off()


## microarray training set ####

group_ord <- c(bacteria , viral , inflammatory , "KD" ,"TB","malaria")
sams.loc <- phenotypes[phenotypes[,"trte"] == "train",1]
sams.ord <- c()
for(i in group_ord){
  sams.ord <- c(sams.ord,sams.loc[phenotypes[sams.loc,"group"] == i])
}

expr.base <- expression_filt[rownames(coef(ridge_fit_final.low)[[1]])[2:162],sams.ord]
d <- dist(scale(expr.base))
clust <- hclust(d)
probe_order <- clust$labels[clust$order]
dend.genes <- as.dendrogram(clust)

# cluster samples
d <- dist(scale(t(expr.base)))
clust <- hclust(d)
sample_order.clust <- clust$labels[clust$order]
dend.samples <- as.dendrogram(clust)

pdf("figures/microarray_heatmap_genedend_train.pdf" ,width = 20,height = 10)
layout(matrix(1:2, byrow = T, nrow = 1), width = c(2,10))
par(xpd = T)
par(mar = c(0,4,4,0))
plot(dend.genes , horiz = T, axes = F, leaflab  = "none", frame.plot = F,
     ylim = c(3.2,166.6), xlim = c(75,2))
expr.loc <- expr.base[probe_order , sams.ord]
expr.loc <- t(scale(t(expr.loc)))
expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
expr.loc <- expr.loc[probe_order,sams.ord]
pretty_plot_area(ytick = c(0,nrow(expr.loc)+4,5),cols = c("white","lightgray"),
                 xtick = c(0,ncol(expr.loc) + 20,10),show_x_tick = F,show_y_tick = F,
                 show_x = F,show_y = F,plot_outside_margin = T,margins = c(1,3.5,6,1))
for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = nrow(expr.loc) + 5 ,border = NA,lwd = 0,
                                ybottom = nrow(expr.loc), col = colour_table[phenotypes[colnames(expr.loc)[j],"group"],1])
for(i in 1:nrow(expr.loc)){
  for(j in 1:ncol(expr.loc)){
    rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = NA, lwd = 0)
  }
}
text( rep(-12,161) , (1:161)-0.5 ,labels = probe_labels, cex = 0.35 ,adj = 0.5)

for(group in group_ord){
  ind <- phenotypes[colnames(expr.loc),"group"] == group
  xrange <- (1:length(ind))[ind]
  xpos <- mean(xrange)
  text(xpos,nrow(expr.loc)+ 6,colour_table[group,"name"], srt = 60,adj = 0)
  rect(xleft = min(xrange) -1 , xright = max(xrange), ytop = nrow(expr.loc) + 5 ,border = T,lwd = 1,
       ybottom = nrow(expr.loc))
}
dev.off()

## RNA-Seq ####

expr.base <-  expr.voom[,]
probe_labels <- probe_conv.tmp[probe_order,c("Symbol")]
probe_labels[is.na(probe_labels)] <- ""
probe_labels[probe_labels == ""] <- probe_order[probe_labels == ""]


sams.loc <- names(rnaseq_validation_pheno)
sams.ord <- c()
for(i in group_ord){
  sams.ord <- c(sams.ord,sams.loc[rnaseq_validation_pheno[sams.loc] == i])
}

pdf("figures/expr_heatmap_rnaseq_top.pdf" ,width = 14,height = 14)
pm <- data.frame()
for(disease in names(coef(rnaseq_pred.allprobe.top))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(rnaseq_pred.allprobe.top)[[disease]][rownames(expr.loc),1]
  expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
  expr.loc <- as.data.frame(expr.loc)
  expr.loc[probe_order[probe_order%ni%rownames(expr.loc)],] <- NA
  expr.loc <- expr.loc[probe_order,sams.ord]
  
  pretty_plot_area(ytick = c(-10,nrow(expr.loc),10),cols = c("white","lightgray"),
                   xtick = c(0,ncol(expr.loc)+20,10),show_x_tick = F,show_y_tick = F,
                   show_x = F,show_y = F,plot_outside_margin = T)
  for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = 0 , ybottom = -10,
                                  col = colour_table[rnaseq_validation_pheno[colnames(expr.loc)[j]],1])
  for(i in 1:nrow(expr.loc)){
    if(is.na(expr.loc[i,1])){
      rect(xleft = 0 , xright = ncol(expr.loc), ybottom = i - 1 , ytop = i, col = "grey", border = F,lwd = 0)
    }else{
      for(j in 1:ncol(expr.loc)){
        rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F,lwd = 0)
      }
    }
    
  }
  title(main = colour_table[disease,"name"])
  
  
  # trueclass red box
  tsams <- (1:ncol(expr.loc))[hier.ma[rnaseq_validation_pheno[colnames(expr.loc)]] == disease]
  rect(xleft = min(tsams), xright = max(tsams), ybottom = 0, ytop = 161,density = 0 , col = "red")
  points(1:ncol(expr.loc),(rnaseq_pred.allprobe.top[sams.ord,disease]*10) - 20,type = "l", col = "grey")
  
  
  # colsums line
  ind <-  hier.ma[rnaseq_validation_pheno[colnames(expr.loc)]] == disease
  vals <- colSums(expr.loc, na.rm = T)
  vals <- vals - min(vals)
  vals <- vals / max(vals)
  points(1:ncol(expr.loc),(vals*10) -20,type = "l")
  pm[sams.ord,disease] <- vals
  # gene labels
  text( rep(0,161) , (1:161)-0.5 ,labels = probe_labels,cex = 0.4 ,adj = 1)
  
  
  # coefficient line
  probecoef <- coef(rnaseq_pred.allprobe.top)[[disease]][,1]
  probecoef[rownames(expr.loc)[rownames(expr.loc) %ni% names(probecoef)]] <- 0
  probecoef <- probecoef[probe_order]
  probecoef <- probecoef * 20
  points(ncol(expr.loc) + 10 + probecoef,1:nrow(expr.loc),type = "l")
  points(c(ncol(expr.loc) + 10,ncol(expr.loc) + 10 ), c(0,nrow(expr.loc))-0.5, typ= "l", col = "grey")
  
}
dev.off()




pdf("figures/expr_heatmap_rnaseq_low.pdf" ,width = 14,height = 14)
pm <- data.frame()
for(disease in names(coef(rnaseq_pred.allprobe.low))){
  expr.loc <- expr.base
  expr.loc <- t(scale(t(expr.loc)))
  expr.loc <- expr.loc * coef(rnaseq_pred.allprobe.low)[[disease]][rownames(expr.loc),1]
  expr.loc <- expr.loc / (2*max(abs(expr.loc))) + 0.5
  expr.loc <- as.data.frame(expr.loc)
  expr.loc[probe_order[probe_order%ni%rownames(expr.loc)],] <- NA
  expr.loc <- expr.loc[probe_order,sams.ord]
  
  pretty_plot_area(ytick = c(-10,nrow(expr.loc),10),cols = c("white","lightgray"),
                   xtick = c(0,ncol(expr.loc)+20,10),show_x_tick = F,show_y_tick = F,
                   show_x = F,show_y = F,plot_outside_margin = T)
  for(j in 1:ncol(expr.loc)) rect(xleft = j - 1 , xright = j, ytop = 0 , ybottom = -10,
                                  col = colour_table[rnaseq_validation_pheno[colnames(expr.loc)[j]],1])
  for(i in 1:nrow(expr.loc)){
    if(is.na(expr.loc[i,1])){
      rect(xleft = 0 , xright = ncol(expr.loc), ybottom = i - 1 , ytop = i, col = "grey", border = F,lwd = 0)
    }else{
      for(j in 1:ncol(expr.loc)){
        rect(xleft = j - 1 , xright = j, ybottom = i - 1 , ytop = i, col = pal(expr.loc[i,j]), border = F,lwd = 0)
      }
    }
    
  }
  title(main = colour_table[disease,"name"])
  
  # trueclass red box
  tsams <- (1:ncol(expr.loc))[rnaseq_validation_pheno[colnames(expr.loc)] == disease]
  rect(xleft = min(tsams), xright = max(tsams), ybottom = 0, ytop = 161,density = 0 , col = "red")
  points(1:ncol(expr.loc),(rnaseq_pred.allprobe.low[sams.ord,disease]*10) - 20,type = "l", col = "grey")
  
  # colsums line
  vals <- colSums(expr.loc, na.rm = T)
  vals <- vals - min(vals)
  vals <- vals / max(vals)
  points(1:ncol(expr.loc),(vals*10) -20,type = "l")
  pm[sams.ord, disease] <- vals
  
  # gene labels
  text( rep(0,161) , (1:161)-0.5 ,labels = probe_labels,cex = 0.4 ,adj = 1)
  
  # coefficient line
  probecoef <- coef(rnaseq_pred.allprobe.low)[[disease]][,1]
  probecoef[rownames(expr.loc)[rownames(expr.loc) %ni% names(probecoef)]] <- 0
  probecoef <- probecoef[probe_order]
  probecoef <- probecoef * 20
  points(ncol(expr.loc) + 10 + probecoef,1:nrow(expr.loc),type = "l")
  points(c(ncol(expr.loc) + 10,ncol(expr.loc) + 10 ), c(0,nrow(expr.loc))-0.5, typ= "l", col = "grey")
}
dev.off()





# prevalence adjustment ####
plot(0,0, ylim = c(0,50), xlim=c(0,50))
for(i in 1:6){
  points(rnasesq_test_confusion.top[[1]][,i],rnasesq_test_confusion.top[[2]][,i])
}

calc_perf_metrics <- function(confusion){
  perfs <- data.frame()
  for(i in rownames(confusion)){
    tp <- confusion[i,i]
    fn <- sum(confusion[i,]) - tp
    fp <- sum(confusion[,i]) - tp
    tn <- sum(confusion) - tp - fp - fn
    
    perfs[i,"ppv"] = tp /(tp + fp)
    perfs[i,"npv"] = tn /(tn + fn)
    perfs[i,"sens"] = tp / (tp + fn)
    perfs[i,"spec"] = tn / (fp + tn)
  }
  perfs
}


prevalence <- data.frame(matrix(nrow = 0, ncol = 6))
prevalence["balanced",] <- rep(0.166,6)
prevalence["EUR_ED",] <- c( 0.2, 0.7,0.03, 0.01, 0.01,0.05)
prevalence["EUR_ED_FWS",] <- c( 0.36, 0.15,0.03, 0.01, 0.01,0.44)
prevalence["EUR_ED_FUO_PF",] <- c( 0.4, 0.16,0.28, 0.05, 0.01,0.1)
prevalence["AFR_OUT",] <- c( 0.15, 0.7,0.01, 0.03, 0.1,0.01)
colnames(prevalence) <- c("bacterial","viral","JIA","TB","malaria","KD")


confusion <- rnasesq_test_confusion.top[[2]]
confusion <- confusion / rowSums(confusion) # normalise for population
perf.uncorrected <- calc_perf_metrics(confusion)

sens <- perf.uncorrected[,"sens"]
spec <- perf.uncorrected[,"spec"]

for(scenario in rownames(prevalence)){
  prev <- unlist(prevalence[scenario,rownames(perf.uncorrected)])
  perf.uncorrected[,paste0("ppv_naive_",scenario)] <- sens * prev / ((sens * prev) + (1-spec) * (1 - prev))
  perf.uncorrected[,paste0("npv_naive_",scenario)] <- spec * (1-prev) / (spec * (1 - prev) + (1 - sens) * prev)
}

prevalence_corrected_performance <- data.frame()
pdf("figures/confusions_circle_size_top_rna_prev.pdf" ,width = 18,height = 15)
for(scenario in rownames(prevalence)){
  print(scenario)
  confusion.adj <- confusion[colnames(prevalence),colnames(prevalence)] * unlist(prevalence[scenario,]) * 100   
  perf <- calc_perf_metrics(confusion.adj)
  perf[colnames(prevalence),"prev"]  <- unlist(prevalence[scenario,])
  
  layout(matrix(c(1,1,1,2) ,
                nrow = 1, ncol = 4, byrow = T))
  
  colord <- c("bacterial","TB","viral","inflammatory","KD","malaria")
  colnames(confusion.adj)[3] <- "inflammatory"
  rownames(confusion.adj)[3] <- "inflammatory"
  rownames(perf)[3] <- "inflammatory"
  perf <- round(perf,digits = 4)
  confusion_heatmap_sizemod(confusion.adj[colord,colord],square_lab = "number",
                            margins = c(21,22,4,2),sf = 7, level = "top",textsize = 3,sf.ax = 3)
  
  title(xlab = "Predicted class",ylab = "True class",line = 17,cex.lab = 3.75)
  title( main =list(prev_explanation[[scenario]] ,cex =3),
         line = 0.5,adj = 0)
  
  pretty_plot_area(ytick = c(0,length(colord),1),cols = c("white","lightgray"),
                   xtick = c(-0.05,0.65,0.01),show_x_tick = F,show_y_tick = F,
                   show_x = F,show_y = F,margins = c(21,0.5,4,2))
  
  text(rep(0 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
       perf[colord,"prev"], cex = 2.5,adj = 0,srt = 45)
  text(rep(0.1 ,length(colord)),seq(length(colord),1,by = -1) - 0.5, ## These are going top to bottom 
       perf[colord,"sens"], cex = 2.5,adj = 0,srt = 45)
  text(rep(0.2 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
       perf[colord,"spec"], cex = 2.5,adj = 0,srt = 45)
  text(rep(0.3 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
       perf[colord,"ppv"], cex = 2.5,adj = 0,srt = 45)
  text(rep(0.4 ,length(colord)),seq(length(colord),1,by = -1) - 0.5,
       perf[colord,"npv"], cex = 2.5,adj = 0,srt = 45)
  par(xpd = T)
  text(c(0,0.1,0.2,0.3,0.4)+0.02,rep(0,5), c("Prevalence","Sensitivity","Specificity","PPV","NPV"),
       cex  = 3*0.7,srt = 90 , adj = 1)
  par(xpd = F)
  
  prevalence_corrected_performance[rownames(perf),paste0(scenario, colnames(perf))] <- perf
  
}
dev.off()

prevalence_corrected_performance <- as.data.frame(t(prevalence_corrected_performance))


prevalence_corrected_performance[c("balancedprev",
                                   "EUR_EDprev",
                                   "EUR_ED_FWSprev",
                                   "EUR_ED_FUO_PFprev",
                                   "AFR_OUTprev"),"Explanation"] <- c("All diseases equal prevalence",
                                                                       "Febrile children undergoing blood tests in European emergency departments",
                                                                       "Hospitalized children with fever without source (<7 days duration)",
                                                                       "Hospitalized children with fever of unknown origin >7 days duration",
                                                                       "African Outpatient children with fever")


prevalence_corrected_performance[c("EUR_EDprev",
                                   "EUR_ED_FWSprev",
                                   "EUR_ED_FUO_PFprev",
                                   "AFR_OUTprev"),"Citation"] <- c("doi.org/10.3389/fped.2021.688272",
                                                                      "doi.org/10.1016/j.jiac.2019.09.015",
                                                                      "doi.org/10.1542/hpeds.2017-0098",
                                                                      "doi.org/10.1056/NEJMoa1214482")


write.csv(prevalence_corrected_performance, file = "figures/prevalence_scenario_performance.csv")

# DE broad category ####



rnaseq_meta.loc <- as.data.frame(rnaseq_a_summary[names(rnaseq_broad_validation_pheno),])
rnaseq_meta.loc[names(rnaseq_broad_validation_pheno),"broadgroup"] <- as.vector(rnaseq_broad_validation_pheno)
rnaseq_meta.loc[,"age"] <- round(as.numeric(rnaseq_meta.loc[,"age"]),digits = 1)
rnaseq_meta.loc[,"sex"] <- rnaseq_meta.loc[,"chromosomal_sex"]
rnaseq_meta.loc[rnaseq_meta.loc[,"sex"] == "XXY","sex"] <- "XY"


library(DESeq2)


gene_counts.filt <- gene_count_matrix[,rownames(rnaseq_meta.loc)]
leg <- rowSums(gene_counts.filt>5)
gene_counts.filt <- gene_counts.filt[leg > 6,] # half smallest group

# disease x vs all 
dds_list <- list()
for(disease in unique(rnaseq_meta.loc[,"broadgroup"])){
  rnaseq_meta.loc[,"bg"] <- rnaseq_meta.loc[,"broadgroup"]
  rnaseq_meta.loc[rnaseq_meta.loc[,"bg"] != disease,"bg"] <- "OD"
  dds_list[[disease]] <- DESeqDataSetFromMatrix(gene_counts.filt , colData = rnaseq_meta.loc,
                                                design = ~ age + sex + bg)  
  dds_list[[disease]]$bg <- relevel(dds_list[[disease]]$bg, ref = "OD")
}


for(disease in names(dds_list)){
  print(disease)
  dds_list[[disease]] <- DESeq(dds_list[[disease]], parallel = T)
}

plot_volcano <- function(res.loc,geneid = F, main = ""){
  pal <- RColorBrewer::brewer.pal(5,"Dark2")
  res.loc <- res.loc[!is.na(res.loc[,"log2FoldChange"]) & !is.na(res.loc[,"padj"]),]
  zerop <- rownames(res.loc)[res.loc[,"padj"] == 0]
  res.loc[0 == res.loc[,"padj"],"padj"] <- min(res.loc[res.loc[,"padj"] > 0,"padj"])
  xvals <- res.loc[,"log2FoldChange"]
  yvals <- -log10(res.loc[,"padj"])
  for( ybr in c(0.1,0.5,1,5,10,15,20,50,100,250,1000)){
    if((max(yvals) - min(yvals) ) / ybr < 15){
      break
    }
  }
  pretty_plot_area(cols = c("white","grey90"),text_col = "grey30",
                   ytick = c(min(yvals),max(yvals)+5,ybr),x_lab = "log2FoldChange",y_lab = "-log10(p(adjusted))",
                   xtick = c(min(xvals)-0.5,max(xvals)+0.5,0.5), ybuffer = 0,main = main)
  sig <- rownames(res.loc)[res.loc[, "padj"]< 0.01 & abs(res.loc[, "log2FoldChange"]) > 1.5]
  sig <- sig[!is.na(sig)]
  title(main = main)
  nsig <- rownames(res.loc)[! rownames(res.loc) %in% sig]
  points(res.loc[nsig,"log2FoldChange"] , -log10(res.loc[nsig,"padj"]) , pch = 16, col =  pal[4])
  
  points(res.loc[sig,"log2FoldChange"] , -log10(res.loc[sig,"padj"]) , pch = 16, col = pal[1])
  if(geneid){
    text(res.loc[sig,"log2FoldChange"] , -log10(res.loc[sig,"padj"]) , labels = sig, pos = 1, col =pal[1])
  }else{
    text(res.loc[sig,"log2FoldChange"] , -log10(res.loc[sig,"padj"]) , labels = genenames[sig,2], pos = 1, col = pal[1])
    
  }
  
}

pdf("figures/volcano_plots.pdf" ,width = 14,height = 14)
plot_volcano( results(dds_list$viral, name=  "bg_viral_vs_OD" ), main = "Viral vs other")
plot_volcano( results(dds_list$bacterial, name=  "bg_bacterial_vs_OD" ), main = "Bacterial vs other")
plot_volcano( results(dds_list$JIA, name=  "bg_JIA_vs_OD" ), main = "JIA vs other")
plot_volcano( results(dds_list$KD, name=  "bg_KD_vs_OD" ), main = "KD vs other")
plot_volcano( results(dds_list$malaria, name=  "bg_malaria_vs_OD" ), main = "Malaria vs other")
plot_volcano( results(dds_list$TB, name=  "bg_TB_vs_OD" ), main = "TB vs other")
dev.off()


lfcThreshold <- 1
baseMeanThreshold <- 10
significant_genes <- list()
res <- results(dds_list$viral, name=  "bg_viral_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_viral_other.csv")
write.csv(res[1:50,], file = "figures/DE_viral_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_viral_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_viral_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_viral_down_1.5_0.01.txt")
significant_genes[["viral"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["viral_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["viral_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold ,])

res <- results(dds_list$bacterial, name=  "bg_bacterial_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_bacterial_other.csv")
write.csv(res[1:50,], file = "figures/DE_bacterial_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_bacterial_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_bacterial_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_bacterial_down_1.5_0.01.txt")
significant_genes[["bacterial"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["bacterial_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["bacterial_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])

res <- results(dds_list$JIA, name=  "bg_JIA_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_JIA_other.csv")
write.csv(res[1:50,], file = "figures/DE_JIA_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_JIA_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_JIA_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_JIA_down_1.5_0.01.txt")
significant_genes[["JIA"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["JIA_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["JIA_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])

res <- results(dds_list$KD, name=  "bg_KD_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_KD_other.csv")
write.csv(res[1:50,], file = "figures/DE_KD_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_KD_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_KD_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_KD_down_1.5_0.01.txt")
significant_genes[["KD"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["KD_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["KD_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])


res <- results(dds_list$malaria, name=  "bg_malaria_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_malaria_other.csv")
write.csv(res[1:50,], file = "figures/DE_malaria_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_malaria_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_malaria_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_malaria_down_1.5_0.01.txt")
significant_genes[["malaria"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["malaria_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["malaria_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])



res <- results(dds_list$TB, name=  "bg_TB_vs_OD" )
res <- res[order(res$padj),]
res[,"Symbol"] <- genenames[rownames(res),2]
write.csv(res, file = "figures/DE_TB_other.csv")
write.csv(res[1:50,], file = "figures/DE_TB_other_top50.csv")
res <- res[!is.na(res[,"padj"]),]
write(paste(rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_TB_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_TB_up_1.5_0.01.txt")
write(paste(rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,]), collapse = " "),
      file = "figures/DE_TB_down_1.5_0.01.txt")
significant_genes[["TB"]] <- rownames(res[abs(res[,"log2FoldChange"]) > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["TB_up"]] <- rownames(res[res[,"log2FoldChange"] > lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])
significant_genes[["TB_down"]] <- rownames(res[res[,"log2FoldChange"] < -lfcThreshold & res[,"padj"] < 0.01 & res[,"baseMean"] > baseMeanThreshold,])

save(significant_genes, file = "figures/significant_genes.RData")


if(F){ # run in version 4 R 
  names(significant_genes)
  library("gprofiler2")
  profiler_res <- list()
  for(i in names(significant_genes)){
    print(i)
    profiler_res[[i]] <-  gost(query = significant_genes[[i]],organism = "hsapiens")
  }
  save(profiler_res, file = "gprofiler2_res.RData", version = "2")
}

load("figures/gprofiler2_res.RData")

pal <- add.alpha(c(RColorBrewer::brewer.pal(8, "Dark2"),"black","grey","lightgrey"),0.7)
names(pal) <- sources

pretty_plot_area <- function(x_lab ="" , y_lab ="" ,main = "",
                             ytick = c(0,2,0.1) ,ytickn = NULL,
                             xtick= c(0,1,0.1) , xtickn = NULL,
                             show_x_tick = T,show_y_tick = T,
                             show_x = T,show_y = T,background = "white",
                             cols = c("lightgray","white"),text_col = "grey20",
                             xbuffer = 0,ybuffer = 0,margins = c(8, 8, 2, 1),
                             add_rect = T,plot_outside_margin = F, axcex = 1,
                             absolute_x = F,absolute_y = F ){
  par(xpd=F)
  xtick2 <- NULL
  ytick2 <- NULL
  if(!is.null(ytickn))  ytick2 <- seq(ytick[1],ytick[2],by = ytickn )
  if(!is.null(xtickn))  xtick2 <- seq(xtick[1],xtick[2],by = xtickn )
  
  # ytick <- seq(ytick[1],ytick[2],by = ytick[3])
  # xtick <- seq(xtick[1],xtick[2],by = xtick[3])
  ytick <- seq(floor(ytick[1] / ytick[3]) * ytick[3],ceiling(ytick[2] / ytick[3]) * ytick[3],by = ytick[3])
  xtick <- seq(floor(xtick[1] / xtick[3]) * xtick[3],ceiling(xtick[2] / xtick[3]) * xtick[3],by = xtick[3])
  
  par(mar = margins, # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i",# Remove plot padding
      bg = background)
  
  # if need to reverse the axes
  if(xtick[1]>xtick[2]){
    xlims <- c(max(xtick)+xbuffer,min(xtick)-xbuffer)
  }else{
    xlims <- c(min(xtick)-xbuffer,max(xtick)+xbuffer)
  }
  if(ytick[1]>ytick[2]){
    ylims <- c(max(ytick)+ybuffer,min(ytick)-ybuffer)
  }else{
    ylims <- c(min(ytick)-ybuffer,max(ytick)+ybuffer)
  }
  
  plot(0,0,cex = 0 , ylim=ylims,
       xlim= xlims,axes =F,xlab = "" ,ylab = "")
  
  title(main =main , xlab = x_lab ,ylab = y_lab,line = 4)
  if(show_x){
    if(absolute_x){
      mtext(side = 1, text = abs(xtick), at = xtick, 
            col = text_col, line = 1, cex = 1.2,las =1) 
    }else{
      mtext(side = 1, text = xtick, at = xtick, 
            col = text_col, line = 1, cex = 1.2,las =1)
    }
  }  
  if(show_y){
    if(absolute_y){
      mtext(side = 2, text = abs(round(ytick, digits = 10)), at = ytick, 
            col = text_col, line = 1, cex = 1.2,las =2)
    }else{
      mtext(side = 2, text = round(ytick, digits = 10), at = ytick, 
            col = text_col, line = 1, cex = 1.2,las =2)  
    }
  } 
  
  if(add_rect){
    rect(xleft = min(xtick)-100,xright = max(xtick)+100 ,
         ybottom =min(ytick)-100 , ytop = max(ytick)+100,col = cols[1],border = F)
  }
  if(show_y_tick)       abline( h = ytick,col =cols[2],lwd =3)
  if(show_x_tick)       abline( v = xtick,col =cols[2],lwd =3)
  if(!is.null(ytick2))  abline( h = ytick2,col =cols[2],lwd =1)
  if(!is.null(xtick2))  abline( v = xtick2,col =cols[2],lwd =1)
  par(xpd=plot_outside_margin)
}






enrichment_plot <- function(src, ylim = c(-20,20), ytextbuffer =0,
                            textjump = 0.4,textsize = 0.6, yjump = 5,rotate = 0){
  pretty_plot_area(ytick = c(ylim,yjump),ybuffer = 5,cols = c("white","lightgray"),margins = rep(6,4),
                   xtick = c(0.8,7,0.2),show_x_tick = F,show_y_tick = T,absolute_y = T,
                   show_x = F,show_y = T,plot_outside_margin = T,y_lab = "-log10(padj)", main = src,axcex = 2)
  for(i in 1:6){
    disease <- c("viral","bacterial","JIA","KD","TB","malaria")[i]
    rect(xleft = i-0.15,xright = i-0.05, ytop = 2, ybottom = -2, col = "white", border = NA)
    text(i-0.1,0, colour_table[disease,"name"], col = colour_table[disease,"col"],srt = 90,cex = 2)
    
    res <- profiler_res[[paste0(disease,"_up")]]$result
    res <- res[res[,"source"] %in% c(src),]
    if(!is.null(res)){
      if(nrow(res) > 0){
        text_lim <- 20
        newypos <- rev(seq(min(-log10(res[,"p_value"])), 2000, by = textjump)[1:nrow(res)])
        newypos <- newypos + ytextbuffer
        if(length(newypos) > 20) newypos <- newypos - newypos[text_lim] + 3
        for( j in 1:text_lim){
          points(c(i,i+0.06),c(-log10(res[j,"p_value"]),newypos[j]), type = "l", col = "lightgrey")
          terms <- res[j,"term_name"]
          terms <- gsub("biological process involved in ","",terms)
          terms <- gsub("regulation","reg",terms)
          terms <- gsub("response","resp",terms)
          if(!is.na(terms)){
            terms <- firstup(terms)
            str_len <- stri_length(terms)
            if(str_len > 50){
              text(i + 0.05, newypos[j]
                   ,terms ,adj = 1,pos = 4,cex = textsize * 0.8,srt = rotate)
              
              
            }else{
              text(i + 0.05, newypos[j]
                   ,terms ,adj = 1,pos = 4,cex = textsize,srt = rotate)
            }
          }
        }
        points(rep(i,nrow(res)),-log10(res[,"p_value"]), col = "red" 
               , pch = 16, cex = 1.5)
      }
    }
    res <- profiler_res[[paste0(disease,"_down")]]$result
    res <- res[res[,"source"] == c(src),]
    if(!is.null(res)){
      if(nrow(res) > 0){
        text_lim <- 20
        newypos <- rev(seq(max(log10(res[,"p_value"])), -2000, by = -1 * textjump)[1:nrow(res)])
        newypos <- newypos - ytextbuffer
        if(length(newypos) > 20) newypos <- newypos - newypos[text_lim] - 3
        for( j in 1:text_lim){
          points(c(i,i+0.06),c(log10(res[j,"p_value"]),newypos[j]), type = "l", col = "lightgrey")
          terms <- res[j,"term_name"]
          terms <- gsub("biological process involved in ","",terms)
          terms <- gsub("regulation","reg",terms)
          terms <- gsub("response","resp",terms)
          if(!is.na(terms)){
            terms <- firstup(terms)
            str_len <- stri_length(terms)
            if(str_len > 50){
              text(i + 0.05, newypos[j]
                   , terms,adj = 1,pos = 4,cex = textsize*0.8,srt = rotate)  
            }else{
              text(i + 0.05, newypos[j]
                   , terms,adj = 1,pos = 4,cex = textsize,srt = rotate)
            }
          }
        }
        points(rep(i,nrow(res)),log10(res[,"p_value"])
               , pch = 16, cex = 1.5, col ="blue")
      }
    }
  }
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

library(stringi)

# save as tables 

for(i in names(profiler_res)){
  res <- profiler_res[[i]]$result
  for(j in 1:nrow(res)) res[j,"Parents"] <- paste0(res[j,"parents"][[1]],collapse = "_")
  res <- res[,c(1:13,15)]
  write.csv(res, file = paste0("figures/broadgroup_DE/gprofiler_",i,".csv"))
}




sources <- c()
for(i in names(profiler_res)){
  sources <- unique(c(sources,profiler_res[[i]]$result[,"source"]))
}


pdf("figures/enrichment.pdf" ,width = 22,height = 10)
for(i in sources){
  ymaxs <- 0
  ymins <- 0
  for(j in names(profiler_res)[grepl("down",names(profiler_res))]){
    res <- profiler_res[[j]]$result
    res <- res[res[,"source"] == i,]
    ymins <- min(c(ymins , log10(res[,"p_value"])))
  }
  for(j in names(profiler_res)[grepl("up",names(profiler_res))]){
    res <- profiler_res[[j]]$result
    res <- res[res[,"source"] == i,]
    ymaxs <- max(c(ymaxs , -log10(res[,"p_value"])))
  }
  
  if(i == "TF") ymaxs <- 10
  
  if(i %in% c("GO:MF","KEGG")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 0.7)
  }else if(i %in% c("GO:CC","REAC")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 2,yjump = 10,textsize = 0.5)
  }else if(i %in% c("GO:BP")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 2.2,yjump = 10)
  }else{
    enrichment_plot(i,ylim = c(ymins,ymaxs))  
  }
  
}
dev.off()



pdf("figures/enrichment.rotated.pdf" ,width = 20,height = 10)
for(i in sources){
  ymaxs <- 0
  ymins <- 0
  for(j in c("viral","bacterial","JIA","KD","TB","malaria")){
    res <- profiler_res[[j]]$result
    res <- res[res[,"source"] == i,]
    ymaxs <- max(c(ymaxs , -log10(res[,"p_value"])))
    ymins <- min(c(ymins , log10(res[,"p_value"])))
  }
  if(i == "TF") ymaxs <- 10
  
  if(i %in% c("GO:MF","KEGG")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 0.9,textsize = 0.8, rotate = 45)
  }else if(i %in% c("REAC")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 2,textsize = 0.8, rotate = 45)
  }else if(i %in% c("GO:CC")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 3,textsize = 0.9,yjump = 10,
                    rotate = 45)
  }else if(i %in% c("GO:BP")){
    enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 3.5,textsize = 0.85,yjump = 10, 
                    rotate = 45)
  }else{
    enrichment_plot(i,ylim = c(ymins,ymaxs),textsize = 0.8)  
  }
  
}
dev.off()




pdf("~/Dropbox/my_files/tmp/enrichment.GO:BP.pdf" ,width = 35,height = 10)
i <- "GO:BP"
ymaxs <- 0
ymins <- 0
for(j in names(profiler_res)[grepl("down",names(profiler_res))]){
  res <- profiler_res[[j]]$result
  res <- res[res[,"source"] == i,]
  ymins <- min(c(ymins , log10(res[,"p_value"])))
}
for(j in names(profiler_res)[grepl("up",names(profiler_res))]){
  res <- profiler_res[[j]]$result
  res <- res[res[,"source"] == i,]
  ymaxs <- max(c(ymaxs , -log10(res[,"p_value"])))
}
enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 3,textsize = 1.4)
dev.off()


pdf("~/Dropbox/my_files/tmp/enrichment.reac.pdf" ,width = 35,height = 10)
i <- "REAC"
ymaxs <- 0
ymins <- 0
for(j in names(profiler_res)[grepl("down",names(profiler_res))]){
  res <- profiler_res[[j]]$result
  res <- res[res[,"source"] == i,]
  ymins <- min(c(ymins , log10(res[,"p_value"])))
}
for(j in names(profiler_res)[grepl("up",names(profiler_res))]){
  res <- profiler_res[[j]]$result
  res <- res[res[,"source"] == i,]
  ymaxs <- max(c(ymaxs , -log10(res[,"p_value"])))
}
enrichment_plot(i,ylim = c(ymins,ymaxs) ,textjump = 2.5,yjump = 10,textsize = 1.2, ytextbuffer = 5)

dev.off()
