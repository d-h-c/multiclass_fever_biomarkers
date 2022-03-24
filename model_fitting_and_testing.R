# initial setup ####
load("data_objects/pre-processed-fold-assigned.RData")
source('scripts/functions.R')
loadlibs()

bacteria <- c("ecoli","gas","gbs","meningococcal","pneumo","staph")
inflammatory <- c("HSP","JIA","SLE")
viral <- c("adeno","enterovirus","flu","HHV6","rhino","RSV")

probemapping <- read.table("dataset/illumina_probes.txt",skip = 16 , header = F)
rownames(probemapping) <- probemapping[,1]
mappingInfo["9qsESLlfOCi_iNeQoI","Symbol"] <- "IFI44L"
mappingInfo["0k_XHUkluCh3ee.9d4","Symbol"] <- "LOC101928092"
mappingInfo["f3y1RvSuIIXQKxXe1I","Symbol"] <- "CAMK1D"
mappingInfo["Knkg_KupfIXrO6KIg4","Symbol"] <- "CDCA2"
mappingInfo["ZkikkSorAkOgCNkqLU","Symbol"] <- "PACRG"

hier.rna <- c(rep("bacterial",length(bacteria)),"TB",rep("viral",length(viral)),"JIA","KD","malaria")
names(hier.rna) <- c(bacteria,"TB",viral,"JIA","KD","malaria")

hier.ma <- c(rep("bacterial",length(bacteria)+1),rep("viral",length(viral)),rep("inflammatory",length(inflammatory)+1),"malaria")
names(hier.ma) <- c(bacteria,"TB",viral,inflammatory,"KD","malaria") 


# run differential expression pre-filter ####

comparisons <- combn(training_groups,2)
training_phenotypes <- training_phenotypes[!is.na(training_phenotypes[,1]),]
training_samples_all <- rownames(training_phenotypes)
de_tables.fulltrain <- mclapply(1:ncol(comparisons), all_de_comparisons ,
                          training_sams.local = training_samples_all,comparisons = comparisons,mc.cores =  10)

input_genes.ma <- pool_de_lists(target_gene_no = 2000,de_tables = de_tables.fulltrain , comparisons.loc = comparisons)
training_genes <- unique(unlist(input_genes.ma))

# Ridge+lasso model fitting ####

## initial fit of lasso ####
ridge_lambdas <- exp(seq(-8,2,by = 0.1))
lasso_lambdas <- rev(exp(seq(-6 , -1.5,by = 0.05)))
lasso_fit_cv.low <- cv.glmnet(t(expression_filt[training_genes , training_samples_all]) ,
                                  training_phenotypes[training_samples_all , "group"],    
                                  weights =  class_weights[training_phenotypes[training_samples_all,"group"] , "cost+imbalance"],
                                  family="multinomial", type.multinomial = "grouped",parallel = T,
                                  alpha=1 ,type.measure = "mse", standardize = TRUE ,lambda = lasso_lambdas,
                                  foldid = training_phenotypes[training_samples_all,"fold"])

folds <- training_phenotypes[training_samples_all,"fold"]

## fit ridges on the lasso ####
fold_predictions <-  mclapply(unique(folds),run_ridges_on_lasso ,mc.cores=5,
                              training_samples = training_samples_all , folds = folds, input_genes = training_genes,
                              lasso_lambdas = lasso_lambdas,ridge_lambdas = ridge_lambdas) 


ridgeylasso_mse.low <- mclapply(lasso_lambdas, get_mse_for_lambda_ridgeylasso , mc.cores = 2 ,
                            fold_predictions = fold_predictions ,training_samples = training_samples_all)

# choose the lasso lambda
mse_plot_dat.low <- select_ridge_lambda(ridgeylasso_mse.low ,lasso_lambdas ,nse = 2)

plot_dat.low <-  plot_ridge_cv(mse_plot_dat.low,xrange = c(-6,-1.5),
                               lasso_fit_cv = lasso_fit_cv.low,mode = "LOW")

## fit second LASSO for relaxed LASSO ######
fold_predictions.relaxed_lasso <-  mclapply(unique(folds),run_ridges_on_lasso ,mc.cores=5,alphaval = 1 ,
                              training_samples = training_samples_all , folds = folds, input_genes = training_genes,
                              lasso_lambdas = lasso_lambdas,ridge_lambdas = ridge_lambdas) 

ridgeylasso_mse.low.relaxed_lasso <- mclapply(lasso_lambdas, get_mse_for_lambda_ridgeylasso , mc.cores = 2 ,
                                fold_predictions = fold_predictions.relaxed_lasso ,training_samples = training_samples_all)

# choose the lasso lambda
mse_plot_dat.low.relaxed_lasso <- select_ridge_lambda(ridgeylasso_mse.low.relaxed_lasso ,lasso_lambdas )

plot_dat.low.relaxed_lasso <-  plot_ridge_cv(mse_plot_dat.low.relaxed_lasso,xrange = c(-6,-1.5),
                                             lasso_fit_cv = lasso_fit_cv.low,mode = "LOW") 


## cv tuning plot over lambda ####
plot_ridge_cv(mse_plot_dat.low,xrange = c(-6,-1.5),lasso_fit_cv = lasso_fit_cv.low,mode = "LOW")
errbar(add = T ,log(lasso_lambdas) , mse_plot_dat.low.relaxed_lasso$mse.r
       ,yplus = mse_plot_dat.low.relaxed_lasso$mse.rp
       ,yminus = mse_plot_dat.low.relaxed_lasso$mse.rm,col = "black")
abline(v = log(mse_plot_dat.low.relaxed_lasso$selected_lambda),col = "black",lty = 2)
abline(v = log(lasso_lambdas)[mse_plot_dat.low.relaxed_lasso$mse.r == min(mse_plot_dat.low.relaxed_lasso$mse.r)],col = "black",lty = 2)

# get selected genes from lasso
sel_lambda <- mse_plot_dat.low$selected_lambda
coefs <- coef(lasso_fit_cv.low,s =sel_lambda)[[1]]
selgenes <- rownames(coefs)[as.matrix(coefs)!=0]
selgenes <- selgenes[selgenes!="(Intercept)"]

# FINAL ridge fit over all training data ####

weights.local <- class_weights[training_phenotypes[training_samples_all,"group"],"cost+imbalance"]
inner_folds <- as.numeric(factor(training_phenotypes[training_samples_all,"fold"]))
ridge_fit_final.low <- cv.glmnet(t(expression_filt[selgenes ,training_samples_all]),
                       y = training_phenotypes_indicator[training_samples_all,rownames(class_weights)]
                       ,foldid = inner_folds,weights = weights.local ,lambda = ridge_lambdas
                       ,type.measure = "mse",alpha = 0,family ="multinomial")

# broad disease classifier fit ##################
disease_set <- c("bacterial","viral","inflammatory","malaria","KD","TB")
groups <- c()
for(sam in training_samples_all){
  groups <- c(groups, disease_set[training_phenotypes_indicator[sam ,disease_set] == 1] )  
}
weights.local <- class_weights_toplevel[groups,"cost+imbalance"]
inner_folds <- as.numeric(factor(training_phenotypes[training_samples_all,"fold"]))
ridge_fit_final.top <- cv.glmnet(t(expression_filt[selgenes  , training_samples_all]),
                                  y = training_phenotypes_indicator[training_samples_all,disease_set]
                                  ,family = "multinomial",type.multinomial = "grouped" ,weights = weights.local,
                                  foldid = inner_folds  ,type = "mse",lambda = lasso_lambdas
                                  ,alpha = 0,standardize = T,parallel = T)

# predict on test set (and training) ####

## test groups
general_test_samples <- rownames(phenotypes)[phenotypes[,"trte"]=="test"]
general_test_samples <- general_test_samples[general_test_samples %in% colnames(expression_filt)]
table(phenotypes[general_test_samples , "group"])

# make test set indicator matrix
test_phenotypes_indicator <- class.ind(factor(phenotypes[general_test_samples,"group"])) 
rownames(test_phenotypes_indicator) <- general_test_samples
specific <- colnames(test_phenotypes_indicator)
test_phenotypes_indicator <- cbind(test_phenotypes_indicator,rep(0,nrow(test_phenotypes_indicator)),
                                   rep(0,nrow(test_phenotypes_indicator)),
                                   rep(0,nrow(test_phenotypes_indicator)))
colnames(test_phenotypes_indicator) <- c(specific , "bacterial","viral","inflammatory")
test_phenotypes_indicator[rowSums(test_phenotypes_indicator[,bacteria]) > 0,"bacterial"] <- 1
test_phenotypes_indicator[rowSums(test_phenotypes_indicator[,inflammatory]) > 0,"inflammatory"] <- 1
test_phenotypes_indicator[rowSums(test_phenotypes_indicator[,viral]) > 0,"viral"] <- 1


prediction_matrices <- list(training = list(),test = list())

# full 
prediction_matrices$test[["full"]] <- list()
prediction_matrices$test$full[["low"]] <- predict(ridge_fit_final.low,
                                newx = t(expression_filt[rownames(coef(ridge_fit_final.low)[[1]])[-1],general_test_samples]) , 
                                s = "lambda.min",type = "response")[,,1]
prediction_matrices$test$full[["top"]] <- predict(ridge_fit_final.top,
                                      newx = t(expression_filt[rownames(coef(ridge_fit_final.top)[[1]])[-1],general_test_samples]) ,
                                    s = "lambda.min",type = "response")[,,1]

# decoupled probabilities (no softmax function to make probabilities sum to 1)
prediction_matrices$test$full[["low.dec"]] <- predict_modified(ridge_fit_final.low,
                                             test_sams.local = general_test_samples ,
                                             lambda = "lambda.min"
                                             ,noprior = F,decouple = T)
prediction_matrices$test$full[["top.dec"]] <- predict_modified(ridge_fit_final.top,
                                                      test_sams.local = general_test_samples , 
                                                      lambda = "lambda.min"
                                                 ,noprior = F,decouple = T)

# predict on training data
prediction_matrices$training[["full"]] <- list()
prediction_matrices$training$full[["low"]] <- predict(ridge_fit_final.low,
                                                  newx = t(expression_filt[rownames(coef(ridge_fit_final.low)[[1]])[-1], rownames(training_phenotypes) ]) , 
                                                  s = "lambda.min",type = "response")[,,1]
prediction_matrices$training$full[["top"]] <- predict(ridge_fit_final.top,
                                                  newx = t(expression_filt[rownames(coef(ridge_fit_final.top)[[1]])[-1], rownames(training_phenotypes) ]) ,
                                                  s = "lambda.min",type = "response")[,,1]
prediction_matrices$training$full[["low.dec"]] <- predict_modified(ridge_fit_final.low,
                                                                      test_sams.local = rownames(training_phenotypes) ,
                                                                      lambda = "lambda.min",noprior = F,decouple = T)
prediction_matrices$training$full[["top.dec"]] <- predict_modified(ridge_fit_final.top,test_sams.local = rownames(training_phenotypes) ,
                                                          lambda = "lambda.min"
                                                      ,noprior = F,decouple = T)
# performance in microarray ####
prediction_matrices <- list_mod(prediction_matrices) # make them all data frames
performance_list <- list_mod.perf(prediction_matrices) # get performance metrics

confusion_microarray_test <- make_confusion(prediction_matrices$test$full$low,
                                            indicator_matrix = test_phenotypes_indicator[rownames(prediction_matrices$test$full$low),
                                                                                         colnames(prediction_matrices$test$full$low)])
confusion_microarray_train <- make_confusion(prediction_matrices$training$full$low,
                                             indicator_matrix = training_phenotypes_indicator[rownames(prediction_matrices$training$full$low),
                                                                                              colnames(prediction_matrices$training$full$low)])


confusion_microarray_test_top <- make_confusion(prediction_matrices$test$full$top,
                                                indicator_matrix = test_phenotypes_indicator[rownames(prediction_matrices$test$full$top),
                                                                                             colnames(prediction_matrices$test$full$top)])
confusion_microarray_train_top <- make_confusion(prediction_matrices$training$full$top,
                                                 indicator_matrix = training_phenotypes_indicator[rownames(prediction_matrices$training$full$top),
                                                                                                  colnames(prediction_matrices$training$full$top)])

ordered_training_groups.top <- c("bacterial","viral","inflammatory","TB","KD","malaria")

# RNA-Seq Data prep ####
get_ensembl <- function(probe,gene = F){
  genesymb <-  as.character(mappingInfo[probe,3])
  genesymb[genesymb=="JARID1A"] <- "KDM5A"
  genesymb[genesymb=="CUGBP1"] <- "CELF1"
  genesymb[genesymb=="LOC641823"] <- "CR16" # WIPF3
  genesymb[genesymb=="LOC648526"] <- "EPPK1"
  genesymb[genesymb=="BHLHB2"] <- "BHLHE40"
  genesymb[genesymb=="TOMM70A"] <- "TOMM70"
  genesymb[genesymb=="LOC644162"] <- "SEPT7"
  genesymb[genesymb=="IL8"] <- "CXCL8"
  genesymb[genesymb=="C14orf131"] <- "ZNF839"
  genesymb[genesymb=="SNX26"] <- "ARHGAP33"
  genesymb[genesymb=="LOC644739"] <- "WASF4"
  genesymb[genesymb=="WASF4"] <- "WASF4P"
  
  genesymb[genesymb=="NGFRAP1"] <- "BEX3"
  genesymb[genesymb=="MGC29506"] <- "MZB1"
  genesymb[genesymb=="SELM"] <- "SELENOM"
  genesymb[genesymb=="SFRS5"] <- "SRSF5"
  genesymb[genesymb=="C15orf28"] <- "RP11-696L21.2"
  genesymb[genesymb=="C9orf21"] <- "AAED1"
  genesymb[genesymb=="LOC644063"] <- "HNRNPKP4"
  genesymb[genesymb=="LOC648000"] <- "RPL7"
  genesymb[genesymb=="SGK"] <- "SGK1"
  genesymb[genesymb=="DEM1"] <- "EXO5"
  genesymb[genesymb=="C1orf183"] <- "FAM212B"
  genesymb[genesymb=="CR16"] <- "WIPF3"
  genesymb[genesymb=="LOC653972"] <- "CBX3"
  genesymb[genesymb=="FAM65A"] <- "RIPOR1"
  
  if(probe == "9qsESLlfOCi_iNeQoI") genesymb <- "IFI44L"
  if(probe == "ZkikkSorAkOgCNkqLU") genesymb <- "PACRG-AS2"
  if(probe == "Knkg_KupfIXrO6KIg4") genesymb <- "CDCA2"
  if(probe == "f3y1RvSuIIXQKxXe1I") genesymb <- "CAMK1D"
  if(probe == "cntX1SX_KLyoJMCqn4") genesymb <- "AAED1"
  if(probe == "HdefXV5T01eVd1OCfg") genesymb <- "OAS1"
  if(probe == "c6dUHXVzXApUek8C54") genesymb <- "CETP"
  if(probe == "NVE8uXlR53H3BSoRo0") genesymb <- "WDFY2"
  if(probe == "W3W6R_5.Xp4id1XvjE") genesymb <- "UBASH3B"  
  if(probe == "Nk5OooPTfS3FJVZ7to") genesymb <- "CTSL"  
  if(probe == "HwlRdVxeiWgRJNNUsI") genesymb <- "TNRC6C"  
  if(probe == "fXUfe5Cfue0wkiNElU") genesymb <- "LINC01089"  
  if(probe == "uR5E7IDqfperSpPiPo") genesymb <- "HIST2H2BF"  
  if(probe == "TQSWUUSSXjlTSSfVVI") genesymb <- "RFX7"  
  if(probe == "TepXsAy7x7DHe.5deU") genesymb <- "FAM160B1"
  if(probe == "uimrUjYJJ0O1I77rJ0") genesymb <- "GEMIN2"  
  if(probe == "NeI0lTt3rZ3fl7o17c") genesymb <- "KIAA1211L"
  if(probe == "T0zNFbuqF9akRJJXiM") genesymb <- "CFAP157"
  if(probe == "06Xqf3il83jl3kvTtU") genesymb <- "SELENON"
  if(probe == "0TnQRAonqQYtQfS45E") genesymb <- "SH3GL1P2"
  if(probe == "TFITfhT0nIFekkEic8") genesymb <- "EBF2"
  if(probe == "l1ILUieKHl7_ScBCR4") genesymb <- "SEC22B"
  if(probe == "WXler3i_eLohjjoiS4") genesymb <- "AC008464.1" ### prev CTC-353G13.1
  if(probe == "xbRJ.Xq__Jdwjrf9d4") genesymb <- "TYMSOS"
  if(probe == "l5f5XV.Ee3lptBCYRI") genesymb <- "LRRC26"
  if(probe == "0Qkr_3iKjrnqoKf6Pk") genesymb <- "AC004080.2" # prev  "HOXA13" but antissense
  if(probe == "okEmkBx6lCOJBVE._o") genesymb <- "EBLN2"
  if(probe == "WUlRXKLv6P406OyTEU") genesymb <- "CDYL2"
  if(probe == "9l5aOODpOCeLpYiPhs") genesymb <- "HNRNPM"
  if(probe == "WQjCT9.Xi9K4eKCJRA") genesymb <- "HNRNPUL1"
  if(probe == "Q6epTsh5ddTi.lJ._I") genesymb <- "PLD1"
  if(probe == "Hay0Gref1XajQd9ShQ") genesymb <- "TESMIN"
  if(probe == "9kAqUtm.c.d2NUp95o") genesymb <- "SPECC1L"
  if(probe == "EaUBNXocgvkqgNASVI") genesymb <- "RAP2C-AS1"
  if(probe == "irAqP6o5unq697dSt0") genesymb <- "" # antisense of FBX07 - no annotated gene
  if(probe == "") genesymb <- ""

  ensembls <- genenames[genenames[,2] %in% genesymb,1]
  ensembls <- ensembls[ensembls%in% rownames(batcha_expr_sets$dn.ln)]
  if(probe =="r731U7Xo5R6C5Eefs4" ){
    genesymb <- "ZCCHC2"
    ensembls <- "ENSG00000141664" 
    # maps to antisense gene  AC064801.1 on the antisense strand -- on ZZCHC2 it is between exons but may detect pre-mrna
  } 
  if(probe =="rWiFXoPwNPJFmV9Vgg" ){
    genesymb <- "RNPC3"
    ensembls <- "ENSG00000185946"
  } 
  if(probe == "x64jX3XE4pS3hC3qIo" ){
    genesymb <- "" # SCMH1"
  } 
  if(probe == "0k_XHUkluCh3ee.9d4" ){ # will remove AS NOT IN ANNOTATION 
    genesymb <- ""
    
  } 
  if( probe ==  "95Eq6PMq6OIpOgVges" ){
    genesymb <-  "95Eq6PMq6OIpOgVges"
    
  } 
  if(probe ==  "0W6nuULuVKju.kCoTc"){
    ensembls <- "ENSG00000275691"  
    genesymb <- "MT1IP"
  } 
  
  if(probe ==  "0yleB.dSUJ3ipUXXl0"){ # removing as multiple potential genes
    ensembls <- c() # c("ENSG00000230724"  ,"ENSG00000236438")
    genesymb <- c() #c("LINC01001","FAM157A")
  } 
  if(probe ==  "ZV51LRFaCeHrUZEngk"){
    ensembls <- c() #c("ENSG00000164180"  ,"ENSG00000089472")
    genesymb <- c() #c("TMEM161B","HEPH")
  } 
  if(gene){
    return(genesymb)
  }else{
    return(ensembls)  
  }
}

#  http://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz
genome_annotation <- read.csv("dataset/Homo_sapiens.GRCh38.89.gtf",
                              sep = "\t", skip = 5,header = F,stringsAsFactors = F) 

genenames <- cbind(gsub(";.*","",gsub(".*gene_id ","",genome_annotation[,"V9"])),
                     gsub(";.*","",gsub(".*gene_name ","",genome_annotation[,"V9"])))
genenames <- unique(genenames)
rownames(genenames)  <- genenames[,1] 
selgenes.ensembl <- c() 
for(i in selgenes){
  ens <- get_ensembl(i)
  if(length(ens)>=1){
    selgenes.ensembl  <- c(selgenes.ensembl , ens)  
  }else{
    print( probemapping[i,])
  }
}

rnaseq_meta <- read.csv("dataset/RNA-Seq_metadata.csv",stringsAsFactors = F,row.names = 1,check.names = F)
gene_count_matrix <- read.csv("dataset/RNA-Seq_genecounts.csv",row.names = 1, check.names = F)


expressed_genes <- rownames(gene_count_matrix[rowSums(gene_count_matrix > 5) > 10, ])
length(expressed_genes)
filtered_genes <- read.csv("dataset/filtered_genes.txt", header = F,stringsAsFactors = F)[,1] # removed globin genes and sex chromosomes
length(filtered_genes)
expressed_genes <- intersect(expressed_genes,filtered_genes)
length(expressed_genes)

countmat.filt <- gene_count_matrix[,rownames(rnaseq_meta)]
sf <- estimateSizeFactorsForMatrix(countmat.filt)

batcha_expr.dn <- t(t(countmat.filt)/sf[colnames(countmat.filt)])


batcha_expr_sets  <- list(dn.log = log(batcha_expr.dn+1),
                          dn.voom = voom(batcha_expr.dn,normalize.method = "none")$E,
                          sf = sf)

rnaseq_validation_pheno <- rnaseq_meta[,"disease_specific"]
names(rnaseq_validation_pheno) <- rownames(rnaseq_meta)

# training test split for rnaseq ####
all_rnaseq_samples <- rownames(rnaseq_meta)
folds_rnaseq <- list()
train_rnaseq <- c()
test_rnaseq <- c()
groups <- unique(rnaseq_meta[,"disease_specific"])
for(class in groups){
  set.seed(4321)
  grp_sams <- names(rnaseq_validation_pheno)[rnaseq_validation_pheno == class] 
  grp_sams <- sample(grp_sams,length(grp_sams))
  train.loc <- grp_sams[1:(table(rnaseq_validation_pheno[grp_sams]) /2)[class]]
  train_rnaseq <- c(train_rnaseq ,train.loc)
  test_rnaseq <- c(test_rnaseq ,grp_sams[!grp_sams %in%train_rnaseq ])
  folds.loc <-  rep(1:10,100)[1:length(train.loc)]
  names(folds.loc) <- train.loc
  folds_rnaseq <- c(folds_rnaseq , folds.loc)
}


# pca #####
pca.uncorrected.rnaseq <- prcomp(log(t(countmat.filt[expressed_genes,]+1)),center=T,scale.=T)     
pca.corrected.rnaseq <- prcomp(t(batcha_expr_sets$dn.log[expressed_genes,colnames(countmat.filt)]),
                               center=TRUE,scale.=TRUE)    
pca.sigspace.rnaseq <- prcomp(t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered,colnames(countmat.filt)]),
                               center=TRUE,scale.=TRUE)    

# genes microarray to rnaseq filtering ####

selgenes.ensembl <- unique(selgenes.ensembl)
present <- rowSums(batcha_expr_sets$dn[selgenes.ensembl,c(train_rnaseq,test_rnaseq)] > 5 ) 
filtered_genes <- names(present[present<10])
selgenes.ensembl.filtered <- names(present[present>10])


## probe->gene conversion table ####

probe_conv <- data.frame(row.names = selgenes)
probe_conv[,1] <- NA
probe_conv[,2] <- F
probe_conv[,3] <- NA
for( i in selgenes){
  ens <- get_ensembl(i,gene = F)
  if(i %in%  c("0k_XHUkluCh3ee.9d4","0yleB.dSUJ3ipUXXl0",
               "x64jX3XE4pS3hC3qIo","ZV51LRFaCeHrUZEngk",
               "irAqP6o5unq697dSt0")){
    probe_conv[i,1] <- ""
    ens <- ""
  }else{
    probe_conv[i,1] <- paste(ens,collapse = "+")
    probe_conv[i,3] <- paste(genenames[ens,2]  ,collapse = "+")
  }
  if(ens %in% selgenes.ensembl.filtered){
    probe_conv[i,2] <- T
  }else{
    probe_conv[i,2] <- F
  }
}
probe_conv["irAqP6o5unq697dSt0" , "ensembl id"] <- "No annotated gene"
probe_conv["0k_XHUkluCh3ee.9d4" , "ensembl id"] <- "No annotated gene"
probe_conv["0yleB.dSUJ3ipUXXl0" , "ensembl id"] <- "Multiple potential genes (ENSG00000230724/ENSG00000236438)"
probe_conv["ZV51LRFaCeHrUZEngk" , "ensembl id"] <- "Multiple potential genes (ENSG00000164180/ENSG00000089472)"
probe_conv["x64jX3XE4pS3hC3qIo" , "ensembl id"] <- "No annotated gene"
colnames(probe_conv) <- c("ensembl id", "use for rnaseq","Symbol")
write.csv(probe_conv , file = "figures/probe_conversion.csv")

# RNA-Seq translation ####
transferprobes <- rownames(probe_conv[probe_conv[,2],])

## refit coefficients on the whole microarray dataset ####
phenotypes_indicator <- as.data.frame(rbind(test_phenotypes_indicator, training_phenotypes_indicator))
sams.tmp <- rownames(phenotypes_indicator)
set.seed(123)
phenotypes_indicator[sample(sams.tmp), "fold"] <- rep(1:10,200)[1:1212]
inner_folds <- phenotypes_indicator[,"fold"]


## refit on only transferred 145 probes
weights.local <- class_weights[phenotypes[sams.tmp,"group"],"cost+imbalance"]
trainmat.loc <- t(expression_filt[transferprobes ,sams.tmp])
ridge_refit.low <- cv.glmnet(trainmat.loc,
                             y = phenotypes[sams.tmp,"group"]
                             ,foldid = inner_folds,weights = weights.local, lambda = exp(seq(-10,5,by = 0.2))
                             ,type.measure = "mse",alpha = 0,family ="multinomial")

weights.local <- class_weights_toplevel[hier.ma[phenotypes[sams.tmp,"group"]],"cost+imbalance"]
ridge_refit.top <- cv.glmnet(trainmat.loc,
                             y = hier.ma[phenotypes[sams.tmp,"group"]]
                             ,foldid = inner_folds,weights = weights.local , lambda = exp(seq(-10,5,by = 0.2))
                             ,type.measure = "mse",alpha = 0,family ="multinomial")


# on all 161 discovered probes
weights.local <- class_weights[phenotypes[sams.tmp,"group"],"cost+imbalance"]
trainmat.loc <- t(expression_filt[selgenes,sams.tmp])
ridge_refit.allprobe.low <- cv.glmnet(trainmat.loc,
                                      y = phenotypes[sams.tmp,"group"]
                                      ,foldid = inner_folds,weights = weights.local, lambda = exp(seq(-10,5,by = 0.2))
                                      ,type.measure = "mse",alpha = 0,family ="multinomial")

weights.local <- class_weights_toplevel[hier.ma[phenotypes[sams.tmp,"group"]],"cost+imbalance"]
ridge_refit.allprobe.top <- cv.glmnet(trainmat.loc,
                                      y = hier.ma[phenotypes[sams.tmp,"group"]]
                                      ,foldid = inner_folds,weights = weights.local , lambda = exp(seq(-10,5,by = 0.2))
                                      ,type.measure = "mse",alpha = 0,family ="multinomial")


# No intercept
weights.local <- class_weights[phenotypes[sams.tmp,"group"],"cost+imbalance"]
ridge_refit.noi.low <- cv.glmnet(trainmat.loc,
                                 y = phenotypes[sams.tmp,"group"]
                                 ,foldid = inner_folds,weights = weights.local, lambda = exp(seq(-10,5,by = 0.2))
                                 ,type.measure = "mse",alpha = 0,family ="multinomial", intercept = F)

weights.local <- class_weights_toplevel[hier.ma[phenotypes[sams.tmp,"group"]],"cost+imbalance"]
ridge_refit.noi.top <- cv.glmnet(trainmat.loc,
                                 y = hier.ma[phenotypes[sams.tmp,"group"]]
                                 ,foldid = inner_folds,weights = weights.local , lambda = exp(seq(-10,5,by = 0.2))
                                 ,type.measure = "mse",alpha = 0,family ="multinomial", intercept = F)

## Predict on RNA-Seq (voom) with these coeficients ####


expr.allprobe.voom <- batcha_expr_sets$dn.voom[probe_conv[transferprobes,1],names(rnaseq_broad_validation_pheno)]
rownames(expr.allprobe.voom) <- transferprobes
expr.allprobe.voom <- as.data.frame(expr.allprobe.voom)
expr.allprobe.voom[selgenes[selgenes %ni% transferprobes],] <- 0 # replace filtered probes with 0

rnaseq_pred.allprobe.low <-  glmnet::predict.cv.glmnet( ridge_refit.allprobe.low ,
                                                        newx = t(expr.allprobe.voom),
                                                        s = "lambda.1se",type ="response")[,,1]

rnaseq_pred.allprobe.top <-  glmnet::predict.cv.glmnet( ridge_refit.allprobe.top ,
                                                        newx = t(expr.allprobe.voom),
                                                        s = "lambda.1se",type ="response")[,,1]


# refit coefficients in rnaseq 50/50 ####

rnaseq_validation_pheno.ind <- class.ind(rnaseq_validation_pheno)
rnaseq_validation_pheno.ind.loc <- as.data.frame(rnaseq_validation_pheno.ind)
rnaseq_validation_pheno.ind.loc[rowSums(rnaseq_validation_pheno.ind[,bacteria[-3]]) == 1,"bacterial"] <- 1 
rnaseq_validation_pheno.ind.loc[rowSums(rnaseq_validation_pheno.ind[,viral[-c(2,4)]]) == 1,"viral"] <- 1 
rnaseq_validation_pheno.ind.loc[is.na(rnaseq_validation_pheno.ind.loc)] <- 0

refit_class_weights <- data.frame() #class_weights
refit_class_weights[colnames(rnaseq_validation_pheno.ind.loc),"freq"] <- colSums(rnaseq_validation_pheno.ind.loc)
refit_class_weights[rownames(class_weights),"cost" ] <- class_weights[,"cost"]
refit_class_weights[c("bacterial","viral"),"cost" ] <- c(10,2)
refit_class_weights[,"cost+imbalance"] <- refit_class_weights[,"cost"] / refit_class_weights[,"freq"]
refit_class_weights[,"imbalance"] <- 1  / refit_class_weights[,"freq"]


rnaseq_groups_low <- colnames(rnaseq_validation_pheno.ind)
refit_ridge.rna.mn__WCI.low <- cv.glmnet( x =  t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered , train_rnaseq  ]),standardize = T,#intercept = F,
                                      y = as.matrix(rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_low]),type.measure = "mse",
                                      alpha = 0,family ="multinomial",foldid = as.numeric(folds_rnaseq[train_rnaseq]), lambda = exp(seq(2,-4,by = -0.1)),
                                      weights = refit_class_weights[rnaseq_validation_pheno[train_rnaseq ],"cost+imbalance"] )


rnaseq_groups_top <- c("bacterial","viral","JIA","TB","KD","malaria")
refit_ridge.rna.mn__WCI.top <- cv.glmnet( x = t( batcha_expr_sets$dn.log[selgenes.ensembl.filtered , train_rnaseq  ]),standardize = T,# intercept = F,
                                          y = as.matrix(rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_top]),type.measure = "mse",
                                          alpha = 0,family ="multinomial", foldid = as.numeric(folds_rnaseq[train_rnaseq]), 
                                          lambda = exp(seq(-6,2,by = 0.2)),
                                          weights = refit_class_weights[hier.rna[rnaseq_validation_pheno[train_rnaseq ]],"cost+imbalance"] )




pred.resp.refit_mn_WCI.low.dec  <- pred_refit_rnaseq(samples = c(test_rnaseq),
                                                fit = refit_ridge.rna.mn__WCI.low,
                                                 expr.loc =  batcha_expr_sets$dn.log,
                                                 lambda = "lambda.min",type = "decoupled.logit")
pred.resp.refit_mn_WCI.top.dec  <- pred_refit_rnaseq(samples = c(test_rnaseq),
                                                  fit = refit_ridge.rna.mn__WCI.top,
                                                  expr.loc =  batcha_expr_sets$dn.log,
                                                  lambda = "lambda.min",type = "decoupled.logit")


pred.resp.refit_mn_WCI.low <- glmnet::predict.cv.glmnet( refit_ridge.rna.mn__WCI.low ,
                                  newx = t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered, c(test_rnaseq)]),
                                  s = "lambda.min",type ="response")[,,1]
pred.resp.refit_mn_WCI.top <- glmnet::predict.cv.glmnet( refit_ridge.rna.mn__WCI.top ,
                                  newx = t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered, c(test_rnaseq)]),
                                  s = "lambda.min",type ="response")[,,1]

pred.resp.refit_mn_WCI.low.train <- glmnet::predict.cv.glmnet( refit_ridge.rna.mn__WCI.low ,
                                                         newx = t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered, c(train_rnaseq)]),
                                                         s = "lambda.min",type ="response")[,,1]
pred.resp.refit_mn_WCI.top.train <- glmnet::predict.cv.glmnet( refit_ridge.rna.mn__WCI.top ,
                                                         newx = t(batcha_expr_sets$dn.log[selgenes.ensembl.filtered, c(train_rnaseq)]),
                                                         s = "lambda.min",type ="response")[,,1]

pred.resp.refit_mn_WCI.low.dec.train  <- pred_refit_rnaseq(samples = c(train_rnaseq),
                                                     fit = refit_ridge.rna.mn__WCI.low,
                                                     expr.loc =  batcha_expr_sets$dn.log,
                                                     lambda = "lambda.min",type = "decoupled.logit")
pred.resp.refit_mn_WCI.top.dec.train  <- pred_refit_rnaseq(samples = c(train_rnaseq),
                                                     fit = refit_ridge.rna.mn__WCI.top,
                                                     expr.loc =  batcha_expr_sets$dn.log,
                                                     lambda = "lambda.min",type = "decoupled.logit")


make_01 <- function(mat){
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat
}

colord <- c(bacteria,viral,inflammatory,"KD","malaria","TB","bacterial" , "viral")
colord <- colord[colord%in%colnames(rnaseq_validation_pheno.ind.loc)]



rnaseq_validation_performance_mn_low <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.low[test_rnaseq,rnaseq_groups_low] ,
                                                        indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,rnaseq_groups_low])
rnaseq_validation_performance_mn_top <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.top[test_rnaseq,rnaseq_groups_top] ,
                                                        indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,rnaseq_groups_top])
rnaseq_validation_performance_mn_low.dec <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.low.dec[test_rnaseq,rnaseq_groups_low] ,
                                                               indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,rnaseq_groups_low])
rnaseq_validation_performance_mn_top.dec <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.top.dec[test_rnaseq,rnaseq_groups_top] ,
                                                               indicator = rnaseq_validation_pheno.ind.loc[test_rnaseq,rnaseq_groups_top])
rnaseq_validation_performance_mn_low.train <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.low.train[train_rnaseq,rnaseq_groups_low] ,
                                                               indicator = rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_low])
rnaseq_validation_performance_mn_top.train <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.top.train[train_rnaseq,rnaseq_groups_top] ,
                                                               indicator = rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_top])
rnaseq_validation_performance_mn_low.dec.train <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.low.dec.train[train_rnaseq,rnaseq_groups_low] ,
                                                                     indicator = rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_low])
rnaseq_validation_performance_mn_top.dec.train <- get_performance.rnaseq(pred_resp = pred.resp.refit_mn_WCI.top.dec.train[train_rnaseq,rnaseq_groups_top] ,
                                                                     indicator = rnaseq_validation_pheno.ind.loc[train_rnaseq,rnaseq_groups_top])

# rnaseq confusion ####

rnasesq_test_confusion.low <- make_confusion(pred.resp.refit_mn_WCI.low[test_rnaseq,],
                                         rnaseq_validation_pheno.ind.loc[test_rnaseq,colnames(pred.resp.refit_mn_WCI.low)])

rnasesq_test_confusion.top <- make_confusion(pred.resp.refit_mn_WCI.top[test_rnaseq,],
                                         rnaseq_validation_pheno.ind.loc[test_rnaseq,colnames(pred.resp.refit_mn_WCI.top)])

# in training set

rnasesq_train_confusion.low <- make_confusion(pred.resp.refit_mn_WCI.low.train[train_rnaseq,],
                                             rnaseq_validation_pheno.ind.loc[train_rnaseq,colnames(pred.resp.refit_mn_WCI.low)])
rnasesq_train_confusion.top <- make_confusion(pred.resp.refit_mn_WCI.top.train[train_rnaseq,],
                                             rnaseq_validation_pheno.ind.loc[train_rnaseq,colnames(pred.resp.refit_mn_WCI.top)])


# level disagreement ####

level_disagreement <- function(pred.resp.low,pred.resp.top,indicator,test_samples,mod = "rnaseq"){
  disagreement_matrix <- data.frame(matrix(nrow = ncol(pred.resp.low),ncol = ncol(pred.resp.top)))
  rownames(disagreement_matrix) <- colnames(pred.resp.low)
  colnames(disagreement_matrix) <- colnames(pred.resp.top)
  disagreement_matrix[] <- 0
  
  confusion_matrix.low <- data.frame(matrix(nrow = ncol(pred.resp.low),ncol = ncol(pred.resp.low)))
  rownames(confusion_matrix.low) <- colnames(pred.resp.low)
  colnames(confusion_matrix.low) <- colnames(pred.resp.low)
  confusion_matrix.low[] <- 0
  

  if(mod == "rnaseq"){
    confusion_matrix.top <- data.frame(matrix(nrow = ncol(pred.resp.top)+1,ncol = ncol(pred.resp.top)))
    rownames(confusion_matrix.top) <- c(colnames(pred.resp.top),"inflammatory")  
    hier.loc <- hier.rna
    hier.loc[hier.loc=="bacteria"] <- "bacterial"
    hier.loc["KD"] <- "KD"
    hier.loc["TB"] <- "TB"
    hier.loc["JIA"] <- "JIA"
  }else{
    hier.loc <- hier.ma
    hier.loc[hier.loc=="bacteria"] <- "bacterial"
    hier.loc["KD"] <- "KD"
    hier.loc["TB"] <- "TB"
    confusion_matrix.top <- data.frame(matrix(nrow = ncol(pred.resp.top),ncol = ncol(pred.resp.top)))
    rownames(confusion_matrix.top) <- c(colnames(pred.resp.top))
  }
  colnames(confusion_matrix.top) <- colnames(pred.resp.top)
  confusion_matrix.top[] <- 0
  disagreement_matrices <- list(confident.low = confusion_matrix.low,
                              confident.top = confusion_matrix.top,
                              ambiguous.low = confusion_matrix.low,
                              ambiguous.top = confusion_matrix.top)
  
  disagreement_sample_list <- list(agree = c(),disagree = c())
  disagreement_matrix.oneright <-  disagreement_matrix
  disagreement_matrix.bothwrong <- disagreement_matrix
  for(sam in test_samples){
    tc <- colnames(pred.resp.low)[indicator[sam,colnames(pred.resp.low)]==1]
    pvl <-pred.resp.low[sam,colnames(pred.resp.low)]
    pvt <-pred.resp.top[sam,colnames(pred.resp.top)]
    pl <- names(pvl)[pvl==max(pvl)]
    pt <- names(pvt)[pvt==max(pvt)]
    
    if(hier.loc[pl] == pt){ # agree
      disagreement_matrices$confident.low[tc,pl] <- disagreement_matrices$confident.low[tc,pl] +1
      disagreement_matrices$confident.top[as.character(hier.loc[tc]),pt] <- disagreement_matrices$confident.top[as.character(hier.loc[tc]),pt] +1
      disagreement_sample_list$agree <- c(disagreement_sample_list$agree ,sam)
    }else{ # disagree
      disagreement_matrices$ambiguous.low[tc,pl] <- disagreement_matrices$ambiguous.low[tc,pl] +1
      disagreement_matrices$ambiguous.top[as.character(hier.loc[tc]),pt] <- disagreement_matrices$ambiguous.top[as.character(hier.loc[tc]),pt] +1
      disagreement_matrix[pl,pt] <- disagreement_matrix[pl,pt] + 1
      disagreement_sample_list$disagree <- c(disagreement_sample_list$disagree ,sam)
      if(tc == pl | hier.loc[tc] == pt){ # of one is right
        disagreement_matrix.oneright[pl,pt] <- disagreement_matrix.oneright[pl,pt] + 1
      }else{
        disagreement_matrix.bothwrong[pl,pt] <- disagreement_matrix.bothwrong[pl,pt] + 1
      }
    }
  }
  disagreement_matrices[["leveldisagreement"]] <- disagreement_matrix
  disagreement_matrices[["leveldisagreement.oneright"]] <- disagreement_matrix.oneright
  disagreement_matrices[["leveldisagreement.bothwrong"]] <- disagreement_matrix.bothwrong
  disagreement_matrices[["samples"]] <- disagreement_sample_list
  disagreement_matrices
}





# disagreement if true and predicted class are different  --(only look at the max(p) prediction)

rnaseq_disagreement <- level_disagreement(pred.resp.low = pred.resp.refit_mn_WCI.low,
                                           pred.resp.top = pred.resp.refit_mn_WCI.top,
                                           indicator = rnaseq_validation_pheno.ind.loc,
                                           test_samples = test_rnaseq,mod ="rnaseq")



microarray_disagreement <- level_disagreement(pred.resp.low = prediction_matrices$test$full$low,
                                          pred.resp.top = prediction_matrices$test$full$top,
                                          indicator = test_phenotypes_indicator,
                                          test_samples = general_test_samples ,mod ="micro")


sum(microarray_disagreement$ambiguous.low) / sum(microarray_disagreement$confident.low)
sum(rnaseq_disagreement$ambiguous.low) / sum(rnaseq_disagreement$confident.low)


# comparison with existing signatures ########################################################

# sweeney TB
# Lancet Respir Med . 2016 March ; 4(3): 213–224. doi:10.1016/S2213-2600(16)00048-5.
# (GBP5 + DUSP3) / 2 – KLF2
gbp5 <- genenames[genenames[,2] == "GBP5",1]
dusp3 <- genenames[genenames[,2] == "DUSP3",1]
klf2 <- genenames[genenames[,2] == "KLF2",1]

sweeney_genes.ens <- c(gbp5=gbp5,dusp3=dusp3,klf2=klf2)
sweeney_coef <- c(0,0.5,0.5 , -1)
names(sweeney_coef) <- c("(Intercept)",sweeney_genes.ens)

# Twogene
ifi44l <- genenames[genenames[,2] == "IFI44L",1]
fam89a <- genenames[genenames[,2] == "FAM89A",1]

twogene_genes.ens <- c(ifi44l,fam89a)
names(twogene_genes.ens) <- c("ifi44l","fam89a")
twogene_coef <- c(0,1,-1)
names(twogene_coef) <- c("(Intercept)" , twogene_genes.ens)


# KD13
library(illuminaHumanv4.db)
adrToIllumina = toTable(illuminaHumanv4ARRAYADDRESS)
kawasaki_genes <- c("CACNA1E","DDIAS","KLHL2",
                    "PYROXD2","SMOX","ZNF185",
                    "LINC02035","CLIC3","S100P",
                    "IFI27", "TIGIT","CD163","RTN1")    # tigit == HS.553068
kawasaki_genes.ens <- c()
for(G in kawasaki_genes){
  if(G =="IFI27"){
    #  "ENSG00000275214"  is on an extra chromosome
    kawasaki_genes.ens <- c(kawasaki_genes.ens ,"ENSG00000165949" )    
  }else if(G =="HS.553068"){
    kawasaki_genes.ens <- c(kawasaki_genes.ens , "")  
  }else{
    kawasaki_genes.ens <- c(kawasaki_genes.ens , rownames(genenames)[genenames[,2]==G])  
  }
  
  print(c(G, rownames(genenames)[genenames[,2]==G]))
}
names(kawasaki_genes.ens) <- kawasaki_genes

kawasaki_coef <- c(0, # intercept
                   0.955,
                   0.844,
                   0.789,
                   0.727,
                   0.675,
                   0.646,
                   0.561,
                   0.464,
                   -0.405,
                   -0.426,
                   -0.599,
                   -0.638,
                   -0.690) 

names(kawasaki_coef) <- c("(Intercept)",kawasaki_genes)


# probe : 1470450  HS.553068 is TIGIT -- (stranding)

# NEJM TB
NEJM_TB_SIG <- data.frame(stringsAsFactors = F,matrix(ncol = 4,byrow = T,c(
  "4180768", "ALAS2", "ILMN_13644", "UP",
  "1070477", "ALDH1A1", "ILMN_177898", "UP",
  "5910019", "C1QB", "ILMN_36274", "UP",
  "4290026", "C20ORF103", "ILMN_165304", "DOWN",
  "2600634", "C3HC4", "ILMN_6980", "DOWN",
  "1580048", "CAST", "ILMN_163108", "UP",
  "3390564", "CCDC52", "ILMN_23129", "UP",
  "3940754", "CD226", "ILMN_3877", "UP",
  "1780440", "CD79A", "ILMN_37614", "UP",
  "5890653", "CDKN1C", "ILMN_20689", "DOWN",
  "5340767", "CEACAM1", "ILMN_21651", "DOWN",
  "130086", "CYB561", "ILMN_8373", "UP",
  "840446", "CYB561", "ILMN_20474", "UP",
  "4540239", "DEFA1", "ILMN_29692", "UP",
  "1050068", "F2RL1", "ILMN_176188", "UP",
  "6510707", "FER1L3", "ILMN_18562", "UP",
  "6840767", "FRMD3", "ILMN_11826", "DOWN",
  "2350189", "GBP3", "ILMN_3653", "UP",
  "1510364", "GBP5", "ILMN_24462", "UP",
  "3780047", "GBP6", "ILMN_1956", "UP",
  "6220739", "GRAMD1B", "ILMN_308544", "DOWN",
  "5260484", "HLA-DRB1", "ILMN_20550", "UP",
  "6370315", "HLA-DRB5", "ILMN_3178", "UP",
  "620544", "HLA-DRB6", "ILMN_5312", "UP",
  "630619", "HPSE", "ILMN_165418", "DOWN",
  "5340762", "HS.106234", "ILMN_74965", "UP",
  "7320678", "HS.171481", "ILMN_80341", "UP",
  "4880370", "JUP", "ILMN_3789", "DOWN",
  "1050215", "KCNJ15", "ILMN_164363", "DOWN",
  "2570438", "KIFC3", "ILMN_4695", "UP",
  "7560114", "KLHDC8B", "ILMN_6513", "UP",
  "5310445", "KREMEN1", "ILMN_41914", "DOWN",
  "4570164", "LOC389386", "ILMN_165610", "UP",
  "4780044", "LOC389386", "ILMN_352098", "UP",
  "2350121", "LOC642678", "ILMN_38908", "UP",
  "6480364", "LOC647460", "ILMN_38026", "DOWN",
  "6900291", "LOC649210", "ILMN_33006", "DOWN",
  "830639", "LOC653778", "ILMN_32201", "DOWN",
  "2260349", "MIR1974", "ILMN_388657", "DOWN",
  "830750", "NCF1B", "ILMN_168368", "UP",
  "6760593", "OSBPL10", "ILMN_11112", "UP",
  "3170246", "PDCD1LG2", "ILMN_3561", "UP",
  "2000292", "SCGB3A1", "ILMN_23096", "DOWN",
  "160368", "SEMA6B", "ILMN_21277", "DOWN",
  "1400593", "SIGLEC14", "ILMN_309673", "UP",
  "460463", "SMARCD3", "ILMN_19301", "UP",
  "540520", "SNORD8", "ILMN_366693", "UP",
  "1240554", "TNFRSF17", "ILMN_17574", "UP",
  "4760747", "TPST1", "ILMN_174128", "UP",
  "2630195", "VAMP5", "ILMN_20179", "DOWN",
  "3940088", "ZBED2", "ILMN_4927", "DOWN")))


for(i in 1:nrow(NEJM_TB_SIG)){
  print("=========================================")
  probe <- IlluminaID2nuID(IlluminaID = adrToIllumina[adrToIllumina[,2] == NEJM_TB_SIG[i,1],1])
  if(NEJM_TB_SIG[i,2] %in% genenames[,2]){
    ens <- rownames(genenames)[genenames[,2]==NEJM_TB_SIG[i,2]] 
    if(length(ens)==1){
      print(c(probe,ens))
      NEJM_TB_SIG[i,5] <- ens
    }else{
      print(c(NEJM_TB_SIG[i,2] , ens))
    }
  }else{
    print("---------------------------------------------------")
    print(NEJM_TB_SIG[i,2])
  }
}

for(i in 1:nrow(NEJM_TB_SIG)){
  print("=========================================")
  probe <- IlluminaID2nuID(IlluminaID = adrToIllumina[adrToIllumina[,2] == NEJM_TB_SIG[i,1],1])
  NEJM_TB_SIG[i,6] <- probe[,"nuID"]
  NEJM_TB_SIG[i,7] <- probemapping[probe[,"nuID"],"V3"]
}  


NEJM_TB_SIG <- NEJM_TB_SIG[!duplicated(NEJM_TB_SIG[,2]),]
rownames(NEJM_TB_SIG) <- NEJM_TB_SIG[,2]
NEJM_TB_SIG[is.na(NEJM_TB_SIG[,5]),]
NEJM_TB_SIG["C20ORF103",5] <- "ENSG00000125869" # LAMP5
NEJM_TB_SIG["C3HC4",5] <- "ENSG00000165406" # MARCHF8 
NEJM_TB_SIG["CCDC52",5] <- "ENSG00000163611" # spice1
NEJM_TB_SIG["CDKN1C",5] <- "ENSG00000129757" # two genes , this one in primary assembly 
NEJM_TB_SIG["FER1L3",5] <- "ENSG00000138119" #  myof
NEJM_TB_SIG["HS.106234",5] <- "ENSG00000226496" # LINC00323 --filtered out
NEJM_TB_SIG["LOC647460",5] <- "ENSG00000239819" # IGKV1D-8
NEJM_TB_SIG["LOC389386",5] <- "ENSG00000124762" # not PANDAR wrong strand  -- ENSG00000213500 and ENSG00000124762 are on right strand --taking the gene (is up in TB) rather than the pseudogene
NEJM_TB_SIG["LOC642678",5] <- "ENSG00000055609" #  KMT2C and pseudogene - pseudogene has 1 extra base matching but is v poorly expressed
NEJM_TB_SIG["LOC653778",5] <- "ENSG00000147454" # 
NEJM_TB_SIG["MIR1974",5] <- "Multiple mapped sites" #  maps to chromosome 11 and mitochondrial genome perfectly
NEJM_TB_SIG["HLA-DRB1",5] <- "ENSG00000196126" # 
NEJM_TB_SIG["LOC649210",5] <- "ENSG00000211638" # iglv8-61
NEJM_TB_SIG["HS.171481",5] <- "ENSG00000173559" # 

colnames(NEJM_TB_SIG)[c(2,4:6)] <- c("Symbol","Direction","Ensembl gene id","Illumina nuID") 
NEJM_TB_SIG[,"expressed"] <-  NEJM_TB_SIG[, "Ensembl gene id"]%in% expressed_genes
write.csv(NEJM_TB_SIG[c(-1,-3,-7)] , file = "figures/NEJM_TB_SIG_genes.csv", row.names = F)

NEJM_TB_SIG.filt <- NEJM_TB_SIG[grep("ENSG",NEJM_TB_SIG[,"Ensembl gene id"]),]
NEJM_TB_SIG.filt <- NEJM_TB_SIG.filt[NEJM_TB_SIG.filt[,"Ensembl gene id"] %in% expressed_genes,]
NEJMTB <- NEJM_TB_SIG.filt[,"Ensembl gene id"]
names(NEJMTB) <- NEJM_TB_SIG.filt[,"Symbol"]

NEJMTB.coef <- c(0,as.numeric(NEJM_TB_SIG.filt[,"Direction"] == "UP"))
NEJMTB.coef[NEJMTB.coef==0] <- -1
names(NEJMTB.coef) <- c("(Intercept)",NEJM_TB_SIG.filt[,"Ensembl gene id"])
NEJMTB.coef[1] <- 0
  

external_perfs <- list()
for(i in 1:(length(batcha_expr_sets)-1)){
  name <- paste("SWEENEY", names(batcha_expr_sets)[i],"retrained")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                             coef = NULL,
                                             genes =sweeney_genes.ens,
                                             train_comparison = list(c("TB"),c("ALL")),
                                             test_comparisons = list(TBvALL =list(c("TB"),c("ALL"))))
  name <- paste("SWEENEY", names(batcha_expr_sets)[i],"givencoefs")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                 coef = sweeney_coef,
                                                 genes =sweeney_genes.ens,
                                                 train_comparison = NULL,
                                                 test_comparisons = list(TBvALL =list(c("TB"),c("ALL"))))


}

for(i in 1:(length(batcha_expr_sets)-1)){
  print(i)
  name <- paste("NEJMTB", names(batcha_expr_sets)[i],"retrained")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = NULL,
                                                     genes =NEJMTB,
                                                     train_comparison = list(c("TB"),c("ALL")),
                                                     test_comparisons = list(TBvALL =list(c("TB"),c("ALL"))))
  name <- paste("NEJMTB", names(batcha_expr_sets)[i],"givencoefs")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = NEJMTB.coef,
                                                     genes =NEJMTB,
                                                     train_comparison = NULL,
                                                     test_comparisons = list(TBvALL =list(c("TB"),c("ALL"))))
}



for(i in 1:(length(batcha_expr_sets)-1)){
  ##
  print(c(i,length(batcha_expr_sets)-1))
  name <- paste("TWOGENE", names(batcha_expr_sets)[i],"retrained")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = NULL,
                                                     genes =twogene_genes.ens,
                                                     train_comparison = list(c(bacteria),c(viral)),
                                                     test_comparisons = list(BvALL =list(c(bacteria),c("ALL")),
                                                                             VvALL = list(c(viral),c("ALL")),
                                                                             BvV = list(c(bacteria),c(viral))))
  name <- paste("TWOGENE", names(batcha_expr_sets)[i],"givencoefs")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = twogene_coef,
                                                     genes =twogene_genes.ens,
                                                     train_comparison = NULL,
                                                     test_comparisons = list(BvALL =list(c(bacteria),c("ALL")),
                                                                              VvALL = list(c(viral),c("ALL")),
                                                                               BvV = list(c(bacteria),c(viral))))
  names(kawasaki_coef) <- c("(Intercept)",kawasaki_genes.ens)
  name <- paste("KD13", names(batcha_expr_sets)[i],"retrained")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = NULL,
                                                     genes =kawasaki_genes.ens,
                                                     train_comparison = list(c("KD"),c("ALL")),
                                                     test_comparisons = list(KDvALL =list(c("KD"),c("ALL"))))
  
  name <- paste("KD13", names(batcha_expr_sets)[i],"givencoefs")
  external_perfs[[name]] <- external_sig_performance(expr.loc =  batcha_expr_sets[[i]] ,
                                                     coef = kawasaki_coef,
                                                     genes =kawasaki_genes.ens,
                                                     train_comparison = NULL,
                                                     test_comparisons = list(KDvALL =list(c("KD"),c("ALL"))))
  
  

}





# Differential expression and enrichment ####

rnaseq_meta.loc <- as.data.frame(meta_data[names(rnaseq_broad_validation_pheno),])
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
names(profiler_res)




pal <- add.alpha(c(RColorBrewer::brewer.pal(8, "Dark2"),"black","grey","lightgrey"),0.7)
names(pal) <- sources


# save as tables 
for(i in names(profiler_res)){
  res <- profiler_res[[i]]$result
  for(j in 1:nrow(res)) res[j,"Parents"] <- paste0(res[j,"parents"][[1]],collapse = "_")
  res <- res[,c(1:13,15)]
  write.csv(res, file = paste0("figures/gprofiler_",i,".csv"))
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

pdf("figures/enrichment.bubble.pdf" ,width = 10,height = 15)
for(i in sources){
  enrichment_table <- data.frame()
  for(j in names(profiler_res)[grep("_",names(profiler_res))]){
    print(j)
    res <- profiler_res[[j]]$result
    res <- res[res[,"source"] == i,]
    res <- res[res[,"p_value"] < 0.01,]
    if(nrow(res) > 0){
      print(all(res[,"significant"]))
      print(min(res[,"intersection_size"],na.rm = T))
      res <- head(res,20)
      enrichment_table[res[,"term_id"],"term_name"] <- res[,"term_name"]
      enrichment_table[res[,"term_id"],j] <- res[,"p_value"]
    }
  }
  if(nrow(enrichment_table) > 0){
    enrichment_plot_2(enrichment_table,i)
  }
}
dev.off()

pdf("figures/enrichment.bubble_plusdirectionless.pdf" ,width = 10,height = 17)
for(i in sources){
  enrichment_table <- data.frame()
  for(j in names(profiler_res)){
    print(j)
    res <- profiler_res[[j]]$result
    res <- res[res[,"source"] == i,]
    res <- res[res[,"p_value"] < 0.05,]
    if(nrow(res) > 0){
      print(all(res[,"significant"]))
      print(min(res[,"intersection_size"],na.rm = T))
      res <- head(res,20)
      enrichment_table[res[,"term_id"],"term_name"] <- res[,"term_name"]
      enrichment_table[res[,"term_id"],j] <- res[,"p_value"]
    }
  }
  if(nrow(enrichment_table) > 0){
    enrichment_plot_2(enrichment_table,i)
  }
}
dev.off()

