

# init ####
source("scripts/functions.R")
loadlibs()
plotdir <- "plots/"
load("data_objects/pre-processed-fold-assigned.RData")

indicator_matrix  <- class.ind(factor(training_phenotypes[,"group"]))
rownames(indicator_matrix) <- rownames(training_phenotypes)

cost_base <- as.data.frame(matrix(nrow = ncol(indicator_matrix),ncol = ncol(indicator_matrix)))
colnames(cost_base) <- colnames(indicator_matrix)
rownames(cost_base) <- colnames(indicator_matrix)
cost_base[] <- 0
cost_matrix.uniform <- cost_base
cost_matrix.uniform[] <- 1
class_weighted_cost <- cost_base
for(i in 1:ncol(class_weighted_cost)){
  class_weighted_cost[,i] <- class_weights[rownames(class_weighted_cost),"cost"]
  class_weighted_cost[i,i] <- 0
  cost_matrix.uniform[i,i] <- 0
}




# SEED CV ##########

## multinomial LASSO and OVA ####
# for combinations of cost and imbalance weighted learning 
# saves to load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))

base_training_phenotypes <- training_phenotypes
mclapply(1:10 , run_seed_cv,mc.cores = 10)


## performance over all seeds  ####

# extract performance  , confusions etc from the seed cv run data files:
runseeds <-   c(1:10)
seed_performance <- c( lapply(runseeds ,calc_performances))# , mc.cores = 3))

seed_modelsize <- c( lapply(runseeds,extract_seed_model_sizes))#,mc.cores = 5))

seed_performance.merged <- merge_seed_runs(seed_performance)
tmp <- seed_performance.merged




## seed rdata -> roc curves   ####

seedvalue <- 1
load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))

all_rocs <- list()

# one versus all
all_rocs[["OVA"]] <- get_rocs(prediction_matrix.onevsall,training_phenotypes_indicator)
# multinomial models (with cost and/or class imbalance respectively (T/F))
all_rocs[["M_FF"]] <- get_rocs( prediction_matrix.multinomial_FF ,training_phenotypes_indicator)
all_rocs[["M_TF"]] <- get_rocs( prediction_matrix.multinomial_TF ,training_phenotypes_indicator)
all_rocs[["M_FT"]] <- get_rocs(prediction_matrix.multinomial_FT  ,training_phenotypes_indicator)
all_rocs[["M_TT"]] <- get_rocs( prediction_matrix.multinomial_TT ,training_phenotypes_indicator)

seedvalue <- 1
shrunk_predmats <- shrink_MTT(seedvalue)

for(i in names(shrunk_predmats)){
  all_rocs[[paste("M_TT_",i,sep="")]] <- get_rocs(shrunk_predmats[[i]]  ,training_phenotypes_indicator)  
}

training_phenotypes_indicator <- training_phenotypes_indicator[rownames(training_phenotypes),]
tmp <- seed_performance[[1]]
for(i in names(seed_performance[[1]])){
  if(i!="seed"){
    for(j in 2:length(seed_performance)){
      tmp[[i]]$meansquared <- c(tmp[[i]]$meansquared , seed_performance[[j]][[i]]$meansquared) 
      tmp[[i]]$fp <- rbind(tmp[[i]]$fp , seed_performance[[j]][[i]]$fp)
      tmp[[i]]$fn <- rbind(tmp[[i]]$fn , seed_performance[[j]][[i]]$fn)
    } 
  }
}

## to show that ova gives more correlated genes : #####
cvfit <- all_fits[[1]]$mn_TT
path <- "flu"
coefficients <- coef(cvfit,s="lambda.1se")

probes <- rownames(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
probes <- probes[probes!="(Intercept)"]

pathogens <- names(coefficients)
mean_exprs <- data.frame(matrix(nrow = length(probes) , ncol = length(pathogens)))
colnames(mean_exprs) <- pathogens
rownames(mean_exprs) <- probes

for(path in names(coefficients)){
  cases <- rownames(training_phenotypes)[training_phenotypes[,"group"] == path]
  mean_exprs[probes , path ] <- rowMeans(expression_filt[probes,cases])
}


corel <-  cor(mean_exprs)
png("plots/path_cor_heatmap_M_TT.png",width = 800 , height = 800)
superheat(corel,
          X.text = round(corel*100,digits = 0),
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,scale = F,
          row.dendrogram = TRUE,   bottom.label.text.angle = 90)
dev.off()

corel.gene <-  cor(t(expression_filt[probes,]))
png("plots/signature_gene_cor_heat_multinom.png",width = 1500 , height = 1500)
superheat(corel.gene,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,scale = F,
          bottom.label.text.angle = 90)
dev.off()

ova_probes <- c()
for(path in pathogens){
  coefs <- as.matrix(coef(all_fits[[1]]$ova[[path]],s = "lambda.1se"))
  probes <- rownames(coefs)[coefs!=0]
  probes <- probes[probes!="(Intercept)"]
  ova_probes  <- c(ova_probes , probes)
}

ova_probes <- unique(ova_probes )

corel.gene2 <-  cor(t(expression_filt[ova_probes,]))
png("plots/signature_gene_cor_heat_ova.png",width = 2000 , height = 2000)
superheat(corel.gene2,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          bottom.label.text.angle = 90)
dev.off()


corel.vec.gene <- as.vector(corel.gene)
corel.vec.gene <- corel.vec.gene[corel.vec.gene!= 1]
corel.vec.gene2 <- as.vector(corel.gene2)
corel.vec.gene2 <- corel.vec.gene2[corel.vec.gene2!=1]

dat <- data.frame(xx = c( corel.vec.gene,corel.vec.gene2 ),
                  yy = c(rep("mn",length(corel.vec.gene)) , rep("ova",length(corel.vec.gene2 ))))


png("plots/signature_gene_correlation_density.png",width = 500 , height = 500)
ggplot(dat,aes(x=xx)) + xlab("correlation") +
  geom_density(data=subset(dat,yy == 'mn'),fill = "red" , adjust=0.5,alpha = 0.5) +
  geom_density(data=subset(dat,yy == 'ova'),fill = "blue", adjust=0.5,alpha = 0.5)
dev.off()


# Seed runs for Hybrid / two-stage methods ##############

## ridge on lassos ####
for(seedvalue in 1:10){
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  ridge_out <- mclapply(1:10 , cv_hybrid , mc.cores = 5,refit_method = "ridge")
  save(ridge_out , file = paste("data_objects/seed_cv/lasso_ridge_seed",seedvalue,".RData",sep=""))
}


## get mse + ms 
pooled_lasso_ridges <- list()
for(seedvalue in 1:10){
  print(seedvalue)
  load(paste("data_objects/seed_cv/lasso_ridge_seed",seedvalue,".RData",sep=""))
  lasso_ridge_pm <- prediction_matrix.multinomial_TT # prediction matrix over all folds
  lasso_ridge_pm[] <- 0
  ms <- c() # model size
  selected_genes <- list()
  
  for(i in 1:10){ # internal folds
    pm <- ridge_out[[i]]$ridge_pred     
    genes <- ridge_out[[i]]$selected_genes     
    lasso_ridge_pm[rownames(pm) , colnames(pm)] <- pm
    ms <- c(ms , length(genes))
    selected_genes[[i]] <- genes
  }
  
  confusion <- make_confusion(lasso_ridge_pm ,
                              training_phenotypes_indicator[rownames(lasso_ridge_pm),
                                                            colnames(lasso_ridge_pm)],mod = "3dconfusion") 
  # 3d confusion is a 3 dimensional array where the intersection of predicted and true values are split by sample along the 3rd dimension
  # confusion$values is for predicted probabilities / confusion$maxvals is where the predicted class is 1 and all others are 0
  
  ew.cost <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost"]
  ew.cost_imb <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost+imbalance"]
  perfs <- list()
  nam <- "ridge_lasso"
  
  perfs[[paste(nam,".u.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ) # predicted probabilities
  perfs[[paste(nam,".u.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ) # discrete predictions
  perfs[[paste(nam,".ew.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".ew.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".cw.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost)
  perfs[[paste(nam,".cw.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost)
      
  errs <- performance_weighted(confusion,mod ="val" ,example_weights = ew.cost)
  
  pooled_lasso_ridges[[seedvalue]] <- c(list(
    genes = selected_genes,
    pm = lasso_ridge_pm,
    ms = mean(ms)
  ), perfs )
  
}


## Relaxed lasso ####
for(seedvalue in 1:10){
  print(seedvalue)
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  hybrid_cv_obj <- mclapply(1:10 , cv_hybrid , mc.cores = 5,refit_method = "lasso")
  save(hybrid_cv_obj , file = paste("data_objects/seed_cv/relaxed_lasso_fit_seed",seedvalue,"_TT.RData",sep=""))
}

## relaxed lasso - get mse + ms

pooled_relaxed_lasso <- list()
for(seedvalue in 1:10){
  print(seedvalue)
  load(paste("data_objects/seed_cv/relaxed_lasso_fit_seed",seedvalue,"_TT.RData",sep=""))
  lasso_ridge_pm <- prediction_matrix.multinomial_TT
  lasso_ridge_pm[] <- 0
  ms <- c()
  selected_genes <- list()
  for(i in 1:10){
    pm <- hybrid_cv_obj[[i]]$ridge_pred     
    genes <- hybrid_cv_obj[[i]]$selected_genes     
    lasso_ridge_pm[rownames(pm) , colnames(pm)] <- pm
    ms <- c(ms , length(genes))
    selected_genes[[i]] <- genes
  }
  
  confusion <- make_confusion(lasso_ridge_pm , training_phenotypes_indicator[rownames(lasso_ridge_pm),colnames(lasso_ridge_pm)],mod = "3dconfusion") 
  
  
  ew.cost <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost"]
  ew.cost_imb <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost+imbalance"]
  perfs <- list()
  nam <- "relaxed_lasso"
  
  perfs[[paste(nam,".u.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" )
  perfs[[paste(nam,".u.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" )
  perfs[[paste(nam,".ew.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".ew.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".cw.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost)
  perfs[[paste(nam,".cw.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost)
  
  errs <- performance_weighted(confusion,mod ="val" ,example_weights = ew.cost)
  
  pooled_relaxed_lasso[[seedvalue]] <- c(list(
    genes = selected_genes,
    pm = lasso_ridge_pm,
    ms = mean(ms)
  ), perfs )
  
}



# check for null mses
for(j in 1:10){
  print(j)
  load(paste("data_objects/seed_cv/lasso_ridge_seed",j,".RData",sep=""))
  for(i in 1:10){
    if(is.null(ridge_out[[i]]$mse_plot_dat$mse.l)) print(c(i,j,"l"))
    if(is.null(ridge_out[[i]]$mse_plot_dat$mse.r)) print(c(i,j,"r"))
  }
}

for(j in 1:10){
  print(j)
  load(paste("data_objects/seed_cv/run_number_",j,".RData",sep=""))
  for(i in 1:10){
    if(is.null(hybrid_cv_obj[[i]]$mse_plot_dat$mse.l)) print(c(i,j,"l"))
    if(is.null(hybrid_cv_obj[[i]]$mse_plot_dat$mse.r)) print(c(i,j,"r"))
  }
}



# ridge on ova ####

ova_ridge <- list()
for(seedvalue in 1:10){
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  for(i in 1:10) training_phenotypes[all_fits[[i]]$test_samples,"fold"] <- i # fix folds to what was selected be seed
  ova_ridge_pm <- prediction_matrix.multinomial_TT
  siglist <- list()
  fits <- list()
  for(foldn in 1:10){
    print(paste("seed =",seedvalue,"fold =", foldn))
    training_samples <- all_fits[[foldn]]$train_samples
    test_samples <- all_fits[[foldn]]$test_samples
    genes_in_ova <- unique(as.character(all_fits[[foldn]]$ova$sig_genes))
    siglist[[foldn]] <- genes_in_ova
    # do fit
    disease_set <- rownames(class_weights)
    weights.local <- class_weights[training_phenotypes[training_samples,"group"],"cost+imbalance"]  
    inner_folds <- training_phenotypes[training_samples,"fold"]
    inner_folds <- as.numeric(as.factor(inner_folds)) # so that 1: ... 
    print(c(length(inner_folds) , length(weights.local), length(training_samples)))
    ridge_fit <- cv.glmnet(t(expression_filt[genes_in_ova ,training_samples]),
                           y = training_phenotypes_indicator[training_samples,disease_set], lambda = exp(seq(-8,2,by = 0.1)) 
                           ,foldid = inner_folds,weights = weights.local 
                           ,type.measure = "mse",alpha = 0,family ="multinomial")
    fits[[foldn]] <- ridge_fit
    print(c("fitted ridge"))
    p2 <- predict(ridge_fit,newx = t(expression_filt[genes_in_ova,test_samples]) , s = "lambda.min",type = "response")[,,1]
    ova_ridge_pm[rownames(p2),colnames(p2)] <- p2
  }
  ova_ridge[[seedvalue]]   <- list( pm = ova_ridge_pm,
                                    fits = fits,
                                    siglist = siglist)
}

save( ova_ridge
     ,file ="data_objects/seed_cv/ridge_on_ova.seed.RData")


pooled_ridgeova <- list()
for(seedvalue in 1:10){
  print(seedvalue)
  ridge_pm <- prediction_matrix.multinomial_TT
  ridge_pm[] <- 0
  ms <- c()
  selected_genes <- list()

  ridge_pm <- ova_ridge[[seedvalue]]$pm     
  genes <- ova_ridge[[seedvalue]]$siglist     
    ms <- lapply(genes ,  FUN = length) ## ?? 

  confusion <- make_confusion(ridge_pm , training_phenotypes_indicator[rownames(ridge_pm),colnames(ridge_pm)],mod = "3dconfusion") 
  
  ew.cost <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost"]
  ew.cost_imb <- class_weights[training_phenotypes[dimnames(confusion$values)[[3]] , "group"] ,"cost+imbalance"]
  perfs <- list()
  nam <- "ova_ridge"
  
  # unweighted performance
  perfs[[paste(nam,".u.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" )
  perfs[[paste(nam,".u.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" )
  # weighted performance  (cost+imbalance or cost)
  perfs[[paste(nam,".ew.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".ew.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost_imb)
  perfs[[paste(nam,".cw.v",sep = "")]] <- performance_weighted(confusion  , mod ="val" ,example_weights = ew.cost)
  perfs[[paste(nam,".cw.m",sep = "")]] <- performance_weighted(confusion  , mod ="max" ,example_weights = ew.cost)
  
  errs <- performance_weighted(confusion,mod ="val" ,example_weights = ew.cost)
  pooled_ridgeova[[seedvalue]] <- c(list(
    genes = ova_ridge[[seedvalue]]$siglist,
    pm = ridge_pm,
    ms = mean(unlist(ms))
  ), perfs )
  
}




# mse lambda plots ####

seedvalue <- 1
fold_value  <- 1

load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
plot(all_fits[[fold_value]]$mn_TT, ylim = c(0,1))

load(paste("data_objects/seed_cv/ridge_fit_seed",seedvalue,"_TT.RData",sep=""))
mse_plot_dat.lr <- ridge_out[[fold_value]]$mse_plot_dat

load(paste("data_objects/seed_cv/relaxed_lasso_fit_seed",seedvalue,"_TT.RData",sep=""))
mse_plot_dat.relaxedlasso <- hybrid_cv_obj[[fold_value]]$mse_plot_dat


plot(mse_plot_dat.lr$mse.l , all_fits[[fold_value]]$mn_TT$cvm, ylim = c(0,1),xlim = c(0,1))
points(mse_plot_dat.lr$mse.r , lasso_fit_cv$cvm)
points(mse_plot_dat.relaxedlasso$mse.l , lasso_fit_cv$cvm)
points(mse_plot_dat.relaxedlasso$mse.r , lasso_fit_cv$cvm)
points(c(0,1),c(0,1), type = "l")

plot(log(lasso_fit_cv$lambda) ,mse_plot_dat.lr$mse.l, ylim = c(0,1),xlim = c(-6,-1.5))
points(log(lasso_fit_cv$lambda) , lasso_fit_cv$cvm)
points(log(lasso_fit_cv$lambda) , mse_plot_dat.lr$mse.r)
points(log(lasso_fit_cv$lambda) , mse_plot_dat.low.relaxed_lasso$mse.l)

plot_ridge_cv(mse_pltdat = mse_plot_dat.lr ,lasso_fit_cv = all_fits[[fold_value]]$mn_TT,xvar = "lambda",
                      mode = "LOW",xrange = c(-6,-1.5),colvec = brewer.pal(2,"Dark2"),lw = 2)

plot(all_fits[[fold_value]]$mn_TT)
df <- coef(all_fits[[fold_value]]$mn_TT)[[1]]
length(df[df!=0])

points(log(lasso_fit_cv$lambda),mse_plot_dat.lr$mse.l)
thresh = mse_plot_dat.lr$mse.lp[mse_plot_dat.lr$mse.l == min(mse_plot_dat.lr$mse.l)]
abline(v=max(log(lasso_fit_cv$lambda)[mse_plot_dat.lr$mse.l <= thresh]))

for(outer_fold in 1:10){
  for(i in 2:9){
      print(which(ridge_out[[outer_fold]]$ridge_sheets[[i]] == min(ridge_out[[outer_fold]]$ridge_sheets[[i]],na.rm = T) , arr.ind = T ))
  }}


# contour plot
for(outer_fold in 1:10){
  mat <-ridge_out[[outer_fold]]$ridge_sheets[[1]]
  for(i in 2:9){
  
    mat <- mat + ridge_out[[outer_fold]]$ridge_sheets[[i]]
  }
  mat <- mat / 10
  p2 <- plot_ly(y=log(as.numeric(rownames( mat))) , x = log(as.numeric(colnames( mat)))) %>%
    add_contour(z = ~as.matrix( mat)) %>%
    #  add_markers(x=X,y=log(Y),z=W,size=(Z), marker = list(symbol = 'circle', sizemode = 'diameter'),sizes=c(1,50)) %>%  
    layout(
      title = " mse",
      scene = list(
        xaxis = list(title = "lambda-lasso"),
        yaxis = list(title = "lambda-ridge"),
        zaxis = list(title = "MSE")
      )) 
  export(p2,paste("plots/ridge_sheet_fold_",outer_fold,".png" ,sep =""))
}


par(mar=c(5.1,4.1,4.1,2.1))
for(outer_fold in 1:10){
  png(paste("plots/ridge_cv_fold_",outer_fold,".png" , sep = "" ),height = 1000,width = 1600)
  mse_pltdat <- ridge_out[[outer_fold]]$mse_plot_dat
  lasso_fit_cv <- all_fits[[outer_fold]]$mn_TT
  plot_ridge_cv(mse_pltdat ,lasso_fit_cv ,ridge_lambdas = as.numeric(rownames(ridge_out[[1]]$ridge_sheets[[1]])))
  dev.off()
}



sig_overlaps <- matrix(nrow = 10 , ncol = 10)
for(i in 1:10){
  for(j in 1:10){
    sig_overlaps[i,j] <- length(intersect(ridge_out[[i]]$selected_genes , ridge_out[[j]]$selected_genes))    
  }
}
fc <- c()
for( i in 1:10){
  fc <- c(fc, ridge_out[[i]]$selected_genes)  
}

superheat(sig_overlaps,
          X.text = sig_overlaps,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,scale = F,
          row.dendrogram = TRUE,   bottom.label.text.angle = 90)


