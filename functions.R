
# handy functions ####

loadlibs <- function(){
  library(magrittr,quietly = T)
  library(pracma , quietly = T)
  library(Hmisc,quietly = T)
  library(plyr,quietly = T)
  library(tidyr,quietly = T)
  library(rsvd,quietly = T)
  library(stats,quietly = T)
  library(limSolve,quietly = T)
  library(abind,quietly = T)
  library(beeswarm,quietly = T)
  library(DESeq2,quietly = T)
  library(factoextra,quietly = T)
  library(FactoMineR,quietly = T)
  library(limma,quietly = T)
  library(htmlwidgets,quietly = T)
  library(cluster,quietly = T)
  library(ggfortify,quietly = T)
  library(PRROC,quietly = T)
  library(gplots,quietly = T)
  library(glmnet,quietly = T)
  library(beeswarm,quietly = T)
  library(circlize,quietly = T)
  library(parallel,quietly = T)
  library(nnet,quietly = T)
  require(doMC,quietly = T)
  library(reshape,quietly = T)
  library(plotly,quietly = T)
  library(dplyr,quietly = T)
  library(ROCR,quietly = T)
  library(pROC,quietly = T)
  library(plotly,quietly = T)
  library(GGally,quietly = T)
  registerDoMC(cores=10)
  options(pROCProgress = list(name = "none"))
  require("VennDiagram")
  
}


`%ni%` <- Negate(`%in%`)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


get_coef<- function(cvfit,  path = "TB", nse = 1 , plot = F,inlambda = NULL){ # nse gives number of standard errors from min to take lambda
  
  se <- cvfit$cvsd[match(cvfit$lambda.min , cvfit$lambda)]
  lambda <- max(cvfit$lambda[cvfit$cvm <= cvfit$cvm[match(cvfit$lambda.min , cvfit$lambda)] + nse * se  ])
  
  if(!is.null(inlambda)){
    lambda <- inlambda
  }
  coefficients <- coef(cvfit,s=lambda )
  
  probes <- rownames(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
  out <- as.numeric(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
  intercept <- out[probes=="(Intercept)"]

  out <- out[probes!="(Intercept)"]
  probes <- probes[probes!="(Intercept)"]
  names(out) <- probes
  if(plot){
    plot(cvfit)
    abline(v = log(lambda))
    
  }
  return(out)
}

rowpaste <- function(mat){
  out <- c()
  for(i in 1:nrow(mat)){
    out <- c(out, paste(mat[i,],collapse = " "))
  }
  out
}


pretty_plot_area <- function(x_lab ="" , y_lab ="" ,main = "",
                             ytick = c(0,2,0.1) ,ytickn = NULL,
                             xtick= c(0,1,0.1) , xtickn = NULL,
                             show_x_tick = T,show_y_tick = T,
                             show_x = T,show_y = T,background = "white",
                             cols = c("lightgray","white"),text_col = "grey20",
                             xbuffer = 0,ybuffer = 0,margins = c(8, 8, 2, 1),
                             add_rect = T,plot_outside_margin = F, axcex = 1){
  par(xpd=F)
  xtick2 <- NULL
  ytick2 <- NULL
  if(!is.null(ytickn))  ytick2 <- seq(ytick[1],ytick[2],by = ytickn )
  if(!is.null(xtickn))  xtick2 <- seq(xtick[1],xtick[2],by = xtickn )
  
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
  if(show_x)  mtext(side = 1, text = xtick, at = xtick, 
                    col = text_col, line = 1, cex = 1.2,las =1)
  if(show_y) mtext(side = 2, text = round(ytick, digits = 10), at = ytick, 
                   col = text_col, line = 1, cex = 1.2,las =2)
  
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


# get performance _ PR, ROC, F1, heierarchical P-R-F1 ####

get_rocs <- function(pred.resp,phenotypes_indicator,toplevel = T){
  rocs <- list()
  colours <- rainbow(ncol(pred.resp))
  names(colours) <- colnames(pred.resp)
  addval <- F
  for(path in colnames(pred.resp)){
    cases <- rownames(phenotypes_indicator)[phenotypes_indicator[,path]==1]
    controls <- rownames(phenotypes_indicator)[phenotypes_indicator[,path]!=1]
    r <- pROC::roc(cases = pred.resp[cases,path],controls = pred.resp[controls,path], ci  = T)
    rocs[[path]] <- r
    if(toplevel){
      for(tl in c("bacterial","viral","inflammatory")){
        controls <- rownames(phenotypes_indicator)[phenotypes_indicator[,tl]==1]
        controls <- controls[!controls%in%cases]
        if(length(controls) > 0){
          r <- pROC::roc(cases = pred.resp[cases,path],controls = pred.resp[controls,path],ci  = T)
          rocs[[paste(path,"-vs-",tl,sep ="")]] <- r
        }
      }
    }
  }
  return(rocs)  
}


pr_curve <- function(cases,controls){
  allvals <-  c(cases,controls)
  allvals <- allvals[order(allvals)]
  curve <- data.frame(matrix(ncol = 4))
  colnames(curve) <- c("tp","fp","tn","fn")  
  for(thresh in allvals){
    tp <- length(cases[cases >= thresh])
    fp <- length(controls[controls >= thresh])
    tn <- length(controls[controls < thresh])
    fn <- length(cases[cases < thresh])
    curve[as.character(thresh),1:4] <- c(tp,fp,tn,fn)  
  }
  curve <- curve[!is.na(curve[,1]),]
  curve[,"precision"] <- curve[,"tp"] / (curve[,"tp"] + curve[,"fp"])
  curve[,"recall"] <- curve[,"tp"] / (curve[,"tp"] + curve[,"fn"])
  curve[,"f1"] <- 2*curve[,"tp"]  / (2*curve[,"tp"]  + curve[,"fp"]  + curve[,"fn"] )
  return(curve)
}

get_PR <- function(pred.resp,phenotypes_indicator,toplevel =T){
  curves <- list()
  colours <- rainbow(ncol(pred.resp))
  names(colours) <- colnames(pred.resp)
  addval <- F
  for(path in colnames(pred.resp)){
    cases <- rownames(phenotypes_indicator)[phenotypes_indicator[,path]==1]
    controls <- rownames(phenotypes_indicator)[phenotypes_indicator[,path]!=1]
    pr <- pr.curve(scores.class0 = pred.resp[cases,path] , scores.class1 =  pred.resp[controls,path] , curve = T)
    curves[[path]] <- pr
    if(toplevel){
      for(tl in c("bacterial","viral","inflammatory")){
        controls <- rownames(phenotypes_indicator)[phenotypes_indicator[,tl]==1]
        controls <- controls[!controls%in%cases]
        if(length(controls) > 0){
          pr <- pr.curve(scores.class0 = pred.resp[cases,path] , scores.class1 =  pred.resp[controls,path] , curve = T)
          curves[[paste(path,"-vs-",tl,sep ="")]] <- pr
        }
      }
    }
  }
  return(curves)  
}

get_auc_table_from_rocs <- function(roc_obj){
  pathogens <- colnames(test_phenotypes_indicator)[!colnames(test_phenotypes_indicator)%in%c("bacterial","viral","inflammatory")]
  colours <- rainbow(length(pathogens))
  auc_table  <- as.data.frame(matrix(rep(0,5)))
  for(i in 1:length(roc_obj)){
    r <- roc_obj[[i]]
    pahthogen_A <- gsub("-vs-.*","",names(roc_obj)[i])
    auc_table[,names(roc_obj)[i]] <- c(length(r$cases),length(r$controls), r$auc,r$ci[1],r$ci[3] )
  }
  rownames(auc_table) <- c("cases","controls","AUC","CI.l","CI.h")
  return(auc_table[,-1])
}


get_performance <- function(pred_resp){
  out <- list()
  if(all(rownames(pred_resp) %in% rownames(training_phenotypes_indicator))){
    indicator <- training_phenotypes_indicator
  }else if(all(rownames(pred_resp) %in% rownames(test_phenotypes_indicator))){
    indicator <- test_phenotypes_indicator
  }else{
    break
  }
  out[["ROCs"]] <- get_rocs(pred_resp ,indicator)
  out[["AUCs"]] <- get_auc_table_from_rocs(out[["ROCs"]])
  out[["PR"]] <- get_PR(pred_resp ,indicator)
  out[["F1"]] <- out[["PR"]] 
  for(i in names(out[["PR"]])){
    pr <- out[["PR"]][[i]]$curve
    out[["F1"]][[i]]   <- 2 * ((pr[,1] * pr[,2]) / (pr[,1] + pr[,2]))
  }
  out[["PR_mod"]] <- get_pr_micro_curves(pred_resp , indicator) 
  # micro averaged PR, class specific PR, multilabel PR, multilabel hierarchical PR
  return(out)  
}


list_mod <- function(list_obj){ # turn list items into dataframes recursively
  list_obj.out <- list_obj
  if(class(list_obj) == "list"){
    for(sublist in names(list_obj)){
      list_obj.out[[sublist]] <- list_mod(list_obj[[sublist]])    
    }
  }else if(class(list_obj) == "data.frame"){
    list_obj.out <- as.data.frame(list_obj)   
  }else if(class(list_obj) == "matrix"){
    list_obj.out <- as.data.frame(list_obj)   
  }
  return(list_obj.out)
}


list_mod.perf <- function(list_obj){ # calc performance from prediction matrix list
  list_obj.out <- list_obj
  print(class(list_obj))
  if(class(list_obj) == "list"){
    print(names(list_obj))
    for(sublist in names(list_obj)){
      if(sublist%ni%c( "leftout","respiratory")){
        list_obj.out[[sublist]] <- list_mod.perf(list_obj[[sublist]])     
      }else{
        list_obj.out[[sublist]] <- NULL
      }
      
    }
  }else if(class(list_obj) == "data.frame"){
    list_obj.out <- get_performance(list_obj)
  }else{
    print("==")
  }
  return(list_obj.out)
}



augment_set <- function(set){ # add gram stain ?? 
  set.aug <- set
  for(class in set){
    if(class %in% bacteria) {
      set.aug <- c(set.aug , "bacterial")
    }else if(class %in% viral){
      set.aug <- c(set.aug , "viral")
    }else if(class %in% inflammatory){
      set.aug <- c(set.aug , "inflammatory")
    }
  }
  set.aug
}


get_pr.ml <- function(y,  y_hat , heierarchical = F){ # multilabel precision recalls
  recalls <- c()
  precisions <- c()
  recalls.aug <- c()
  precisions.aug <- c()
  specificities <- c()
  for(i in rownames(y_hat)){
    y_hat.i <- colnames(y_hat)[y_hat[i,] == 1]
    y.i <- colnames(y)[y[i,]==1]
    if(heierarchical== T){ # augment the sets to contain tree info
      y_hat.i <- augment_set(y_hat.i)
      y.i <- augment_set(y.i)
    }
    precision.i <- length(intersect(y.i , y_hat.i)) / length(y_hat.i) # see  http://arxiv.org/abs/1306.6802v2
    recall.i <- length(intersect(y.i , y_hat.i)) / length(y.i)
    
    minus_y_hat.i <- colnames(y_hat)[y_hat[i,] != 1]
    minus_y.i <- colnames(y)[y[i,]!=1]
    specificity <- length(intersect(minus_y.i, minus_y_hat.i)) / length(minus_y.i)
    
    recalls <- c(recalls , recall.i)
    precisions <- c(precisions , precision.i)
    specificities <- c(specificities , specificity)
    
  }
  precisions[is.na(precisions)] <- 0
  precision <- mean(precisions)
  recall <- mean(recalls)
  specificity <- mean(specificities)
  return(list(precision = precision,  recall = recall , f = 2*precision*recall / (precision + recall), specificity = specificity))
}


get_pr_micro <- function(y,  y_hat , heierarchical = F){ # multilabel precision recalls
  
  tp <- colSums(y == 1 & y_hat == 1 )
  fp <- colSums(y == 0 & y_hat == 1 )
  tn <- colSums(y == 0 & y_hat == 0 )
  fn <- colSums(y == 1 & y_hat == 0 )
  
  precision.class  = tp / (tp + fp)
  recall.class  = tp / (tp + fn)
  f1.class = 2*tp / (2*tp + fp + fn)
  
  sensitivity.class <- tp / (tp + fn)
  specificity.class <- tn / (tn + fp)
  
  return(list(  precision.class  = precision.class,
                recall.class  = recall.class,
                f1.class = f1.class,
                sensitivity.class = sensitivity.class,
                specificity.class = specificity.class,
                
                precision.macro = mean(precision.class,na.rm =T),
                recall.macro = mean(recall.class),
                f1.macro = mean(f1.class),
                sensitivity.macro = mean(sensitivity.class ),
                specificity.macro = mean(specificity.class,na.rm =T),
                
                precision.micro  = sum(tp) / (sum(tp) + sum(fp)),
                recall.micro  = sum(tp) / (sum(tp) + sum(fn)),
                f1.micro = 2*sum(tp) / (2*sum(tp) + sum(fp) + sum(fn)),
                sensitivity.micro = sum(tp) / sum(c(tp, fn)),
                specificity.micro = sum(tn) / sum(c(tn, fp))
                
  ))
}


get_pr_micro_curves <- function(pred_resp , indicator){
  out <- list()
  pred_resp.max <- pred_resp
  pred_resp.max[pred_resp == rowMaxs(as.matrix(pred_resp))] <- 1
  pred_resp.max[pred_resp != rowMaxs(as.matrix(pred_resp))] <- 0
  y_hat <- pred_resp.max
  y <- indicator[rownames(y_hat),colnames(y_hat)] 
  out[["single_class_hierarchical"]] <- get_pr.ml(y,y_hat,heierarchical = T)
  # multi-class
  out[["single_class"]] <- get_pr_micro(y,y_hat) # 
  thresh_range <- seq(min(pred_resp)-0.01,max(pred_resp)+0.01,by = 0.01)
  
  for(thresh in thresh_range){ # single threshold
    y_hat <- pred_resp
    y_hat[pred_resp >= thresh ] <- 1
    y_hat[pred_resp < thresh ] <- 0
    y <- indicator[rownames(y_hat),colnames(y_hat)] 
    if(thresh == thresh_range[1]){
      perf.list <- get_pr_micro(y,y_hat)
      tmp <- get_pr.ml(y,y_hat)
      perf.list[["multilabel.precision"]] <- tmp$precision
      perf.list[["multilabel.recall"]] <- tmp$recall
      perf.list[["multilabel.specificity"]] <- tmp$specificity
      perf.list[["multilabel.f"]] <- tmp$f
      
      tmp <- get_pr.ml(y,y_hat,heierarchical = T)
      perf.list[["multilabel_heierarchical.precision"]] <- tmp$precision
      perf.list[["multilabel_heierarchical.recall"]] <- tmp$recall
      perf.list[["multilabel_heierarchical.specificity"]] <- tmp$specificity
      perf.list[["multilabel_heierarchical.f"]] <- tmp$f
    }
    
    perf <- get_pr_micro(y,y_hat)
    tmp <- get_pr.ml(y,y_hat)
    perf[["multilabel.precision"]] <- tmp$precision
    perf[["multilabel.recall"]] <- tmp$recall
    perf[["multilabel.specificity"]] <- tmp$specificity
    perf[["multilabel.f"]] <- tmp$f
    
    tmp <- get_pr.ml(y,y_hat,heierarchical = T)
    perf[["multilabel_heierarchical.precision"]] <- tmp$precision
    perf[["multilabel_heierarchical.recall"]] <- tmp$recall
    perf[["multilabel_heierarchical.specificity"]] <- tmp$specificity
    perf[["multilabel_heierarchical.f"]] <- tmp$f
    for(i in names(perf)){
      perf.list[[i]] <- rbind(perf.list[[i]] ,perf[[i]])
    } 
  }
  return(perf.list)
}

get_performance.rnaseq <- function(pred_resp,indicator){
  out <- list()
  out[["ROCs"]] <- get_rocs(pred_resp ,indicator,toplevel = F)
  out[["AUCs"]] <- get_auc_table_from_rocs(out[["ROCs"]])
  out[["PR"]] <- get_PR(pred_resp ,indicator,toplevel = F)
  out[["F1"]] <- out[["PR"]] 
  for(i in names(out[["PR"]])){
    pr <- out[["PR"]][[i]]$curve
    out[["F1"]][[i]]   <- 2 * ((pr[,1] * pr[,2]) / (pr[,1] + pr[,2]))
  }
  out[["PR_mod"]] <- get_pr_micro_curves(pred_resp , indicator) 
  # micro averaged PR, class specific PR, multilabel PR, multilabel hierarchical PR
  out
}

external_sig_performance <- function(expr.loc ,genes ,coef = NULL ,train_comparison = list(c("TB"),c("ALL")),
                                     test_comparisons = list()){
  if(is.null(coef)){ # refit coefficients
    expr <- as.data.frame(t(expr.loc[genes,train_rnaseq] ))
    expr[,"class"] <- rnaseq_validation_pheno[train_rnaseq] %in% train_comparison[[1]]
    expr[,"class"] <- as.numeric(expr[,"class"])
    if(train_comparison[[2]]!="ALL"){ # if not ova is a vs b
      expr <- expr[rnaseq_validation_pheno[train_rnaseq] %in% unlist(train_comparison),]
    }

    fit  <- glm(class ~ . , data = expr,maxit = 10000)
    coef <- coef(fit)
  }
  pred <- coef["(Intercept)"] + colSums(expr.loc[names(coef)[-1],] * coef[-1])
  out <- list()
  for(i in 1:length(test_comparisons)){
    comparison <-  test_comparisons[[i]]
    # test
    cases <- test_rnaseq[rnaseq_validation_pheno[test_rnaseq] %in% comparison[[1]]]
    if(comparison[[2]]=="ALL"){
      controls <- test_rnaseq[!rnaseq_validation_pheno[test_rnaseq] %in% comparison[[1]]]  
      controls <- controls[rnaseq_validation_pheno[controls]!="HC"]
    }else if(comparison[[2]]=="OD"){
      controls <- rnaseq_ODs
    }else{
      controls <- test_rnaseq[rnaseq_validation_pheno[test_rnaseq] %in% comparison[[2]]]
    }
    
    out[[names(test_comparisons)[[i]]]] <- list(comparison = comparison,
                                                roc  = roc(cases = pred[cases] , controls = pred[controls],ci =T),
                                                roc.smooth  = roc(cases = pred[cases] , controls = pred[controls],
                                                                  ci =T,smooth = T,smooth.method = "binormal"),
                                                curve = pr_curve(cases = pred[cases] , controls = pred[controls]),
                                                vals.cases = pred[cases],
                                                vals.controls = pred[controls])
  }
  return(out)
}




# plotting functions #####################################################################


bees <- function(pltdat,main =""){
  boxplot(vals ~ xvals  , data = pltdat, main = main, # mappingInfo[i,"Symbol"], 
          outline = FALSE ,ylab = "expression" ,at = unique(pltdat[,"xvals"]),xaxt = 'n')
  beeswarm(vals ~ xvals ,data = pltdat, corral = "wrap" ,at = unique(pltdat[,"xvals"]),add =T,
           pwcol = pltdat[,"cols"] , pch = 16 ,corralWidth=.5) # , data = pltdat , xaxt = 'n',add = F )#,)
  axis(side = 1 , at = unique(pltdat[,"xvals"]) , labels = unique(pltdat[,"labs"]) ,las = 2) # rep(comp_classes,each =1)
  
}


plot_radar <- function(prediction.matrix,sam,pm2, title ="" ){
  library(fmsb)
  training_groups <- training_groups[order(training_groups)]
  colours <- rainbow(length(training_groups))
  
  base <- t(cbind(rep(1,ncol(prediction.matrix)), rep(0 , ncol(prediction.matrix))))
  colnames(base) <- colnames(prediction.matrix)
  
  phecol <- colours[match(training_phenotypes[sam,"group"] , training_groups)]
  
  linecol <- "black"
  sample_data <-rbind(base, prediction.matrix[sam,] )
  
  if("ecoli" %in%  colnames(prediction.matrix)){
    sample_data <- sample_data[,ordered_training_groups]
  }else{
    sample_data <- sample_data[,ordered_training_groups_toplevel]
  }
  
  if(!is.null(pm2)){
    linecol <- c("grey","black")
    phecol <- c(gsub("FF$","99",phecol) , phecol  )
    sample_data <- rbind(sample_data, pm2[sam , ])
  }

  radarchart( sample_data  , axistype=1 , 
              title = title ,
              pcol=linecol , pfcol=phecol , plwd=3 , 
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1), cglwd=0.8,
              vlcex=0.8 
  )
}



plot_radar_multi <- function(prediction.matrix,sam,pm2, title ="" ){
  library(fmsb)
  training_groups <- training_groups[order(training_groups)]
  colours <- rainbow(length(training_groups))
  
  base <- t(cbind(rep(1,ncol(prediction.matrix)), rep(0 , ncol(prediction.matrix))))
  colnames(base) <- colnames(prediction.matrix)
  
  phecol <- colours[match(phenotypes[sam,"group"] , training_groups)]
  linecol <- "black"
  sample_data <-base
  if("ecoli" %in%  colnames(prediction.matrix)){
    sample_data <- sample_data[,ordered_training_groups]
  }else{
    sample_data <- sample_data[,ordered_training_groups_toplevel]
  }
  
  if(!is.null(pm2)){
    linecol <- c(rep("black",length(sam)),rep("black",length(sam)))
    phecol <- c(gsub("FF$","48",phecol) , gsub("FF$","48",phecol)   )
    sample_data <- rbind(sample_data, pm2[sam , ])
  }
  
  radarchart( sample_data  , axistype=1 , 
              title = title ,
              pcol=linecol , pfcol=phecol , plwd=3 , 
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1), cglwd=0.8,
              vlcex=0.8 
  )
}




plot_radar.2 <- function(prediction.matrix,sam, title ="", level = "low", plot_dec = T ){
  library(fmsb)
  if(level == "low"){
    groups <- c(bacteria,viral,inflammatory,"TB","KD","malaria")  
  }else{
    groups <- c("bacterial","viral","inflammatory","TB","KD","malaria")
  }
  colour_table[groups,"col"]
  pm <- prediction.matrix[[level]]
  pm2 <- prediction.matrix[[paste0(level,".dec")]]
  base <- t(cbind(rep(1,ncol(pm)), rep(0 , ncol(pm))))
  colnames(base) <- colnames(pm)
  phecol <- colour_table[phenotypes[sam,"group"],"col"]
  linecol <- "black"
  sample_data <-rbind(base, pm[sam,] )
  sample_data <- sample_data[,groups]
  if(plot_dec){
    linecol <- c("grey","black")
    phecol <- c(gsub("FF$","99",phecol) , phecol  )
    sample_data <- rbind(sample_data, pm2[sam , ])
    sample_data <- sample_data[c(1,2,4,3),]
  }
  colnames(sample_data) <- colour_table[colnames(sample_data) , "name"]
  radarchart( sample_data  , axistype=1 , 
              title = title ,
              pcol=linecol , pfcol=phecol , plwd=3 , 
              cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
              vlcex=0.8)
}

plot_parallel <- function(prediction.matrix,sam_list, main = "", plot_dec = T , multiplot = F,rnaseq = F,axis_buffer = 0.5){
  groups.low <- c(bacteria,"TB",viral,"malaria",inflammatory,"KD")  
  groups <- c("bacterial","TB","viral","malaria","inflammatory","KD")
  if("JIA" %in% colnames(prediction.matrix$top)){
    groups[5] <- "JIA"
    rnaseq <- T
  } 
  groups.low <- groups.low[groups.low %in% colnames(prediction.matrix$low)]  
  groups <- groups[groups %in% colnames(prediction.matrix$top)]  
  phecol <- colour_table[phenotypes[sam_list,"group"],"col"]
  ng <- length(groups.low) + length(groups)
  if(!multiplot){
    pretty_plot_area(cols = c("white","white"),text_col = "grey10",
                     ytick = c(0,length(sam_list),1),x_lab = "Disease",y_lab = "Sample predicted pribability",
                     xtick = c(1,ng,1), show_x = F, show_x_tick = T,show_y = F, show_y_tick = F,xbuffer = 0.4)
  }else{
    pretty_plot_area(cols = c("white","white"),text_col = "grey10",
                     ytick = c(0,length(sam_list),1),x_lab = "",y_lab = "",
                     xtick = c(1,ng,1), show_x = F, show_x_tick = T,show_y = F, show_y_tick = F,
                     xbuffer = 0.4,margins = c(4,1,5,0.2))
  }
  title(main = main,line = 3.5)

  # group highlight  
  if(rnaseq){
    rect(xleft = 0.5,ybottom = 0,xright = 5.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["bacterial","col"],0.3),border = NA)
    rect(xleft = 5.5+1,ybottom = 0,xright = 9.5+1,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["viral","col"],0.3),border = NA)
    rect(xleft = 9.5+2,ybottom = 0,xright = 10.5+2,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["inflammatory","col"],0.3),border = NA)
    rect(xleft = 13.5,ybottom = 0,xright = 14.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["bacterial","col"],0.3),border = NA)
    rect(xleft = 15.5,ybottom = 0,xright = 16.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["viral","col"],0.3),border = NA)
    rect(xleft = 17.5,ybottom = 0,xright = 18.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["inflammatory","col"],0.3),border = NA)
    abline(v = 13.5)
  }else{
    rect(xleft = 0.5,ybottom = 0,xright = 6.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["bacterial","col"],0.3),border = NA)
    rect(xleft = 6.5+1,ybottom = 0,xright = 12.5+1,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["viral","col"],0.3),border = NA)
    rect(xleft = 12.5+2,ybottom = 0,xright = 15.5+2,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["inflammatory","col"],0.3),border = NA)
    rect(xleft = 18.5,ybottom = 0,xright = 19.5,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["bacterial","col"],0.3),border = NA)
    rect(xleft = 19.5+1,ybottom = 0,xright = 20.5+1,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["viral","col"],0.3),border = NA)
    rect(xleft = 20.5+2,ybottom = 0,xright = 21.5+2,ytop = length(sam_list) + 1,lwd = 0,col = add.alpha(colour_table["inflammatory","col"],0.3),border = NA)
    abline(v = 18.5)
  }    
  # disease labels
  n_dis <- which(!is.na(match(sam_list, "N")))
  if(length(n_dis) > 0){
    n_dis  <- c(0,n_dis,length(sam_list)+1)
    for( i in 2:length(n_dis)){
      if(rnaseq){
        dis <- rnaseq_validation_pheno[sam_list[n_dis[i]-1]]  
      }else{
        dis <- phenotypes[sam_list[n_dis[i]-1],"group"]  
      }
      j <- match(dis,groups.low)
      rect(ybottom = n_dis[i-1],ytop = n_dis[i]-1, xleft = j-0.5 ,xright = j+0.5, col = "#66A61E", border = NA)
    }
  }else{
    if(rnaseq){
      dis <- unique(rnaseq_validation_pheno[sam_list]  )
    }else{
      dis <- unique(phenotypes[sam_list,"group"])
    }
    j <- match(dis,groups.low)
    rect(ybottom = 0,ytop = length(sam_list), xleft = j-0.5 ,xright = j+0.5, col = "#66A61E", border = NA)
  }
  
  if(rnaseq){ # FYI this will not let multiple broad disease groups be plotted together
    dis <- unique(rnaseq_broad_validation_pheno[sam_list]  )
  }else{
    dis <- unique(phenotypes[sam_list,"group"])
    if(any(dis %in% bacteria)){
      dis <- "bacterial"
    }else if(any(dis %in% viral)){
      dis <- "viral"
    }else if(any(dis %in% inflammatory)){
      dis <- "inflammatory"
    }
  }
  j <- match(dis,groups) + length(groups.low)
  rect(ybottom = 0,ytop = length(sam_list), xleft = j-0.5 ,xright = j+0.5, col = "#66A61E", border = NA)
  # the ridgeline
  for(i in 1:length(sam_list)){
    if(sam_list[i] != "N"){
      polygon(c(0,1:ng,ng + 1), c(i-1,i - 1 + as.numeric(c(prediction.matrix$low.dec[sam_list[i],groups.low],prediction.matrix$top.dec[sam_list[i],groups])),i-1),col = "grey",border = "grey", lwd = 0.5)
      polygon(c(0,1:ng,ng + 1), c(i-1,i - 1 + as.numeric(c(prediction.matrix$low[sam_list[i],groups.low],prediction.matrix$top[sam_list[i],groups])),i-1),  col = "black",border = "black", lwd = 0.5)      
    }
  }
  
  text(1:ng, par("usr")[3]-axis_buffer, 
       srt = 60, adj = 1, xpd = TRUE,
       labels = c(colour_table[groups.low,"name"],colour_table[groups,"name"]), cex = 0.8)
  text(1:ng, par("usr")[4]+axis_buffer, 
       srt = 300, adj = 1, xpd = TRUE,
       labels = c(colour_table[groups.low,"name"],colour_table[groups,"name"]), cex = 0.8)
}





plot_all_pr <- function(pr_obj,solidclass =c(),plotclasses,indicator,lwd_sf = 1 ){
  prcols <- add.alpha(rainbow(length(plotclasses)),0.6)
  plot(pr_obj[[1]]$curve[,1] , pr_obj[[1]]$curve[,2],type = "l",xlab ="Recall",ylab ="Precision",xlim = c(0,1),ylim= c(0,1),col ="white")
  leg <- c()
  auc_tab <- data.frame(matrix(nrow = length(plotclasses),ncol = 3))
  colnames(auc_tab) <- c("disease","AUC","n_samples")
  for(i in 1:length(pr_obj)){
    classname <- names(pr_obj)[i]
    if(classname %in% plotclasses){
      points(pr_obj[[i]]$curve[,1] , pr_obj[[i]]$curve[,2]
             ,type = "l",col = prcols[match(classname ,plotclasses)],lty =1,
             lwd = lwd_sf * colSums(indicator)[classname] * (1/max(colSums(indicator))))
      auc_tab[i, ] <- c(classname,round(digits = 3,pr_obj[[i]]$auc.integral),colSums(indicator)[classname])
      leg <- c(leg, paste(classname,round(digits = 3,pr_obj[[i]]$auc.integral)))
    }
  }
  auc_tab[!is.na(auc_tab[,1]),]
}


plot_all_roc <- function(roc_obj,plotclasses,main ="",indicator, lwd_sf = 1,split = F){
  par(mar=c(3.5,3.5,2,2))
  if(split){
    par(mfrow = c(ceiling(length(plotclasses)/4),4))
  }else{
    plot(roc_obj[[1]], main =main,col="white")
  }
  leg <- c()
  auc_tab <- data.frame(matrix(nrow = length(plotclasses),ncol = 3))
  colnames(auc_tab) <- c("disease","AUC","n_samples")
  for(i in 1:length(roc_obj)){
    classname <- names(roc_obj)[i]
    if(classname %in% plotclasses){
      if(split){
        plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(1,0),axes =F,xlab = "",ylab = "")
        title(main =  colour_table[classname,"name"],line =1)
        title(xlab = "Specificity",ylab = "Sensitivity",line =lin)
        axis(1,pos = 0)
        axis(2,las =2,pos = 1)
        points(c(1,0),c(1,1),type = "l")
        points(c(1,0),c(0,1),type = "l",col ="grey")
        par(xpd = T)
        points(roc_obj[[i]]$specificities,
               roc_obj[[i]]$sensitivities,
               col = colour_table[classname , "col"],
               type = "l",lwd =2)
        
      }else{
        plot(roc_obj[[i]],add =T,col = colour_table[classname,"col"],
             lwd = lwd_sf * colSums(indicator)[classname] * (1/max(colSums(indicator))))
      }
      auc_tab[i,] <- c(classname, round(digits = 3,pROC::auc(roc_obj[[i]])) , colSums(indicator)[classname] )
      leg <- c(leg ,  paste(classname, round(digits = 3,pROC::auc(roc_obj[[i]]))) )
    }
  }
  auc_tab[!is.na(auc_tab[,1]),]
  par(mfrow = c(1,1))
  par(xpd = F)
}



plot_R <- function(rocs_obj ,group="all",toplevel = F,main = "", rnaseq = F,plot_legend = T,only_legned = F){
  if(only_legned)plot(0,0,cex =0 ,ylim =c(-10,1),xlim = c(0,1.5))
  pathogens <- colnames(test_phenotypes_indicator)[!colnames(test_phenotypes_indicator)%in%c("bacterial","viral","inflammatory")]
  pathogens <- c(bacteria,viral,inflammatory, "malaria","TB","KD" )
  if(toplevel == T){
    pathogens <- c("bacterial","viral","inflammatory", "malaria","TB","KD"   )
  }
  colours <- rainbow(length(pathogens))
  names(colours) <- pathogens
  addval <- F
  newcolours <- c()
  labels_auc <- c()
  for(path in pathogens){
    cases <- rownames(test_phenotypes_indicator)[test_phenotypes_indicator[,path]==1]
    r <- rocs_obj[[path]]
    if(group==path | group=="all"){
      colour <- colours[path]# rainbow
      if(path %in% bacteria){
        fc <- colorRampPalette(c("pink", "darkred"))
        colour <- fc(6)[match(path,bacteria)]
      }else if(path %in% viral){
        fc <- colorRampPalette(c("cyan", "darkblue"))
        colour <- fc(6)[match(path,viral)]
      }else if(path %in% inflammatory){
        fc <- colorRampPalette(c("orange", "yellow"))
        colour <- fc(3)[match(path,inflammatory)]
      }else{
        colour <- col2hex(c("black","purple","grey"))[match(path,c("malaria","TB","KD"))]
      }
      newcolours <- c(newcolours ,colour)  
      sf <- length(cases)/52
      if(!only_legned) plot(r,add=addval , col = colour,main=main,lwd = 2,cex.lab =2) #,lwd =sf*10) 
      youd <- coords(r, "b", best.method="youden")
      if(path %in%c("gas","gbs")){
        path.u <- toupper(path)
      }else{
        path.u <- capitalize(path)
        if(path=="ecoli"){
          path.u <- "E.coli"  
        }else if(path=="pneumo"){
          path.u <- "Pneumococcal"  
        }else if(path=="staph"){
          path.u <- "Staphylococcal"  
        }else if(path=="adeno"){
          path.u <- "Adenovirus"  
        }else if(path=="flu"){
          path.u <- "Influenza"  
        }
        
      }
      pathperf <- paste(path.u ," AUC = ",round(pROC::auc(r)*100,digits = 1),"","n(cases) = ",length(cases))
      labels_auc <- c(labels_auc ,pathperf)
      ypos <- 1 - match(path,pathogens)*0.05
      if(group==path){
        for(vsgroup in c("bacterial","viral","inflammatory")){
          addval <- T
          r <- rocs_obj[[paste(path,"-vs-",vsgroup,sep = "")]]
          if(!only_legned) plot(r,add=addval , col = colours[path],main=main,lty = 3)
          labels_auc <- c(labels_auc ,paste(paste(path,"-vs-",vsgroup,sep = "") ," AUC = ",
                                            round(pROC::auc(r)*100,digits = 4),"","n(cases) = ",length(cases)))
          if(plot_legend){
            if(!only_legned) text(0.2,0.95,cex =0.5,paste(path," vs ",vsgroup," (using ",path," response) ",
                                                          sep = "","AUC = ",round(pROC::auc(r)*100,digits = 4),"\n","n(cases) = ",length(cases)))
          }
          
        }
      }
      addval <- T
    }

  }
  if(only_legned)  legend(0.2,0.8,legend = labels_auc,col = newcolours,
                          fill = colours,pt.cex = 0.6,cex = 0.9,y.intersp = 0.2)   
  if(group=="all" ){
    if(plot_legend ){
      if(!only_legned ) legend(0.75,0.9,legend = labels_auc,col = newcolours,fill = newcolours,bty = "n",
                               pt.cex = 0.4,cex = 1.3,x.intersp=0.5,y.intersp = 0.5)   
    }
  }
}


plot_pr_micro <- function(perf.list,mod = "classes",add = F){
  sf <- 2 # scale pointsize
  if(!add){
    plot(0,0,cex = 0 , ylim = c(0,1),xlim = c(0,1), xlab = "Recall",ylab = "Precision")
  }
  if(mod == "classes"){
    for(i  in 1:18){
      class <- colnames(perf.list$recall.class)[i]
      points(perf.list$recall.class[,class] , perf.list$precision.class[,class], cex = sf * perf.list$f1.class[,class] ,
             pch = 16,type = "b",col = rainbow(18)[i])  
    }
  }else if(mod == "micro"){
    points(perf.list$recall.micro , perf.list$precision.micro, cex = sf * perf.list$f1.micro , # do better in the micro averaged because do badly on rare classes
           pch = 16,type = "b",col = "black")
  }else if(mod == "ml.hierarchy"){
    points(perf.list$multilabel_heierarchical.recall,perf.list$multilabel_heierarchical.precision, 
           cex = sf * perf.list$multilabel_heierarchical.f,pch = 16,type = "b",col = "black")
  }else if(mod == "ml"){
    points(perf.list$multilabel.recall,perf.list$multilabel.precision,
           cex = sf * perf.list$multilabel.f,pch = 16,type = "b",col = "black")
  }
}



plot_PR <- function(curve_obj ,group="all",toplevel = F,main = "", rnaseq = F,legendpos = c(0,0.8)){
  
  pathogens <- colnames(test_phenotypes_indicator)[!colnames(test_phenotypes_indicator)%in%c("bacterial","viral","inflammatory")]
  if(toplevel == T){
    pathogens <- c("bacterial","viral","inflammatory", "malaria","TB","KD"   )
  }else if(rnaseq == T){
    pathogens <- c("bacterial","viral","inflammatory","TB","KD"   )
  }
  colours <- rainbow(length(pathogens))
  names(colours) <- pathogens
  labels_auc <- c()
  auclist <- data.frame()
  plot(0,0,cex = 0 ,ylim = c(0,1),xlim = c(0,1),xlab = "Recall",ylab = "Precision",main = main)
  for(path in pathogens){
    cases <- rownames(test_phenotypes_indicator)[test_phenotypes_indicator[,path]==1]
    pr <- curve_obj[[path]]
    auclist[path,1]  <- pr$auc.integral
    if(group==path | group=="all"){
      points(pr$curve[,1],pr$curve[,2] , col = colours[path],type = "l") 
      labels_auc <- c(labels_auc ,paste(path ," AUC = ",round(pr$auc.integral*100,digits = 4),"",
                                        "n(cases) = ",length(cases)))
      if(group==path){
        text(0.2,0.8,paste("AUC = ",round(pr$auc.integral*100,digits = 4),"\n","n(cases) = ",length(cases)))
        for(vsgroup in c("bacterial","viral","inflammatory")){
          addval <- T
          pr <- curve_obj[[paste(path,"-vs-",vsgroup,sep = "")]]
          points(pr$curve[,1],pr$curve[,2] , col = colours[path],lty = 3,type = "l") 
          ypos <- 0.8 - (0.1 * match(vsgroup , c("bacterial","viral","inflammatory")))
          text(0.2,ypos,paste(path," vs ",vsgroup," (using ",path," response) ",
                              sep = "","AUC = ",round(pr$auc.integral*100,digits = 4),"\n","n(cases) = ",length(cases)))
        }
      }
      addval <- T
    }
  }
  if(group=="all"){
    legend(legendpos[1],legendpos[2],legend = labels_auc,col = colours,fill = colours,cex = 0.8)    
  }
  return(auclist)
}




plot_PR_vs_nsam <- function(curve_obj ,group="all",toplevel = F,main = "", rnaseq = F,legendpos = c(0,0.8),sf =1){
  pathogens <- colnames(test_phenotypes_indicator)[!colnames(test_phenotypes_indicator)%in%c("bacterial","viral","inflammatory")]
  if(toplevel == T){
    pathogens <- c("bacterial","viral","inflammatory", "malaria","TB","KD"   )
  }else if(rnaseq == T){
    pathogens <- c("bacterial","viral","inflammatory","TB","KD"   )
  }
  par(mar = c(5,5,1,1))
  colours <- rainbow(length(pathogens))
  names(colours) <- pathogens
  labels_auc <- c()
  auclist <- data.frame()
  plot(0,0,cex = 0 ,ylim = c(0,1.1),xlim = c(0,210),
       xlab = "number of samples",ylab = "AUPRC"
       ,main = main,cex.axis = 1,las = 2,cex.lab = 1)
  for(path in pathogens){
    cases <- rownames(phenotypes)[phenotypes[,"group"] == path]
    pr <- curve_obj[[path]]
    auclist[path,1]  <- pr$auc.integral
    points(length(cases), pr$auc.integral ,pch = 16, col = colours[path],cex = sf)
    if(path %in% c("TB","HSP")){
      text(length(cases), pr$auc.integral,pos = 1, labels = path,cex = sf)
    }else if(path == "malaria"){
      text(length(cases), pr$auc.integral,pos =2 , labels = path,cex = sf)
    }else if(path %in% c("HHV6","adeno","staph","SLE","JIA")){
      text(length(cases), pr$auc.integral,pos =3, labels = path,cex = sf)
    }else{
      text(length(cases), pr$auc.integral,pos =4, labels = path,cex = sf)
    }
  }
  
}




plot_curve_loc <- function(curve , type = "PR", sf = 1,more_pr_loc = NULL,
                           plotclasses =NULL,linethickness =T,dstype ="micro"){
  if(type =="PR"){
    pretty_plot_area(cols = c("grey90","white"),text_col = "grey10",
                     ytick = c(0,1,0.2),x_lab = "Recall",y_lab = "Precision",
                     xtick = c(0,1,0.2))
  }else if(type == "ROC"){
    pretty_plot_area(cols = c("white","grey90"),text_col = "grey30",
                     ytick = c(0,1,0.2),x_lab = "Specificity",y_lab = "Sensitivity",
                     xtick = c(1,0,-0.2),xbuffer = 0.05,ybuffer = 0.05)
  }
  if(dstype == "micro"){
    indicator <- test_phenotypes_indicator
  }else if(dstype == "rnaseq"){
    indicator <- rnaseq_validation_pheno.ind.loc[test_rnaseq,]
  }
  class_abundance <- colSums(indicator)
  class_abundance <- log(sf * 10 * class_abundance / max(class_abundance))  
  
  if(!linethickness){
    class_abundance[] <- sf
  }  
  if(is.null(plotclasses)){
    grps_loc <- names(curve)[grep(invert = T ,"vs",names(curve))]  
  }else{
    grps_loc <- plotclasses
  }
  auc_mat <- data.frame(matrix(nrow = length(grps_loc),ncol = 5))
  for(i in grps_loc){
    if(type =="PR"){
      auc_mat[match(i , grps_loc),] <- c(i,
                                         round( curve[[i]]$auc.integral,digits = 4),
                                         as.numeric(colSums(indicator)[i]))
      points(curve[[i]]$curve[,1],curve[[i]]$curve[,2] , lwd = (class_abundance[i]),
             col = add.alpha(colour_table[i,"col"],alpha = 0.9),lty = as.numeric(colour_table[i,"line"]),type = "l") 
    }else if(type == "ROC"){
      auc_mat[match(i , grps_loc),] <- c(i,
                                         round( curve[[i]]$auc,digits = 4) ,
                                         as.numeric(colSums(indicator)[i]),
                                         round( curve[[i]]$ci,digits = 4)[c(1,3)]) # lower and upper ci
      points(curve[[i]]$specificities,curve[[i]]$sensitivities , lwd = class_abundance[i],
             col = add.alpha(colour_table[i,"col"],alpha = 0.9),lty =colour_table[i,"line"],type = "l") 
    }
    
    
    # AUCS by different methods ;; but  auc(r) == DescTools::AUC(r$specificities,r$sensitivities) to ~7dp
  }
  if(!is.null(more_pr_loc)){
    if(type =="PR"){
      points(more_pr_loc$recall.micro , more_pr_loc$precision.micro,type = "l",lty = 1,lwd =2)
      points(more_pr_loc$recall.macro , more_pr_loc$precision.macro,type = "l",lty = 3,lwd =2)
      
      aucval <-  trapz(x = more_pr_loc$recall.micro ,y = more_pr_loc$precision.micro )
      auc_mat[nrow(auc_mat)+1,] <- c("micro average",
                                     round( aucval,digits = 4),
                                     dim(indicator)[[1]],NA,NA)
      
      aucval <-  trapz(x = more_pr_loc$recall.macro ,y = more_pr_loc$precision.macro)
      auc_mat[nrow(auc_mat)+1,] <- c("macro average",
                                     round( aucval,digits = 4) ,
                                     dim(indicator)[[1]],NA,NA)
    }else if(type == "ROC"){
      points(more_pr_loc$specificity.micro, more_pr_loc$sensitivity.micro,type = "l",lty = 1,lwd =2)
      points(more_pr_loc$specificity.macro, more_pr_loc$sensitivity.macro,type = "l",lty = 3,lwd =2)
      
      aucval <-  trapz(x = more_pr_loc$specificity.micro ,y = more_pr_loc$sensitivity.micro)
      auc_mat[nrow(auc_mat)+1,] <- c("micro average",
                                     round( aucval,digits = 4) ,
                                     dim(indicator)[[1]],NA,NA)
      
      aucval <-  trapz(x = more_pr_loc$specificity.macro ,y = more_pr_loc$sensitivity.macro)
      auc_mat[nrow(auc_mat)+1,] <- c("macro average",
                                     round( aucval,digits = 4) ,
                                     dim(indicator)[[1]],NA,NA)
      
    }
  }
  colnames(auc_mat) <- c("disease" , "AUC","#cases","CI.l","CI.h")
  auc_mat[,2] <- as.numeric(auc_mat[,2])
  auc_mat[,3] <- as.numeric(auc_mat[,3])
  auc_mat <- auc_mat[order(auc_mat[,3],decreasing = T),]
  return(auc_mat)
}


confusion_heatmap <- function(confusion.loc,main ="",col_r = c(0,5),
                              scale = "proportion",level="low",square_lab="fraction",
                              class_abundance_vector =NULL,margins = c(10, 10, 2, 1),textsize = 1){

  pretty_plot_area(ytick = c(0,nrow(confusion.loc),1),
                   xtick = c(0,ncol(confusion.loc),1),
                   show_x = F,show_y = F,margins = margins)
  
  if(scale =="raw"){
    col_range_wrong <- colorRamp2(col_r,colors = c("grey90","red"))
    col_range_true <- colorRamp2(col_r,colors = c("grey90","green"))
    col_range_close <- colorRamp2(col_r,colors = c("grey90","yellow"))  
  }else{
    col_r <- c(0,1)
    col_range_wrong <- colorRamp2(col_r,colors = c("grey90","red"))
    col_range_true <- colorRamp2(col_r,colors = c("grey90","green"))
    col_range_close <- colorRamp2(col_r,colors = c("grey90","yellow"))  
    
  }
  
  axis(side = 1,at = 1:ncol(confusion.loc) - 0.5,tick = F,cex.axis = 1.7,labels = colour_table[colnames(confusion.loc),"name"],las =2)
  axis(side = 2,at = 1:nrow(confusion.loc) - 0.5,tick = F,cex.axis = 1.7,labels = colour_table[colnames(confusion.loc),"name"],las =2)
  
  if(level =="top"){
    hier <- c("bacterial","viral","inflammatory","KD","TB","malaria")
    names(hier) <- c("bacterial","viral","inflammatory","KD","TB","malaria")
  }else if(level == "top_rna"){
    hier <- c("bacterial","viral","inflammatory","KD","TB","malaria")
    names(hier) <- c("bacterial","viral","JIA","KD","TB","malaria")
  }else if(level == "low"){
    hier <- c(rep("bacterial",length(bacteria)+1),rep("viral",length(viral)),rep("inflammatory",length(inflammatory)+1),"malaria")
    names(hier) <- c(bacteria,"TB",viral,inflammatory,"KD","malaria")
  }else{
    hier <- c(rep("bacterial",length(bacteria)+2),rep("viral",length(viral)+1),
              rep("inflammatory",length(inflammatory)+2),"malaria")
    names(hier) <- c(bacteria,"TB","bacterial",viral,"viral",inflammatory,"KD","inflammatory","malaria")
  }
  
  hier["RF"] <- "RF"
  
  if(is.null( class_abundance_vector)){
    class_abundance_vector <- rowSums(confusion.loc)
  }else{
    class_abundance_vector <- class_abundance_vector[colnames(confusion.loc)]
  }
  
  
  for(i in 1:nrow(confusion.loc)){
    for(j in 1:ncol(confusion.loc)){
      if(scale =="raw"){
        if(i==j){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_true(confusion.loc[i,j]),border =F)
        }else if(hier[colnames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_close(confusion.loc[i,j]),border =F)
        }else{
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_wrong(confusion.loc[i,j]),border =F) 
        }
      }else if(scale =="proportion"){
        if(i==j){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_true(confusion.loc[i,j] /  class_abundance_vector[i]),border =F)
        }else if(hier[colnames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_close(confusion.loc[i,j] / class_abundance_vector[i]),border =F)
        }else{
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_wrong(confusion.loc[i,j] /  class_abundance_vector[i]),border =F) 
        }
      }else if(!is.null(class_abundance_vector)){
        if(i==j){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_true(confusion.loc[i,j] / class_abundance_vector[i]),border =F)
        }else if(hier[colnames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_close(confusion.loc[i,j] / class_abundance_vector[i]),border =F)
        }else{
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_wrong(confusion.loc[i,j] / class_abundance_vector[i]),border =F) 
        }
      }
      
      if(confusion.loc[i,j]!=0){
        if(square_lab =="number"){
          text(j-0.5,i-0.5,round(digits = 2,confusion.loc[i,j]),cex = textsize)  
        }else if(square_lab == "fraction"){
          text(j-0.5,i-0.5,paste(confusion.loc[i,j] ,"/", class_abundance_vector[i],sep =" " ),cex = textsize)  
        }else if(square_lab == "percent"){
          text(j-0.5,i-0.5,paste(round(100*confusion.loc[i,j]/ class_abundance_vector[i], digits = 1),"%") ,cex = textsize)  
        }else if(square_lab == "both"){
          text(j-0.5,i-0.5,paste(confusion.loc[i,j] ,"/", class_abundance_vector[i] ,"\n"
                                 ,round(100*confusion.loc[i,j]/ class_abundance_vector[i] , digits = 1),"%") ,cex = textsize)  
        }
      }
    }
  }
  abline(h=0:20,col="white",lwd =2)
  abline(v=0:20,col="white",lwd =2)
  title(xlab = "Predicted",ylab = "True",line = 8.8,cex.lab =1.7)
}



confusion_heatmap_sizemod <- function(confusion.loc,main ="",col_r = c(0,5),
                                      scale = "proportion",level="low",square_lab="fraction",
                                      class_abundance_vector =NULL,margins = c(10, 10, 2, 1),textsize = 1,sf = 1,addrects =NULL,sf.ax = 1){
  
  confusion.loc <- apply(confusion.loc,2,rev) 
  pretty_plot_area(ytick = c(0,nrow(confusion.loc),1),cols = c("white","white"),
                   xtick = c(0,ncol(confusion.loc),1),
                   show_x = F,show_y = F,margins = margins)
  abline(v = seq(0.5,20.5, by = 1), col = "lightgray", lwd = 2)
  abline(h = seq(0.5,20.5, by = 1), col = "lightgray", lwd = 2)
  if(scale =="raw"){
    col_range_wrong <- colorRamp2(col_r,colors = c("grey90","red"))
    col_range_true <- colorRamp2(col_r,colors = c("grey90","green"))
    col_range_close <- colorRamp2(col_r,colors = c("grey90","yellow"))  
  }else{
    col_r <- c(0,1)
    col_range_wrong <- colorRamp2(col_r,colors = c("#E7298A"  ,"#E7298A"  ))
    col_range_true <- colorRamp2(col_r,colors = c("#66A61E" ,"#66A61E" ))
    col_range_close <- colorRamp2(col_r,colors = c("#A6761D" ,"#A6761D" ))  
  }
  if(!is.null(addrects)){
    addrects()
  }
  axis(side = 2,at = 1:ncol(confusion.loc) - 0.5,tick = F,cex.axis = sf.ax,
       labels = colour_table[rownames(confusion.loc),"name"] ,las =2)
  axis(side = 1,at = 1:nrow(confusion.loc) - 0.5,tick = F,cex.axis = sf.ax,
       labels = colour_table[colnames(confusion.loc),"name"],las =2)
  
  if(level =="top"){
    hier <- c("bacterial","viral","inflammatory","KD","TB","malaria")
    names(hier) <- c("bacterial","viral","inflammatory","KD","TB","malaria")
  }else if(level == "top_rna"){
    hier <- c("bacterial","viral","inflammatory","KD","TB","malaria")
    names(hier) <- c("bacterial","viral","JIA","KD","TB","malaria")
  }else if(level == "low"){
    hier <- c(rep("bacterial",length(bacteria)+1),rep("viral",length(viral)),rep("inflammatory",length(inflammatory)+1),"malaria")
    names(hier) <- c(bacteria,"TB",viral,inflammatory,"KD","malaria")
    
  }else{
    hier <- c(rep("bacterial",length(bacteria)+2),rep("viral",length(viral)+1),
              rep("inflammatory",length(inflammatory)+2),"malaria")
    names(hier) <- c(bacteria,"TB","bacterial",viral,"viral",inflammatory,"KD","inflammatory","malaria")
  }
  hier["RF"] <- "RF"
  if(is.null( class_abundance_vector)){
    class_abundance_vector <- rowSums(confusion.loc)
  }else{
    class_abundance_vector <- class_abundance_vector[colnames(confusion.loc)]
  }
  
  for(i in 1:nrow(confusion.loc)){
    for(j in 1:ncol(confusion.loc)){
      if(scale =="raw"){
        if(rownames(confusion.loc)[i]==colnames(confusion.loc)[j]){
          points(j-0.5,i-0.5,col = col_range_true(confusion.loc[i,j]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }else if(hier[rownames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          points(j-0.5,i-0.5,col = col_range_close(confusion.loc[i,j]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }else{
          points(j-0.5,i-0.5,col = col_range_wrong(confusion.loc[i,j]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }
      }else if(scale =="proportion"){
        if(rownames(confusion.loc)[i]==colnames(confusion.loc)[j]){
          points(j-0.5,i-0.5,col = col_range_true(confusion.loc[i,j] /  class_abundance_vector[i]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }else if(hier[rownames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          points(j-0.5,i-0.5,col = col_range_close(confusion.loc[i,j] /  class_abundance_vector[i]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }else{
          points(j-0.5,i-0.5,col = col_range_wrong(confusion.loc[i,j] /  class_abundance_vector[i]) , cex =  sf * sqrt(confusion.loc[i,j]/ pi),pch =16)
        }
      }else if(!is.null(class_abundance_vector)){
        if(rownames(confusion.loc)[i]==colnames(confusion.loc)[j]){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_true(confusion.loc[i,j] / class_abundance_vector[i]),border =F)
        }else if(hier[colnames(confusion.loc)[i]]==hier[colnames(confusion.loc)[j]] ){
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_close(confusion.loc[i,j] / class_abundance_vector[i]),border =F)
        }else{
          rect(xleft = j-1,xright = j,ybottom = i-1,ytop = i,col = col_range_wrong(confusion.loc[i,j] / class_abundance_vector[i]),border =F) 
        }
      }
      
      if(confusion.loc[i,j]!=0){
        if(square_lab =="number"){
          if(confusion.loc[i,j]<5){
            text(j-0.3,i-0.3,round(digits = 2,confusion.loc[i,j]),cex = textsize)  
          }else{
            text(j-0.5,i-0.5,round(digits = 2,confusion.loc[i,j]),cex = textsize)  
          }
        }else if(square_lab == "fraction"){
          text(j-0.5,i-0.5,paste(confusion.loc[i,j] ,"/", class_abundance_vector[i],sep =" " ),cex = textsize)  
        }else if(square_lab == "percent"){
          text(j-0.5,i-0.5,paste(round(100*confusion.loc[i,j]/ class_abundance_vector[i], digits = 1),"%") ,cex = textsize)  
        }else if(square_lab == "both"){
          text(j-0.5,i-0.5,paste(confusion.loc[i,j] ,"/", class_abundance_vector[i] ,"\n"
                                 ,round(100*confusion.loc[i,j]/ class_abundance_vector[i] , digits = 1),"%") ,cex = textsize)  
        }
      }
    }
  }
}




plot_prediction_circs <- function(pred.resp,trueclasses = NULL ,main = "",scale_fact = 2.5,
                                  no_true = F){
  par(mar = c(10,10,1,1))
  pred.resp <- reorder_predmat(pred.resp)
  # reorder to group paths 
  rs <- rownames(pred.resp)[order(phenotypes[rownames(pred.resp),"group"])]
  pred.resp <- pred.resp[rs,]
  
  grps <- colnames(pred.resp)[order( colnames(pred.resp))]
  sams <- rownames(pred.resp)
  cols <- rainbow(length(grps))  
  layout(matrix(1, nrow=1)) 
  plot(-10,-10,ylim=c(0,length(grps)),xlim=c(0,length(sams)),yaxt = 'n',xaxt = 'n',xlab="",ylab = "",main=main)
  axis(2, at=1:length(grps), labels=grps,las = 2)

  for(i in 1:length(sams)){
    sam <- sams[i]
    for(j in 1:length(grps)){
      if(!no_true){
        points(i,j,cex= sqrt(pred.resp[i,j]) * scale_fact,pch=16,col= cols[match(phenotypes[sam,"group"] , grps)])
      }else{
        points(i,j,cex= sqrt(pred.resp[i,j]) * scale_fact,pch=16,col= "red")  
      }
      if(pred.resp[i,j] == max(pred.resp[i,1:(length(grps))])){
        points(i,j,cex= sqrt(pred.resp[i,j]) * scale_fact,pch=8,lwd = 2,col="black")
        if(no_true){
          points(i,j,cex=sqrt(1) * scale_fact  ,pch=1,lwd = 2,col="orange")  
        }
      }
     if(!no_true){ 
        if(pred.resp[i,j] == max(pred.resp)){
          points(i,j,cex= pred.resp[i,j]*5,pch=8,lwd = 2,col="black")
        }
        if(grps[j] == phenotypes[sam,"group"]){
          points(i,j,cex= sqrt(pred.resp[i,j])*scale_fact,pch=1,lwd = 2,col="black")
        }
    }
    }
  }
  axis(1, at=1:length(sams), labels=sams,las = 2,cex = 0.75)
}



plot_prediction_circs_toplevel <- function(pred.resp,disease_set ,main = "",scale_fact = 2.5,no_true = F, 
                                           indicator = training_phenotypes_indicator){
  par(mar = c(10,10,1,1))
  pred.resp <- reorder_predmat(pred.resp)
  # reorder to group paths 
  rs <- rownames(pred.resp)[order(phenotypes[rownames(pred.resp),"group"])]
  pred.resp <- pred.resp[rs,]
  
  grps <- colnames(pred.resp)[order( colnames(pred.resp))]
  sams <- rownames(pred.resp)
  training_groups <- colnames(pred.low)
  cols <- rainbow(length(training_groups))  
  layout(matrix(1, nrow=1)) 
  plot(-10,-10,ylim=c(0,length(grps)+1),xlim=c(0,length(sams)),yaxt = 'n',xaxt = 'n',xlab="",ylab = "",main=main)
  axis(2, at=1:length(grps), labels=grps,las = 2)
  for(i in 1:length(sams)){
    sam <- sams[i]
    for(j in 1:length(grps)){
      points(i,j,cex= sqrt(pred.resp[i,j]) * scale_fact,pch=16,col= cols[match(phenotypes[sam,"group"] , training_groups)])
      text(i,j+0.5 , srt = 90 , labels = phenotypes[sam,"group"])
      if(pred.resp[i,j] == max(pred.resp[i,1:(length(grps))])){
        points(i,j,cex= sqrt(pred.resp[i,j]) * scale_fact,pch=8,lwd = 2,col="black")

      }
        if(pred.resp[i,j] == max(pred.resp[i,])){
          points(i,j,cex= pred.resp[i,j]*5,pch=8,lwd = 2,col="black")
        }
        if(indicator[sam , grps[j]]==1){
          points(i,j,cex= sqrt(pred.resp[i,j])*scale_fact,pch=1,lwd = 2,col="black")
        }
      }
    }
  axis(1, at=1:length(sams), labels=sams,las = 2,cex = 0.75)
}




plot_signature_beeswarm <- function(prediction_matrix,path , indicator = NULL,pcex = 1,toplevel = F,main = ""){
  ordered_training_groups <- c("ecoli","gas","gbs","meningococcal","pneumo","staph","adeno","enterovirus","flu","HHV6","rhino","RSV","HSP","JIA","SLE","KD","malaria" ,"TB" )
  ordered_training_groups_toplevel <- c("bacterial","viral","inflammatory","KD","malaria" ,"TB" )
  training_groups <- ordered_training_groups[order(ordered_training_groups)]
  if(toplevel){
    training_groups <- ordered_training_groups_toplevel
    ordered_training_groups <- ordered_training_groups_toplevel
  } 
  pltdat <- data.frame()
  colours <- rainbow(length(training_groups))
  colours <- gsub("FF$","80",colours)
  for(i in 1: length(training_groups)){
    sams <- rownames(indicator)[indicator[,training_groups[i]] == 1]
    pltdat[sams , 1] <- prediction_matrix[sams,path]
    pltdat[sams , 2] <- match(training_groups[i] , ordered_training_groups)
    pltdat[sams , 3] <- colours[i]
  }
  colnames(pltdat) <- c("p","pathogen","col")
  boxplot(p ~ pathogen , cols = col , data = pltdat , outline = FALSE ,xaxt = 'n',ylab = "p",ylim = c(0,1.1),main = main)
  beeswarm(p ~ pathogen , corral = "wrap",
           pwcol = pltdat[,"col"] , pch = 16 , data = pltdat ,cex = pcex, xaxt = 'n',add = T)
  axis(side = 1 , at = 1:(length(training_groups)) , labels = c(ordered_training_groups),las = 2)
  if(!toplevel){
    rect(xleft = 0.7,xright = 6.3 , ybottom = 1,ytop = 1.05,col = rgb(1, 0, 0,alpha = 0.5),border = NA)
    text(3,1.07 ,"bacterial")
    rect(xleft = 6.7,xright = 12.3 , ybottom = 1,ytop = 1.05,col = rgb(0, 1, 0,alpha = 0.5),border = NA)
    text(9,1.07 ,"viral")
    rect(xleft = 12.7,xright = 15.3 , ybottom = 1,ytop = 1.05,col = rgb(0, 0, 1,alpha = 0.5),border = NA)
    text(14,1.07 ,"inflammatory")
  }
}



plot_scree <- function(pca,main ="",tickint =10){
  par(mar = c(5,6,2,2))
  eigen <- get_eigenvalue(pca)[1:10,2]
  plot(eigen,cex =0,ylab = "Percentage of\n variance explained",xlab = "Dimensions",
       axes = F,type ="b",pch =1,xlim =c(0,11),ylim = c(0,max(eigen)),main = main,xaxs="r", yaxs="r")
  for(i in 1:10){
    rect(xleft = i-0.3,xright = i + 0.3,ybottom = 0,ytop = eigen[i],col = "lightgrey",border = F)
  }
  abline(h=seq(0,100,by = tickint) , col = "white",lwd =2)
  points(eigen,type ="b")
  axis(side = 1,1:10,1:10,tick = F)
  axis(side = 2,tick=T,lwd =1,lwd.ticks = 1,las = 2)
}



signature_set_heatmap <- function(signatures,sigmethod = "--", fold = "--",mat =F){
  overlaps <- data.frame(matrix(ncol = length(signatures),nrow = length(signatures)))
  colnames(overlaps) <- names(signatures)
  rownames(overlaps) <- names(signatures)
  overlaps[] <- 0
  comparisons <- combn(names(signatures),2)
  for(i in 1:ncol(comparisons)){
    path_a <- comparisons[1,i]
    path_b <- comparisons[2,i]
    shared <- length(intersect(signatures[[path_a]] , signatures[[path_b]]))
    overlaps[path_a,path_b] <- shared
    overlaps[path_b,path_a] <- shared
    overlaps[path_a,path_a] <- length(signatures[[path_a]])
    overlaps[path_b,path_b] <- length(signatures[[path_b]])
  } 
  if(mat == F){
    superheat(overlaps,
              X.text = round(as.matrix(overlaps), 1),
              X.text.size = 4)
  }else{
    return(overlaps)
  }
}



plot_errchange <- function(opt,m1,m2,pn,fold){
  lim <- 0.02
  # fp
  plot(0,0,cex = 0 , ylim=c(0,lim) , xlim=c(-0.5,1.5),xaxt ="n",xlab = "",ylab = "false positive error")
  pathcols <- rainbow(length(training_groups))
  for(path in training_groups){
    points(c(0,1),c(err_list[[m1]][[pn]][fold,path],
                    err_list[[m2]][[pn]][fold,path]),
           lwd = class_weights[path,"cost"]/4 ,type = "b",
           col = pathcols[match(path, training_groups)])
    text(-0.1 , err_list[[m1]][[pn]][fold,path] ,pos = 2 ,labels = path)
    text(1.1 , err_list[[m1]][[pn]][fold,path] ,pos = 4 ,labels = path)
    axis(1,at = c(0,1) , labels = c(m1,m2))
  }
}

plot_prediction_squares <- function(pred.top , pred.low  ,main = "" ,annot = NULL  ,mar =NULL,
                                    indicator  = NULL, cutoff = NA ,batches = F ,shift_text = 10,
                                    colmat = NULL,cluster_sams = F , axis_midpoints = F,annot_2 = NULL,
                                    vertlines = T,mark_prob =F){
  bacteria <- c("ecoli","gas","gbs","meningococcal","pneumo","staph")
  inflammatory <- c("HSP","JIA","SLE")
  viral <- c("adeno","enterovirus","flu","HHV6","rhino","RSV")
  
  if(is.null(annot)){
    par(mar = c(5,8,1,1) )  
  }else{
    par(mar = c(10,8,1,1) )
  }
  if(!is.null(mar)){
    par(mar = mar )
  }
  pred.top <- pred.top[,c("bacterial","viral","inflammatory","malaria","TB","KD")]
  top_order <- c("bacterial","viral","inflammatory","malaria_top","TB_top","KD_top")
  colnames(pred.top) <-  top_order
  low_order <- c(bacteria,viral,inflammatory,"malaria","TB","KD")
  pred.resp <- cbind(pred.low[,low_order], 0 ,pred.top[,top_order])
  colnames(pred.resp)[19] <- ""
  if(is.null(indicator)){
    indicator <- matrix(nrow = nrow(pred.top), ncol = ncol(pred.top) + ncol(pred.low) + 1)
    rownames(indicator) <- rownames(pred.top)
    colnames(indicator) <- c(colnames(pred.low),"",colnames(pred.top))
    indicator[] <- 0
    notrue <- T
  }else{
    notrue <- F
    reordered_sams <- c()
    for(i in low_order){
      reordered_sams <- c(reordered_sams , rownames(pred.resp)[indicator[rownames(pred.resp),i] ==1  ])
    }
    pred.resp <- pred.resp[reordered_sams,]
  }
  if(cluster_sams == T){
    for(class in colnames(indicator)[1:18]){
      samples <- rownames(pred.top)
      samples <- samples[indicator[samples,class] == 1]
      dist <- dist(cbind(pred.low[samples,] , pred.top[samples,]))
      positions <- match(samples ,reordered_sams  )
      reordered_sams[positions] <- samples[hclust(dist)$order]
    }
    pred.resp <- pred.resp[reordered_sams,]
  }
  sams <- rownames(pred.resp)
  col_range <- colorRampPalette(c('white','green'))
  col_range_true <- col_range(250)
  col_range <- colorRampPalette(c('white','blue'))
  col_range_false <- col_range(250)
  col_range <- colorRampPalette(c('white','orange'))
  col_range_close <- col_range(250)
  
  layout(matrix(1, nrow=1)) 
  plot(-10,-10,ylim=c(-1,ncol(pred.resp))+0.5,xlim=c(0,length(sams))+0.5 ,yaxt = 'n',xaxt = 'n',xlab="",ylab = "",xaxs = "i" , yaxs = "i")
  
  if(batches){
    axis(2, at=0, labels="Batch",las = 2,tick = F)
  }
  axis(2, at=1:ncol(pred.resp), labels=gsub("_.*" ,"" ,colnames(pred.resp)),las = 2,tick = F)
  if(!is.null(annot)){
    axis(1, at=1:nrow(pred.resp), labels=annot,las = 2,tick = F)
  }
  rect(xleft = 0.5 , xright = length(sams) + 0.5, ytop = 19.5 , ybottom = 18.5,
       col = "white" , lwd = 0)
  phem1 <- ""
  phe_X <- 0.5
  for(i in 1:length(sams)){
    sam <- sams[i]
    for(j in 1:ncol(pred.resp)){
      rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
           col = col_range_false[round(249 * pred.resp[i,j])+1 ],border = NA) 
      # close calls 
      if(colnames(pred.resp)[j] %in% low_order){
        pred <- colnames(pred.resp)[j]
        if(!notrue){
          truth <- colnames(indicator[,low_order])[indicator[sam, low_order] == 1]
          if(pred %in% bacteria & truth %in% bacteria 
             | pred %in% viral & truth %in% viral 
             | pred %in% inflammatory & truth %in% inflammatory ){
            rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
                 col = col_range_close[round(249 * pred.resp[i,j]) + 1] , lwd = 0) 
          }
        }
      }
      # true classes
      if(colnames(pred.resp)[j]!= "" & !notrue){
        if(indicator[sam , gsub("_.*" ,"" ,colnames(pred.resp))[j]] == 1 ){
          rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
               col = col_range_true[round(249 * pred.resp[i,j]) + 1] , lwd =0) 
        } 
      }
      
      # colouring for uncertain / co
      if(!is.null(colmat)){
        col_range <- colorRampPalette(c('white',colmat[i,j]))
        col_range.loc <- col_range(250)
        rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
             col = col_range.loc[round(249 * pred.resp[i,j]) + 1] , lwd = 0) 
        if(colmat[i,j] == "green"){
          rect(xleft = i-0.05 , xright = i + 0.05, ytop = j + 0.25 , ybottom = j - 0.25,
               col = "black" , lwd = 0)
        }
        
      }
      # annotate with probabilities
      if(mark_prob == T){
        if(pred.resp[i,j] > 0){
          text(i,j , round(pred.resp[i,j],digits = 2),pos = 2, cex =0.5,offset = 4)
        }
      }
      # max marks
      if(colnames(pred.resp)[j] %in% low_order){
        if(is.na(cutoff)){
          if(pred.resp[i,j] == max(pred.resp[i,low_order])){  
            rect(xleft = i-0.25 , xright = i + 0.25, ytop = j + 0.05 , ybottom = j - 0.05,
                 col = "black" , lwd = 0)
          }  
        }else{
          if(pred.resp[i,j]  >= cutoff){  
            rect(xleft = i-0.25 , xright = i + 0.25, ytop = j + 0.05 , ybottom = j - 0.05,
                 col = "black" , lwd = 0)
          }
        }
      }else{
        if(pred.resp[i,j] == max(pred.resp[i,top_order])){
          rect(xleft = i-0.25 , xright = i + 0.25, ytop = j + 0.05 , ybottom = j - 0.05,
               col = "black" , lwd = 0)
        }  
      }
    }
    # group separators
    if(phem1 != phenotypes[sam , "group"]){
      if(vertlines){
        abline( v = i -0.5 )  
      }
      midpoint <- (phe_X + i -0.5 )  / 2 
      if(axis_midpoints){
        axis(1,at = midpoint , las = 2 , labels = phem1,cex = 0.75)
      }
      
      phe_X <- i -0.5
    } 
    phem1 <- phenotypes[sam , "group"]
    # batch marks 
    if(batches == T){
      rect(xleft = i - 0.5 , xright = i + 0.5, ytop = -0.5 , ybottom = 0.5,
           col = rainbow(20)[phenotypes[sam,"batch"]], lwd = 0)
    }
  }
  # final x lab
  midpoint <- (phe_X + i -0.5 )  / 2 
  if(axis_midpoints){
    axis(1,at = midpoint , las = 2 , labels = phem1,cex = 0.75)  
  }
  abline(h = seq(0.5,ncol(pred.resp),by = 1),lwd  = 1,col = "grey")
  abline(h = 18.5,lwd  = 2.5)
  abline(h = 19.5,lwd  = 2.5)
  rect(xleft = 0.5,xright = nrow(pred.resp)+0.5 , ytop = 19.5 , ybottom = 18.5,lwd = 2,col = "white")
  abline(h = 6.5,lwd  = 2)
  abline(h = 12.5,lwd  = 2)
  abline(h = 15.5,lwd  = 2)
  abline(h = 16.5,lwd  = 2)
  abline(h = 17.5,lwd  = 2)
  for(i in 1:length(sams)){
    # unseen pathogens 
    if(!is.null(annot_2)){
      text(i-0.2 ,shift_text ,labels = annot_2[i]  ,pos = 3,srt = 90)
    }
  }
}


plot_rnaseq_squares <- function(pred.top  ,pred.low ,main = "" ,annot = NULL  ,plot_order.low, plot_order.top, # = c("bacterial","viral","inflammatory","TB","KD","malaria"),
                                indicator  = NULL, cutoff = NA ,batches = F , colmat = NULL,notrue =F, control =F){
  if(is.null(annot)){
    par(mar = c(5,8,1,1) )  
  }else{
    par(mar = c(10,8,1,1) )
  }
  reordered_sams <- c()
  for(i in plot_order.low){
    reordered_sams <- c(reordered_sams , rownames(pred.low)[indicator[rownames(pred.low),i] ==1  ])
  }
  pred.top <- pred.top[reordered_sams,plot_order.top]
  pred.low <- pred.low[reordered_sams,plot_order.low]
  sams <- reordered_sams
  col_range <- colorRampPalette(c('blue','white','green'))
  col_range_true <- colorRamp2(breaks = c(-1,0,1),colors = col_range(3)) 
  col_range <- colorRampPalette(c('blue','white','red'))
  col_range_false <- colorRamp2(breaks = c(-1,0,1),colors = col_range(3)) 
  col_range <- colorRampPalette(c('blue','white','yellow'))
  col_range_close <- colorRamp2(breaks = c(-1,0,1),colors = col_range(3))  # blue for negative vals
  layout(matrix(1, nrow=1)) 
  plot(-10,-10,ylim=c(0,ncol(pred.low) + ncol(pred.top) + 1 )+0.5,xlim=c(0,length(sams))+0.5,
       yaxt = 'n',xaxt = 'n',xlab="",ylab = "",xaxs = "i" , yaxs = "i")
  axis(2, at=1:ncol(pred.low), labels=gsub("_.*" ,"" ,colnames(pred.low)),las = 2,tick = F)  
  axis(2, at=(ncol(pred.low) + 2):(ncol(pred.low)+ncol(pred.top)+1), labels=gsub("_.*" ,"" ,colnames(pred.top)),las = 2,tick = F)  
  if(!is.null(annot)){
    names(annot) <- rownames(pred.top)
    axis(1, at=1:nrow(pred.top), labels=annot[reordered_sams],las = 2,tick = F)
  }
  rect(xleft = 0.5 , xright = length(sams) + 0.5, ytop = 19.5 , ybottom = 18.5,
       col = "white" , lwd = 0)
  phem1 <- ""
  for(i in 1:length(sams)){
    sam <- sams[i]
    for(j in 1:ncol(pred.low)){
      rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
           col = col_range_false(pred.low[i,j]),border = NA) 
      # true classes
      if(indicator[sam , gsub("_.*" ,"" ,colnames(pred.low))[j]] == 1 ){
        rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
             col = col_range_true( pred.low[i,j]) , lwd = 0) 
        rect(xleft = i-0.5 , xright = i + 0.5, ytop = j + 0.5 , ybottom = j - 0.5,
             lwd = 1) 
      } 
      if(pred.low[i,j] == max(pred.low[i,])){
        rect(xleft = i-0.25 , xright = i + 0.25, ytop = j + 0.05 , ybottom = j - 0.05,
             col = "black" , lwd = 0)
      }  
    }
    # toplevel 
    j_base <- ncol(pred.low)+1
    for(j in 1:ncol(pred.top)){
      rect(xleft = i-0.5 , xright = i + 0.5, ytop = j_base + j + 0.5 , ybottom = j_base + j - 0.5,
           col = col_range_false(pred.top[i,j]),border = NA) 
      if(indicator[sam , gsub("_.*" ,"" ,colnames(pred.top))[j]] == 1 ){
        rect(xleft = i-0.5 , xright = i + 0.5, ytop = j_base + j + 0.5 , ybottom = j_base + j - 0.5,
             col = col_range_true( pred.top[i,j]) , lwd = 0) 
        rect(xleft = i-0.5 , xright = i + 0.5, ytop = j_base + j + 0.5 , ybottom = j_base + j - 0.5,
             lwd = 1) 
      } 
      if(pred.top[i,j] == max(pred.top[i,])){
        rect(xleft = i-0.25 , xright = i + 0.25, ytop = j_base + j + 0.05 , ybottom = j_base + j - 0.05,
             col = "black" , lwd = 0)
      }  
    }
  }
}



plot_ridgeline <- function(pm,im,rl,sf = (1/20)){
  colord <- c(bacteria,"bacterial",viral,"viral",inflammatory,"inflammatory","TB","KD","malaria")
  colord <- colord[colord%in% colnames(pm)]
  layout(matrix(c(rep(seq(1,length(colord)),4), seq(length(colord)+1,length(colord)*2)),
                nrow = 5, ncol = length(colord), byrow = TRUE))
  for(class in colord){
    pretty_plot_area(ytick = c(0,ncol(pm)+2,1),cols = c("white","lightgray"),
                     xtick = c(0,1,0.2),margins = c(2,0.5,3,0.5),
                     show_x = T,show_y = F,show_y_tick = F)
    if(class == colord[1]) mtext(side = 2 , colord, at = 1:length(colord))
    title(colour_table[class,"name"],line = 0.5)
    df <- data.frame()
    for(i in rev(1:length(colord))){
      pclass <-  colord[i]
      cases <- rownames(im)[im[,pclass]==1]
      d <- density(pm[cases,class])
      df[cases,1] <- pm[cases,class]
      df[cases,2] <- pclass
      df[cases,3] <- colour_table[pclass, "col"]
      if(class == pclass){
        sf.m <-  length(cases) * sf
        polygon(d$x,((sf * d$y )/ (max(d$y)*1.2))+i - 0.05 ,
                col = add.alpha(colour_table[pclass,"col"],0.7),lwd = 2, border = T)
      }else{
        sf.m <-  length(cases) * sf
        
        polygon(d$x,( (sf.m * d$y )/ (max(d$y)*1.2))+i - 0.05 ,
                col = add.alpha(colour_table[pclass,"col"],0.7),lwd = 2, border = F)
      }
    }
    colnames(df) <- c("val","class","col")
    df$class <- factor(df$class,levels = colord)
    beeswarm(val ~ class, data = df , horizontal = T,
             add = T, side = 1, spacing = 0.5, pch = 16, corral = "wrap" ) 
  }
  for(class in colord){
    pretty_plot_area(ytick = c(0,1,0.5),cols = c("white","lightgray"),
                     xtick = c(1,0,-0.5),margins = c(2.5,2.5,3,0.5),
                     show_x = T,show_y = F,show_y_tick = T)
    plot(rl$ROCs[[class]], add = T, col =colour_table[class, "col"], lwd = 2)
    points(c(1,0), c(0,1), col = "lightgray", type = "l")
  }
  par(mfrow = c(1,1))
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



enrichment_plot <- function(src, ylim = c(-20,20), ytextbuffer =0,
                            textjump = 0.4,textsize = 0.6, yjump = 5,rotate = 0){
  pretty_plot_area(ytick = c(ylim,yjump),ybuffer = 5,cols = c("white","lightgray"),margins = rep(6,4),
                   xtick = c(0.8,7,0.2),show_x_tick = F,show_y_tick = T,absolute_y = T,
                   show_x = F,show_y = T,plot_outside_margin = T,y_lab = "-log10(padj)", main = src)
  for(i in 1:6){
    disease <- c("viral","bacterial","JIA","KD","TB","malaria")[i]
    rect(xleft = i-0.15,xright = i-0.05, ytop = 2, ybottom = -2, col = "white", border = NA)
    text(i-0.1,0, colour_table[disease,"name"], col = colour_table[disease,"col"],srt = 90)
    
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
          text(i + 0.05, newypos[j]
               ,terms ,adj = 1,pos = 4,cex = textsize,srt = rotate)
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
        newypos <- newypos + ytextbuffer
        if(length(newypos) > 20) newypos <- newypos - newypos[text_lim] - 3
        for( j in 1:text_lim){
          points(c(i,i+0.06),c(log10(res[j,"p_value"]),newypos[j]), type = "l", col = "lightgrey")
          terms <- res[j,"term_name"]
          terms <- gsub("biological process involved in ","",terms)
          text(i + 0.05, newypos[j]
               , terms,adj = 1,pos = 4,cex = textsize,srt = rotate)
        }
        points(rep(i,nrow(res)),log10(res[,"p_value"])
               , pch = 16, cex = 1.5, col ="blue")
      }
    }
  }
}



enrichment_plot_2 <- function(enrichment_table,src){
  pretty_plot_area(ytick = c(1,nrow(enrichment_table),1),cols = c("white","lightgray"),margins = c(10,35,3,2),
                   xtick = c(1,6,1),show_x_tick = T,show_y_tick = T,ybuffer = 0.4,xbuffer = 0.4,
                   show_x = F,show_y = F,plot_outside_margin = T)
  text(3, nrow(enrichment_table)+2,src)
  sf <- 0.5
  for(i in 1:6){
    dis <- c("viral","bacterial","JIA","KD","TB","malaria")[i]
    
    if(dis %in% colnames(enrichment_table)){
      ind <- (1:nrow(enrichment_table))[!is.na(enrichment_table[,dis])]
      points(rep(i,length(ind)),ind,cex = sf  * sqrt(-log10(enrichment_table[ind,dis])), col = "black")    
    }
    if(paste0(dis,"_up") %in% colnames(enrichment_table)){
      ind <- (1:nrow(enrichment_table))[!is.na(enrichment_table[,paste0(dis,"_up")])]
      points(rep(i,length(ind)),ind,cex = sf * sqrt(-log10(enrichment_table[ind,paste0(dis,"_up")])), col = "red")    
    }
    
    if(paste0(dis,"_down") %in% colnames(enrichment_table)){
      ind <- (1:nrow(enrichment_table))[!is.na(enrichment_table[,paste0(dis,"_down")])]
      points(rep(i,length(ind)),ind,cex = sf  * sqrt(-log10(enrichment_table[ind,paste0(dis,"_down")])), col = "blue")    
    }
  }
  par(xpd = T)
  text(rep(0.5, nrow(enrichment_table)) , 1:nrow(enrichment_table),enrichment_table[,"term_name"],adj = 1)
  text(1:6,-0.5,srt = 45,labels =c("Viral vs Other","Bacterial vs Other","JIA vs Other",
                                   "KD vs Other","TB vs Other","Malaria vs Other") , adj = 1)
  points(rep(-2,4) ,(-1:-4)*1.5,cex = sf  * sqrt(seq(10,40,by =10)))
  text(rep(-2.5,4) ,(-1:-4)*1.5,labels = seq(10,40,by =10))
  text(-2,-0.5,"log10(p value)")
  
  points(rep(-4,3), lwd = 2 ,(-1:-3)*1.3,cex = sf  * sqrt(rep(10,3)), col = c("red","blue","black"))
  text( rep(-5,3),(-1:-3)*1.3,labels = c("Upregulated genes","Downregulated genes",
                                         "Differentially expressed genes"),adj = 1)
  
}


# pairs plot ####

null_finc <- function(...){
  
}

upper.panel <- function(x, y,samples,indicator.loc ,sf = 1,inmatrix,...){
  cases <- samples[indicator.loc[samples,x]==1]
  cases1 <- samples[indicator.loc[samples,y]==1]
  r1 <- roc(cases = inmatrix[cases,x],controls = inmatrix[cases1,x],direction = "<")
  r2 <- roc(cases = inmatrix[cases,y],controls = inmatrix[cases1,y],direction = ">")
  points(c(0,1),c(1,0),type = "l",col = "white")
  plot(r1,add = T,col = colour_table[x,"col"],lwd = 2)
  plot(r2,add = T,col = colour_table[y,"col"],lwd = 2)
  text(0.45,0.4, format(c(round(auc(r1),digits = 2),0.12))[1],cex = 3*sf,col = colour_table[x,"col"],adj = 0)
  text(0.45,0.18, format(c(round(auc(r2),digits = 2),0.12))[1],cex = 3*sf,col = colour_table[y,"col"],adj = 0)
  
}

upper.panel_one_thresh <- function(x, y,samples,indicator.loc ,sf = 1,
                                   inmatrix,conf = F, linecol = "white",...){
  cases <- samples[indicator.loc[samples,x]==1]
  cases1 <- samples[indicator.loc[samples,y]==1]
  r <- roc_diagonal(x = inmatrix[samples,x] ,
                    y = inmatrix[samples,y],
                    casesx = cases,
                    casesy = cases1, samples , conf = conf)
  points(c(0,1),c(1,0),type = "l",col = linecol)
  points(r$spec,r$sens,type = "l",col = "black",lwd = 2)
  if(conf) plot(ci.sp(r$rocobj, sensitivities=seq(0, 1, .01)), type="shape")
  if(conf){
    string <- paste("    ",round(r$auc,digits = 2),"\n(",round(ci(r$rocobj$auc)[1], digits = 2),"-",round(ci(r$rocobj$auc)[2], digits = 2),")",sep = "")
    if(ci(r$rocobj$auc)[1] == 1)  string <- gsub("\n\\(1-1\\)","",string)
    text(0.6,0.3, string,
         cex = 2*sf,adj = 0)
  }else{
    text(0.45,0.4, format(c(round(r$auc,digits = 2),0.12))[1],cex = 3*sf,adj = 0)
  }
}


roc_diagonal <- function(x,y,casesx , casesy, samples,conf){
  names(x) <- samples
  names(y) <- samples
  ns <- y / x # gradients
  ns <- ns[order(ns)] 
  sens <- c()
  spec <- c()
  for(n in ns){ # for each gradient
    pos <- names(y)[y >= n*x] #pos if above line 
    neg <- names(y)[y < n*x]
    tp <- length(intersect(casesy,pos))
    tn <- length(intersect(casesx,neg))
    fp <- length(intersect(casesx,pos))
    fn <- length(intersect(casesy,neg))
    sens <- c(sens , tp / (tp + fn))
    spec <- c(spec , tn / (tn + fp))
  }
  r.obj <- NULL
  if(conf)  r.obj <- roc(controls = ns[casesx], cases = ns[casesy],ci = T)
  return(list(sens = sens , spec = spec , auc = trapz(spec,sens), rocobj = r.obj))
} 

lower.panel<-function(x, y, samples,sf = 1, indicator.loc,classes,annot){
  if(!is.null(indicator.loc)){
    phe <- c()
    for(sam in samples){
      phe <- c(phe,classes[indicator.loc[sam,classes]==1])  
    }
    points(x,y, pch = 19, col = colour_table[phe,"col"],cex = sf)
  }else{
    points(x,y, pch = 19, col = annot[samples,"col"],cex = sf)
  }
}

bottom.panel <- function(x,topclass,samples,indicator.loc,sf = 1,scale_to_max,inmatrix){
  dens_list <- list()
  ymax <- 0
  for(class in colnames(inmatrix)){
    d <- density(x[indicator.loc[samples,class]==1])
    dens_list[[class]] <- d
    ymax <- max(c(ymax,d$y))
  }
  if(scale_to_max){
    pretty_plot_area(ytick = c(0,1,1),
                     xtick = c(0,1,1),show_x_tick = F, show_y_tick = F,plot_outside_margin = F,
                     show_x = F,show_y = F,margins = c(0.5,0.5,0.5,0.5)) # ,xbuffer = 0.01,ybuffer = 0.01)
  }else{
    pretty_plot_area(ytick = c(0,ceiling(ymax),1),
                     xtick = c(0,1,1),show_x_tick = F, show_y_tick = F,plot_outside_margin = F,
                     show_x = F,show_y = F,margins = c(0.5,0.5,0.5,0.5)) # ,xbuffer = 0.01,ybuffer = 0.01)
  }
  for(class in colnames(inmatrix)){
    if(scale_to_max){
      points(dens_list[[class]]$x,dens_list[[class]]$y / max(dens_list[[class]]$y) ,
             col = add.alpha(colour_table[class, "col"],0.6),
             type = "l",lwd = 2 * sf)   
    }else{
      points(dens_list[[class]]$x,dens_list[[class]]$y , col = add.alpha(colour_table[class, "col"],0.6),
             type = "l",lwd = 2 * sf)   
    }
  } 
  
  
}

middle.panel <- function(x,topclass,samples,sf = 1,indicator.loc,...){
  rect(0,0,1,1,col = colour_table[topclass,"col"])
  if(!is.null(indicator.loc)){
    text(0.5,0.5,paste(colour_table[topclass,"name"],"\n n =",colSums(indicator.loc)[topclass]),cex = 2 * sf)
  }else{
    text(0.5,0.5,paste(colour_table[topclass,"name"]),cex = 2 * sf)
  }
}

right1.panel <- function(topclass,samples,indicator.loc,sf = 1,inmatrix,...){ # boxplot / beeswarm
  classes <- colnames(inmatrix)
  df <- data.frame()
  par(xpd = F)
  for(k in 1:length(classes)){
    samples <- rownames(indicator.loc)[indicator.loc[,classes[k]]==1]
    df[samples,1] <- inmatrix[samples,topclass]
    df[samples,2] <- classes[k]
    df[samples,3] <- colour_table[classes[k],"col"]
    
    d <- density(inmatrix[samples,topclass])
    outline <- list( y = c(d$x, rev(d$x)),
                     x = c((d$y/(2.2*max(d$y)))+k , rev(-(d$y/(2.2*max(d$y)))+k))) 
    polygon(outline$x ,outline$y, col = colour_table[classes[k],"col"],border = F)
  }
  par(xpd = T)
  colnames(df) <- c("expr","class","cols")
  df[,2] <- factor(df[,2] , levels = classes)
  beeswarm(expr ~ class ,data = df,add = T,pch = 16,corral = "wrap",spacing =0.5,cex = 0.8 * sf,
           col= add.alpha("black",0.5))#,pwcol = df[,"cols"])  
}

right1.panel.small <- function(topclass,samples,indicator.loc,sf = 1,inmatrix,...){ # boxplot / beeswarm
  classes <- colnames(inmatrix)
  df <- data.frame()
  par(xpd = F)
  for(k in 1:length(classes)){
    samples <- rownames(indicator.loc)[indicator.loc[,classes[k]]==1]
    df[samples,1] <- inmatrix[samples,topclass]
    df[samples,2] <- classes[k]
    df[samples,3] <- colour_table[classes[k],"col"]
    points(rep(k,2),summary(inmatrix[samples,topclass])[c( "1st Qu.","3rd Qu.")],
           type = "l",lwd = sf*10, col = colour_table[classes[k],"col"])
    points(c(k-0.45,k+0.45),rep(summary(inmatrix[samples,topclass])[c( "Median")],2),
           type = "l",lwd = sf * 5, col = colour_table[classes[k],"col"])
    
    
  }
  par(xpd = T)
  colnames(df) <- c("expr","class","cols")
  df[,2] <- factor(df[,2] , levels = classes)
  beeswarm(expr ~ class ,data = df,add = T,pch = 16,corral = "wrap",spacing =0.5,cex = 0.8 * sf,
           col= add.alpha("black",0.5))#,pwcol = df[,"cols"])  
}


right2.panel <- function(topclass,samples,indicator.loc,sf = 1,inmatrix,
                         linecol = "white",conf,...){ # OVA roc
  par(xpd = T)
  points(c(0,1),c(1,0),type = "l",col = linecol)
  cases <- samples[indicator.loc[samples,topclass]==1]
  cases1 <- samples[indicator.loc[samples,topclass]!=1]
  r <- roc(cases = inmatrix[cases,topclass],controls = inmatrix[cases1,topclass],ci = conf)
  plot(r,add = T, col = colour_table[topclass,"col"],lwd = 2)
  if(conf) plot(ci.sp(r, sensitivities=seq(0, 1, .01)), type="shape", col = colour_table[topclass,"col"])
  if(conf){
    string <- paste("    ",round(r$auc,digits = 2)," \n(",round(ci(r$auc)[1], digits = 2),"-",round(ci(r$auc)[2], digits = 2),")",sep = "")
    if(ci(r$auc)[1] == 1) string <- gsub("\n\\(.*\\)","",string)
    text(0.6,0.3, string,
         cex = 2*sf,adj = 0, col = colour_table[topclass,"col"])
  }else{
    text(0.45,0.4, format(c(round(r$auc,digits = 2),0.12))[1],
         cex = 3*sf,adj = 0, col = colour_table[topclass,"col"])
  }
  par(xpd = F)
}


plot_custom_pairs <- function(inmatrix ,upper_panel ,middle_panel, lower_panel ,scatter_pair = T,
                              bottom_panel,right1_panel, right2_panel, indicator,sf = 1,margin_val = 0.5,
                              scale_to_max = F,plotcols = c("white","lightgrey"), axcex = 2,n_wide_boxplot = 2,regular = T,annot = NULL, conf = F){
  endind <- ncol(inmatrix)
  par(xpd = T)
  par(mfrow = c(endind,endind+2),oma=c(5,5,5,5)) # + 1
  layoutmat <- matrix(1:((endind+2)*endind), byrow = T,nrow = endind, ncol = endind + 2)
  if(regular){
    layout(layoutmat)
  }else{
    layout(layoutmat[,c(1:endind,rep(endind+1,n_wide_boxplot),endind + 2)])
  }
  
  for(i in 1:(endind)){ # i == row
    for(j in 1:(endind+2)){ # j  = column
      if(i > j & i <= endind){# lower SCATTER
        sams <- rownames(inmatrix)
        # filter to just the two groups involved:
        if(scatter_pair) sams <- sams[indicator[sams, colnames(inmatrix)[i]] == 1 | indicator[sams, colnames(inmatrix)[j]] == 1]
        pretty_plot_area(ytick = c(0,1,1),cols = plotcols,
                         xtick = c(0,1,1),show_x_tick = T, show_y_tick = T,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        lower_panel(inmatrix[sams,j],inmatrix[sams,i], sams ,sf = sf,annot=annot,
                    indicator.loc = indicator,classes = colnames(inmatrix)) # flip axes
        
      }else if(j > i & j <= endind ){ # upper ROC
        pretty_plot_area(ytick = c(0,1,1),cols = plotcols,
                         xtick = c(1,0,-1),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        upper_panel(colnames(inmatrix)[i],colnames(inmatrix)[j], rownames(inmatrix),inmatrix = inmatrix
                    ,indicator.loc = indicator,sf = sf,annot=annot,linecol = plotcols[2],conf)
        
      }else if(j == i & j!= endind +1){ # middle label
        pretty_plot_area(ytick = c(0,1,1),cols = plotcols,
                         xtick = c(0,1,1),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        middle_panel(inmatrix[,i] ,rownames(inmatrix), topclass = colnames(inmatrix)[i]
                     ,sf = sf,indicator.loc = indicator ,annot=annot)
        
      }else if(j == endind +1 & i<= endind ){ # boxplot - right 1
        pretty_plot_area(ytick = c(0,1,1),cols = plotcols,
                         xtick = c(0.5,ncol(inmatrix)+0.5,1),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        right1_panel(rownames(inmatrix), topclass = colnames(inmatrix)[i],inmatrix = inmatrix
                     ,indicator.loc = indicator, sf = sf ,annot=annot)
        
      }else if(j == endind +2 & i<= endind){ # ROC - right 2
        pretty_plot_area(ytick = c(0,1,1),cols = plotcols,
                         xtick = c(1,0,-1),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        right2_panel(samples = rownames(inmatrix), topclass = colnames(inmatrix)[i],inmatrix = inmatrix
                     ,indicator.loc = indicator ,sf = sf ,linecol = plotcols[2],annot=annot, conf)
        
      }else if(i == endind +1  & j<= endind){ # bottom - density
        bottom_panel(inmatrix[,j] ,rownames(inmatrix),inmatrix = inmatrix,
                     topclass = colnames(inmatrix)[j], indicator.loc = indicator,sf = sf,scale_to_max = scale_to_max)
        
      }else if(i == endind +1  & j == i){ # bottom right - class proportions
        class_prop <- colSums(indicator[rownames(inmatrix),colnames(inmatrix)])
        pretty_plot_area(ytick = c(0,max(class_prop),1),cols = plotcols,
                         xtick = c(0,ncol(inmatrix)+0.2,0.2),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        for( k in 1:length(class_prop)){
          rect(xleft = k-0.8,xright = k,
               ybottom = 0,ytop = class_prop[k],col = colour_table[names(class_prop)[k],"col"],border = F)
          text(k-0.4,class_prop[k]-5,class_prop[k],col = "white", cex = 2*sf)
        }
      }else if(i == endind +1  & j == endind + 2){ # bottom right
        pretty_plot_area(ytick = c(0,max(class_prop),1),cols = plotcols,
                         xtick = c(0,ncol(inmatrix)+0.2,0.2),show_x_tick = F, show_y_tick = F,plot_outside_margin = T,
                         show_x = F,show_y = F,margins = rep(margin_val,4)) # ,xbuffer = 0.01,ybuffer = 0.01)
        
      }
      
      
      if(i ==1 & j ==1 | i ==endind+1 & j == endind+2 | i ==1 & j == endind + 1){
        
      }else{
        if(i == endind & j <= endind) axis(1,las = 2,cex = sf,cex.axis = axcex,line = 0.5)
        if(i == 1) axis(3,las = 2,cex = sf,cex.axis = axcex,line = 0.5)
        if(j == 1) axis(2,cex = sf,cex.axis = axcex,line = 0.5)
        if(j == endind+2 ) axis(4,cex = sf,cex.axis = axcex,line = 0.5)
      }
      
    } 
  }
  
  par(mfrow = c(1,1))
}



# running functions #####################################################################

# MULTI #### 

run_fit_multinomial <- function(target_gene_no,de_tables,training_sams.local,parallel = F,modifiers = list(cost = T, rebalance = T),
                                left_out_fold.loc=left_out_fold,fit_type = "multinomial",comp){
  
  target_gene_no.loc <- target_gene_no
  de_tables.loc <- de_tables
  genes <- pool_de_lists(target_gene_no = target_gene_no.loc,de_tables = de_tables.loc , comparisons.loc = comp)
  training_genes <- unique(unlist(genes))
  tmp_folds <- training_phenotypes[training_sams.local,"fold"]
  tmp_folds[tmp_folds > left_out_fold.loc]  <- tmp_folds[tmp_folds > left_out_fold.loc] - 1 # remove gap in seq

  if(modifiers$cost == T & modifiers$rebalance == T){
    weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"cost+imbalance"]
  }else if(modifiers$cost == T & modifiers$rebalance == F){
    weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"cost"]
  }else if(modifiers$cost == F & modifiers$rebalance == T){
    weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"imbalance"]
  }else if(modifiers$cost == F & modifiers$rebalance == F){
    weights.local <- rep(1,length(training_sams.local))
  }
  if(fit_type == "multiresponse"){
    cvfit <- cv.glmnet(t(expression_filt[training_genes , training_sams.local]) ,
                       training_phenotypes_indicator[training_sams.local , ],    
                       weights = weights.local, 
                       family="mgaussian",
                       alpha=1 ,type.measure = "mse", standardize = TRUE ,
                       foldid = tmp_folds,parallel = parallel)    
  }else if(fit_type == "multinomial"){
    cvfit <- cv.glmnet(t(expression_filt[training_genes , training_sams.local]) ,
                       training_phenotypes[training_sams.local , "group"],    
                       weights = weights.local, 
                       family="multinomial", type.multinomial = "grouped",
                       alpha=1 ,type.measure = "mse", standardize = TRUE ,
                       foldid = tmp_folds,parallel = parallel)
  }else{
    print("no fit type")
    break
  }
  return(cvfit)
}






choose_ngenes <- function(ngenes_cv,ngenes_range){ # select number of genes to prefilter to .. by min mse
  mses <-  unlist(lapply(ngenes_cv , function(i) i$cvm[match(i$lambda.1se,i$lambda)]))
  nzeros <-  unlist(lapply(ngenes_cv , function(i) i$nzero[match(i$lambda.1se,i$lambda)]))
  out <- list()
  out[["ngenes"]] <- ngenes_range[match(min(mses),mses)]
  out[["mse"]] <- mses
  out[["nzeros"]] <- nzeros
  out[["ngenes_range"]] <- ngenes_range
  return(out)
}



train_multinomial <- function(train_sams_local,ncores,modifiers , de_tables.loc = NULL,ngenes_range =NULL
                              ,left_out_fold = left_out_fold){
  print(training_groups)
  comparisons <- combn(training_groups,2)
  if(is.null(de_tables.loc)){
    print("running DE")
    de_tables.loc <- mclapply(1:ncol(comparisons), all_de_comparisons ,
                              training_sams.local = train_sams_local,comparisons.loc = comparisons,mc.cores =  ncores)
  }
  if(is.null(ngenes_range)){
    ngenes_range <- rev(seq(1000,10000,by=500))  # bin packing   
  }
  print("cv - tune prefilter ")
  cvfits <- mclapply(ngenes_range  , run_fit_multinomial , training_sams.local = train_sams_local,
                     de_tables = de_tables.loc , mc.cores = ncores,modifiers = modifiers,
                     left_out_fold.loc = left_out_fold, comp = comparisons)
  prefilt_tune <- choose_ngenes(cvfits , ngenes_range)
  ngenes_no_select <- prefilt_tune$ngenes
  print(c(ngenes_no_select,"genes sel by prefilt"))
  final_fit <- run_fit_multinomial(ngenes_no_select , de_tables.loc , 
                                   train_sams_local , parallel = ncores > 1 ,modifiers = modifiers,
                                   left_out_fold.loc = left_out_fold, comp = comparisons)
  genes <- pool_de_lists(target_gene_no = ngenes_no_select,de_tables = de_tables.loc ,
                         comparisons.loc = comparisons)
  training_genes <- unique(unlist(genes)) # will need in predictions !!
  print("done")
  return(list(fit = final_fit , input_genes = training_genes ,DE = de_tables.loc,
              prefilt_tune = prefilt_tune))
}


predict_modified <- function(cvfit, test_sams.local ,silence_features = NULL , lambda = "lambda.1se" ,
                             class_priors = NULL ,noprior = T,decouple = F){
  coefficients <- coef(cvfit,s=lambda)
   if(is.null(class_priors)){
    training_groups <- names(coefficients)
    class_priors.uniform <- data.frame()
    class_priors.uniform[ training_groups,1] <- 1/length(training_groups)
    class_priors <- class_priors.uniform
  }
  pred.resp <- data.frame(matrix(nrow = length(test_sams.local) , ncol = length(names(coefficients))))
  colnames(pred.resp) <- names(coefficients)
  rownames(pred.resp) <- test_sams.local
  exponents <- pred.resp
  pred.decoupled <- pred.resp
  for(path in names(coefficients)){
    betas <- coefficients[[path]][which(as.numeric(coefficients[[path]])!=0)] 
    probes  <- rownames(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
    coef.loc <- data.frame(matrix(nrow = length(betas)))
    rownames(coef.loc) <- probes
    coef.loc[probes,] <- betas
    probes <- probes[probes!="(Intercept)"]
    xnew <- expression_filt[probes,test_sams.local]
    var1 <-  exp(coef.loc["(Intercept)",] + colSums( coef.loc[probes,] * xnew[probes,]  )  # e ^ ( B0 + B1 + B2 + ...) 
                 + log(class_priors[path,1] / (1 - class_priors[path,1] ) )) # bias to include prior
    
    if(noprior ==T){
      var1 <-  exp(coef.loc["(Intercept)",] + colSums( coef.loc[probes,] * xnew[probes,]  ) )
      # no prior included
    }
    exponents[names(var1),path] <- var1
    pred.decoupled[names(var1),path] <- var1 / (1 + var1) # each output is an independent logit
  }
  pred.resp <-  exponents / (1 + rowSums(exponents))
  if(decouple){
    pred.resp <-  pred.decoupled
  }
  return(pred.resp)
  
}





predict.mod2 <- function(cvfit, test_expr ,silence_features = NULL , lambda = "lambda.min" ,
                         class_priors = NULL ,noprior = T,decouple = F){
  coefficients <- coef(cvfit,s=lambda)
  if(is.null(class_priors)){
    training_groups <- names(coefficients)
    class_priors.uniform <- data.frame()
    class_priors.uniform[ training_groups,1] <- 1/length(training_groups)
    class_priors <- class_priors.uniform
  }
  pred.resp <- data.frame(matrix(nrow = ncol(test_expr) , ncol = length(names(coefficients))))
  colnames(pred.resp) <- names(coefficients)
  rownames(pred.resp) <- colnames(test_expr)
  exponents <- pred.resp
  pred.decoupled <- pred.resp
  for(path in names(coefficients)){
    betas <- coefficients[[path]][which(as.numeric(coefficients[[path]])!=0)] 
    probes  <- rownames(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
    coef.loc <- data.frame(matrix(nrow = length(betas)))
    rownames(coef.loc) <- probes
    coef.loc[probes,] <- betas
    probes <- probes[probes!="(Intercept)"]
    xnew <- t(scale(t(test_expr[probes,])))
    xnew[is.na(xnew)] <- 0
    var1 <-  exp(coef.loc["(Intercept)",] + colSums( coef.loc[probes,] * xnew[probes,]  )  # e ^ ( B0 + B1 + B2 + ...) 
                 + log(class_priors[path,1] / (1 - class_priors[path,1] ) )) # bias to include prior
    if(noprior ==T){
      var1 <-  exp( colSums( coef.loc[probes,] * xnew[probes,]  ) )
      # no prior included + no intercept
    }
    exponents[names(var1),path] <- var1
    pred.decoupled[names(var1),path] <- var1 / (1 + var1) # each output is an independent logit
  }
  pred.resp <-  exponents / (1 + rowSums(exponents))
  if(decouple){
    pred.resp <-  pred.decoupled
  }
  return(pred.resp)
  
}


pred_refit_rnaseq <- function(samples,fit,type = c("decoupled.logit","logit","response") ,
                              lambda = "lambda.min",expr.loc = expr.mod.raw,scale = T){
  coefficients <- coef(fit,s=lambda)
  pred.resp <- data.frame(matrix(nrow = length(samples) , ncol = length(names(coefficients))))
  colnames(pred.resp) <- names(coefficients)
  rownames(pred.resp) <- samples
  exponents <- pred.resp
  pred.decoupled <- pred.resp
  pred.linmod <- pred.resp
  nullprior <- 1 / length(coefficients)
  probes  <- rownames(coefficients[[1]])[as.numeric(coefficients[[1]])!=0]
  probes <- probes[probes!="(Intercept)"]
  xnew <- expr.loc[probes,samples] # not scaling (as for microarrays)
  xnew[is.na(xnew)] <- 0
  for(path in names(coefficients)){
    betas <- coefficients[[path]][which(as.numeric(coefficients[[path]])!=0)] 
    probes  <- rownames(coefficients[[path]])[as.numeric(coefficients[[path]])!=0]
    coef.loc <- data.frame(matrix(nrow = length(betas)))
    rownames(coef.loc) <- probes
    coef.loc[probes,] <- betas
    probes <- probes[probes!="(Intercept)"]
    var1 <-  exp( colSums( coef.loc[probes,] * xnew[probes,]  ) + coef.loc["(Intercept)",] )
    exponents[names(var1),path] <- var1
    pred.decoupled[names(var1),path] <- var1 / (1 + var1) # each output is an independent logit
    pred.linmod[names(var1),path] <- colSums(coef.loc[probes,] * xnew[probes,] )
  }
  pred.resp <-  exponents / (1 + rowSums(exponents))
  if(type == "decoupled.logit"){
    return(pred.decoupled)
  }else if(type == "logit"){
    return(pred.resp)  
  }else if(type == "response"){
    return(pred.linmod)
  }
  
}




# OVA #########################################################################

# run a single fit with a given number of genes as input :
run_fit_one_vs_all <- function(target_gene_no,de_tables,training_sams.local,parallel = F,comp,
                               training_groups.loc = training_groups,modifiers = NULL,left_out_fold = left_out_fold){
  
  genes <- pool_de_lists(target_gene_no = target_gene_no,de_tables = de_tables , comparisons.loc = comp)
  training_genes <- unique(unlist(genes))
  
  tmp_folds <- training_phenotypes[training_sams.local,"fold"]
  tmp_folds[tmp_folds > left_out_fold]  <- tmp_folds[tmp_folds > left_out_fold] - 1 # remove gap in seq
  onevall_fit <- list()
  for( path in training_groups.loc ){
    cases <- training_sams.local[training_phenotypes[training_sams.local,"group"]==path]
    controls <- training_sams.local[training_phenotypes[training_sams.local,"group"]!=path]
    sams_tmp <- c(cases,controls)
    phenos <- c(rep(1,length(cases)),rep(0,length(controls)))
    weights.tmp <- c(rep(1/length(cases) , length(cases)),rep(1/length(controls),length(controls))) # only rebalancing.3
    
    cvfit <- cv.glmnet(t(expression_filt[training_genes , sams_tmp] ),phenos ,
                       weights = weights.tmp ,alpha = 1  ,
                       type.measure = "mse",standardize =T , foldid = tmp_folds , parallel = parallel)
    onevall_fit[[path]] <- cvfit
  }
  
  onevall_fit[["input_genes"]] <- training_genes
  sigs <- extract_geneno_onevall(onevall_fit,training_groups.loc)
  train_dat <- t(expression_filt[unique(unlist(sigs)) ,training_sams.local ])
  resp <- as.factor(training_phenotypes[training_sams.local,"group"])
  train_dat2 <- as.data.frame(cbind(train_dat , resp))
  resp <-  training_phenotypes[training_sams.local,"group"]
  num_iterations <- 1000  # <- might slow by having too high but poss necessary
  ## FF
  onevall_fit[["ridge_FF"]] <- cv.glmnet(train_dat , resp ,alpha = 0  ,family = "multinomial",
                         type.measure = "mse",standardize =T , foldid = tmp_folds )
  onevall_fit[["mn_FF"]] <- multinom(resp ~ . , train_dat2 , MaxNWts = 100000,maxit = num_iterations)
  ##TF
  weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"cost"]
  onevall_fit[["ridge_TF"]] <- cv.glmnet(train_dat , resp , weights = weights.local ,alpha = 0  ,family = "multinomial",
                            type.measure = "mse",standardize =T , foldid = tmp_folds )
  onevall_fit[["mnTF"]] <- multinom(resp ~ . , train_dat2 , MaxNWts = 100000,weights = weights.local,maxit = num_iterations)
  ##FT
  weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"imbalance"]
  onevall_fit[["ridge_FT"]] <- cv.glmnet(train_dat , resp , weights = weights.local ,alpha = 0  ,family = "multinomial",
                            type.measure = "mse",standardize =T , foldid = tmp_folds )
  onevall_fit[["mn_FT"]] <- multinom(resp ~ . , train_dat2 , MaxNWts = 100000,weights = weights.local,maxit = num_iterations)
  ##TT
  weights.local <- class_weights[training_phenotypes[training_sams.local,"group"],"cost+imbalance"]
  onevall_fit[["ridge_TT"]] <- cv.glmnet(train_dat , resp , weights = weights.local ,alpha = 0  ,family = "multinomial",
                            type.measure = "mse",standardize =T , foldid = tmp_folds )
  onevall_fit[["mn_TT"]] <- multinom(resp ~ . , train_dat2 , MaxNWts = 100000,weights = weights.local,maxit = num_iterations)
  onevall_fit[["sig_genes"]] <- unlist(sigs)
  return(onevall_fit)
}





extract_geneno_onevall <- function(onevall_fit,training_groups.loc){
  signatures <- list()
  for(path in training_groups.loc){
    fit <- onevall_fit[[path]]
    coefs <- coef(fit,s = "lambda.1se")[,1]
    genes <- names(coefs)[coefs!=0]
    genes <- genes[genes != "(Intercept)"]
    signatures[[path]] <- genes
  }
  return(signatures) 
}




choose_ngenes.onevsall <- function(ngenes_cv,ngenes_range){ # select number of genes to prefilter to .. by min mse
  mses <- c()
  nzeros <- c()
  for(i in 1: length(ngenes_cv)){
    cvfit <- ngenes_cv[[i]]    
    meanerrs <- c()
    pooled_sigs <- c()
    for(path in training_groups){
      fit <- cvfit[[path]]
      meanerr <- fit$cvm[fit$lambda == fit$lambda.1se] # mean error at 1se for cross validation of lambda in pathogen specific signature
      allcoef <- coef(fit,s = "lambda.1se")[,1]
      sig <- names(allcoef)[allcoef!=0]
      pooled_sigs <- c(pooled_sigs, sig )
      meanerrs <- c(meanerrs , meanerr)
    }
    print(sum(meanerrs)) ## sum mean errors ????
    print(length(unique(pooled_sigs)))
    
    mses <- c(mses , sum(meanerrs))
    nzeros <- c(nzeros , length(unique(pooled_sigs)))
  }
  out <- list()
  out[["ngenes"]] <- ngenes_range[match(min(mses),mses)]
  out[["mse"]] <- mses
  out[["nzeros"]] <- nzeros
  out[["ngenes_range"]] <- ngenes_range
  return(out)
}


train_onevsall <- function(train_sams_local,ncores = 1,modifiers , de_tables.loc = NULL,ngenes_range  = NULL
                           ,left_out_fold = left_out_fold){ # all these samples are used in the CV to find the prameters
  print(training_groups)
  comparisons <- combn(training_groups,2)
  if(is.null(de_tables.loc)){
    print("running DE")
    de_tables.loc  <- mclapply(1:ncol(comparisons), all_de_comparisons ,
                               training_sams.local = train_sams_local,comparisons.loc = comparisons,mc.cores =  ncores)
  }
  if(is.null(ngenes_range)){
    ngenes_range <- rev(seq(1000,10000,by=500))
  }
  print("cv - tune prefilter ")
  cvfits <- mclapply(ngenes_range  , run_fit_one_vs_all , training_sams.local = train_sams_local,
                     de_tables = de_tables.loc , mc.cores = ncores,modifiers = modifiers,
                     training_groups.loc = training_groups,left_out_fold = left_out_fold)
  
  prefilt_tune <- choose_ngenes.onevsall(cvfits , ngenes_range)
  ngenes_no_select <- prefilt_tune$ngenes
  print(c(ngenes_no_select,"genes selected by pre-filter"))
  final_fit <- run_fit_one_vs_all(ngenes_no_select , de_tables.loc , train_sams_local , parallel = ncores > 1 )
  pooled_sigs <- c()
  for(path in training_groups){
    fit <- final_fit[[path]]
    allcoef <- coef(fit,s = "lambda.1se")[,1]
    sig <- names(allcoef)[allcoef!=0]
    pooled_sigs <- c(pooled_sigs, sig )
  }
  print(length(unique(pooled_sigs)))
  genes <- pool_de_lists(target_gene_no = ngenes_no_select,de_tables = de_tables.loc , comparisons.loc = comparisons)
  training_genes <- unique(unlist(genes)) # will need in predictions !!
  print("done")
  return(list(fit = final_fit , input_genes = training_genes ,
              DE = de_tables.loc, prefilt_tune = prefilt_tune))
}




# predict ova  ####

predict_onevsall <- function(fit,test_sams.local){
  newx <-  t(as.matrix(expression_filt[fit$input_genes,test_sams.local]))
  pred.resp <- matrix(nrow =  nrow(newx) ,ncol = length(training_groups))
  rownames(pred.resp) <- rownames(newx)
  colnames(pred.resp) <- training_groups
  for(path in training_groups){
    fit.loc <- fit[[path]]
    pred.resp[,path] <- predict.cv.glmnet(fit.loc,newx = newx,s="lambda.1se")
  }
  return(pred.resp)
}

predict_onevsall_multinom <- function(fit,test_sams.local){
  newx <-  t(as.matrix(expression_filt[fit$sig_genes,test_sams.local]))
  ova.multinom.pred <- list()
  ova.multinom.pred[["TF"]] <- predict(fit$mn_TF, newx,type = "probs")
  ova.multinom.pred[["FT"]] <- predict(fit$mn_FT, newx,type = "probs")
  ova.multinom.pred[["TT"]] <- predict(fit$mn_TT, newx,type = "probs")
  ova.multinom.pred[["FF"]] <- predict(fit$mn_FF, newx,type = "probs")
  return(ova.multinom.pred)
}

reorder_predmat <- function(predmat){ predmat[,order(colnames(predmat))] }

# performance #####

pr_curve <- function(cases ,controls,invert = F){
  if(invert){ # in case cases not up from controls
    cases <- -cases
    controls <- -controls 
  }
  pnt <- c(cases, controls)
  pnt <- pnt[order(pnt)]
  out <- data.frame(matrix(nrow = length(pnt),ncol = 6))
  for(j in 1:length(pnt)){
    i <- pnt[j]
    tp <- length(cases[cases >= i])
    fn <- length(cases[cases < i])
    fp <- length(controls[controls >= i])
    tn <- length(controls[controls < i])
    precision <- tp/(tp + fp)
    recall <- tp/(tp+fn)
    accuracy <- (tp+tn)/(tp+fp+tn+fn)
    specificity <- tn/(tn+fp)
    f1 <- (2*tp) /((2*tp)+fp+fn) 
    out[j,] <- c(i , precision, recall,accuracy,specificity,f1)
  }
  colnames(out) <- c("threshold" , "precision", "recall","accuracy","specificity","f1")
  return(out)
}



make_confusion <- function(prediction_matrix, indicator_matrix,mod=NULL){
  confusion_base <- as.data.frame(matrix(nrow = ncol(indicator_matrix),ncol = ncol(indicator_matrix)))
  colnames(confusion_base) <- colnames(indicator_matrix)
  rownames(confusion_base) <- colnames(indicator_matrix)
  confusion_base[] <- 0
  samples <- rownames(indicator_matrix)
  confusion <- array(dim = c(nrow(confusion_base), ncol(confusion_base) , length(samples)))
  colnames(confusion) <- colnames(confusion_base)
  rownames(confusion) <- rownames(confusion_base)
  dimnames(confusion)[[3]] <- samples
  if(is.null(mod)){
    confusion[] <- 0 # for summing need non NA ;; for errr need to ignore
  }
  for(sam in samples){
    actual <- indicator_matrix[sam , colnames(confusion)]
    confusion[actual == 1 ,,sam] <-  as.numeric(prediction_matrix[sam,colnames(confusion)])
  }
  confusion.max <- confusion
  for( sam in samples){
    confusion.max[,,sam][ confusion[,,sam] == max(confusion[,,sam],na.rm = T) ] <- 1
    confusion.max[,,sam][ confusion[,,sam] != max(confusion[,,sam],na.rm = T) ] <- 0
  }
  confusion_sum <- rowSums(confusion,dims = 2)
  confusion_sum.max <- rowSums(confusion.max,dims = 2)
  confusions <- list(confusion_sum , confusion_sum.max)
  if(is.null(mod)){
    return(confusions)
  }else{
    return(list(values = confusion , maxvals = confusion.max))
  }
}





performance_calc <- function(resp_true,resp_pred){
  performance <- data.frame(matrix(nrow = ncol(resp_true) + 2 ,ncol = 9))
  rownames(performance) <- c(colnames(resp_true),"micro","macro")
  colnames(performance) <- c("precision","recall","accuracy","specificity","f1","tp","tn","fn","fp")
  for(i in 1:ncol(resp_true)){
    path <- colnames(resp_true)[i]
    group <- rownames(resp_true)[resp_true[,path]==1]
    tp <- sum(  resp_pred[group,path])
    fn <- - sum(  resp_pred[group,path] - 1 )
    ngroup <- rownames(resp_pred)[!rownames(resp_pred)%in%group]
    fp  <- sum(resp_pred[ngroup,path])
    tn <- - sum(  resp_pred[ngroup,path] - 1 )
    
    precision <- tp/(tp + fp)
    recall <- tp/(tp+fn)
    accuracy <- (tp+tn)/(tp+fp+tn+fn)
    specificity <- tn/(tn+fp)
    f1 <- (2*tp) /((2*tp)+fp+fn) 
    
    performance[path, ] <- c(precision,recall,accuracy,specificity,f1,tp,tn,fn,fp)
  }
  
  performance["macro",] <- colMeans(performance[colnames(resp_true),])
  performance["micro",c("tp","tn","fn","fp")] <- colSums(performance[colnames(resp_true),c("tp","tn","fn","fp")])
  tp <- performance["micro","tp"]
  tn <- performance["micro","tn"]
  fp <- performance["micro","fp"]
  fn <- performance["micro","fn"]
  precision <- tp/(tp + fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  specificity <- tn/(tn+fp)
  f1 <- (2*tp) /((2*tp)+fp+fn) 
  performance["micro", ] <- c(precision,recall,accuracy,specificity,f1,tp,tn,fn,fp)
  return(performance)
}




merge_seed_runs <- function(seed_performance){
  tmp <- seed_performance[[1]]
  for(i in names(seed_performance[[1]])){
    if(i!="seed"){
      for(j in 2:length(seed_performance)){
        tmp[[i]]$meansquared <- c(tmp[[i]]$meansquared , seed_performance[[j]][[i]]$meansquared) 
        tmp[[i]]$fp <- rbind(tmp[[i]]$fp , seed_performance[[j]][[i]]$fp)
        tmp[[i]]$fn <- rbind(tmp[[i]]$fn , seed_performance[[j]][[i]]$fn)
        tmp[[i]]$macromwse <- c(tmp[[i]]$macromwse , seed_performance[[j]][[i]]$macromwse) 
      } 
    }
  }
  return(tmp)
}


performance_weighted <- function(confusion,mod = "val",example_weights = NULL){
  if(mod == "val"){
    confusion_loc <- confusion$values  
  }else if(mod =="max"){
    confusion_loc <- confusion$maxvals
  }else{break}
  if(is.null(example_weights)){example_weights <- rep(1,length(dimnames(confusion_loc)[[3]]))}
  names(example_weights) <- dimnames(confusion_loc)[[3]]
  square_errs <- get_square_errors(confusion_loc )
  mse.all <- sum(square_errs * example_weights)
  se <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  errs <- get_fpfn_errors(confusion_loc,example_weights)
  # class errs 
  per_class_mse <- as.data.frame(matrix(nrow = nrow(confusion_loc[,,1]),ncol = 1))
  rownames(per_class_mse) <- dimnames(confusion_loc)[[1]]
  for(class in dimnames(confusion_loc)[[1]]){
    sams.loc <- c()
    for(sam in dimnames(confusion_loc)[[3]]){
      tc <- colnames(confusion_loc[,,sam])[!is.na(rowSums(confusion_loc[,,sam]))] # anchor
      if(tc == class) sams.loc <- c(sams.loc , sam)
    }
    mse <- sum(square_errs[sams.loc] * example_weights[sams.loc])
    per_class_mse[class,1] <-  mse
  }
  errors <- list(meansquared =  mse.all,
                 se = se,
                 fp  = errs$fp ,
                 fn = errs$fn,
                 allerr = square_errs,
                 per_class_mse = per_class_mse,
                 macromwse = mean(per_class_mse[,1]))
  return(errors)
}





get_square_errors <- function(confusion_loc){
  truecall <- diag(ncol(confusion_loc))
  for(sam in dimnames(confusion_loc)[[3]]){
    tc <- confusion_loc[,,sam][truecall ==1 ]
    confusion_loc[,,sam][truecall ==1 ] <- 1- tc
    confusion_loc[,,sam] <- confusion_loc[,,sam]^2
  }
  sample_errs <- aaply(confusion_loc , 3, sum ,na.rm = T) # sum the weighted square errors 
  return(sample_errs)
}

get_fpfn_errors <- function(confusion_loc,example_weights){
  truecall <- diag(ncol(confusion_loc))
  for(sam in dimnames(confusion_loc)[[3]]){
    tc <- confusion_loc[,,sam][truecall ==1 ]
    confusion_loc[,,sam][truecall ==1 ] <- 0 # if on diag is true 
    confusion_loc[,,sam] <- ( confusion_loc[,,sam]^2 ) * example_weights[sam]
  }
  fp  <- aaply(confusion_loc , 2, sum ,na.rm = T) 
  fn  <- aaply(confusion_loc , 1, sum ,na.rm = T) 
  errs <- list(fp = fp , fn = fn)
  return(errs)
}



shrink_MTT <- function(seedvalue){
  print(seedvalue)
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  indicator_matrix  <- class.ind(factor(training_phenotypes[rownames(prediction_matrix.multinomial_TT),"group"]))
  rownames(indicator_matrix) <- rownames(prediction_matrix.multinomial_TT)
  confusion_multinom_TT <- make_confusion(prediction_matrix.multinomial_TT , indicator_matrix,mod = "3dconfusion") 
  pred_blank <- prediction_matrix.multinomial_TT[,training_groups[order(training_groups)]]
  pred_blank[] <- 0 
  shrunk_predmats <- list(pred_blank , pred_blank ,pred_blank , pred_blank ,pred_blank , pred_blank )
  nfolds <- 10
  set.seed(seedvalue)
  for(grp in training_groups){
    grp_sams <- rownames(training_phenotypes)[training_phenotypes["group"]==grp]
    grp_sams <- sample(grp_sams,length(grp_sams))
    training_phenotypes[grp_sams,"fold"]   <- rep(1:nfolds,100)[1:length(grp_sams)]
    print(grp_sams )
  }
  for(left_out_fold in 1:10 ){
    training_sams.local <- rownames(training_phenotypes)[training_phenotypes[,"fold"] !=  left_out_fold]
    test_sams.local <- rownames(training_phenotypes)[training_phenotypes[,"fold"] ==  left_out_fold]
    cvfit <- all_fits[[left_out_fold]]$mn_TT
    msizes <- c(50,100,150,200,250,300)
    shrink_lambdas <- unlist(lapply(msizes , function(l) cvfit$lambda[cvfit$cvm==min(cvfit$cvm[cvfit$nzero < l])]))
    for( i in  1:length(msizes) ){
      lambda_val <- shrink_lambdas[i]
      pred.resp <- stats::predict(cvfit, 
                                  newx = t(as.matrix(expression_filt[all_fits[[left_out_fold]]$input_genes,test_sams.local])),
                                  s = lambda_val, type = "response")[,,1]
      shrunk_predmats[[i]][test_sams.local , colnames(pred.resp)] <- pred.resp
    }
  }
  names(shrunk_predmats) <- as.character(msizes)
  return(shrunk_predmats)
}



calc_performances <- function(seedvalue){
  print(seedvalue)
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  indicator_matrix  <- class.ind(factor(training_phenotypes[rownames(prediction_matrix.multinomial_TT),"group"]))
  rownames(indicator_matrix) <- rownames(prediction_matrix.multinomial_TT)
  confusion_multinom_TT <- make_confusion(prediction_matrix.multinomial_TT , indicator_matrix,mod = "3dconfusion") 
  ## shrink modelsizes with lambda 
  pred_blank <- prediction_matrix.multinomial_TT[,training_groups[order(training_groups)]]
  pred_blank[] <- 0 
  shrunk_predmats <- list(pred_blank , pred_blank ,pred_blank , pred_blank ,pred_blank , pred_blank )
  nfolds <- 10
  set.seed(seedvalue)
  for(grp in training_groups){
    grp_sams <- rownames(training_phenotypes)[training_phenotypes["group"]==grp]
    grp_sams <- sample(grp_sams,length(grp_sams))
    training_phenotypes[grp_sams,"fold"]   <- rep(1:nfolds,100)[1:length(grp_sams)]
  }
  confusion_list <- list(confusions = list() , names = c())
  # DIFF weights LASSO
  confusion_multinom_TF <- make_confusion(prediction_matrix.multinomial_TF , indicator_matrix,mod = "3dconfusion")   
  confusion_multinom_FT <- make_confusion(prediction_matrix.multinomial_FT , indicator_matrix,mod = "3dconfusion") 
  confusion_multinom_FF <- make_confusion(prediction_matrix.multinomial_FF , indicator_matrix,mod = "3dconfusion") 
  for(cm in list(confusion_multinom_TT,confusion_multinom_TF,
                 confusion_multinom_FT,confusion_multinom_FF)){
    confusion_list$confusions[[length(confusion_list$confusions) + 1]] <- cm 
  }
  confusion_list$names <- c(confusion_list$names , c("M_TT", 
                                                     "M_TF",
                                                     "M_FT"
                                                     ,"M_FF") )
  
  # OVA
  confusion_ova <- make_confusion(prediction_matrix.onevsall , indicator_matrix,mod = "3dconfusion") 
  confusion_list$confusions[[length(confusion_list$confusions) + 1]] <- cm 
  confusion_list$names <- c(confusion_list$names , c("OVA"))
  
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
  
  perfs <- list(seed = seedvalue)
  
  ew.cost <- class_weights[training_phenotypes[dimnames(confusion_multinom_TT$values)[[3]] , "group"] ,"cost"]
  ew.cost_imb <- class_weights[training_phenotypes[dimnames(confusion_multinom_TT$values)[[3]] , "group"] ,"cost+imbalance"]
  
  for(i in 1:length(confusion_list$confusions)){
    nam <- confusion_list$names[i]
    
    perfs[[paste(nam,".u.v",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="val" )
    perfs[[paste(nam,".u.m",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="max" )
    perfs[[paste(nam,".ew.v",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="val" ,example_weights = ew.cost_imb)
    perfs[[paste(nam,".ew.m",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="max" ,example_weights = ew.cost_imb)

    # cost weighting (no imbalance in mwse calculation)
    perfs[[paste(nam,".cw.v",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="val" ,example_weights = ew.cost)
    perfs[[paste(nam,".cw.m",sep = "")]] <- performance_weighted(confusion_list$confusions[[i]]  , mod ="max" ,example_weights = ew.cost)
    
  }
  
  if(length(  confusion_list$names ) != length(  confusion_list$confusions )){
    break # if names dont line up
  }
  return( perfs)
}






extract_seed_model_sizes <- function(seedvalue){
  print(seedvalue)
  load(paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  model_size_matrix <- matrix(nrow = 10, ncol = length(names(all_fits[[1]])) - 4)
  colnames(model_size_matrix) <- names(all_fits[[1]])[!names(all_fits[[1]])%in% c("train_samples","test_samples","input_genes","DE")]
  for(fold in 1:10){
    for(method in colnames(model_size_matrix)){
      if(method!="ova"){
        coefmat <- coef(all_fits[[fold]][[method]],s="lambda.1se")
        ms <- rownames(coefmat[[1]])[as.numeric(coefmat[[1]])!=0]
      }else{
        ms <- unique(unlist(extract_geneno_onevall(all_fits[[fold]][[method]],training_groups)))
      }
      model_size_matrix[fold,method] <- length(ms)
    }
  }
  return(model_size_matrix)
}


# pre-filtering #### 

all_de_comparisons <- function(i,comparisons.loc,training_sams.local){
  tab <- run_de(comparisons.loc[1,i],comparisons.loc[2,i] , training_sams.local )
  return(tab)
}


pool_de_lists <- function(target_gene_no,comparisons.loc, de_tables){
  tmp <- t(comparisons.loc)
  tmp <- cbind(tmp,target_gene_no / nrow(tmp))
  init <- target_gene_no / ncol(comparisons.loc)
  needed <- target_gene_no
  while(needed > 0){
    full_list <- lapply(1:ncol(comparisons.loc) ,function(i) rownames(de_tables[[i]])[1:as.numeric(tmp[i,3])])
    full_list <- lapply(1:ncol(comparisons.loc) , function(i) unique(full_list[[i]][!is.na(full_list[[i]])]))
    gene_weights <- 1 / table(unlist(full_list))
    reweighted <- unlist(lapply(full_list , function(i) sum(gene_weights[i])))
    extra <- init - reweighted
    extra[extra < 0 ] <- 0
    tmp[,3] <-  ceiling(as.numeric(tmp[,3]) + extra)
    needed <- sum(round(init - reweighted))
  }
  return(full_list)
}



run_de <- function(m1,m2,training_sams.local){
  a <- rownames(training_phenotypes)[training_phenotypes[,"group"]==m1]
  b <- rownames(training_phenotypes)[training_phenotypes[,"group"]==m2]
  a <- a[a%in%training_sams.local]
  b <- b[b%in%training_sams.local]
  group <- as.factor(c(rep(1,length(a)),rep(2,length(b))))
  age <- as.numeric(training_phenotypes[c(a,b),"age"])
  gender <- training_phenotypes[c(a,b),"gender"]
  design <- model.matrix(~ group + age  + gender)
  counts.filt <- expression_filt[,c(a,b)]
  fit <- lmFit(counts.filt,design)
  fit2 <- eBayes(fit,trend=TRUE)
  tt <- topTable(fit2,number = Inf, coef = 2)
  tt <- tt[abs(tt[,1])>0.5,]
  tt <- tt[order(tt[,"adj.P.Val"]),]
  tt <- cbind(tt,mappingInfo[rownames(tt),"Symbol"])
  return(tt)
}


# seed cross validation #####

fill_missing_age_mean <- function(training_phenotypes){
  missing_age <- rownames(training_phenotypes)[training_phenotypes[,"age"] == -40]
  notmissing_age <- rownames(training_phenotypes)[training_phenotypes[,"age"] != -40]
  for(path in unique(training_phenotypes[,"group"])){
    pathsams <- rownames(training_phenotypes)[training_phenotypes[,"group"] == path]
    training_phenotypes[intersect(missing_age , pathsams),"age"] <- mean(training_phenotypes[intersect(notmissing_age , pathsams),"age"])
  }
  return(training_phenotypes)
}

run_seed_cv <- function(seedvalue){
  set.seed(seedvalue)
  training_phenotypes <- base_training_phenotypes
  para <- F
  nfolds <- 10
  for(grp in training_groups){
    grp_sams <- rownames(training_phenotypes)[training_phenotypes["group"]==grp]
    grp_sams <- sample(grp_sams,length(grp_sams))
    training_phenotypes[grp_sams,"fold"]   <- rep(1:nfolds,100)[1:length(grp_sams)]
  }
  
  prediction_matrix <- data.frame(matrix(nrow = nrow(training_phenotypes),
                                         ncol = length(training_groups)))
  rownames(prediction_matrix) <- rownames(training_phenotypes)
  colnames(prediction_matrix) <- training_groups
  
  prediction_matrix.multinomial <- prediction_matrix
  prediction_matrix.multinomial_nnet <- prediction_matrix
  prediction_matrix.onevsall <- prediction_matrix
  prediction_matrix.multiresponse <- prediction_matrix
  prediction_matrix.multinomial_TT <- prediction_matrix
  prediction_matrix.multinomial_TF <- prediction_matrix
  prediction_matrix.multinomial_FT <- prediction_matrix
  prediction_matrix.multinomial_FF <- prediction_matrix
  prediction_matrix.onevsall.mn_TT <- prediction_matrix
  prediction_matrix.onevsall.mn_TF <- prediction_matrix
  prediction_matrix.onevsall.mn_FT <- prediction_matrix
  prediction_matrix.onevsall.mn_FF <- prediction_matrix
  
  fit_list_TT <- list()
  fit_list_TF <- list()
  fit_list_FT <- list()
  fit_list_FF <- list()
  fit_list_onevall <- list()
  comparisons <- combn(training_groups,2)
  ngenes_no_select <- 2000
  all_fits <- list() 
  
  for(left_out_fold in 1:max(training_phenotypes[,"fold"] ) ){
    all_fits[[left_out_fold]] <- list() 
    ncores <- 1
    print("-----------------------------------")
    print(paste("seed =" ,seedvalue , "fold =",left_out_fold))
    training_sams.local <- rownames(training_phenotypes)[training_phenotypes[,"fold"] !=  left_out_fold]
    test_sams.local <- rownames(training_phenotypes)[training_phenotypes[,"fold"] ==  left_out_fold]
    # DE
    print("run differential expression")
    de_tables.loc <- lapply(1:ncol(comparisons), all_de_comparisons ,
                            training_sams.local = training_sams.local,comparisons.loc = comparisons)
    # use same gene list for all methods
    print("pool these lists")
    genes <- pool_de_lists(target_gene_no = ngenes_no_select,de_tables = de_tables.loc , comparisons.loc = comparisons)
    training_genes <- unique(unlist(genes)) # will need in predictions !!
    all_fits[[left_out_fold]][["train_samples"]] <- training_sams.local
    all_fits[[left_out_fold]][["test_samples"]] <- test_sams.local
    # fit multinomial lasso
    print("mn TT")
    modifiers <-  list(cost = T, rebalance = T) 
    final_fit <- run_fit_multinomial(ngenes_no_select , de_tables.loc , training_sams.local , parallel = para ,
                                     modifiers = modifiers,left_out_fold = left_out_fold, comp = comparisons)
    pred.resp <- stats::predict(final_fit, 
                                newx = t(as.matrix(expression_filt[training_genes,test_sams.local])),
                                s = "lambda.1se", type = "response")[,,1]
    prediction_matrix.multinomial_TT[rownames(pred.resp),colnames(pred.resp)] <- pred.resp
    all_fits[[left_out_fold]][["mn_TT"]] <-   final_fit
    
    
    print("mn TF")
    modifiers <-  list(cost = T, rebalance = F)
    final_fit <- run_fit_multinomial(ngenes_no_select , de_tables.loc , training_sams.local , parallel = para ,
                                     modifiers = modifiers,left_out_fold = left_out_fold, comp = comparisons)
    pred.resp <- stats::predict(final_fit,
                                newx = t(as.matrix(expression_filt[training_genes,test_sams.local])),
                                s = "lambda.1se", type = "response")[,,1]
    prediction_matrix.multinomial_TF[rownames(pred.resp),colnames(pred.resp)] <- pred.resp
    all_fits[[left_out_fold]][["mn_TF"]] <-   final_fit
     
    print("mn FT")
    modifiers <-  list(cost = F, rebalance = T) 
    final_fit <- run_fit_multinomial(ngenes_no_select , de_tables.loc , training_sams.local , parallel = para ,
                                     modifiers = modifiers,left_out_fold = left_out_fold, comp = comparisons)
    pred.resp <- stats::predict(final_fit, 
                                newx = t(as.matrix(expression_filt[training_genes,test_sams.local])),
                                s = "lambda.1se", type = "response")[,,1]
    prediction_matrix.multinomial_FT[rownames(pred.resp),colnames(pred.resp)] <- pred.resp
    all_fits[[left_out_fold]][["mn_FT"]] <-   final_fit
    
    
    print("mn FF")
    modifiers <-  list(cost = F, rebalance = F)
    final_fit <- run_fit_multinomial(ngenes_no_select , de_tables.loc , training_sams.local , parallel = para
                                     ,modifiers = modifiers,left_out_fold = left_out_fold, comp = comparisons)
    pred.resp <- stats::predict(final_fit,
                                newx = t(as.matrix(expression_filt[training_genes,test_sams.local])),
                                s = "lambda.1se", type = "response")[,,1]
    prediction_matrix.multinomial_FF[rownames(pred.resp),colnames(pred.resp)] <- pred.resp
    all_fits[[left_out_fold]][["mn_FF"]] <-   final_fit
    
    ## onevsall 
    # train
    print("OVA")
    final_fit <- run_fit_one_vs_all(ngenes_no_select , de_tables.loc , training_sams.local , parallel = para ,
                                    left_out_fold = left_out_fold, comp = comparisons)
    pred.resp <- predict_onevsall(fit =  final_fit , test_sams.local = test_sams.local)
    prediction_matrix.onevsall[rownames(pred.resp),colnames(pred.resp)] <- pred.resp
    all_fits[[left_out_fold]][["ova"]] <- final_fit
    
    all_fits[[left_out_fold]][["input_genes"]] <- training_genes
    all_fits[[left_out_fold]][["DE"]] <- de_tables.loc
    
    save(prediction_matrix.onevsall , prediction_matrix.multinomial_FF, prediction_matrix.multinomial_TF,
         prediction_matrix.multinomial_FT , prediction_matrix.multinomial_TT,all_fits
         ,file = paste("data_objects/seed_cv/run_number_",seedvalue,".RData",sep=""))
  }
}



## ridge+lasso ########################################


run_ridges_on_lasso <- function(testfold,training_samples,folds ,
                                input_genes, lasso_lambdas,weight_vector=NULL,
                                alphaval = 0 ,ridge_lambdas ,mode ="LOW"){ # test fold is for 9 fold cv to tune the lambdas 
  print("running test fold")
  print(testfold)
  train_l1 <- training_samples[training_phenotypes[training_samples,"fold"]!=testfold]
  test_l1 <- training_samples[training_phenotypes[training_samples,"fold"]==testfold]
  if(is.null(weight_vector)){
    if(mode == "LOW"){
      weights.local <- class_weights[training_phenotypes[train_l1,"group"],"cost+imbalance"]  
      disease_set <- rownames(class_weights)
    }else if(mode == "TOP"){
      disease_set <- c("bacterial","viral","inflammatory","malaria","TB","KD")
      group <- c()
      for(sam in train_l1) group <- c(group , disease_set[training_phenotypes_indicator[sam , disease_set] == 1 ])
      weights.local <- class_weights_toplevel[group,"cost+imbalance"]
    }else{break}
  }else{
    weights.loc <- weight_vector[train_l1]
  }
  
  # the fit for each inner fold .. so no mse
  lasso_fit <- glmnet(t(expression_filt[input_genes , train_l1]),
                      y = training_phenotypes_indicator[train_l1,disease_set]
          ,family = "multinomial",lambda = lasso_lambdas,type.multinomial = "grouped"
          ,alpha = 1,standardize = T,weights = weights.local)

  print("fitted lasso")
  fold_prediction <- list()
  fold_prediction$lassofit <- lasso_fit
  for(l_lambda in lasso_lambdas){
    # lasso prediction which am comparing to : 
    p1 <- predict(lasso_fit,t(expression_filt[input_genes,test_l1]) , s = l_lambda,type="response")[,,1]
    coefs <- coef(lasso_fit,s=l_lambda)[[1]]
    selgenes <- rownames(coefs)[as.matrix(coefs)!=0]
    selgenes <- selgenes[selgenes!="(Intercept)"]
    p2 <- NULL
    ridge_fit <- NULL
    if(length(selgenes) > 2){
      print(length(selgenes))
      inner_folds <- training_phenotypes[train_l1,"fold"]
      inner_folds <- as.numeric(as.factor(inner_folds)) # so that 1:8 ... 
      ridge_fit <- NULL
      ridge_fit <- try(cv.glmnet(t(expression_filt[selgenes ,train_l1]),# try is here to escape a weird case where dim(nz) fails
                             y = training_phenotypes_indicator[train_l1,disease_set]
                             ,foldid = inner_folds,weights = as.numeric(weights.local)
                             ,type.measure = "mse",alpha = alphaval,family ="multinomial"),silent = T)
      if(!is.null(names(ridge_fit))){
        print(c("fitted ridge",l_lambda))
        p2 <- predict(ridge_fit,newx = t(expression_filt[selgenes,test_l1]) , s = "lambda.min",type = "response")[,,1]
      }
    }
    if(!is.null(names(ridge_fit))){
      fold_prediction[[as.character(l_lambda)]] <- list(lasso = p1,ridge = p2,genes = selgenes, fit = ridge_fit)
    }else{
      p2 <- p1
      p2[] <- 0
      fold_prediction[[as.character(l_lambda)]] <- list(lasso = p1,ridge = p2,genes = selgenes, fit = NULL)
    }
  }
  return(fold_prediction)
}





get_mse_for_lambda_ridgeylasso <- function(l_lambda,fold_predictions,
                                           weight_vector.loc = NULL,training_samples ,mode = "LOW"){ 
  # in is lasso lambda used for gene selection , ridge lambda is taken as min
  training_groups <- names(coef(fold_predictions[[1]]$lassofit))
  lasso.pm <- data.frame(matrix(nrow = length(training_samples), ncol = length(training_groups)))
  rownames(lasso.pm) <- training_samples
  colnames(lasso.pm) <- training_groups
  ridge.pm <- lasso.pm   
  ridge.pm.1se <- lasso.pm   
  ridge.pm.decoupled <- lasso.pm
  failed <- F
  for(foldn in 1:length(fold_predictions)){ # merge prediction matrices over folds
    pm <- fold_predictions[[foldn]][[as.character(l_lambda)]]$lasso
    lasso.pm[rownames(pm),colnames(pm)] <- pm
    pm <- fold_predictions[[foldn]][[as.character(l_lambda)]]$ridge
    ridge.pm[rownames(pm),colnames(pm)] <- pm
    selgenes <- fold_predictions[[foldn]][[as.character(l_lambda)]]$genes
    if(length(selgenes) > 2){
      # 1se for ridge lambda
      test_l1 <- rownames(fold_predictions[[foldn]][[as.character(l_lambda)]]$lasso)
      ridge_fit <- fold_predictions[[foldn]][[as.character(l_lambda)]]$fit
      if(!is.null(ridge_fit)){
        pm <-predict(ridge_fit,newx = t(expression_filt[selgenes,test_l1]) , s = "lambda.1se",type = "response")[,,1] 
        ridge.pm.1se[rownames(pm),colnames(pm)] <- pm
        ## (the lambda here is the ridge lambda)
        # with decoupled probabilities
        pm <- predict_modified(ridge_fit ,test_sams.local = test_l1, lambda = "lambda.min",decouple =T )
        ridge.pm.decoupled[rownames(pm),colnames(pm)] <- pm
      }else{ # if null fit error for any fold completely disregard the entire lambda
        failed <- T
      }
    }
    
  }
  if(failed){
    ridge.pm[] <- 0
    ridge.pm.1se[] <- 0   
    ridge.pm.decoupled[] <- 0 
  }
  
  out <- list()
  if(mode == "LOW"){
    disease_set <- training_groups
    example_weights <- class_weights[training_phenotypes[rownames(lasso.pm),"group"] , "cost+imbalance"] 
  }else if(mode =="TOP"){
    disease_set <- c("bacterial","viral","inflammatory","malaria","TB","KD")
    group <- c()
    for(sam in rownames(lasso.pm)) group <- c(group , disease_set[training_phenotypes_indicator[sam , disease_set] == 1 ])
    example_weights <- class_weights_toplevel[group,"cost+imbalance"]
  }else return(NULL)

  if(!is.null(weight_vector.loc)) example_weights <- weight_vector.loc[rownames(lasso.pm)]
  
  example_weights <-  example_weights/sum(example_weights)
  
  confusion.l <- make_confusion(lasso.pm , training_phenotypes_indicator[rownames(lasso.pm),disease_set] ,mod = "3dconfusion") ##
  square_errs <- get_square_errors(confusion.l$values)
  out$lasso <-  sum(square_errs * example_weights)
  out$lassose <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  
  confusion.r <- make_confusion(ridge.pm , training_phenotypes_indicator[rownames(ridge.pm),disease_set] ,mod = "3dconfusion")
  square_errs <- get_square_errors(confusion.r$values)
  out$ridge <-  sum(square_errs * example_weights)
  out$ridgese <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  
  confusion.rd <- make_confusion(ridge.pm.decoupled , training_phenotypes_indicator[rownames(ridge.pm.decoupled),disease_set] ,mod = "3dconfusion")
  square_errs <- get_square_errors(confusion.rd$values)
  out$ridge_dec <-  sum(square_errs * example_weights)
  out$ridge_decse <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  
  out$partition_ridge <- list()
  conf <- confusion.r$values
  conf[is.na(conf)] <- 0
  for(path in training_groups){
    if(mode == "LOW"){
      samples <- training_samples[training_phenotypes[training_samples,"group"] == path]
      example_weights.loc <- class_weights[training_phenotypes[samples,"group"] , "cost+imbalance"]
      example_weights.loc <-  example_weights.loc/sum(example_weights.loc)
    }else if(mode =="TOP"){
      samples <- rownames(training_phenotypes_indicator)[training_phenotypes_indicator[,path] ==1 ]
      disease_set <- c("bacterial","viral","inflammatory","malaria","TB","KD")
      group <- c()
      for(sam in samples) group <- c(group , disease_set[training_phenotypes_indicator[sam , disease_set] == 1 ])
      example_weights.loc <- class_weights_toplevel[group,"cost+imbalance"]
      example_weights.loc <-  example_weights.loc/sum(example_weights.loc)
    }
    out$partition_ridge[[path]] <- list()
    square_errs <- get_square_errors(confusion.r$values[,,samples])
    out$partition_ridge[[path]]$ridge <-  sum(square_errs * example_weights.loc)
    out$partition_ridge[[path]]$ridgese <- sqrt(var(square_errs) * sum(example_weights.loc^2) / sum(example_weights.loc)^2 )
  }
  confusion.r.1se <- make_confusion(ridge.pm.1se , training_phenotypes_indicator[rownames(ridge.pm.1se),disease_set] ,mod = "3dconfusion")
  square_errs <- get_square_errors(confusion.r.1se$values)
  out$ridge.1se <-  sum(square_errs * example_weights)
  out$ridgese.1se <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  return(out)
}



select_ridge_lambda <- function(ridgeylasso_mse,lasso_lambdas,nse = 1){
  
  mse.l <- c() # lasso
  mse.r <- c() # ridge (min)
  mse.r.1se <- c() # ridge (1se)
  mse.lp <- c() #(+1se)
  mse.rp <- c()
  mse.rp.1se <- c()
  mse.lm <- c()
  mse.rm <- c()
  mse.rm.1se <- c()
  # decoupled probabilities
  mse.rp.dec <- c()
  mse.rm.dec <- c()
  mse.r.dec <- c()
  for(lambda in lasso_lambdas){
    mse.l <- c(mse.l , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$lasso )
    mse.lp <- c(mse.lp, ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$lasso + ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$lassose)
    mse.lm <- c(mse.lm , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$lasso - ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$lassose)
    
    mse.r <- c(mse.r , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge )
    mse.rp <- c(mse.rp, ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge + ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridgese)
    mse.rm <- c(mse.rm , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge - ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridgese)
    
    mse.r.1se <- c(mse.r.1se , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge.1se )
    mse.rp.1se <- c(mse.rp.1se,
                    ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge.1se + ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridgese.1se)
    mse.rm.1se <- c(mse.rm.1se ,
                    ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge.1se - ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridgese.1se)
    
    # decoupled
    MSE <- ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge_dec
    SE <- ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$ridge_decse
    mse.rp.dec <- c(mse.rp.dec ,  MSE + SE )
    mse.rm.dec <- c(mse.rm.dec , MSE - SE )
    mse.r.dec <- c(mse.r.dec , MSE )
  }
  
  pathogen_mses <- list()
  for(path in names(ridgeylasso_mse[[1]]$partition_ridge)){
    mse.m <- c()
    mse.p <- c()
    mse <- c()
    for(lambda in lasso_lambdas){
      mse <- c(mse , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$partition_ridge[[path]]$ridge )
      mse.m <- c(mse.p, ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$partition_ridge[[path]]$ridge +
                   ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$partition_ridge[[path]]$ridgese)
      mse.p <- c(mse.m , ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$partition_ridge[[path]]$ridge -
                   ridgeylasso_mse[[match(lambda,lasso_lambdas)]]$partition_ridge[[path]]$ridgese )
    }
    pathogen_mses[[path]] <- list(mse = mse , mse.m = mse.m , mse.p = mse.p) 
  }

  mse.r[1:3] <- 1.1 # the solution of not predicting anything
  
  ind <- mse.r == min(mse.r)
  se <- mse.rp[ind] -  mse.r[ind]
  mse_cutoff <- mse.r[ind] + ( nse * se )
  # minimise model size given 1se from min mse (over ridge and lasso)
  ridgey_min_lasso_1se <- max(lasso_lambdas[mse.r < mse_cutoff & log(lasso_lambdas) < -1.7]) 
  # the -2 is a fudge to stop the selection of super small erroneous signatures
  # these might be the "predict nothing" trap
  
  out <- list(mse.l = mse.l,
              mse.lp = mse.lp,
              mse.lm = mse.lm,
              mse.r = mse.r,
              mse.rp = mse.rp,
              mse.rm = mse.rm,
              mse.r.1se = mse.r.1se,
              mse.rp.1se = mse.rp.1se,
              mse.rm.1se = mse.rm.1se,
              mse.r.dec = mse.r.dec,
              mse.rp.dec = mse.rp.dec,
              mse.rm.dec = mse.rm.dec,
              selected_lambda = ridgey_min_lasso_1se,
              pathogen_mses)
  return(out)
}

cv_hybrid <- function(outer_fold,prerun = NULL,refit_method = "ridge"){
  training_samples <- all_fits[[outer_fold]]$train_samples
  test_samples <- all_fits[[outer_fold]]$test_samples
  if(any(training_samples%ni%rownames(training_phenotypes))) break
  if(any(test_samples%ni%rownames(training_phenotypes))) break
  lasso_fit_cv <- all_fits[[outer_fold]]$mn_TT
  lambdas <- lasso_fit_cv$lambda
  # cant shortcut need to include the minimum of the ridge cv mse
  # reassign folds according to the selection in vanilla lasso seed run
  for(i in 1:10){ training_phenotypes[all_fits[[i]]$test_samples , "fold"] <- i}
  folds <- training_phenotypes[training_samples,"fold"]
  input_genes <- all_fits[[outer_fold]]$input_genes
  if(refit_method=="ridge"){
    phis <- exp(seq(-8,2,by = 0.1))  
  }else if(refit_method == "lasso"){
    phis <-  exp(seq(-10,-1.5,by = 0.1))
  }
  weights.local <- class_weights[training_phenotypes[training_samples,"group"],"cost+imbalance"]
  tr <- training_samples
  fld <- folds
  inp_gen <- input_genes 
  ll <- lambdas
  rl <- phis
  wt <- weights.local
  if(is.null(prerun)){
    # here is 9 fold cross validation as in case of vanilla lasso but fit for each lambda
    if(refit_method=="ridge"){ # lasso ridge hybrid
    fold_predictions <-  lapply(unique(folds),run_ridges_on_lasso ,
                                  training_samples = tr , folds = fld, input_genes = inp_gen,
                                  lasso_lambdas = ll,ridge_lambdas = rl) 
    }else if(refit_method == "lasso"){ # relaxed lasso
      fold_predictions <-  lapply(unique(folds),run_ridges_on_lasso ,alphaval = 1,
                                  training_samples = tr , folds = fld, input_genes = inp_gen,
                                  lasso_lambdas = ll,ridge_lambdas = rl)
    }
    ridgeylasso_mse <- lapply(ll, get_mse_for_lambda_ridgeylasso , 
                                fold_predictions = fold_predictions ,training_samples = tr)
    
    # select lambda based on the refitted mses
    mse_plot_dat <- select_ridge_lambda(ridgeylasso_mse ,lasso_lambdas = ll)
    selected_lambda <- mse_plot_dat$selected_lambda
  }else{
    print("skip cv")
    # have already run the cv just want to select lambda and fit overall models
    mse_cutoff <- lasso_fit_cv$cvup[lasso_fit_cv$lambda == lasso_fit_cv$lambda.min]
    selected_lambda <- max(lambdas[prerun[[outer_fold]]$mse_plot_dat$mse.r < mse_cutoff &
                                     log(lambdas) < -2])
    mse_plot_dat <- prerun[[outer_fold]]$mse_plot_dat
  }
  plotind <- mse_plot_dat$mse.r > 0
  # refit lasso on whole set
  # get selected genes from lasso
  sel_lambda <- selected_lambda
  lasso_fit <- all_fits[[outer_fold]]$mn_TT
  coefs <- coef(lasso_fit,s=sel_lambda)[[1]]
  selgenes <- rownames(coefs)[as.matrix(coefs)!=0]
  selgenes <- selgenes[selgenes!="(Intercept)"]
  # fit ridge over all data
  inner_folds <- as.numeric(factor(training_phenotypes[training_samples,"fold"]))
  ridge_fit <- cv.glmnet(t(expression_filt[selgenes ,training_samples]),
                         y = training_phenotypes_indicator[training_samples,training_groups]
                         ,foldid = inner_folds,weights = weights.local # ,lambda = phis ## let phi be selected by glmnet !!!!!??
                         ,type.measure = "mse",alpha = 0,family ="multinomial")
  p2 <- predict(ridge_fit,newx = t(expression_filt[selgenes,test_samples]) , 
                s = "lambda.min",type = "response")[,,1]
  out <- list(
    sel_lambda = sel_lambda ,  
    selected_ridge = ridge_fit,
    selected_genes = selgenes,
    ridge_pred = p2,
    mse_plot_dat = mse_plot_dat
  )
  return(out)
}



plot_ridge_cv <- function(mse_pltdat,lasso_fit_cv = NULL,aux_mse_pltdat = NULL,mode = "LOW",
                          path = F,
                          ridge_lambdas ,xrange = NULL, colvec = c("blue","red","blue","purple"),cexval = 1,
                          xvar ="lambda",add_ = F,lw = 1,capval = 0.001){
  ind <- round(lasso_fit_cv$cvm,digits = 8) %in% round(mse_pltdat$mse.l,digits = 8) # round to get around floating point errs 
  modelsizes <- as.numeric(lasso_fit_cv$nzero[ind])
  mod_pltdat <- cbind(modelsizes,
                      log(lasso_fit_cv$lambda[ind]) 
                    ,mse_pltdat$mse.l ,mse_pltdat$mse.lp,mse_pltdat$mse.lm 
                    ,mse_pltdat$mse.r ,mse_pltdat$mse.rp,mse_pltdat$mse.rm )
  colnames(mod_pltdat) <- c("ms","lambda","l","lp","lm","r","rp","rm")
  # vanilla lasso
  lasso_min <- mod_pltdat[mod_pltdat[,"l"] == min(mod_pltdat[,"l"]) ,]
  lasso_1se <- mod_pltdat[mod_pltdat[,"l"] < lasso_min["lp"],][1,]
  se <- lasso_min["lp"] - lasso_min["l"]
  lasso_2se <- mod_pltdat[mod_pltdat[,"l"] < (lasso_min["l"] + (2*se)),][1,]
  if(add_){
    errbar(mod_pltdat[,"l"], x = mod_pltdat[,xvar],cex =cexval,yplus = mod_pltdat[,"lp"] , yminus = mod_pltdat[,"lm"],
           col = colvec[1] ,add = T,cap =capval)
  }else{
    pretty_plot_area(cols = c("grey95","white"),text_col = "grey30",
                     ytick = c(0,0.8,0.1),x_lab = xvar,y_lab = "MWSE",
                     xtick = c(min(xrange),max(xrange),0.5), margins = c(5,5,2,2))
    
    points( mod_pltdat[,xvar], mod_pltdat[,"l"] ,col = colvec[1],cex =cexval,pch =16)
    points( mod_pltdat[,xvar], mod_pltdat[,"lp"] ,type="l",col = colvec[1],cex =cexval,lwd = 0.7)
    points( mod_pltdat[,xvar], mod_pltdat[,"lm"] ,type="l",col = colvec[1],cex =cexval,lwd = 0.7)
  }
  points(lasso_min[xvar],lasso_min["l"],col = "black",pch = 16,cex = cexval)
  abline(v = lasso_min[xvar],lty =1,col =colvec[1],lwd = lw) # min vert 
  text(lasso_min[xvar] , 0.7, labels = lasso_min["ms"],col =colvec[1])
  
  abline(h = lasso_min["lp"],lty =1,col =add.alpha(colvec[1],alpha = 0.5),lwd = lw) # min +se horiz
  points(lasso_1se[xvar],lasso_1se["l"],col = "black",pch = 16)
  abline(v = lasso_1se[xvar],lty =3,col = colvec[1],lwd = lw) # 1se vert 
  text(lasso_1se[xvar] , 0.7, labels = lasso_1se["ms"],col =colvec[1])
  
  abline(h = (lasso_min["l"] + (2*se)),lty =2,col =add.alpha(colvec[1],alpha = 0.5),lwd = lw) # min +se horiz
  abline(v = lasso_2se[xvar],lty =3,col = colvec[1],lwd = lw) # 1se vert 
  text(lasso_2se[xvar] , 0.7, labels = lasso_2se["ms"],col =colvec[1])
  
  if(is.null(aux_mse_pltdat)){
    mse_pltdat <- list(mse_pltdat)
  }else{
    mse_pltdat <- c(list(mse_pltdat),aux_mse_pltdat)
  }
  
  for(i in 1:length(mse_pltdat)){
    mse_pltdat_obj <- mse_pltdat[[i]]
    mod_pltdat <- cbind(modelsizes,
                        log(lasso_fit_cv$lambda[ind])
                        ,mse_pltdat_obj$mse.r ,mse_pltdat_obj$mse.rp,mse_pltdat_obj$mse.rm )
    colnames(mod_pltdat) <- c("ms","lambda","r","rp","rm")
    # hybrid lasso methods
    hybrid_min <- mod_pltdat[mod_pltdat[,"r"] == min(mod_pltdat[,"r"]) ,]
    hybrid_1se <- mod_pltdat[mod_pltdat[,"r"] < hybrid_min["rp"],][1,]
    se <- hybrid_min["rp"] - hybrid_min["r"]
    hybrid_2se <- mod_pltdat[mod_pltdat[,"r"] < (hybrid_min["r"] + (2* se)),][1,]
    
    points( mod_pltdat[,xvar], mod_pltdat[,"r"] ,col = colvec[i+1],cex =cexval,pch =16)
    points( mod_pltdat[,xvar], mod_pltdat[,"rp"] ,type="l",col = colvec[i+1],cex =cexval,lwd = 0.7)
    points( mod_pltdat[,xvar], mod_pltdat[,"rm"] ,type="l",col = colvec[i+1],cex =cexval,lwd = 0.7)
    
    points(hybrid_min[xvar],hybrid_min["r"],col = "black",pch = 16,cex = cexval)
    abline(v = hybrid_min[xvar],lty =1 ,col =colvec[i+1],lwd = lw)
    text(hybrid_min[xvar] , 0.7, labels = hybrid_min["ms"],col =colvec[i + 1])
    
    abline(h = hybrid_min["rp"],lty = 1 ,col = add.alpha(colvec[i+1],alpha = 0.5),lwd = lw)
    points(hybrid_1se[xvar],hybrid_1se["r"],col = "black",pch = 16)
    abline(v = hybrid_1se[xvar],lty =3,col =colvec[i+1],lwd = lw)
    text(hybrid_1se[xvar] , 0.7, labels = hybrid_1se["ms"],col =colvec[i + 1])
    points(hybrid_2se[xvar],hybrid_2se["r"],col = "black",pch = 16)
    abline(v = hybrid_2se[xvar],lty =3,col =colvec[i+1],lwd = lw)
    text(hybrid_2se[xvar] , 0.7, labels = hybrid_2se["ms"],col =colvec[i + 1])
    
  }
}








.run_hybrid_fold <- function(foldn,relax_alpha = 0 , lambdas , phis ){
  weights <-   class_weights[training_phenotypes[training_samples_all,"group"] , "cost+imbalance"]
  
  weights.loc <- weights[training_phenotypes[training_samples_all , "fold"] != foldn]
  training_samples.loc <- training_samples_all[training_phenotypes[training_samples_all , "fold"] != foldn]
  test_samples.loc <- training_samples_all[training_phenotypes[training_samples_all , "fold"] == foldn]
  
  
  lasso_fit <- glmnet(t(expression_filt[training_genes , training_samples.loc]),
                      training_phenotypes[training_samples.loc , "group"],
                      weights = weights.loc,   family="multinomial", type.multinomial = "grouped",
                      alpha=1 , standardize = TRUE ,lambda = lambdas)
  print("lasso fit")
  coef_mat <- as.matrix(coef(lasso_fit)[[1]])
  colnames(coef_mat) <- lambdas
  coef_mat <- coef_mat[-1,] #remove intercept
  M_lambda_list <- lapply(lambdas , function(l)  rownames(coef_mat)[abs(coef_mat[,as.character(l)]) > 0])
  names(M_lambda_list) <- as.character(lambdas)

  .fit_single_ridge <- function(M_lambda){
    if(length(M_lambda) > 1){
      ridge_fit <- glmnet(t(expression_filt[M_lambda , training_samples.loc]),
                          training_phenotypes[training_samples.loc , "group"],
                          weights = weights.loc,family="multinomial", type.multinomial = "grouped",
                          alpha= relax_alpha , standardize = TRUE ,lambda = phis)
      return(ridge_fit)
    }else{
      return(NULL)
    }
  }
  
  ridge_fits <- mclapply(M_lambda_list , .fit_single_ridge ,mc.cores = 10)
  names(ridge_fits) < lambdas
  print("ridge fit")
  .get_predmats <- function(i){
    ridge_fit<- ridge_fits[[i]]
    M_lambda <- M_lambda_list[[i]]
    predmats <- list()
    for(lambda in phis){
      if(is.null(ridge_fit)){
        return(NULL)
      }else{
        pm <- predict(ridge_fit , s = lambda,newx = t(expression_filt[M_lambda , test_samples.loc])
                      ,type = "response")[,,1]
      }
      predmats[[as.character(lambda)]] <- pm  
    }
    return(predmats)
  }
  
  .get_predmats_lasso <- function(i){
    lambda <- lambdas[i]
    pm <- predict(lasso_fit , s = lambda,newx = t(expression_filt[training_genes, test_samples.loc])
                  ,type = "response")[,,1]
    return(pm)
  }
  i <- 1
  while(is.null(ridge_fits[[i]])){
    i <- i + 1# skip the first few as null
    
  } 
  tmp <- lapply(1:length(lambdas) , .get_predmats_lasso ) # the lasso prediction matrices
  predmat_list <- mclapply(i:length(lambdas) , .get_predmats , mc.cores = 10)
  print("predictions made")
  return(list(ridge = predmat_list,lasso = tmp,fits_ridge  = ridge_fits,m_l = M_lambda_list ,
              lambdas_new = lambdas[i:length(lambdas)]))
}


.get_hybrid_mse <- function(invar ,predictions,lambdas , phis){
  i <- invar[1]
  j <- invar[2]
  out <- list()
  full_pm <- training_phenotypes_indicator[training_samples_all , ]
  full_pm <- full_pm[,!colnames(full_pm) %in% c("bacterial" , "viral","inflammatory")]
  for(foldn in 1:10){
    if(is.na(j)){
      pm <- predictions[[foldn]]$lasso[[i]]
    }else{
      pm <- predictions[[foldn]]$ridge[[i]][[j]]  
    }
    full_pm[rownames(pm),colnames(pm)] <- pm
  }
  example_weights <- class_weights[training_phenotypes[training_samples_all,"group"] , "cost+imbalance"] 
  example_weights <-  example_weights/sum(example_weights)
  confusion <- make_confusion(full_pm , training_phenotypes_indicator[rownames(full_pm),colnames(full_pm)] ,mod = "3dconfusion") ##
  square_errs <- get_square_errors(confusion$values)
  out$mse <-  sum(square_errs * example_weights)
  out$se <- sqrt(var(square_errs) * sum(example_weights^2) / sum(example_weights)^2 )
  
  out$lasso_lambda <- lambdas[i]
  if(!is.na(j)){
    out$ridge_lambda <- phis[j]
  }
  out
}


two_stage_bridge <- function(lambdas , phis,relax_alpha,label = "0"){
  predictions_cv_lasso_ridge <- list()
  for(i in 1:10){
    print(i)
    predictions_cv_lasso_ridge[[i]]  <- .run_hybrid_fold(i,relax_alpha,lambdas , phis) # unsuprisingly is a big ol' object
  }
  # makes a simple list of the lambdas at whcih to look at mse
  ridge_lasso_hybrid_mse <- list()
  for(i in 1:length(predictions_cv_lasso_ridge[[1]]$ridge) ){
    lambda <- lambdas[i]
    for(j in 1:length(predictions_cv_lasso_ridge[[1]]$ridge[[1]]) ){
      ridge_lasso_hybrid_mse[[length(ridge_lasso_hybrid_mse) + 1]] <-   c(i,j)
    }
  }
  
  new_lambdas <- predictions_cv_lasso_ridge[[1]]$lambdas_new # to ignore those where no gene
  
  ridge_lasso_hybrid_mse.run <- mclapply(ridge_lasso_hybrid_mse ,
                                         .get_hybrid_mse ,mc.cores = 5, 
                                         predictions = predictions_cv_lasso_ridge,
                                         lambdas= new_lambdas, 
                                         phis=phis)
  
  lasso_mse.run <- mclapply(1:length(new_lambdas) ,
                            .get_hybrid_mse ,mc.cores = 5, 
                            predictions = predictions_cv_lasso_ridge,lambdas=new_lambdas , phis=phis)
  
  ridge_lasso_hybrid_mse_mat <- list(mse = data.frame() , se = data.frame() )
  for(pnt in ridge_lasso_hybrid_mse.run){
    ridge_lasso_hybrid_mse_mat$mse[as.character(pnt$lasso_lambda),
                                   as.character(pnt$ridge_lambda)] <- pnt$mse
    ridge_lasso_hybrid_mse_mat$se[as.character(pnt$lasso_lambda),
                                  as.character(pnt$ridge_lambda)] <- pnt$se
  }
  return(list(pred_mats = predictions_cv_lasso_ridge , 
              mse = ridge_lasso_hybrid_mse.run , 
              mse.lasso = lasso_mse.run ,
              new_lambdas=new_lambdas,
              mse_mat = ridge_lasso_hybrid_mse_mat))
}


plot_cv_twostage_bridge <- function(fit_obj,type = "contour", yaxl ="lambda"){
  
  mse_mat <- fit_obj$mse_mat$mse
  se_mat <- fit_obj$mse_mat$se
  mse.lasso <- fit_obj$mse.lasso
  
  # find ridge minima
  min_points <- as.data.frame(matrix(nrow = length(rownames(mse_mat)),ncol = 4))
  rownames(min_points) <- rownames(mse_mat)
  colnames(min_points) <- c("x","y","z","se")
  for(lasso_lambda in rownames(mse_mat)){
    r_l <- colnames(mse_mat)[mse_mat[lasso_lambda,] == min(mse_mat[lasso_lambda,]) ]
    r_l <- as.character(max(as.numeric(r_l)))
    if(length(r_l)>1){break}
    min_points[lasso_lambda,] <- c(as.numeric(lasso_lambda), as.numeric(r_l) ,
                                   mse_mat[lasso_lambda , r_l],  se_mat[lasso_lambda , r_l])
    
  }
  
  # cv plot
  mse_mat <- mse_mat[!is.na(mse_mat[,1]),]
  mse.lasso.2 <- matrix(unlist(mse.lasso),ncol =3,byrow = T)
  rownames(mse.lasso.2) <- as.character(mse.lasso.2[,3])
  mse_mat <- mse_mat[order(as.numeric(rownames(mse_mat)),decreasing = F)
                     ,order(as.numeric(colnames(mse_mat)),decreasing = F)]
  
  if(type == "contour"){ # for a contour plot 
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,xlim = c(-6,-2),
            y = log(as.numeric(colnames(mse_mat))),xlab = "lambda (LASSO)",ylab = yaxl, levels = seq(0,0.8,by = 0.01)) 
    
    # for the 2d tuning
    min_lambda <- min_points[min_points[,"z"] == min(min_points[,"z"]),1]
    se_cutoff <- min_points[as.character(min_lambda),"z"] + min_points[as.character(min_lambda),"se"]
    se2_cutoff <- min_points[as.character(min_lambda),"z"] + (2 * min_points[as.character(min_lambda),"se"])
    # min
    points(log(min_points[as.character(min_lambda),"x"]) , log(min_points[as.character(min_lambda),"y"]),
           col = "black", cex = 1.5, pch = 16)
    # se line
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,
            y = log(as.numeric(colnames(mse_mat))),levels = round(se_cutoff,digits = 3),add = T,col ="red")
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,
            y = log(as.numeric(colnames(mse_mat))),levels = round(se2_cutoff, digits = 3),add = T,col ="blue")
    
    
  }else if(type == "heat"){ # heatmap inc lasso
    
    col_fun = colorRamp2(c(log(0.1),log(0.4)), c("red","black") )
    
    
    lambdas <- log(as.numeric(rownames(mse_mat)))
    phis <- log(as.numeric(colnames(mse_mat)))
    
    plot(0,0,cex = 0 , xlim = c(min(lambdas),max(lambdas)),ylim = c(min(phis)-2,max(phis)),
         xlab = "lambda",ylab ="phi")
    lambda_jump <- lambdas[2] - lambdas[1]
    phi_jump <- phis[2] - phis[1]
    for( i in 1:ncol(mse_mat)){
      for(j in 1:nrow(mse_mat)){
        val <- mse_mat[j,i]
        rect(xleft = lambdas[j] - lambda_jump / 2,
             xright = lambdas[j] + lambda_jump / 2
             ,ybottom = phis[i] - phi_jump /2
             ,ytop = phis[i] + phi_jump /2,col = col_fun(log(val)),border =  col_fun(log(val))) 
      }
    }  
    # vanilla lasso -- ignore / remove the -10 and put label on y
    lambdas <- log(as.numeric(rownames(mse.lasso.2)))
    lambda_jump <- lambdas[2] - lambdas[1]
    for(j in 1:nrow(mse_mat)){
      val <- mse.lasso.2[j,1]
      rect(ybottom = min(phis)-2 , ytop = min(phis)-1
           ,xleft = lambdas[j] - lambda_jump / 2,
           xright = lambdas[j] + lambda_jump / 2, col = col_fun(log(val)) ,border =  col_fun(log(val)) )
    }
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,
            y = log(as.numeric(colnames(mse_mat))),add = T,levels = seq(0,0.8,by = 0.01)) 
    points(log(min_points[,"x"]) , log(min_points[,"y"]))
    
    # for the 2d tuning
    min_lambda <- min_points[min_points[,"z"] == min(min_points[,"z"]),1]
    se_cutoff <- min_points[as.character(min_lambda),"z"] + min_points[as.character(min_lambda),"se"]
    # min
    points(log(min_points[as.character(min_lambda),"x"]) , log(min_points[as.character(min_lambda),"y"]),col = "red")
    # se line
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,
            y = log(as.numeric(colnames(mse_mat))),levels = se_cutoff,add = T,col ="red")
    contour(z = as.matrix(mse_mat) ,x = log(as.numeric(rownames(mse_mat))) ,
            y = log(as.numeric(colnames(mse_mat))),levels = 0.25,add = T,col ="blue")
  }
}








