

load("dataset/mega/base_microarray.RData")
load("dataset/status.RData")
dim(e.set)

library("BiocParallel") 
BPPARAM = MulticoreParam(workers=8)
library(as.color)
library(lumi)
library(limma)
library(ggplot2)
library(ggfortify)
library(beeswarm)
library("sva")
library(GEOquery)
source("scripts/functions.R")
library(factoextra)
library(stringr)
library(COCONUT)
library(xml2)
# library(superheat) 

plotdir <- "plots/"

# GEOpatch

getAndParseGSEMatrices=function (GEO, destdir, AnnotGPL, getGPL = TRUE) 
{
  GEO <- toupper(GEO)
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
  b = getDirListing(sprintf(gdsurl, stub, GEO))
  b=b[-1]
  message(sprintf("Found %d file(s)", length(b)))
  ret <- list()
  for (i in 1:length(b)) {
    message(b[i])
    destfile = file.path(destdir, b[i])
    if (file.exists(destfile)) {
      message(sprintf("Using locally cached version: %s", 
                      destfile))
    }
    else {
      print(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                    stub, GEO, b[i]))
      download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                            stub, GEO, b[i]), destfile = destfile, mode = "wb", 
                    method = getOption("download.file.method.GEOquery"))
    }
    ret[[b[i]]] <- parseGSEMatrix(destfile, destdir = destdir, 
                                  AnnotGPL = AnnotGPL, getGPL = getGPL)$eset
  }
  return(ret)
}

environment(getAndParseGSEMatrices) <- asNamespace("GEOquery")
assignInNamespace("getAndParseGSEMatrices", getAndParseGSEMatrices, ns="GEOquery")

# Load data from GEO ####

get_dataset <- function(expt_id,log = T,local=F,M=F){
  if(local==T){
    eset_tmp <-  getGEO(filename = expt_id)
  }else{
    eset_tmp <-  getGEO(expt_id)
  }
  get_subbatch <- function(eset_tmp){
    # filter any rows with NA 
    rmv <- which(apply(exprs(eset_tmp), 1, function(x) any ( is.na(x))))
    if(length(rmv) >0 ){
      eset_tmp <- eset_tmp[-rmv,]  
    }
    png(paste(plotdir,expt_id,"_box_prenorm.png",sep=""),height = 500,width = 5000)
    boxplot(exprs(eset_tmp))
    dev.off()
    if(log==T){
      eset_tmp.N<-lumiN(lumiQ(lumiT(eset_tmp, method='log')),method='rsn') 
    }else{
      eset_tmp.N<-lumiN(lumiQ(eset_tmp),method='rsn') 
    }
    png(paste(plotdir,expt_id,"_box_postnorm.png",sep=""),height = 500,width = 5000)
    boxplot(eset_tmp.N)
    dev.off()
    return( eset_tmp.N)
  }
  print(length(eset_tmp))
  if(M){
    return(get_subbatch(eset_tmp )  )
  }else{
    return(lapply(eset_tmp , get_subbatch)  )  
  }
}

## ramilo viral ###################################################

expt_id <- "GSE38900"

file1 <- "dataset/external_datasets/GSE38900/GSE38900-GPL10558_series_matrix.txt"
file2 <- "dataset/external_datasets/GSE38900/GSE38900-GPL6884_series_matrix.txt"


tab1 <- read.table(file1,sep="\t",quote ="\"" , stringsAsFactors = F,skip=73,header=T,fill=NA)
tab2 <- read.table(file2,sep="\t",quote ="" , stringsAsFactors = F,skip=73,header=T,fill=NA)
#"nas are from single line"

rownames(tab1) <- tab1[,1]
tab1 <- tab1[,2:ncol(tab1)]
tab1 <- tab1[rowSums(is.na(tab1))==0,]

rownames(tab2) <- tab2[,1]
tab2 <- tab2[,2:ncol(tab2)]
tab2 <- tab2[rowSums(is.na(tab2))==0,]


pheno1 <- read.table(file1,sep="\t",quote ="\"" , stringsAsFactors = F,skip=grep("Sample_title",readLines(file1)) 
                     ,nrows = 29,header=T,fill=NA)
rownames(pheno1) <- make.names(pheno1[,1],unique =T)
pheno1 <- pheno1[,-1]
pheno1 <- as.data.frame(t(pheno1))
pheno1 <- new("AnnotatedDataFrame" , data =   pheno1 ) #, varMetadata = meta)


pheno2 <- read.table(file2,sep="\t",quote ="" , stringsAsFactors = F,skip=grep("Sample_title",readLines(file2)) 
                     ,nrows = 29,header=T,fill=NA)
rownames(pheno2) <- make.names(pheno2[,1],unique =T)
pheno2 <- pheno2[,-1]
pheno2 <- as.data.frame(t(pheno2))
pheno2 <- new("AnnotatedDataFrame" , data =   pheno2 ) #, varMetadata = meta)

eset1 <- ExpressionSet(assayData = as.matrix(tab1) , phenoData = pheno1)
png(paste(plotdir,expt_id,"_box_prenorm_b1.png",sep=""),height = 500,width = 5000)
boxplot(eset1)
dev.off()

eset1.N<-lumiN(lumiQ(lumiT(eset1, method='log')),method='rsn') 
png(paste(plotdir,expt_id,"_box_postnorm_b1.png",sep=""),height = 500,width = 5000)
boxplot(eset1.N)
dev.off()


eset2 <- ExpressionSet(assayData = as.matrix(tab2) , phenoData = pheno2)
png(paste(plotdir,expt_id,"_box_prenorm_b2.png",sep=""),height = 500,width = 5000)
boxplot(eset2)
dev.off()

eset2.N<-lumiN(lumiQ(lumiT(eset2, method='log')),method='rsn') 
png(paste(plotdir,expt_id,"_box_postnorm_b2.png",sep=""),height = 500,width = 5000)
boxplot(eset2.N)
dev.off()
# second has many failed probes (all exactly 1 in input)


## read iris from file ####
filename<-"dataset/iris/IRIS\ Apr-Nov\ 2010\ background\ subtraction_Sample_Probe_Profile_150413.txt"
object<-lumiR(filename, detectionTh = 0.01, checkDupId = TRUE, QC = TRUE, columnNameGrepPattern =
                list('AVG_SIGNAL', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS'), verbose = TRUE)

lumi.x.N <- lumiN(lumiQ(lumiT(object , method='log')) ,method='rsn') 
GSE72810  <- lumi.x.N


## download from GEO ####
options('download.file.method'='curl')

#GSE97741 rsv no controls just convalescent
#GSE25504  neonatal sepsis (#4)

GSE68004 <- get_dataset("GSE68004") # bv
GSE22098 <- get_dataset("GSE22098") # sle
GSE40396 <- get_dataset("GSE40396" ,log = F) # bv
GSE65391 <- get_dataset("GSE65391",log = F)# sle  # 158 patients multiple timepoints
GSE38900 <- get_dataset("GSE38900") # ramilo viral 
GSE29366 <- get_dataset("GSE29366") # flu
GSE63881 <- get_dataset("GSE63881",log = F)   # kd
GSE72810 <- get_dataset("GSE72810",log=F) # iris validation set
GSE30119 <- get_dataset("GSE30119") # ramilo staph aureus
GSE39941 <- get_dataset("GSE39941") # anderson tb # no healthy controls  -- use ltbi
GSE42026 <- get_dataset("GSE42026")
GSE64456 <- get_dataset("GSE64456")
GSE34404 <- get_dataset("GSE34404")


dsdir <- "dataset/external_datasets/"


datasets <- list( GSE29366[[1]] ,  # flu
                  GSE38900[[1]], # rsv
                  GSE38900[[2]], # rsv + rhino
                  GSE68004 , # kd gas adenovirus
                  GSE22098  , # sle staph strep jia/still
                  GSE40396 , # spread
                  GSE65391  # longditudinal sle
                  #     ,GSE63881[[1]] # kd no controls
                  ,GSE72810 # iris
                  ,GSE30119[[1]] # ramil staph
                  ,GSE25504[[4]] #neonate sepsis
                  ,GSE39941[[1]] # anderson tb
                  ,GSE42026[[1]] #l. Transcriptomic profiling in childhood H1N1/09 influenza reveals reduced expression of protein synthesis genes. 
                  ,GSE64456[[1]] # majahan neonates
                  ,GSE34404[[1]] # malaria
                  
)


names(datasets) <- c( "GSE29366",  # flu
                      "GSE38900[[1]]", # rsv
                      "GSE38900[[2]]", # rsv + rhino
                      "GSE68004" , # kd gas adenovirus
                      "GSE22098"  , # sle staph strep jia/still
                      "GSE40396" , # spread
                      "GSE65391"  # longditudinal sle
                      #  ,"GSE63881" # kd no controls
                      ,"GSE72810" # iris
                      ,"GSE30119" # ramil staph
                      ,"GSE25504[[4]]" #neonate sepsis
                      ,"GSE39941"
                      ,"GSE42026" 
                      ,"GSE64456"
                      ,"GSE34404"
)

# not included
# GSE40586 # not paediatric
#GSE66099 # affy
#GSE11755
#EMEXP3567
#E-GEOD-11907 - Transcription profiling of human stemic juvenile idiopathic arthritis, , systemic lupus erythematosus,, type I diabetes, metastatic melanoma, acute infections, or liver-transplant recipients undergoing immunosuppressive therapy (n = 37)


# combine datasets into one matrix ####

mega_expr <- e.set
rownames(mega_expr) <- unlist(lapply(rownames(mega_expr),function(i) strsplit(i,split = "_")[[1]][2]))
illumina_to_nuid <- cbind(IlluminaID2nuID(rownames(mega_expr),species = "Human")[,"nuID"] ,  rownames(mega_expr))
save(illumina_to_nuid , file = "data_objects/illumina_to_nuid.RData")
rownames(mega_expr) <- IlluminaID2nuID(rownames(mega_expr),species = "Human")[,"nuID"] # nuids not probe ids
present_ids <- list()
present_ids[[1]] <- rownames(mega_expr)
# filter samples  -- some adults 
colours <- rep(length(datasets)+1,ncol(mega_expr))
merged_expr <- mega_expr
for(i in 1:length(datasets) ){
  print("------------------")
  print(names(datasets)[i])
  if(names(datasets)[i] %in% c("GSE25504[[4]]","GSE38900[[2]]")){ # excluding datasets to ensure missing probes arent excluded unnecessarily
    print("EXCLUDE")
  }else{
    array <- datasets[[i]]
    print(c(dim(exprs(array)),dim(pData(array))))
    probes <- rownames(exprs(array))
    nuids <- IlluminaID2nuID(probes , species = "Human")[,"nuID"]
    present_ids[[i+1]] <- nuids
    exprs_tmp <- exprs(array)
    rownames(exprs_tmp) <- nuids
    shared <- intersect(rownames(merged_expr) , rownames(exprs_tmp))
    print(length(shared))
    print(dim(merged_expr[shared , ]))
    merged_expr <- cbind(merged_expr[shared , ],exprs_tmp[shared , ])
    colours <- c(colours , rep(i , ncol(exprs_tmp)))  
    print(dim(exprs_tmp[shared , ]))
    print(dim(merged_expr[shared , ]))  
  }
  
}


# merge phenotype data ####


mega_link <- read.table("dataset/sample_details.txt",sep = "\t",head = T,quote="",fill=NA,stringsAsFactors = F)

mega_linkertab <- cbind(mega_link[,8],mega_link[,"site.number.episode"])
rownames(mega_linkertab) <- mega_link[,8]

mega_phenos2 <- read.table("dataset/mega/mega_pheno.tab",stringsAsFactors = F, sep = "\t",
                           quote ="",fill=NA,skip=1,header = T)


merged_phenos <- data.frame(samples = c() ,age=c(),gender = c(), ethnicity = c(),phenotype=c())
samples_to_exclude <- c()
samples_complex <-  c()

num <- "GSE29366"
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
age_tmp <- gsub("age: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
age_tmp[grep("y",age_tmp)] <- as.numeric(gsub("y","",age_tmp[grep("y",age_tmp)])) *12
age_tmp <- as.numeric(gsub("[a-z]","",age_tmp))
merged_phenos[samples,"age"] <-  age_tmp
merged_phenos[samples,"gender"] <- gsub("gender: ","",pData(datasets[[num]])[,"characteristics_ch1.3"])
merged_phenos[samples,"ethnicity"] <- gsub("race: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
phenos <- as.character(pData(datasets[[num]])[,"characteristics_ch1"])
phenos[grep("FLU",phenos)] <- "flu"
phenos[grep("control",phenos)] <- "HC"
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- 1
merged_phenos[samples,"experiment"] <- num


num <- "GSE38900[[1]]"
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
# all cildren below 2 years so months
merged_phenos[samples,"age"] <-  as.numeric(gsub("age: ","",pData(datasets[[num]])[,"characteristics_ch1"]))
merged_phenos[samples,"gender"] <- gsub("gender: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
merged_phenos[samples,"ethnicity"] <- "?" #gsub("race: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub("","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos[grep("health",phenos)] <- "HC"
phenos[grep("RSV",phenos)] <- "RSV"
merged_phenos[samples,"phenotype"] <- phenos 
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num
#merged_phenos[samples,] 

num <- "GSE38900[[2]]" #  not used -- serious problems with the data 
#  many probes are set to zero / a single value 
# excluding these probes would reduce scope in the other datasets
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub("age [(]months[)]: ","",as.character(pData(datasets[[num]])[,"characteristics_ch1"])))
merged_phenos[samples,"gender"] <- gsub("gender: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
merged_phenos[samples,"ethnicity"] <- "?" #gsub("race: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub("","",pData(datasets[[num]])[,"characteristics_ch1.2"])
samples_to_exclude <- c(samples_to_exclude ,samples[gsub("","",pData(datasets[[num]])[,"characteristics_ch1.2"]) ==  "sample group: 1-2 month after the acute hospitalization in children with RSV"])
phenos[grep("acute human Influenza A",phenos)] <- "flu"
phenos[grep("RSV",phenos)] <- "RSV"
phenos[grep("health",phenos)] <- "HC"
phenos[grep("inovirus",phenos)] <- "rhino"
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num

num <- "GSE68004"
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub("age .mos..: ","",pData(datasets[[num]])[,"characteristics_ch1"]))
merged_phenos[samples,"gender"] <- gsub("gender: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
merged_phenos[samples,"ethnicity"] <- gsub("race: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.3"])
phenos[grep("health",phenos)] <- "HC"
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num
# inKD == incomplete KD ,, cKD  == complete

num <- "GSE22098"
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1"])) * 12
merged_phenos[samples,"gender"] <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
merged_phenos[samples,"ethnicity"] <-  gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub("","",pData(datasets[[num]])[,"characteristics_ch1.3"])
phenos[grep("health",phenos)] <- "HC"
phenos <- gsub("illness: ","",phenos)
phenos[grep("till",phenos)] <- "JIA" #most of stills disease cases are too old
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num


num <- "GSE40396"
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.1"]))
merged_phenos[samples,"gender"] <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1"])
merged_phenos[samples,"ethnicity"] <-  gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <-gsub("pathogen: ","",pData(datasets[[num]])[,"characteristics_ch1.4"])
phenos[grep("None",phenos)] <- "HC"
phenos[grep("Adenovirus",phenos)] <- "adeno"
phenos[grep("Rhinovirus",phenos)] <- "rhino"
merged_phenos[samples,"phenotype"] <-  phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num


num <- "GSE65391"
samples <- rownames(pData(datasets[[num]]))
samples <- rownames(pData(datasets[[num]]))[pData(datasets[[num]])[samples,"characteristics_ch1.3"]=="visit: 1"]
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.13"])) * 12
merged_phenos[samples,"gender"] <- gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.11"])
merged_phenos[samples,"ethnicity"] <-  gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.12"])
phenos <-  gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.10"])
phenos[grep("Health",phenos)] <- "HC"
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + as.numeric(gsub("batch: ","",pData(datasets[[num]])[samples,10]))
merged_phenos[samples,"experiment"] <- num
# col 11 batch replicate  


num <- "GSE72810"
# IRIS 
pheno_file <- "dataset/iris/patients_in_the_paper.tab"
iris_data1 <- read.table(pheno_file,stringsAsFactors = F,sep="\t",header = T,quote = "")
samples <- rownames(pData(datasets[[num]]))
rownames(iris_data1) <- iris_data1[,1]

samples <- samples[samples%in%rownames(iris_data1)]

merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  iris_data1[samples,"Age"] *12
merged_phenos[samples,"gender"] <- iris_data1[samples ,"Sex"]
merged_phenos[samples,"ethnicity"] <- iris_data1[samples,"Ethnicity"]
merged_phenos[samples,"phenotype"] <-  iris_data1[samples,15]

merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm =T) +1
merged_phenos[samples,"experiment"] <- num
merged_phenos[samples,][merged_phenos[samples,"phenotype"] == "R","phenotype"]  <- "RSV"
merged_phenos[samples,][merged_phenos[samples,"phenotype"] == "H","phenotype"]  <- "flu"
merged_phenos[samples,][merged_phenos[samples,"phenotype"] == "H R","phenotype"]  <- "flu RSV"
merged_phenos[samples,][merged_phenos[samples,"phenotype"] == "C","phenotype"]  <- "HC"

bact <- samples[merged_phenos[samples,"phenotype"] == "bact"]
oth <- samples[merged_phenos[samples,"phenotype"] %in% c("H other","H other ","H other bact","R bact","R other bact","R other")]
oth2 <- samples[merged_phenos[samples,"phenotype"] %in% c("U","other")]

merged_phenos["99_11ep1_1U_1","phenotype"] <- "gas"
merged_phenos["102_23ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["117_91ep1_1P_1","phenotype"] <- "gas piv" # parainfluenza
merged_phenos["120_116ep1_1M_1","phenotype"] <- "pneumo mpv" # metapneumo
merged_phenos["145_111ep1_1U_1","phenotype"] <- "gas"
merged_phenos["164_163ep1_1M_1","phenotype"] <- "pneumo mpv"
merged_phenos["173_176ep1_1A_1","phenotype"] <- "adeno staph pseudomonas"
merged_phenos["177_188ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["179_288ep1_1U_1","phenotype"] <- "saureus"
merged_phenos["180_290ep1_1F_1","phenotype"] <- "flu pneumo"
merged_phenos["191_291ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["192_292ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["216_301ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["226_293ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["228_302ep1_1U_1","phenotype"] <- "gas cdiff"
merged_phenos["236_222ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["238_294ep1_1U_1","phenotype"] <- "pneumo"
merged_phenos["239_295ep1_1U_1","phenotype"] <- "pneumo"


merged_phenos[oth,"phenotype"]  <- gsub("/","" ,iris_data1[oth,c("combined.Virology")] )
samples_to_exclude <- c(samples_to_exclude,samples[merged_phenos[samples,"phenotype"] == "excl"])
merged_phenos["21_122ep1_1RBMN_1","phenotype"] <- paste(merged_phenos["21_122ep1_1RBMN_1","phenotype"] , "pneumo")
merged_phenos["53_61ep1_1HPN_1","phenotype"] <- paste(merged_phenos["53_61ep1_1HPN_1","phenotype"] , "pneumo")
merged_phenos["6_110ep1_1RA_1","phenotype"] <- paste(merged_phenos["6_110ep1_1RA_1","phenotype"] , "pneumo")

merged_phenos["81_101ep1_1HB_1","phenotype"] <- paste(merged_phenos["81_101ep1_1HB_1","phenotype"] , "b?")
merged_phenos["6_110ep1_1RA_1","phenotype"] <- paste(merged_phenos["6_110ep1_1RA_1","phenotype"] , "b?")
merged_phenos["83_102ep1_1RB_1","phenotype"] <- paste(merged_phenos["83_102ep1_1RB_1","phenotype"] , "b?")
merged_phenos["74_67ep1_1HB_1","phenotype"] <- paste(merged_phenos["74_67ep1_1HB_1","phenotype"] , "b?")
merged_phenos["86_81ep1_1HN_1","phenotype"] <- paste(merged_phenos["86_81ep1_1HN_1","phenotype"] , "b?")
merged_phenos["43_58ep1_1HB_1","phenotype"] <- paste(merged_phenos["43_58ep1_1HB_1","phenotype"] , "b?")
merged_phenos["4_109ep1_1RA_1","phenotype"] <- paste(merged_phenos["4_109ep1_1RA_1","phenotype"] , "b?")
merged_phenos[oth2,"phenotype"] <- unlist(lapply(oth2 , function(l) paste(iris_data1[l,c("combined.Virology","Bacteria","diag.group.for.paper")],collapse = "_")))
samples_complex <-  c(samples_complex , oth2,
                      "4_109ep1_1RA_1","43_58ep1_1HB_1","86_81ep1_1HN_1",
                      "74_67ep1_1HB_1","83_102ep1_1RB_1","6_110ep1_1RA_1","81_101ep1_1HB_1",
                      "117_91ep1_1P_1","120_116ep1_1M_1","164_163ep1_1M_1",
                      "173_176ep1_1A_1","180_290ep1_1F_1","228_302ep1_1U_1")

num <- "GSE30119"
samples <- rownames(pData(datasets[[num]]))
pData(datasets[[num]])[1,]
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.4"])) * 12
merged_phenos[samples,"gender"] <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.6"])
merged_phenos[samples,"ethnicity"] <-  gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.5"])
phenos <-gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.3"])
phenos[grep("Health",phenos)] <- "HC"
phenos[grep("aureus",phenos)] <- "saureus"
merged_phenos[samples,"phenotype"] <-   phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) +1
merged_phenos[samples,"experiment"] <- num


num <-  "GSE39941"
samples <- rownames(pData(datasets[[num]]))
alt_name <- unlist(lapply(pData(datasets[[num]])[,"title"],function(l) strsplit(as.character(l) , split = "_")[[1]][4]))
pData(datasets[[num]])[1,]

anderson_pheno1 <- read.table("dataset/mlw.csv",sep="\t",stringsAsFactors = F,header= T , fill = NA,quote = "\"")
anderson_pheno2 <- read.table("dataset/rxh.csv",sep="\t",stringsAsFactors = F,header= T , fill = NA,quote= "\"")
anderson_pheno3 <- read.table("dataset/KDH.csv",sep="\t",stringsAsFactors = F,header= T , fill = NA,quote= "\"")

rownames(anderson_pheno1) <- unlist(lapply(anderson_pheno1[,1],function(l) strsplit(l,"-")[[1]][2]))
rownames(anderson_pheno2) <- unlist(lapply(anderson_pheno2[,1],function(l) strsplit(l,"-")[[1]][2]))
rownames(anderson_pheno3) <- anderson_pheno3[,"tbid"]
anderson_pheno <- rbind(anderson_pheno1 , anderson_pheno2)
kenyan.alt <-  alt_name[!alt_name%in% rownames(anderson_pheno)] 
kenyan.sam <- samples[!alt_name%in%rownames(anderson_pheno)]


merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"alt_id"] <- alt_name
merged_phenos[samples,"age"] <- anderson_pheno[alt_name,"ageint1"] # as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.4"])) * 12
merged_phenos[kenyan.sam,"age"] <- anderson_pheno3[kenyan.alt,"agem"]


merged_phenos[samples,"gender"] <- anderson_pheno[alt_name,"gend"] # gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.6"])
merged_phenos[kenyan.sam,"gender"] <- anderson_pheno3[kenyan.alt,"male"]
merged_phenos[kenyan.sam,"gender"][merged_phenos[kenyan.sam,"gender"] ==1] <- "m"
merged_phenos[kenyan.sam,"gender"][merged_phenos[kenyan.sam,"gender"] ==0] <- "f"
merged_phenos[kenyan.sam,"gender"][is.na(merged_phenos[kenyan.sam,"gender"])] <- "?"

merged_phenos[samples,"ethnicity"] <-  gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1"])
phenos[phenos == "other disease"] <- "od"
phenos[phenos ==    "latent TB infection" ] <- "HC (LTBI)"
phenos[phenos ==   "active tuberculosis (culture confirmed)"] <- "TBpos"
phenos[phenos ==   "active tuberculosis (culture negative)"] <- "TBneg"
phenos[phenos ==  "active tuberculosis"  ] <- "TB"


merged_phenos[samples,"phenotype"] <-  phenos
od_ids <- merged_phenos[samples,][merged_phenos[samples,"phenotype"] == "od",1]
od_alt_ids <- merged_phenos[od_ids , "alt_id"]
tmp <- od_alt_ids[od_alt_ids%in%kenyan.alt]
tmp2 <- od_ids[od_alt_ids%in%kenyan.alt]
merged_phenos[tmp2,"phenotype"][anderson_pheno3[tmp,"malaria"] == 1] <- "malaria"

merged_phenos[samples[alt_name %in% rownames(anderson_pheno1)],"batch"] <- max(merged_phenos[,"batch"],na.rm =T ) +1
merged_phenos[samples[alt_name %in% rownames(anderson_pheno2)],"batch"] <- max(merged_phenos[,"batch"],na.rm =T ) +1
merged_phenos[samples[alt_name %in% rownames(anderson_pheno3)],"batch"] <- max(merged_phenos[,"batch"],na.rm =T ) +1

merged_phenos[samples,"experiment"] <- num
samples_to_exclude <- c(samples_to_exclude ,samples[pData(datasets[[num]])[,"characteristics_ch1.1"]=="hiv status: HIV positive"])

num <- "GSE42026"
samples <- rownames(pData(datasets[[num]]))

h1n1_pheno <- read.table("dataset/patients_in_the_h1n1_paper.csv",sep="\t",stringsAsFactors = F,header= T , fill = NA,quote= "\"")
h1n1_pheno <- h1n1_pheno[!is.na(h1n1_pheno[,"microarray.sample.number"]),]
tmp <- unlist(lapply(as.character(pData(datasets[[num]])[,1]) , function(l) strsplit(l , "m")[[1]][2]))


merged_phenos[samples,"samples"] <- samples
for(j  in 1:length(samples)){
  dat <- h1n1_pheno[h1n1_pheno[,"microarray.sample.number"] == tmp[j],]
  sam <- samples[j]
  merged_phenos[sam,"age"] <- (12 * dat[,"Age"]) 
  merged_phenos[sam , "gender"] <- substring(dat[,"Sex"] , 1 , 1)
  merged_phenos[sam , "ethnicity"] <- dat[,"Ethnicity"] 
  
  if(dat[,"diag.group.for.paper"]=="R"){
    merged_phenos[sam , "phenotype"] <- "RSV" 
  }else if(dat[,"diag.group.for.paper"]=="H"){
    merged_phenos[sam , "phenotype"] <- "flu"
  }else if(dat[,"diag.group.for.paper"]=="C"){
    merged_phenos[sam , "phenotype"] <- "HC"
  }else if(dat[,"diag.group.for.paper"]=="bact"){
    bact <- dat[,"Bacteria"]
    if(str_detect(bact,"neumo")){
      merged_phenos[sam , "phenotype"] <- "pneumo"  
    }else if(bact %in% c("GrpA Strep Blood" ,
                         "Gp A Strep blood cx",
                         "GAS wound" )){
      merged_phenos[sam , "phenotype"] <- "gas"  
    }else if(bact %in% c("PVL staph in pleural fluid, pseudomonas on NPA & throat swab; BC neg" ,
                         "Gp A Strep, C. difficile")){
      merged_phenos[sam , "phenotype"] <- bact  
      merged_phenos[sam , "group"] <- "complex/uncertain"  
    }else if( bact == "Staph aureus blood") {
      merged_phenos[sam , "phenotype"] <- "staph"  
    }else{
      print(bact)
    }
  }
}

merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm =T ) +1
merged_phenos[samples,"experiment"] <- num


num <- "GSE64456"
path_table <- read.table("dataset/GSE64456-species_list_mod.txt",header = T )
samples <- rownames(pData(datasets[[num]]))
pData(datasets[[num]])[1,]
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1"])) / 30 # age in days -> months
merged_phenos[samples,"gender"] <-gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.1"])
merged_phenos[samples,"ethnicity"] <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.2"])
phenos <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.5"])

phenos[phenos=="SBI-(ENTERO)"] <- "enterovirus"
phenos[phenos=="SBI-(ENTERO) / RHINO"] <- "enterovirus + rhino"
phenos[phenos=="SBI- ENTEROVIRUS_NSBI"] <- "enterovirus"
phenos[phenos=="SBI- ENTEROVIRUS_NSBI / RHINO"] <- "enterovirus + rhino"
phenos[phenos%in%c("SBI- flu+  H1N1-_NSBI","SBI- flu+  H1N1+_NSBI")] <- "flu"
phenos[phenos=="SBI-(No viral test)"] <- "U"
phenos[phenos=="SBI- VIRUS+(PARAIN)_NSBI"] <- "parainfluenza"
phenos[phenos=="CTRL"] <- "HC"
merged_phenos[samples,"phenotype"] <-  phenos

phenos2 <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.4"])  
tmp <- pData(datasets[[num]])

bactage <- gsub(".*: ","",tmp[phenos2=="SBI","characteristics_ch1"])
bactgend <- gsub(".*: ","",tmp[phenos2=="SBI","characteristics_ch1.1"])
bactsam <- as.character(tmp[phenos2=="SBI",2])
bactsite <- gsub(".*: ","",tmp[phenos2=="SBI","characteristics_ch1.5"])
bactsite[bactsite=="Bacteremia" ] <- "Blood?"
bactsite[bactsite=="UTI" ] <- "Urine"
bactsite[bactsite=="SBI+urine_SBI" ] <- "Urine"
bactsite[bactsite=="SBI+blood/CSF_SBI" ] <- "Blood/CSF"
bactsite[bactsite=="SBI+urine/blood/CSF_SBI" ] <- "Blood/Urine/CSF"
bactsite[bactsite=="SBI+blood_SBI" ] <- "Blood"
bactsite[bactsite=="SBI+urine/blood_SBI" ] <- "Blood/Urine"

path <- "--"
path_table_assembled <- data.frame(cbind(bactsam,bactage,bactgend,path),stringsAsFactors = F)
rownames(path_table_assembled) <- path_table_assembled[,1]
#  path_table_assembled[1,"path"] <- 
for(j in 1:length(bactsam)){
  print(c(bactsam[j],bactsite[j]))
  potential.t <- path_table[path_table[,"Age.days."] == bactage[j] & path_table[,"Gender"] == bactgend[j] ,]
  potential <- as.character(unique(potential.t[,"Organism"]))
  if(path_table_assembled[j,4]=="--"){
    if(length(potential)==1){
      path_table_assembled[j,"path"] <- potential
      
    }else{
      print(potential.t)
    }
  }
}

path_table_assembled["GSM1571387","path"] <-  "Klebsiella-pneumoniae"
path_table_assembled["GSM1571390","path"] <-  "E.coli"
path_table_assembled["GSM1571563","path"] <-  "S.aureus"
path_table_assembled["GSM1571567","path"] <-  "GBS"
path_table_assembled["GSM1571586","path"] <-  "E.coli"
path_table_assembled["GSM1571556","path"] <-  "GBS"
path_table_assembled["GSM1571395","path"] <-  "E.coli"
path_table_assembled["GSM1571394","path"] <-  "E.coli"
path_table_assembled["GSM1571577","path"] <-  "E.coli"
path_table_assembled["GSM1571409","path"] <-  "E.coli" # 53 m 
path_table_assembled["GSM1571562","path"] <-  "GBS"
path_table_assembled["GSM1571583","path"] <-  "E.coli"
path_table_assembled["GSM1571396","path"] <-  "E.coli"
path_table_assembled["GSM1571573","path"] <-  "Viridans-streptococcus"
path_table_assembled["GSM1571568","path"] <-  "S.aureus"
path_table_assembled["GSM1571582","path"] <-  "gneg"
path_table_assembled["GSM1571578","path"] <-  "gneg"
path_table_assembled["GSM1571388","path"] <-  "gneg"
path_table_assembled["GSM1571557","path"] <-  "gneg"
path_table_assembled["GSM1571412","path"] <-  "gneg"
path_table_assembled["GSM1571397","path"] <-  "gneg"


unassigned <- rownames(path_table_assembled)[path_table_assembled[,4]=="--"]

unassigned <- c(unassigned,"GSM1571397","GSM1571412","GSM1571388","GSM1571557","GSM1571578","GSM1571582")
path_table_assembled[path_table_assembled[,4]=="--","path"] <- "bact"
merged_phenos[rownames(path_table_assembled),"phenotype"] <- path_table_assembled[,"path"]

batchno <-  max(merged_phenos[,"batch"],na.rm =T ) + 1
tmp <- pData(datasets[[num]])
merged_phenos[as.character(tmp[tmp[,27]=="Illumina HT12 V4 R1",2]),"batch"] <- batchno
merged_phenos[as.character(tmp[tmp[,27]=="Illumina HT12 V4 R2",2]),"batch"] <-  batchno + 1

merged_phenos[samples,"experiment"] <- num

num <- "GSE34404" # controls are healthy -- siblings of sickle cell cases
samples <- rownames(pData(datasets[[num]]))
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  as.numeric(gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.6"])) * 12
merged_phenos[samples,"gender"] <- gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.5"])
merged_phenos[samples,"ethnicity"] <- "?"
phenos <-  gsub(".*: ","",pData(datasets[[num]])[samples,"characteristics_ch1.7"])
phenos[grep("POSITIF",phenos)] <- "malaria"
phenos[grep("NEGATIF",phenos)] <- "HC"
merged_phenos[samples,"phenotype"] <- phenos
merged_phenos[samples,"batch"] <- max(merged_phenos[,"batch"],na.rm = T) + 1
merged_phenos[samples,"experiment"] <- num

samples <- colnames(mega_expr)
merged_phenos[samples,"samples"] <- samples
merged_phenos[samples,"age"] <-  status[samples,"Age..months."]
merged_phenos[samples,"gender"] <- as.character(status[samples,"Sex"])
merged_phenos[samples,"ethnicity"] <-  as.character(status[samples,"ethnicity"])
# get more detailed phenos 
conv_mega_ids <- as.character(mega_linkertab[samples,2])
tmp_tab <- mega_phenos2[mega_phenos2[,"site.number.episode"]%in%conv_mega_ids,c("site.number.episode","X","BorV.coding")]
duplicates <-  mega_phenos2[mega_phenos2[,"site.number.episode"]%in%names(table(tmp_tab[,1])[table(tmp_tab[,1])>1]),1]
samples_to_exclude <- c(samples_to_exclude,"greyu_14_SMH","greyu_15_SMH",rownames(status)[ "exclude" == status[,"general"]],
                        "greyu_16_SMH","JIAact_3_AMC","JIAexa_16_AMC","JIAexa_18_AMC" )# for duplicate jia take younger
tmp_tab <- unique(tmp_tab)
tmp_tab <- tmp_tab[tmp_tab[,1]!= "SMH-389-1", ]
rownames(tmp_tab) <- tmp_tab[,1]
tmp_tab <- cbind(labels[samples,],tmp_tab[conv_mega_ids,],samples)
rownames(tmp_tab) <- samples
tmp_tab <- tmp_tab[!is.na(tmp_tab[,4]),]
tmp_tab[tmp_tab[,4]=="",4] <- as.character(tmp_tab[tmp_tab[,4]=="",1])
modsamples <- samples[samples%in%rownames(tmp_tab)]
merged_phenos[samples,"phenotype"] <- as.character(status[samples,"general"])
tmp_tab2 <- cbind(as.character(status[samples[samples%in%rownames(tmp_tab)],"general"]),unlist(lapply(modsamples, function(l) if(tmp_tab[l,"X"]!=""){
  return(paste(tmp_tab[l,c("BorV.coding","X")],collapse = ","))
}else{return(tmp_tab[l,"BorV.coding"])})))
rownames(tmp_tab2) <- samples[samples%in%rownames(tmp_tab)]
tmp_tab2[tmp_tab2[,2] == "t",2] <- "TBneg"
tmp_tab2[tmp_tab2[,2] == "T",2] <- "TBpos"
tmp_tab2[tmp_tab2[,1] == "KDconv" & tmp_tab2[,2] == "C",2] <- "KDconv"
tmp_tab2[tmp_tab2[,1] == "adeno" & tmp_tab2[,2] == "V",2] <- "V,adeno"
tmp_tab2[tmp_tab2[,1] == "RSV" & tmp_tab2[,2] == "V",2] <- "V,RSV"
tmp_tab2[tmp_tab2[,1] == "HSP" & tmp_tab2[,2] == "J",2] <- "HSP"

merged_phenos[modsamples,"phenotype"]  <- tmp_tab2[,2]
merged_phenos[modsamples,"phenotype"] [merged_phenos[modsamples,"phenotype"] =="C"] <- "HC"

complex <- modsamples[merged_phenos[modsamples,"phenotype"] %in% c("B","B,ecoli pseudomonas")]
complex <- c(complex ,modsamples[grep("grey",merged_phenos[modsamples,"phenotype"])])
complex <- c(complex ,modsamples[grep("prob",merged_phenos[modsamples,"phenotype"])])
complex <- c(complex ,modsamples[grep("VB",merged_phenos[modsamples,"phenotype"])])
complex <- c(complex ,modsamples[grep("BorV",merged_phenos[modsamples,"phenotype"])])
complex <- c(complex ,modsamples[grep("BT",merged_phenos[modsamples,"phenotype"])])
complex <- c(complex ,modsamples[grep("V,.* .*",merged_phenos[modsamples,"phenotype"])])
samples_complex <- c(samples_complex , complex)
table(merged_phenos[modsamples[!modsamples%in%complex],"phenotype"])

merged_phenos[modsamples[!modsamples%in%complex],"phenotype"] <- gsub(".,","",merged_phenos[modsamples[!modsamples%in%complex],"phenotype"])
kenyan_mega <- samples[merged_phenos[samples,"phenotype"]%in%c("TB","OD","HC (LTBI)")]
not_ken <- samples[!samples%in%kenyan_mega]

merged_phenos[kenyan_mega,"batch"] <- max(merged_phenos[,"batch"],na.rm=T) +1
merged_phenos[not_ken,"batch"] <- max(merged_phenos[,"batch"],na.rm=T) +1
merged_phenos[samples,"experiment"] <- "mega"


# filter
merged_phenos[,2] <- as.numeric(merged_phenos[,2])
merged_phenos[is.na(merged_phenos[,2]),2] <- -40
paediatric_sams <- rownames(merged_phenos)[as.numeric(merged_phenos[,2])<(12 * 18)] # AGE
length(paediatric_sams)
dim(merged_phenos)
merged_phenos[,3] <- tolower(substring(merged_phenos[,3],1,1))

merged_phenos <- merged_phenos[paediatric_sams,]
merged_phenos <- merged_phenos[!rownames(merged_phenos)%in%samples_to_exclude,]


simple_phenos <- paediatric_sams[!paediatric_sams%in%c(samples_complex,samples_to_exclude)]

table(merged_phenos[simple_phenos,5])

merged_phenos[merged_phenos[,5]=="E.coli",5] <- "ecoli"
merged_phenos[merged_phenos[,5]=="control",5] <- "HC"
merged_phenos[merged_phenos[,5]=="cKD",5] <- "KD" # complete KD
merged_phenos[merged_phenos[,5]=="inKD",5] <- "KD" # incomplete KD
merged_phenos[merged_phenos[,5]=="GAS",5] <- "gas"
merged_phenos[merged_phenos[,5]=="GAS/SF",5] <- "gas" # strep faecalis
merged_phenos[merged_phenos[,5]=="Staph",5] <- "staph"
merged_phenos[merged_phenos[,5]=="rsv",5] <- "RSV"
merged_phenos[merged_phenos[,5]=="GBS",5] <- "gbs"
merged_phenos[merged_phenos[,5]=="pSLE",5] <- "SLE"
merged_phenos[merged_phenos[,5]=="PSLE",5] <- "SLE"
merged_phenos[merged_phenos[,5]=="Adenovirus",5] <- "adeno"
merged_phenos[ merged_phenos[,5]%in%c("JIAact" , "JIAexa","illness: Still") ,5] <- "JIA"
merged_phenos[ merged_phenos[,5]%in%c(  "TBpos") ,5] <- "TB" # only culture positives
merged_phenos[ grep("HAdV",merged_phenos[,5]) ,5] <- "adeno"

samples_complex <- c(samples_complex ,
                     rownames(merged_phenos)[merged_phenos[,5]%in%c("enterococcus staph","flu RSV","RSV MPV Boca RV pneumo","V",
                                                                    "Bacteria","infected","hsv1","nonJIA","parecho","Salmonella",
                                                                    "kingella","posJIA","Gp A Strep, C. difficile","PVL staph in pleural fluid, pseudomonas on NPA & throat swab; BC neg",
                                                                    "H1N1 boca","Tbneg",  "H1N1 RV boca", "H1N1 RV PIV 4 pneumo","RSV RV"
                                                                    ,"FUO","greyu","RSV adeno boca","RSV boca","RSVFluA")])
samples_complex <- unique(samples_complex)

small_groups <- names(table(merged_phenos[simple_phenos,5]))[table(merged_phenos[simple_phenos,5])<10]
small_group_samples <-rownames(merged_phenos)[merged_phenos[,"phenotype"]%in%small_groups]

samples_to_exclude <- c(samples_to_exclude , "TBneg_36_SMH","TBneg_35_SMH","exclude_3_SMH",
                        rownames(merged_phenos)[merged_phenos[,5]%in%c("KDconv","Rexcl - technical","exclude","FPIESp","FPIESf")]) 

big_groups <- names(table(merged_phenos[simple_phenos,5]))[table(merged_phenos[simple_phenos,5])>=10]

tmp <- merged_phenos[rownames(merged_phenos)[merged_phenos[,"phenotype"]%in%big_groups],]
ordids <- rownames(tmp)[order(tmp[,"phenotype"])]

merged_phenos[ordids,"group"] <- merged_phenos[ordids,"phenotype"]
merged_phenos[small_group_samples,"group"] <- "smallgroup"
merged_phenos[samples_complex,"group"] <- "complex/uncertain"

sams_use <- rownames(merged_phenos)
b9sams <- rownames(merged_phenos)[merged_phenos[,"batch"]=="9"]
rownames(merged_phenos)[rownames(merged_phenos)%in%b9sams] <-  paste("X",rownames(merged_phenos)[rownames(merged_phenos)%in%b9sams],sep="")
colnames(merged_expr)[colnames(merged_expr)%in%b9sams] <- paste("X",colnames(merged_expr)[colnames(merged_expr)%in%b9sams],sep="")




merged_phenos <- merged_phenos[!is.na(merged_phenos[,"experiment"]),]


# premature baby exclusion of GSE25504[[4]]
num <-  "GSE25504[[4]]"
pData(datasets[[num]])[1,]
samples <- rownames(pData(datasets[[num]]))

age_tmp <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.6"])
age_tmp <- as.numeric(gsub("\\+.*","",age_tmp))

birthweight <- as.numeric(gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1.7"]))
phe <- gsub(".*: ","",pData(datasets[[num]])[,"characteristics_ch1"])


plot(age_tmp,birthweight, col = c("black","red")[as.numeric(phe=="infected")+1])
# the curve of age vs birthweight shows that these are mostly quite premature and small

# fix mega phenos
mega_sams <- rownames(merged_phenos)[merged_phenos[,"experiment"]=="mega"] 
tvec <- merged_phenos[mega_sams,"phenotype"] %in% c("B","K","T","V","t")
merged_phenos[mega_sams,"phenotype"][tvec]  <- labels[mega_sams[tvec],]

sams_use <- rownames(merged_phenos)[rownames(merged_phenos) %in% colnames(merged_expr)]

# link GSE42026 to the iris 
tmp_pheno <- pData(datasets$GSE42026)
tmpsams <- merged_phenos[merged_phenos[,"experiment"]=="GSE42026","samples"]
barcode <- as.character(tmp_pheno[tmpsams , "description"])

tmp_pheno2 <- mega_phenos2[mega_phenos2[,"barcode.if.run.on.IRIS1..v3."]%in%barcode,]
rownames(tmp_pheno2) <- tmp_pheno2[,"barcode.if.run.on.IRIS1..v3."]

merged_phenos[tmpsams,"alt_id"] <- tmp_pheno2[barcode,"array.name.as.used.in.IRIS1"]
# "5147466019_H" missing but control

tmpsams <- tmpsams[barcode%in%rownames(tmp_pheno2)]
barcode <- barcode[barcode%in%rownames(tmp_pheno2)]
names(tmpsams) <- barcode
rownames(tmp_pheno2) <- tmpsams[rownames(tmp_pheno2)] 

excl_sams <- tmpsams[tmp_pheno2[tmpsams,"BorV.coding"]%in%c("BorV","INForINF","VB","VprobB")]
keep_sams <- tmpsams[!tmp_pheno2[tmpsams,"BorV.coding"]%in%c("BorV","INForINF","VB","VprobB")]
merged_phenos[excl_sams,"group"] <- "complex/uncertain"
merged_phenos["GSM1030796","group"] <- "complex/uncertain" # osteopetrosis


# combat coconut #########

# san diego batch effect :
merged_phenos.filt <- merged_phenos
mega_sams <- merged_phenos.filt[merged_phenos.filt[,"experiment"]=="mega",1]
mega_sams <- mega_sams[merged_phenos.filt[mega_sams , "batch"] != 18 ] 
# batch 18 is kenyan tbs with no latent infections to use in coconut


ucsd <- mega_sams[grep("UCSD",mega_sams)]
batch <- rep(1,length(mega_sams))
batch[mega_sams%in% ucsd] <- 2

mod <- model.matrix( ~ as.factor(group), data = merged_phenos.filt[mega_sams,])
mod0 <- model.matrix( ~ 1 , data = merged_phenos.filt[mega_sams,])

mod_combat  <- model.matrix( ~ 1  + phenotype , data = merged_phenos[mega_sams,])
combat_edata <- ComBat(dat = merged_expr[,mega_sams] , batch = batch, mod = mod_combat ,   par.prior=TRUE, prior.plots=FALSE)
merged_expr[,mega_sams] <- combat_edata[,mega_sams]


# coconut ######################

merged_phenos[merged_phenos[,"group"]%in%c("HC","HC (LTBI)"),"batchcontrol"] <- 0
merged_phenos[!merged_phenos[,"group"]%in%c("HC","HC (LTBI)"),"batchcontrol"] <- 1


b9sams <- rownames(merged_phenos)[merged_phenos[,"batch"]=="9"]
b9sams <- c(b9sams,gsub("^X+","",rownames(merged_phenos)[rownames(merged_phenos)%in%b9sams]))
rownames(merged_phenos)[rownames(merged_phenos)%in%b9sams] <-  gsub("^X*","X",rownames(merged_phenos)[rownames(merged_phenos)%in%b9sams])
colnames(merged_expr)[colnames(merged_expr)%in%b9sams] <- gsub("^X*","X",colnames(merged_expr)[colnames(merged_expr)%in%b9sams])

totsams <- c()
coconut_input <- list()
for(batchno in unique(merged_phenos[,"batch"])){
  if(batchno!=18){ # leave out mega kenyan tbs
    print("------------")
    samples  <- rownames(merged_phenos)[merged_phenos[,"batch"]==batchno]
    samples <-samples[!is.na(samples)] 
    print(head(merged_phenos[samples,]))
    print(table(merged_phenos[samples,"batchcontrol"]))
    
    print(table(merged_phenos[samples,"batch"]))
    batchid <- paste(unique(merged_phenos[samples,"experiment"]),batchno,sep="_")
    if(!unique(merged_phenos[samples,"experiment"]) %in% c("GSE38900[[2]]","GSE25504[[4]]")){
      totsams <- c(totsams , samples)
      coconut_input[[batchid]] <-  list(genes = merged_expr[,samples],
                                        pheno = merged_phenos[samples,])
    }
  }
}


coconut_out <- COCONUT(GSEs = coconut_input ,
                       control.0.col = "batchcontrol",
                       byPlatform = FALSE)


## make gene matrices
COCONUTgenes <- Reduce(cbind, lapply(coconut_out$COCONUTList, function(x) x$genes))
rawgenes <- Reduce(cbind, lapply(coconut_out$rawDiseaseList, function(x) x$genes))

if (require(lumiHumanIDMapping)) lumiHumanIDMapping_nuID()#
mappingInfo <- nuID2RefSeqID(rownames(merged_expr), lib.mapping='lumiHumanIDMapping', returnAllInfo =TRUE)

gene_probes <- rownames(mappingInfo[mappingInfo[,"Symbol"] == "ACTB",])

batches <-merged_phenos[colnames(COCONUTgenes),"batch"]
colours <- rainbow(max(batches))

for(probe in gene_probes){
  plot(x=1:ncol(rawgenes), y=rawgenes[probe, ], pch=20, col="grey" , ylim = c(11,16))
  points(x=1:ncol(COCONUTgenes), y=COCONUTgenes[probe, ], pch=20, col = colours[batches])
  points(x=1:ncol(COCONUTgenes), y=COCONUTgenes[probe, ], pch=1, col = "black")
  runn <- 0
  runx <- 1
  for(i in 1:length(batches)){
    if(batches[i]!=runn){
      
      if(runx!=i){
        text(mean(c(runx , i)),min(COCONUTgenes[probe, ]) ,labels = batches[i -1 ],cex = 1)
        runx <- i
      }
      runn <- batches[i]
      abline(v = i - 0.5)
    }
  }
  abline(v = i + 0.5)
  text(mean(c(runx , i)),min(COCONUTgenes[probe, ]) ,labels = batches[i  ],cex = 1)
}









# checking for samples shared between datasets ####

sams.tmp <- c(merged_phenos[merged_phenos[,"experiment"]%in%c("GSE72810"),1])
sams.tmp <- gsub("^","X",sams.tmp)
sams.tmp <- c(sams.tmp,merged_phenos[merged_phenos[,"experiment"]%in%c("GSE42026"),1])
sams.tmp <- sams.tmp[sams.tmp%in% colnames(COCONUTgenes)]
pca <- prcomp(t(COCONUTgenes[,sams.tmp]))
plot(pca$x[,c(1,2)], col = c("red","green")[1+as.numeric(merged_phenos[rownames(pca$x),"experiment"]=="GSE42026")])


cor_mat <- cor(COCONUTgenes[,sams.tmp])

plot(COCONUTgenes[,"X51_50ep1_0H_1"] , COCONUTgenes[,"GSM1030845"])
plot(cor_mat["X51_50ep1_0H_1",], col = c("red","green")[1+as.numeric(colnames(cor_mat)%in%c("GSM1030845","X51_50ep1_0H_1"))], ylim = c(0.8,1))
plot(cor_mat["GSM1030845",], col = c("red","green")[1+as.numeric(colnames(cor_mat)%in%c("GSM1030845","X51_50ep1_0H_1"))], ylim = c(0.8,1))
# pca showing these two samples


cor_mat <- cor(COCONUTgenes)

filt_experiments <-unique(merged_phenos[,"experiment"])
filt_experiments <- filt_experiments[!filt_experiments%in%c(  "GSE38900[[2]]")]
experiment_pair <- combn(filt_experiments,2)


dup <- c()
for(i in 1:ncol(experiment_pair)){ # pairwise dataset correlation of samples
  sams_a <- rownames(merged_phenos)[merged_phenos[,"experiment"]== experiment_pair[1,i] & merged_phenos[,"group"]!="HC"]
  sams_b <- rownames(merged_phenos)[merged_phenos[,"experiment"]== experiment_pair[2,i] & merged_phenos[,"group"]!="HC"]
  
  sams.tmp <- c(sams_a,sams_b)
  superheat(cor_mat[sams.tmp,sams.tmp], heat.pal =  c("green","white","red","black"), heat.pal.values = c(0,0.9,0.95,1),title =paste(experiment_pair[,i], collapse = "-"))
  dup <- c(dup,readline(prompt="has duplicates? ")  )
}
inter_batch_duplicates <- experiment_pair[,dup==1]

dup <- c()
for(i in unique(phenotypes[,"experiment"])){ # within dataset duplicates
  sams.tmp <- rownames(phenotypes)[phenotypes[,"experiment"]==i ]
  
  superheat(cor_mat[sams.tmp,sams.tmp], heat.pal =  c("green","white","red","black"), 
            heat.pal.values = c(0,0.9,0.99,1),title =i)
  dup <- c(dup,readline(prompt="has duplicates? ")  )
  
}

# possibly internal
# "GSE22098" "mega" 

tmp <- cor_mat[sams.tmp,sams.tmp]
tmp[tmp == 1] <- 0
phenotypes[rownames(tmp[rowMaxs(tmp)==max(tmp),colMaxs(tmp)==max(tmp)]),]

# fixing inter batch duplicates ####

inter_batch_duplicates <- experiment_pair[,c(37,60)]
#  "GSE22098" "GSE72810"
#  "GSE30119" "GSE42026"


sams_a <- rownames(merged_phenos)[merged_phenos[,"experiment"]== inter_batch_duplicates[1,1]] 
sams_b <- rownames(merged_phenos)[merged_phenos[,"experiment"]== inter_batch_duplicates[2,1]]

sams_a <- sams_a[sams_a%in%rownames(cor_mat)]
sams_b <- sams_b[sams_b%in%rownames(cor_mat)]

# not any besides staph
sams.tmp <- c(sams_a,sams_b)
sams.tmp <- sams.tmp[merged_phenos[sams.tmp,"group"]!= "staph"]
superheat(cor_mat[sams.tmp,sams.tmp], heat.pal =  c("green","white","red","black"), heat.pal.values = c(0,0.9,0.95,1),
          title =paste(experiment_pair[,37], collapse = "-"))


# staph
sams.tmp <- c(sams_a,sams_b)
sams.tmp <- sams.tmp[merged_phenos[sams.tmp,"group"]== "staph"]
sams_a.loc <- sams_a[merged_phenos[sams_a,"group"]== "staph"]
sams_b.loc <- sams_b[merged_phenos[sams_b,"group"]== "staph"]

tmp <- cor_mat[sams.tmp,sams.tmp]
tmp <- tmp[sams_a.loc,]
for(i in 1:length(sams_a.loc)){
  tmp[i,i] <- 0.6
}

matchup <- c()
for(i in 1:length(sams_a.loc)){
  tmp[i,tmp[i,]==max(tmp[i,])] <- 1
  matchup <- c(matchup, colnames(tmp)[tmp[i,]==max(tmp[i,])])
}
names(matchup) <- sams_a.loc
for( i in sams_a.loc) tmp[i,matchup[i]] <-  2
superheat(tmp, heat.pal =  c("green","white","red","black","blue"), heat.pal.values = c(0,0.9,0.95,1,2))




for( i in sams_a.loc){ # perfect match up  for ethnicity and gender
  if((merged_phenos[i,"age"]!= merged_phenos[matchup[i],"age"])) print(c(merged_phenos[i,"age"], 
                                                                         merged_phenos[matchup[i],"age"]))
}
# difference is rounding to years in sams_a dataset so want to keep the  sams_b

merged_phenos[names(matchup),"note"] <- "exclude" # have removed the s. aureus cases from 22098,

# iris overlap
sams_a <- rownames(merged_phenos)[merged_phenos[,"experiment"]== inter_batch_duplicates[1,2]]
sams_b <- rownames(merged_phenos)[merged_phenos[,"experiment"]== inter_batch_duplicates[2,2]]
sams_a <- sams_a[sams_a%in%rownames(cor_mat)]
sams_b <- sams_b[sams_b%in%rownames(cor_mat)]


sams.tmp <- c(sams_a,sams_b)
tmp <- cor_mat[sams.tmp,sams.tmp]
tmp <- tmp[sams_a,] 
matchup  <- c()
for(i in sams_a)tmp[i,i] <- 0.8
for(i in sams_a){
  tosam <- colnames(tmp)[tmp[i,]==max(tmp[i,])]
  matchup  <- c(matchup ,tosam)
  tmp[i,tosam] <- 2
}
names(matchup) <- sams_a

superheat(tmp, heat.pal =  c("green","white","red","black","blue"), heat.pal.values = c(0,0.9,0.95,1,2))


all(merged_phenos[sams_a,"gender"] == merged_phenos[matchup,"gender"])
cbind(merged_phenos[sams_a,"group"], merged_phenos[matchup,"group"])  


all(merged_phenos[sams_a,"age"]== merged_phenos[matchup,"age"])  

tmp <- cor_mat[sams_b[sams_b %ni% matchup],]
for(i in rownames(tmp)){
  tmp[i,i] <- 0
  tosam <- colnames(tmp)[tmp[i,]==max(tmp[i,])]
  print(c(i, tosam))
} 


# there are some samples in b which arent in a
sams_a <- rownames(merged_phenos)[merged_phenos[,"experiment"]== inter_batch_duplicates[1,2]]
merged_phenos[sams_a,"note"] <- "exclude"

merged_phenos[is.na(merged_phenos[,"note"]),"note"] <- "-"

tmp <- cor_mat[merged_phenos[rownames(cor_mat),"note"]!= "exclude",
               merged_phenos[rownames(cor_mat),"note"]!= "exclude"]
for(i in 1:ncol(tmp)){
  tmp[i,i] <- 0
}

for(i in as.vector(tmp)[order(as.vector(tmp), decreasing = T)][1:20]){
  print(i)
  print(merged_phenos[rownames(tmp[rowMaxs(tmp)==i,colMaxs(tmp)==i]),])  
}




# PCA ####

# add back 38900[[2]] to show on pca- ignore in main pipeline

array <- datasets[["GSE38900[[2]]"]]
print(c(dim(exprs(array)),dim(pData(array))))
probes <- rownames(exprs(array))
nuids <- IlluminaID2nuID(probes , species = "Human")[,"nuID"]
present_ids[[i+1]] <- nuids
exprs_tmp <- exprs(array)
rownames(exprs_tmp) <- nuids
shared <- intersect(rownames(merged_expr) , rownames(exprs_tmp))
merged_expr <- cbind(merged_expr[shared , ],exprs_tmp[shared , ])

sams_use <- rownames(merged_phenos)[rownames(merged_phenos) %in% colnames(merged_expr)]

# before coconut
pca.uncorrected <- prcomp(t(merged_expr[,sams_use]),center=TRUE,scale.=TRUE)     

# filtering out the odd ramilo batch 
colnames(merged_expr_coconut) <- gsub("^X+","X",colnames(merged_expr_coconut))
sams_use <- sams_use[sams_use %in% colnames(merged_expr_coconut)]
pca.uncorrected_filtered <- prcomp(t(merged_expr[,sams_use]),center=TRUE,scale.=TRUE)     

# after coconut
sams_use <- sams_use[sams_use %in% colnames(merged_expr_coconut)]
pca.coconut <- prcomp(t(merged_expr_coconut[,sams_use]),center=TRUE,scale.=TRUE)     


## tidy phenotypes ########

phenotypes <- merged_phenos
expression_matrix <- merged_expr_coconut
colnames(expression_matrix) <- gsub("^X+","X",colnames(expression_matrix))
samples  <- colnames(expression_matrix)
phenotypes <- phenotypes[rownames(phenotypes) %in% samples,]


sams_noncomplex <- samples[phenotypes[samples,"group"]!= "complex/uncertain"] 
# complex/uncertain is samples with multiple detected pathogens or a diagnosis which is uncertain
unique(phenotypes[sams_noncomplex,c("phenotype","group")])

#two malaria and two tb ages set to mean for experiment, ages are just used for DE prefilter
phenotypes[phenotypes[,"age"]==-40 & phenotypes[,"experiment"]=="GSE34404" ,"age"] <- mean(phenotypes[phenotypes[,"age"]!=-40 & phenotypes[,"experiment"]=="GSE34404" ,"age"])
phenotypes[phenotypes[,"age"]==-40 & phenotypes[,"experiment"]=="GSE39941" ,"age"] <- mean(phenotypes[phenotypes[,"age"]!=-40 & phenotypes[,"experiment"]=="GSE39941" ,"age"])



png(paste(plotdir,"beeswarm_age-bypheno.png",sep=""),height = 1000,width = 2000)
pltdat <- phenotypes[,c("age","group") ]
phenos <- unique(phenotypes[,"group"])
pltdat[,2] <- match(pltdat[,2],phenos)
phecols <- rainbow(length(phenos ))
beeswarm(age ~ group ,data = pltdat,pch=16,labels = phenos ,corral ="wrap"
         , pwcol = phecols[match(phenotypes[,"group"],phenos)],cex = 1,las =2)
dev.off()

tmp <- phenotypes[phenotypes[,"experiment"]=="GSE22098" & phenotypes[,"phenotype"] == "Strep" & phenotypes[,"group"]!="complex/uncertain","samples"] 
phenotypes[tmp,c("group","phenotype")] <- "gas" # in the paper "Strep"s are all Group A

tmp <- phenotypes[phenotypes[,"phenotype"] == "V,rhino Asthma" ,"samples"] # are down as complex but need to include
phenotypes[tmp,c("group","phenotype")] <- "rhino"

table(phenotypes[phenotypes[,"group"]=="complex/uncertain","phenotype"])

phenotypes[phenotypes[,"group"]=="complex/uncertain" & phenotypes[,"phenotype"]%in%c("viral","staph","RSV","pneumo"),]
# iris RSVs contain some with potential bacterial infection and high crp

# Check 
phenotypes[phenotypes[,"group"]=="complex/uncertain" & grepl("flu",ignore.case = T,phenotypes[,"phenotype"]),]
# iris GSE42026 --BorVs and BVs see tmp_pheno2 in run_inc...
phenotypes[phenotypes[,"group"]=="complex/uncertain" & grepl("staph",ignore.case = T,phenotypes[,"phenotype"]),]
# GSM1030796 => immunodeficiency
phenotypes[phenotypes[,"group"]=="complex/uncertain" & grepl("rsv",ignore.case = T,phenotypes[,"phenotype"]),]
# IRIS BorVs
phenotypes[phenotypes[,"group"]=="complex/uncertain" & grepl("pneumo",ignore.case = T,phenotypes[,"phenotype"]),]
# IRIS VBs


# groups with too few samples
phenotypes[phenotypes[,"phenotype"]%in%c("Salmonella","kingella","klebsiella","Klebsiella-pneumoniae"
                                         ,"enterobacter cloacae","gneg","Enterobacter-cloacae","hinfluenzae",
                                         "gram positive bacterial infection","E faecalis",
                                         "E faecium","Enterococcus-faecalis","bacterialgpos")
                                                      ,"group"] <- "smallgroup"


phenotypes[phenotypes[,"phenotype"]%in%c("enterovirus + rhino")
           ,"group"] <- "complex/uncertain"
phenotypes[phenotypes[,"phenotype"]%in%c("SBI-(Viral -)","SBI- VIRUS-_NSBI","SBI- VTN_NSBI"),c("phenotype","group")] <- c("U","U")
phenotypes[phenotypes[,"phenotype"]%in%c("Enterovirus") & phenotypes[,"group"]=="smallgroup",c("phenotype","group")] <- c("enterovirus","enterovirus")



# collapse a couple of the groups:
mv_to_grp <- function(pheno,newgroup){
  phenotypes[phenotypes[,"phenotype"] %in% pheno & phenotypes[,"group"]!="complex/uncertain","group"] <- newgroup
  return(phenotypes)
}

phenotypes <- mv_to_grp(c("saureus","S.aureus","MSSA","MRSA","staph_aureus"),"staph")
phenotypes <- mv_to_grp(c("od","other disease (IGRA +)"),"OD")

phenotypes <- mv_to_grp(c("Rexcl - technical", "SBI- VIRUS+_NSBI",
                          "SBI- VIRUS+(RSV)_NSBI" ,"J","one"
                          ,"KDconv","exclude","FPIESf","FPIESp"),"exclude")
phenotypes <- mv_to_grp(c("hsv1"),"smallgroup")
phenotypes <- mv_to_grp(c("mening"),"meningococcal")
phenotypes <- mv_to_grp(c("gneg","gram positive bacterial infection","V","bacterialgpos","bact"),"uncertain")
phenotypes[phenotypes[,"group"] %in% c("U") , "group"] <- "uncertain"
phenotypes[grep("IVIG",rownames(phenotypes)),"group"] <- "exclude"


## which groups used in training ##########################################################

leaveout <- rownames(phenotypes[phenotypes[,"group"]%in%c("smallgroup","exclude","uncertain","complex/uncertain","TBneg"
                                                          ,"OD"),] )
leaveout <- c(leaveout , rownames(phenotypes)[ grep("exclude",rownames(phenotypes)) ])
leaveout <- c(leaveout , rownames(phenotypes)[ grep("TBneg",rownames(phenotypes)) ])
leaveout <- c(leaveout , rownames(phenotypes)[ phenotypes[,"note"]=="exclude"])
leavein <- rownames(phenotypes)[!rownames(phenotypes) %in% leaveout]

# training test split ##############

training_groups <- unique(phenotypes[leavein,"group"])

for(grp in training_groups){
  set.seed(112358)
  grp_sams <- rownames(phenotypes)[phenotypes["group"]==grp]
  grp_sams <- grp_sams[grp_sams%in%leavein]
  grp_sams <- sample(grp_sams,size = length(grp_sams))
  test_sams <- grp_sams[1:ceiling(length(grp_sams)*0.25)]
  training_sams <- grp_sams[!grp_sams %in% test_sams]
  phenotypes[test_sams,"trte"] <- "test"
  phenotypes[training_sams,"trte"] <- "train"#
}

phenotypes[is.na(phenotypes[,"trte"]),"trte"] <- "leftout"
save(phenotypes,file = "data_objects/phenotypes.RData")

expression_filt <- merged_expr_coconut[rowSums(merged_expr_coconut > 6) > 50,]

## set up nested cv structure  within 1st level training set #####

training_phenotypes <- phenotypes[phenotypes[,"trte"]=="train",] # just level 2 samples
training_phenotypes <- training_phenotypes[training_phenotypes[,"group"]%in%training_groups,]

table(phenotypes[phenotypes[,"trte" ]=="train","group"])
table(phenotypes[phenotypes[,"trte" ]=="test","group"])


# set up folds 
nfolds <- 10
for(grp in training_groups){
  set.seed(12345)
  grp_sams <- rownames(training_phenotypes)[training_phenotypes["group"]==grp]
  grp_sams <- sample(grp_sams,length(grp_sams))
  training_phenotypes[grp_sams,"fold"]   <- rep(1:nfolds,100)[1:length(grp_sams)]
  
}


# make indicator matrix
training_phenotypes_indicator <- class.ind(factor(training_phenotypes[,"group"])) 
rownames(training_phenotypes_indicator) <- rownames(training_phenotypes)

bacteria <- c("ecoli","gas","gbs","meningococcal","pneumo","staph")
inflammatory <- c("HSP","JIA","SLE")
viral <- c("adeno","enterovirus","flu","HHV6","rhino","RSV")


specific <- colnames(training_phenotypes_indicator)
training_phenotypes_indicator <- cbind(training_phenotypes_indicator,rep(0,nrow(training_phenotypes_indicator)),
                                       rep(0,nrow(training_phenotypes_indicator)),
                                       rep(0,nrow(training_phenotypes_indicator)))
colnames(training_phenotypes_indicator) <- c(specific , "bacterial","viral","inflammatory")


training_phenotypes_indicator[rowSums(training_phenotypes_indicator[,bacteria]) > 0,"bacterial"] <- 1
training_phenotypes_indicator[rowSums(training_phenotypes_indicator[,inflammatory]) > 0,"inflammatory"] <- 1
training_phenotypes_indicator[rowSums(training_phenotypes_indicator[,viral]) > 0,"viral"] <- 1



comparisons <- combn(training_groups,2)

## class weights ####

class_weights <- cbind(table(phenotypes[leavein,"group"]),NA)
class_weights[,2] <- c(2,10,2,3,10,10,2,1,5,9,10,10,10,2,2,8,10,8)
colnames(class_weights) <- c("freq" , "cost")
class_weights <- cbind(class_weights , class_weights[,2] / class_weights[,1])
colnames(class_weights)[3] <- "cost+imbalance"
imbalance <- 1/ class_weights[,1]
class_weights <- cbind(class_weights , imbalance)

class_weights_toplevel <- data.frame()
class_weights_toplevel[ "bacterial" ,"freq"] <- sum(class_weights[bacteria , "freq"])
class_weights_toplevel[ "viral" ,"freq"] <- sum(class_weights[viral , "freq"])
class_weights_toplevel[ "inflammatory" ,"freq"] <- sum(class_weights[inflammatory , "freq"])
class_weights_toplevel[ "malaria" ,"freq"] <- class_weights["malaria" , "freq"]
class_weights_toplevel[ "TB" ,"freq"] <- class_weights["TB" , "freq"]
class_weights_toplevel[ "KD" ,"freq"] <- class_weights["KD" , "freq"]

class_weights_toplevel[  ,"cost"] <- c(10,2,5,10,8,9)
class_weights_toplevel <- cbind(class_weights_toplevel , class_weights_toplevel[,2] / class_weights_toplevel[,1])
imbalance <- 1/ class_weights_toplevel[,1]
colnames(class_weights_toplevel)[3] <- "cost+imbalance"
class_weights_toplevel <- cbind(class_weights_toplevel , imbalance)

phenotypes <- phenotypes[!is.na(phenotypes[,1]),]
training_phenotypes <- training_phenotypes[training_phenotypes[,"experiment"] != "GSE25504[[4]]",] #  remove premature babies
phenotypes <- phenotypes[phenotypes[,"experiment"] != "GSE25504[[4]]",]
training_phenotypes_indicator <- training_phenotypes_indicator[rownames(training_phenotypes),]


# saving ####


save(expression_filt ,merged_phenos,
     phenotypes,training_phenotypes,training_phenotypes_indicator,
     test_sams,training_sams,mappingInfo,
     class_weights_toplevel,class_weights,training_groups,bacteria,viral,inflammatory,
     training_groups,file = "data_objects/pre-processed-fold-assigned.RData")

