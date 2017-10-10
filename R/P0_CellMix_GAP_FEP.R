
#CellMix approach:



# Install and load Libraries: ---------------------------


library(CellMix)
library(lumi)
library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)






# Set directories: ---------------------------

setwd("/home/daniel/Documents/A_Year_4_PHD/Paper_1_Diff_ex")
getwd()

data_dir <-"./data/"
P0_output_dir <-"./P0_Characterise/output/"
P0_figs_dir <-"./P0_Characterise/figs/"




# Load gene expression data (Should be LumiBatch object): ---------------------------

lumidata<-"GAP_FEP_Full_Gene_Expression_Data_Linear.RData"
load(file=paste(data_dir,lumidata,sep=""))



# Subset to good probes: ---------------------------

dim(eset_bg_log2_rsn_SVA)
eset_bg_log2_rsn_SVA_Good<-eset_bg_log2_rsn_SVA

#make exprs data with good probes
exprs_data<-exprs(eset_bg_log2_rsn_SVA)
feature_data<-fData(eset_bg_log2_rsn_SVA)

#get good probe
good_fdata<-filter(feature_data, good_probe =="TRUE")

exprs_into_lumibatch<-exprs_data[rownames(exprs_data)%in%good_fdata$nuID,]

################ REMOVE DUPLICATES BY SELCTING HIGHEST AVERAGE EXPRESSED PROBE ###############################
#rowmeans
exprsmean<-rowMeans(exprs_into_lumibatch)
#add id
rownames(good_fdata)<-good_fdata$nuID
exprsmean_ids<-merge(as.data.frame(exprsmean),good_fdata[,c("nuID","TargetID")],by="row.names")
# order data frame by truncated probe id and then expression level
exprsmean_ids<-exprsmean_ids[order(exprsmean_ids$TargetID, -exprsmean_ids$exprsmean), ]

# remove all duplicate probe id - keep one with highest mean expression
exprsmean_ids_unique<-exprsmean_ids[!duplicated(exprsmean_ids$TargetID),]

#reduce exprs set again
exprs_into_lumibatch2<-exprs_into_lumibatch[exprsmean_ids_unique$Row.names,]

good_fdata_small<-good_fdata[exprsmean_ids_unique$Row.names,]

#check they are same order
all.equal(good_fdata_small$nuID,rownames(exprs_into_lumibatch2))



#LUMIBATCH made
eset_bg_log2_rsn_SVA_Good<-eset_bg_log2_rsn_SVA[rownames(exprs_into_lumibatch2),]
fData(eset_bg_log2_rsn_SVA_Good)<-good_fdata_small

save(eset_bg_log2_rsn_SVA_Good,file=paste(P0_output_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))

