# Load Libaries and clean -------------------------------------------------

##cleanup
rm(list = ls())
dev.off()
##Libaries
library(reshape)
library(WGCNA)
library(flashClust)
library(data.table)
library(dplyr)
library(tidyr)


# Set directories ---------------------------------------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()
top_dir <-getwd()
data_dir <-"./data/Enrichment_libaries/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P2_output_dir <-"./P2_Enrichment/output/"



# Functions ---------------------------------------------------------------
#function to remove uneccesary lists, by overlap with background. 
smaller_enrichmentlists<-function(enrichmentlist,Background,overlap){
  print("number of categories at start")
  print(length(unique(enrichmentlist$Eclass)))
  #add column indicating presence in background.
  print("col1")
  enrichmentlist$InGAP<-enrichmentlist[,1]%in%Background
  #Select only True rows
  print("True")
  dataTrue<-enrichmentlist%>%group_by(Eclass)%>%filter(InGAP==TRUE)
  #Table of Eclass with False and True numbers
  print("Eclass")
  dataTrue<-as.data.frame(table(dataTrue$Eclass))
  #find Eclass list with at least X (overlap) True probes.
  print("Eclass2")
  smalllist<-droplevels(dataTrue[dataTrue$Freq>=overlap,1])
  print("number of categories at end")
  print(length(smalllist))
  #Make reduced list
  return(enrichmentlist[enrichmentlist$Eclass%in%smalllist,1:2])
}


##Function for adding column with list type. Only works with data input I create myself. Userlistenrichment internal lists probably wont work. 
New_Types_danlists<-function(Data){
  #turn to factor
  Data$UserDefinedCategories<-as.character(Data$UserDefinedCategories)
  #split out category
  splitD<-strsplit(Data$UserDefinedCategories,"__")
  new_Category<-unlist(lapply(splitD, `[[`, 1))
  new_Type<-unlist(lapply(splitD, `[[`, 2))
  Data$UserDefinedCategories<-new_Category
  Data$Type<-new_Type
  return(Data)
}

#Round function
round_df_fun <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


# Enrichment data modification: Adjust to data that is available.  --------

#load background data
all_probes_file<-"FEP vs Control_all_probes_limma_results.tsv"

#all_probes_file<-"FEP_vs_Control_all_probes_limma.tsv"
all_probes<-read.csv(paste(P1_output_dir,all_probes_file,sep=""),sep="\t",header=TRUE)


#load(paste(data_dir,"BrainRegionMarkers.rda",sep=""))
#colnames(BrainRegionMarkers)<-c("Gene_names","Eclass")
#write.table(BrainRegionMarkers,paste(data_dir,"BrainReg_WGCNA.csv",sep=""),sep=",",row.names=F)
#write.table(BloodLists,"Blood_WGCNA.csv",sep="",row.names=F) 
#write.table(BrainLists,"Brain_WGCNA.csv",sep="",row.names=F)
#head(BloodLists)

#define lists
lists<-c("GO_Biological_Process_2015.csv","GO_Cellular_Component_2015.csv","GO_Molecular_Function_2015.csv","Kegg2016.csv","GeneSetsMental_Pirooznia_2016final2.csv","Brain_WGCNA.csv","Blood_WGCNA.csv","BrainReg_WGCNA.csv")

#read data and perform functions and write to output dir.  
for(listname in lists){
  inputdata<-read.csv(paste(data_dir,listname,sep=""))
  listdata<-smaller_enrichmentlists(inputdata,all_probes$TargetID,20)
  write.csv(listdata,file=paste(P2_output_dir,"GAP_reduced_",listname,sep=""),row.names=F)  
}


# User List Enrichement for Diffferentially expressed genes ---------------

#load background data
analysis_names<-c("FEP vs Control","Scz vs Con","OP vs Con","Scz vs OP")


for (a in 1:length(analysis_names)){
  setwd(top_dir)  
  #analysis name
  print(analysis_names[a])
  analysis_name<-analysis_names[a]
  
  #load data
  all_probes_file<-paste(analysis_names[a],"_all_probes_limma_results.tsv",sep="")
  all_probes<-read.csv(paste(P1_output_dir,all_probes_file,sep=""),sep="\t",header=TRUE)
  
  #clean data
  all_probes$Groups<-replace(all_probes$Groups, all_probes$Sig_LogFC_probes=="BACKGROUND", "background")
  all_probes$Groups<-replace(all_probes$Groups, all_probes$Sig_LogFC_probes=="Diffexprs", "FDR_Pass")
  
  #enrichment
  setwd(top_dir)
  setwd(P2_output_dir)
  enrichments = userListEnrichment(all_probes$TargetID,all_probes$Groups,fnIn=c("GAP_reduced_GO_Biological_Process_2015.csv","GAP_reduced_GO_Cellular_Component_2015.csv","GAP_reduced_GO_Molecular_Function_2015.csv","GAP_reduced_Kegg2016.csv","GAP_reduced_GeneSetsMental_Pirooznia_2016final2.csv","GAP_reduced_Blood_WGCNA.csv","GAP_reduced_Brain_WGCNA.csv","GAP_reduced_BrainReg_WGCNA.csv"),catNmIn=c("GO_BP","GO_CC","GO_MF","KEGG_2016","Pirooznia","Blood","Brain","BrainReg"),minGenesInCategory = 5)
  setwd(top_dir)
  
  #Extract categories
  enpv<-enrichments$pValue
  enrichment_probes<-unlist(lapply(enrichments$ovGenes,paste, collapse =";"))
  enpv$Genes<-enrichment_probes
  enpv<-enpv[order(enpv$InputCategories,enpv$Pvalues),]
  
  #data table and subset to overlap and p-val
  enpvDT<-data.table(enpv)
  enpvDT<-enpvDT[Pvalues < 0.05 & NumOverlap > 5,.SD[],by=InputCategories]
  enpv_DF<-as.data.frame(enpvDT,stringsAsFactors=FALSE)
  
  #Cange types
  enpv_DF_type<-New_Types_danlists(enpv_DF)
  
  #round
  enpv_DF_type_r<-as.data.frame(round_df_fun(enpv_DF_type,4))
  
  #Save
  write.csv(enpv_DF_type_r,file=paste(P2_output_dir,"1_out_",analysis_name,"_diff_expression_results.csv",sep=""),row.names=F)

}  
  
  
  

  
