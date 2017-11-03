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
full_data_dir <-"./data/"
data_dir <-"./data/Enrichment_libaries/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P2_output_dir <-"./P2_Enrichment/output/"
P4_output_dir <-"./P4_Correlation/output/"



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




# Load Feature data -------------------------------------------------------------

# Full background
lumidata<-"GAP_FEP_Full_Gene_Expression_Data_Linear.RData"
load(file=paste(full_data_dir,lumidata,sep=""))
#Get fdata
Feature_data_full<-fData(eset_bg_log2_rsn_SVA)

#Reduce Feature Data
Full_probe_list<-Feature_data_full[!duplicated(Feature_data_full$TargetID),] %>% 
      mutate(PROBE_KEEP=ifelse( grepl("^LOC",TargetID),"DROP",ifelse( grepl("^HS\\.",TargetID), "DROP","KEEP"))) %>% 
      filter(PROBE_KEEP == "KEEP")



# Enrichment data modification: Adjust to data that is available.  --------


#load background data
#all_probes_file<-"FEP vs Control_all_probes_limma_results.tsv"
#all_probes<-read.csv(paste(P1_output_dir,all_probes_file,sep=""),sep="\t",header=TRUE)
#all_probes<-Feature_data_full
#load(paste(data_dir,"BrainRegionMarkers.rda",sep=""))
#colnames(BrainRegionMarkers)<-c("Gene_names","Eclass")
#write.table(BrainRegionMarkers,paste(data_dir,"BrainReg_WGCNA.csv",sep=""),sep=",",row.names=F)
#write.table(BloodLists,"Blood_WGCNA.csv",sep="",row.names=F) 
#write.table(BrainLists,"Brain_WGCNA.csv",sep="",row.names=F)


setwd(top_dir)
#define lists
lists<-c("GO_Biological_Process_2015.csv","GO_Cellular_Component_2015.csv","GO_Molecular_Function_2015.csv","Kegg2016.csv","GeneSetsMental_Pirooznia_2016final2.csv","Brain_WGCNA.csv","Blood_WGCNA.csv","BrainReg_WGCNA.csv")

#read data and perform functions and write to output dir.  
for(listname in lists){
  inputdata<-read.csv(paste(data_dir,listname,sep=""))
  listdata<-smaller_enrichmentlists(inputdata,Full_probe_list$TargetID,20)
  write.csv(listdata,file=paste(P2_output_dir,"GAP_reduced_",listname,sep=""),row.names=F)  
}


# User List Enrichement for Diffferentially expressed genes ---------------

#load background data
analysis_names<-c("FEP vs Control","Scz vs Con","OP vs Con")
analysis_names_save<-c("FEP_vs_Control","Scz_vs_Con","OP_vs_Con")

#define IDS and files
category_ids<-c("GO_BP","GO_CC","GO_MF","KEGG_2016","Pirooznia","Blood","Brain","BrainReg")
category_file_names<-c("GAP_reduced_GO_Biological_Process_2015.csv","GAP_reduced_GO_Cellular_Component_2015.csv","GAP_reduced_GO_Molecular_Function_2015.csv","GAP_reduced_Kegg2016.csv","GAP_reduced_GeneSetsMental_Pirooznia_2016final2.csv","GAP_reduced_Blood_WGCNA.csv","GAP_reduced_Brain_WGCNA.csv","GAP_reduced_BrainReg_WGCNA.csv")
#select lists to tests
cat_id=c(1:8)


for (a in 1:length(analysis_names)){
  
  setwd(top_dir)  
  #analysis name
  print(analysis_names[a])
  analysis_name<-analysis_names[a]
  
  
  #load data from Limma results
  all_probes_file<-paste(analysis_names[a],"_all_probes_limma_results.tsv",sep="")
  all_probes<-read.csv(paste(P1_output_dir,all_probes_file,sep=""),sep="\t",header=TRUE)
  
  #Filter to diff exprs
  all_probes_reduced<-filter(all_probes,Sig_LogFC_probes == "Diffexprs")
  
  #Make copy of Full probe list
  FData_reduced<-Full_probe_list
  
  #adjust full probe list adding background using diff ex probes
  FData_reduced$Groups<-FData_reduced$TargetID%in%as.character(all_probes_reduced$TargetID)
  FData_reduced$Groups<-replace(FData_reduced$Groups,FData_reduced$Groups == "FALSE", "background")
  FData_reduced$Groups<-replace(FData_reduced$Groups,FData_reduced$Groups == "TRUE", "FDR_Pass")        

  #enrichment
  setwd(top_dir)
  setwd(P2_output_dir)
  enrichments = userListEnrichment(FData_reduced$TargetID,FData_reduced$Groups,fnIn=category_file_names[cat_id],catNmIn=category_ids[cat_id],minGenesInCategory = 10)
  setwd(top_dir)  
  
  
  #Extract categories
  enpv<-enrichments$pValue
  enrichment_probes<-unlist(lapply(enrichments$ovGenes,paste, collapse =";"))
  enpv$Genes<-enrichment_probes
  enpv<-enpv[order(enpv$InputCategories,enpv$CorrectedPvalues),]

  #data table and subset to overlap and p-val
  enpvDT<-data.table(enpv)
  enpvDT<-enpvDT[CorrectedPvalues < 0.05 & NumOverlap > 5,.SD[],by=InputCategories]
  enpv_DF<-as.data.frame(enpvDT,stringsAsFactors=FALSE)
  
  #Change types
  enpv_DF_type<-New_Types_danlists(enpv_DF)
  
  #round
  enpv_DF_type_r<-as.data.frame(round_df_fun(enpv_DF_type,4))
  
  #Save
  write.csv(enpv_DF_type_r,file=paste(P2_output_dir,"1_out_",analysis_names_save[a],"_diff_expression_results.csv",sep=""),row.names=F)
  
}  



  
  



# Correlation results PANSS enrichment ------------------------------------


load(file=paste(P4_output_dir,"correlation_results.Rdata",sep=""))




#load background data
analysis_names_save<-c("FEP_vs_Control","Scz_vs_Con","OP_vs_Con")
analysis_names<-c("PanssScore","PanssPositive","PanssNegative","PRS_0.1_adj","PanssPsycho")
#analysis_names_save<-c("OP_vs_Con")
#analysis_names<-c("PRS_0.1_adj")

#define IDS and files
category_ids<-c("GO_BP","GO_CC","GO_MF","KEGG_2016","Pirooznia","Blood","Brain","BrainReg")
category_file_names<-c("GAP_reduced_GO_Biological_Process_2015.csv","GAP_reduced_GO_Cellular_Component_2015.csv","GAP_reduced_GO_Molecular_Function_2015.csv","GAP_reduced_Kegg2016.csv","GAP_reduced_GeneSetsMental_Pirooznia_2016final2.csv","GAP_reduced_Blood_WGCNA.csv","GAP_reduced_Brain_WGCNA.csv","GAP_reduced_BrainReg_WGCNA.csv")
#select lists to tests
cat_id=c(1:8)

str(Final_results_list_gene)
enrichment_results<-list()
for (i in 1:length(analysis_names_save)){
  analysis_name_save<-analysis_names_save[i]
  temp_list<-Final_results_list_gene[[analysis_name_save]]
  #enrichment_list
  enrichment_var_list_temp<-list()
  
  for (a in 1:length(analysis_names)){

    
    setwd(top_dir)  
    #analysis name
    analysis_name<-analysis_names[a]
    print(paste("Analysis:",analysis_name_save,"Variable:",analysis_name,sep=" "))

    
    
    #Get Data from List
    corr_probes<-temp_list[[analysis_name]]
    #print to check
    print(temp_list[analysis_name])
    #Make copy of Full probe list
    FData_reduced<-Full_probe_list
    
    #adjust full probe list adding background using diff ex probes
    FData_reduced$Groups<-FData_reduced$TargetID%in%corr_probes
    FData_reduced$Groups<-replace(FData_reduced$Groups,FData_reduced$Groups == "FALSE", "background")
    FData_reduced$Groups<-replace(FData_reduced$Groups,FData_reduced$Groups == "TRUE", "FDR_Pass")        
    
    #enrichment
    setwd(top_dir)
    setwd(P2_output_dir)
    enrichments = userListEnrichment(FData_reduced$TargetID,FData_reduced$Groups,fnIn=category_file_names[cat_id],catNmIn=category_ids[cat_id],minGenesInCategory = 10)
    setwd(top_dir)  
    
    
    #Extract categories
    enpv<-enrichments$pValue
    enrichment_probes<-unlist(lapply(enrichments$ovGenes,paste, collapse =";"))
    enpv$Genes<-enrichment_probes
    enpv<-enpv[order(enpv$InputCategories,enpv$CorrectedPvalues),]
    
    
    
    #data table and subset to overlap and p-val
    enpvDT<-data.table(enpv)
    enpvDT<-enpvDT[CorrectedPvalues < 0.05 & NumOverlap > 5,.SD[],by=InputCategories]
    enpv_DF<-as.data.frame(enpvDT,stringsAsFactors=FALSE)
    
    #Change types
    enpv_DF_type<-New_Types_danlists(enpv_DF)
    
    #round
    enpv_DF_type_r<-as.data.frame(round_df_fun(enpv_DF_type,4))
    

    enrichment_var_list_temp[[analysis_name]]<-enpv_DF_type_r
    #Save
    write.csv(enpv_DF_type_r,file=paste(P4_output_dir,"enrichment",analysis_name_save,analysis_name,"_diff_expression_results.csv",sep=""),row.names=F)
    
  }    

  enrichment_results[[analysis_name_save]]<-enrichment_var_list_temp
}

enrichment_results


# Make Enrichment Table Meta ----------------------------------------------


# Load data FEP vs Con, SCZ vs con, OP vs COn, PANSS * 3
#load background data
analysis_names<-c("FEP_vs_Control","Scz_vs_Con","OP_vs_Con")
Panss_names<-c("PanssScore","PanssPositive","PanssNegative")


# dput(names(FEP_PosPan))
#DFcolnames<-c("InputCategories", "UserDefinedCategories", "Type", "NumOverlap","Pvalues", "CorrectedPvalues", "Genes")
# for loop for the range of documents
Enrich_Meta_DF<-data.frame(InputCategories = character(),
                           UserDefinedCategories = character(),
                           Type = character(),
                           NumOverlap = numeric(),
                           Pvalues = numeric(),
                           CorrectedPvalues = numeric(),
                           Genes = character())
Enrich_Meta_DF[1:4]
for (i in 1:3){
  eTempDF<-read.csv(paste(P2_output_dir,"1_out_",analysis_names[i],"_diff_expression_results.csv",sep=""),stringsAsFactors = FALSE)
  eTempDF$InputCategories<-paste(analysis_names[i])
  Enrich_Meta_DF <- rbind(Enrich_Meta_DF, eTempDF)
}


for (i in 1:2){
  for(v in 1:3){
    eTempDF<-read.csv(paste(P4_output_dir,"enrichment",analysis_names[i],variable_names[v],"_diff_expression_results.csv",sep=""),stringsAsFactors = FALSE)
    eTempDF$InputCategories<-paste(analysis_names[i],variable_names[v],sep="")
    Enrich_Meta_DF <- rbind(Enrich_Meta_DF, eTempDF)
  }
}

str(Enrich_Meta_DF)
table(Enrich_Meta_DF[,1])

#Make table for table
Paper_table_colnames<-c("Module","Library","Enriched Categories")
Table_paper <- data.frame(matrix(NA,ncol=3))
names(Table_paper)<-Paper_table_colnames
rowindex = 1
Table_paper
#group by color, > group by library >
#First for loop selects module color. The second for loop selects 
for (module in unique(Enrich_Meta_DF$InputCategories)){
  modulesub<-filter(Enrich_Meta_DF,InputCategories==module)
  print(module)
  #loop starts for module through all categories. Categories being GO, KEGG, Piroznia, ETC.
  for(Type_cat in unique(modulesub$Type)){
    typesub<-filter(modulesub,Type==Type_cat)
    Table_paper[rowindex,1:2] <- typesub[1,c(1,3)]
    Table_paper[rowindex,3]<-paste(head(typesub$UserDefinedCategories), sep="", collapse="; ")
    rowindex = rowindex + 1
  }
  
}

Table_paper[,1:2]
tablename<-"Meta_"
write.csv(Table_paper,file=paste(P2_output_dir,tablename,"Enrichment_results.csv",sep=""),row.names=F)

