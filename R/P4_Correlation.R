# Load Libaries and clean -------------------------------------------------

##cleanup
rm(list = ls())
dev.off()
##Libaries
library(lumi)
library(Biobase)
library(tidyr)
library(dplyr)
library(Hmisc)
library(reshape)
#library(flashClust)
library(ComplexHeatmap)
library(circlize)




# Functions ---------------------------------------------------------------

# list of Probes function
load_probe_lists_fun<-function(analysis_names,analysis_names_save){
  list_probes<-list()
  for (a in 1:length(analysis_names)){
    #analysis name
    print(analysis_names[a])
    
    #load data
    all_probes_file<-paste(analysis_names[a],"_all_probes_limma_results.tsv",sep="")
    all_probes<-read.csv(paste(P1_output_dir,all_probes_file,sep=""),sep="\t",header=TRUE)
    
    #subset to vector
    print("dimension of data")
    print(dim(all_probes))
    list_probes[[analysis_names_save[a]]]<-as.character(filter(all_probes,Sig_LogFC_probes == "Diffexprs")$TargetID)
    
  }
  return(list_probes)
}

#Function for correlation of gene expressiond data.
corr_GX_fun<-function(input_Data,Variable,Probe_index=c(41:4770)){
  
  # check if Character or Integer 
  if (class(Probe_index)=="character"){
    Reduced_Data<-input_Data[,c(Variable,Probe_index)]  
  }else if (class(Probe_index)=="integer"){
    Reduced_Data<-input_Data[,c(Variable,colnames(input_Data)[Probe_index])]  
  }
  
  #return(Reduced_Data)
  return(results<-rcorr(as.matrix(Reduced_Data)))
  
}



# Set directories ---------------------------------------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()
top_dir <-getwd()
data_dir <-"./P0_Characterise/output/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P4_output_dir <-"./P4_Correlation/output/"




# Load Data ---------------------------------------------------------------


#Columns 1-40 are Demographics the rest is Gene Expression. Adjusted for Sex age ethnicity, Cellmix
load(paste(data_dir,"GX_DF_adj_data.Rdata",sep=""))


#Cor PRS, White

#Cor PANSS subscales




# Correlation -------------------------------------------------------------




## Load diff ex list: FEP-Con, SCZ-Con, OP-Con
# Make list
analysis_names<-c("FEP vs Control","Scz vs Con","OP vs Con")
analysis_names_save<-c("FEP_vs_Control","Scz_vs_Con","OP_vs_Con")

# Create sample index
sample_index<-list()
sample_index[[analysis_names_save[1]]]<-GX_DF_adj[,c("sampleID","ICD_DSM")]
sample_index[[analysis_names_save[2]]]<-filter(GX_DF_adj,ICD_DSM != "Other_Psychosis")[,c("sampleID","ICD_DSM")]
sample_index[[analysis_names_save[3]]]<-filter(GX_DF_adj,ICD_DSM != "Schizophrenia")[,c("sampleID","ICD_DSM")]




# dput(names(GX_DF_adj[,1:50])) # Variables to use in correlation.
Variables_for_corr<-c("BMI","PanssScore", "PanssPositive", "PanssNegative", 
                      "PanssPsycho","PRS_5e08_adj", "PRS_1e05_adj", "PRS_1e04_adj", "PRS_0.001_adj", 
                      "PRS_0.01_adj", "PRS_0.05_adj", "PRS_0.1_adj", "PRS_0.2_adj", 
                      "PRS_0.5_adj", "PRS_1_adj")


# Get list of diff ex probes from limma
index_of_probes<-load_probe_lists_fun(analysis_names,analysis_names_save)


# define output list
output<-list()
# loop through analysis
for (a2 in 1:length(analysis_names)){

  # get sample IDs
  sample_ids<-sample_index[[analysis_names_save[a2]]]$sampleID  

  # get probe IDs
  probes_temp<-index_of_probes[[analysis_names_save[a2]]]
  
  
  # analysis name and sample size output
  print(paste("Anlysis name = ",analysis_names[a2]," | number of samples = ",length(sample_ids)," | number of probes = ",length(probes_temp),sep=""))
  
  # subset GX data using sample ids
  GX_DF_temp<-GX_DF_adj[sample_ids,]
  
  # correlation function called
  output[[analysis_names_save[a2]]]<-corr_GX_fun(input_Data=GX_DF_temp,Variable=Variables_for_corr,Probe_index=probes_temp)
  
  #Save
  
  #write.csv(enpv_DF_type_r,file=paste(P2_output_dir,"1_out_",analysis_name,"_diff_expression_results.csv",sep=""),row.names=F)
}  

#check output
output$FEP_vs_Control$n[c(1:20),c(200:210)]
output$Scz_vs_Con$n[c(10:20),c(200:210)]
output$OP_vs_Con$n[c(10:20),c(200:210)]




# Get interesting data ----------------------------------------------------


#DF cols GeneSymbol, Cor, N, Pval


temp_output<-output[temp_name]$FEP_vs_Control

Final_results_list_full<-list()
Final_results_list_gene<-list()
for (a4 in 1:length(analysis_names_save)){
  
  #define analysis
  temp_name<-analysis_names_save[a4]
  print(temp_name)
  #get df
  temp_output<-output[[temp_name]]
  
  #prep list
  results_by_var <- list()
  gene_list <- list()

  #start loop for vars
  for (i in 1:length(Variables_for_corr)){
    #ignore rows
    row_x=c(-1:-15)
    #select column
    col_x=i
    df1<-as.data.frame(cbind(temp_output$r[row_x,col_x],temp_output$n[row_x,col_x],temp_output$P[row_x,col_x]))
    df1$GeneSymbol<-rownames(df1)
    colnames(df1) <-c("cor","N","p_val","GeneSymbol")
    #df1<-filter(df1,p_val < 0.05)
    df1 <- df1 %>% filter(p_val < 0.05) %>% arrange(p_val)
    print(paste ("Significant Probes for ",Variables_for_corr[i]," = ", length(df1$p_val),sep=""))
    results_by_var[[Variables_for_corr[i]]] <- df1
    gene_list[[Variables_for_corr[i]]]<-df1$GeneSymbol
  }

  Final_results_list_full[[temp_name]]<-results_by_var
  Final_results_list_gene[[temp_name]]<-gene_list
  print(paste(temp_name," Done",sep=""))
}
save(Final_results_list_full,Final_results_list_gene,Variables_for_corr,file=paste(P4_output_dir,"correlation_results.Rdata",sep=""), compress = T)






# Correlation + Heatmap Scripts ---------------------------------------------




getwd()
wgcna_dir = "/home/daniel/Documents/A_Year_1_PhD/PhD_Final_Year_projects_2016_to_2017/2016_11_03_Diff_Expression_Paper_clean/P2_WGCNA/output/"
WGCNA_data<-read.csv(paste(wgcna_dir,"Supplementary_table_6_WGCNA_geneInfo.csv",sep=""))
data_dir = "./data/"
load(paste(data_dir,"GX_DF_adj_data.Rdata",sep=""))


GreenYellow_WGCNA<-filter(WGCNA_data,moduleColor == "greenyellow")
GX_DF_adj$PanssPositive
GX_DF_adj$PanssNegative
?rcorr
names()
str(GX_DF_adj[14:47])
Filer_GX_DF<-filter(GX_DF_adj,Medication =="Olanzapine")


str(GX_DF_adj[14:47])



results<-rcorr(as.matrix(Filer_GX_DF[14:4767]))


dim(Filer_GX_DF)
################################ Correlations
results<-rcorr(as.matrix(Filer_GX_DF[14:4767]))
cor_allPanss<-as.data.frame(results$r)
pval_allPanss<-as.data.frame(results$P)
table(cor_allPanss$PanssPositive > 0.3)
Filer_GX_DF[14:20]
table(WGCNA_data[4])
#positive
panss_pos_topprobes[c(20:100),1:4]
panss_pos_topprobes<-cor_allPanss[abs(cor_allPanss$PanssPositive) >0.3, c(1:4,21)]
panss_pos_wgcna<-WGCNA_data[WGCNA_data$nuID%in%rownames(panss_pos_topprobes)[],]
table(panss_pos_wgcna[4])
dim(panss_pos_wgcna)
#negative
panss_neg_topprobes<-cor_allPanss[abs(cor_allPanss$PanssNegative)>0.7, 1:4]
panss_neg_wgcna<-WGCNA_data[WGCNA_data$nuID%in%rownames(panss_neg_topprobes)[-c(1:3)],]
table(panss_neg_wgcna[4])
dim(panss_neg_wgcna)
#PRS
panss_PRS_topprobes<-cor_allPanss[abs(cor_allPanss$PRS_0.1_adj) >0.17,c(1:3, 21)]
panss_PRS_wgcna<-WGCNA_data[WGCNA_data$nuID%in%rownames(panss_PRS_topprobes)[-c(1:3)],]

pval_allPanss[rownames(pval_allPanss)%in%rownames(panss_PRS_topprobes),c(1:3, 21)]


table(panss_PRS_wgcna[4])
dim(panss_PRS_wgcna[1:4])
#combined
panss_neg_wgcna[panss_neg_wgcna$nuID%in%panss_PRS_wgcna$nuID,1]
panss_neg_wgcna[panss_neg_wgcna$nuID%in%panss_pos_wgcna$nuID,1]
panss_pos_wgcna[panss_pos_wgcna$nuID%in%panss_PRS_wgcna$nuID,1]


library(Hmisc)
rcorr(x, type="pearson") 

Probes<-colnames(GX_DF_adj[,c(1600:1633,1750:1875)])
Probes<-colnames(GX_DF_adj[,-c(1:37)])
Probes<-as.character(GreenYellow_WGCNA$nuID)
#Probes<-c(as.character(GreenYellow_WGCNA$nuID),"RBCK1","RNF31")
#colnames(GX_DF_adj)
#select probes from ML
#Probes<-filter(glm_var_list,Model=="2_Gx_Scz")[1:20,2]
#Probes<-glm_var_list[1:20,2]
#Probes<-glm_var_list[1:15,2]
#Probes<-GX_DF_adj[,Probes]
#scale and modifiy
#Probes_scale<-Probes
##########################################
#positive Symptoms Heatmap

PhenoProbes2<-GX_DF_adj
PhenoProbes2$panssposcat<-PhenoProbes2$PanssPositive
table(PhenoProbes2$panssposcat)
PhenoProbes2$panssposcat[GX_DF_adj$Phenotype == "Control"]<-"Control"
PhenoProbes2$panssposcat[GX_DF_adj$PanssPositive >= 20]<-"Panss3_20plus"
PhenoProbes2$panssposcat[GX_DF_adj$PanssPositive >= 15 & GX_DF_adj$PanssPositive < 20]<-"Panss2_15to20"
PhenoProbes2$panssposcat[GX_DF_adj$PanssPositive  < 15]<-"Panss1_below_15"
PhenoProbes2<-PhenoProbes2[!is.na(PhenoProbes2$panssposcat),]


Probes_scale<-scale(PhenoProbes2[,Probes], center = TRUE, scale = TRUE)
Probes_scale<-as.data.frame(Probes_scale)
Probes_scale$sampleID<-rownames(Probes_scale)

PhenoProbes2<-PhenoProbes2[c(1,length(PhenoProbes2))]
#combine data
PhenoProbes2<-inner_join(PhenoProbes2,Probes_scale, by = "sampleID")


dim(PhenoProbes2)
categories<-unique(PhenoProbes2$panssposcat)
output_listp<-list()
for(i in 1:length(categories)){
  group_temp<-categories[i]
  temp_DF<-PhenoProbes2[PhenoProbes2$panssposcat==group_temp,]
  temp_matrix<-as.data.frame(apply(temp_DF[,3:length(names(temp_DF))],2,mean))
  colnames(temp_matrix)<-group_temp
  #rownames(temp_matrix)<-group_temp
  output_listp[[group_temp]] <-temp_matrix
  #final_matrix<-rbind(final_matrix,temp_matrix)
}
new_matrix<-data.frame(output_listp)
new_matrix<-new_matrix[,c(1,2,4,3)]

Heatmap(new_matrix[abs(new_matrix$Control+new_matrix$Panss3_20plus)>0.3,], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = FALSE,row_hclust_side = "right")
Heatmap(new_matrix, name = "Expression", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right",km=3)
Heatmap(new_matrix, name = "Expression", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right")

Heatmap(new_matrix[,c(1,2,4,3)], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right",km = 2)

?Heatmap
colnames(new_matrix) = c("Control","Low","Medium","High")
Heatmap(new_matrix[c(40:60),c(1,7,6,2:5)], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = FALSE,row_hclust_side = "right")
?Heatmap
dim(expr)
dim(mat)
heatmap_filename = paste("/home/daniel/Documents/A_Year_1_PhD/PhD_Final_Year_projects_2016_to_2017/2016_11_03_Diff_Expression_Paper_clean/","GreenYellow","_Panss_Positive_Gx_Heatmap.jpeg",sep="")
jpeg(file=heatmap_filename,width=800,height=1000,pointsize = 20)
Heatmap(new_matrix, name = "Expression", row_title="Probes in Greenyellow Module",column_title="Grouped by Positive Symptom Severity",column_title_side="top",cluster_rows = T, cluster_columns = F,row_dend_side = "right",column_title_rot=0)

dev.off()
##########################################
#Negative Symptoms Heatmap
PhenoProbes2<-GX_DF_adj
PhenoProbes2$panssposcat<-PhenoProbes2$panssposcat
table(PhenoProbes2$panssposcat)
PhenoProbes2$panssposcat[GX_DF_adj$Phenotype == "Control"]<-"Control"
PhenoProbes2$panssposcat[GX_DF_adj$PanssNegative >= 20]<-"Panss3_20plus"
PhenoProbes2$panssposcat[GX_DF_adj$PanssNegative >= 15 & GX_DF_adj$PanssNegative < 20]<-"Panss2_15to20"
PhenoProbes2$panssposcat[GX_DF_adj$PanssNegative  < 15]<-"Panss1_below_15"
PhenoProbes2<-PhenoProbes2[!is.na(PhenoProbes2$panssposcat),]


Probes_scale<-scale(PhenoProbes2[,Probes], center = TRUE, scale = TRUE)
Probes_scale<-as.data.frame(Probes_scale)
Probes_scale$sampleID<-rownames(Probes_scale)

PhenoProbes2<-PhenoProbes2[c(1,length(PhenoProbes2))]
#combine data
PhenoProbes2<-inner_join(PhenoProbes2,Probes_scale, by = "sampleID")


dim(PhenoProbes2)
categories<-unique(PhenoProbes2$panssposcat)
output_listp<-list()
for(i in 1:length(categories)){
  group_temp<-categories[i]
  temp_DF<-PhenoProbes2[PhenoProbes2$panssposcat==group_temp,]
  temp_matrix<-as.data.frame(apply(temp_DF[,3:length(names(temp_DF))],2,mean))
  colnames(temp_matrix)<-group_temp
  #rownames(temp_matrix)<-group_temp
  output_listp[[group_temp]] <-temp_matrix
  #final_matrix<-rbind(final_matrix,temp_matrix)
}
new_matrix<-data.frame(output_listp)
new_matrix<-new_matrix[,c(1,2,4,3)]

Heatmap(new_matrix[abs(new_matrix$Control+new_matrix$Panss3_20plus)>0.3,], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = FALSE,row_hclust_side = "right")
Heatmap(new_matrix, name = "Expression", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right",km=3)
Heatmap(new_matrix, name = "Expression", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right")
Heatmap(new_matrix[,c(1,2,4,3)], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right",km = 2)

?Heatmap
colnames(new_matrix) = c("Control","Low","Medium","High")
Heatmap(new_matrix[c(40:60),c(1,7,6,2:5)], name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = FALSE,row_hclust_side = "right")
?Heatmap
dim(expr)
dim(mat)
heatmap_filename = paste("/home/daniel/Documents/A_Year_1_PhD/PhD_Final_Year_projects_2016_to_2017/2016_11_03_Diff_Expression_Paper_clean/","GreenYellow","_Panss_Negative_Gx_Heatmap.jpeg",sep="")
jpeg(file=heatmap_filename,width=800,height=1000,pointsize = 20)
Heatmap(new_matrix, name = "Expression", row_title="Probes in Greenyellow Module",column_title="Grouped by Negative Symptom Severity",column_title_side="top",cluster_rows = T, cluster_columns = F,row_dend_side = "right",column_title_rot=0)

dev.off()



######################################
expr = readRDS(paste0(system.file(package = "ComplexHeatmap"), "/extdata/gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))

type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(df = data.frame(type = type))

Heatmap(mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(4, "mm"), 
        show_row_names = FALSE, show_column_names = FALSE) +
  Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm")) +
  #Heatmap(expr$length, name = "length", col = colorRamp2(c(0, 1000000), c("white", "orange")),
  #        heatmap_legend_param = list(at = c(0, 200000, 400000, 60000, 800000, 1000000), 
  #                                    labels = c("0kb", "200kb", "400kb", "600kb", "800kb", "1mb")),
  #        width = unit(5, "mm")) +
  Heatmap(expr$type, name = "type", width = unit(5, "mm"))

#Heatmap(t(PhenoProbes2[,48:length(names(PhenoProbes2))]), name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = F,row_hclust_side = "right")
Probes<-GX_DF_adj[,Probes]
#scale and modifiy
#Probes_scale<-Probes
PhenoProbes2<-PhenoProbes2[!is.na(PhenoProbes2$PanssPositive),]

probe_names=colnames(Probes)
variables<-names(PhenoProbes2[c(1:20,46)])
expr$PanssPositive
expr[,10:20]
mat[1:2,1:10]
mat_scaled[1:2,1:10]
str(PhenoProbes2)
expr<-PhenoProbes2[,c(variables,probe_names)]
expr<-expr[order(expr$PanssPositive),]
rownames(expr)<-PhenoProbes2$gap_id
mat = as.matrix(expr[,probe_names])
base_mean = rowMeans(mat)
mat_scaled = apply(mat, 2, scale)
dim(mat_scaled)
dim(mat)

Genes = 
  Type_Panss = expr[,c(5:7,13,15)]
ha = HeatmapAnnotation(df = data.frame(Type = Type_Panss))
?HeatmapAnnotation

Heatmap(t(mat), name = "expression", km = 3, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(10, "mm"), 
        show_row_names = T, show_column_names = FALSE,cluster_columns = F,cluster_rows = T,row_hclust_side = "left") +
  Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm")) +
  Heatmap(expr$length, name = "length", col = colorRamp2(c(0, 1000000), c("white", "orange")),
          heatmap_legend_param = list(at = c(0, 200000, 400000, 60000, 800000, 1000000), 
                                      labels = c("0kb", "200kb", "400kb", "600kb", "800kb", "1mb")),
          width = unit(5, "mm")) +
  Heatmap(expr$type, name = "type", width = unit(5, "mm"))



jpeg(file=paste(figs_dir,project_id,project_name,"_1_Gx_Heatmap.jpeg",sep=""),width=800,height=800,pointsize = 15)
Heatmap(new_matrix4, name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = FALSE,row_hclust_side = "right")
dev.off()

PhenoProbes2$panssposcat
panssposcat
str(PhenoProbes2)
ht_list = NULL
for(s in levels(factor(PhenoProbes2$panssposcat))) {
  #print(filter(PhenoProbes2[46:50],panssposcat == s))
  df_temp<-filter(PhenoProbes2,panssposcat == s)
  ht_list = ht_list + Heatmap(t(df_temp[,48:length(names(PhenoProbes2))]), name = "Relative Expression Level", row_title="Genes",cluster_rows = T, cluster_columns = T,column_title = s)
}

mat1 = matrix(rnorm(80, 2), 8, 10)
mat1 = rbind(mat1, matrix(rnorm(40, -2), 4, 10))
rownames(mat1) = paste0("R", 1:12)
colnames(mat1) = paste0("C", 1:10)

mat2 = matrix(rnorm(60, 2), 6, 10)
mat2 = rbind(mat2, matrix(rnorm(60, -2), 6, 10))
rownames(mat2) = paste0("R", 1:12)
colnames(mat2) = paste0("C", 1:10)

ht1 = Heatmap(mat1, name = "ht1")
ht2 = Heatmap(mat2, name = "ht2")
class(ht1)


ht1 + ht2


