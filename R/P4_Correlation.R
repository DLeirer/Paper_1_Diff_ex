# Load Libaries and clean -------------------------------------------------

##cleanup
rm(list = ls())
dev.off()
##Libaries
library(lumi)
library(Biobase)
library(tidyr)
library(dplyr)



library(reshape)
library(WGCNA)
library(flashClust)
library(ComplexHeatmap)
library(circlize)


# Functions ---------------------------------------------------------------




# Set directories ---------------------------------------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()
top_dir <-getwd()
data_dir <-"./P0_Characterise/output/"
P1_limma_dir <-"./P1_Diff_Ex/output/"
P4_output_dir <-"./P4_Correlation/output/"




# Load Data ---------------------------------------------------------------
getwd()

load(paste(data_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))



## Load diff ex list: FEP-Con, SCZ-Con, OP-Con
# Make list
GX_DF[1:10,1:10]
## Load Gene ex matrix
GX_DF<-exprs(eset_bg_log2_rsn_SVA_Good)
# Adjust for Cellmix etc.  
# Sex+Age+Ethnicity+Tc+neutro
GX_DF_adj[1:10,1:40]
## Load P-Data PANSS, PRS
Pheno_data<-pData(eset_bg_log2_rsn_SVA_Good)
#Cor PRS, White
#Cor PANSS subscales






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


