## Overview
#Differential Expression on GAP FEP Data. We are removing all batch and covariate effects, (technical variables, age, gender, ethnicity). We can't control effectively for Medication, Tobacco and BMI or anything else. This will have to be done in a posthoc analysis, so these results should be interpreted cautiously.
#That being said we think this is still appropiate and usefull, since First episode psychosis patients will at most only be exposed to medication for a short while, and we regressed statistically significant cell effects out (TC and Neutro) at a previous step. 

## Aim
#Overall we hope to see fundamental gene expression difference between cases and controls, regardless of age, gender or ethnicity. These can be investigated elswhere.
#It is useful to find biological differences independent of the above mentioned phenotypes, since this might elucidate fundamental changes in peripheral blood of psychosis patients that could be incorporated in clinical tests, or taken further to investigate biological underpinnings of the disease. 

# Load Libaries and clean -------------------------------------------------

##cleanup
rm(list = ls())
dev.off()
library(lumi)
library(limma)
library(dplyr)
library(Tmisc)
library(calibrate)
library(ggplot2)



# Functions ---------------------------------------------------------------

# Limma Cleaner Function
Limma_Datatoplot_fun <- function(limma_matrix){
  limma_matrix$logFC <- round(limma_matrix$logFC,2)
  limma_matrix$CI.L <- round(limma_matrix$CI.L,2)
  limma_matrix$CI.R <- round(limma_matrix$CI.R,2)
  limma_matrix$AveExpr <- round(limma_matrix$AveExpr,2)
  limma_matrix$t <- round(limma_matrix$t,2)
  limma_matrix$P.Value <- signif(limma_matrix$P.Value,3)
  limma_matrix$adj.P.Val <- signif(limma_matrix$adj.P.Val,3)
  limma_matrix$B <- round(limma_matrix$B,2)
  small_matrix <- limma_matrix[,c("TargetID","logFC","P.Value","adj.P.Val","CHROMOSOME","DEFINITION")]
  
  small_matrix <- small_matrix %>% 
    #define new columns
    mutate(SIG_DE=adj.P.Val <=0.05, 
           LogFC_DIRECTION=ifelse(logFC >= 0, "up-regulated", ifelse(logFC < 0, "down-regulated", "no-change")),
           LogFC_BIOLOCICAL=ifelse(logFC >= 0.1, "up-regulated", ifelse(logFC <= -0.1, "down-regulated", "no-sig-change")),
           PROBE_KEEP=ifelse( grepl("^LOC",TargetID),"DROP",ifelse( grepl("^HS\\.",TargetID), "DROP","KEEP"))) %>% 
    mutate(Sig_LogFC_probes = ifelse(SIG_DE==TRUE & LogFC_BIOLOCICAL != "no-sig-change","Diffexprs","BACKGROUND"))
  
  
  ## Significantly DE Probes
  return(small_matrix %>% filter(PROBE_KEEP=="KEEP")) 
}

#Volcano Plot function
Volcano_plot_fun <- function(datatoplot,n_up_probes=10,n_down_probes=10,n_sig_probes = 20){
  x_max = max(datatoplot$logFC)+0.5
  x_min = min(datatoplot$logFC)-0.5
  with(datatoplot, plot(logFC, -log10(P.Value), pch=20, main="Volcano Plot of Differentially Expressed Probes", xlim=c(x_min,x_max)))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(datatoplot, logFC >= 0 ), points(logFC, -log10(P.Value), pch=20, col="navy"))
  with(subset(datatoplot, logFC <= 0 ), points(logFC, -log10(P.Value), pch=20, col="darkgreen"))
  with(subset(datatoplot, Sig_LogFC_probes =="BACKGROUND" ), points(logFC, -log10(P.Value), pch=20, col="firebrick"))
  
  # Find top probes
  upreg_logFC<-filter(datatoplot,rank(desc(logFC))<=n_up_probes)$TargetID
  downreg_logFC<-filter(datatoplot,rank(logFC)<=n_down_probes)$TargetID
  probes_adj_p_val<-filter(datatoplot,rank(adj.P.Val)<=n_sig_probes)$TargetID
  # make list of unique probes
  Probes_plot<-unique(c(probes_adj_p_val,downreg_logFC,upreg_logFC))
  # Label points with the textxy function from the calibrate plot
  with(datatoplot[datatoplot$TargetID%in%Probes_plot,], textxy(logFC, -log10(P.Value), labs=TargetID, cex=.5))
}


# Set directories ---------------------------------------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()
data_dir <-"./P0_Characterise/output/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P1_figs_dir <-"./P1_Diff_Ex/figs/"


# Load eset ---------------------------------------------------------------

load(paste(data_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))



# Linear model to remove Age Ethnicity Sex and Blood differences: ---------


# Limma case control ------------------------------------------------------
#get pheno data
pData_rAESB<-pData(eset_bg_log2_rsn_SVA_Good)
names(pData_rAESB)

#make phenodata objects for model matrix
design=model.matrix(~0+Phenotype+Sex+Age+Ethnicity+Tc+neutro,data=pData_rAESB)

#lmFite
fit <- lmFit(eset_bg_log2_rsn_SVA_Good, design)

## Limma makeContrasts 
contrasts <- makeContrasts(PhenotypeFEP-PhenotypeControl, levels=design)

## Limma eBayes on makeContrasts
contrast.fit <- contrasts.fit(fit, contrasts)
contrast.fit <- eBayes(contrast.fit)


# Top Table ---------------------------------------------------------------

top_de_genes <- topTable(contrast.fit, coef=1, number=5000,adjust.method="fdr",p.value=1,confint=TRUE)

head(top_de_genes)



# tidy up
top_de_genes$logFC <- round(top_de_genes$logFC,2)
top_de_genes$CI.L <- round(top_de_genes$CI.L,2)
top_de_genes$CI.R <- round(top_de_genes$CI.R,2)
top_de_genes$AveExpr <- round(top_de_genes$AveExpr,2)
top_de_genes$t <- round(top_de_genes$t,2)
top_de_genes$P.Value <- signif(top_de_genes$P.Value,3)
top_de_genes$adj.P.Val <- signif(top_de_genes$adj.P.Val,3)
top_de_genes$B <- round(top_de_genes$B,2)




##Save full dataframe
full_limma_results_AgeSexEthnicity_Adj <- top_de_genes
write.table(full_limma_results_AgeSexEthnicity_Adj, file=paste(P1_output_dir,"full_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")

##Save Reduced dataframe
small_limma_results_AgeSexEthnicity_Adj <- top_de_genes[,c("TargetID","logFC","P.Value","adj.P.Val","CHROMOSOME","DEFINITION")]
write.table(small_limma_results_AgeSexEthnicity_Adj, file=paste(P1_output_dir,"small_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")



# Significantly DE Probes -------------------------------------------------

de_res <- top_de_genes[,c("TargetID","logFC","P.Value","adj.P.Val","CHROMOSOME","DEFINITION")]

de_res <- de_res %>% 
  mutate(SIG_DE=adj.P.Val <=0.05, 
         LogFC_DIRECTION=ifelse(logFC >= 0, "up-regulated",
                                ifelse(logFC < 0, "down-regulated", "no-change")),
         LogFC_BIOLOCICAL=ifelse(logFC >= 0.1, "up-regulated",
                                 ifelse(logFC <= -0.1, "down-regulated", "no-sig-change")),
         PROBE_KEEP=ifelse( grepl("^LOC",TargetID),"DROP",
                            ifelse( grepl("^HS\\.",TargetID), "DROP","KEEP"))
  ) %>% mutate(Sig_LogFC_probes = ifelse(SIG_DE==TRUE & LogFC_BIOLOCICAL != "no-sig-change","Diffexprs","BACKGROUND"))


## Significantly DE Probes
datatoplot<- de_res %>% filter(PROBE_KEEP=="KEEP")

write.table(datatoplot, file=paste(P1_output_dir,"Datatoplot.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")


## Number of probes
dim(de_res)
dim(datatoplot)
table(datatoplot$PROBE_KEEP)
4730-4063

## Number of significant probes
table(datatoplot$Sig_LogFC_probes)


##All probes
write.table(de_res, file=paste(P1_output_dir,"Supplementary_Table_1_all_probes_4730_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")
#write.table(datatoplot, file=paste(P1_output_dir,"Supplementary_Table_2_all_probes_4063_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")

##Significant probes
#datatoplot_sig<-datatoplot[datatoplot$adj.P.Val <= 0.05,]
#datatoplot_sig<-datatoplot[datatoplot$LogFC_BIOLOCICAL != "no-sig-change" & datatoplot$SIG_DE == TRUE,]
datatoplot_sig<-datatoplot[datatoplot$Sig_LogFC_probes == "Diffexprs",]
dim(datatoplot_sig)
write.table(datatoplot_sig, file=paste(P1_output_dir,"Supplementary_Table_2_significant_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")




##upregulated
datatoplot_sig_up<-datatoplot_sig[datatoplot_sig$logFC >= 0.1,]
dim(datatoplot_sig_up)
write.table(filter(datatoplot_sig_up,Sig_LogFC_probes == "Diffexprs"), file=paste(P1_output_dir,"Supplementary_Table_3_upregulated_significant_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")

table1up<-filter(datatoplot_sig_up,Sig_LogFC_probes == "Diffexprs")[1:50,c("TargetID","CHROMOSOME","DEFINITION","logFC","adj.P.Val")]


table1up$logFC<-signif(table1up$logFC,2)
table1up$adj.P.Val<-signif(table1up$adj.P.Val,2)

write.table(table1up, file=paste(P1_output_dir,"Table_1_upregulated_significant_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")


##downregulated
datatoplot_sig_down<-datatoplot_sig[datatoplot_sig$logFC <= -0.1,]
dim(datatoplot_sig_down)
write.table(datatoplot_sig_down, file=paste(P1_output_dir,"Supplementary_Table_5_downregulated_significant_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")

dim(filter(datatoplot_sig_down,Sig_LogFC_probes == "Diffexprs"))

table2down<-filter(datatoplot_sig_down,Sig_LogFC_probes == "Diffexprs")[1:50,c("TargetID","CHROMOSOME","DEFINITION","logFC","adj.P.Val")]

table2down$logFC<-signif(table2down$logFC,2)
table2down$adj.P.Val<-signif(table2down$adj.P.Val,2)

write.table(table2down, file=paste(P1_output_dir,"Table_2_downregulated_significant_limma_results_AgeSexEthnicity_Adj.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")




# Plots: Volcano and Chromosome -------------------------------------------

# Make a basic volcano plot

filename <- "Volcanoplot.jpeg"
jpeg(file = paste(P1_figs_dir,filename,sep=""), pointsize = 20, width = 1500, height = 1300)
with(datatoplot, plot(logFC, -log10(P.Value), pch=20, main="Volcano Plot of Differentially Expressed Probes", xlim=c(-1,1)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(datatoplot, logFC >= 0 ), points(logFC, -log10(P.Value), pch=20, col="navy"))
with(subset(datatoplot, logFC <= 0 ), points(logFC, -log10(P.Value), pch=20, col="darkgreen"))
with(subset(datatoplot, Sig_LogFC_probes =="BACKGROUND" ), points(logFC, -log10(P.Value), pch=20, col="firebrick"))


# Label points with the textxy function from the calibrate plot
with(datatoplot [c(1,4,8,9,10,12),], textxy(logFC, -log10(P.Value), labs=TargetID, cex=.5))
with(subset(datatoplot, logFC < -0.25 & adj.P.Val <0.002 ), textxy(logFC, -log10(P.Value), labs=TargetID, cex=.5,offset = .6))
with(subset(datatoplot, logFC < -0.33 & adj.P.Val <0.05 ), textxy(logFC, -log10(P.Value), labs=TargetID, cex=.5,offset = .6))
with(subset(datatoplot, logFC >= 0.38 & adj.P.Val <0.05 ), textxy(logFC, -log10(P.Value), labs=TargetID, cex=.5,offset = .6))
dev.off()


#LogFC by chromosome
title = "LogFC by Chromosome"
ggplot(data = datatoplot, 
       aes(x = CHROMOSOME, y=logFC, color=CHROMOSOME)) +
  geom_boxplot(alpha = 0)+
  geom_text(data=filter(datatoplot, logFC >= 0.3 |logFC <= -0.27),check_overlap = TRUE,angle=45,
            aes(CHROMOSOME,logFC,label=TargetID),size=3)+
  ggtitle(title)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 1))
ggsave(paste(P1_figs_dir,title,".png",sep=""))




# Limma Diagnosis ------------------------------------------------------

names(pData_rAESB)

#make model Matrix
design_diagnosis=model.matrix(~0+ICD_DSM+Sex+Age+Ethnicity+Tc+neutro,data=pData_rAESB)

#lmFit
fit_diagnosis <- lmFit(eset_bg_log2_rsn_SVA_Good, design_diagnosis)

## Limma makeContrasts 
contrast.matrix <- makeContrasts(ICD_DSMSchizophrenia-ICD_DSMControl,ICD_DSMOther_Psychosis-ICD_DSMControl,ICD_DSMSchizophrenia-ICD_DSMOther_Psychosis, levels=design_diagnosis)

## Limma eBayes on makeContrasts
fit_diagnosis2 <- contrasts.fit(fit_diagnosis, contrast.matrix)
fit_diagnosis2 <- eBayes(fit_diagnosis2)

results <- decideTests(fit_diagnosis2)
results

VennDia<-"VennDiagram_Compare_diagnosis"
jpeg(paste(P1_figs_dir,VennDia,".jpeg",sep=""))
vennDiagram(results,names=c("Scz vs Con", "OP vs Con", "Sch vs OP"),cex=1)
dev.off()


Dia_Scz_Con<-topTable(fit_diagnosis2, coef=1, number=5000,adjust.method="fdr",p.value=1,confint=TRUE)
topTable(fit_diagnosis2, coef=2, adjust="BH")
topTable(fit_diagnosis2, coef=3, adjust="BH")


top_de_genes <- topTable(contrast.fit, coef=1, number=5000,adjust.method="fdr",p.value=1,confint=TRUE)



Dia_Scz_Con<-Limma_Datatoplot_fun(Dia_Scz_Con)
head(Dia_Scz_Con)



filename <- "Test.jpeg"
jpeg(file = paste(P1_figs_dir,filename,sep=""), pointsize = 20, width = 1500, height = 1300)
Volcano_plot_fun(Dia_Scz_Con)
dev.off()
