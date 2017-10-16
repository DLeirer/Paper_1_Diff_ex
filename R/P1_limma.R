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



# Set directories ---------------------------------------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()
data_dir <-"./P0_Characterise/output/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P1_figs_dir <-"./P1_Diff_Ex/figs/"


# Load eset ---------------------------------------------------------------

load(paste(data_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))



# Linear model to remove Age Ethnicity Sex and Blood differences: ---------
pData_rAESB<-pData(eset_bg_log2_rsn_SVA_Good)
exprs_rAESB<-exprs(eset_bg_log2_rsn_SVA_Good)

mod = model.matrix(~Tc+neutro+as.factor(Sex)+Age+as.factor(Ethnicity),data=pData_rAESB)
fit = lm.fit(mod,t(exprs_rAESB))
#add residuals to average expression data
exprs_res_cor_adj<-t(fit$residuals)+apply(exprs_rAESB,1,mean)
exprs_res_cor_adj[1:10,1:10]
exprs_rAESB[1:10,1:10]

all.equal(colnames(exprs_res_cor_adj),pData_rAESB)

#new("LumiBatch", exprs = exprs_into_lumibatch, se.exprs = [matrix], beadNum = [matrix], detection = [matrix], phenoData = [AnnotatedDataFrame], history = [data.frame], ...)



# Limma case control ------------------------------------------------------

GX_lumi_object<-eset_linear_adj
exprs_rAESB
pData_rAESB
names(pData_rAESB)
p_Phenotype<-factor(pData_rAESB$Phenotype)
p_Gender<-factor(pData_rAESB$Sex)
p_Age<-pData_rAESB$Age
p_Ethnicity<-factor(pData_rAESB$Ethnicity)
p_Tc<-pData_rAESB$Tc
p_neutro<-pData_rAESB$neutro

design = model.matrix(~0+p_Phenotype+p_Gender+p_Age+p_Ethnicity+p_Tc+p_neutro)  


fit <- lmFit(exprs_rAESB, design)

## Limma makeContrasts 
contrasts <- makeContrasts(p_PhenotypeFEP-p_PhenotypeControl, levels=design)

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



#Add fData probe names etc.
fData_rAESB<-fData(eset_bg_log2_rsn_SVA_Good)
fData_rAESB_order<-fData_rAESB[rownames(top_de_genes),]

#check equal
all.equal(rownames(top_de_genes),fData_rAESB_order$nuID)
#if equal combine dfs
top_de_genes_fdata<-cbind(top_de_genes,fData_rAESB_order)
#check equal again
all.equal(rownames(top_de_genes_fdata),top_de_genes_fdata$nuID)

#back to top_de_genes
top_de_genes<-top_de_genes_fdata

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

