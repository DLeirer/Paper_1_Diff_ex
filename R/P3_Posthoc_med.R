

# Load Libraries ----------------------------------------------------------

rm(list=ls())
gc()
library(lumi)
library(limma)
library(plyr)
library(dplyr)
library(Tmisc)
library(calibrate)
library(ggplot2)


#library(ggbiplot)
#library(tableone)


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



#Save results in TSV
save_limma_results_fun<-function(limma_df,analysis_name,save_dir){
  
  ##All probes
  write.table(limma_df, file=paste(save_dir,analysis_name,"_all_probes_limma_results.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")
  
  #Make smaller and clean
  limma_df_small<-filter(limma_df,Sig_LogFC_probes == "Diffexprs")[,c("TargetID","CHROMOSOME","DEFINITION","logFC","adj.P.Val")]
  limma_df_small$logFC<-signif(limma_df_small$logFC,2)
  limma_df_small$adj.P.Val<-signif(limma_df_small$adj.P.Val,2)
  
  ##upregulated
  table_sig_up<-limma_df_small[limma_df_small$logFC >= 0.1,]
  write.table(table_sig_up, file=paste(save_dir,analysis_name,"_upregulated_significant_limma_results.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")
  
  
  ##downregulated
  table_sig_down<-limma_df_small[limma_df_small$logFC <= -0.1,]
  write.table(table_sig_down, file=paste(save_dir,analysis_name,"_downregulated_significant_limma_results.tsv",sep=""),row.names=FALSE,quote=FALSE,sep = "\t")
}


#Volcano Plot function
Volcano_plot_fun <- function(datatoplot,title_comp,n_up_probes=10,n_down_probes=10,n_sig_probes = 20){
  x_max = max(datatoplot$logFC)+0.5
  x_min = min(datatoplot$logFC)-0.5
  with(datatoplot, plot(logFC, -log10(P.Value), pch=20, main=paste("Volcano Plot of",title_comp,"comparison",sep=" "), xlim=c(x_min,x_max)))
  
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

#Chromosome plot function
Chromosome_plot_fun<-function(data_limma,title,logFC_up= 0.3,logFC_down = -0.27){
  ggplot(data = data_limma, 
         aes(x = CHROMOSOME, y=logFC, color=CHROMOSOME)) +
    geom_boxplot(alpha = 0)+
    geom_text(data=filter(data_limma, logFC >= logFC_up |logFC <= logFC_down),check_overlap = TRUE,angle=45,
              aes(CHROMOSOME,logFC,label=TargetID),size=3)+
    ggtitle(title)+ 
    theme_bw(base_size = 10) + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 1))
  
}




# Define Directories ------------------------------------------------------

data_dir <-"./P0_Characterise/output/"
P1_output_dir <-"./P1_Diff_Ex/output/"
P3_output_dir <-"./P3_Posthoc/output/"
P3_figs_dir <-"./P3_Posthoc/figs/"


# Load eset ---------------------------------------------------------------

load(paste(data_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))

#load previous limma results
significant_probes<-read.csv(paste(P1_output_dir,"FEP vs Control_all_probes_limma_results.tsv",sep=""),sep="\t",header=TRUE)



# Modify Phenodata --------------------------------------------------------

Pheno_data_temp<-pData(eset_bg_log2_rsn_SVA_Good)

# 1 Med vs control, Ola vs Control, Ris vs Control, AF vs Control, AF vs Ola, Af vs Ris. 
#change af to antipsychotic_free
#get rid of NA
New_Pheno_Data<-Pheno_data_temp %>% mutate(Medication=replace(Medication, Medication=="Antipsychotic free", "AntipsychoticFree")) %>% mutate(Medication=replace(Medication, is.na(Medication), "Unknown")) 
New_Pheno_Data<-New_Pheno_Data %>% mutate(AllMed=ifelse(Medication == "CONTROL","CONTROL",
                                                        ifelse(Medication == "Unknown" | Medication == "AntipsychoticFree","Notused","Med")))

table(New_Pheno_Data$Medication)
table(New_Pheno_Data$AllMed)


#Maybe do later. 
#filter(significant_probes$Sig_LogFC_probes,)


#this is bad and needs to be changed
l<-eset_bg_log2_rsn_SVA_Good
rownames(l@phenoData@data)<-l@phenoData@data$sampleID
lf<-l@featureData
lf@data<-l@featureData@data[,]
le<-exprs(l)[,]
NewExprs<-ExpressionSet(assayData=le,phenoData=l@phenoData,featureData = lf)


#Back into eset

Modified_eset<-NewExprs
#Modified_eset<-eset_bg_log2_rsn_SVA_Good
pData(Modified_eset)<-New_Pheno_Data
#exprs(Modified_eset)<-exprs(Modified_eset)
Modified_eset
Modified_eset@phenoData@data$Medication


# Limma Medication ------------------------------------------------------



#make model Matrix
design_diagnosis=model.matrix(~0+Medication+Sex+Age+Ethnicity+Tc+neutro,data=New_Pheno_Data)



#lmFit
fit_diagnosis <- lmFit(Modified_eset, design_diagnosis)
#noquote(colnames(fit_diagnosis$design))



## Limma makeContrasts 1 = AF - Con, 2 Ola - Con,  3Ris - Con, 4AF - Ola, AF - Ris
contrast.matrix <- makeContrasts(MedicationAntipsychoticFree-MedicationCONTROL,MedicationOlanzapine-MedicationCONTROL,MedicationRisperidone-MedicationCONTROL,MedicationAntipsychoticFree-MedicationOlanzapine,MedicationAntipsychoticFree-MedicationRisperidone, levels=design_diagnosis)


## Limma eBayes on makeContrasts
fit_diagnosis2 <- contrasts.fit(fit_diagnosis, contrast.matrix)
fit_diagnosis2 <- eBayes(fit_diagnosis2)



# diagnosis outputs (tables and figures) ---------------------------------------------------------------

#vennDiagram plot output
results <- decideTests(fit_diagnosis2)
results

analysis_names = c("AFvsCon", "OlavsCon","RisvsCon","AFvsOla","AFvsRis")
VennDia<-"VennDiagram_Compare_Med"
jpeg(paste(P3_figs_dir,VennDia,".jpeg",sep=""))
vennDiagram(results,names=c(analysis_names),cex=1)
dev.off()


for (i in 1:length(analysis_names)){
  print(analysis_names[i])
  analysis_name_temp<-analysis_names[i]
  
  #Generate top list for each comparision and save
  Top_De_Genes_temp <- topTable(fit_diagnosis2, coef=i, number=5000,adjust.method="fdr",p.value=1,confint=TRUE)
  
  Top_De_Genes_temp<-Limma_Datatoplot_fun(Top_De_Genes_temp)
  #Save Limma Results up down and full.
  save_limma_results_fun(Top_De_Genes_temp,analysis_name_temp,save_dir=P3_output_dir)  
  
  ## Plots
  #Volcano Plot
  filename <- paste(analysis_name_temp,"_Volcanoplot.jpeg",sep="")
  jpeg(file = paste(P3_figs_dir,filename,sep=""), pointsize = 20, width = 1500, height = 1300)
  Volcano_plot_fun(Top_De_Genes_temp,analysis_name_temp,10,10,15)
  dev.off()
  
  
  #Chromosome Plot
  title = paste(analysis_name_temp," LogFC by Chromosome",sep="")
  Chromosome_plot_fun(Top_De_Genes_temp,title = title)
  ggsave(paste(P3_figs_dir,title,".png",sep=""))
  
}
