
#CellMix approach:



# Install and load Libraries: ---------------------------

rm(list = ls())
dev.off()
library(CellMix)
library(lumi)
library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)
library(tableone)
library(stargazer)


# Functions: --------------------------------------------------------------

## Table function 1 demographics
table_dem_fun<-function(pdata,table_name,stratify,listVars,catVars){
  table1 <- CreateTableOne(vars = listVars, data = pdata, factorVars = catVars,strata=c(stratify),includeNA = T)
  table1print<-print(table1)
  table1print<-table1print[,-length(names(data.frame(table1print)))]
  write.csv(table1print, file=paste(P0_output_dir,"no_latex_",table_name,sep=""),row.names = TRUE, col.names = TRUE,quote=FALSE,sep = ",")
  write.csv(stargazer(table1print,summary=FALSE), file=paste(P0_output_dir,"latex_",table_name,sep=""),row.names=F,col.names = F,sep="")
}




# Set directories: ---------------------------

setwd("/home/daniel/Documents/Post_PhD_Papers/Paper_1_Diff_ex")
getwd()

data_dir <-"./data/"
P0_output_dir <-"./P0_Characterise/output/"
P0_figs_dir <-"./P0_Characterise/figs/"






# Load gene expression data (Should be LumiBatch object): ---------------------------

lumidata<-"GAP_FEP_Full_Gene_Expression_Data_Linear.RData"
load(file=paste(data_dir,lumidata,sep=""))
PRS_data<-read.csv(paste(data_dir,"GAP_FEP_Polygenic_Risk_scores.csv",sep=""))
PRS_eigenvalues_data<-read.csv(paste(data_dir,"GAP_eigenvectors_from_Vangelis_04_03_2015.csv",sep=""),na.strings=c("#N/A","notfound"))



# Subset to good probes: ---------------------------

dim(eset_bg_log2_rsn_SVA)
eset_bg_log2_rsn_SVA_Good<-eset_bg_log2_rsn_SVA

#make exprs data with good probes
exprs_data<-exprs(eset_bg_log2_rsn_SVA)
feature_data<-fData(eset_bg_log2_rsn_SVA)

#get good probe
good_fdata<-filter(feature_data, good_probe =="TRUE")

exprs_into_lumibatch<-exprs_data[rownames(exprs_data)%in%good_fdata$nuID,]

# REMOVE DUPLICATES BY SELCTING HIGHEST AVERAGE EXPRESSED PROBE 
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




##GED BLOOD ABBAS WHOLE BLOOD:
# gedBlood is a meta function. It uses a standard set of functions included in CellMix to generate results. It seems to be able to handle nuIDs. gedBlood automatically adjusts gene IDs so I suspect it searches the LumiBatch file for IDs it recognises and matches them with whatever dataset it uses. 
# In this case the ABBAS blood atlas is used. 
# CLsubset allows you to choose between WB for whole blood and PBMCs. 
# verbose simply tells you the steps it takes which would usually have to be performed manually in the package. 
# These settings just use a linear model based on gene expression values for different cell types from the ABBAS blood atlas (Abbas et al. 2009)* to estimate cell proportions in each sample.
# This is after preprocessing, so we only use the approx. 5000 probes that passed QC.
# Cellmix has a lot more functions and ways to estimate cell proportions, but after looking through all of it, you could probably make it into a study by itself. Some of the methods seem to be quite computationally intensive. See table 2 on page 26 [here](http://web.cbio.uct.ac.za/~renaud/CRAN/web/CellMix/vignettes/Introduction.pdf) for list of approaches. 
# 
# ^*Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark HF (2009). "Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus."^
#   




# GEDBlood Function: ---------------------------

res_all <- gedBlood(eset_bg_log2_rsn_SVA_Good, CLsubset = "WB", verbose = TRUE)


#Extract cell proportions from res_all
wbloodprop<-coef(res_all)

#Remove rows with sum of 0. Rows represent cell types. If a cell type has a sum of 0 across all samples I exclude it. 

reduced_props<-wbloodprop[apply(wbloodprop, 1, function(x) !all(x==0)),]





# Check for differences between cases and controls: -----------------------
# Here I split the blood proportion data in case and control so I can plot it using ggplot2. 


#Get Case Control Status
pheno_data <- pData(eset_bg_log2_rsn_SVA_Good)

#transform data
reduced_props<-data.frame(t(reduced_props))
#check samples in correct order
all.equal(rownames(reduced_props),rownames(pheno_data))
#add pheno data
reduced_probs_pheno<-cbind(reduced_props,pheno_data)
names(reduced_probs_pheno)

#Melt Data for ggplot
mwbdata <- melt(reduced_probs_pheno, id=c(11:71))
head(mwbdata)



#Graph data
bloodgraph<-ggplot(mwbdata,aes(x=variable,y=value,fill=Gender)) +  
  geom_boxplot()+
  ggtitle("Blood Cell Proportions")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  xlab("Cell Type")+ylab("Percentage for each individual")
bloodgraph


#select plots to make
plotnames<-names(reduced_probs_pheno[c("Phenotype","Gender","Tobacco","Ethnicity","Medication","tech.Sentrix.Barcode","tech.Date_Washing_and_scanning","tech.Date_Quantitation_by_RiboGreen","tech.SampleSection")])

names(reduced_probs_pheno)

#make function
f <- function(mwbdata, fill_name) {
  bloodgraph<-ggplot(mwbdata,aes_string(x="variable",y="value",fill=fill_name) ) +  
    geom_boxplot()+
    ggtitle("Blood Cell Proportions")+
    theme(plot.title = element_text(lineheight=.8, face="bold"))+
    xlab("Cell Type")+ylab("Percentage for each individual")
  print(bloodgraph)
}

#Create loop for all variables and save
CellMixgraph<-"Cell_Mix_Boxplots.pdf"
pdf(paste(P0_figs_dir,CellMixgraph,sep=""))
for (fill_number in 1:length(plotnames)){
  fill_name<-plotnames[fill_number]
  f(mwbdata,fill_name)  
  
}
dev.off()






# Statistics: -------------------------------------------------------------
names(mwbdata)
stats_pheno<-ddply(mwbdata,"variable",
                   function(x) {
                     w <- wilcox.test(value~Phenotype,data=x)
                     with(w,data.frame(statistic,p.value))
                   })

stats_pheno

# Add PRS, ICD_DSM, TC and neutro columns to phenodata: -----------------------------------------

cor(reduced_probs_pheno[1:10])
#reduced_probs_pheno<-reduced_probs_pheno[,c(11:71,1,10)]
#pData(eset_bg_log2_rsn_SVA_Good)<-reduced_probs_pheno[,c(11:71,1,10)]
#names(pData(eset_bg_log2_rsn_SVA_Good))




## add PRS
reduced_probs_pheno_prs=merge(reduced_probs_pheno,PRS_data[,-c(2:3)],by="sampleID")
names(reduced_probs_pheno_prs)


##Define Columns to keep.
#dput(names(reduced_probs_pheno_prs))
DataColGAP_keep<-c("sampleID","Study_ID","gap_id","Phenotype","SEX", "Age", "Ethnicity", "BMI", "Tobacco", 
                   "Medication", "dsmiv.opcrit", "icd10.opcrit", "panss.date", "PanssScore", 
                   "PanssPositive", "PanssNegative", "PanssPsycho","Tc", "neutro", 
                   "Pol_5e08_GAP_all_strict_excl_WTCCC2", "Pol_1e05_GAP_all_strict_excl_WTCCC2", 
                   "Pol_1e04_GAP_all_strict_excl_WTCCC2", "Pol_0.001_GAP_all_strict_excl_WTCCC2", 
                   "Pol_0.01_GAP_all_strict_excl_WTCCC2", "Pol_0.05_GAP_all_strict_excl_WTCCC2", 
                   "Pol_0.1_GAP_all_strict_excl_WTCCC2", "Pol_0.2_GAP_all_strict_excl_WTCCC2", 
                   "Pol_0.5_GAP_all_strict_excl_WTCCC2", "Pol_1_GAP_all_strict_excl_WTCCC2"
)
names(reduced_probs_pheno_prs[,DataColGAP_keep])

reduced_probs_pheno_order<-reduced_probs_pheno_prs[,DataColGAP_keep]
names(reduced_probs_pheno_order)


##Change colnames
New_col_names<-c("sampleID","Study_ID","gap_id","Phenotype","Sex", "Age", "Ethnicity", "BMI", "Tobacco", 
                   "Medication", "dsmiv", "icd10", "panss_date", "PanssScore", 
                   "PanssPositive", "PanssNegative", "PanssPsycho","Tc", "neutro", 
                   "PRS_5e08", "PRS_1e05", "PRS_1e04", "PRS_0.001", "PRS_0.01", "PRS_0.05", "PRS_0.1", "PRS_0.2", "PRS_0.5", "PRS_1")


colnames(reduced_probs_pheno_order)<-New_col_names
names(reduced_probs_pheno_order)

## add control to dsmiv and icd.

reduced_probs_pheno_order$Phenotype[reduced_probs_pheno_order$Phenotype == "control"] <- "Control"
reduced_probs_pheno_order$dsmiv[reduced_probs_pheno_order$Phenotype == "Control"] <- "Control"
reduced_probs_pheno_order$icd10[reduced_probs_pheno_order$Phenotype == "Control"] <- "Control"

#remove NAs make no criteria met.
reduced_probs_pheno_order$icd10[is.na(reduced_probs_pheno_order$icd10)] <-"No criteria Met" 
reduced_probs_pheno_order$dsmiv[is.na(reduced_probs_pheno_order$dsmiv)] <-"No criteria Met" 
reduced_probs_pheno_order$dsmiv[reduced_probs_pheno_order$dsmiv == "No Criteria Met"] <-"No criteria Met" 

reduced_pheno_final<-reduced_probs_pheno_order %>% mutate(ICD_DSM=ifelse(dsmiv=="Schizophrenia" | icd10 == "Schizophrenia","Schizophrenia",
                                                                             ifelse(dsmiv =="Control", "Control", "Other_Psychosis")))

##Final order
reduced_pheno_final<-reduced_pheno_final[,c(1:12,30,13:29)]
table(reduced_pheno_final$icd10)
table(reduced_pheno_final$dsmiv)
table(reduced_pheno_final$ICD_DSM)


# Adjust PRS for Principle components -------------------------------------

#New matrix Remove mising values
names(reduced_pheno_final)
head(reduced_pheno_final)
PRS_to_adjust<-reduced_pheno_final[c(1:3,21:30)]
rownames(PRS_to_adjust)<-reduced_pheno_final$sampleID

#remove columns not needed
PRS_to_adjust<-PRS_to_adjust[complete.cases(PRS_to_adjust),]


##get principle components to use for adjustment
#Reduce
PRS_eigenvalue_adjust<-PRS_eigenvalues_data[c(2,3,9:18)]
PRS_eigenvalue_adjust<-PRS_eigenvalue_adjust[complete.cases(PRS_eigenvalue_adjust),]
#Reduce to relevant Samples
PRS_eigenvalue_adjust<-PRS_eigenvalue_adjust[PRS_eigenvalue_adjust$GAP_ID%in%PRS_to_adjust$gap_id,]
colnames(PRS_eigenvalue_adjust)[2]<-"gap_id"
#change to integer
PRS_eigenvalue_adjust$gap_id<-as.integer(as.character(PRS_eigenvalue_adjust$gap_id))
#combine PRS values and eigenvalues.
Combined_PRS<-inner_join(PRS_eigenvalue_adjust,PRS_to_adjust,by="gap_id")

#add rownames
rownames(Combined_PRS)<-Combined_PRS$sampleID

#Make adjustment data and pdata
names(Combined_PRS)
PRS_matrix<-t(Combined_PRS[,15:24])


#Model For matrix, have samples as columns. variables as rows.
#mod = model.matrix(~eigen_1+eigen_2,data=Combined_PRS)
formula_formod<-as.formula(paste("~", paste(paste("eigen_",1:10,sep=""), collapse="+")))
mod = model.matrix(formula_formod,data=Combined_PRS)
fit = lm.fit(mod,t(PRS_matrix))
#add residuals to average expression data
PRS_adjusted<-t(fit$residuals)+apply(PRS_matrix,1,mean)

#clean PRS_adjusted
PRS_adjusted<-t(PRS_adjusted)
colnames(PRS_adjusted)<-paste(colnames(PRS_adjusted),"_adj",sep="")
PRS_adjusted<-as.data.frame(PRS_adjusted)
PRS_adjusted$sampleID<-rownames(PRS_adjusted)

PRS_adjusted_full<-cbind(Combined_PRS,PRS_adjusted)

#combine with orginal data
reduced_pheno_final_adj_PRS<-full_join(reduced_pheno_final,PRS_adjusted,by="sampleID")


title = "PRS 0.1 adjusted for 10 principle components"
file_name = "PRS_01_adjusted_for_10_principle_components"
ggplot(data = reduced_pheno_final_adj_PRS, 
       aes(x = ICD_DSM, y=PRS_0.1_adj,colour = ICD_DSM)) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~ Ethnicity)+
  ggtitle(title)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(P0_figs_dir,title,".png",sep=""))


title = "PRS 0.1 before adjustment"
file_name = "PRS_01_before_adjustment"
ggplot(data = reduced_pheno_final_adj_PRS, 
       aes(x = ICD_DSM, y=PRS_0.1,colour = ICD_DSM)) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~ Ethnicity) +
  ggtitle(title)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(P0_figs_dir,file_name,".png",sep=""))



# Save with good probes and prs ------------------------------

#check pheno df
head(reduced_pheno_final_adj_PRS)
names(reduced_pheno_final_adj_PRS)

#add to lumibatch
pData(eset_bg_log2_rsn_SVA_Good)<-reduced_pheno_final_adj_PRS
names(pData(eset_bg_log2_rsn_SVA_Good))
save(eset_bg_log2_rsn_SVA_Good,file=paste(P0_output_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))



# Demographics Table ------------------------------------------------------


load(file=paste(P0_output_dir,"GAP_FEP_small_Gene_Expression_Data.RData",sep=""))

PRS_index<-c(31:40)

phenodat<-pData(eset_bg_log2_rsn_SVA_Good)
phenodat[PRS_index]<-scale(phenodat[PRS_index])

names(phenodat)

#######Table1: Demographics
##variables
lVars <- c("Sex","Age", "Ethnicity","BMI","Tobacco",names(phenodat)[PRS_index])
cVars <- c("Ethnicity","Tobacco")
lVars2 <- c("Medication","dsmiv","icd10","ICD_DSM","PanssScore","PanssPositive","PanssNegative","PanssPsycho")
cVars2 <- c("Medication","dsmiv","icd10","ICD_DSM")


##table statified by case control
stratVar <- "Phenotype" 
table_name <-"Table_1_Demographics.csv"
table_dem_fun(phenodat,table_name,stratVar,lVars,cVars)
table_name <-"Table_2_Clinicalinformation.csv"
table_dem_fun(phenodat,table_name,stratVar,lVars2,cVars2)

##table statified by Diagnosis consensus.
stratVar <- "ICD_DSM" 
table_name <-"Table_1_SczCat_Demographics.csv"
table_dem_fun(phenodat,table_name,stratVar,lVars,cVars)
table_name <-"Table_2_SczCat_Clinicalinformation.csv"
table_dem_fun(filter(phenodat,Phenotype=="FEP"),table_name,stratVar,lVars2,cVars2)

##table statified by Diagnosis consensus.
stratVar <- "ICD_DSM" 
table_name <-"Table_1_SczCat_white_Demographics.csv"
table_dem_fun(filter(phenodat,Ethnicity =="White"),table_name,stratVar,lVars,cVars)
table_name <-"Table_2_SczCat_white_Clinicalinformation.csv"
table_dem_fun(filter(phenodat,Phenotype=="FEP" & Ethnicity =="White"),table_name,stratVar,lVars2,cVars2)


# Regress out CellMix and Demographics ------------------------------------

# Make pheno_Data + expresison data
pData_rAESB<-pData(eset_bg_log2_rsn_SVA_Good)
exprs_rAESB<-exprs(eset_bg_log2_rsn_SVA_Good)
names(pData_rAESB)

# Regress out variables
mod = model.matrix(~Tc+neutro+Sex+Age+Ethnicity,data=pData_rAESB)
fit = lm.fit(mod,t(exprs_rAESB))
#add residuals to average expression data
exprs_res_cor_adj<-t(fit$residuals)+apply(exprs_rAESB,1,mean)
exprs_res_cor_adj[1:10,1:10] #new
exprs_rAESB[1:10,1:10] # old



# Combine GX data and Pheno data ------------------------------------------

# remove duplicate probes and use gene symbols
gene<-t(exprs_res_cor_adj)
probenames<-eset_bg_log2_rsn_SVA_Good@featureData@data

#change to gene symbol
all.equal(colnames(gene),probenames$nuID) #Must be true
colnames(gene) <- probenames$TargetID

gene_expression_matrix<-t(gene)#delete?
exprs_data<-gene
dim(exprs_data)

# pheno data
names(pData(eset_bg_log2_rsn_SVA_Good))
pheno_data<-pData(eset_bg_log2_rsn_SVA_Good)
names(pheno_data)


#Add phenotype to Gene expression. 
GX_DF_adj<-droplevels(merge(pheno_data,exprs_data, by.x="sampleID",by.y="row.names"))
row.names(GX_DF_adj)<-GX_DF_adj[,1]
GX_DF_adj[1:10,1:50]


#Save results
save(GX_DF_adj,file=paste(P0_output_dir,"GX_DF_adj_data.Rdata",sep=""), compress = T)

