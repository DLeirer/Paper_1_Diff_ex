
##### NOT NEEDED!!!! ####

# Install and load Libraries: ---------------------------------------------
rm(list = ls())
dev.off()
library(lumi)
library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)
library(tableone)
library(stargazer)

# Set directories: --------------------------------------------------------

### DONT RUN!!! Incorporated into Cellmix script.
#setwd("/home/daniel/Documents/A_Year_4_PHD/Paper_1_Diff_ex")
#getwd()

#data_dir <-"./data/"
#P0_output_dir <-"./P0_Characterise/output/"
#P0_figs_dir <-"./P0_Characterise/figs/"



# Load gene expression data (Should be LumiBatch object): -----------------

lumidata<-"GAP_FEP_eset_linear_adj_Data.RData"
load(paste(P0_output_dir,lumidata,sep=""))



# Functions: --------------------------------------------------------------

## Table function 1 demographics
table_dem_fun<-function(pdata,table_name,stratify,listVars,catVars){
  table1 <- CreateTableOne(vars = listVars, data = pdata, factorVars = catVars,strata=c(stratify),includeNA = T)
  table1print<-print(table1)
  table1print<-table1print[,-length(names(data.frame(table1print)))]
  write.csv(table1print, file=paste(P0_output_dir,"no_latex_",table_name,sep=""),row.names = TRUE, col.names = TRUE,quote=FALSE,sep = ",")
  write.csv(stargazer(table1print,summary=FALSE), file=paste(P0_output_dir,"latex_",table_name,sep=""),row.names=F,col.names = F,sep="")
}






# Demographics: -----------------------------------------------------------



phenodat<-pData(eset_linear_adj)[1:24]
phenodat[64:73]<-scale(phenodat[64:73])

phenodat
names(phenodat)

#######Table1: Demographics
##variables
lVars <- c("Gender","Age", "Ethnicity","BMI","Tobacco",names(phenodat)[65:73])
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





