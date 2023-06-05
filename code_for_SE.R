# FILE:    lipidomics_analysis_pipeline_psoriasis.R
# DATE:    11 October 2022
# AUTHOR:  Nisha Stephan
# PURPOSE: reading metabolon lipidomics data and processing it in autonomics and/or maplet
# This is setup under github project https://github.com/wcmq-abmc/Psoriasis         

# cleanup
rm(list=ls())

# libraries
library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(cowplot)
library(plotly)
#library(xlsx)

source("mt_plots_volcano_log_pvalue.R")
source("mt_plots_pls.R")
source("mt_plots_boxplot.R")

source("mt_reporting_html_resize.R")
source("mt_plots_volcano_log_pvalue.R")
source("mt_plots_box_scatter_test.R")



old <- theme_set(theme_bw())
theme_set(theme_bw())

# define basename for output
basename = "lipidomics_analysis_psoriasis"
subgroupvar = 'GROUP_ID'

# open a log file
logfile = paste(basename, ".log", sep = "")
options(warn=-1); sink(); options(warn=0)
sink(file = logfile, split = TRUE)
cat("Start logging to", logfile, "\n")
print(Sys.time())

######################################################################################################
# # load the data from excel file  "WCQA-01-22MLCLP PLASMA CLP DATA TABLES.xlsx"
######################################################################################################

file_lipidomics_data <- "WCQA-01-22MLCLP-PLASMA-CLP-DATA-TABLES.xlsx"


pheno_data <- "pheno_data.xlsx"
# To get the checksum
# library(openssl)
# as.character(md5(file(file_lipidomics_data, open = "rb")))



D_CLP <-
  # validate checksum
  mt_load_checksum(file=file_lipidomics_data, checksum = "83c4beae42b164f609faf6218485c7ec") %>%
  # load data
  mt_load_metabolon_lipidomics(file=file_lipidomics_data, sheet_list = c("Lipid Class Concentrations", "Species Concentrations","Fatty Acid Concentrations")) %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# the make.names function adds X to the rownames



rowData(D_CLP)$CHEMICAL_NAME =  make.names(rowData(D_CLP)$name)
rownames(D_CLP) = rowData(D_CLP)$CHEMICAL_NAME

# Requirements for autonomics

rowData(D_CLP)$name = rownames(D_CLP)
rowData(D_CLP)$BIOCHEMICAL = rownames(D_CLP)
rowData(D_CLP)$feature_id = rownames(D_CLP)
colData(D_CLP)$sample_id = colData(D_CLP)$CLIENT_SAMPLE_ID 
colnames(D_CLP) = colData(D_CLP)$CLIENT_SAMPLE_ID 



# Read additional phenotypes

pheno = readxl::read_excel(path = pheno_data, sheet = "Sheet1", col_names = TRUE,
                           na = c("", "NA"))

pheno$Age = pheno$`Age        (in years)` 
pheno$`Age        (in years)` = NULL
pheno$CLIENT_SAMPLE_ID = pheno$`Sample Code/Client sample ID          (NP11-XXX)` 

pheno$`Sample Code/Client sample ID          (NP11-XXX)` = NULL


pheno$Ethnicity = pheno$`Ethnicity                            (Arab, Asian, Caucasian)`
pheno$`Ethnicity                            (Arab, Asian, Caucasian)`= NULL

pheno$Gender = pheno$"Gender              (male, female)"
pheno$"Gender              (male, female)" = NULL

pheno$BMI = pheno$"BMI             (kg*m-2)"
pheno$"BMI             (kg*m-2)" = NULL

pheno_add = pheno[,c("CLIENT_SAMPLE_ID","Gender", "BMI", "Age","Ethnicity", "HbA1C")]





data_join <- merge(as.data.frame(colData(D_CLP)), pheno_add,by.x = "CLIENT_SAMPLE_ID", by.y = "CLIENT_SAMPLE_ID" ,all = TRUE)


data_join_select <- data_join %>% select(CLIENT_SAMPLE_ID, Gender, GENDER, BMI.x,BMI.y,Age,Ethnicity,RACE_ETHNICITY,HbA1C) %>% 
                    mutate(gen = ifelse(Gender == GENDER, TRUE, FALSE)) %>% 
                    mutate(bm = ifelse(BMI.x == BMI.y, TRUE, FALSE)) %>%
                    mutate(eth = ifelse(Ethnicity == RACE_ETHNICITY, TRUE, FALSE))

data_join_select %>% filter(gen== FALSE| bm == FALSE | eth == FALSE |(is.na(Gender) & !is.na(GENDER)))
#colData(D_CLP)$AGE =  NULL
colData(D_CLP)$AGE =  data_join$Age
colData(D_CLP)$HbA1C = as.numeric(data_join$HbA1C)
#colData(D_CLP)$GENDER = NULL
colData(D_CLP)$GENDER = data_join$Gender
colData(D_CLP)$RACE_ETHNICITY = NULL
colData(D_CLP)$ETHNICITY = data_join$Ethnicity
#colData(D_CLP)$BMI = NULL
colData(D_CLP)$BMI = data_join$BMI.y
colData(D_CLP)$GROUP_ID = colData(D_CLP)$Group


#############################################################################

# PART 2 - DATA CLEANING ------------plot and decide on sample and feature missingness % to filter later

##############################################################################
source("mt_plots_volcano_log_pvalue.R")
D <- D_CLP


D %>% mt_reporting_heading(heading = "Data Clean-up", lvl = 1) %>%
  # filter samples
  mt_modify_filter_samples(filter = !is.na(GROUP_ID)) %>%
  # modify variable to factor
  mt_anno_apply(anno_type = "samples", col_name = "GROUP_ID", fun = as.factor) %>%
  mt_reporting_data() 

# PART 3.1 - PREPROCESSING: ANALYSIS MISSING VALUES ----------------------------------------------------





colData(D)$Diabetic_status <- case_when(
  (colData(D)$GROUP_ID == "Healthy" | colData(D)$GROUP_ID == "Psoriasis_Naïve" |colData(D)$GROUP_ID == "Psoriasis_FollowUp") ~ "Non-diabetic",
  (colData(D)$GROUP_ID == "Psoriasis_Prediab_Naïve" | colData(D)$GROUP_ID == "Psoriasis_Prediab_FollowUp")  ~ "Prediab",
  TRUE ~ "diabetic"
)

save(D, file="SE.Rdata")