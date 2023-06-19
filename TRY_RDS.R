# devtools::install_github(repo="krumsieklab/maplet@v1.0.1", subdir="maplet")
library(shiny)
library(shinyWidgets)
library(maplet)
library(shinydashboard)
library(shinyjs)
library(tidyverse)
library(DT)
library(plotly)
library(openxlsx)
library(readxl)
library(RColorBrewer)
library(grid)
library(graphics)
# refer help functions
#source("help_functions.R")

#file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")

D <- mt_load_se_xls(file='D_Olink_new.xlsx', sheet_names = c("data","proteins","samples")) %>%
  # # log assay dimensions and number of columns for both metabolite and clinical annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# Save the SE object as an RDS file
saveRDS(D, file = "se_object.rds")
