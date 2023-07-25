#All 3 comments are working
#interactiveness for statistical analysis tab Working
#Differential expression separated, diff exp analaysis, correlation analysis working
#Partial Correlation working
# Pathway enrichment also working: tables for both above are working not plots
#Loading animations working

rm(list=ls())

################################################################################
####################load packages, functions and result SE object ##############
################################################################################

# load packages
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
library(DT)
library(shinymaterial)
library(shinycssloaders)


# refer help functions
source("help_functions.R")

# load SE with fixed name
#load("SE.Rdata")

# extract object names from result SE 'D'

options(shiny.maxRequestSize = 50*1024^2)  # Set maximum file size to 50MB (50 * 1024^2 bytes)


################################################################################
########################## Define UI for Shiny application #####################
################################################################################

ui <- fluidPage(
  
  # set appearance customization -------------------------------------------------
  
  theme = "bootstrap.css",
  includeCSS("www/style.css"),
  setBackgroundColor("#FFFFFF"),# set canvas background color
  div(style = "padding: 1px 0px; width: '100%'",
      titlePanel(
        title = "",
        windowTitle = "Maplet"
      )
  ),
  # remove shiny "red" warning messages on GUI
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  # adjust tab height
  tags$head(
    tags$style(HTML(' .navbar {
                          height: 80px;
                          min-height:25px !important;
                        }
                      .navbar-nav > li > a, .navbar-brand {
                            padding-top:1px !important; 
                            padding-bottom:1px !important;
                            height: 80px;
                            }'))),
  navbarPage(
    # embed Maplet logo and title
    title = div(img(src='logo.png',
                    style="float:left; margin-top: 5px; padding-right:20px;padding-bottom:5px",
                    height = 60),
                tags$a("Krumsiek Lab", href="https://github.com/krumsieklab/maplet-shiny", style="color: White"),
                tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://weill.cornell.edu\"><img src=\"wcm2.png\" alt=\"logo\" style=\"float:right;height:50px;margin-top: 10px; padding-right:1px; \"> </a></div>');console.log(header)")),
                windowTitle = "Maplet"),
    # sticky tabs while scrolling main panel
    position = c("fixed-top"), 
    
    # Define layout of Uploading Data(coded as mod0) ----------------------------------------------------
    
    tabPanel(HTML(paste("Upload", "Data", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod0_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            
                            tags$p(
                              HTML("<b>Hint:<br></b>Please upload an excel file which with 3 sheets containing assay data, 
                              rowData and colData each OR Please upload an rds file which contains an SE object to be used in the 
                              rest of the pipeline."
                                   
                              )),
                            br(),
                            
                            fileInput("file", "Upload Excel File", accept = c(".xlsx", ".xls")), 
                            # delay the output
                            actionButton("mod0_go", "Upload Excel File"), br(),
                            br(),
                            br(),
                            
                            fileInput("se_file", "Upload rds File with an SE object", accept = ".rds"), 
                            # delay the output for R file
                            actionButton("mod0_go_R", "Upload R file"), 
                            
                            br(),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'Upload' button after selection."
                              )),
                            br()
                            
               ), 
               mainPanel(id = "mod0_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         uiOutput("mod0_sheet_dropdowns"),
                         
                         br(), 
                         br(), 
                         br(),
                         #verbatimTextOutput("mod0_assay_display")
                         uiOutput("mod0_sheet_viewer"),
                         DTOutput("excel_table")
                         
               )
             )
    ),
    
    
    # Define layout of Module-Real-Time Pipeline(coded as mod6) ----------------------------------------------------
    tabPanel(HTML(paste("Real-Time", "Pipeline", sep = "<br/>")), 
             # Sidebar layout with input and output definitions ----
             dashboardPage(
               dashboardHeader(disable = TRUE),
               dashboardSidebar(disable = TRUE),
               dashboardBody(
                 sidebarLayout(
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     id = "mod6_panel1",
                     style = "overflow-y: scroll; max-height: 700px; width: 80%; position:relative; margin-left: -5px; margin-top: 45px; margin-bottom: 5px;",
                     tags$p(
                       HTML("<b>Real-Time Pipeline Module</b> starts with original data, creates a pipeline and download it to local."
                       )),
                     tags$p(
                       HTML("Pipeline is constrained to run <b>Data Loading->Preprocessing->Differential Analysis</b> and no section should be skipped."
                       )),
                     tags$p(
                       HTML("The result SE object is dependent on the <b>instant parameters</b>."
                       )),
                     
                     
                     tags$hr(),
                     
                     
                     
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Preprocessing", width = "220px",
                         tags$p(HTML("Max % missingness per feature (Filter):")),
                         numericInput("mod6_filter_feat_max", label = NULL,
                                      value = 50,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Max % missingness per feature (normalization):")),
                         numericInput("mod6_feat_max_norm", label = NULL,
                                      value = 50,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Max % missingness per sample:")),
                         numericInput("mod6_filter_sample_max", label = NULL,
                                      value = 100,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Sample coloring column:")),
                         uiOutput("mod6_pre_sample_color_column"),
                         
                         tags$p(HTML("Reference sample column:")),
                         uiOutput("mod6_norm"),
                         
                         tags$p(HTML("Reference column value:")),
                         uiOutput("mod6_norm_value"),
                         
                         tags$p(HTML("PCA/UMAP coloring column:")),
                         uiOutput("mod6_pre_pca_color_column"),
                         tags$p(HTML("Heatmap annotation column:")),
                         uiOutput("mod6_pre_heatmap_anno_column"),
                         tags$p(HTML("Heatmap annotation row:")),
                         uiOutput("mod6_pre_heatmap_anno_row"),
                         tags$p(HTML("Run to see log text of data loading and preprocessing. This step may cost a few seconds to run.")),
                         actionButton("mod6_go_preprocess", "Run", width = "110px"),
                         
                         tags$p(HTML("Reference sample column:")),
                         
                         
                         renderUI({selectInput("mod6_reference_samp",
                                               label = NULL,
                                               selected = "GROUP_ID",
                                               choices = names(colData(D())),
                                               width = "220px")
                         }),
                         
                         tags$p(HTML("Reference column value:")),
                         renderUI({selectInput("mod6_reference_val", label = NULL,
                                               choices = colData(D()) %>% as.data.frame() %>%
                                                 dplyr::select(input$mod2_reference_samp) %>% unique() %>%
                                                 unlist() %>% as.character(),width = "220px" )
                         })
                         
                         
                     ),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Differential Analysis", width = "220px",
                         tags$p(HTML("Outcome variable:")),
                         uiOutput("mod6_outcome"),
                         checkboxInput("mod6_outcome_binary", "Binary outcome?", FALSE),
                         tags$p(HTML("Type of analysis:")),
                         selectInput("mod6_analysis_type", label = NULL,
                                     width = "220px",
                                     choices = c("lm","pearson","spearman","kendall"),
                                     selected = "lm"),
                         tags$p(HTML("Multiple testing correction:")),
                         selectInput("mod6_mult_test_method", label = NULL,
                                     width = "220px",
                                     choices = c("BH","bonferroni","BY"),
                                     selected = "BH"),
                         tags$p(HTML("Significance threshold:")),
                         numericInput("mod6_sig_threshold", label = NULL,
                                      value = 0.05,
                                      min = 0,
                                      max = 1,
                                      step = 0.01,
                                      width = "220px"),
                         tags$p(HTML("Pathway aggregation in barplot:")),
                         uiOutput("mod6_group_col_barplot"),
                         tags$p(HTML("Barplot coloring column:")),
                         uiOutput("mod6_color_col_barplot"),
                         tags$p(HTML("Run to see log text of data loading, preprocessing and differential analysis. This step may cost a few seconds to run.")),
                         actionButton("mod6_go_differ", "Run", width = "110px")
                     )
                   ),
                   
                   # Main panel for displaying outputs ----
                   mainPanel(
                     id = "mod6_panel2", 
                     style = "overflow-y: auto; max-height: 85vh; position: absolute; left: 28%",
                     br(), 
                     br(), 
                     br(), 
                     # Output: Data file ----
                     #tags$p(HTML("Downloading SE.Rdata may cost more than one minute. Please wait for the prompt.")),
                     #downloadButton("download_se", "Download result SE .Rdata"),
                     br(),
                     br(),
                     #uiOutput("mod6_main_panel1"),
                     # br(),
                     #br(),
                     #fluidRow(column(width = 12, uiOutput("mod6_main_panel")))
                     uiOutput("mod6_main_panel")
                     
                   )
                 )
               )
             )
    ),
    
    # Define layout of Pre-processing (coded as mod2) ----------------------------------------------------
    tabPanel(HTML(paste("Pre-processing")), 
             # Sidebar layout with input and output definitions ----
             dashboardPage(
               dashboardHeader(disable = TRUE),
               dashboardSidebar(disable = TRUE),
               dashboardBody(
                 sidebarLayout(
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     id = "mod2_panel1",
                     style = "overflow-y: scroll; max-height: 700px; width: 80%; position:relative; margin-left: -5px; margin-top: 45px; margin-bottom: 5px;",
                     tags$p(
                       HTML("<b>Pre-processing Module</b> starts with original data, creates a pipeline and download it to local."
                       )),
                     tags$p(
                       HTML("Pipeline is constrained to run <b>Data Loading->Preprocessing->Differential Analysis</b> and no section should be skipped."
                       )),
                     tags$p(
                       HTML("The result SE object is dependent on the <b>instant parameters</b>."
                       )),
                     
                     
                     tags$hr(),
                     
                     
                     
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Missingness Analysis", width = "220px",
                         tags$p(HTML("Max % missingness per feature (Filter):")),
                         numericInput("mod2_filter_feat_max", label = NULL,
                                      value = 50,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Max % missingness per feature (normalization):")),
                         numericInput("mod2_feat_max_norm", label = NULL,
                                      value = 50,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Max % missingness per sample:")),
                         numericInput("mod2_filter_sample_max", label = NULL,
                                      value = 100,
                                      min = 0,
                                      max = 100,
                                      step = 5,
                                      width = "220px"),
                         tags$p(HTML("Sample coloring column:")),
                         uiOutput("mod2_pre_sample_color_column"),
                         
                         
                         tags$p(HTML("Run to see the plots. This step may cost a few seconds to run.")),
                         actionButton("mod2_go_missingness", "Run", width = "110px")
                     ),
                     
                     tags$hr(),
                     
                     box(
                       solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                       title = "Filtration", width = "220px",
                       tags$p(HTML("Select dimension to filter on:")),
                       radioButtons("mod2_choice_filter", "Select one dimension:", 
                                    choices = list("Column Data" = "mod2_col", 
                                                   "Row Data" = "mod2_row")
                       ),
                       br(),
                       uiOutput("mod2_dimension_filter"),
                       br(),
                       uiOutput("mod2_dimension_filter_multi"),
                       br(),
                       tags$p(
                         HTML("Run to filter the data. This step may cost a few seconds to run.")
                       ),
                       br(),
                       # delay the output
                       actionButton("mod2_go_filter", "Run")
                     ),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Imputation", width = "220px",
                         
                         
                         tags$p(HTML("Type of imputation:")),
                         selectInput("mod2_impute_type", label = NULL,
                                     width = "220px",
                                     choices = c("Min value","knn"),
                                     selected = "Min value"),
                         uiOutput("mod2_dimension_ui"),
                         tags$p(HTML("Run to see plots, Imputation. This step may cost a few seconds to run.")),
                         actionButton("mod2_go_impute", "Run", width = "110px")
                     ),
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Normalization", width = "220px",
                         
                         checkboxInput("mod2_trans_log", "Log transform the data?", TRUE),
                         tags$p(HTML("Type of normalization:")),
                         selectInput("mod2_norm_type", label = NULL,
                                     width = "220px",
                                     choices = c("probabilistic quotient","external column"),
                                     selected = "probabilistic quotient"),
                         uiOutput("mod2_dimension_ui_norm"),
                         
                         tags$p(HTML("Run to see plots, after normalization This step may cost a few seconds to run.")),
                         actionButton("mod2_go_norm", "Run", width = "110px")
                     )
                   ),
                   
                   # Main panel for displaying outputs ----
                   mainPanel(
                     id = "mod2_panel2", 
                     style = "overflow-y: auto; max-height: 85vh; position: absolute; left: 28%",
                     br(), 
                     br(), 
                     br(), 
                     # Output: Data file ----
                     #tags$p(HTML("Downloading SE.Rdata may cost more than one minute. Please wait for the prompt.")),
                     #downloadButton("download_se", "Download result SE .Rdata"),
                     br(),
                     br(),
                     #uiOutput("mod6_main_panel1"),
                     # br(),
                     #br(),
                     #fluidRow(column(width = 12, uiOutput("mod6_main_panel")))
                     uiOutput("mod2_main_panel"),
                     br(),
                     br(),
                     uiOutput("mod2_main_panel_part2"),
                     
                   )
                 )
               )
             )
    ),
    
    
    
    
    # Define layout of Module-Annotations Explorer(coded as mod5) ----------------------------------------------------
    tabPanel(HTML(paste("Annotations", "Explorer", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod5_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Annotations Explorer Module</b> creates tables, distribution plots, or other graphics to explore the SE object."
                              )),
                            radioButtons("mod5_dimension", "Select one dimension:", 
                                         choices = list("Column Data" = "col", 
                                                        "Row Data" = "row")
                            ),
                            br(),
                            uiOutput("mod5_dimension_ui"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                              )),
                            br(),
                            # delay the output
                            actionButton("mod5_go", "Update")
               ), 
               mainPanel(id = "mod5_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         uiOutput("mod5_output_ui")
               )
             )
    ),
    
    
    
    ###################################################
    
    
    # Define layout of Module-2D Projection(coded as mod3) ----------------------------------------------------
    
    tabPanel(HTML(paste("2D", "Projection", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod3_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>2D Projection Module</b> generates an interactive 2D projection of PCA/UMAP."
                              )),
                            tags$p(
                              HTML("It displays a drop-down menu of all colData columns for coloring."
                              )),
                            # select one plot type
                            radioButtons("mod3_select_plot", "Select one plot type:", 
                                         choices = list("PCA" = "pca", 
                                                        "UMAP" = "umap",
                                                        "PLS" = "pls")
                            ),
                            # function argument
                            uiOutput("mod3_pca_data"),
                            # select coloring colData and factor it
                            uiOutput("mod3_plot_argument"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                              )),
                            br(),
                            # delay the output
                            actionButton("mod3_go", "Update")
               ), 
               mainPanel(id = "mod3_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         # plotly
                         downloadButton("mod3_download_plotly", "download plotly"),
                         plotlyOutput('mod3_plot', height = 700)
               )
             )
    ),
    
    
    
    # Define layout of Module-Differential Expression Analysis(coded as mod7) ----------------------------------------------------
    
    
    tabPanel(HTML(paste("Statistical", "Analyses", sep = "<br/>")),
             # Sidebar layout with input and output definitions ----
             dashboardPage(
               dashboardHeader(disable = TRUE),
               dashboardSidebar(disable = TRUE),
               dashboardBody(
                 
                 sidebarLayout(
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     
                     id = "mod7_panel1",
                     style = "overflow-y: scroll; max-height: 700px; width: 80%; position:relative; margin-left: -5px; margin-top: 45px; margin-bottom: 5px;",
                     tags$p(
                       HTML("<b>Real-Time Pipeline Module</b> starts with original data, creates a pipeline and download it to local."
                       )),
                     tags$p(
                       HTML("Pipeline is constrained to run <b>Data Loading->Preprocessing->Differential Analysis</b> and no section should be skipped."
                       )),
                     tags$p(
                       HTML("The result SE object is dependent on the <b>instant parameters</b>."
                       )),
                     
                     
                     tags$hr(),
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Differential Analysis", width = "220px",
                         tags$p(HTML("Outcome variable:")),
                         uiOutput("mod7_outcome_diff"),
                         checkboxInput("mod7_outcome_binary_diff", "Binary outcome?", FALSE),
                         checkboxInput("mod7_filt_samp_diff", "Filter Samples?", FALSE),
                         
                         uiOutput("mod7_filt_samp_dim_diff"),
                         
                         tags$p(HTML("Type of analysis:")),
                         selectInput("mod7_analysis_type_diff", label = NULL,
                                     width = "220px",
                                     choices = c("lm"),
                                     selected = "lm"),
                         
                         uiOutput("mod7_dimension_ui_diff"),
                         
                         
                         tags$p(HTML("Multiple testing correction:")),
                         selectInput("mod7_mult_test_method_diff", label = NULL,
                                     width = "220px",
                                     choices = c("BH","bonferroni","BY"),
                                     selected = "BH"),
                         tags$p(HTML("Significance threshold:")),
                         numericInput("mod7_sig_threshold_diff", label = NULL,
                                      value = 0.05,
                                      min = 0,
                                      max = 1,
                                      width = "220px"),
                         
                         tags$p(HTML("Pathway aggregation in barplot:")),
                         uiOutput("mod7_group_col_barplot_diff"),
                         actionButton("mod7_go_diff", "Run", width = "110px")
                         
                     ),
                     tags$hr(),
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Correlation Analysis", width = "220px",
                         tags$p(HTML("Outcome variable:")),
                         uiOutput("mod7_outcome_corr"),
                         checkboxInput("mod7_outcome_binary_corr", "Binary outcome?", FALSE),
                         
                         
                         tags$p(HTML("Type of analysis:")),
                         selectInput("mod7_analysis_type_corr", label = NULL,
                                     width = "220px",
                                     choices = c("pearson","spearman","kendall"),
                                     selected = "pearson"),
                         
                         uiOutput("mod7_dimension_ui_corr"),
                         
                         
                         tags$p(HTML("Multiple testing correction:")),
                         selectInput("mod7_mult_test_method_corr", label = NULL,
                                     width = "220px",
                                     choices = c("BH","bonferroni","BY"),
                                     selected = "BH"),
                         tags$p(HTML("Significance threshold:")),
                         numericInput("mod7_sig_threshold_corr", label = NULL,
                                      value = 0.05,
                                      min = 0,
                                      max = 1,
                                      width = "220px"),
                         tags$p(HTML("Pathway aggregation in barplot:")),
                         uiOutput("mod7_group_col_barplot_corr"),
                         
                         actionButton("mod7_go_corr", "Run", width = "110px")
                     ),
                     tags$hr(),
                     # box(solidHeader = T, collapsible = T, collapsed = T,
                     #     title="Pathway Analysis", width = "220px",
                     #     tags$p(HTML("Outcome variable:")),
                     #     uiOutput("mod7_outcome_path"),
                     #     
                     #     tags$p(HTML("Significance threshold:")),
                     #     #change below to accept very small values
                     #     numericInput("mod7_sig_threshold_path", label = NULL,
                     #                  value = 0.05,
                     #                  min = 0,
                     #                  max = 1,
                     #                  step = 0.01, #remove step argument
                     #                  width = "220px"),
                     # 
                     #     tags$p(HTML("Pathway aggregation in barplot:")),
                     #     uiOutput("mod7_group_col_barplot_path"),
                     #     tags$p(HTML("Barplot coloring column:")),
                     #     uiOutput("mod7_color_col_barplot"),
                     #     tags$p(HTML("Run to see log text of data loading, preprocessing and differential analysis. This step may cost a few seconds to run.")),
                     #     actionButton("mod7_go_path", "Run", width = "110px")
                     # ),
                     tags$hr(),
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Partial Correlation Network", width = "220px",
                         tags$p(HTML("Stat Name")),
                         uiOutput("mod7_stat_name"),
                         actionButton("mod7_go_partial_corr", "Run", width = "110px")
                     ),
                     tags$hr(),
                     box(solidHeader = T, collapsible = T, collapsed = T,
                         title="Pathway Enrichment", width = "220px",
                         tags$p(HTML("Outcome variable:")),
                         uiOutput("mod7_outcome_path_enrich"),
                         tags$p(HTML("rowData column to use for pathway fetching. Hint: The selected column must contain HMDB metabolite
                      identifiers")),
                         uiOutput("mod7_in_col"),
                         actionButton("mod7_go_path_enrich", "Run", width = "110px")
                     ),
                     tags$hr(),
                     downloadLink("download_link", "Download Report")
                     
                   ),
                   
                   # Main panel for displaying outputs ----
                   mainPanel(
                     id = "mod7_panel2",
                     style = "overflow-y: auto; max-height: 85vh; position: absolute; left: 28%",
                     br(),
                     br(),
                     br(),
                     # Output: Data file ----
                     #tags$p(HTML("Downloading SE.Rdata may cost more than one minute. Please wait for the prompt.")),
                     #downloadButton("download_se", "Download result SE .Rdata"),
                     br(),
                     br(),
                     #fluidRow(column(width = 12, uiOutput("mod7_main_panel")))
                     uiOutput("mod7_main_panel_1"),
                     
                     withSpinner(uiOutput("mod7_main_panel_graph1")),
                     br(),
                     withSpinner(uiOutput("mod7_main_panel_graph2")),
                     
                     br(),
                     withSpinner(uiOutput("mod7_main_panel_graph3")),
                    
                     DTOutput("mod7_main_panel_2")
                   )
                 )
                 
               )
             )
    ),
    
    # # Define layout of Module-All Results Explorer(coded as mod1) ----------------------------------------------------
    
    tabPanel(HTML(paste("All Results", "Explorer", sep = "<br/>")),
             sidebarLayout(
               sidebarPanel(id = "mod1_panel1",
                            # sidebar auto-scrolling with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>All Results Explorer Module</b> extracts all the result objects one at a time."
                              )),
                            tags$p(
                              HTML("Users can assess results in a drop-down menu that offers a list of a stat_name and a plot type (e.g. “missingness”, “pval”)."
                              )),
                            br(),
                            # select plot type or stats table
                            radioButtons("mod1_radio", "Select output type:",
                                         choices = list("Plot" = "plots",
                                                        "Table" = "stats"),
                                         selected = "stats"
                            ),
                            br(),
                            # define one UI object to select stat_name
                            uiOutput("mod1_select_statname_ui"),
                            br(),
                            # define one UI object to select output type
                            uiOutput("mod1_select_object_ui"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection. Some plots such as box plot or multiple plots may take dozens of seconds to show up."
                              )),
                            # delay the output
                            actionButton("mod1_go", "Update")
               ),
               mainPanel(id = "mod1_panel2",
                         # scrollable panel
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         br(),
                         br(),
                         br(),
                         # dynamic number of plots
                         uiOutput('mod1_output')
               )
             )
    )
    
    
    
  )
)

################################################################################
################ Define server logic required to draw outputs ##################
################################################################################

server <- function(input, output,session) {
  
  # Declare global variables
  
  pl_list <- reactiveValues()
  pdiffer_list <- reactiveValues()
  pmissing_list <- reactiveValues()
  pimpute_list <- reactiveValues()
  pnorm_list <- reactiveValues()
  pheatmap_list <- reactiveValues()
  
  # Define rendering logic of uploading data in Module-Upload Data (coded as mod0)
  
  categories <- reactiveValues(sheet1 = "assay", sheet2 = "assay", sheet3 = "assay")
  sheet_list <- reactiveVal()
  D_excel <- reactiveVal(NULL)
  se_r <- reactiveVal(NULL)                     
  D <- reactiveVal(NULL)
  
  
  observeEvent(input$file, { 
    print("happening")
    # Read the uploaded Excel file
    
    file <- input$file$datapath
    sheets <- excel_sheets(file)
    updateSelectInput(session, "sheet_viewer", choices = sheets)
    
    # Create a named list to store sheet names and their corresponding category
    sheet_list(lapply(sheets, function(sheet) {
      name <- paste0("category_", sheet)
      dropdown_id <- paste0("dropdown_", sheet)
      input_id <- paste0("input_", sheet)
      output_id <- paste0("output_", sheet)
      
      list(
        sheet_name = sheet,
        #dropdown_id = dropdown_id,
        input_id = input_id,
        output_id = output_id,
        category = reactiveVal("data")  # Initialize with a default value
      )
    }))
    
    # Reset the selected categories when a new file is uploaded
    categories$sheet1 <- "data"
    categories$sheet2 <- "data"
    categories$sheet3 <- "data"
    
    # If the file has less than 3 sheets, show an error message and hide the dropdowns
    if (length(sheets) < 3) {
      showModal(
        modalDialog(
          title = "Error",
          "The file must have at least three sheets.",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    
  })
  
  
  
  output$mod0_sheet_dropdowns <- renderUI({
    dropdowns <- c("assay", "row", "col")
    file <- input$file$datapath
    if (!is.null(file)) {
      sheets <- excel_sheets(file)
      
      # Only render the dropdowns if the file has at least 3 sheets
      if (length(sheets) >= 3) {  
        # Generate dropdown menus for each category
        ui_list <- lapply(dropdowns, function(dropdown) {
          selectInput(
            inputId = paste0("dropdown_", dropdown),
            label = paste("Choose", dropdown),
            choices = sheets
          )
        })
        
        do.call(tagList, ui_list)
      }
      
      
      
      
    }
  })
  
  output$mod0_sheet_viewer <- renderUI({
    file <- input$file$datapath
    if (!is.null(file)) {
      sheets <- excel_sheets(file)
      selectInput("sheet_viewer", "Select Sheet", choices = NULL)
    }
  })
  
  
  observeEvent(input$mod0_go, {
    print("excel buton pressed")
    file <- input$file$datapath
    if (!is.null(file)) {
      sheet_list <- sheet_list()
      
      selected_data_sheet <- input$dropdown_assay
      
      
      selected_row_data_sheet <- input$dropdown_row
      
      selected_col_data_sheet <- input$dropdown_col
      
      print(selected_row_data_sheet)
      
      # Check if any duplicate selections exist
      if (length(unique(c(selected_data_sheet, selected_row_data_sheet, selected_col_data_sheet))) < 3) {
        showModal(
          modalDialog(
            title = "Error",
            "Please select only one unique sheet per category.",
            easyClose = TRUE,
            footer = NULL
          )
        )
        return()  # Stop execution if there are duplicate selections
      }
      
      print(sheet_list)
      
      tryCatch(
        {
          D_excel(mt_load_se_xls(file=input$file$datapath, sheet_names=c(selected_data_sheet, selected_row_data_sheet, selected_col_data_sheet)))
          # Print the coldata
          print(names(colData(D_excel())))
          #output$mod0_assay_display <- renderPrint(assay(D_excel()))
          file <- input$file$datapath
          #sheets = excel_sheets(file)
          
          
          output$excel_table = DT::renderDataTable({
            sheet_v = input$sheet_viewer
            df <- read.xlsx(file, sheet_v)  # Adjust the sheet number as needed
            DT::datatable(df)
          })
          
          showModal(
            modalDialog(
              title = "Success",
              "SE object was created successfully. The assay data is printed for your reference. Please proceed to the next steps of the pipeline in the Navigation Bar.",
              easyClose = TRUE,
              footer = NULL
            )
          )
          
        },
        error = function(e) {
          showModal(
            modalDialog(
              title = "Error",
              "SE object was not created successfully. Please select the sheets from the dropdown again.",
              easyClose = TRUE,
              footer = NULL
            )
          )
        }
      )
      
    }
    
  })
  
  
  # Function to generate SE object from the uploaded R file
  #generateSEObject <- function(file_path) {
  #  file_content <- readLines(file_path)
  #  print("whyyyy")
  #  se <- readRDS(file_path)
  
  #  return(se)
  #}
  tryCatch({
    # Reactive expression for storing the SE object
    se_r <- eventReactive(input$mod0_go_R, {
      #req(input$se_file)
      file <- input$se_file
      file_path <- file$datapath
      se <- readRDS(file_path)
      #se2 <- readRDS("se_object.rds")
      #output$mod0_assay_display <- renderPrint(assay(se2))
      
      showModal(
        modalDialog(
          title = "Success",
          "SE object was created successfully. The assay data is printed for your reference. Please proceed to the next steps of the pipeline in the Navigation Bar.",
          easyClose = TRUE,
          footer = NULL
        )
      )
      return(se)
    })
  },
  error = function(e) {
    showModal(
      modalDialog(
        title = "Error",
        "SE object was not created successfully. Please select the sheets from the dropdown again.",
        easyClose = TRUE,
        footer = NULL
      )
    )
  }
  )
  
  # Print the dimensions of the SE object
  #observe({
  #  print(dim(se_r()))
  #})
  
  # Final variable to store the SE object
  D <- reactive({
    if (!is.null(D_excel())) {
      print("excel uploaded")
      D_excel()  # Use SE object from Excel file if it exists
    } else if (!is.null(se_r())) {
      print("R uploaded")
      se_r()  # Use SE object from R file if it exists
    } else {
      NULL  # Default value if no SE object is available
    }
  })
  
  # Print the dimensions of the final SE object
  observe({
    if (!is.null(D())) {
      print(dim(D()))
    }
  })
  
  
  
  # Define rendering logic of control widgets in Module-Real-Time Pipeline (coded as mod6)----------------------
  # control widget of selecting file
  
  # control widget of preprocessing
  
  
  # Define rendering logic of outputs in Module-Real-Time Pipeline(coded as mod6) ------------------------------
  
  
  
  output$mod6_pre_sample_color_column <- renderUI({
    selectInput("pre_sample_color_column", label = NULL,
                width = "220px",
                choices = names(colData(D()))
    )
  })
  
  #Newly added normalization column
  output$mod6_norm <- renderUI({selectInput("mod6_reference_samp",
                                            label = NULL,
                                            selected = "GROUP_ID",
                                            choices = names(colData(D())),
                                            width = "220px"
  )
  })
  
  
  
  #Newly added normalization value
  output$mod6_norm_value <- renderUI({selectInput("mod6_reference_val", 
                                                  label = NULL,
                                                  choices = colData(D()) %>% as.data.frame() %>%
                                                    dplyr::select(input$mod6_reference_samp) %>% unique() %>%
                                                    unlist() %>% as.character(),width = "220px" )
  })
  
  
  
  output$mod6_pre_pca_color_column <- renderUI({
    selectInput("pre_pca_color_column", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = names(colData(D()))
    )
  })
  
  output$mod6_pre_heatmap_anno_column <- renderUI({
    selectInput("pre_heatmap_anno_column", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = names(colData(D()))
    )
  })
  
  output$mod6_pre_heatmap_anno_row <- renderUI({
    selectInput("pre_heatmap_anno_row", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  # control widget of differential analysis
  output$mod6_outcome <- renderUI({
    selectInput("outcome_mod6", label = NULL,
                width = "220px",
                choices = names(colData(D()))
    )
  })
  
  output$mod6_group_col_barplot <- renderUI({
    selectInput("group_col_barplot_mod6", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod6_color_col_barplot <- renderUI({
    selectInput("color_col_barplot_mod6", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  # Define rendering logic of outputs in Module-Real-Time Pipeline(coded as mod6) ------------------------------
  
  
  
  
  observeEvent(input$mod6_go_preprocess,{
    # define main panel for preprocessing section
    # output$mod6_main_panel1 <- renderUI({
    #   tagAppendAttributes(verbatimTextOutput("log_preprocess"), 
    #                        style="white-space:pre-wrap;")
    # })
    pl_list$value <- get_plots_SE(D_preprocess())
    pl_list$length <-  pl_list$value %>% length()
    
    
    output$mod6_main_panel  <- renderUI({
      
      mod6_output_plotlist <-   lapply(1: pl_list$length, function(i){
        local({
          len_j <- length(pl_list$value[[i]])
          lapply(1:(len_j), function(j) {
            
            plotname <- paste("Plot", i,j, sep="")
            #print(paste0("The value of my variable from UI is ", plotname))
            
            plotOutput(plotname)
            
          })
        })
        
        
      })
      
      do.call(tagList, mod6_output_plotlist)
      
      
    })
    
    lapply(1: pl_list$length, function(i){
      local({
        
        len_j <- length(pl_list$value[[i]])
        
        lapply(1:(len_j), function(j) {
          
          plotname <- paste("Plot", i,j, sep="")
          
          output[[paste("Plot", i,j, sep="")]] <-
            renderPlot({
              #grid.force()
              pl_list$value[[i]][j]
              
            })
        })
      })
      
    })
    
  }) 
  
  observeEvent(input$mod6_go_differ,{
    # define main panel of differential analysis
    # output$mod6_main_panel <- renderUI({
    #   tagAppendAttributes(verbatimTextOutput("log_differ"), 
    #                       style="white-space:pre-wrap;")
    # })
    
    pdiffer_list$value <- get_plots_SE(D_differ())
    pdiffer_list$length <-  pdiffer_list$value %>% length()
    
    
    output$mod6_main_panel  <- renderUI({
      
      
      
      mod6_output_plotlist <-   lapply(1: pdiffer_list$length, function(i){
        local({
          len_j <- length(pdiffer_list$value[[i]])
          lapply(1:(len_j), function(j) {
            
            plotname <- paste("Plot", i,j, sep="")
            #print(paste0("The value of my variable from UI is ", plotname))
            
            plotOutput(plotname)
            
          })
        })
        
        
      })
      
      do.call(tagList, mod6_output_plotlist)
      
    })
    
    lapply(1: pdiffer_list$length, function(i){
      local({
        
        len_j <- length(pdiffer_list$value[[i]])
        
        lapply(1:(len_j), function(j) {
          
          plotname <- paste("Plot", i,j, sep="")
          
          output[[paste("Plot", i,j, sep="")]] <-
            renderPlot({
              #grid.force()
              pdiffer_list$value[[i]][j]
              
            })
        })
      })
      
    })
    
  })
  
  # get proprocessing SE
  D_preprocess <- reactive({
    ## preprocessing D
    D <- D() %>%
      mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
      mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
      mt_plots_missingness(feat_max=(input$mod6_filter_feat_max)/100,samp_max = (input$mod6_filter_sample_max)/100) %>%
      mt_pre_filter_missingness(feat_max = (input$mod6_filter_feat_max)/100, samp_max = (input$mod6_filter_sample_max)/100) %>%
      mt_plots_missingness(feat_max=(input$mod6_filter_feat_max)/100, samp_max = (input$mod6_filter_sample_max)/100) %>%
      mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
      mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
      mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "Original", plot_logged = T) %>%
      {.}
    if(!is.null(input$pre_batch_column)){
      D %<>%
        mt_pre_batch_median(batch_col = input$pre_batch_column)
    }
    D <- D %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After batch correction", plot_logged = T) %>%
      #mt_pre_trans_exp() %>%
      mt_pre_norm_quot(feat_max = (input$mod6_feat_max_norm)/100, ref_samples = (!!sym(input$mod6_reference_samp) == input$mod6_reference_val)) %>%
      #mt_pre_norm_quot(feat_max = (input$mod6_feat_max_norm)/100, ref_samples = GROUP_ID=="Healthy") %>%
      mt_plots_dilution_factor(in_col=input$pre_sample_color_column) %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After normalization", plot_logged = T) %>%
      
      mt_pre_trans_log() %>%
      mt_pre_impute_min() %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After imputation", plot_logged = T) %>%
      mt_pre_outlier_detection_univariate() %>%
      mt_reporting_data() %>%
      mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
      {.}
    ## add PCA/UMAP plots
    lapply(input$pre_pca_color_column, function(x){
      D <<- D %>%
        mt_plots_pca(scale_data = T, title = sprintf("scaled PCA - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
        mt_plots_umap(scale_data = T, title = sprintf("scaled UMAP - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
        {.}
    }) %>% invisible
    ## add heatmap
    D %<>%
      mt_plots_heatmap(scale_data = T, annotation_col = input$pre_heatmap_anno_column, annotation_row = input$pre_heatmap_anno_row,
                       clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
      {.}
    ## return D
    D
  })
  
  ## get differential analysis SE
  D_differ <- reactive({
    # Differential analysis D
    D <- D_preprocess()  %>%
      mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
      diff_analysis_func(var=input$outcome_mod6,
                         binary=input$mod6_outcome_binary,
                         analysis_type=input$mod6_analysis_type,
                         mult_test_method=input$mod6_mult_test_method,
                         alpha=input$mod6_sig_threshold,
                         group_col_barplot=input$group_col_barplot_mod6,
                         color_col_barplot=input$color_col_barplot_mod6) %>%
      {.}
    ## return D
    D
  })
  
  # render logic of the log text of preprocessing
  output$log_preprocess <- renderPrint({
    # loading log
    
    # preprocessing log
    text_preprocess <- get_log_text(D_preprocess())
    # paste log text
    
    cat(text_preprocess)
  })
  
  
  session_store <- reactiveValues()
  
  
  # render logic of the log text of differential analysis
  output$log_differ <- renderPrint({
    # loading log
    
    # preprocessing log
    text_preprocess <- get_log_text(D_preprocess())
    # differential analysis log
    text_differ <- get_log_text(D_differ())
    # paste log text
    str <- paste( text_preprocess, text_differ, sep = "\n")
    cat(str)
  })
  
  
  # Define rendering logic of outputs in Module-Pre-processing (coded as mod2) ------------------------------
  
  
  
  output$mod2_pre_sample_color_column <- renderUI({
    selectInput("pre_sample_color_column_mod2", label = NULL,
                width = "220px",
                selected = "GROUP_ID",
                choices = names(colData(D()))
    )
  })
  
  
  
  
  
  output$mod2_group_col_barplot <- renderUI({
    selectInput("group_col_barplot_mod2", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod2_color_col_barplot <- renderUI({
    selectInput("color_col_barplot_mod2", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod2_dimension_filter <- renderUI({
    switch(input$mod2_choice_filter,
           "mod2_col" = {
             selectInput("mod2_samp1_select",
                         "Select the sample to filter on",
                         choices = names(colData(D())),
                         width = "220px")
           },
           "mod2_row" = {
             selectInput("mod2_feat1_select",
                         "Select the feature to filter on",
                         choices = names(rowData(D())),
                         width = "220px")
           }
    )
  })
  
  output$mod2_dimension_filter_multi <- renderUI({
    switch(input$mod2_choice_filter,
           "mod2_col" = {
             multiInput("mod2_filter_col",
                        label = input$mod2_samp1_select,
                        choices = colData(D()) %>%
                          as.data.frame() %>%
                          dplyr::select(input$mod2_samp1_select) %>%
                          unique() %>%
                          unlist() %>%
                          as.character(),
                        width = "350px"
             )
           },
           "mod2_row" = {
             multiInput("mod2_filter_row",
                        label = input$mod2_feat1_select,
                        choices = rowData(D()) %>%
                          as.data.frame() %>%
                          dplyr::select(input$mod2_feat1_select) %>%
                          unique() %>%
                          unlist() %>%
                          as.character(),
                        width = "350px"
             )
           }
    )
  })
  
  
  
  output$mod2_dimension_ui <- renderUI({
    switch(input$mod2_impute_type,
           "knn"=numericInput("mod2_k_knn", 
                              "Number of nearest neighbors for KNN:", 
                              value = 10,
                              width = "220px"
           )
           
           
    )
    
  })
  
  
  output$mod2_dimension_ui_norm <- renderUI({
    switch(input$mod2_norm_type,
           "probabilistic quotient"= list(tags$p(HTML("Max % missingness per feature (normalization):")),
                                          numericInput("mod2_norm_feat_max", label = NULL,
                                                       value = 50,
                                                       min = 0,
                                                       max = 100,
                                                       step = 5,
                                                       width = "220px"),
                                          
                                          tags$p(HTML("Reference sample column:")),
                                          
                                          
                                          renderUI({selectInput("mod2_reference_samp",
                                                                label = NULL,
                                                                selected = "GROUP_ID",
                                                                choices = names(colData(D())),
                                                                width = "220px"
                                          )
                                          }),
                                          
                                          tags$p(HTML("Reference column value:")),
                                          renderUI({selectInput("mod2_reference_val", label = NULL,
                                                                choices = colData(D()) %>% as.data.frame() %>%
                                                                  dplyr::select(input$mod2_reference_samp) %>% unique() %>%
                                                                  unlist() %>% as.character(),width = "220px" )
                                          })
                                          
           ),
           "external column" = list(tags$p(HTML("Numeric-value column in colData:")),
                                    selectInput("mod2_ext_norm",
                                                label = NULL,
                                                choices = colData(D()) %>% as.data.frame() %>% dplyr::select(where(is.numeric)) %>% names(),
                                                selected = "GROUP_ID",
                                                width = "220px"
                                    )
                                    
                                    
           )
    )
    
  })
  
  
  
  # Define rendering logic of outputs in Module-Pre-processing(coded as mod2) ------------------------------
  
  
  # Define reactive values to store SE objects
  D_missingness <- reactiveVal(NULL)
  D_filter <- reactiveVal(NULL)
  D_impute <- reactiveVal(NULL)
  D_norm <- reactiveVal(NULL)
  
  D_for_analysis <- reactiveVal()
  
  
  
  
  #original version with all graphs - working
  # observeEvent(input$mod2_go_missingness,{
  #   
  #   pmissing_list$value <- get_plots_SE(D_missingness())
  #   pmissing_list$length <-  pmissing_list$value %>% length()
  #   
  #   # Update D_for_analysis
  #   D_for_analysis(D_missingness())
  #   
  #   
  #   output$mod2_main_panel  <- renderUI({
  #     
  #     
  #     
  #     mod2_output_plotlist <-   lapply(1: pmissing_list$length, function(i){
  #       local({
  #         len_j <- length(pmissing_list$value[[i]])
  #         lapply(1:(len_j), function(j) {
  #           
  #           plotname <- paste("Plot", i,j, sep="")
  #           #print(paste0("The value of my variable from UI is ", plotname))
  #           
  #           plotOutput(plotname)
  #           
  #         })
  #       })
  #       
  #       
  #     })
  #     
  #     do.call(tagList, mod2_output_plotlist)
  #     
  #     
  #   })
  #   
  #   lapply(1: pmissing_list$length, function(i){
  #     local({
  #       
  #       len_j <- length(pmissing_list$value[[i]])
  #       
  #       lapply(1:(len_j), function(j) {
  #         
  #         plotname <- paste("Plot", i,j, sep="")
  #         
  #         output[[paste("Plot", i,j, sep="")]] <-
  #           renderPlot({
  #             #grid.force()
  #             pmissing_list$value[[i]][j]
  #             
  #           })
  #       })
  #     })
  #     
  #   })
  #   
  # })
  
  #progress bar version
  observeEvent(input$mod2_go_missingness, {
    
    tryCatch({
    
      pmissing_list$value <- get_plots_SE(D_missingness())
      pmissing_list$length <-  pmissing_list$value %>% length()
      
      # Update D_for_analysis
      D_for_analysis(D_missingness())
      
      # Calculate the total number of plots
      total_plots <- pmissing_list$length
      
      # Initialize a variable to track the progress
      progress <- 0
      
      # Function to increment progress
      updateProgress <- function() {
        progress <<- progress + 1
      }
      
      # Render the plots
      output$mod2_main_panel <- renderUI({
        mod2_output_plotlist <- lapply(1:total_plots, function(i) {
          local({
            len_j <- length(pmissing_list$value[[i]])
            lapply(1:len_j, function(j) {
              plotname <- paste("Plot", i, j, sep = "")
              plotOutput(plotname) %>% withSpinner()
            })
          })
        })
        
        do.call(tagList, mod2_output_plotlist)
      })
      
      # Render individual plots
      lapply(1:total_plots, function(i) {
        local({
          len_j <- length(pmissing_list$value[[i]])
          lapply(1:len_j, function(j) {
            plotname <- paste("Plot", i, j, sep = "")
            output[[plotname]] <- renderPlot({
              updateProgress()
              pmissing_list$value[[i]][j]
            })
          })
        })
      })
      
    }, error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        glue::glue("An error occurred: {e$message}"),
        type = "error"
      )
    })
    

  })
  
  
  
  
  
  
  
  
  
  
  
  
  # pheatmap_list <- reactiveValues()
  
  observeEvent(input$mod2_go_impute,{
    # define main panel of differential analysis
    # output$mod6_main_panel <- renderUI({
    #   tagAppendAttributes(verbatimTextOutput("log_differ"), 
    #                       style="white-space:pre-wrap;")
    # })
    
    tryCatch({
    
      pimpute_list$value <- get_plots_SE_preprocess(D_impute(), title = c("Original","After imputation"))
      pimpute_list$length <-  pimpute_list$value %>% length()
      
      # Update D_for_analysis
      D_for_analysis(D_impute())
      
      # Initialize a variable to track the progress
      progress <- 0
      
      # Function to increment progress
      updateProgress <- function() {
        progress <<- progress + 1
      }
      
      output$mod2_main_panel  <- renderUI({
        mod2_output_plotlist <-   lapply(1: pimpute_list$length, function(i){
          local({
            len_j <- length(pimpute_list$value[[i]])
            lapply(1:(len_j), function(j) {
              
              plotname <- paste("Plot", i,j, sep="")
              #print(paste0("The value of my variable from UI is ", plotname))
              
              plotOutput(plotname) %>% withSpinner()
              
            })
          })
          
          
        })
        
        do.call(tagList, mod2_output_plotlist)
        
      })
      
      
      lapply(1: pimpute_list$length, function(i){
        local({
          
          len_j <- length(pimpute_list$value[[i]])
          
          lapply(1:(len_j), function(j) {
            
            plotname <- paste("Plot", i,j, sep="")
            
            output[[paste("Plot", i,j, sep="")]] <-
              renderPlot({
                #grid.force()
                updateProgress()
                pimpute_list$value[[i]][j]
                
              })
          })
        })
        
      })
      
    },error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        glue::glue("An error occurred: {e$message}"),
        type = "error"
      )
    })
    
  })
  
  
  
  
  observeEvent(input$mod2_go_norm,{
    # define main panel of differential analysis
    # output$mod6_main_panel <- renderUI({
    #   tagAppendAttributes(verbatimTextOutput("log_differ"), 
    #                       style="white-space:pre-wrap;")
    # })
    
    tryCatch({
    
      pnorm_list$value <- get_plots_SE_preprocess(D_norm(), title = c("Original","After normalization"))
      pnorm_list$length <-  pnorm_list$value %>% length()
      
      # Update D_for_analysis
      D_for_analysis(D_norm())
      
      
      output$mod2_main_panel  <- renderUI({
        
        
        
        mod2_output_plotlist <-   lapply(1: pnorm_list$length, function(i){
          local({
            len_j <- length(pnorm_list$value[[i]])
            lapply(1:(len_j), function(j) {
              
              plotname <- paste("Plot", i,j, sep="")
              #print(paste0("The value of my variable from UI is ", plotname))
              
              plotOutput(plotname)
              
            })
          })
          
          
        })
        
        do.call(tagList, mod2_output_plotlist)
        
      })
      
      lapply(1: pnorm_list$length, function(i){
        local({
          
          len_j <- length(pnorm_list$value[[i]])
          
          lapply(1:(len_j), function(j) {
            
            plotname <- paste("Plot", i,j, sep="")
            
            output[[paste("Plot", i,j, sep="")]] <-
              renderPlot({
                #grid.force()
                pnorm_list$value[[i]][j]
                
              })
          })
        })
        
      })
    },error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        glue::glue("An error occurred: {e$message}"),
        type = "error"
      )
    })
    
  })
  
  # get proprocessing SE
  D_missingness <- reactive({
    # Initialize an empty result
    D <- NULL
    
    # Use tryCatch to handle errors
    tryCatch({
      # Your original code block
      D <- D() %>%
        mt_reporting_heading(heading = "Preprocessing", lvl = 1) %>%
        mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
        mt_plots_missingness(feat_max = (input$mod2_filter_feat_max) / 100, samp_max = (input$mod2_filter_sample_max) / 100) %>%
        mt_pre_filter_missingness(feat_max = (input$mod2_filter_feat_max) / 100, samp_max = (input$mod2_filter_sample_max) / 100) %>%
        mt_plots_missingness(feat_max = (input$mod2_filter_feat_max) / 100, samp_max = (input$mod2_filter_sample_max) / 100) %>%
        mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
        mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
        mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
        mt_plots_sample_boxplot(color = !!sym(input$pre_sample_color_column_mod2), title = "Original", plot_logged = T)
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging
      print(paste("Error:", e))
      # Return NULL to indicate no valid result
      NULL
    })
    
    # Return the resulting D (could be NULL if an error occurred)
    D
  })
  
  
  
  
  
  
  D_impute <- reactive({
    # Initialize an empty result
    D <- NULL
    
    # Use tryCatch to handle errors
    tryCatch({
      # Your original code block
      D <- D_for_analysis()
      
      if (input$mod2_impute_type == "Min value") {
        D <- D %>% mt_pre_impute_min() 
      } else {
        D <- D %>% mt_pre_impute_knn(k = input$mod2_k_knn) 
      }
      
      D <- D %>% mt_plots_sample_boxplot(color = !!sym(input$pre_sample_color_column_mod2),
                                         title = "After imputation", plot_logged = T) %>%
        mt_pre_outlier_detection_univariate() %>%
        mt_reporting_data()
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging
      print(paste("Error:", e))
      # Return NULL to indicate no valid result
      NULL
    })
    
    # Return the resulting D (could be NULL if an error occurred)
    D
  })
  
  
  
  D_norm <- reactive({
    # Initialize an empty result
    D <- NULL
    
    # Use tryCatch to handle errors
    tryCatch({
      # Your original code block
      D <- D_for_analysis() 
      
      if (input$mod2_norm_type == "probabilistic quotient") {
        D <- D %>%
          #mt_pre_trans_exp() %>% commenting for simulated_data
          #mt_reporting_heading(heading = "STEP1") %>%
          mt_pre_norm_quot(feat_max = (input$mod2_norm_feat_max) / 100,
                           ref_samples = (!!sym(input$mod2_reference_samp) == input$mod2_reference_val)) %>%
          #mt_reporting_heading(heading = "STEP2") %>%
          mt_plots_dilution_factor(in_col = input$pre_sample_color_column_mod2) 
      } else {
        D <- D %>%
          #mt_pre_trans_exp() %>% simulated data
          mt_pre_norm_external(col_name = input$mod2_ext_norm)
      }
      
      if (input$mod2_trans_log) {
        D <- D %>%
          # mt_pre_trans_exp() %>% simulated_data
          mt_pre_trans_log()
      }
      
      D <- D %>%     
        mt_plots_sample_boxplot(color = !!sym(input$pre_sample_color_column_mod2),
                                title = "After normalization", plot_logged = T)
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging
      print(paste("Error:", e))
      # Return NULL to indicate no valid result
      NULL
    })
    
    # Return the resulting D (could be NULL if an error occurred)
    D
  })
  
  
  
  D_heatmap <- reactive({ 
    D <- D_norm() %>%
      
      mt_plots_heatmap(scale_data = T, annotation_col = input$pre_heatmap_anno_column, annotation_row = input$pre_heatmap_anno_row,
                       clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
      
      
      {.}
    
    
    D
  }) 
  
  #Filtration section
  
  observeEvent(input$mod2_go_filter, {
    if (input$mod2_choice_filter == "mod2_row") {
      selected_values <- input$mod2_filter_row
      
      # Construct logical expression based on selected values
      filter_expr <- paste0("rowData(D())$", input$mod2_feat1_select, " %in% c(", paste0("'", selected_values, "'", collapse = ","), ")")
      
      # Evaluate the filter expression
      D_filter(mt_modify_filter_features(D = D_for_analysis(),
                                         filter = eval(parse(text = filter_expr))))
      
      if (identical(D_filter, D_for_analysis())) {
        # Filtration unsuccessful
        showModal(modalDialog(
          title = "Filtration Unsuccessful",
          "No features were filtered based on the selected criteria.",
          easyClose = TRUE
        ))
      } else {
        # Filtration successful
        showModal(modalDialog(
          title = "Filtration Successful",
          "Features were successfully filtered based on the selected criteria.",
          easyClose = TRUE
        ))
      }
    } else if (input$mod2_choice_filter == "mod2_col") {
      selected_values <- input$mod2_filter_col
      
      # Construct logical expression based on selected values
      filter_expr <- paste0("colData(D())$", input$mod2_samp1_select, " %in% c(", paste0("'", selected_values, "'", collapse = ","), ")")
      
      
      # Evaluate the filter expression
      D_filter(mt_modify_filter_samples(D = D_for_analysis(),
                                        filter = eval(parse(text = filter_expr))))
      
      if (identical(D_filter, D_for_analysis())) {
        # Filtration unsuccessful
        showModal(modalDialog(
          title = "Filtration Unsuccessful",
          "No samples were filtered based on the selected criteria.",
          easyClose = TRUE
        ))
      } else {
        # Filtration successful
        showModal(modalDialog(
          title = "Filtration Successful",
          "Samples were successfully filtered based on the selected criteria.",
          easyClose = TRUE
        ))
      }
    }
    
  })
  
  # Retrieve the last available SE object for further analysis
  #D_for_analysis <- reactive({
  #  if (!is.null(D_norm())) {
  #    D_norm()  # Use D_norm if available
  #  } else if (!is.null(D_impute())) {
  #    D_impute()  # Use D_impute if available
  #  } else if (!is.null(D_filter())) {
  #    D_filter()  # Use D_filter if available
  #  } else if (!is.null(D_missingness())) {
  #    D_missingness()  # Use D_missingness if available (fallback option)
  #  } else {
  #    D()
  #  }
  
  #})
  
  
  
  
  
  session_store <- reactiveValues()
  
  
  
  
  
  
  
  # Define rendering logic of control widgets in Module-Annotations Explorer(coded as mod5) ----------------------
  output$mod5_dimension_ui <- renderUI({
    switch(input$mod5_dimension,
           "col"=list(selectInput("mod5_var1_select", 
                                  "Select the primary variable:", 
                                  choices = names(colData(D())),
                                  selected = "Age",
                                  width = "220px"),
                      checkboxInput("mod5_var1_type", 
                                    "Continuous", 
                                    value = TRUE),
                      tags$hr(),
                      selectInput("mod5_var2_select", 
                                  "Select the secondary variable:", 
                                  choices = names(colData(D())),
                                  selected = "sample",
                                  width = "220px"),
                      checkboxInput("mod5_var2_type", 
                                    "Continuous", 
                                    value = TRUE),
                      tags$hr(),
                      selectInput("mod5_select_hover", 
                                  "Select hovering text:", 
                                  choices = names(colData(D())),
                                  selected = names(colData(D()))[1],
                                  width = "220px",
                                  multiple=TRUE)
           ),
           "row"=selectInput("mod5_rowdata_plot", 
                             "Select one plot for row data:", 
                             choices = names(rowData(D())),
                             width = "220px")
    )
  })
  
  output$mod5_output_ui <- renderUI({
    switch(input$mod5_dimension,
           "col"=list(downloadButton("mod5_download_plotly", "download plotly"),
                      plotlyOutput('mod5_plot', height = 600)),
           "row"=list(fluidRow(
             splitLayout(style = "border: 1px", cellWidths = c(1000, 1000), 
                         downloadButton("mod5_download_plotly", "download plotly"), 
                         downloadButton("mod5_download_plotly2", "download plotly")
             )
           ),
           fluidRow(
             splitLayout(style = "height:600px; border: 1px", cellWidths = c(1000, 1000), 
                         plotlyOutput('mod5_plot', height = 600), 
                         plotlyOutput('mod5_plot2', height = 600)
             )
           ))
    )
  })
  
  # Define rendering logic of outputs in Module-Annotations Explorer(coded as mod5) ------------------------------
  mod5_input <- eventReactive(input$mod5_go,{
    c(input$mod5_var1_select,
      input$mod5_var1_type,
      input$mod5_var2_select,
      input$mod5_var2_type,
      input$mod5_rowdata_plot)
  })
  
  
  output$mod5_plot <- renderPlotly({
    tryCatch({
      session_store$mod5_plotly <- switch(input$mod5_dimension,
                                          "col"=
                                            if(mod5_input()[2]==TRUE & mod5_input()[4]==TRUE){
                                              mod5_scatter(D_for_analysis(), x=mod5_input()[3],
                                                           y=mod5_input()[1],
                                                           hover = input$mod5_select_hover)
                                              
                                              
                                            } else if(mod5_input()[2]==TRUE & mod5_input()[4]==FALSE) {
                                              mod5_boxplot(D_for_analysis(), x=mod5_input()[3],
                                                           x_cate = FALSE,
                                                           y=mod5_input()[1],
                                                           y_cate = TRUE,
                                                           fill=mod5_input()[3],
                                                           hover=input$mod5_select_hover)
                                            } else if(mod5_input()[2]==FALSE & mod5_input()[4]==TRUE) {
                                              mod5_boxplot(D_for_analysis(), x=mod5_input()[1],
                                                           x_cate = FALSE,
                                                           y=mod5_input()[3],
                                                           y_cate = TRUE,
                                                           fill=mod5_input()[1],
                                                           hover=input$mod5_select_hover)
                                            } else {
                                              mod5_barplot(D_for_analysis(), x=mod5_input()[3],
                                                           fill=mod5_input()[1],
                                                           hover = input$mod5_select_hover)
                                            },
                                          "row"=
                                            rowData(D()) %>%
                                            data.frame %>%
                                            dplyr::rename(var=mod5_input()[5]) %>%
                                            dplyr::group_by(var) %>%
                                            dplyr::summarise(count=n()) %>%
  
  
  
                                            plot_ly(labels = ~var,
                                                    values = ~count,
                                                    type = 'pie',
                                                    textposition = 'inside',
                                                    source="mod5-click",
                                                    title= sprintf("<b>Distribution of %s</b>",mod5_input()[5])),
                                          layout(autosize = F, width = 1000, height = 500,
                                                 uniformtext=list(minsize=12, mode='hide'),
                                                 legend = list(x = 1,
                                                               y = .5,
                                                               tracegroupgap = 5)
                                          )
      )
      session_store$mod5_plotly
  }, error = function(e) {
    # Handle the error here
    # You can show a notification or error message to the user
    showNotification(e$message, type = "error")
    # You can also print the error message to the console for debugging purposes
    print(paste("Error in renderPlotly:", e))
  })
})  
  
 
  
  # download button
  output$mod5_download_plotly <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      tryCatch({
        saveWidget(as_widget(session_store$mod5_plotly), file, selfcontained = TRUE)
      }, error = function(e) {
        # Handle the error here
        # You can show a notification or error message to the user
        showNotification(e$message, type = "error")
        # You can also print the error message to the console for debugging purposes
        print(paste("Error in downloadHandler:", e))
      })
    }
  )
  
  ## to see the stored data of clicking
  # output$info <- renderPrint({
  #   d5 <- event_data("plotly_click", source = "mod5-click")
  #   if(!is.null(d5)){
  #     d5
  #   }
  # })
  
  output$mod5_plot2 <- renderPlotly({
    
    tryCatch({
      d5 <- event_data("plotly_click", source = "mod5-click")
      pie_dat <- as.data.frame(rowData(D()))
      
      if (!is.null(d5)){
        lvls <- rev(pie_dat$SUPER_PATHWAY)
        label <- lvls[round(as.numeric(d5$pointNumber))+1]
        
        session_store$mod5_plot2 <- 
          pie_dat[pie_dat$SUPER_PATHWAY == label, ] %>%
          dplyr::rename(var="SUB_PATHWAY") %>%
          dplyr::group_by(var) %>%
          dplyr::summarise(count=n()) %>%
          plot_ly(labels = ~var, 
                  values = ~count, 
                  type = 'pie', 
                  textposition = 'inside',
                  title=paste0("<b>Distribution of Sub Pathway in Specified Super Pathway - </b>", label)
          ) %>%
          layout(autosize = F, width = 1000, height = 500,
                 uniformtext=list(minsize=12, mode='hide'),
                 legend = list(x = 1,
                               y = .5,
                               tracegroupgap = 5)
          )
        session_store$mod5_plot2
      }
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in renderPlotly (mod5_plot2):", e))
    })
    
  })
  
  
  # download button
  output$mod5_download_plotly2 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      tryCatch({
        saveWidget(as_widget(session_store$mod5_plotly2), file, selfcontained = TRUE)
      }, error = function(e) {
        # Handle the error here
        # You can show a notification or error message to the user
        showNotification(e$message, type = "error")
        # You can also print the error message to the console for debugging purposes
        print(paste("Error in downloadHandler (mod5_plot2):", e))
      })
    }
  )
  
  
  # Define rendering logic of control widgets in Module-2D projection(coded as mod3) ------------------------
  output$mod3_pca_data <- renderUI({
    if(input$mod3_select_plot=="pca"){
      selectInput("mod3_pca_data_type", "Select data type for PCA:",
                  width = "220px",
                  choices = c("scores", "loadings"),
                  selected = "scores"
      )
    } else {
      NULL
    }
  })
  
  # create intermediate var to indicate coloring widgets
  inter_var <- reactive({
    if (input$mod3_select_plot=="pca" & input$mod3_pca_data_type=="scores") {
      
      "pca-scores"
    } else if(input$mod3_select_plot=="pca" & input$mod3_pca_data_type=="loadings"){
      "pca-loadings"
    } else if(input$mod3_select_plot=="pls"){
      "pls"
    } else {
      "umap"
    }
  })
  
  
  # create reactive plotting argument for PCA/UMAP
  output$mod3_plot_argument <- renderUI({
    switch(
      inter_var(),
      "pca-scores"=list(
        checkboxInput("mod3_scale_data", "Scaled data", 
                      value = TRUE
        ),
        selectInput("mod3_select_colData", 
                    "Select one coloring variable:", 
                    choices = names(colData(D())),
                    selected = "BOX.NUMBER",
                    width = "220px"
        ),
        checkboxInput("mod3_checkbox_factor", 
                      "Categorical Coloring", 
                      value = FALSE
        ),
        selectInput("mod3_select_hover", 
                    "Select hovering text:", 
                    # selectInput coerces its output to character
                    # https://github.com/rstudio/shiny/issues/2367
                    # choices = setNames(seq_along(colData(D)), names(colData(D))),
                    choices = names(colData(D())),
                    selected = "sample",
                    width = "220px",
                    multiple=TRUE
        )
      ),
      "pca-loadings"=list(
        checkboxInput("mod3_scale_data", "Scaled data", 
                      value = TRUE
        ),
        selectInput("mod3_select_colData", 
                    "Select one coloring variable:", 
                    choices = names(rowData(D())),
                    selected = "SUPER_PATHWAY",
                    width = "220px"
        ),
        checkboxInput("mod3_checkbox_factor", 
                      "Categorical Coloring", 
                      value = FALSE
        ),
        selectInput("mod3_select_hover", 
                    "Select hovering text:", 
                    # choices = setNames(seq_along(rowData(D)), names(rowData(D))),
                    choices = names(rowData(D())),
                    selected = "name",
                    width = "220px",
                    multiple=TRUE
        )
      ),
      "umap"=list(numericInput("mod3_umap_n_neighbors", 
                               "Number of neighbors for UMAP:", 
                               value = 15,
                               width = "220px"
      ),
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
      ),
      selectInput("mod3_select_colData", 
                  "Select one coloring variable:", 
                  choices = names(colData(D())),
                  selected = "BOX.NUMBER",
                  width = "220px"
      ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
      ),
      selectInput("mod3_select_hover", 
                  "Select hovering text:", 
                  # choices = setNames(seq_along(colData(D)), names(colData(D))),
                  choices = names(colData(D())),
                  selected = "sample",
                  width = "220px",
                  multiple=TRUE
      )
      ),
      
      "pls"=list(
        #   numericInput("mod3_pls_n_dim",
        #                         "Number of dimensions for PLS:",
        #                         value = 2,
        #                         width = "220px"
        # ),
        
        selectInput("mod3_select_subgroup", 
                    "Select one subgroup variable:", 
                    choices = names(colData(D())),
                    selected = "GROUP_ID",
                    width = "220px"
        ),
        
        selectInput("mod3_select_hover", 
                    "Select hovering text:", 
                    # choices = setNames(seq_along(colData(D)), names(colData(D))),
                    choices = names(colData(D())),
                    selected = "sample",
                    width = "220px",
                    multiple=TRUE
        )
      )
    )
  })
  
  # create reactive inputs list
  mod3_input_object <- eventReactive(input$mod3_go, 
                                     {c(input$mod3_select_plot, 
                                        input$mod3_select_colData,
                                        input$mod3_scale_data,
                                        input$mod3_checkbox_factor,
                                        input$mod3_pca_data_type,
                                        input$mod3_umap_n_neighbors
                                     )}
  )
  
  
  # Define rendering logic of outputs in Module-2D projection(coded as mod3) --------------------------------
  
  # render pca/umap of mod3
  output$mod3_plot <- renderPlotly({
    
    tryCatch({
      session_store$mod3_plotly <- if (mod3_input_object()[1]=="pca"){
        mod3_plots_pca(D = D_for_analysis(),
                       scale_data = mod3_input_object()[3],
                       color = mod3_input_object()[2],
                       categorizing=mod3_input_object()[4],
                       data_type = mod3_input_object()[5],
                       hover = input$mod3_select_hover
        )
      } else if (mod3_input_object()[1]=="umap"){
        print("true")
        mod3_plots_umap(D = D_for_analysis(),
                        scale_data = mod3_input_object()[3],  
                        color = mod3_input_object()[2],
                        categorizing=mod3_input_object()[4],
                        n_neighbors = as.numeric(mod3_input_object()[6]),
                        hover = input$mod3_select_hover
        )
      } else
        mod3_plots_pls(D = D_for_analysis(),
                       subgroupvar  = input$mod3_select_subgroup,
                       # ndim = input$mod3_pls_n_dim,
                       hover = input$mod3_select_hover
                       
                       
        )
      
      session_store$mod3_plotly
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in renderPlotly (mod3_plot):", e))
    })
  })
  
  
  # download button
  output$mod3_download_plotly <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      tryCatch({
        saveWidget(as_widget(session_store$mod3_plotly), file, selfcontained = TRUE)
      }, error = function(e) {
        # Handle the error here
        # You can show a notification or error message to the user
        showNotification(e$message, type = "error")
        # You can also print the error message to the console for debugging purposes
        print(paste("Error in downloadHandler (mod3_plotly):", e))
      })
    }
  )
  
  
  # Define rendering logic of control widgets in Module-All Results Explorer(coded as mod1) ------------------------
  
  obj_name <- reactive({ obj <- get_obj_name(D = D_differ())
  
  obj})
  
  # create stat_name list dependent on radio button
  output$mod1_select_statname_ui <- renderUI({
    selectInput("mod1_select_statname", "Select one stat name:",
                width = "220px",
                
                #choices = obj_name() %>% .[.$V1 == input$mod1_radio, ] %>% dplyr::distinct(stat_name) %>% .$stat_name
                choices = obj_name()$V1 %>% dplyr::distinct()
                
    )
  })
  
  # create object list dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    if (input$mod1_radio=="stats"){
      NULL
    } else {
      selectInput("mod1_select_object", "Select one object:",
                  width = "220px",
                  choices = obj_name() %>% .[.$stat_name == input$mod1_select_statname&.$V1==input$mod1_radio, ] %>% dplyr::distinct(V2) %>% .$V2
                  
      )
    }
    
  })
  
  
  # create indicator of box plot output
  box_switch <- reactive({
    if (input$mod1_select_object=="box"){
      "box_plot"
    } else {
      "non_box_plot"
    }
  })
  
  
  
  ## get the order of selected stat_name
  
  ord <- reactive({
    # assign a data frame of all the object names of box plots
    # filter() cannot run in Shiny, use subset() instead
    box_obj_name <- subset(obj_name(), V1=="plots"&V2=="box")
    box_output_order <- box_obj_name %>%
      dplyr::mutate(order=seq(from=1, to=n()))
    
    if(input$mod1_select_statname %in% box_output_order$stat_name){
      box_output_order[box_output_order$stat_name==input$mod1_select_statname, ]$order
    } else {
      1
    }
  })
  
  # create reactive inputs list
  mod1_input_object <- eventReactive(input$mod1_go, ## delayed output
                                     {c(input$mod1_radio,
                                        input$mod1_select_statname,
                                        input$mod1_select_object)}
  )
  
  
  
  
  # Define rendering logic of outputs in Module-Differential Expression Analysis(coded as mod7) ------------------------------
  
  
  
  
  # control widget of differential analysis
  
  
  
  
  output$mod7_outcome_diff <- renderUI({
    selectInput("outcome_mod7_diff", label = NULL,
                width = "220px",
                choices = names(colData(D())),
                selected = "GROUP_ID"
    )
  })
  
  output$mod7_outcome_corr <- renderUI({
    selectInput("outcome_mod7_corr", label = NULL,
                width = "220px",
                choices = names(colData(D())),
                selected = "GROUP_ID"
    )
  })
  
  output$mod7_outcome_path <- renderUI({
    selectInput("outcome_mod7_path", label = NULL,
                width = "220px",
                choices = names(colData(D())),
                selected = "GROUP_ID"
    )
  })
  
  output$mod7_outcome_path_enrich <- renderUI({
    selectInput("outcome_mod7_path_enrich", label = NULL,
                width = "220px",
                choices = names(colData(D())),
                selected = "GROUP_ID"
    )
  })
  
  output$mod7_stat_name <- renderUI({
    selectInput("stat_name_mod7", label = NULL,
                width = "220px",
                choices = c("GGM","pcor"),
                selected = NULL
    )
  })
  
  
  
  
  output$mod7_filt_samp_dim_diff <- renderUI({
    
    if(input$mod7_filt_samp_diff)
    {
      list( 
        renderUI({selectInput("mod7_samp_filter_diff", label = "Sample Filter Variable:",
                              width = "220px",
                              choices = names(colData(D())),
                              selected = "GROUP_ID"
        )
        }),
        
        renderUI({multiInput("mod7_filter_val_diff", label = input$mod7_samp_filter_diff,
                             choices = colData(D()) %>% as.data.frame() %>% dplyr::select(input$mod7_samp_filter_diff) %>% unique() %>% unlist() %>% as.character(),width = "350px"
                             
        )
        })
        
      )
    }
    
    
  })
  
  samp_filter_diff <- reactive({
    if(input$mod7_filt_samp_diff)
    {
      input$mod7_samp_filter_diff
    }
    else
      NULL
  })
  
  
  
  
  output$mod7_dimension_ui_diff <- renderUI({
    switch(input$mod7_analysis_type_diff,
           "lm"=list(multiInput("mod7_covar_col_select_diff", 
                                label = "Select Covariates from Column Data:", 
                                choices = names(colData(D())),width = "350px"),
                     tags$hr(),
                     multiInput("mod7_covar_row_select_diff", 
                                label = "Select Covariates from Row Data:", 
                                choices = names(rowData(D())),width = "350px")
           )
           
           
    )
    
  })
  
  
  
  
  
  output$mod7_group_col_barplot_diff <- renderUI({
    selectInput("group_col_barplot_mod7_diff", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod7_group_col_barplot_corr <- renderUI({
    selectInput("group_col_barplot_mod7_corr", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod7_group_col_barplot_path <- renderUI({
    selectInput("group_col_barplot_mod7_path", label = NULL,
                width = "220px",
                selected=NULL,
                choices = names(rowData(D()))
    )
  })
  
  output$mod7_color_col_barplot <- renderUI({
    selectInput("color_col_barplot_mod7", label = NULL,
                width = "220px",
                selected="BIOCHEMICAL",
                choices = names(rowData(D()))
    )
  })
  
  output$mod7_in_col <- renderUI({
    selectInput("in_col_mod7", label = NULL,
                width = "220px",
                choices = names(rowData(D())),
                selected = "GROUP_ID"
    )
  })
  
  # Define rendering logic of outputs in Module-Differential Expression Analysis(coded as mod7) ------------------------------
    D_download <- reactiveVal()
  
    observeEvent(input$mod7_go_diff,{

      tryCatch({

        print("done")
        D_download(D_differ_tab_diff())
        print(class(D_differ_tab_diff()))
        
  
        covar_col_select <- input$mod7_covar_col_select_diff
        covar_row_select <- input$mod7_covar_row_select_diff
        var <- input$outcome_mod7_diff
        if(is.null(covar_col_select) && is.null(covar_row_select))
        {
          covar = NULL
        } else
        {
          covar <- paste("+", paste(c(covar_row_select,covar_col_select), collapse = "+"), sep = "")
        }
        
        x <- if (input$mod7_outcome_binary_diff) {
          "fc"
        } else {
          "statistic"
        }
        
        p1 <- mt_plots_volcano_new(
          D_differ_tab_diff(),
          stat_name = sprintf("~  %s%s Analysis", var, replace(covar, is.null(covar), "")),
          x = !!sym(x),
          feat_filter = p.adj < input$mod7_sig_threshold_diff,
          color = p.adj < input$mod7_sig_threshold_diff
        )
        
        interactive_plot1 <- plotly::ggplotly(p1)
        
        interactive_plot1 <- plotly::style(interactive_plot1, tooltip = input$mod5_select_hover)
        interactive_plot1 <- plotly::event_register(interactive_plot1, "plotly_hover") 
        
   
        
        output$mod7_main_panel_graph1 <- renderUI({
          plotly::plotlyOutput("interactive_plot1") 
        })
        
        output$interactive_plot1 <- plotly::renderPlotly({
          interactive_plot1
        })
        
        
        plot_type <- if (input$mod7_outcome_binary_diff) {
          "box"
        } else {
          "scatter"
        }
        
        p2 <- mt_plots_box_scatter_new(D_differ_tab_diff(),
                                       stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                                       x = !!sym(var),
                                       plot_type = plot_type,
                                       feat_filter = p.adj < input$mod7_sig_threshold_diff,
                                       feat_sort = p.value,
                                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}")
        
        interactive_plot2 <- plotly::ggplotly(p2)
        
        interactive_plot2 <- plotly::style(interactive_plot2, tooltip = input$mod5_select_hover)
        interactive_plot2 <- plotly::event_register(interactive_plot2, "plotly_hover")
        
  
        output$mod7_main_panel_graph2 <- renderUI({
          plotly::plotlyOutput("interactive_plot2") 
        })
        
        output$interactive_plot2 <- plotly::renderPlotly({
          interactive_plot2
        })
        
       
        
        p3 <- mt_plots_stats_pathway_bar_new(D_differ_tab_diff(),
                                             stat_list = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                                             y_scale = "count",
                                             feat_filter = p.adj < input$mod7_sig_threshold_diff,
                                             group_col = input$group_col_barplot_mod7_diff)
        
        interactive_plot3 <- plotly::ggplotly(p3)
        
        interactive_plot3 <- plotly::style(interactive_plot3, tooltip = input$mod5_select_hover)
        interactive_plot3 <- plotly::event_register(interactive_plot3, "plotly_hover")
        
        output$mod7_main_panel_graph3 <- renderUI({
          plotlyOutput("interactive_plot3") 
        })
        
        output$interactive_plot3 <- renderPlotly({
          interactive_plot3
        })
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in observeEvent (mod7_go_diff):", e))
    })
     
  
  })
  
  observeEvent(input$mod7_go_corr,{
    
    tryCatch({
    
      print("done")
      print(class(D_differ_tab_corr()))
      D_download(D_differ_tab_corr())
      
      #covar_col_select <- input$mod7_covar_col_select_diff
      #covar_row_select <- input$mod7_covar_row_select_diff
      var <- input$outcome_mod7_corr
      covar = NULL
       
      
      p1 <- mt_plots_volcano_new(D_differ_tab_corr(),
                                 stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                                 x = !!sym(ifelse(input$mod7_outcome_binary_corr,"fc","statistic")),
                                 feat_filter = p.adj < input$mod7_sig_threshold_corr,
                                 color = p.adj < input$mod7_sig_threshold_corr)
  
  
  
  
      interactive_plot1 <- plotly::ggplotly(p1)
  
      interactive_plot1 <- plotly::style(interactive_plot1, tooltip = input$mod5_select_hover)
      interactive_plot1 <- plotly::event_register(interactive_plot1, "plotly_hover")
  
      output$mod7_main_panel_graph1 <- renderUI({
        plotly::plotlyOutput("interactive_plot1")
      })
  
      output$interactive_plot1 <- plotly::renderPlotly({
        interactive_plot1
      })
      
      
      p2 <- mt_plots_box_scatter_new(D_differ_tab_corr(),
                                     stat_name = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                                     x = !!sym(var),
                                     plot_type = ifelse(input$mod7_outcome_binary_corr,"box","scatter"),
                                     feat_filter = p.adj < input$mod7_sig_threshold_corr,
                                     feat_sort = p.value,
                                     annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}")
      
      interactive_plot2 <- plotly::ggplotly(p2)
      
      interactive_plot2 <- plotly::style(interactive_plot2, tooltip = input$mod5_select_hover)
      interactive_plot2 <- plotly::event_register(interactive_plot2, "plotly_hover")
      
      output$mod7_main_panel_graph2 <- renderUI({
        plotly::plotlyOutput("interactive_plot2")
      })
      
      output$interactive_plot2 <- plotly::renderPlotly({
        interactive_plot2
      })
      
      p3 <- mt_plots_stats_pathway_bar_new(D_differ_tab_corr(),
                                           stat_list = sprintf("~  %s%s Analysis",var, replace(covar, is.null(covar),"")),
                                           y_scale = "count",
                                           feat_filter = p.adj < input$mod7_sig_threshold_corr,
                                           group_col = input$group_col_barplot_mod7_corr)
      
      interactive_plot3 <- plotly::ggplotly(p3)
      
      interactive_plot3 <- plotly::style(interactive_plot3, tooltip = input$mod5_select_hover)
      interactive_plot3 <- plotly::event_register(interactive_plot3, "plotly_hover")
      
      output$mod7_main_panel_graph3 <- renderUI({
        plotly::plotlyOutput("interactive_plot3")
      })
      
      output$interactive_plot3 <- plotly::renderPlotly({
        interactive_plot3
      })
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in observeEvent (mod7_go_corr):", e))
    })
      
    
  })
  
  
  #actions for correlation analysis:
  # observeEvent(input$mod7_go_corr,{
  #   
  #   
  #   pdiffer_list$value <- get_plots_SE_differ(D_differ_tab_corr())
  #   pdiffer_list$length <-  pdiffer_list$value %>% length()
  #   
  #   output$mod7_main_panel  <- renderUI({
  #     
  #     mod7_output_plotlist <-   lapply(1: pdiffer_list$length, function(i){
  #       local({
  #         len_j <- length(pdiffer_list$value[[i]])
  #         lapply(1:(len_j), function(j) {
  #           
  #           plotname <- paste("Plot_differ", i,j, sep="")
  #           
  #           plotOutput(plotname)
  #           
  #         })
  #       })
  #       
  #       
  #     })
  #     
  #     do.call(tagList, mod7_output_plotlist)
  #     
  #     
  #   })
  #   
  #   lapply(1: pdiffer_list$length, function(i){
  #     local({
  #       
  #       len_j <- length(pdiffer_list$value[[i]])
  #       
  #       lapply(1:(len_j), function(j) {
  #         
  #         plotname <- paste("Plot_differ", i,j, sep="")
  #         
  #         output[[plotname]] <-
  #           renderPlot({
  #             #grid.force()
  #             pdiffer_list$value[[i]][j]
  #             
  #             #print(pdiffer_list$value[[i]][j])
  #             
  #           })
  #       })
  #     })
  #     
  #   })
  #   
  #   
  # })
  
  
  
  
  
  #actions for partial correlation
  observeEvent(input$mod7_go_partial_corr,{
    
   
    
    # Wrap the code in tryCatch to handle errors gracefully
    tryCatch({
    
      # Get the pathway enrichment results
      output_table <- mt_stats_cormat_genenet_new(D_for_analysis(),
                                                  stat_name = input$stat_name_mod7)
      # Display the table
      print(output_table)
      
      output$mod7_main_panel_2 <- renderDT({
        datatable(output_table)
      })
    }, error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        glue::glue("An error occurred: {e$message}"),
        type = "error"
      )
    })
    
    
    # pdiffer_list$value <- get_plots_SE_differ(D_differ_partial_corr())
    # print(pdiffer_list$value)
    # pdiffer_list$length <-  pdiffer_list$value %>% length()
    # 
    # output$mod7_main_panel  <- renderUI({
    #   
    #   mod7_output_plotlist <-   lapply(1: pdiffer_list$length, function(i){
    #     local({
    #       len_j <- length(pdiffer_list$value[[i]])
    #       lapply(1:(len_j), function(j) {
    #         
    #         plotname <- paste("Plot_differ", i,j, sep="")
    #         
    #         plotOutput(plotname)
    #         
    #       })
    #     })
    #     
    #     
    #   })
    #   
    #   do.call(tagList, mod7_output_plotlist)
    #   
    #   
    # })
    # 
    # lapply(1: pdiffer_list$length, function(i){
    #   local({
    #     
    #     len_j <- length(pdiffer_list$value[[i]])
    #     
    #     lapply(1:(len_j), function(j) {
    #       
    #       plotname <- paste("Plot_differ", i,j, sep="")
    #       
    #       output[[plotname]] <-
    #         renderPlot({
    #           #grid.force()
    #           pdiffer_list$value[[i]][j]
    #           
    #           #print(pdiffer_list$value[[i]][j])
    #           
    #         })
    #     })
    #   })
    #   
    # })
    
  })
  
  
  #actions for pathway enrichment
  observeEvent(input$mod7_go_path_enrich, {
    tryCatch({
      # Get the pathway enrichment results
      D_download(D_final_enrich())
      print(D_final_enrich())
      enrichment_results <- metadata(D_final_enrich())$pathways$enrichment_results
      
      print(class(enrichment_results))
      # Display the results in the desired output element
      output$mod7_main_panel_2 <- renderDT({
        datatable(enrichment_results)
      })
    }, error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        glue::glue("An error occurred: {e$message}"),
        type = "error"
      )
    })
    
    
    #want barplot for D_final_enrich
    # Generate the stat_list based on input values
    #covar_col_select <- input$mod7_covar_col_select_diff
    #covar_row_select <- input$mod7_covar_row_select_diff
    var <- input$outcome_mod7_path_enrich
    #if(is.null(covar_col_select) && is.null(covar_row_select))
    #{
    #  covar = NULL
    #} else
    #{
    #  covar <- paste("+", paste(c(covar_row_select,covar_col_select), collapse = "+"), sep = "")
    #}
    
    
    # p1 <- mt_plots_stats_pathway_bar_new(
    #   D_final_enrich(),
    #   stat_list = stat_list,
    #   feat_filter = p.value < 1,
    #   group_col = input$in_col_mod7
    # )
    # result <- get_stat_names(D_final_enrich)
    # 
    # # Print the result
    # print(result)
    # 
    # 
    # 
    # 
    # p1 <- mt_plots_stats_pathway_bar_new(D_final_enrich(),
    #                                      stat_list = stat_list, #sprintf("~  %s%s",input$outcome_mod7_path_enrich, ""),
    #                                      feat_filter = p.adj < input$mod7_sig_threshold_diff)
    # 
    # 
    # 
    # 
    # 
    # interactive_plot1 <- plotly::style(interactive_plot1, tooltip = input$mod5_select_hover)
    # interactive_plot1 <- plotly::event_register(interactive_plot1, "plotly_hover")
    # 
    # output$mod7_main_panel_1 <- renderUI({
    #   plotly::plotlyOutput("interactive_plot1")
    # })
    # 
    # output$interactive_plot1 <- plotly::renderPlotly({
    #   interactive_plot1
    # })
  })
  
  
  # Differential analysis D
  # D_differ_tab_diff <- reactive({
  #   
  #   
  #   D <- D_for_analysis()  %>%
  #     mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
  #     diff_analysis_func_tab(var=input$outcome_mod7_diff,
  #                            binary=input$mod7_outcome_binary_diff,
  #                            sample_filter = samp_filter_diff(),
  #                            filter_val = input$mod7_filter_val_diff,
  #                            covar_col_select = input$mod7_covar_col_select_diff,
  #                            covar_row_select = input$mod7_covar_row_select_diff,
  #                            analysis_type=input$mod7_analysis_type_diff,
  #                            mult_test_method=input$mod7_mult_test_method_diff,
  #                            alpha=input$mod7_sig_threshold_diff,
  #                            group_col_barplot = input$group_col_barplot_mod7_diff) %>%
  #     {.}
  #   ## return D
  #   D
  # }) 
  
  #SE object withSpinner:
  D_differ_tab_diff <- reactive({
    tryCatch({
      D_for_analysis()  %>%
        mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
        diff_analysis_func_tab(var = input$outcome_mod7_diff,
                               binary = input$mod7_outcome_binary_diff,
                               sample_filter = samp_filter_diff(),
                               filter_val = input$mod7_filter_val_diff,
                               covar_col_select = input$mod7_covar_col_select_diff,
                               covar_row_select = input$mod7_covar_row_select_diff,
                               analysis_type = input$mod7_analysis_type_diff,
                               mult_test_method = input$mod7_mult_test_method_diff,
                               alpha = input$mod7_sig_threshold_diff,
                               group_col_barplot = input$group_col_barplot_mod7_diff)
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in D_differ_tab_diff:", e))
      # Return an appropriate value in case of an error, for example:
      # return(data.frame()) # Return an empty data frame
      return(NULL) # Return NULL
    })
  })
  
  
  
  
  # Correlation analysis D
  D_differ_tab_corr <- reactive({
    tryCatch({
      # Differential analysis D
      D_for_analysis()  %>%
        mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
        diff_analysis_func_tab(var = input$outcome_mod7_corr,
                               binary = input$mod7_outcome_binary_corr,
                               analysis_type = input$mod7_analysis_type_corr,
                               mult_test_method = input$mod7_mult_test_method_corr,
                               alpha = input$mod7_sig_threshold_corr,
                               group_col_barplot = input$group_col_barplot_mod7_corr) %>%
        {.}
    }, error = function(e) {
      # Handle the error here
      # You can show a notification or error message to the user
      showNotification(e$message, type = "error")
      # You can also print the error message to the console for debugging purposes
      print(paste("Error in D_differ_tab_corr:", e))
      # Return an appropriate value in case of an error, for example:
      # return(data.frame()) # Return an empty data frame
      return(NULL) # Return NULL
    })
  })
  
  
  # Pathway analysis D
  # D_differ_tab_path <- reactive({
  #   
  #   
  #   D <- D_for_analysis()  %>%
  #     mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
  #     diff_analysis_func_tab(var=input$outcome_mod7_path,
  #                            alpha=input$mod7_sig_threshold_path,
  #                            group_col_barplot=input$group_col_barplot_mod7_path,
  #                            color_col_barplot = input$color_col_barplot_mod7) %>%
  #     {.}
  #   ## return D
  #   D
  # })
  
  # Partial Correlation Network D
  D_differ_partial_corr <- reactive({
    result <- tryCatch({
      D <- D_for_analysis() %>%
        mt_reporting_heading(heading = "Partial Correlation Network", lvl = 2) %>%
        mt_stats_cormat_genenet(stat_name = input$stat_name_mod7) %>%
        mt_post_multtest(stat_name = input$stat_name_mod7, method = "BH") %>%
        {.}
      ## return D
      D
    }, error = function(e) {
      # Handle the error here
      showNotification(e$message, type = "error")
      message(paste("Error in D_differ_partial_corr:", conditionMessage(e))) # Print error message to console
      return(NULL) # Return NULL
    })
    
    return(result)
  })
  
  
  
  #Pathway enrichment D
  D_path_enrich <- reactive({
    
    
    D <- D_for_analysis() %>%
      mt_anno_pathways_hmdb_new(in_col = input$in_col_mod7,
                                out_col = "kegg_db",
                                pwdb_name = "KEGG",
                                db_dir = system.file("extdata", "precalc/hmdb", package = "maplet")) %>%
      mt_anno_pathways_remove_redundant(feat_col = "HMDb", pw_col = "kegg_db") %>%
    D
  })
  
  D_final_enrich <- reactive({
    
    tryCatch({
      D_path_enrich() %>%
        mt_stats_univ_lm(formula = as.formula(sprintf("~ %s%s",input$outcome_mod7_path_enrich, "")), stat_name = sprintf("~  %s%s",input$outcome_mod7_path_enrich, "")) %>%
        mt_post_multtest(stat_name = sprintf("~  %s%s",input$outcome_mod7_path_enrich, ""), method = "BH") %>%
        mt_post_fold_change(stat_name = sprintf("~  %s%s",input$outcome_mod7_path_enrich, "")) %>%
        mt_stats_pathway_enrichment(pw_col = "kegg_db",
                                        stat_name = sprintf("~  %s%s",input$outcome_mod7_path_enrich, ""),
                                        cutoff = 0.4)
    }, error = function(e) {
      # Show the error message to the user using showNotification
      showNotification(
        e$message,
        type = "error"
      )
      # Return an empty data frame or NULL to avoid other issues in the app
      NULL
    })
  })
  
  
  generate_report <- reactive({
    
    D_differ_tab_diff() %>% mt_reporting_html(file = "Example_Pipeline_Metabolite_Analysis.html",
                                              title = "Example Pipeline - Statistical Analysis")
    
    return("Report generation successful!")
  })
  
  # Provide a link to download the report (Code for downloadHandler)
  # output$download_link <- downloadHandler(
  #   filename = function() {
  #     "Example_Pipeline_Metabolite_Analysis.html"
  #   },
  #   content = function(file) {
  #     # Create a list to hold the SE objects that have been created and are not empty
  #     valid_D_list <- list()
  #     
  #     # Check and include only the valid SE objects in the list
  #     if (!is.null(D_differ_tab_diff()) && nrow(D_differ_tab_diff()) > 0) {
  #       valid_D_list$D_differ_tab_diff <- D_differ_tab_diff()
  #     }
  #     
  #     if (!is.null(D_differ_tab_corr()) && nrow(D_differ_tab_corr()) > 0) {
  #       valid_D_list$D_differ_tab_corr <- D_differ_tab_corr()
  #     }
  #     
  #     if (!is.null(D_final_enrich()) && nrow(D_final_enrich()) > 0) {
  #       valid_D_list$D_final_enrich <- D_final_enrich()
  #     }
  #     
  #     # Check if the valid_D_list is not empty before generating the report
  #     if (length(valid_D_list) > 0) {
  #       # Generate the HTML content using mt_reporting_html_nonlinear
  #       html_content <- mt_reporting_html_nonlinear(D_list = valid_D_list, title = "Example Pipeline - Statistical Analysis")
  #       
  #       # Write the HTML content to the file
  #       cat(html_content, file = file)
  #     } else {
  #       # If valid_D_list is empty, show an error message to the user
  #       showNotification("No valid data found for generating the report.", type = "error")
  #     }
  #   }
  # )
  
  
  output$download_link <- downloadHandler(
    filename = function() {
      "Example_Pipeline_Metabolite_Analysis.html"
    },
    content = function(file) {
      
      D_diff <- D_differ_tab_diff()
      D_corr <- D_differ_tab_corr()
      D_partial_corr <- D_differ_partial_corr()
      
      print(class(D_diff))
      print(class(D_corr))
      print(class(D_partial_corr))
      
      # Call mt_reporting_html and save the report with the specified filename
      mt_reporting_html_nonlinear(D_list = list(D_diff, D_corr, D_partial_corr), file = file, 
                            title = "Example Pipeline - Statistical Analysis")
      # Note: Replace D with your actual SummarizedExperiment object.
      
    }
  )
  
  
  
  # output$download_link <- downloadHandler(
  #   filename = function() {
  #     "Example_Pipeline_Metabolite_Analysis.html"
  #   },
  #   content = function(file) {
  #     # Call mt_reporting_html and save the report with the specified filename
  #     mt_reporting_html_nonlinear(D_list = list(D_differ_tab_diff(),D_differ_tab_corr(),D_final_enrich()), file = file, title = "Example Pipeline - Statistical Analysis")
  #     # Note: Replace D with your actual SummarizedExperiment object.
  #   }
  # ) 
  
  
  
  # Define rendering logic of outputs in Module-All Results Explorer(coded as mod1) --------------------------------
  
  # Insert the right number of plot output objects into UI
  
  
  # # define pathway annotation column (extracted from corresponding stat_bar
  # pwvar <- mtm_res_get_entries(D_results, c("plots", "stats"))[[1]]$args$group_col
  # # define threshold for significance (extracted from corresponding stat_bar plot)
  # alpha <- mtm_res_get_entries(D_results, c("plots", "stats"))[[1]]$args$feat_filter[[3]]
  # # get pathway annotations
  # rd <- get_pathway_annotations(D_results, pwvar)
  
  output$mod1_output_plot <- renderUI({
    ## limit plots to specified stat_name
    obj_name <- subset(obj_name(), V1==mod1_input_object()[1])
    obj_name <- subset(obj_name(), V2==mod1_input_object()[3])
    output_order <- obj_name() %>%
      dplyr::mutate(order=seq(from=1, to=n()))
    output_order <- subset(output_order, stat_name==mod1_input_object()[2])
    plots <- list()
    for(plot_i in seq_along(output_order$order)){
      plots[[plot_i]] <- mtm_res_get_entries(D = D_differ(), c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]
    }
    # there are multiple plots
    len_i <- length(plots)
    # some plots have multiple objects
    len_j <- length(plots[[1]]$output)
    # name every plot object in UI
    mod1_plot_output_list <- lapply(1:(len_i*len_j), function(i) {
      plotname <- paste("Plot", i, sep="")
      # locate the row in the `plots`
      row_n <- ceiling(i/len_j)
      ## set dynamic height of box scatter plots based on output2
      height <- if(plots[[1]]$fun[2]=="box"&plots[[1]]$fun[3]=="scatter"&!is.null(plots[[row_n]]$output2)){
        as.numeric(plots[[row_n]]$output2)*150
      } else {
        560
      }
      plotOutput(plotname, height = height, width = 850)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, mod1_plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  # get the max number of objects
  
  
  max_plot <- eventReactive(input$mod1_go,{
    
    num_df <-  obj_name() %>% subset(V1=="plots")
    num_df <- num_df %>%
      dplyr::group_by(V2, stat_name) %>%
      dplyr::summarise(cnt_sum=sum(cnt))
    max_plot <- max(num_df$cnt_sum)
    max_plot
    
    
    
    
    for (i in 1:max_plot) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i
        plotname <- paste("Plot", my_i, sep="")
        output[[plotname]] <- renderPlot({
          ## limit plots to specified stat_name
          obj_name <-reactive({ subset(obj_name(), V1==mod1_input_object()[1])})
          obj_name <- reactive({ subset(obj_name(), V2==mod1_input_object()[3])})
          output_order <- obj_name %>%
            dplyr::mutate(order=seq(from=1, to=n()))
          output_order <- subset(output_order, stat_name==mod1_input_object()[2])
          plots <- list()
          for(plot_i in seq_along(output_order$order)){
            plots[[plot_i]] <- mtm_res_get_entries(D = D_differ(), c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]$output
          }
          # there are multiple plots
          len_i <- length(plots)
          # some plots have multiple objects
          len_j <- length(plots[[1]])
          # locate the row in the `plots`
          row_n <- ceiling(my_i/len_j)
          # locate the column in the `plots`
          col_n <- ifelse((my_i %% len_j)==0, len_j, (my_i %% len_j))
          # render the plot object in each loop
          plots[[row_n]][col_n]
          
        })
      })
    }
    
  })
  # render stats table of Mod1
  output$mod1_output_table <- renderDataTable({
    table <- data.frame(var=row.names(rowData( D_differ())), rowData(D_differ())) %>%
      left_join(mtm_get_stat_by_name(D = D_differ(), mod1_input_object()[2]),
                by=c("var"="var")
      ) %>%
      dplyr::select(c(2, 20:26))
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        dplyr::select(name, statistic, p.value, p.adj, term, dplyr::everything()) %>%
        ## scientific notation
        dplyr::mutate(statistic=formatC(statistic, format = "E", digits = 2),
                      p.value=formatC(p.value, format = "E", digits = 2),
                      p.adj=formatC(p.adj, format = "E", digits = 2),
                      estimate=formatC(estimate, format = "E", digits = 2),
                      std.error=formatC(std.error, format = "E", digits = 2)
        )
    } else {
      table %>%
        dplyr::select(name, statistic, p.value, p.adj, dplyr::everything())
    }
    datatable(table,
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                ## set column width
                autoWidth = TRUE,
                columnDefs = list(list(width = '100px', targets = c(2:4))),
                scrollX = TRUE
              ))
  })
  # render plots or table
  output$mod1_output <- renderUI({
    switch(
      mod1_input_object()[1],
      "plots" = uiOutput("mod1_output_plot"),
      "stats" = dataTableOutput("mod1_output_table")
    )
  })
  
  
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)


