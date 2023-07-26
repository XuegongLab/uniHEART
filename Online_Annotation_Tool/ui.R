library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(markdown)
library(shinyalert)
library(igraph)
library(reshape2)
library(circlize)
library(Seurat)
library(SingleR)
library(future)


shinyUI(fluidPage(theme = shinytheme("cerulean"),
                  
                  navbarPage("Human Heart Cell Atlas",
                             
                             #################### Annotation ####################
                             tabPanel(p(icon("feather-alt",lib = "font-awesome"),em("Annotation"),style="text-align:center"),
                                      br(),
                                      h1("Online analysis",style="text-align:center"),
                                      p("You can upload your data to conduct online analysis and annotation. The uploaded file should be a expression matrix in .csv or seurat object file in .rds. The uploaded file should be <= 200MB an <= 10,000 cells."),
                                      fluidRow(
                                        column(4,radioButtons("PA_stage", "Choose the development stage of your data:",c("Adult", "Fetal")),offset=2),
                                        column(6,
                                               fileInput("data_upload", "Choose your data file:", multiple = FALSE, accept = c(".csv",".rds",'.tsv')),
                                               actionButton(inputId="BTN_Demo_Dataset",label= "Using Demo Dataset",icon=icon("book-open",lib = "font-awesome")),
                                               checkboxInput("Standardize_Gene_CKB", "Do Gene Standardization", FALSE),
                                               p("[Warning]: The gene standardization may take about one hour, we suggest do it offline with our package",a("hECA_GeneSymbolUniform_Rtoolkit",href="https://github.com/XuegongLab/hECA_GeneSymbolUniform_Rtoolkit"))
                                        )
                                      ),
                                      fluidRow(column(verbatimTextOutput("TO_Standardize"),style="text-align:center",width = 12)),
                                      fluidRow(column(verbatimTextOutput("TO_Load"),style="text-align:center",width = 12)),
                                      hr(),
                                      
                                      #################### QC ####################
                                      wellPanel(
                                        fluidRow(
                                          column(12,h3("Quality control"),style="text-align:center")
                                        ),
                                        fluidRow(
                                          column(3,h4("QC parameters:"),
                                                 p("nFeature range:"),
                                                 fluidRow(
                                                    column(numericInput("nF_min", label = NULL, value = 200),width = 5),
                                                    column(p("~"),width = 1),
                                                    column(numericInput("nF_max", label = NULL,value = 7500),width = 5)
                                                 ),
                                                 p("nCount range:"),
                                                 fluidRow(
                                                   column(numericInput("nC_min", label = NULL, value = 0),width = 5),
                                                   column(p("~"),width = 1),
                                                   column(numericInput("nC_max", label = NULL,value = 500000),width = 5)
                                                 ),
                                                 p("MT percentage range:"),
                                                 fluidRow(
                                                   column(numericInput("MT_min", label = NULL, value = 0),width = 5),
                                                   column(p("~"),width = 1),
                                                   column(numericInput("MT_max", label = NULL,value = 15),width = 5)
                                                 ),
                                                 actionButton("BTN_MAP","Start",icon=icon("clipboard-check",lib = "font-awesome"))
                                                 ,style="text-align:center"),
                                          column(3,plotOutput("PO_nFeature")),
                                          column(3,plotOutput("PO_nCount")),
                                          column(3,plotOutput("PO_MT"))
                                        ),
                                        br(),
                                        fluidRow(column(verbatimTextOutput("TO_QCResult"),style="text-align:center",width = 12))
                                      ),
                                      hr(),

                                      #################### Visualization ####################
                                      wellPanel(
                                        fluidRow(
                                          column(12,h3("Visualization"),style="text-align:center")
                                        ),
                                        fluidRow(column(verbatimTextOutput("TO_Preprocess1"),style="text-align:center",width = 12)),
                                        fluidRow(column(verbatimTextOutput("TO_Preprocess2"),style="text-align:center",width = 12)),
                                        fluidRow(
                                          column(6,plotOutput("PO_UMAP1",height=800),style="text-align:center",offset=3)
                                        )
                                      ),
                                      hr(),

                                      #################### Annotation ####################
                                      wellPanel(
                                        fluidRow(
                                          column(12,h3("Annotation"),style="text-align:center")
                                        ),
                                        fluidRow(
                                          column(4,"Annoation method:",style="text-align:right"),
                                          column(4,selectInput("select_AM", label = NULL, 
                                                                choices = list("singleR" = 1,"Seurat.Mapping"=2), 
                                                                selected = 2)
                                                ),
                                          column(4,actionButton("BTN_ANNO","Start",icon=icon("clipboard-check",lib = "font-awesome")))
                                        ),
                                        fluidRow(column(verbatimTextOutput("TO_Annotation"),style="text-align:center",width = 12)),
                                        fluidRow(
                                          column(8,plotOutput('PO_UMAP2',height=800)),
                                          column(4,plotOutput('PO_hist',height=800))
                                        ),
                                        br(),
                                        fluidRow(
                                          column(4,verbatimTextOutput('TO_Proportion')),
                                          column(8,plotOutput('PO_Pie',height=800))
                                        )
                                      ),
                                      p(downloadButton("BTN_DOWNLOAD1","Download result",icon = shiny::icon("download")),style="text-align:right"),
                                      hr(),

                                      #################### Feature Plot ####################
                                      wellPanel(
                                        fluidRow(
                                          column(12,h3("Feature Plot"),style="text-align:center")
                                        ),
                                        fluidRow(
                                          column(4,"Gene to view:"),
                                          column(4,textInput("select_FP", label = NULL,value = "ACTB")),
                                          column(4,actionButton("BTN_FP","Plot",icon=icon("clipboard-check",lib = "font-awesome")))
                                        ),
                                        fluidRow(
                                          column(6,plotOutput("PO_FeaturePlot",height=800),style="text-align:center",offset=3)
                                        ),
                                        fluidRow(
                                          column(6,plotOutput("PO_VlnPlot",height=800),style="text-align:center",offset=3)
                                        )
                                      ),
                                      hr(),

                                      #################### DEgenes ####################
                                      wellPanel(
                                        fluidRow(
                                          column(12,h3("DEGs of the annotation results"),style="text-align:center")
                                        ),
                                        fluidRow(
                                          DT::dataTableOutput("DEG_table")
                                        )
                                      ),
                                      p(downloadButton("BTN_DOWNLOAD2","Download DEGs",icon = shiny::icon("download")),style="text-align:right"),
                             )),
                  hr(),
                  p(em("Developed by Yixin Chen"),br("XGlab, Tsinghua University, Beijing, China"),style="text-align:center; font-family: times")
))
