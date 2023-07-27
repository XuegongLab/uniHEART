library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(shinyalert)
library(igraph)
library(reshape2)
library(Seurat)
library(SingleR)
library(future)

library(stringr)
library(dplyr)
library(data.table)
library(rlist)

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# gene uniform tool
Gene_Uniform <- function(query_data){
  # load reference gene
  ref_table_raw <- read.csv("www/GeneSymbolRef_SelectAll_upd0731.csv", header=TRUE, na.strings=TRUE, stringsAsFactors=FALSE)
  # Separate dataframe to previous and alias symbol sub-dataframe
  # Remove duplicates and empty values
  ref_table_raw <- ref_table_raw[,c("Approved.symbol","Alias.symbol","Previous.symbol")]
  ref_table <- ref_table_raw[ref_table_raw[,"Previous.symbol"]!="" | ref_table_raw[,"Alias.symbol"]!="",]
  # Seurat changes all "_" to "-".
  ref_table$Previous.symbol <- str_replace(ref_table$Previous.symbol, "_", "-")
  ref_table$Alias.symbol <- str_replace(ref_table$Alias.symbol, "_", "-")
  ref_table$Approved.symbol <- str_replace(ref_table$Approved.symbol, "_", "-")
  # print(dim(ref_table))
  ref_table_prev <- unique(ref_table[,c("Approved.symbol","Previous.symbol")])
  ref_table_prev <- ref_table_prev[ref_table_prev[,"Previous.symbol"]!="",]
  ref_table_alia <- unique(ref_table[,c("Approved.symbol","Alias.symbol")])
  ref_table_alia <- ref_table_alia[ref_table_alia[,"Alias.symbol"]!="",]
  
  query_gene_list <- rownames(query_data)

  total_gene_list_raw = read.table("www/total_gene_list_43878.txt", 
                     header=TRUE, sep='\t', fill=TRUE, stringsAsFactors=FALSE)
  total_gene_list = total_gene_list_raw[,1]
  # Seurat changes all "_" to "-".
  total_gene_list <- str_replace(total_gene_list, "_", "-")


  setDT(result_data)
  colstoavg <- names(result_data)[1:(dim(result_data)[2]-1)]
  result_data_grouped <- result_data[,lapply(.SD, mean, na.rm=TRUE),by=genenames,.SDcols=colstoavg]

  result_data_sub <- as.data.frame(result_data_grouped)[which(!result_data_grouped$genenames %in% outlier_gene_list),]
  result_data_out <- subset(result_data_sub, select = -genenames )
  rownames(result_data_out) <- result_data_sub$genenames
  print("Shape of processed query data: ")
  print(dim(result_data_out))

  add_df <- data.frame(matrix(nrow=sum(gene_appearance_list$appearance==FALSE),ncol=dim(result_data)[2]-1, 0))
  rownames(add_df) <- gene_appearance_list$gene_name[!gene_appearance_list$appearance]
  colnames(add_df) <- colnames(result_data_out)

  result_data_out <- rbind(result_data_out, add_df)
  return(result_data_out)
}

# limit the upload data size to 200M
options(shiny.maxRequestSize=200*1024^2)
plan("multicore", workers = 20)

shinyServer(function(input, output){
    
    
    # global_variables
    app.env <- reactiveValues(
      path = NULL,
      Query.raw = NULL,
      Query.processed = NULL,
      Pre.result = NULL,
      DEGs = NULL
    )

    # load data
    observeEvent(
      eventExpr = input$data_upload,
      handlerExpr = {
        if (nchar(x = input$data_upload$datapath)) {
          app.env$path <- input$data_upload$datapath
        }
      }
    )

    observeEvent(
      eventExpr = app.env$path,
      handlerExpr = {
        file_name <- app.env$path
        if(strsplit(file_name,"[.]")[[1]][-1]=="rds"){
          dataobj <- readRDS(file_name)
          output$TO_Load <- renderPrint({
              print("A .rds file loaded.")
              print(paste0(dim(dataobj)[2]," cells and ",dim(dataobj)[1]," genes loaded."))
              print("Print out first 5 genes in query data, in case something wrong happens in data loading: ")
              print(row.names(dataobj)[1:5])
          })

        }else{
          if(strsplit(file_name,"[.]")[[1]][-1]=="tsv"){
            df.raw <- read.csv(file_name,row.names = 1,sep="\t")
            # print data info
            output$TO_Load <- renderPrint({
                print("A .tsv file loaded.")
                print(paste0(dim(df.raw)[2]," cells and ",dim(df.raw)[1]," genes loaded."))
                print("Print out first 5 genes in query data, in case something wrong happens in data loading: ")
                print(row.names(df.raw)[1:5])
            })
          }
          if(strsplit(file_name,"[.]")[[1]][-1]=="csv"){
            df.raw <- read.csv(file_name,row.names = 1)
            # print data info
            output$TO_Load <- renderPrint({
                print("A .csv file loaded.")
                print(paste0(dim(df.raw)[2]," cells and ",dim(df.raw)[1]," genes loaded."))
                print("Print out first 5 genes in query data, in case something wrong happens in data loading: ")
                print(row.names(df.raw)[1:5])
            })
          }
          
        }
        # standardize
        if(input$Standardize_Gene_CKB){
          # Do gene standardize
          output$TO_Standardize <- renderPrint({
              print("################# Start standardizing genes #################")
              print("This step may take minutes, please wait patiently.")
              if(strsplit(file_name,"[.]")[[1]][-1]=="rds"){
                query_data <- as.data.frame(as_matrix(dataobj@assays$RNA@data))
              }else{
                query_data <- df.raw
              }
              uniform_result <- Gene_Uniform(query_data)
              dataobj <-CreateSeuratObject(uniform_result, min.cells = 0, min.features = 0)
              dataobj[["percent.mt"]] <- PercentageFeatureSet(dataobj, pattern = "^MT-")
              Idents(dataobj) <- "data"
              app.env$Query.raw <- dataobj
          })
        }else{
          output$TO_Standardize <- renderPrint({
              print("################# Skip Gene Standardization #################")
          })
          if(strsplit(file_name,"[.]")[[1]][-1]=="rds"){
            dataobj[["percent.mt"]] <- PercentageFeatureSet(dataobj, pattern = "^MT-")
            Idents(dataobj) <- "data"
            app.env$Query.raw <- dataobj
          }else{
            dataobj <-CreateSeuratObject(df.raw, min.cells = 0, min.features = 0)
            dataobj[["percent.mt"]] <- PercentageFeatureSet(dataobj, pattern = "^MT-")
            Idents(dataobj) <- "data"
            app.env$Query.raw <- dataobj
          }
        }
      }
    )

    observeEvent(input$BTN_Demo_Dataset,{
      output$TO_Load <- renderPrint({
        output$TO_Load <- renderPrint({
          print("Demo data loaded.")
        })
        app.env$path <- "www/DemoDataobj.rds"
      })
    })
    
    observeEvent(eventExpr = app.env$Query.raw,
        handlerExpr = {
          output$PO_nFeature <- renderPlot({
              VlnPlot(app.env$Query.raw, features = "nFeature_RNA", pt.size = 0) + NoLegend() + geom_hline(aes(yintercept = input$nF_max), colour = 'blue') + geom_hline(aes(yintercept = input$nF_min), colour = 'blue')
          })
        
          output$PO_nCount <- renderPlot({
              VlnPlot(app.env$Query.raw, features = "nCount_RNA", pt.size = 0) + NoLegend() + geom_hline(aes(yintercept = input$nC_max), colour = 'blue') + geom_hline(aes(yintercept = input$nC_min), colour = 'blue')
          })
        
          output$PO_MT <- renderPlot({
              VlnPlot(app.env$Query.raw, features = "percent.mt", pt.size = 0) + NoLegend() + geom_hline(aes(yintercept = input$MT_max), colour = 'blue') + geom_hline(aes(yintercept = input$MT_min), colour = 'blue')
          })
        }
    )

    # QC
    observeEvent(eventExpr = input$BTN_MAP,
        handlerExpr = {
          if(is.null(app.env$Query.raw)){
              shinyalert(
                title = "Error",
                text = "Please upload your data before annotation.",
                size = "xs", 
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = FALSE,
                type = "error",
                showConfirmButton = TRUE,
                showCancelButton = FALSE,
                confirmButtonText = "ok",
                confirmButtonCol = "#AEDEF4",
                timer = 0,
                imageUrl = "",
                animation = TRUE
              )
          }else{
              # QC
              dataobj <- app.env$Query.raw
              dataobj <- subset(dataobj,nFeature_RNA>=input$nF_min & nFeature_RNA<=input$nF_max & nCount_RNA>=input$nC_min & nCount_RNA<=input$nC_max & percent.mt>=input$MT_min & percent.mt<=input$MT_max)
              output$TO_QCResult <- renderPrint({
                  print("Quality control finished.")
                  print(paste0(dim(dataobj)[2]," cells and ",dim(dataobj)[1]," genes after QC."))
              })
              # preprocessing
              output$TO_Preprocess1 <- renderPrint({
                print("Preprocess started.")
              })
              output$TO_Preprocess2 <- renderPrint({
                  dataobj <- NormalizeData(dataobj)
                  dataobj <- FindVariableFeatures(dataobj)
                  dataobj <- ScaleData(dataobj,features = VariableFeatures(dataobj))
                  print("Preprocess finished.")
                  dataobj <- RunPCA(dataobj, features = VariableFeatures(object = dataobj))
                  dataobj <- FindNeighbors(dataobj, reduction = "pca", dims = 1:30, nn.eps = 0.5)
                  dataobj <- RunUMAP(dataobj, dims = 1:30)
                  print("Embdedding calculated.")
                  app.env$Query.processed <- dataobj
              })
              output$PO_UMAP1 <- renderPlot({
                  dataobj <- app.env$Query.processed
                  Idents(dataobj) <- "Unannoated"
                  DimPlot(dataobj,label=T)+ NoLegend()
              })
          }
       })
    # Annotation
    observeEvent(eventExpr = input$BTN_ANNO,
        handlerExpr = {
            message(input$select_AM)
            # singleR Method
            if(input$select_AM==1){
                if(input$PA_stage=="Adult"){
                    load("www/trainedModel.Adult.Rdata")
                }else{
                    load("www/trainedModel.Fetal.Rdata")
                }
                message("Model loaded.")
                message("Annotation started.")
                dataobj <- app.env$Query.processed
                predict <- classifySingleR(dataobj@assays$RNA@data,trainedR)
                app.env$Pre.result <- predict
                message("Annotation finished.")
                dataobj$predict.result <- predict$pruned.labels
                dataobj$predict.score <- predict$tuning.scores$first
                app.env$Query.processed <- dataobj
                
                output$PO_UMAP2 <- renderPlot({
                  DimPlot(app.env$Query.processed,group.by = "predict.result",label=T)+ NoLegend()
                })
                output$PO_hist <- renderPlot({
                  hist(app.env$Query.processed$predict.score)
                })
                output$TO_Proportion <- renderPrint({
                  print(table(app.env$Query.processed$predict.result ))
                })
                output$PO_Pie <- renderPlot({
                  df <- data.frame(table(app.env$Query.processed$predict.result))
                  colnames(df) <- c("cell_type","Freq")
                  label_value <- paste0("(",signif(df$Freq/sum(df$Freq),3)*100,"%)")
                  label <-paste(sort(df$cell_type), label_value[order(df$cell_type)], sep = ' ')
                  
                  ggplot(data = df, mapping = aes(x = 'Content', y = Freq, fill = cell_type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)+ coord_polar(theta = 'y') +
                    labs(x = '', y = '', title = 'The cell-type proportion in query data') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + 
                    scale_fill_discrete(labels = label) + theme_bw() + theme(panel.grid=element_blank()) + theme(legend.title = element_text(size=15)) + theme(legend.text = element_text(size=14)) +theme(
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      panel.border = element_blank(),
                      panel.grid=element_blank(),
                      axis.ticks = element_blank(),
                      plot.title=element_text(size=15, face="bold")
                    )
                })
            }
            # Seurat Method
            if(input$select_AM==2){
                if(input$PA_stage=="Adult"){
                    dataobj.reference <- readRDS("www/Ref.Adult.rds")
                }else{
                    dataobj.reference <- readRDS("www/Ref.Fetal.rds")
                }
                message("Reference loaded.")
                message("Mapping started.")
                dataobj <- app.env$Query.processed
                anchors <- FindTransferAnchors(reference = dataobj.reference, query = dataobj,
                    dims = 1:30, reference.reduction = "pca",)
                predictions <- TransferData(anchorset = anchors, refdata = dataobj.reference$cell_type,dims = 1:30,k.weight = 10)
                message("Mapping finished.")
                dataobj$predict.result <- predictions$predicted.id
                dataobj$predict.score <- predictions$prediction.score.max
                app.env$Query.processed <- dataobj
                
                output$PO_UMAP2 <- renderPlot({
                  DimPlot(app.env$Query.processed,group.by = "predict.result",label=T)+ NoLegend()
                })
                
                output$PO_hist <- renderPlot({
                  hist(app.env$Query.processed$predict.score)
                })
                
                output$TO_Proportion <- renderPrint({
                  print("The annotated cells number of different cell-types:")
                  print(table(app.env$Query.processed$predict.result ))
                })
                
                output$PO_Pie <- renderPlot({
                  df <- data.frame(table(app.env$Query.processed$predict.result))
                  colnames(df) <- c("cell_type","Freq")
                  label_value <- paste0("(",signif(df$Freq/sum(df$Freq),3)*100,"%)")
                  label <-paste(sort(df$cell_type), label_value[order(df$cell_type)], sep = ' ')
                  
                  ggplot(data = df, mapping = aes(x = 'Content', y = Freq, fill = cell_type)) + geom_bar(stat = 'identity', position = 'stack', width = 1)+ coord_polar(theta = 'y') +
                    labs(x = '', y = '', title = 'The cell-type proportion in query data') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + 
                    scale_fill_discrete(labels = label) + theme_bw() + theme(panel.grid=element_blank()) + theme(legend.title = element_text(size=15)) + theme(legend.text = element_text(size=14)) +theme(
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      panel.border = element_blank(),
                      panel.grid=element_blank(),
                      axis.ticks = element_blank(),
                      plot.title=element_text(size=15, face="bold")
                    )
                })
            }

            # DEgene
            Idents(dataobj) <- "predict.result"
            DEGs <- FindAllMarkers(dataobj,features= VariableFeatures(dataobj),only.pos = T)
            DEGs <- DEGs[DEGs$p_val_adj <= 0.05,c(2:7)]
            app.env$DEGs <- DEGs
            output$DEG_table <- DT::renderDataTable(DT::datatable({
                DEGs
            }))
      })

      observeEvent(input$BTN_FP, {
      if(! input$select_FP %in% rownames(app.env$Query.processed)){
          shinyalert(
            title = "Error",
            text = "This gene is not found.",
            size = "xs", 
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = FALSE,
            type = "error",
            showConfirmButton = TRUE,
            showCancelButton = FALSE,
            confirmButtonText = "ok",
            confirmButtonCol = "#AEDEF4",
            timer = 0,
            imageUrl = "",
            animation = TRUE
          )
      }else{
            output$PO_FeaturePlot <- renderPlot({
                FeaturePlot(app.env$Query.processed, features=c(input$select_FP))
            })
            output$PO_VlnPlot <- renderPlot({
                options(repr.plot.width=12, repr.plot.height=6, repr.plot.res = 200) 
                VlnPlot(app.env$Query.processed, features=c(input$select_FP),group.by = "predict.result",pt.size = 0)
            })
        }
      })

    output$BTN_DOWNLOAD1 <- downloadHandler(
      filename = function() {
        "Query.processed.rds"
      },
      content = function(file) {
        saveRDS(app.env$Query.processed, file)
      }
    )

    
    output$BTN_DOWNLOAD2 <- downloadHandler(
      filename = function() {
        if(is.null(app.env$DEGs)){
          shinyalert(
            title = "Error",
            text = "Please calculate the DEGs first.",
            size = "xs", 
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = FALSE,
            type = "error",
            showConfirmButton = TRUE,
            showCancelButton = FALSE,
            confirmButtonText = "ok",
            confirmButtonCol = "#AEDEF4",
            timer = 0,
            imageUrl = "",
            animation = TRUE
          )
        }else{
          return("PredictResult.DEGs.csv")
        }
      },
      content = function(file) {
        write.csv(app.env$DEGs, file, row.names = FALSE)
      }
    )


})