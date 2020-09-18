
# Human Infertility Single_cell Testis Atlas (HISTA)


# .libPaths(c("/home/groups/monkeydo/R_LIBS/3.6.1_em"))

.libPaths()
library(BiocManager)
options(repos = BiocManager::repositories())

# InstallPkg <- function(pkg, lib="/home/groups/monkeydo/R_LIBS/3.6.1_em") {
#   if(!(pkg %in% rownames(installed.packages()))){
#     BiocManager::install(pkg , ask =F, lib=lib)
#     print(paste0(pkg, " installed"))
#   } else print(paste0(pkg, " prev. installed"))
# }
# 
# reqPkgs <- c("shiny", "rclipboard", "shinydashboard", "ggplot2", "data.table", "ggrepel", "viridis", "ggnewscale",
#              "RColorBrewer", "RColorBrewer", "grid", "gridExtra", "dplyr", "ggpubr", "BiocParallel")
# 
# for(PK in reqPkgs){
#   InstallPkg(PK)
# }

# setwd("/home/groups/monkeydo/acc/shiny/apps/HISTA")

library(shiny)
library(rclipboard)
library(shinydashboard)

# library(SDAtools)
library(shiny)
library(ggplot2)
library(data.table)
library(ggrepel)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(grid)
library(gridExtra) 
require(dplyr)
library(ggrastr)


library(ggpubr)

library("BiocParallel")
register(MulticoreParam(4))
LocalRun=T


source(system.file('app/Fxs.R', package = 'HISTA', mustWork = TRUE), local = TRUE)

if (Sys.getenv("SCRATCH_DIR") != "") {
  init.path = paste0(Sys.getenv("SCRATCH_DIR"), "/data")
  load.data.path = paste0(init.path, "/ConradLab/HISTA/ShinyServerLSV3_Sep2020.rds" )
}  else {
  if(LocalRun) init.path = "/Volumes/Maggie/Work/OHSU/Conrad/R/TestisII/HISTA_orig/data" else init.path = getwd()
  
  load.data.path = paste0(init.path, "/ShinyServerLSV3_Sep2020.rds" )
  
}

list2env(readRDS(load.data.path), envir = globalenv())

# datat <- as.data.frame(datat)
# rownames(datat) <- datat$barcode
# 
# datat$donor <- factor(datat$donor,levels=c("UtahD1", "UtahD2", "UtahD3", "AdHu173", "AdHu174", "AdHu175", "UtahI1",  "UtahI2",  "UtahK1",  "UtahK2", "Juv1", "Juv2"))
# 
# datat$donor2 <- as.character(datat$donor)
# 
# datat$donor2[which(datat$donor %in% c("UtahD1", "UtahD2", "UtahD3", "AdHu173", "AdHu174", "AdHu175"))] <-"CNT"
# datat$donor2[which(datat$donor %in% c("Juv1", "Juv2"))] <-"JUV"
# datat$donor2[which(datat$donor %in% c("UtahK1",  "UtahK2"))] <-"KS"
# datat$donor2[which(datat$donor %in% c("UtahI1"))] <-"INF1"
# datat$donor2[which(datat$donor %in% c("UtahI2"))] <-"INF2"
# datat$donor2 <- factor(datat$donor2,
#                               levels=c("CNT", "INF1", "INF2", "KS", "JUV"))
# 
# datat$experiment <- datat$donor2


# 
# names(orig)
# names(new)


col_vector <- col_vector[c(1:8, 12, 16:19, 20:26, sample(setdiff(1:length(col_vector), c(1:8, 12, 16:19, 20:26)), length(setdiff(1:length(col_vector), c(1:8, 12, 16:19, 20:26))), replace = F))]


SDA_Top100pos <- (as.data.frame(lapply(1:150, function(xN){
  as.data.frame(print_gene_list(xN, PosOnly = T, NegOnly = F))[1:150,1]
})))
colnames(SDA_Top100pos) <- paste0("SDAV", 1:150)
colnames(SDA_Top100pos) <- paste0(colnames(SDA_Top100pos), "_" , StatFac$Lab)

SDA_Top100neg <- (as.data.frame(lapply(1:150, function(xN){
  as.data.frame(print_gene_list(xN, PosOnly = F, NegOnly = T))[1:150,1]
})))
colnames(SDA_Top100neg) <- paste0("SDAV", 1:150)
colnames(SDA_Top100neg) <- paste0(colnames(SDA_Top100neg), "_" , StatFac$Lab)




ui <- dashboardPage(
  dashboardHeader(title = "HISTA"
  ),
  
  
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Germ + Soma", tabName = "combodash", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Gene expr stat. sig", tabName = "geneexprstatsig", icon = icon("affiliatetheme"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Cell type stat. sig", tabName = "celltypestatsig", icon = icon("affiliatetheme"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      # menuItem("Germ Only", tabName = "germdash", icon = icon("affiliatetheme"),
      #          badgeLabel = "soon", badgeColor = "red"),
      # menuItem("Soma Only", tabName = "somadash", icon = icon("allergies"),
      #          badgeLabel = "soon", badgeColor = "red"),
      menuItem("Score order per. Comp", tabName = "CellScoreOrderingSDA", icon = icon("arrows-alt"),
               badgeLabel = "underconst.", badgeColor = "red"),
      menuItem("Pseudotime (Germ-Only)", tabName = "pseudotimeSDA", icon = icon("arrows-alt"),
               badgeLabel = "underconst.", badgeColor = "red"),
      menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("dashboard"),
               badgeLabel = "underconst.", badgeColor = "yellow"),
      menuItem("Conrad Lab", icon = icon("file-code-o"), 
               href = "https://conradlab.org")
    )
  ),
  
  dashboardBody(
    tabItems(
      # First
      tabItem(tabName = "combodash",
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  
                  textInput("ComponentNtext", "Numerical input:", "1"),
                  textInput("Genetext", "Text input:", "PRM1"),
                  radioButtons("data", "Data origin:",
                               c("tSNE - Batch-removed SDA-CellScores" = "tsneBrSDACS",
                                 "tSNE - DropSim Norm. DGE" = "tsneDSNDGE",
                                 #"Batch-removed-Imputed-DGE" = "DGEimp",
                                 #"DromSim-Normalized-DGE" = "DGEnorm",
                                 #"tSNE - Seurat-Norm. DGE" = "tsneSNDGE",
                                 "tSNE - Batch-removed-Imputed-DGE" = "tsneImpDGE"
                               )),
                  radioButtons("metaselect", "Metadata Selection:",
                               c("Cell types" = "celltype",
                                 "Donor-replicates" = "donrep",
                                 "Donors" = "donor",
                                 "Condition" = "COND.ID",
                                 "Experiments" = "experiment"
                               )),
                  textInput("NoOfGenes", "No. of Genes to output:", "20"),
                  actionButton("C2Cpos", "Copy2ClipPosGenes"),
                  actionButton("C2Cneg", "Copy2ClipNegGenes"),
                  downloadButton("TXTall", "Copy2TxtAll"),
                  actionButton("PrevSDA", "Prev SDA"),
                  actionButton("NextSDA", "Next SDA"),
                  width = 3
                ),
                
                
                valueBoxOutput("CellType1", width = 3),
                
                valueBoxOutput("GeneName", width = 3)
                
                
              ), 
              
              box(
                #actionButton("button", "Next"),
                title = "Cell score", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot1"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Gene expr", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot4"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Pheno - tSNE", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot3"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Pheno Legend", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("plot2"),
                width = 5
              ),
              
              
              
              box(
                title = "Cell Scores Across", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("SDAScoresAcross", height = 400),
                width = 10
              ),
              
              box(
                title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("SDAScoresChiPos", height = 400), 
                width = 10, background = "black"
              ),
              box(
                title = "ChiSqrRes Scores Neg cellscores", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("SDAScoresChiNeg", height = 400), 
                width = 10, background = "black"
              ),
              
              box(
                title = "Pos. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("GOpos"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Neg. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("GOneg"),
                width = 5
              ),
              
              box(
                title = "Chrom. Location", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("ChrLoc"),
                width = 10
              ),
              
              box(title = "Pos. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                  tableOutput("packageTablePos")
              ),
              box(title = "Neg. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                  tableOutput("packageTableNeg")
              )
              
              
              
      ),
      
      # Gene expression
      tabItem(tabName = "geneexprstatsig",
              h2("Gene Expression Stat. Sig."),
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  textInput("Genetext2", "Text input:", "STAR"),
                  radioButtons("metaselect2", "Metadata Selection:",
                               c("Experiments" = "experiment"
                               ))
                ),
                box(
                  title = "Inputs 2", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  radioButtons("celltypeselect2", "Celltype Selection:",
                               c("Leydig" = "leydig",
                                 "Sertoli" = "sertoli",
                                 "Germ-All" = "germ",
                                 "Germ-UndiffSg" = "germ_UnDiffSgSct",
                                 "Germ-DiffSg" = "germ_DiffSgSct",
                                 "Germ-PrePachSct" = "germ_PrePachSct",
                                 "Germ-PachSct" = "germ_PachSct",
                                 "Germ_Std" = "germ_Std",
                                 "Endothelial" = "endothelial",
                                 "Myeloid" = "myeloid",
                                 "Adaptive" = "adaptive",
                                 "All" = "all"
                               ))
                ),
                box(
                  title = "Gene Expression Stat. Sig. Meta", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("geneexprstatsig_download"),
                  plotOutput("GeneExprSigMeta"),
                  width = 10
                )
                
                
                
              )
      ),
      
      
      
      
      # Cell type
      tabItem(tabName = "celltypestatsig",
              h2("Cell type Expression Stat. Sig."),
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  textInput("Genetext3", "Text input:", "STAR"),
                  radioButtons("metaselect3", "Metadata Selection:",
                               c("Pheno" = "FinalFinalPheno_old"
                               ))
                ),
                box(
                  title = "Gene Expression Stat. Sig. Meta", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("GeneExprSigMeta2"),
                  width = 10
                )
                
                
                
              )
      ),
      
      
      
      # # Germ
      # tabItem(tabName = "germdash",
      #         h2("Germ only")
      # ),
      # 
      # # Soma
      # tabItem(tabName = "somadash",
      #         h2("Soma only")
      # ),
      
      # CellScoreOrdering
      tabItem(tabName = "CellScoreOrderingSDA",
              h2("Cell Score Ordering"),
              fluidRow(
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  
                  textInput("ComponentNtext2", "Numerical input:", "1"),
                  # radioButtons("metaselect3", "Metadata Selection:",
                  #              c("Experiments" = "experiment"
                  #              )),
                  radioButtons("celltypeselect3", "Celltype Selection:",
                               c("Leydig" = "leydig",
                                 "Sertoli" = "sertoli",
                                 "Germ-All" = "germ",
                                 "Germ-UndiffSg" = "germ_UnDiffSgSct",
                                 "Germ-DiffSg" = "germ_DiffSgSct",
                                 "Germ-PrePachSct" = "germ_PrePachSct",
                                 "Germ-PachSct" = "germ_PachSct",
                                 "Germ_Std" = "germ_Std",
                                 "Endothelial" = "endothelial",
                                 "Myeloid" = "myeloid",
                                 "Adaptive" = "adaptive",
                                 "All" = "all"
                               ))
                ),
                box(title = "SDA comp scores ordered", status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,
                    plotOutput("CellScoreOrderSDA", height = 800),
                    width = 10
                )
                
              )
      ),
      
      # pseudotime
      tabItem(tabName = "pseudotimeSDA",
              h2("Pseudotime of Germ Cells Only"),
              fluidRow(
                box(title = "Inputs", status = "warning", solidHeader = TRUE,
                    textInput("ComponentNtext3", "Numerical input:", "1"),
                    radioButtons("metaselect4", "Metadata Selection:",
                                 c("Pseudotime" = "pseudotime",
                                   "Cell types" = "celltype",
                                   "Donor-replicates" = "donrep",
                                   "Donors" = "donor",
                                   "Condition" = "COND.ID",
                                   "Experiments" = "experiment"
                                 )),
                    width = 3),
                box(title = "Germ Cell Only tSNE", status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,
                    plotOutput("tSNEPseudoSDA", height = 400),
                ),
                box(title = "Pseudotime", status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,
                    downloadButton("PseudotimeSDA_download"),
                    plotOutput("PseudotimeSDA", height = 800),
                    width = 10
                )
                
              )
      ),
      
      # Enrichment
      tabItem(tabName = "Enrichment",
              fluidRow(
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  "Multiple formatting of gene sets accepted", 
                  br(), "List can be seperated by comma e.g. from ", 
                  br(), "   or spaces e.g. from Excel", 
                  br(), "Also, single or double quotes or not",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  textInput("GeneSet", "A set of genes", "'PRM1', 'SPATA42', 'SPRR4', 'NUPR2', 'HBZ', 'DYNLL2'"),
                  width = 8
                ),
                
                box(
                  title = "Positive Loadings", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot5"),
                  width = 10
                ),
                
                box(
                  title = "Negative Loadings", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("plot6"),
                  width = 10
                )
              )
      )
      
    )
  )
  
  
)




server <- function(input, output, session) {
  
  datat <- as.data.frame(cbind(datat, results$scores[rownames(datat),])); 
  rownames(datat) <- datat$barcode 
  
  ColFac_DONR.ID  <- as.data.frame(datat)[rownames(results$scores), ]$donor

  datat <- data.table(datat)
  
  # dl.y <- callModule(dlmodule, "input1")
  
  # SDALoadings <- results$loadings[[1]]
  
  # print(input$data)
  
  output$messageMenu <- renderMenu({
    # Code to generate each of the messageItems here, in a list. This assumes
    # that messageData is a data frame with two columns, 'from' and 'message'.
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })
    
    # This is equivalent to calling:
    #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = "messages", .list = msgs)
  })
  
  
  ################################ Reactive sections
  
  tDF <- reactive({
    # data, tsneBrSDACS, tsneDSNDGE, tsneSNDGE
    # print(input$data)
    
    
    if(input$data == "tsneBrSDACS") {
      tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC4b", "Tsne2_SDAQC4b")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
    } else {
      if(input$data == "tsneDSNDGE") {
        tempDF <- as.data.frame(datat)[, c("Tsne1", "Tsne2")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
      } else {
        if(input$data == "tsneImpDGE"){
          tempDF <- as.data.frame(datat)[, c("Tsne1_imputed", "Tsne2_imputed")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
          tempDF$tSNE1 <- as.numeric(tempDF$tSNE1)
          tempDF$tSNE2 <- as.numeric(tempDF$tSNE2)
          
        }
        
      }
    }
    
    
    
    rownames(tempDF) <- datat$barcode
    
    
    tempDF$GeneExpr <- rep(0, nrow(tempDF))
    
    return(tempDF)
  })
  
  #germ soma
  tDF_GS <- reactive({
    # data, tsneBrSDACS, tsneDSNDGE, tsneSNDGE
    # print(input$data)
    
    tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC4sperm", "Tsne2_SDAQC4sperm")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
    
    
    rownames(tempDF) <- datat$barcode
    
    # tempDF <- tempDF[!is.na(tempDF$tSNE1), ]
    
    tempDF$GeneExpr <- rep(0, nrow(tempDF))
    
    return(tempDF)
  })
  
  GEx <- reactive({
    
    if(input$celltypeselect2 == "leydig"){
      MyCells <- datat[datat$FinalFinalPheno == "Leydig",]$barcode
    } else {
      if(input$celltypeselect2 == "sertoli"){
        MyCells <- datat[datat$FinalFinalPheno == "Sertoli",]$barcode
      } else {
        if(input$celltypeselect2 == "endothelial"){
          MyCells <- datat[datat$FinalFinalPheno == "Endothelial",]$barcode
        } else {
          if(input$celltypeselect2 == "myeloid"){
            MyCells <- datat[datat$FinalFinalPheno %in% c("Macrophage-M2", "Macrophage-M1"),]$barcode
          } else {
            if(input$celltypeselect2 == "adaptive"){
              MyCells <- datat[datat$FinalFinalPheno %in% c("Lymphoid-Bcell", "Lymphoid-Tcell"),]$barcode
            } else {
              if(input$celltypeselect2 == "germ"){
                MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct", 
                                                              "Gamete_Meiotic_Pach_Dip_2nd_Scts",
                                                              "Gamete_Meiotic_preLep_Lep_Zyg_Scts",
                                                              "Gamete_RoundSpermatid",
                                                              "Gamete_UndiffSg"),]$barcode
              } else {
                if(input$celltypeselect2 == "all"){
                  MyCells <- datat$barcode
                } else {
                  if(input$celltypeselect2 == "germ_DiffSgSct"){
                    MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct"),]$barcode
                  } else {
                    if(input$celltypeselect2 == "germ_PrePachSct"){
                      MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_preLep_Lep_Zyg_Scts"),]$barcode
                    } else {
                      if(input$celltypeselect2 == "germ_PachSct"){
                        MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_Pach_Dip_2nd_Scts"),]$barcode
                      }  else {
                        if(input$celltypeselect2 == "germ_Std"){
                          MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_RoundSpermatid"),]$barcode
                        } else {
                          if(input$celltypeselect2 == "germ_UnDiffSgSct"){
                            MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_UndiffSg"),]$barcode
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
    }
    
    
    
    
    
    
    
    
    if(input$metaselect2 == "experiment") {
      MetaFac <- (datat$experiment)
    } else  if(input$metaselect2 == "FinalFinalPheno_old") {
      MetaFac <- (datat$experiment)
    } else {
      
      
    }
    
    
    
    
    if(input$Genetext2 %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneExpr <- results$scores %*% results$loadings[[1]][,as.character(input$Genetext2)]
    } else {
      GeneExpr <- results$scores %*% rep(0, nrow(results$loadings[[1]]))
      
    }
    GeneExpr <- as.data.frame(GeneExpr)
    GeneExpr$barcode <- rownames(GeneExpr)
    GeneExpr<- GeneExpr[datat$barcode, ]
    GeneExpr$MetaFac <- MetaFac
    GeneExpr <- GeneExpr[MyCells,-2]
    colnames(GeneExpr) <- c("gene", "meta")
    
    tempCom <- combn(levels(factor(GeneExpr$meta)),2)
    
    my_comparisons <- lapply(1:ncol(tempCom), function(xN){
      c(tempCom[1,xN], tempCom[2,xN])
    })
    
    return(list(GeneExpr = GeneExpr, my_comparisons = my_comparisons))
    
  })
  
  GEx2 <- reactive({
    MyCells <- datat$barcode
    
    
    
    
    
    if(input$metaselect3 == "experiment") {
      MetaFac <- (datat$experiment)
    } else  if(input$metaselect3 == "FinalFinalPheno_old") {
      MetaFac <- (datat$FinalFinalPheno_old)
    } else {
      
      
    }
    
    
    
    
    if(input$Genetext3 %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneExpr <- results$scores %*% results$loadings[[1]][,as.character(input$Genetext3)]
    } else {
      GeneExpr <- results$scores %*% rep(0, nrow(results$loadings[[1]]))
      
    }
    
    GeneExpr <- as.data.frame(GeneExpr)
    GeneExpr$barcode <- rownames(GeneExpr)
    GeneExpr<- GeneExpr[datat$barcode, ]
    GeneExpr$MetaFac <- MetaFac
    GeneExpr <- GeneExpr[MyCells,-2]
    colnames(GeneExpr) <- c("gene", "meta")
    
    tempCom <- combn(levels(factor(GeneExpr$meta)),2)
    
    my_comparisons <- lapply(1:ncol(tempCom), function(xN){
      c(tempCom[1,xN], tempCom[2,xN])
    })
    
    return(list(GeneExpr = GeneExpr, my_comparisons = my_comparisons))
    
  })
  
  GEx3 <- reactive({
    
    if(input$celltypeselect3 == "leydig"){
      MyCells <- datat[datat$FinalFinalPheno == "Leydig",]$barcode
    } else {
      if(input$celltypeselect3 == "sertoli"){
        MyCells <- datat[datat$FinalFinalPheno == "Sertoli",]$barcode
      } else {
        if(input$celltypeselect3 == "endothelial"){
          MyCells <- datat[datat$FinalFinalPheno == "Endothelial",]$barcode
        } else {
          if(input$celltypeselect3 == "myeloid"){
            MyCells <- datat[datat$FinalFinalPheno %in% c("Macrophage-M2", "Macrophage-M1"),]$barcode
          } else {
            if(input$celltypeselect3 == "adaptive"){
              MyCells <- datat[datat$FinalFinalPheno %in% c("Lymphoid-Bcell", "Lymphoid-Tcell"),]$barcode
            } else {
              if(input$celltypeselect3 == "germ"){
                MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct", 
                                                              "Gamete_Meiotic_Pach_Dip_2nd_Scts",
                                                              "Gamete_Meiotic_preLep_Lep_Zyg_Scts",
                                                              "Gamete_RoundSpermatid",
                                                              "Gamete_UndiffSg"),]$barcode
              } else {
                if(input$celltypeselect3 == "all"){
                  MyCells <- datat$barcode
                } else {
                  if(input$celltypeselect3 == "germ_DiffSgSct"){
                    MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct"),]$barcode
                  } else {
                    if(input$celltypeselect3 == "germ_PrePachSct"){
                      MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_preLep_Lep_Zyg_Scts"),]$barcode
                    } else {
                      if(input$celltypeselect3 == "germ_PachSct"){
                        MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_Pach_Dip_2nd_Scts"),]$barcode
                      }  else {
                        if(input$celltypeselect3 == "germ_Std"){
                          MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_RoundSpermatid"),]$barcode
                        } else {
                          if(input$celltypeselect3 == "germ_UnDiffSgSct"){
                            MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_UndiffSg"),]$barcode
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
    }
    
    
    
    
    
    
    
    
    # if(input$metaselect3 == "experiment") {
    #   MetaFac <- (datat$experiment)
    # } else  if(input$metaselect3 == "FinalFinalPheno_old") {
    #   MetaFac <- (datat$FinalFinalPheno_old)
    # } else {
    #   
    #   
    # }
    
    
    # print(head(results$scores))
    
    
    # GeneExpr <- as.data.frame(GeneExpr)
    # GeneExpr$barcode <- rownames(GeneExpr)
    # GeneExpr<- GeneExpr[datat$barcode, ]
    # GeneExpr$MetaFac <- MetaFac
    # GeneExpr <- GeneExpr[MyCells,-2]
    # colnames(GeneExpr) <- c("gene", "meta")
    # 
    # tempCom <- combn(levels(factor(GeneExpr$meta)),2)
    # 
    # my_comparisons <- lapply(1:ncol(tempCom), function(xN){
    #   c(tempCom[1,xN], tempCom[2,xN])
    # })
    # 
    
    
    
    return(list(Scores = results$scores[MyCells,], Meta=as.data.frame( subset(datat, barcode %in% MyCells))$experiment))
    
  })
  
  PseudotimeSDA_Rx <- reactive({
    
    Scores <- results$scores
    tempDF <- tDF_GS()
    
    PT <- datat$PseudoTime
    
    if(input$metaselect4 == "pseudotime") {
      MetaFac <- datat$PseudoTime
    } else{
      
      if(input$metaselect4 == "celltype") {
        MetaFac <- (datat$FinalFinalPheno_old)
      } else {
        if(input$metaselect4 == "donrep"){
          MetaFac <- (datat$DonRep)
        } else {
          if(input$metaselect4 == "donor"){
            MetaFac <- (datat$donor)
          } else {
            if(input$metaselect4 == "COND.ID"){
              MetaFac <- (datat$COND.ID)
            } else {
              if(input$metaselect4 == "experiment"){
                MetaFac <- (datat$experiment)
              } else {
                
              }
            }
          }
        }
      }
      
    }
    
    
    tempDF$MetFacZ <- MetaFac
    tempDF$PT <- PT
    
    tempDF <- tempDF[!is.na(tempDF$tSNE1),]
    
    tempDF$Scores <- Scores[rownames(tempDF), paste0("SDAV", as.numeric(input$ComponentNtext3))]
    tempDF$barcode <- rownames(tempDF)
    
    if(input$metaselect4 == "pseudotime") {
      tempDF$MetFacZ <- as.numeric(as.character(tempDF$MetFacZ)) 
    } else {
      tempDF$MetFacZ <- factor(as.character(tempDF$MetFacZ))
    }
    
    merge_sda_melt <- reshape2::melt(tempDF, id.vars = c("barcode","tSNE1", "tSNE2", "GeneExpr", "MetFacZ", "PT"))
    # print(head(rownames(tempDF)))
    # print(head(rownames(Scores)))
    
    # tempDF <- tempDF[!is.na(tempDF$tSNE1),]
    # Scores <-Scores[rownames(tempDF),]
    print(head(merge_sda_melt))
    # plot(merge_sda_melt$PT, 
    #      merge_sda_melt$value)
    
    
    ggpp <- ggplot(merge_sda_melt, aes(PT, value, colour=(MetFacZ))) +
      geom_point_rast(alpha=1, size=1.2, stroke=0) +
      geom_smooth(aes(PT, value), size=1, alpha = 0.6, method = "gam", formula = y ~ s(x, k = 20), se = F) +#colour="black",
      ylab("Cell Component Score") +
      xlab("Pseudotime") +
      ggtitle(paste0("SDA Comp: ", as.numeric(input$ComponentNtext3)))+
      theme_bw() +
      # scale_colour_manual(values=col_vector)+ #RColorBrewer::brewer.pal(9,"Set1")[-6]
      theme(legend.position = "none") +
      ylim(-8,8)
    
    if(input$metaselect4 == "pseudotime") {
      ggpp <-  ggpp +  scale_color_viridis()
    } else {
      ggpp <-  ggpp +  scale_colour_manual(values=col_vector)  + facet_wrap(~MetFacZ, 
                                                                            ncol=3, scales = "fixed")
    }
    ggpp
    
  })
  
  ComboTopSDAgenes <- reactive({
    Out1 <- print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>% head(as.numeric(input$NoOfGenes))
    Out1 <- Out1$Gene.Name
    
    Out2 <- print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>% head(as.numeric(input$NoOfGenes))
    Out2 <- Out2$Gene.Name
    
    data.frame(Pos=Out1, Neg=Out2)
  })
  
  GeneExprSigMeta_Rx <- reactive({
    
    GeneExpr <- GEx()$GeneExpr
    my_comparisons <- GEx()$my_comparisons
    
    
    CellType = input$celltypeselect2
    
    
    TestName = "Wilcox Rank Sum"
    
    ggboxplot(GeneExpr, x = "meta", y = "gene", palette = "jco",
              add = "jitter", col="meta") + 
      stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
      theme_bw() + 
      ggtitle( paste0(as.character(input$Genetext2), " expression :: ", TestName, " test :: ", CellType)) + 
      xlab("") + ylab(as.character(input$Genetext2))
    
    
  })
  
  ################################ observeEvent sections
  
  observeEvent(input$NextSDA, {
    Val <- as.character(min(c(150, as.numeric(input$ComponentNtext)+1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$PrevSDA, {
    Val <- as.character(max(c(1, as.numeric(input$ComponentNtext)-1)))
    
    updateTextInput(session, "ComponentNtext", value = Val)
  })
  
  observeEvent(input$C2Cpos, {
    
    
    Out1 <- print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out1 <- Out1$Gene.Name
    
    # print(Out1)
    clipr::write_clip(Out1)
    
  })

  observeEvent(input$C2Cneg, {
    
    
    Out2 <- print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes)) 
    Out2 <- Out2$Gene.Name
    
    # print(Out1)
    clipr::write_clip(Out2)
    
  })
  
  ################################ renderPlot sections
  
  output$CellScoreOrderSDA <- renderPlot({
    
    Scores <- GEx3()$Scores
    Meta <- GEx3()$Meta
    
    # print(head(Scores[,as.numeric(input$ComponentNtext2)]))
    
    SC_SDA_ScoreDF <- data.frame(Score=Scores[,as.numeric(input$ComponentNtext2)],
                                 Rank=rank(Scores[,as.numeric(input$ComponentNtext2)]))
    # print(head(SC_SDA_ScoreDF))
    # print(head(Meta))
    
    SC_SDA_ScoreDF$Meta <- Meta #datat$experiment2[rownames(Scores)]
    
    
    cowplot::plot_grid(ggplot(SC_SDA_ScoreDF ,
                              aes(x = Score, y = Rank, colour = Meta)) +
                         geom_point(alpha=.2) +
                         ggtitle(paste0("Cells ordered by SDA Scores in SDA", input$ComponentNtext2)) +
                         xlab(paste0("Score : SDA", input$ComponentNtext2)) +
                         ylab(paste0("Rank : SDA", input$ComponentNtext2)) +
                         theme_bw() + scale_color_manual(values=col_vector) +
                         theme(legend.position="below",
                               legend.direction="horizontal",
                               legend.title = element_blank(),
                               axis.text.x = element_text(angle = 90)),
                       ggplot(SC_SDA_ScoreDF, 
                              aes(x = Score, 
                                  y = Meta, 
                                  colour = Meta)) +
                         ggbeeswarm::geom_quasirandom(groupOnX = FALSE)  +
                         xlab(paste0("Score : SDA", input$ComponentNtext2)) + 
                         ylab("Categ.")  +
                         theme_bw() + scale_color_manual(values=col_vector) +
                         theme(legend.position="bottom",
                               legend.direction="horizontal",
                               legend.title = element_blank(),
                               axis.text.x = element_text(angle = 90)),
                       ncol=1)
    
    
    
    
    
  })
  
  output$tSNEPseudoSDA <- renderPlot({
    
    
    tempDF <- tDF_GS()
    
    if(input$metaselect4 == "pseudotime") {
      MetaFac <- datat$PseudoTime
    } else{
      
      if(input$metaselect4 == "celltype") {
        MetaFac <- (datat$FinalFinalPheno_old)
      } else {
        if(input$metaselect4 == "donrep"){
          MetaFac <- (datat$DonRep)
        } else {
          if(input$metaselect4 == "donor"){
            MetaFac <- (datat$donor)
          } else {
            if(input$metaselect4 == "COND.ID"){
              MetaFac <- (datat$COND.ID)
            } else {
              if(input$metaselect4 == "experiment"){
                MetaFac <- (datat$experiment)
              } else {
                
              }
            }
          }
        }
      }
      
    }
    
    
    tempDF$MetFacZ <- MetaFac
    
    tempDF <- tempDF[!is.na(tempDF$tSNE1),]
    
    
    if(input$metaselect4 == "pseudotime") {
      
      ggplot(tempDF, aes(tSNE1, tSNE2, color=MetFacZ)) +
        geom_point(size=0.1) + theme_bw() +
        scale_color_viridis() +
        theme(legend.position = "bottom", aspect.ratio=1,
              legend.title = element_blank())  +
        ggtitle("SDA t-SNE - cells coloured by PseudoTime")  +
        coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
      
    } else {
      
      #ggplotly
      ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(as.character(MetFacZ)))) +
        geom_point(size=0.1)+ theme_bw() +
        theme(legend.position = "bottom", aspect.ratio=1,
              legend.title = element_blank()) +
        ggtitle("Germ-cell Only t-SNE") +
        scale_color_manual(values=(col_vector)) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1), nrow =3))  +
        coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
      
    }
    
    
    
    
    
    
    
  })
  
  output$PseudotimeSDA <- renderPlot({
    
    PseudotimeSDA_Rx()
    
  })

  output$SDAScoresChiPos <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    
    SDAScores <- results$scores
    ComponentN <- 1:ncol(SDAScores)
    MetaDF <- as.data.frame(datat)
    rownames(MetaDF) <- datat$barcode
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    
    
    PosCompsDF <- as.data.frame(lapply(levels(factor(MetaDF$experiment)), function(CondX){
      
      apply(SDAScores[rownames(MetaDF)[which(MetaDF$experiment == CondX)], ], 2, 
            function(x){
              round(sum(x>0)/nrow(SDAScores)*100, 2)
            })
    }))
    
    colnames(PosCompsDF) <- levels(factor(MetaDF$experiment))
    
    # print(head(PosCompsDF))
    # print(min(PosCompsDF))
    
    
    PosCompsDF <- PosCompsDF[gtools::mixedsort(rownames(PosCompsDF)),]
    PosCompsDF$SDA <- factor(rownames(PosCompsDF), levels=rownames(PosCompsDF))
    
    print(head(PosCompsDF))
    
    ChiT <- chisq.test(PosCompsDF[,1:(ncol(PosCompsDF)-1)])
    
    ChiTres <- ChiT$residuals
    ChiTres[which(is.na(ChiTres))] = 0
    
    ChiResSD <- round(apply(ChiTres, 1, sd),2)
    ChiResSD[which(is.na(ChiResSD))] <- 0
    ChiResSD[ChiResSD < 0.2] <- ""
    
    # if(is.null(envv$SDAScoresChi_clusBTN)) {
    #   clustStat = F
    # } else {
    #   clustStat <- ifelse(envv$SDAScoresChi_clusBTN=="ON", T, F)
    # }
    
    clustStat = T
    
    pheatmap::pheatmap((t(ChiT$residuals)),
                       cluster_cols = clustStat, cluster_rows = clustStat,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                       labels_col = paste0(rownames(PosCompsDF), " sd_", ChiResSD)
    )
    
    
    
    
    
    
    
    
    
    
    
  })
  
  output$SDAScoresChiNeg <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    
    SDAScores <- results$scores
    ComponentN <- 1:ncol(SDAScores)
    MetaDF <- as.data.frame(datat)
    rownames(MetaDF) <- datat$barcode
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    
    
    NegCompsDF <- as.data.frame(lapply(levels(factor(MetaDF$experiment)), function(CondX){
      
      apply(SDAScores[rownames(MetaDF)[which(MetaDF$experiment == CondX)], ], 2, 
            function(x){
              round(sum(x<0)/nrow(SDAScores)*100, 2)
            })
    }))
    
    colnames(NegCompsDF) <- levels(factor(MetaDF$experiment))
    
    # print(head(NegCompsDF))
    # print(min(NegCompsDF))
    
    
    NegCompsDF <- NegCompsDF[gtools::mixedsort(rownames(NegCompsDF)),]
    NegCompsDF$SDA <- factor(rownames(NegCompsDF), levels=rownames(NegCompsDF))
    
    print(head(NegCompsDF))
    
    ChiT <- chisq.test(NegCompsDF[,1:(ncol(NegCompsDF)-1)])
    
    ChiTres <- ChiT$residuals
    ChiTres[which(is.na(ChiTres))] = 0
    
    ChiResSD <- round(apply(ChiTres, 1, sd),2)
    ChiResSD[which(is.na(ChiResSD))] <- 0
    ChiResSD[ChiResSD < 0.2] <- ""
    
    # if(is.null(envv$SDAScoresChi_clusBTN)) {
    #   clustStat = F
    # } else {
    #   clustStat <- ifelse(envv$SDAScoresChi_clusBTN=="ON", T, F)
    # }
    
    clustStat = T
    
    pheatmap::pheatmap((t(ChiT$residuals)),
                       cluster_cols = clustStat, cluster_rows = clustStat,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                       labels_col = paste0(rownames(NegCompsDF), " sd_", ChiResSD)
    )
    
    
    
    
    
    
    
    
    
    
    
  })
  
  output$GeneExprSigMeta <- renderPlot({
    (GeneExprSigMeta_Rx())
  })

  output$GeneExprSigMeta2 <- renderPlot({
    
    GeneExpr <- GEx2()$GeneExpr
    my_comparisons <- GEx2()$my_comparisons
    
    # print(head(GeneExpr))
    # print(head(my_comparisons))
    
    # CellType = input$celltypeselect
    
    
    
    GeneExpr$meta <- factor(GeneExpr$meta)
    GeneExpr$meta <- factor(GeneExpr$meta, levels = gtools::mixedsort(levels(GeneExpr$meta)) )
    
    
    TestName = "Wilcox Rank Sum"
    
    
    
    ggboxplot(GeneExpr, x = "meta", y = "gene", palette = col_vector,
              add = "jitter", col="meta") +
      #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
      theme_bw() +
      ggtitle( paste0(as.character(input$Genetext3), " expression :: "
                      #, TestName, " test "
      )) +
      xlab("") + ylab(as.character(input$Genetext3))  +
      theme(legend.position="bottom",
            legend.direction="horizontal",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    
    
  })
  
  output$packageTablePos <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), PosOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    print_gene_list(as.numeric(input$ComponentNtext), NegOnly = T) %>%
      #group_by(package) %>%
      #tally() %>%
      #arrange(desc(n), tolower(package)) %>%
      #mutate(percentage = n / nrow(pkgData()) * 100) %>%
      #select("Package name" = package, "% of downloads" = percentage) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)

  #renderPlotly
  output$plot1 <- renderPlot({
    #ggplotly
    tempDF <- tDF()
    (ggplot(cbind(tempDF, SDAComp=datat[,get(paste0("SDAV", input$ComponentNtext, sep=""))]), 
            aes(tSNE1, tSNE2, color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
        geom_point(size=0.1) + theme_bw() +
        scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        theme(legend.position = "bottom", aspect.ratio=1) + 
        ggtitle(paste0("SDAV", input$ComponentNtext, " \n", 
                       StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], sep="")) + 
        simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE))
    
  })
  
  output$plot2 <- renderPlot({
    tempDF <- tDF()
    
    if(input$metaselect == "celltype") {
      MetaFac <- (datat$FinalFinalPheno_old)
    } else {
      if(input$metaselect == "donrep"){
        MetaFac <- (datat$DonRep)
      } else {
        if(input$metaselect == "donor"){
          MetaFac <- (datat$donor)
        } else {
          if(input$metaselect == "COND.ID"){
            MetaFac <- (datat$COND.ID)
          } else {
            if(input$metaselect == "experiment"){
              MetaFac <- (datat$experiment)
            } else {
              
            }
          }
        }
      }
    }
    
    
    
    ggFph <- ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(MetaFac))) +
      geom_point(size=0.1)+ theme_bw() +
      theme(legend.position = "bottom", aspect.ratio=1,
            legend.title = element_blank()) +
      ggtitle("t-SNE - Final Pheno") +
      scale_color_manual(values=(col_vector)) + guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
      coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
    legend <- cowplot::get_legend(ggFph)
    
    #grid.newpage()
    grid.draw(legend)
    
  })

  #renderPlotly
  output$plot3 <- renderPlot({
    tempDF <- tDF()
    
    if(input$metaselect == "celltype") {
      MetaFac <- (datat$FinalFinalPheno_old)
    } else {
      if(input$metaselect == "donrep"){
        MetaFac <- (datat$DonRep)
      } else {
        if(input$metaselect == "donor"){
          MetaFac <- (datat$donor)
        } else {
          if(input$metaselect == "COND.ID"){
            MetaFac <- (datat$COND.ID)
          } else {
            if(input$metaselect == "experiment"){
              MetaFac <- (datat$experiment)
            } else {
              
            }
          }
        }
      }
    }
    
    #ggplotly
    (ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(MetaFac))) +
        geom_point(size=0.1)+ theme_bw() +
        theme(legend.position = "none", aspect.ratio=1) +
        ggtitle("t-SNE - Meta") +
        scale_color_manual(values=(col_vector))   +
        coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE))
    
    
  })
  
  output$plot4 <- renderPlot({
    
    tempDF <- as.data.frame(tDF())
    
    if(input$Genetext %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneExpr <- results$scores %*% results$loadings[[1]][,as.character(input$Genetext)]
    } else {
      GeneExpr <- results$scores %*% rep(0, nrow(results$loadings[[1]]))
      
    }
    #ggplotly
    #as.numeric(input$NoOfGenes)
    
    LoadOrdVal <- round(results$loadings[[1]][,as.character(input$Genetext)][order(abs(results$loadings[[1]][,as.character(input$Genetext)]), decreasing = T)], 3)
    
    
    tempDF[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]
    
    (ggplot(tempDF, 
            aes(tSNE1, tSNE2, color=cut(asinh(GeneExpr), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
        geom_point(size=0.1) + theme_bw() +
        scale_color_manual("EX", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        theme(legend.position = "bottom", aspect.ratio=1) + 
        
        simplify2 + 
        coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)) + 
      labs(title = paste("Gene: ", input$Genetext, sep=""), 
           subtitle = paste("Found in comps: \n",
                            paste(names(LoadOrdVal)[1:5], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[1:5], collapse = ", "), 
                            "\n",
                            paste(names(LoadOrdVal)[6:10], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[6:10], collapse = ", "), 
                            "\n"), 
           caption = "Caption here")
    
    
    
    
    
  })
  
  output$GOpos <- renderPlot({
    
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
      print("No GO")
    } else {
      go_volcano_plot(component = paste("V", input$ComponentNtext, "P", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  output$GOneg <- renderPlot({
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
      print("No GO")
    } else {
      go_volcano_plot(component = paste("V", input$ComponentNtext, "N", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  output$ChrLoc <- renderPlot({
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
      print("No Comp")
    } else {
      pgl <- genome_loadings(results$loadings[[1]][as.numeric(input$ComponentNtext),], 
                             label_both = T, 
                             max.items = 10, 
                             gene_locations =   gene_locations,
                             chromosome_lengths = chromosome.lengths)+ theme(aspect.ratio = .5)
      print(pgl)
      
    }
    
  })
  
  output$SDAScoresAcross <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
      print("No Comp")
    } else {
      SDAScores <- results$scores
      ComponentN <- as.numeric(input$ComponentNtext)
      
      pgl <- ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                               score = SDAScores[, paste0("SDAV", ComponentN)], 
                               experiment = gsub("_.*", "", gsub("[A-Z]+\\.", "", rownames(SDAScores))), 
                               ColFac = ColFac_DONR.ID), 
                    aes(cell_index, score, colour = ColFac)) + 
        geom_point(size = 0.5, stroke = 0) + 
        xlab("Cell Index") + ylab("Score") + 
        #scale_color_brewer(palette = "Paired") + 
        theme_bw() + 
        theme(legend.position = "bottom") + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        scale_colour_manual(values =(col_vector),
                            guide = guide_legend(nrow=2)) +
        #guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + 
        ggtitle(paste0("SDAV", ComponentN))
      
      print(pgl)
      
    }
    
  })
  
  output$plot5 <- renderPlot({
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=8025
    # k = number of genes submitted, top 100
    k = 150 #100
    GeneSet <- input$GeneSet
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      #print(GeneSet)
    }else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      #print(GeneSet)
    }
    
    GeneSetNot <- GeneSet[!GeneSet %in% colnames(results$loadings[[1]][,])]
    
    print("length of your genes:")
    print(length(GeneSet))
    GeneSet <- GeneSet[GeneSet %in% colnames(results$loadings[[1]][,])]
    print("length of your genes in this dataset:")
    print(length(GeneSet))
    
    
    
    
    # print("length of your genes in this dataset:")
    # print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_Top100neg, 
               GeneVec = GeneSet, 
               plotTitle= paste0("Gene-set enrichment\n SDA top 150 pos loadings\n Cust. Input. genes \n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                                 paste0(GeneSetNot, collapse = ", ")),
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  output$plot6 <- renderPlot({
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=8025
    # k = number of genes submitted, top 100
    k = 150 #100
    GeneSet <- input$GeneSet
    #GeneSet <- "'PRM1', 'SPATA42', 'SPRR4', 'NUPR2', 'HBZ', 'DYNLL2'"
    
    
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      #print(GeneSet)
    }else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      #print(GeneSet)
    }
    
    # print("length of your genes:")
    # print(length(GeneSet))
    GeneSetNot <- GeneSet[!GeneSet %in% colnames(results$loadings[[1]][,])]
    
    
    
    
    GeneSet <- GeneSet[GeneSet %in% colnames(results$loadings[[1]][,])]
    # print("length of your genes in this dataset:")
    # print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_Top100neg, 
               GeneVec = GeneSet, 
               plotTitle= paste0("Gene-set enrichment\n SDA top 150 neg loadings\n Cust. Input. genes \n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                                 paste0(GeneSetNot, collapse = ", ")),
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  ################################ renderValueBox sections
  
  output$CellType1 <- renderValueBox({
    valueBox(
      value = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),6],
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  
  output$GeneName <- renderValueBox({
    if(input$Genetext %in% colnames(results$loadings[[1]])){
      # results$loadings[[1]][,"PRM1"]
      GeneName <- input$Genetext
    } else {
      GeneName <- paste0(input$Genetext, " Not Found")
      
    }
    
    valueBox(
      value = GeneName, #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = "Gene Name",
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  
  ################################ downloadHandler sections
  
  output$TXTall <- downloadHandler(
    filename = function(){
      paste("data-TopGenes_SDA_negNpos", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # print(ComboTopSDAgenes())
      write.csv(ComboTopSDAgenes(), file, row.names=F)
      # write.table(paste(text,collapse=", "), file,col.names=FALSE)
    }
  )
  
  output$geneexprstatsig_download <- downloadHandler(
    filename = function(){
      paste("geneexprstatsig_download", Sys.Date(), ".pdf", sep = "")
      # "test.pdf"
    },
    content = function(file) {
      pdf(file, width = 12, height =9, compress = T, pointsize = 15)
      # grid.text("This is some initial text",  x=0.5, y=.9, gp=gpar(fontsize=18), check=TRUE)
      # grid::grid.newpage()
      plot(GeneExprSigMeta_Rx())
      # grid::grid.newpage()
      # grid.text("This is some final text",  x=0.5, y=.9, gp=gpar(fontsize=18), check=TRUE)
      # ggsave(file,GeneExprSigMeta_Rx(), width = 12, height =9,  units="in", device = "pdf")
      
      dev.off()
    })
  
  output$PseudotimeSDA_download <- downloadHandler(
    filename = function(){
      paste("PseudotimeSDA_download_15x5", Sys.Date(), ".pdf", sep = "")
      # "test.pdf"
    },
    content = function(file) {
      pdf(file, width = 15, height =5, compress = T, pointsize = 15)
      plot(PseudotimeSDA_Rx())
      dev.off()
    })
  
}


shinyApp(ui, server)

