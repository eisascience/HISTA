
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
  load.data.path = paste0(init.path, "/ConradLab/HISTA/ShinyServerLSV3_Sep2020B.rds" )
}  else {
  if(LocalRun) init.path = "/Volumes/Maggie/Work/OHSU/Conrad/R/TestisII/HISTA_orig/data" else init.path = getwd()
  
  load.data.path = paste0(init.path, "/ShinyServerLSV3_Sep2020B.rds" )
  
}

list2env(readRDS(load.data.path), envir = globalenv())

print(table(datat$FinalFinalPheno))
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
      menuItem("tSNE-SDA score per Celltype", tabName = "tsnepercelltype", icon = icon("arrows-alt"),
               badgeLabel = "soon", badgeColor = "yellow"),
      menuItem("tSNE-Meta per Celltype", tabName = "tsnepercelltype_meta", icon = icon("arrows-alt"),
               badgeLabel = "soon", badgeColor = "yellow"),
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
                downloadButton("tSNEwSDAScoreProj_download"),
                plotOutput("tSNEwSDAScoreProj"), #plotlyOutput
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
                downloadButton("tSNEwMeta_download"),
                plotOutput("tSNEwMeta"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Pheno Legend", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("tSNEwMetaLegend"),
                width = 5
              ),
              
              
              
              box(
                downloadButton("SDAScoresAcross_download"),
                title = "Cell Scores Across", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                plotOutput("SDAScoresAcross", height = 400),
                width = 10
              ),
              
              box(
                title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("SDAScoresChiPos_download"),
                plotOutput("SDAScoresChiPos", height = 400), 
                width = 10, background = "black"
              ),
              box(
                title = "ChiSqrRes Scores Neg cellscores", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("SDAScoresChiNeg_download"),
                plotOutput("SDAScoresChiNeg", height = 400), 
                width = 10, background = "black"
              ),
              
              box(
                title = "Pos. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("SDAGOpos_download"),
                plotOutput("GOpos"), #plotlyOutput
                width = 5
              ),
              box(
                title = "Neg. Loadings GO", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("SDAGOneg_download"),
                plotOutput("GOneg"),
                width = 5
              ),
              
              box(
                title = "Chrom. Location", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("ChrLoc_download"),
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
                                 "Myoid" = "myoid",
                                 "Neuro" = "neuro",
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
      
      
      
     
      tabItem(tabName = "tsnepercelltype",
              h2("tSNE with SDA score by cell type"),
              fluidRow(
                box(
                  title = "Inputs 2", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  textInput("ComponentNtext_tsnepercelltype", "Numerical input:", "1"),
                  radioButtons("tsnepercelltype_ctselect", "Celltype Selection:",
                               c("Leydig" = "leydig",
                                 "Sertoli" = "sertoli",
                                 "Myoid" = "myoid",
                                 "Neuro" = "neuro",
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
                               ), selected = "all")
                ),
                box(
                  title = "tSNE with SDA score Projected per CellType", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnepercelltype_download"),
                  plotOutput("tSNEperCellType"),
                  width = 10
                )
                
                
                
              )
      ),
      
      
      
      
      tabItem(tabName = "tsnepercelltype_meta",
              h2("tSNE with Metada by cell type"),
              fluidRow(
                box(
                  title = "Inputs 2", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  radioButtons("tsnepercelltype_ctselect_meta", "Celltype Selection:",
                               c("Leydig" = "leydig",
                                 "Sertoli" = "sertoli",
                                 "Myoid" = "myoid",
                                 "Neuro" = "neuro",
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
                               ), selected = "all"),
                  radioButtons("ctselect_meta", "Metadata Selection:",
                               c("Cell types" = "celltype",
                                 "Donor-replicates" = "donrep",
                                 "Donors" = "donor",
                                 "Condition" = "COND.ID",
                                 "Experiments" = "experiment"
                               )),
                ),
                box(
                  title = "tSNE with SDA score Projected per CellType", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnepercelltype_meta_download"),
                  plotOutput("tSNEperCellType_meta"),
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
                                 "Myoid" = "myoid",
                                 "Neuro" = "neuro",
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
  
  AddPer <- function(x, perc=0.1){
    x + x *  perc
  }
  
  
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
  source("app_Reactive.R",local = TRUE)
  
################################ observeEvent sections
  source("app_ObserveEvents.R",local = TRUE)
  
################################ renderPlot sections
  source("app_RenderPlots.R",local = TRUE)
  
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
  source("app_DownloadButtons.R",local = TRUE)
  
  
}


shinyApp(ui, server)

