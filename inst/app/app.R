
# Human Infertility Single_cell Testis Atlas (HISTA)



library(shiny)
library(rclipboard)
library(shinydashboard)


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

library(DT)


pathi = getwd()


source("app_Fxs.R",local = TRUE)


# source(paste0(pathi, "inst/app/app_Fxs.R"))

list2env(readRDS( paste0(pathi, "/data/HISTAv1_dataLS_oct2021.rds")), envir = globalenv()) #for deploy on shinyapps



col_vector <- col_vector[c(1:8, 12, 16:19, 20:26, sample(setdiff(1:length(col_vector), c(1:8, 12, 16:19, 20:26)), length(setdiff(1:length(col_vector), c(1:8, 12, 16:19, 20:26))), replace = F))]


# SDA_Top100pos <- (as.data.frame(lapply(1:150, function(xN){
#   as.data.frame(print_gene_list(xN, PosOnly = T, NegOnly = F))[1:150,1]
# })))
# colnames(SDA_Top100pos) <- paste0("SDAV", 1:150)
# colnames(SDA_Top100pos) <- paste0(colnames(SDA_Top100pos), "_" , StatFac$Lab)
# 
# SDA_Top100neg <- (as.data.frame(lapply(1:150, function(xN){
#   as.data.frame(print_gene_list(xN, PosOnly = F, NegOnly = T))[1:150,1]
# })))
# colnames(SDA_Top100neg) <- paste0("SDAV", 1:150)
# colnames(SDA_Top100neg) <- paste0(colnames(SDA_Top100neg), "_" , StatFac$Lab)


## ui ------

ui <- dashboardPage(
  dashboardHeader(title = "HISTA v2.8"
  ),
  
  ## dashboard items ------
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home Page", tabName = "homepage", icon = icon("home"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Main tab", tabName = "combodash", icon = icon("star"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Index of Components", tabName = "sdaindex", icon = icon("tree"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per pathology (boxplot)", tabName = "geneexprstatsig", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per cell type (boxplot)", tabName = "celltypestatsig", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per cell type (tSNE)", tabName = "celltypeGeneExpr2D", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Fingerprinting (heatmap)", tabName = "chisqrSDAenrichHeatmaps", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Cell score per celltype (tSNE)", tabName = "tsnepercelltype", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Metadata per celltype (tSNE)", tabName = "tsnepercelltype_meta", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      #menuItem("Score order per. Comp", tabName = "CellScoreOrderingSDA", icon = icon("plus"),
               #badgeLabel = "final", badgeColor = "green"),
      menuItem("Pseudotime (germ only)", tabName = "pseudotimeSDA", icon = icon("random"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Pseudotime Component Index", tabName = "pseudotimeSDAIndex", icon = icon("tree"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("asterisk"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Soma only W. LN19", tabName = "somaWLN", icon = icon("minus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("LC only W. Zhao21", tabName = "LConly", icon = icon("minus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("User Manual", tabName = "usermanual", icon = icon("book"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Conrad Lab", icon = icon("file-code-o"), 
               href = "https://conradlab.org"),
      menuItem("@eisamahyari", icon = icon("heart"), 
               href = "https://eisascience.github.io")
    )
  ),
  
  # dashboard body ------
  
  dashboardBody(
    tabItems(
      
      ## template tabs -----
      
      # tabItem(tabName = "germdash",
      #         h2("Germ only"),
      #         fluidRow(
      #           
      #         )
      # ),

      
      
      ## homepage -----
      tabItem(tabName = "homepage",
              h2("Welcome to HISTA's Homepage"),
              fluidRow(
                
                imageOutput("homepage", width = 900)
                
                
              )
      ),
      
      
 
     
     
     
      ## sdaindex -----

     tabItem(tabName = "sdaindex",
             h2("Index of SDA components and annotation"),
             fluidRow(
              
               
               DT::dataTableOutput("SDAannotations")
               
             )
     ),

     
     ## pseudotimeSDAIndex -----
     
     tabItem(tabName = "pseudotimeSDAIndex",
             h2("Index of SDA components and their pseudo-order"),
             fluidRow(
               DT::dataTableOutput("SDApseudotime")

             )
     ),
     
      
     ## combo main dashboard ------
      
      # First
      tabItem(tabName = "combodash",
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  
                  textInput("ComponentNtext", "SDA component search (numerical):", "1"),
                  textInput("Genetext", "Gene search (text input):", "PRM1"),
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
                  actionButton("C2Cpos", "Copy2ClipPosGenes (local only)"),
                  actionButton("C2Cneg", "Copy2ClipNegGenes (local only)"),
                  downloadButton("TXTall", "Copy2TxtAll (any browser"),
                  actionButton("PrevSDA", "Prev SDA"),
                  actionButton("NextSDA", "Next SDA"),
                  #actionButton("ScreenShot", "Take a screenshot"),
                  width = 3
                ),
                
                
                valueBoxOutput("CellType1", width = 3),
                
                valueBoxOutput("GeneName", width = 3)
                
                
              ), 
              box(
                title = "Gene expr", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,
                downloadButton("tSNEwGeneExpr_download"),
                plotOutput("tSNE_geneExpr"), #plotlyOutput
                width = 5
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
              
              # box(
              #   title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
              #   collapsible = TRUE,
              #   downloadButton("SDAScoresChiPos_download"),
              #   plotOutput("SDAScoresChiPos", height = 400), 
              #   width = 10, background = "black"
              # ),
              # box(
              #   title = "ChiSqrRes Scores Neg cellscores", status = "primary", solidHeader = TRUE,
              #   collapsible = TRUE,
              #   downloadButton("SDAScoresChiNeg_download"),
              #   plotOutput("SDAScoresChiNeg", height = 400), 
              #   width = 10, background = "black"
              # ),
              
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
      
      
     ##Soma only w LN19 -----
     
      tabItem(tabName = "somaWLN",
              h2("Reprocessing of: \nMahyari-Conrad et al. (2021) + \n Laurentino-Neuhaus (2019)*"),
              fluidRow(
                box(
                  title = "Phenotype", status = "primary", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnesomaonlywln_phenotype_download"),
                  plotOutput("tSNE_somaWLN_Pheno3_Rx"),
                  width = 5,
                  # footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/"
                ),
                box(
                  title = "Condition", status = "primary", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnesomaonlywln_condition_download"),
                  plotOutput("tSNE_somaWLN_COND.ID_Rx"),
                  width = 5,
                  # footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/"
                ),
                box(
                  title = "Donors", status = "primary", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnesomaonlywln_donor_download"),
                  plotOutput("tSNE_somaWLN_DONR.ID_Rx"),
                  width = 5,
                  footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/"
                ),
                box(
                  title = "Total No. of transcripts", status = "primary", 
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  downloadButton("tsnesomaonlywln_ncount_download"),
                  plotOutput("tSNE_somaWLN_nCount_RNA_Rx"),
                  width = 5,
                  # footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/"
                )
                
              )
        
      ),
     
      
     ##LC only w Zhao20 LN19 -----
     
     tabItem(tabName = "LConly",
             h2("Reprocessing of Leydig Cells (LC): \nMahyari-Conrad et al. (2021) + \n Laurentino-Neuhaus (2019) + \n Zhao-Li (2020)**"),
             fluidRow(
               
               box(
                 title = "Mahyari et al. only data", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConly_donors_Rx"),
                 width = 5,
                 #footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               box(
                 title = "Mahyari et al. only data", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConly_phenotype_Rx"),
                 width = 5,
                 #footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               box(
                 title = "LC from Mahyari21 + LN19 + Zhao20", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConlyZhao_donors_Rx"),
                 width = 5,
                 footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               box(
                 title = "LC from Mahyari21 + LN19 + Zhao20", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConlyZhao_phenotype_Rx"),
                 width = 5,
                 #footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               box(
                 title = "LC from Mahyari21 + LN19 + Zhao20", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConlyZhao_phenotypeProp_Rx"),
                 width = 5,
                 #footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               box(
                 title = "LC from Mahyari21 + LN19 + Zhao20", status = "primary", 
                 solidHeader = TRUE,
                 collapsible = TRUE,
                 # downloadButton("tsnesomaonlywln_donor_download"),
                 plotOutput("DimRedux_LConlyZhao_KeyGenesViolin_Rx"),
                 width = 10,
                 #footer = "* Laurentino et al. 2019 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6714305/ \n** Zhao et al. 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7655944/"
               ),
               
               
             )
             
     ),
      
     
     ##Gene expression -----

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
     
     ##tSNE per cell type gene expr-----
     
     
     tabItem(tabName = "celltypeGeneExpr2D",
             h2("tSNE with gene expr by cell type"),
             fluidRow(
               box(
                 title = "Inputs 2", status = "warning", solidHeader = TRUE,
                 #"Box content here", br(), "More box content",
                 #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                 width = 5,
                 textInput("Genetext_celltype_gex", "Gene search (text input):", "PRM1"),
                 radioButtons("tsnepercelltype_ctselect_gex", "Celltype Selection:",
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
                 title = "tSNE with gene expr Projected per CellType", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 downloadButton("tsnepercelltype_gex_download"),
                 plotOutput("tSNEperCellType_GEX"),
                 width = 10
               )
               
               
               
             )
     ),
      
     ##ChisqrHeatmaps -----
     
     tabItem(tabName = "chisqrSDAenrichHeatmaps",
             h2("SDA component pathology fingerprinting"),
             fluidRow(
               
               
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
               
               
             )
     ),
      
     ##tSNE per cell type cell score -----
     
     
      tabItem(tabName = "tsnepercelltype",
              h2("tSNE with SDA score by cell type"),
              fluidRow(
                box(
                  title = "Inputs 2", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  textInput("ComponentNtext_tsnepercelltype", "SDA component search (numerical):", "1"),
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
      
      
     ##tSNE per meta -----
     
      
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
      
      
      
     ##Cell type stat sig -----
     
      tabItem(tabName = "celltypestatsig",
              h2("Cell type Expression Stat. Sig."),
              fluidRow(
                
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  width = 5,
                  textInput("Genetext3", "Text input:", "STAR"),
                  radioButtons("metaselect3", "Phenotype scope:",
                               c("Pheno" = "scope1"
                               )),
                  radioButtons("condSelect1", "Conditions:",
                               c("All donors" = "all",
                                 "Adult controls" = "cnt",
                                 "Inf 1" = "inf1",
                                 "Inf 2" = "inf2",
                                 "KS" = "ks", 
                                 "JUV"= "juv"
                               ), selected = "cnt")
                ),
                box(
                  title = "Gene Expression Stat. Sig. Meta", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("GeneExprSigMeta2"),
                  width = 10
                )
                
                
                
              )
      ),
      
      
      
      
      
     ##Cell score ordering -----
     
      tabItem(tabName = "CellScoreOrderingSDA",
              h2("Cell Score Ordering"),
              fluidRow(
                box(
                  title = "Inputs", status = "warning", solidHeader = TRUE,
                  #"Box content here", br(), "More box content",
                  #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                  
                  textInput("ComponentNtext2", "SDA component search (numerical):", "1"),
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
      
     ##Pseudotime -----
     
      tabItem(tabName = "pseudotimeSDA",
              h2("Pseudotime of Germ Cells Only"),
              fluidRow(
                box(title = "Inputs", status = "warning", solidHeader = TRUE,
                    textInput("ComponentNtext3", "SDA component search (numerical):", "1"),
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
      
     ## Enrichment -----
     
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
  
  
################################ Reactive sections-----
  source("app_Reactive.R",local = TRUE)
  
################################ observeEvent sections-----
  source("app_ObserveEvents.R",local = TRUE)
  
################################ renderPlot sections-----
  source("app_RenderPlots.R",local = TRUE)
  
################################ renderValueBox sections-----
  
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
  
################################ downloadHandler sections ----
  source("app_DownloadButtons.R",local = TRUE)
  
  
  ### table SDA annotations -----
  output$SDAannotations <- DT::renderDataTable(SDAannotation,
                                     options = list(paging = TRUE,    ## paginate the output
                                                    pageLength = 20,  ## number of rows to output for each page
                                                    scrollX = TRUE,   ## enable scrolling on X axis
                                                    scrollY = TRUE,   ## enable scrolling on Y axis
                                                    autoWidth = TRUE, ## use smart column width handling
                                                    server = FALSE,   ## use client-side processing
                                                    dom = 'Bfrtip',
                                                    buttons = c('csv', 'excel')#,
                                                    # columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                                    #                   list(targets = c(0, 8, 9), visible = FALSE))
                                     ),
                                     extensions = 'Buttons',
                                     selection = 'single', ## enable selection of a single row
                                     filter = 'bottom',              ## include column filters at the bottom
                                     rownames = FALSE                ## don't show row numbers/names
  )
  
  ### table SDA pseudotime  -----
  output$SDApseudotime<- DT::renderDataTable(SDApseudotime,
                                             options = list(paging = TRUE,    ## paginate the output
                                                            pageLength = 20,  ## number of rows to output for each page
                                                            scrollX = TRUE,   ## enable scrolling on X axis
                                                            scrollY = TRUE,   ## enable scrolling on Y axis
                                                            autoWidth = TRUE, ## use smart column width handling
                                                            server = FALSE,   ## use client-side processing
                                                            dom = 'Bfrtip',
                                                            buttons = c('csv', 'excel')#,
                                                            # columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                                            #                   list(targets = c(0, 8, 9), visible = FALSE))
                                             ),
                                             extensions = 'Buttons',
                                             selection = 'single', ## enable selection of a single row
                                             filter = 'bottom',              ## include column filters at the bottom
                                             rownames = FALSE                ## don't show row numbers/names
  )
  
}


shinyApp(ui, server)

