
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
library(knitr)


library(HISTA)


pathi = getwd()

source("app_Fxs.R",local = TRUE)


# source(paste0(pathi, "inst/app/app_Fxs.R"))

# list2env(readRDS( paste0(pathi, "/data/HISTAv1_dataLS_feb2022.rds")), envir = globalenv()) #for deploy on shinyapps

#monkey shiny server do has SCRATCH_DIR defined
if (Sys.getenv("SCRATCH_DIR") != "") {
  
  pathi = paste0(Sys.getenv("SCRATCH_DIR"), "/data")
  load.data.path = paste0(pathi, "/ConradLab/HISTA/HISTAv1_dataLS_June2023.rds" )
}  else {
  load.data.path = paste0(pathi, "/HISTAv1_dataLS_June2023.rds" )
  
  if(!file.exists(load.data.path)) {
    load.data.path = paste0(pathi, "/data/HISTAv1_dataLS_June2023.rds" )
    
  }
  
}
list2env(readRDS(load.data.path), envir = globalenv())


# palettes <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

col_vector = c("#7FC97F", "#38170B", "#BEAED4", "#BF1B0B", "#FFC465", "#386CB0", 
  "#66ADE5", "#F0027F", "#252A52", "#BF5B17", "#999999", "#666666", 
  "#E69F00", "#1B9E77", "#7570B3", "#F0E442", "#B07AA1", "#F28E2B", 
  "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
  "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295", 
  "#FABFD2", "#E7298A", "#D4A6C8", "#9D7660", "#D7B5A6","#e6194b", 
  "#666666", "#3cb44b", "#A6CEE3", "#ffe119", "#1F78B4", "#B15928")


## ui ------

ui <- dashboardPage(
  dashboardHeader(title = "HISTA v2.9.5"
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
      menuItem("Fingerprinting (heatmap)", tabName = "chisqrSDAenrichHeatmaps", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per pathology (boxplot)", tabName = "geneexprstatsig", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per cell type (boxplot)", tabName = "celltypestatsig", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Expr per cell type (2D)", tabName = "celltypeGeneExpr2D", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Cell score per celltype (2D)", tabName = "tsnepercelltype", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Metadata per celltype (2D)", tabName = "tsnepercelltype_meta", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Gene Correlations", tabName = "GeneCor", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Component Corrrelations", tabName = "CompCor", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      #menuItem("Score order per. Comp", tabName = "CellScoreOrderingSDA", icon = icon("plus"),
               #badgeLabel = "final", badgeColor = "green"),
      menuItem("Pseudotime Meta (germ only)", tabName = "pseudotimeSDA", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Pseudotime Gene (germ only)", tabName = "pseudotimeSDAgene", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Pseudotime Component Index", tabName = "pseudotimeSDAIndex", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("plus"),#asterisk
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Top Loaded Components", tabName = "TopLoadedComps", icon = icon("plus"),
              badgeLabel = "final", badgeColor = "green"),
      menuItem("lncRNAs", tabName = "lncRNA_expr_toploaded", icon = icon("plus"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Soma only W. LN19", tabName = "somaWLN", icon = icon("asterisk"), #minus
               badgeLabel = "final", badgeColor = "green"),
      menuItem("LC only W. Zhao21", tabName = "LConly", icon = icon("asterisk"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("User Manual", tabName = "usermanual", icon = icon("book"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Version History", tabName = "versionhistory", icon = icon("book"),
               badgeLabel = "final", badgeColor = "green"),
      menuItem("Conrad Lab", icon = icon("book"), 
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
     

     
     ## versionhistory -----
     
     tabItem(tabName = "versionhistory",
             h2("Version History"),
             fluidRow(
               uiOutput('VerHistHTML')
               

             )
     ),
     
     
     ## usermanual -----
     
     tabItem(tabName = "usermanual",
             h2("User Manual"),
             fluidRow(
               uiOutput('UserManualHTML')
               
               
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
                                 "UMAP - Batch-removed SDA-CellScores" = "umapBrSDACS",
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
                                 "Experiments" = "experiment",
                                 "CellCycle" = "Phase"
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
             h2("Gene expr by cell type"),
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
               box(title = "DimRedux", status = "warning", solidHeader = TRUE,
                 radioButtons("data2", "Data origin:",
                              c("tSNE - Batch-removed SDA-CellScores" = "tsneBrSDACS",
                                "UMAP - Batch-removed SDA-CellScores" = "umapBrSDACS",
                                "tSNE - DropSim Norm. DGE" = "tsneDSNDGE",
                                #"Batch-removed-Imputed-DGE" = "DGEimp",
                                #"DromSim-Normalized-DGE" = "DGEnorm",
                                #"tSNE - Seurat-Norm. DGE" = "tsneSNDGE",
                                "tSNE - Batch-removed-Imputed-DGE" = "tsneImpDGE"
                              )),
                 width = 5
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
               
               box(title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                   collapsible = TRUE,
                   radioButtons("clustStat_chisqr", "Hierarchical cluster:",
                                c("Cluster" = "True",
                                  "No Cluster" = "False"
                                )),
                 radioButtons("metaselect_chisqr", "Metadata Selection:",
                              c("Cell types" = "celltype",
                                "Donor-replicates" = "donrep",
                                "Donors" = "donor",
                                "Condition" = "COND.ID",
                                "Experiments" = "experiment",
                                "CellCycle" = "Phase"
                              ), selected = "experiment"),
                 
                 width = 5#, background = "black"
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
               
               
             )
     ),
      
     ##tSNE per cell type cell score -----
     
     
      tabItem(tabName = "tsnepercelltype",
              h2("SDA score by cell type"),
              fluidRow(
                box(title = "DimRedux", status = "warning", solidHeader = TRUE,
                    radioButtons("data3", "Data origin:",
                                 c("tSNE - Batch-removed SDA-CellScores" = "tsneBrSDACS",
                                   "UMAP - Batch-removed SDA-CellScores" = "umapBrSDACS",
                                   "tSNE - DropSim Norm. DGE" = "tsneDSNDGE",
                                   #"Batch-removed-Imputed-DGE" = "DGEimp",
                                   #"DromSim-Normalized-DGE" = "DGEnorm",
                                   #"tSNE - Seurat-Norm. DGE" = "tsneSNDGE",
                                   "tSNE - Batch-removed-Imputed-DGE" = "tsneImpDGE"
                                 )),
                    width = 5
                ),
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
              h2("Metadata by cell type"),
              fluidRow(
                box(title = "DimRedux", status = "warning", solidHeader = TRUE,
                    radioButtons("data4", "Data origin:",
                                 c("tSNE - Batch-removed SDA-CellScores" = "tsneBrSDACS",
                                   "UMAP - Batch-removed SDA-CellScores" = "umapBrSDACS",
                                   "tSNE - DropSim Norm. DGE" = "tsneDSNDGE",
                                   #"Batch-removed-Imputed-DGE" = "DGEimp",
                                   #"DromSim-Normalized-DGE" = "DGEnorm",
                                   #"tSNE - Seurat-Norm. DGE" = "tsneSNDGE",
                                   "tSNE - Batch-removed-Imputed-DGE" = "tsneImpDGE"
                                 )),
                    width = 5
                ),
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
                                 "Experiments" = "experiment",
                                 "CellCycle" = "Phase"
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
                    textInput("ComponentName_pseudo", "SDA component search (numerical):", "1"),
                    radioButtons("metaselect_pseudo", "Metadata Selection:",
                                 c("Pseudotime" = "pseudotime",
                                   "Cell types" = "celltype",
                                   "Donor-replicates" = "donrep",
                                   "Donors" = "donor",
                                   "Condition" = "COND.ID",
                                   "Experiments" = "experiment",
                                   "CellCycle" = "Phase"
                                 )),
                    width = 3),
                box(title = "Germ Cell Only tSNE", status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,
                    plotOutput("tSNEPseudoSDA", height = 400),
                ),
                box(title = "Pseudotime", status = "primary", solidHeader = TRUE,
                    collapsible = TRUE,
                    downloadButton("PseudotimeSDA_download"),
                    plotOutput("PseudotimeSDA", height = 400),
                    width = 10
                )
                
              )
      ),
     
     
     ### Pseudotime gene -----
     
     tabItem(tabName = "pseudotimeSDAgene",
             h2("Gene expr on Pseudotime of Germ Cells Only"),
             fluidRow(
               box(title = "Inputs", status = "warning", solidHeader = TRUE,
                   textInput("Genetext_pseudo", "Gene search (text input):", "PRM1"),
                   radioButtons("metaselect_pseudo_gene", "Metadata Selection:",
                                c("Cell types" = "celltype",
                                  "Donor-replicates" = "donrep",
                                  "Donors" = "donor",
                                  "Condition" = "COND.ID",
                                  "Experiments" = "experiment",
                                  "CellCycle" = "Phase"
                                ), selected = "Phase"),
                   width = 5
               ),
               box(title = "Gene Expr Pseudotime", status = "primary", solidHeader = TRUE,
                   collapsible = TRUE,
                   downloadButton("PseudotimeSDA_gene_download"),
                   plotOutput("PseudotimeSDAgene", height = 400),
                   plotOutput("PseudotimeSDAgeneMeta", height = 400),
                   
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
                  plotOutput("PosEnrichPlot"),
                  width = 10
                ),
                
                box(
                  title = "Negative Loadings", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE,
                  plotOutput("NegEnrichPlot"),
                  width = 10
                )
              )
      ),
     
     
     ## TopLoadedComps -----
     
     tabItem(tabName = "TopLoadedComps",
             fluidRow(
               box(
                 title = "Inputs", status = "warning", solidHeader = TRUE,
                 "Multiple formatting of gene sets accepted", 
                 br(), "List can be seperated by comma e.g. from ", 
                 br(), "   or spaces e.g. from Excel", 
                 br(), "Also, single or double quotes or not",
                 br(), "This works best with larger sets of genes!",
                 br(), "As an example set, here are all the -AS genes",
                 #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                 textInput("GeneSet_TLC", "A set of genes", "'AADACL2-AS1', 'C21orf62-AS1', 'DLGAP2-AS1', 'TRIM31-AS1', 'MAMDC2-AS1', 'KCNQ1-AS1', 
               'MIF-AS1', 'KANSL1-AS1', 'MAPT-AS1', 'LENG8-AS1', 'DDX39B-AS1', 'FAM181A-AS1', 'ITPK1-AS1', 
               'SERTAD4-AS1', 'ZFY-AS1', 'GLIS2-AS1', 'ZNF436-AS1', 'EPHA1-AS1', 'MCPH1-AS1', 'SLC14A2-AS1', 
               'DHRS4-AS1', 'PPP4R1-AS1', 'C1RL-AS1', 'AQP4-AS1', 'MCM8-AS1', 'SIRPG-AS1', 'TTC39C-AS1', 
               'ISM1-AS1', 'MACROD2-AS1', 'TTC3-AS1', 'ASMTL-AS1', 'CACNA1C-AS4', 'CACNA1C-AS2', 
               'KCNK15-AS1', 'USP12-AS2', 'ARHGAP5-AS1', 'STK4-AS1', 'PAN3-AS1', 'HORMAD2-AS1',
               'TEX26-AS1', 'SLC25A21-AS1', 'DLGAP1-AS3', 'URB1-AS1', 'ZNF295-AS1', 'ZNF674-AS1',
               'ITCH-AS1', 'EP300-AS1', 'FUT8-AS1', 'USP12-AS1', 'DIAPH3-AS1', 'MIS18A-AS1', 'PCDH9-AS4', 
               'NDFIP2-AS1', 'BRWD1-AS2', 'BRWD1-AS1', 'BMP7-AS1', 'STX18-AS1', 'PRR34-AS1', 'RBM26-AS1', 
               'FGF13-AS1', 'WASF3-AS1', 'DSG2-AS1', 'FOXN3-AS1', 'PCDH9-AS3', 'RAP2C-AS1', 'UCKL1-AS1', 
               'B4GALT1-AS1', 'TAPT1-AS1', 'ZNF337-AS1', 'CTBP1-AS', 'DIO2-AS1', 'TRPM2-AS', 'TSPEAR-AS1', 
               'TSPEAR-AS2', 'FRMD6-AS1', 'ST8SIA6-AS1', 'EFCAB6-AS1', 'FRMD6-AS2', 'PCCA-AS1', 'CDKN2B-AS1', 
               'SCEL-AS1', 'MCM3AP-AS1', 'ADORA2A-AS1', 'GPC5-AS2', 'LY86-AS1', 'PCDH9-AS1', 'PROSER2-AS1', 
               'TSC22D1-AS1', 'SGMS1-AS1', 'FAM83C-AS1', 'MORC2-AS1', 'PARD3-AS1', 'F10-AS1', 'JARID2-AS1', 
               'LNX1-AS1', 'FAM170B-AS1', 'ELOVL2-AS1', 'IGF2-AS', 'MGAT3-AS1', 'UBOX5-AS1', 'PAXBP1-AS1', 
               'ITGB2-AS1', 'HM13-AS1', 'ZNF503-AS2', 'GRID1-AS1', 'JMJD1C-AS1', 'PRMT5-AS1', 'KIZ-AS1', 
               'EDNRB-AS1', 'MACC1-AS1', 'ZNF341-AS1', 'NAV2-AS5', 'GPC6-AS1', 'FGF14-AS1', 'MORF4L2-AS1', 
               'SYP-AS1', 'FGF14-AS2', 'BDNF-AS', 'PLCG1-AS1', 'A2ML1-AS1', 'GABRG3-AS1', 'BTBD9-AS1', 
               'TRAM2-AS1', 'PPP1R26-AS1', 'SACS-AS1', 'RNASEH2B-AS1', 'PLBD1-AS1', 'SLC39A12-AS1', 'COL5A1-AS1', 
               'CSE1L-AS1', 'DSG1-AS1', 'OGFR-AS1', 'VPS13A-AS1', 'PRRX2-AS1', 'FAM53B-AS1', 'MCHR2-AS1', 'SNCA-AS1', 
               'WDFY3-AS1', 'NRSN2-AS1', 'DBH-AS1', 'IQCH-AS1', 'ARNTL2-AS1', 'THAP9-AS1', 'ADNP-AS1', 'SLC25A25-AS1', 
               'LZTS1-AS1', 'ZNRF3-AS1', 'ODF2-AS1', 'GAS6-AS1', 'HOXC-AS1', 'TEX36-AS1', 'FOCAD-AS1', 'RARA-AS1', 
               'LMO7-AS1', 'EDRF1-AS1', 'ZBTB46-AS1', 'TPT1-AS1', 'USP2-AS1', 'SLC25A30-AS1', 'GSN-AS1', 'IFT74-AS1', 
               'HTR2A-AS1', 'DIAPH3-AS2', 'ARAP1-AS2'"),
                 width = 5
               )),
               
               
        
               
               box(
                 title = "Top positive loaded comps", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 plotOutput("TopLoadComp_Plot"),
                 width = 5
               ),
             
             box(
               title = "Enrichment across components", status = "primary", solidHeader = TRUE,
               collapsible = TRUE,
               plotOutput("TopLoadedBarplot"),
               width = 10
             ),
             
             box(
               title = "Cor(loadings) of SDA components", status = "primary", solidHeader = TRUE,
               collapsible = TRUE,
               downloadButton("CompCorCustPlot_download"),
               plotOutput("CompCorCustPlot"),
               width = 10
             )
             
             
     ),
     
     
     
     ## lncRNA_expr_toploaded -----
     
     tabItem(tabName = "lncRNA_expr_toploaded",
             fluidRow(
               
               
               box(
                 title = "Overlap of lncRNAs (biomaRt 04/2022)", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 plotOutput("lncRNA_Venn"),
                 width = 5
               ),

               box(
                 title = "Distribution of lncRNAs vs Random gene set", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 plotOutput("lncRNA_toploaded"),
                 width = 5
               ),
               
               box(
                 title = "Enrichment across components", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 plotOutput("lncRNA_BarplotSDA"),
                 width = 10
               ),
               
               box(
                 title = "Selection of component", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 textInput("ComponentNtext_lncRNATopLoaded", "SDA component search (numerical):", "1"),
                 # plotOutput("lncRNA_topLoaded"),
                 width = 5
               ),
               box(
                 title = "Top Loaded lncRNAs Neg and Pos", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 textOutput("lncRNAgenes_Neg_txt"),
                 textOutput("lncRNAgenes_Pos_txt"),
                 width = 10
               ),
               
               # valueBoxOutput("lncRNAgenes", width = NULL),
               
               
             )
     ),
     
     
     
     ## CompCor -----
     
     tabItem(tabName = "CompCor",
             fluidRow(
               
               valueBoxOutput("CorPlot_CellType1", width = 3),
               
               box(
                 title = "Cor(loadings) of SDA components", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 textInput("CompCorSDAnum", "Top 150 genes of SDA component (numerical):", "1"),
                 sliderInput(inputId = "CompCorSDAnum_ngene",
                             label = "Number of top genes:",
                             min = 10,
                             max = 150,
                             value = 24),
                 downloadButton("CompCorPlot_download"),
                 plotOutput("CompCorPlot"),
                 width = 10
               ),
               box(
                 title = "Top Loaded genes", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 textOutput("top_Neg_txt"),
                 textOutput("top_Peg_txt"),
                 width = 10
               )
               
               
             )
     ),
     
     ## CompCor -----
     
     tabItem(tabName = "GeneCor",
             fluidRow(
               
               # valueBoxOutput("CorPlot_CellType1", width = 3),
               
               box(
                 radioButtons("celltypeselect_geneCor", "Celltype Selection:",
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
                              ), selected = "germ_UnDiffSgSct"),
                 width = 5
               ),
               
               box(
                 title = "Cor(GEX) gene expression correlation", status = "primary", solidHeader = TRUE,
                 collapsible = TRUE,
                 textInput("GeneSet_geneCor", "A set of genes", "'NANOS1', 'NANOS2', 'NANOS3', 'DMRT1', 'KIT', 'PLPPR3', 'CCK', 'HES5', 'PLPPR5', 'RHCE', 'SCT', 'FAM25G', 'DPPA4', 'TDRD1', 'MAGEA4', 'STRA8', 'DPPA4', 'TDRD1', 'MAGEA4', 'STRA8', 'GFRA1', 'BCL6B', 'ID4', 'PIWIL4', 'DNMT1', 'DMRT1', 'KIT', 'DMC1', 'SYCP3', 'ETV5', 'EGR4'"),
                 
                 downloadButton("GeneCorPlot_download"),
                 plotOutput("GeneCorPlot"),
                 width = 10
               )
               
               
             )
     )
     
     
      
    )#tabitems
  )#dashboardBody
  
  
)#dashboardPage




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
  
################################ lncRNA sections-----  
  source("app_lncRNA.R",local = TRUE)
  
################################ renderValueBox sections-----
  
  output$VerHistHTML <- renderUI({
    # paste0(pathi, "/data/HISTAv1_dataLS_feb2022.rds")
    # HTML(markdown::markdownToHTML(knit('VersionHistory.md', quiet = TRUE)))
    includeHTML('VersionHistory.html')
  })
  
  output$UserManualHTML <- renderUI({
    # paste0(pathi, "/data/HISTAv1_dataLS_feb2022.rds")
    # HTML(markdown::markdownToHTML(knit('VersionHistory.md', quiet = TRUE)))
    includeHTML('UserManual.html')
  })
  
  
  # renderText()
  
  
  
  output$lncRNAgenes_Neg_txt <- renderText({ 
    paste0("neg = ", paste0(lncLS$NegLoaded_top[[ as.numeric(input$ComponentNtext_lncRNATopLoaded)]], collapse = ", "))
  })
  output$lncRNAgenes_Pos_txt <- renderText({ 
   paste0("pos = ", paste0(lncLS$PosLoaded_top[[ as.numeric(input$ComponentNtext_lncRNATopLoaded)]], collapse = ", "))
  })
  
  
  
  output$top_Neg_txt <- renderText({ 
    paste0("neg = ", paste0( SDA_Top100neg[1:input$CompCorSDAnum_ngene, as.numeric(input$CompCorSDAnum)], collapse = ", "))
  })
  
 
  output$top_Peg_txt <- renderText({ 
    paste0("pos = ", paste0( SDA_Top100pos[1:input$CompCorSDAnum_ngene, as.numeric(input$CompCorSDAnum)], collapse = ", "))
  })
  


# output$lncRNAgenes <- renderValueBox({
#   valueBox(
#     value = paste0(paste0("pos = ", paste0(lncLS$PosLoaded_top[[ as.numeric(input$ComponentNtext_lncRNATopLoaded)]], collapse = ", ")), " -- -- -- ",
#                    paste0("neg = ", paste0(lncLS$NegLoaded_top[[ as.numeric(input$ComponentNtext_lncRNATopLoaded)]], collapse = ", "))),
#     subtitle = paste0("SDAV", input$ComponentNtext_lncRNATopLoaded, sep=""),
#     icon = icon("area-chart"),
#     color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
#   )
# })
  
  output$CellType1 <- renderValueBox({
    valueBox(
      value = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = StatFac[paste0("SDAV", input$ComponentNtext, sep=""),6],
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  
  output$CorPlot_CellType1 <- renderValueBox({
    valueBox(
      value = StatFac[paste0("SDAV", input$CompCorSDAnum, sep=""),2], #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = StatFac[paste0("SDAV", input$CompCorSDAnum, sep=""),6],
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

