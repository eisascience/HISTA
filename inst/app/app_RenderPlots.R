output$CellScoreOrderSDA <- renderPlot({
  
  Scores <- ComboTopSDAgenes_Rx()$Scores
  Meta <- ComboTopSDAgenes_Rx()$Meta
  
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
  
  
  tempDF <- tSNE_GermCells_DF_Rx()
  
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
  tempLS <- SDAScoresChiPos_Rx()
  
  pheatmap::pheatmap(tempLS$obj,
                     cluster_cols = tempLS$clustStat, cluster_rows = tempLS$clustStat,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                     labels_col = tempLS$label_col)
})

output$SDAScoresChiNeg <- renderPlot({
  tempLS <- SDAScoresChiNeg_Rx()
  
  pheatmap::pheatmap(tempLS$obj,
                           cluster_cols = tempLS$clustStat, cluster_rows = tempLS$clustStat,
                           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                           labels_col = tempLS$label_col)
})

output$GeneExprSigMeta <- renderPlot({
  (GeneExprSigMeta_Rx())
})

output$GeneExprSigMeta2 <- renderPlot({
  
  GeneExpr <- GeneExprAcroosCellType_DF_Rx()$GeneExpr
  my_comparisons <- GeneExprAcroosCellType_DF_Rx()$my_comparisons
  
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

output$tSNEwSDAScoreProj <- renderPlot({
  tSNEwSDAScoreProj_Rx()
})

output$tSNEperCellType <- renderPlot({
  tSNEwSDAScoreProjPerCT_Rx()
})

# output$tSNE_somaWLN <- renderPlot({
#   cowplot::plot_grid(tSNE_somaWLN_Pheno3_Rx(), 
#                      tSNE_somaWLN_COND.ID_Rx(), 
#                      tSNE_somaWLN_DONR.ID_Rx(), 
#                      tSNE_somaWLN_nCount_RNA_Rx(),
#                      ncol=2)
# })


output$tSNE_somaWLN_Pheno3_Rx <- renderPlot({
  tSNE_somaWLN_Pheno3_Rx()
})
output$tSNE_somaWLN_COND.ID_Rx <- renderPlot({
  tSNE_somaWLN_COND.ID_Rx()
})
output$tSNE_somaWLN_DONR.ID_Rx <- renderPlot({
  tSNE_somaWLN_DONR.ID_Rx()
})
output$tSNE_somaWLN_nCount_RNA_Rx <- renderPlot({
  tSNE_somaWLN_nCount_RNA_Rx()
})



output$tSNEperCellType_meta <- renderPlot({
  tSNEwMetaPerCT_Rx()
})


output$tSNEwMetaLegend <- renderPlot({
  
  legend <- cowplot::get_legend(tSNEwMeta_Rx())
  
  #grid.newpage()
  grid.draw(legend)
})

#renderPlotly
output$tSNEwMeta <- renderPlot({
  tSNEwMeta_Rx()+
    theme(legend.position = "none", aspect.ratio=1,
          legend.title = element_blank())
})

output$plot4 <- renderPlot({
  
  tempDF <- as.data.frame(tSNE_AllCells_DF_Rx())
  
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
  
  SDAGOpos_Rx()

})

output$GOneg <- renderPlot({
  
  SDAGOneg_Rx()
  
})

output$ChrLoc <- renderPlot({
  
  ChrLocLoadings_Rx()
  
})

output$SDAScoresAcross <- renderPlot({
  SDAScoresAcross_Rx()
  
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
  
  plotEnrich(GeneSetsDF=SDA_Top100pos, 
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