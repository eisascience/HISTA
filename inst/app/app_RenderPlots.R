


# image2 sends pre-rendered images
output$homepage <- renderImage({
  return(list(
    src = "images/HISTA_Figs_v2.png",
    filetype = "image/png",
    alt = "Tutorial"
  ))
  
}, deleteFile = FALSE)


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
  
  tempDF <- PseudotimeGeneral_RX()
  
  # 
  # tempDF <- tSNE_GermCells_DF_Rx()
  # 
  # if(input$metaselect4 == "pseudotime") {
  #   MetaFac <- datat$PseudoTime
  # } else{
  #   
  #   if(input$metaselect4 == "celltype") {
  #     MetaFac <- (datat$FinalFinalPheno_old)
  #   } else {
  #     if(input$metaselect4 == "donrep"){
  #       MetaFac <- (datat$DonRep)
  #     } else {
  #       if(input$metaselect4 == "donor"){
  #         MetaFac <- (datat$donor)
  #       } else {
  #         if(input$metaselect4 == "COND.ID"){
  #           MetaFac <- (datat$COND.ID)
  #         } else {
  #           if(input$metaselect4 == "experiment"){
  #             MetaFac <- (datat$experiment)
  #           } else {
  #             
  #           }
  #         }
  #       }
  #     }
  #   }
  #   
  # }
  # 
  # 
  # tempDF$MetFacZ <- MetaFac
  # 
  # tempDF <- tempDF[!is.na(tempDF$tSNE1),]
  # 
  # #reduce memory
  
  # tempDF = tempDF[sample(1:nrow(tempDF), 500, replace = F), ]
  
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
  
  tempDF <- PseudotimeGeneral_RX()
  
  merge_sda_melt <- reshape2::melt(tempDF, id.vars = c("barcode","tSNE1", "tSNE2", "GeneExpr", "MetFacZ", "PT"))
  
  print("melted table")
  
  # print(head(rownames(tempDF)))
  # print(head(rownames(Scores)))
  
  # tempDF <- tempDF[!is.na(tempDF$tSNE1),]
  # Scores <-Scores[rownames(tempDF),]
  # print(head(merge_sda_melt))
  # plot(merge_sda_melt$PT, 
  #      merge_sda_melt$value)
  
  
  
  ggpp = ggplot(merge_sda_melt, aes(PT, value, colour=(MetFacZ))) +
    geom_point(alpha=1, size=.2) +
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 50), se = FALSE) +
    # stat_smooth(aes(PT, value), size=1, alpha = 0.6, method = "gam", formula = y ~ s(x, k = 20), se = F) +#colour="black",
    ylab("Cell Component Score") +
    xlab("Pseudotime") +
    # ggtitle(paste0("SDA Comp: ", as.numeric(input$ComponentNtext3)))+
    theme_bw() +
    theme(legend.position = "none") +
    ylim(-8,8)
  
  print("ggpp made")

  if(input$metaselect4 == "pseudotime") {
    ggpp =  ggpp +  scale_color_viridis()
  } else {
    ggpp =  ggpp +  scale_colour_manual(values=col_vector)  + facet_wrap(~MetFacZ,
                                                                         ncol=3, 
                                                                         scales = "fixed")
  }
  
  print("color type corrected")
  
  ggpp
  
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

output$tSNEperCellType_GEX <- renderPlot({
  tSNEwSDAScoreProjPerCT_GEX_Rx()
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

output$tSNE_geneExpr <- renderPlot({
  tSNE_geneExpr_Rx()
  
})

output$tSNEperCellType_meta <- renderPlot({
  tSNEwMetaPerCT_Rx()
})

# output$tSNE_somaWLN <- renderPlot({
#   cowplot::plot_grid(tSNE_somaWLN_Pheno3_Rx(), 
#                      tSNE_somaWLN_COND.ID_Rx(), 
#                      tSNE_somaWLN_DONR.ID_Rx(), 
#                      tSNE_somaWLN_nCount_RNA_Rx(),
#                      ncol=2)
# })


## Soma only with LN19 ------

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


## LC with Zhao20 and LN19 ------

output$DimRedux_LConly_donors_Rx <- renderPlot({
  DimRedux_LConly_donors_Rx()
})

output$DimRedux_LConly_phenotype_Rx <- renderPlot({
  DimRedux_LConly_phenotype_Rx()
})

output$DimRedux_LConlyZhao_phenotype_Rx <- renderPlot({
  DimRedux_LConlyZhao_phenotype_Rx()
})
output$DimRedux_LConlyZhao_donors_Rx <- renderPlot({
  DimRedux_LConlyZhao_donors_Rx()
})
output$DimRedux_LConlyZhao_phenotypeProp_Rx <- renderPlot({
  DimRedux_LConlyZhao_phenotypeProp_Rx()
})
output$DimRedux_LConlyZhao_KeyGenesViolin_Rx <- renderPlot({
  DimRedux_LConlyZhao_KeyGenesViolin_Rx()
})

## GO enrichment plots ------

output$GOpos <- renderPlot({
  
  SDAGOpos_Rx()

})

output$GOneg <- renderPlot({
  
  SDAGOneg_Rx()
  
})


## chrom loading location -----

output$ChrLoc <- renderPlot({
  
  ChrLocLoadings_Rx()
  
})


## SDA scpres across ------

output$SDAScoresAcross <- renderPlot({
  SDAScoresAcross_Rx()
  
})


## Enrichment ------

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