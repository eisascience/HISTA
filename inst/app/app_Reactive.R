################################ Reactive sections

tSNE_AllCells_DF_Rx <- reactive({
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

tSNE_GermCells_DF_Rx <- reactive({
  # data, tsneBrSDACS, tsneDSNDGE, tsneSNDGE
  # print(input$data)
  
  tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC4sperm", "Tsne2_SDAQC4sperm")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
  
  
  rownames(tempDF) <- datat$barcode
  
  # tempDF <- tempDF[!is.na(tempDF$tSNE1), ]
  
  tempDF$GeneExpr <- rep(0, nrow(tempDF))
  
  return(tempDF)
})

tSNEwSDAScoreProj_Rx <- reactive({
  #ggplotly
  tempDF <- tSNE_AllCells_DF_Rx()
  ggplot(cbind(tempDF, SDAComp=datat[,get(paste0("SDAV", input$ComponentNtext, sep=""))]), 
         aes(tSNE1, tSNE2, color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
    geom_point(size=0.1) + theme_bw() +
    scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
    theme(legend.position = "bottom", aspect.ratio=1) + 
    ggtitle(paste0("SDAV", input$ComponentNtext, " \n", 
                   StatFac[paste0("SDAV", input$ComponentNtext, sep=""),2], sep="")) + 
    simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
  
})

tSNEwSDAScoreProjPerCT_Rx <- reactive({
  #ggplotly
  tempDF <- tSNE_SDA_CT_Rx()
  # print(head(tempDF))
  
  tempMeta <- datat[,get(paste0("SDAV", input$ComponentNtext_tsnepercelltype, sep=""))]
  names(tempMeta) <- datat$barcode
  tempMeta <- tempMeta[rownames(tempDF)]
  print(head(tempMeta))
  
  
  ggplot(cbind(tempDF, SDAComp=tempMeta), 
         aes(tSNE1, tSNE2, color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
    geom_point(size=0.1) + theme_bw() +
    scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
    theme(legend.position = "bottom", aspect.ratio=1) + 
    ggtitle(paste0("SDAV", input$ComponentNtext_tsnepercelltype, " \n", 
                   StatFac[paste0("SDAV", input$ComponentNtext_tsnepercelltype, sep=""),2], sep="")) + 
    simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
  
})

tSNEwMeta_Rx <- reactive({
  tempDF <- tSNE_AllCells_DF_Rx()
  
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
  ggplot(tempDF, aes(tSNE1, tSNE2, color=factor(MetaFac))) +
    geom_point(size=0.1)+ theme_bw() +
    theme(legend.position = "right", aspect.ratio=1,
          legend.title = element_blank()) +
    ggtitle("t-SNE - Final Pheno") +
    scale_color_manual(values=(col_vector)) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol =2))  +
    coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)
  
  
})

tSNE_SDA_CT_Rx <- reactive({
  
  if(input$tsnepercelltype_ctselect == "leydig"){
    MyCells <- datat[datat$FinalFinalPheno == "Leydig",]$barcode
  } else {
    if(input$tsnepercelltype_ctselect == "sertoli"){
      MyCells <- datat[datat$FinalFinalPheno == "Sertoli",]$barcode
    } else {
      if(input$tsnepercelltype_ctselect == "endothelial"){
        MyCells <- datat[datat$FinalFinalPheno == "Endothelial",]$barcode
      } else {
        if(input$tsnepercelltype_ctselect == "myeloid"){
          MyCells <- datat[datat$FinalFinalPheno %in% c("Macrophage-M2", "Macrophage-M1"),]$barcode
        } else {
          if(input$tsnepercelltype_ctselect == "adaptive"){
            MyCells <- datat[datat$FinalFinalPheno %in% c("Lymphoid-Bcell", "Lymphoid-Tcell"),]$barcode
          } else {
            if(input$tsnepercelltype_ctselect == "germ"){
              MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct", 
                                                            "Gamete_Meiotic_Pach_Dip_2nd_Scts",
                                                            "Gamete_Meiotic_preLep_Lep_Zyg_Scts",
                                                            "Gamete_RoundSpermatid",
                                                            "Gamete_UndiffSg"),]$barcode
            } else {
              if(input$tsnepercelltype_ctselect == "all"){
                MyCells <- datat$barcode
              } else {
                if(input$tsnepercelltype_ctselect == "germ_DiffSgSct"){
                  MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_DifferentiatingSgSct"),]$barcode
                } else {
                  if(input$tsnepercelltype_ctselect == "germ_PrePachSct"){
                    MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_preLep_Lep_Zyg_Scts"),]$barcode
                  } else {
                    if(input$tsnepercelltype_ctselect == "germ_PachSct"){
                      MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_Meiotic_Pach_Dip_2nd_Scts"),]$barcode
                    }  else {
                      if(input$tsnepercelltype_ctselect == "germ_Std"){
                        MyCells <- datat[datat$FinalFinalPheno %in% c("Gamete_RoundSpermatid"),]$barcode
                      } else {
                        if(input$tsnepercelltype_ctselect == "germ_UnDiffSgSct"){
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
  
  
  
  
  
  tempDF <- as.data.frame(datat)[, c("Tsne1_SDAQC4b", "Tsne2_SDAQC4b")]; colnames(tempDF) <- c("tSNE1", "tSNE2")
  
  rownames(tempDF) <- datat$barcode
  tempDF <- tempDF[MyCells, ]
  
  tempDF$GeneExpr <- rep(0, nrow(tempDF))
  
  return(tempDF)
  
})

GeneExprPerCellType_DF_Rx <- reactive({
  
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

GeneExprAcroosCellType_DF_Rx <- reactive({
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

GeneExprSigMeta_Rx <- reactive({
  
  GeneExpr <- GeneExprPerCellType_DF_Rx()$GeneExpr
  my_comparisons <- GeneExprPerCellType_DF_Rx()$my_comparisons
  
  
  CellType = input$celltypeselect2
  
  
  TestName = "Wilcox Rank Sum"
  
  ggboxplot(GeneExpr, x = "meta", y = "gene", palette = "jco",
            add = "jitter", col="meta") + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    theme_bw() + 
    ggtitle( paste0(as.character(input$Genetext2), " expression :: ", TestName, " test :: ", CellType)) + 
    xlab("") + ylab(as.character(input$Genetext2))
  
  
})

ComboTopSDAgenes_Rx <- reactive({
  
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

ComboTopSDAgenes_DF_Rx <- reactive({
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

PseudotimeSDA_Rx <- reactive({
  
  Scores <- results$scores
  tempDF <- tSNE_GermCells_DF_Rx()
  
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
  # print(head(merge_sda_melt))
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

SDAScoresChiNeg_Rx <- reactive({
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
  
  # print(head(NegCompsDF))
  
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
  
  # HM <- pheatmap::pheatmap((t(ChiT$residuals)),
  #                    cluster_cols = clustStat, cluster_rows = clustStat,
  #                    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
  #                    labels_col = paste0(rownames(NegCompsDF), " sd_", ChiResSD))
  outLS <- list(obj=(t(ChiT$residuals)),
       clustStat = clustStat,
       label_col = paste0(rownames(NegCompsDF), " sd_", ChiResSD))
  
  # return((HM$gtable))
  return(outLS)
  
  
})

SDAScoresChiPos_Rx <- reactive({
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
  
  # print(head(PosCompsDF))
  
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
  
  # HM <- pheatmap::pheatmap((t(ChiT$residuals)),
  #                    cluster_cols = clustStat, cluster_rows = clustStat,
  #                    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
  #                    labels_col = paste0(rownames(PosCompsDF), " sd_", ChiResSD))
  # 
  # return(grid.draw(HM$gtable))
  
  outLS <- list(obj=(t(ChiT$residuals)),
                clustStat = clustStat,
                label_col = paste0(rownames(PosCompsDF), " sd_", ChiResSD))
  
  # return((HM$gtable))
  return(outLS)
})

SDAScoresAcross_Rx <- reactive({
  
  # ColFac_DONR.ID <- CDID()
  
  if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
    print("No Comp")
  } else {
    SDAScores <- results$scores
    ComponentN <- as.numeric(input$ComponentNtext)
    
    ggplot(data.table(cell_index = 1:nrow(SDAScores), 
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
    
    
  }
})

SDAGOpos_Rx <- reactive({
  if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
    print("No GO")
  } else {
    go_volcano_plot(component = paste("V", input$ComponentNtext, "P", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
    
  }
})

SDAGOneg_Rx <- reactive({
  if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
    print("No GO")
  } else {
    go_volcano_plot(component = paste("V", input$ComponentNtext, "N", sep=""))+ theme_bw()+ theme(aspect.ratio = 1)
    
  }
})

ChrLocLoadings_Rx <- reactive({
  if(! (as.numeric(input$ComponentNtext) %in% 1:150)){
    print("No Comp")
  } else {
    pgl <- genome_loadings(results$loadings[[1]][as.numeric(input$ComponentNtext),], 
                           label_both = T, 
                           max.items = as.numeric(input$NoOfGenes), 
                           gene_locations =   gene_locations,
                           chromosome_lengths = chromosome.lengths) + theme(aspect.ratio = .5)
    return(pgl)
    
  }
})

