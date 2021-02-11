output$TXTall <- downloadHandler(
  filename = function(){
    paste("data-TopGenes_SDA_negNpos", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    # print(ComboTopSDAgenes_DF_Rx())
    write.csv(ComboTopSDAgenes_DF_Rx(), file, row.names=F)
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



output$tSNEwSDAScoreProj_download <- downloadHandler(
  filename = function(){
    paste("tSNEwSDAScoreProj_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9.5, height =9, compress = T, pointsize = 15)
    plot(tSNEwSDAScoreProj_Rx())
    dev.off()
  })

output$tSNEwMeta_download <- downloadHandler(
  filename = function(){
    paste("tSNEwMeta_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9.5, height =9, compress = T, pointsize = 15)
    plot(tSNEwMeta_Rx())
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

output$tsnepercelltype_download <- downloadHandler(
  filename = function(){
    paste("tsnepercelltype_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNEwSDAScoreProjPerCT_Rx())
    dev.off()
  })

output$tsnepercelltype_meta_download <- downloadHandler(
  filename = function(){
    paste("tsnepercelltype_meta_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNEwMetaPerCT_Rx())
    dev.off()
  })

output$SDAScoresChiNeg_download <- downloadHandler(
  filename = function(){
    paste("SDAScoresChiNeg_download_15x6", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 15, height =6, compress = T, pointsize = 15)
    # plot(SDAScoresChiNeg_Rx())
    tempLS <- SDAScoresChiNeg_Rx()
    pdf(file, width = 15, height =6, compress = T, pointsize = 15)
    # plot(SDAScoresChiPos_Rx())
    pheatmap::pheatmap(tempLS$obj,
                       cluster_cols = tempLS$clustStat, cluster_rows = tempLS$clustStat,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                       labels_col = tempLS$label_col)
    dev.off()
  })

output$SDAScoresChiPos_download <- downloadHandler(
  filename = function(){
    paste("SDAScoresChiPos_download_15x6", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    tempLS <- SDAScoresChiPos_Rx()
    pdf(file, width = 15, height =6, compress = T, pointsize = 15)
    # plot(SDAScoresChiPos_Rx())
    pheatmap::pheatmap(tempLS$obj,
                       cluster_cols = tempLS$clustStat, cluster_rows = tempLS$clustStat,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                       labels_col = tempLS$label_col)
    dev.off()
  })

output$SDAScoresAcross_download <- downloadHandler(
  filename = function(){
    paste("SDAScoresAcross_download_15x6", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 15, height =6, compress = T, pointsize = 15)
    plot(SDAScoresAcross_Rx())
    dev.off()
  })


output$SDAGOneg_download <- downloadHandler(
  filename = function(){
    paste("SDAGOneg_download_9x9", Sys.Date(), ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(SDAGOneg_Rx())
    dev.off()
  })
output$SDAGOpos_download <- downloadHandler(
  filename = function(){
    paste("SDAGOpos_download_9x9", Sys.Date(), ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(SDAGOpos_Rx())
    dev.off()
  })


output$ChrLoc_download <- downloadHandler(
  filename = function(){
    paste("ChrLoc_download_12x6", Sys.Date(), ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file, width = 12, height =6, compress = T, pointsize = 15)
    plot(ChrLocLoadings_Rx())
    dev.off()
  })




output$tsnesomaonlywln_phenotype_download <- downloadHandler(
  filename = function(){
    paste("tsnesomaonlywln_phenotype_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNE_somaWLN_Pheno3_Rx())
    dev.off()
  })
output$tsnesomaonlywln_condition_download <- downloadHandler(
  filename = function(){
    paste("tsnesomaonlywln_condition_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNE_somaWLN_COND.ID_Rx())
    dev.off()
  })
output$tsnesomaonlywln_donor_download <- downloadHandler(
  filename = function(){
    paste("tsnesomaonlywln_donor_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNE_somaWLN_DONR.ID_Rx())
    dev.off()
  })
output$tsnesomaonlywln_ncount_download <- downloadHandler(
  filename = function(){
    paste("tsnesomaonlywln_ncount_download_9x9", Sys.Date(), ".pdf", sep = "")
    # "test.pdf"
  },
  content = function(file) {
    pdf(file, width = 9, height =9, compress = T, pointsize = 15)
    plot(tSNE_somaWLN_nCount_RNA_Rx())
    dev.off()
  })