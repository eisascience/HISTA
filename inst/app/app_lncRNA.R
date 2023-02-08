
assign("lincrna", readRDS(paste0(pathi, "/data/lincrnaDF.rds")), envir=globalenv())



# ggvenn::ggvenn(list(names1 = lngrna_in_data1, names2 = lngrna_in_data2))


lncrna_overlap = unique(lincrna$hgnc_symbol)[unique(lincrna$hgnc_symbol)  %in% unique(c(unique(lincrna$hgnc_symbol)[unique(lincrna$hgnc_symbol)  %in% colnames(results$loadings[[1]])], 
                                                                                             c(colnames(results$loadings[[1]])[grep("-AS", colnames(results$loadings[[1]]))],
                                                                                               colnames(results$loadings[[1]])[grep("LINC", colnames(results$loadings[[1]]))]))) ]
assign("lncrna_overlap", lncrna_overlap, envir=globalenv())




EnumSDA = function(geneV = NULL, Ladings = NULL){
  
  NegLoaded_top = lapply(1:150, function(x){
    SDA_comp = Ladings[x, ]
    top_load_genes = sort(SDA_comp, decreasing = F) %>% head(200)
    geneV_in_top = intersect(geneV, names(top_load_genes))
    geneV_in_top
  })
  
  PosLoaded_top = lapply(1:150, function(x){
    SDA_comp = Ladings[x, ]
    top_load_genes = sort(SDA_comp, decreasing = T) %>% head(200)
    geneV_in_top = intersect(geneV, names(top_load_genes))
    geneV_in_top
  })
  
  N_overlap_NegLoaded = unlist(lapply(NegLoaded_top, function(x){
    length(x)
  }))
  
  
  N_overlapp_PosLoaded = unlist(lapply(PosLoaded_top, function(x){
    length(x)
  }))
  
  dfm1 = data.frame(NegN = N_overlap_NegLoaded, 
                    PosN= N_overlapp_PosLoaded, 
                    comp = paste0("SDA", 1:150), 
                    row.names = paste0("SDA", 1:150))
  
  return(list(DF = dfm1,
              PosLoaded_top = PosLoaded_top, 
              NegLoaded_top = NegLoaded_top))
  
}


assign("lncLS", EnumSDA(geneV = lncrna_overlap, Ladings = results$loadings[[1]]), envir=globalenv())



sigComps = apply(lncLS$DF[,1:2], 1, function(x){
  any(x>22.3)
})

assign("sigComps", sigComps, envir=globalenv())


plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.9, position="identity", aes(y = ..density..), color="black", bins = 50) +
    geom_density(alpha=0.4) +
    # geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  
  plt + guides(fill=guide_legend(title=label_column))
}
