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