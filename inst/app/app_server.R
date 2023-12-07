AddPer <- function(x, perc=0.1){
  x + x *  perc
}


datat <- as.data.frame(cbind(datat, results$scores[rownames(datat),]));
rownames(datat) <- datat$barcode

ColFac_DONR.ID  <- as.data.frame(datat)[rownames(results$scores), ]$donor

datat <- data.table(datat)


### render message menu ----
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


### DE genes ----
### dropdowns and datatable
observe({
  DEgenes_selectedItem <- input$DEgenes_selectedItem
  items <- names(DEgenes[[DEgenes_selectedItem]])
  updateSelectInput(session, "DEgenes_selectedDEDF", choices = items)
})

# Generate a dataframe based on the selected DEgenes_selectedItem and item
selected_data <- reactive({
  DEgenes_selectedItem <- input$DEgenes_selectedItem
  item <- input$DEgenes_selectedDEDF
  DEgenes[[DEgenes_selectedItem]][[item]]
})

# Render the resulting dataframe in a DataTable
output$DEgenes_dataTable <- DT::renderDataTable({
  selected_data()
}, options = list(
  paging = TRUE,
  pageLength = 20,
  scrollX = TRUE,
  scrollY = TRUE,
  autoWidth = TRUE,
  server = FALSE,
  dom = 'Bfrtip',
  buttons = c('csv', 'excel')
), extensions = 'Buttons',
selection = 'single',
filter = 'bottom',
rownames = FALSE)




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

################################ downloadHandler sections
source("app_DownloadButtons.R",local = TRUE)


### table SDA annotations
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

### table SDA pseudotime
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