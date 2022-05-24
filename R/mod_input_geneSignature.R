geneSignatureInputModeUI <- function(id, geneSignatures){
  ns <- NS(id)
  
  choices <- geneSignatures$signature
  
  list(
    fluidRow(
      column_8(
        br(),
        containerDT(ns("table"))
      ),
      column_4()
    ),
    brushPlotUI(ns("barplot"))
  )
}

geneSignatureInputMode <- function(input, output, session, TissuePrefilter, geneSignatures){
  output$table <- DT::renderDT(
    expr = {
      table <- geneSignatures
      validate(need(!is.null(table) && nrow(table) > 0, "no gene signatures found"))
      table %>% 
        mutate(description = paste0('<a href="', hyperlink, '" target="_blank">', description, '</a>')) %>%
        select(-hyperlink)
    },
    rownames = FALSE,
    filter = "none",
    options = list(
      pageLength = 5, 
      search = list(regex = TRUE)
    ),
    escape = FALSE,
    selection = list(
      mode = "single", 
      selected = 1, 
      target = "row"
    )
  )
  
  SelectedSignature <- reactive({
    s <- input$table_rows_selected
    validate(need(s, "no gene signature selected"))
    geneSignatures$signature[s]
  })
  
  GeneSignatureData <- reactive({
    signature <- SelectedSignature()
    req(signature)
    getWaterfallDataGeneSignatures(signature, TissuePrefilter())
  })
  
  GeneSignaturePlot <- reactive({
    gsd <- GeneSignatureData()
    validate(need(gsd, "no data found"))
    generateGeneSignatureWaterfallPlot(gsd, fill = TissuePrefilter()$db_col)
  })
  GeneSignaturePlotCheck <- reactive({
    !is.null(SelectedSignature()) && !is.null(TissuePrefilter())
  })
  
  textCallback <-  makeTextCallback(
    "Gene signature values from {ry[2]} to {ry[1]}"
  )
  
  GeneSignature_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = GeneSignaturePlot,
    checkExpr = GeneSignaturePlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving gene signature data",
    value = 0.3
  )
  
  GeneSignature_ti
}
