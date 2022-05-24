geneExpressionInputModeUI <- function(id){
  ns <- NS(id)
  
  list(
    fluidRow(
      column_3(
        br(),
        textInput(
          inputId = ns("symbol"), 
          label = "gene symbol", 
          value = ""
        )
      ),
      column_9(
        br(), 
        h3(textOutput(ns("genename")))
      )
    ), 
    fluidRow_12(brushPlotUI(ns("barplot")))
  )
}

geneExpressionInputMode <- function(input, output, session, species, TissuePrefilter){
  GeneExpression <- reactive({
    symbol <- input$symbol
    if (is.null(symbol)) return()
    if (nchar(symbol) < 3) return()
    
    getGeneFromSymbol(symbol, species)
  })
  
  ExpressionData <- reactive({
    if (is.null(GeneExpression())) return()
    a <- getWaterfallDataGeneExpression(GeneExpression()$ensg, TissuePrefilter())
  })
  
  output$genename <- renderText(GeneExpression()$name)
  
  ExpressionPlot <- reactive({
    validate(need(ExpressionData(), "waiting for expression data..."))
    generateExpressionWaterfallPlot(ExpressionData(), fill = TissuePrefilter()$db_col)
  })
  
  ExpressionPlotCheck <- reactive({
    symbol <- input$symbol
    !is.null(symbol) && nchar(symbol) >= 3 && !is.null(TissuePrefilter()) && !is.null(species)
  })
  
  textCallback <- makeTextCallback(
    "Gene expression from {ry[1]} to {ry[2]} TPM"
  )
  
  Expression_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = ExpressionPlot,
    checkExpr = ExpressionPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving expression data",
    value = 0.3
  )
  
  Expression_ti
}
