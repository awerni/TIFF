copyNumberInputModeUI <- function(id){
  ns <- NS(id)
  
  list(
    fluidRow(
      column_2(
        br(),
        textInput(
          inputId = ns("symbol"), 
          label = "gene symbol", 
          value = ""
        )
      ),
      column_10(
        br(),
        h3(textOutput(ns("genename")))
      )
    ), 
    fluidRow_12(
      radioButtons(
        inputId = ns("plot_type"),
        label = "Select using",
        choices = c(
          "relative DNA copy number" = "score", 
          "copy number class" = "types"
        ),
        selected = "score"
      )
    ),
    fluidRow_12(brushPlotUI(ns("barplot")))
  )
}

copyNumberInputMode <- function(input, output, session, species, TissuePrefilter){
  GeneCopynumber <- reactive({
    symbol <- input$symbol
    if (is.null(symbol)) return()
    if (nchar(symbol) < 3) return()
    
    getGeneFromSymbol(symbol, species)
  })
  
  CopynumberData <- reactive({
    if (is.null(GeneCopynumber())) return()
    getWaterfallDataCopyNumber(GeneCopynumber()$ensg, TissuePrefilter())
  })
  
  output$genename <- renderText(GeneCopynumber()$name)
  
  CopyNumberPlot <- reactive({
    data <- CopynumberData()
    plotType <- input$plot_type
    
    validate(
      need(data, "no copy number data available..."),
      need(plotType, "select plot type")
    )
    
    if (plotType == "score"){
      generateCopynumberWaterfallPlot(data, fill = TissuePrefilter()$db_col)
    } else {
      generateCopynumberBarPlot(data)
    }
    
  })
  
  CopyNumberPlotCheck <- reactive({
    symbol <- input$symbol
    !is.null(symbol) && nchar(symbol) >= 3 && 
      !is.null(TissuePrefilter()) && !is.null(species) && !is.null(input$plot_type)
  })
  
  textCallback <- makeTextCallback(
    "Relative copy number from {ry[1]} to {ry[2]} copies"
  )
  
  CopyNumber_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = CopyNumberPlot,
    checkExpr = CopyNumberPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving copy number data",
    value = 0.3
  )
  
  CopyNumber_ti
}
