immuneCellsInputModeUI <- function(id){
  ns <- NS(id)
  
  fluidRow_12(
    h4(
      "The immune cells score is an estimation of an", 
      a(
        href = "https://doi.org/10.1186/s13059-017-1349-1", 
        target = "_blank",
        "enrichment (and not a deconvolution) algorithm"
      ), 
      "based on whole tissue RNAseq. Use with caution and consider the limitations explained in the paper."
    ),
    br(),
    selectInput(
      inputId = ns("cell_type"), 
      label = "Immune cell type", 
      choices = NULL
    ),
    brushPlotUI(ns("barplot")),
    uiOutput(ns("indicator"))
  )
}

immuneCellsInputMode <- function(input, output, session, TissuePrefilter){
  output$indicator <- renderUI({
    updateSelectInput(
      session = session,
      inputId = "cell_type",
      choices = getAvailableImmuneCellTypes(),
      selected = "CD8+ T-cells"
    )
    
    NULL
  })
  
  ImmuneCellData <- reactive({
    cellType <- input$cell_type
    validate(need(cellType, "select cell type"))
    getWaterfallDataImmuneCells(cellType, TissuePrefilter())
  })
  
  ImmuneCellPlot <- reactive({
    validate(need(ImmuneCellData(), "no immune cell prevalence data available..."))
    generateImmuneCellWaterfallPlot(ImmuneCellData(), fill = TissuePrefilter()$db_col)
  })
  
  ImmuneCellPlotCheck <- reactive({
    !is.null(input$cell_type) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "Score from {ry[1]} to {ry[2]}"
  )
  
  ImmuneCell_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = ImmuneCellPlot,
    checkExpr = ImmuneCellPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving immune cell data",
    value = 0.3
  )
  
  ImmuneCell_ti
}
