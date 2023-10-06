# UI --------------------------------------------------------------------------
hallmarkSetsInputModeUI <- function(id){
  ns <- NS(id)
  
  list(
    fluidRow(
      selectInput(
        inputId = ns("geneset"), 
        label = "Gene set:", 
        multiple = FALSE, 
        choices = NULL
      ),
      uiOutput(ns("indicator"))
    ),
    brushPlotUI(ns("barplot"))
  )
}

# Server ----------------------------------------------------------------------
hallmarkSetsInputMode <- function(input, output, session, TissuePrefilter, 
                                  geneSets, default = "HALLMARK_APOPTOSIS"){
  
  output$indicator <- renderUI({
    choices <- getHallmarkGeneSetChoices(sort(geneSets))
    
    updateSelectInput(
      inputId = "geneset",
      choices = choices,
      selected = default
    )
    
    NULL
  })
  
  HallmarkSetData <- reactive({
    geneSet <- input$geneset
    
    req(geneSet)
    getWaterfallDataHallmarkSet(geneSet, TissuePrefilter())
  })
  
  HallmarkSetPlot <- reactive({
    hsd <- HallmarkSetData()
    validate(need(hsd, "no data found"))
    generateHallmarkSetWaterfallPlot(hsd, fill = TissuePrefilter()$db_col)
  })
  HallmarkSetPlotCheck <- reactive({
    !is.null(input$geneset) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "Score values from {ry[2]} to {ry[1]}"
  )
  
  HallmarkSet_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = HallmarkSetPlot,
    checkExpr = HallmarkSetPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 0,
    message = "Retrieving hallmark set data",
    value = 0.3
  )
  
  HallmarkSet_ti
}
