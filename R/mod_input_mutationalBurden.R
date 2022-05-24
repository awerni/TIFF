mutationalBurdenInputModeUI <- function(id){
  ns <- NS(id)
  brushPlotUI(ns("barplot"))
}

mutationalBurdenInputMode <- function(input, output, session, TissuePrefilter){
  MutationalBurdenData <- reactive({
    getWaterfallDataMutationalBurden(TissuePrefilter())
  })
  
  MutationalBurdenPlot <- reactive({
    validate(need(MutationalBurdenData(), "no mutational burden data available..."))
    generateMutationalBurdenWaterfallPlot(MutationalBurdenData(), fill = TissuePrefilter()$db_col)
  })
  
  MutationalBurdenPlotCheck <- reactive({
    !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "Mutational fraction from {ry[1]} to {ry[2]} %"
  )
  
  MutationalBurden_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = MutationalBurdenPlot,
    checkExpr = MutationalBurdenPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 10,
    message = "Retrieving mutational burden data",
    value = 0.3
  )
  
  MutationalBurden_ti
}
