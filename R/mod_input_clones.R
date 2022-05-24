clonesInputModeUI <- function(id){
  ns <- NS(id)
  
  fluidRow_12(
    h4(
      "The clonal architecture was derived from ", 
      a(
        href = "https://doi.org/10.1371/journal.pgen.1007669", 
        target = "_blank",
        "allelic frequency of mutations and copy number estimations."
      ), 
      "Use with caution and consider the limitations explained in the paper."
    ),
    br(),
    selectInput(
      inputId = ns("score"), 
      label = "Clone score", 
      choices = c(
        "number of clones" = "number_of_clones", 
        "tree score" = "clone_tree_score"
      ), 
      width = "600px", 
      selected = "number of clones"
    ),
    brushPlotUI(ns("barplot"))
  )
}

clonesInputMode <- function(input, output, session, TissuePrefilter){
  CloneData <- reactive({
    score <- input$score
    if (is.null(score)) return()
    
    getWaterfallDataClones(score, TissuePrefilter())
  })
  
  ClonesPlot <- reactive({
    validate(need(CloneData(), "no clone data available..."))
    generateClonesWaterfallPlot(CloneData(), fill = TissuePrefilter()$db_col)
  })
  
  ClonesPlotCheck <- reactive({
    !is.null(input$score) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "Score from {ry[1]} to {ry[2]}"
  )
  
  Clones_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = ClonesPlot,
    checkExpr = ClonesPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving clone data",
    value = 0.3
  )
  
  Clones_ti
}
