metabolicsInputModeUI <- function(id){
  ns <- NS(id)
  
  fluidRow_12(
    h4(
      "The metabolics score is based on ", 
      a(
        href = "https://doi.org/10.1186/s12943-018-0895-9", 
        target = "_blank",
        "single sample GSEA from metabolics pathways"
      ), 
      "based on whole tissue RNAseq. Use with caution and consider the limitations explained in the paper."
    ),
    br(),
    selectInput(
      inputId = ns("pathway"), 
      label = "Metabolic pathway", 
      choices = NULL, 
      width = "600px"
    ),
    brushPlotUI(ns("barplot")),
    uiOutput(ns("indicator"))
  )
}

metabolicsInputMode <- function(input, output, session, TissuePrefilter){
  output$indicator <- renderUI({
    updateSelectInput(
      session = session,
      inputId = "pathway",
      choices = getAvailableMetabolicPathways(),
      selected = "metabolism of carbohydrates"
    )
    
    NULL
  })
  
  MetabolicsData <- reactive({
    pathway <- input$pathway
    validate(need(pathway, "select pathway"))
    getWaterfallDataMetabolics(pathway, TissuePrefilter())
  })
  
  MetabolicsPlot <- reactive({
    validate(need(MetabolicsData(), "no metabolics data available..."))
    generateMetabolicsWaterfallPlot(MetabolicsData(), fill = TissuePrefilter()$db_col)
  })
  MetabolicsPlotCheck <- reactive({
    !is.null(input$pathway) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "Score from {ry[1]} to {ry[2]}"
  )
  
  Metabolics_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = MetabolicsPlot,
    checkExpr = MetabolicsPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving metabolics data",
    value = 0.3
  )
  
  Metabolics_ti
}
