signalingInputModeUI <- function(id){
  ns <- NS(id)
  
  fluidRow_12(
    h4(
      "The signaling pathway data are based on the paper ", 
      a(
        href = "https://doi.org/10.1016/j.cell.2018.03.035", 
        target = "_blank",
        "Oncogenic Signaling Pathways in The Cancer Genome Atlas"
      ), 
      "Use with caution and consider the limitations explained in the paper."
    ),
    br(),
    selectInput(
      inputId = ns("pathway"), 
      label = "Signaling pathway", 
      choices = NULL, 
      width = "300px"
    ),
    brushPlotUI(ns("barplot")),
    uiOutput(ns("indicator"))
  )
}

signalingInputMode <- function(input, output, session, TissuePrefilter){
  output$indicator <- renderUI({
    updateSelectInput(
      session = session,
      inputId = "pathway",
      choices = getAvailableSignalingPathways(),
      selected = "rtk_ras"
    )
    
    NULL
  })
  
  SignalingData <- reactive({
    pathway <- input$pathway
    validate(need(pathway, "select pathway"))
    df <- getInputDataSignaling(pathway, TissuePrefilter())
    
    validate(need(nrow(df) > 0, "no data available"))
    df
  })
  
  SignalingPlot <- reactive({
    validate(need(SignalingData(), "no signaling data available..."))
    generateSignalingInputBarplot(SignalingData(), fill = TissuePrefilter()$db_col)
  })
  
  SignalingPlotCheck <- reactive({
    !is.null(input$pathway) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback("Signaling from {ry[1]} to {ry[2]}")
  
  Signaling_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = SignalingPlot,
    checkExpr = SignalingPlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 1,
    message = "Retrieving signaling data",
    value = 0.3
  )
  
  Signaling_ti
}
