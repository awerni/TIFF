proteinExpressionRppaInputModeUI <- function(id){
  ns <- NS(id)
  
  list(
    fluidRow(
      column_3(
        br(),
        selectizeInput(
          inputId = ns("antibody"), 
          label = "antibody:", 
          choices = NULL
        )
      ),
      column_4(
        br(), 
        uiOutput(ns("antibody_desc"))
      ),
      column_5()
    ), 
    fluidRow_12(brushPlotUI(ns("barplot"))),
    uiOutput(ns("indicator"))
  )
}

proteinExpressionRppaInputMode <- function(input, output, session, TissuePrefilter){
  output$indicator <- renderUI({
    updateSelectizeInput(
      session = session, 
      inputId = "antibody", 
      choices = sort(getAntibodyInformation()$antibody), 
      selected = "MEK1_pS217_S221", 
      server = TRUE
    )
    
    NULL
  })
  
  ProteinExpressionData <- reactive({
    antibody <- input$antibody
    if (is.null(antibody) || !nzchar(antibody)) return()
    getWaterfallDataRppa(antibody, TissuePrefilter())
  })
  
  ProteinExpressionPlot <- reactive({
    d <- ProteinExpressionData()
    if (is.null(d)) return()
    generateProteinWaterfallPlot(d, fill = TissuePrefilter()$db_col)
  })
  
  output$antibody_desc <- renderUI({
    a <- getAntibodyInformation(input$antibody) %>% unlist()
    if (length(a) == 0) return()
    
    tagList(
      h3(paste0(a["vendor"], " (", a["catalog_number"], ")")),
      h4(a["validation_status"])
    )
  })
  
  PlotCheck <- reactive({
    antibody <- input$antibody
    !is.null(antibody) && nzchar(antibody) && !is.null(TissuePrefilter())
  })
  
  textCallback <- makeTextCallback(
    "RPPA Protein expression score from {ry[2]} to {ry[1]}"
  )
  
  ProteinExpression_ti <- callModule(
    module = brushPlot,
    id = "barplot",
    plotExpr = ProteinExpressionPlot,
    checkExpr = PlotCheck,
    textCallback = textCallback,
    defaultCutoff_y = 0,
    message = "Retrieving protein expression data",
    value = 0.3
  )
  
  ProteinExpression_ti
}
