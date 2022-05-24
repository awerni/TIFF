survivalTabUI <- function(id){
  ns <- NS(id)
  fluidRow(
    column_2(survivalTabUI_sidebar(ns)),
    column_10(plotWrapperUI(ns("plot"), height = "800px"))
  )
}

survivalTabUI_sidebar <- function(ns){
  choices <- c("hide", "show")
  list(
    radioButtons(
      inputId = ns("confidence"), 
      label = "Show confidence interval", 
      choices = choices, 
      selected = "hide"
    ),
    radioButtons(
      inputId = ns("pval"), 
      label = "Show p-value", 
      choices = choices, 
      selected = "show"
    ),
    radioButtons(
      inputId = ns("risk"), 
      label = "Show risk table", 
      choices = choices, 
      selected = "show"
    )
  )
}

survivalTab <- function(input, output, session, classSelection, classLabel, PatientAnnotationFuller){
  PatientAnnotationFocus <- reactive({
    PatientAnnotationFuller() %>% filter(tissuename %in% c(classSelection$class1, classSelection$class2))  
  })
  
  SurvivalData <- reactive({
    cs <- reactiveValuesToList(classSelection)
    
    validate(
      need(cs$class1, "goto input tab and select tissues for class1"),
      need(cs$class2, "goto input tab and select tissues for class2"),
      need(input$confidence, "show confidence button is not ready yet"),
      need(input$pval, "show p-value button is not ready yet"),
      need(input$risk, "show risk table button is not ready yet"),
      need(PatientAnnotationFocus(), "no patient data available")
    )
    
    cl <- reactiveValuesToList(classLabel)
    data <- PatientAnnotationFocus() %>% 
      inner_join(stackClasses(cs, cl, return_factor = TRUE), by = "tissuename") %>%
      select(tissuename, class, days_to_last_followup, days_to_death)
    
    validate(need(any(!is.na(data$days_to_last_followup)), "no survival data available"))
    
    survivalAnalysis(data)
  })
  
  KaplanMeierPlotExpr <- reactive({
    analysisData <- SurvivalData()
    req(analysisData)
    
    p <- generateKaplanMeierPlot(
      data = analysisData,
      plot_p = input$pval == "show",
      plot_ci = input$confidence == "show",
      plot_risk = input$risk == "show"
    )
    
    validate(need(p, "no survival curves to be compared"))
    p
  })
  
  callModule(
    module = plotWrapper,
    id = "plot",
    PlotExpr = KaplanMeierPlotExpr
  )
}
