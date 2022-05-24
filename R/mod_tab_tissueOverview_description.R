tissueOverviewDescriptionTabUI_main <- function(id){
  ns <- NS(id)
  
  fluidRow_12(
    br(),
    plotWrapperUI(ns("plot"), height = "800px")
  )
}

tissueOverviewDescriptionTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  barPieInputId <- ns("bar_or_pie")
  conditionId <- paste0("input['", barPieInputId, "']")
  
  div(
    selectInput(
      inputId = ns("property1"), 
      label = "Property 1", 
      choices = NULL, 
      selectize = FALSE
    ),
    selectInput(
      inputId = ns("property2"), 
      label = "Property 2", 
      choices = NULL,
      selectize = FALSE
    ),
    conditionalPanel(
      condition = paste0(conditionId, " && ", conditionId, " == 'bar'"),
      radioButtons(
        inputId = ns("percent"),
        label = "Show Property 1 by",
        choices = c("count", "percent"),
        selected = "count"
      )
    ),
    sliderInput(
      inputId = ns("nclasses"), 
      label = "Lump together if more than n classes",
      min = 1, 
      max = 30, 
      value = 10, 
      step = 1
    ),
    radioButtons(
      inputId = barPieInputId, 
      label = "Plot type", 
      choices = c(
        "Bar chart" = "bar", 
        "Pie chart" = "pie"
      ), 
      selected = "bar"
    ),
    hr(),
    textInput(
      inputId = ns("filename"), 
      label = "Choose file name", 
      value = "sample_description"
    ),
    downloadButton(
      outputId = ns("download"), 
      label = "Download"
    ),
    uiOutput(ns("indicator"))
  )
}

tissueOverviewDescriptionTab <- function(input, output, session, classSelection, classLabel, TissueAnnotation, TissueAnnotationFocus){
  output$indicator <- renderUI({
    ti_anno <- TissueAnnotationFocus()
    if (is.null(ti_anno)) return()
    
    ti_anno <- ti_anno %>% select_if(purrr::negate(is.numeric))
    choices <- setdiff(names(ti_anno), "tissuename")
    
    updateSelectInput(
      session = session,
      inputId = "property1",
      choices = choices,
      selected = isolate(input$property1) %||% choices[1]
    )
    
    updateSelectInput(
      session = session,
      inputId = "property2",
      choices = c("none", choices),
      selected = isolate(input$property2) %||% "none"
    )

    NULL
  })
  
  TissueDescriptionPlotExpr <- reactive({
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    p1 <- input$property1
    p2 <- input$property2
    nClasses <- input$nclasses
    
    validate(
      need(length(cs$class1) + length(cs$class2) > 0, "no tissues selected..."),
      need(cl$class1_name, "class name is missing..."),
      need(cl$class2_name, "class name is missing..."),
      need(p1, "property 1 is not available..."),
      need(p2, "property 2 is not available..."),
      need(nClasses, "number of classes is not ready yet ...")
    )
    
    usePercent <- (input$percent %||% "dummy") == "percent"
    generateClassSelectionPlot(
      sampleClasses = cs, 
      classLabel = cl, 
      prop1 = p1, 
      prop2 = p2, 
      n_classes = nClasses, 
      plot_type = input$bar_or_pie, 
      usePercent = usePercent,
      annotation = TissueAnnotation(), 
      annotationFocus = TissueAnnotationFocus()
    )
  })
  
  callModule(
    module = plotWrapper,
    id = "plot",
    PlotExpr = TissueDescriptionPlotExpr
  )
  
  output$download <- downloadHandler (
    filename <- function() { 
      paste0(input$filename, '.csv') 
    },
    content = function(file) {
      cs <- reactiveValuesToList(classSelection)
      cl <- reactiveValuesToList(classLabel)
      
      assignment <- stackClasses(cs, cl) %>%
        left_join(TissueAnnotation(), by = "tissuename")
      
      readr::write_excel_csv(assignment, file, na = "")
    }
  )
}
