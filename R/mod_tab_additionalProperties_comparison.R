additionalPropertiesClassComparisonTabUI_main <- function(id){
  ns <- NS(id)
  
  list(
    fluidRow_12(plotWrapperUI(ns("plot"))),
    fluidRow_12(containerDT(ns("table")))
  )
}

additionalPropertiesClassComparisonTabUI_sidebar <- function(id){
  ns <- NS(id)
  NULL
}

additionalPropertiesClassComparisonTab <- function(input, output, session, classSelection, classLabel, TissueAnnotationFocus, Result_otherPropStatistic){
  TableData <- reactive({
    validate(need(Result_otherPropStatistic(), "no statistics results available..."))
    res <- Result_otherPropStatistic() %>% 
      filter(! property %in% c("Mutational fraction", "Vendor name"))
    res[,setdiff(colnames(res), "higher")]
  })
  
  SelectedOtherProperty <- reactive({
    s <- input$table_rows_selected
    validate(need(s, "no property selected..."))
    as.list(TableData()[s, c("property", "data")])
  })
  
  OtherPropPlotExpr <- reactive({
    validate(
      need(nrow(TissueAnnotationFocus()) > 0, "no tissue data available..."),
      need(SelectedOtherProperty(), "no property selected")
    )
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    generateOtherPropertiesPlot(
      sampleClasses = cs,
      classLabel = cl,
      sAnno = TissueAnnotationFocus(),
      property = SelectedOtherProperty()
    )
  })
  
  callModule(
    module = plotWrapper,
    id = "plot",
    PlotExpr = OtherPropPlotExpr
  )
  
  output$table <- DT::renderDataTable(
    expr = TableData(),
    rownames = FALSE,
    filter = "none",
    options = list(
      pageLength = 10, 
      search = list(regex = TRUE)
    ),
    escape = FALSE,
    selection = list(
      mode = "single", 
      selected = 1, 
      target = "row"
    )
  )
}
