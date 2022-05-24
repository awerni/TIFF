mutationMutationalFractionTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

mutationMutationalFractionTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    tabLayoutUI_sidebar(
      id = ns("tab"),
      defaults = list(
        left = "roc",
        middle = "point",
        right = "violin"
      )
    ),
    hr(),
    extendedTableDownloadUI(ns("download"), "mutational fraction")
  )
}

mutationMutationalFractionTab <- function(input, output, session, classSelection, classLabel, TissueAnnotationFocus, Result_otherPropStatistic){
  # Reactives -----------------------------------------------------------------
  SelectedProperty <- reactive({
    s <- SelectedRow()
    validate(need(s, "no property selected..."))
    
    df <- TableData()
    as.list(df[s, "property", drop = FALSE])
  })
  
  PropertyData <- reactive({
    property <- as.character(SelectedProperty()$property)
    TissueAnnotationFocus() %>%
      mutate(property = property) %>%
      select(tissuename, tumortype, property, value = !!property)
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType){
    property <- SelectedProperty()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateMutationalFractionPlot(
      df = PropertyData(), 
      ca = makeClassAssignment(cs, cl),
      plotType = plotType,
      title = paste0("\n", property$property)
    )
  }
  
  TableData <- reactive({
    res <- Result_otherPropStatistic()
    validate(need(res, "no statistics results available..."))
    
    cl <- reactiveValuesToList(classLabel)
    res <- res %>% filter(property == "Mutational fraction")
    
    validate(need(nrow(res) > 0, "no mutational fraction data available..."))
    
    res %>%
      mutate(higher = styleHigherCol(higher == "class1", cl)) %>%
      select(property, method, higher, p.value, statistic)
  })
  
  rowInfo <- callModule(
    module = tabLayout,
    id = "tab",
    plotFun = plotFun,
    TableData = TableData
  )
  SelectedRow <- rowInfo$SelectedRow
  
  # Download ------------------------------------------------------------------
  callModule(
    module = extendedTableDownload,
    id = "download",
    Table = TableData,
    Subject = PropertyData,
    Item = SelectedProperty,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "property"
  )
}
