derivedDataMetabolicsTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

derivedDataMetabolicsTabUI_sidebar <- function(id){
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
    extendedTableDownloadUI(ns("download"), "score", "metabolics")
  )
}

derivedDataMetabolicsTab <- function(input, output, session, classSelection, classLabel){
  Result_metabolicsStatistics <- reactive({
    cs <- reactiveValuesToList(classSelection)
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    
    df <- compareMetabolics(cs, p = session)
    validateFunctionResult(df)
    df
  })
  
  ScoreTissue <- reactive({
    cs <- reactiveValuesToList(classSelection)
    pathway <- SelectedMetabolicPathway()[["metabolic pathway"]]
    
    getTissueDataMetabolicsById(pathway, cs)
  })
  
  SelectedMetabolicPathway <- reactive({
    s <- SelectedRow()
    validate(need(s, "no metabolic pathway selected"))
    
    df <- TableData()
    as.list(df[s, "metabolic pathway", drop = FALSE])
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    pathway <- SelectedMetabolicPathway()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateMetabolicsPlot(
      df = ScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title = paste0("\n", pathway[["metabolic pathway"]])
    )
  }
  
  TableData <- reactive({
    res <- Result_metabolicsStatistics()
    validate(need(res, "no statistics results available..."))
    
    cl <- reactiveValuesToList(classLabel)
    
    res %>% 
      mutate_at(c("diff", "p.value", "adj.P.Value"), signif, 3) %>%
      mutate(higher = styleHigherCol(diff > 0, cl)) %>%
      relocate(higher, .after = method)
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
    Subject = ScoreTissue,
    Item = SelectedMetabolicPathway,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "metabolic pathway"
  )
}
