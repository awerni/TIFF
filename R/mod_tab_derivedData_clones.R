derivedDataClonesTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

derivedDataClonesTabUI_sidebar <- function(id){
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
    extendedTableDownloadUI(ns("download"), "score", "clones")
  )
}

derivedDataClonesTab <- function(input, output, session, classSelection, classLabel){
  Result_clonesStatistics <- reactive({
    cs <- reactiveValuesToList(classSelection)
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    
    df <- compareClones(cs, p = session)
    validateFunctionResult(df)
    df
  })
  
  ScoreTissue <- reactive({
    cs <- reactiveValuesToList(classSelection)
    score <- SelectedCloneScore()[["clone score"]]
    
    getTissueDataClones(cs, score) %>%
      mutate(`clone score` = score) %>%
      select(tissuename, `clone score`, score = !!score)
  })
  
  SelectedCloneScore <- reactive({
    s <- SelectedRow()
    validate(need(s, "no clone score selected"))
    
    df <- TableData()
    as.list(df[s, "clone score", drop = FALSE])
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    cloneScore <- SelectedCloneScore()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateClonesPlot(
      df = ScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title = paste0("\n", cloneScore[["clone score"]])
    )
  }
  
  TableData <- reactive({
    res <- Result_clonesStatistics()
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
    Item = SelectedCloneScore,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "clone score"
  )
}
