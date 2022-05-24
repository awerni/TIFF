derivedDataImmuneCellsTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

derivedDataImmuneCellsTabUI_sidebar <- function(id){
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
    extendedTableDownloadUI("download", "score", "immune_cells")
  )
}

derivedDataImmuneCellsTab <- function(input, output, session, classSelection, classLabel){
  # Reactives -----------------------------------------------------------------
  Result_immuneCellStatistics <- reactive({
    cs <- reactiveValuesToList(classSelection)
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    
    df <- compareImmuneCells(cs, p = session)
    validateFunctionResult(df)
    df
  })
  
  ScoreTissue <- reactive({
    cs <- reactiveValuesToList(classSelection)
    cellType <- SelectedImmuneCellType()[["cell type"]]
    
    getTissueDataImmuneCellsById(cellType, cs)
  })
  
  SelectedImmuneCellType <- reactive({
    s <- SelectedRow()
    validate(need(s, "no immune cell selected"))
    
    df <- TableData()
    as.list(df[s, "cell type", drop = FALSE])
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    cellType <- SelectedImmuneCellType()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateImmuneCellsPlot(
      df = ScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title =  paste0("\n", cellType[["cell type"]])
    )
  }
  
  TableData <- reactive({
    res <- Result_immuneCellStatistics()
    validate(need(res, "no statistics results available..."))
    
    cl <- reactiveValuesToList(classLabel)
    
    res %>% 
      mutate_at(c("log2FC", "p.value", "adj.P.Value"), signif, 3) %>%
      mutate(higher = styleHigherCol(log2FC > 0, cl)) %>%
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
    Item = SelectedImmuneCellType,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "cell type"
  )
}
