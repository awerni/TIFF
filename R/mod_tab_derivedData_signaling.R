derivedDataSignalingTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"), FALSE)
}

derivedDataSignalingTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    tabLayoutUI_sidebar(
      id = ns("tab"),
      defaults = list(
        left = "bar",
        right = "pie"
      ),
      additionalChoices = c("bar plot" = "bar", "pie chart" = "pie"),
      hidden = c("roc", "point", "violin", "box"),
      useMiddlePlot = FALSE
    ),
    hr(),
    extendedTableDownloadUI(ns("download"), "state", "signaling")
  )
}

derivedDataSignalingTab <- function(input, output, session, classSelection, classLabel){
  Result_signalingStatistics <- reactive({
    cs <- reactiveValuesToList(classSelection)
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    
    df <- compareSignaling(cs, p = session)
    validateFunctionResult(df)
    df
  })
  
  ScoreTissue <- reactive({
    cs <- reactiveValuesToList(classSelection)
    pathway <- SelectedSignalingPathway()[["signaling pathway"]]
    
    getTissueDataSignalingById(pathway, cs) %>%
      mutate(active = ifelse(active, "active", "inactive"))
  })
  
  SelectedSignalingPathway <- reactive({
    s <- SelectedRow()
    validate(need(s, "no signaling pathway selected"))
    
    df <- TableData()
    as.list(df[s, "signaling pathway", drop = FALSE])
  })
  
  plotFun <- function(plotType){
    pathway <- SelectedSignalingPathway()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateSignalingPlot(
      df = ScoreTissue(),
      ca = makeClassAssignment(cs, cl),
      plotType = plotType,
      title = paste("\n", pathway, "pathway")
    )
  }
  
  TableData <- reactive({
    res_sig <- Result_signalingStatistics()
    validate(need(res_sig, "no statistics results available..."))
    
    cl <- reactiveValuesToList(classLabel)
    
    res_sig %>%
      mutate(
        class1.active = glue::glue(
          "{signif(class1.percentage * 100, 3)}% ({class1.n_active}/{class1.n_total})"
        ),
        class2.active = glue::glue(
          "{signif(class2.percentage * 100, 3)}% ({class2.n_active}/{class2.n_total})"
        ),
        pathway = toupper(pathway)
      ) %>%
      select(
        `signaling pathway` = pathway,
        !!rlang::sym(paste0(cl$class1_name, ".active")) := class1.active,
        !!rlang::sym(paste0(cl$class2_name, ".active")) := class2.active, 
        higher, p.value, adj.P.Value
      ) %>%
      mutate_at(c("p.value", "adj.P.Value"), signif, 3) %>%
      mutate(higher = styleHigherCol(higher == "class1", cl))
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
    Item = SelectedSignalingPathway,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "signaling pathway"
  )
}
