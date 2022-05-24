additionalPropertiesHallmarkComparisonTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

additionalPropertiesHallmarkComparisonTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    tabLayoutUI_sidebar(ns("tab")),
    hr(),
    extendedTableDownloadUI(ns("download"), "score", "gene set")
  )
}

additionalPropertiesHallmarkComparisonTab <- function(input, output, session, classSelection, classLabel, 
                                                      TissueAnnotationFocus, msigDBLink){
  # Reactives -----------------------------------------------------------------
  TableData <- reactive({
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    cs <- reactiveValuesToList(classSelection)
    df <- differentialHallmarkSets(cs, p = session) 
    
    validateFunctionResult(df)
    
    df %>%
      mutate(
        raw_gene_set = gene_set,
        gene_set = names(getHallmarkGeneSetChoices(raw_gene_set))
      ) %>%
      mutate_at(c("logFC", "P.Value", "adj.p.val"), signif, 3) %>% 
      select(gene_set, higher, logFC, P.Value, adj.p.val, raw_gene_set)
  })
  
  GeneSet <- reactive({
    s <- SelectedRow()
    validate(need(s, "no gene set selected"))
    as.list(TableData()[s, c("gene_set", "raw_gene_set")])
  })
  
  GeneSetScoreTissue <- reactive({
    geneSet <- GeneSet()$raw_gene_set
    cs <- reactiveValuesToList(classSelection)
    
    getTissueDataGeneSetById(geneSet, cs)
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    geneSet <- GeneSet()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateGeneSetPlot(
      df = GeneSetScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title = paste0("\n", geneSet$gene_set)
    )
  }
  
  PrettyTableData <- reactive({
    df <- TableData()
    req(df)
    
    cl <- reactiveValuesToList(classLabel)
    
    df %>%
      mutate(
        gene_set = getGeneSetLink(gene_set, msigDBLink),
        higher = styleHigherCol(higher == "class1", cl)
      ) %>%
      select(gene_set, logFC, higher, everything(), -raw_gene_set)
  })
  
  rowInfo <- callModule(
    module = tabLayout,
    id = "tab",
    plotFun = plotFun, 
    TableData = PrettyTableData
  )
  SelectedRow <- rowInfo$SelectedRow
  
  # Download ------------------------------------------------------------------
  RawGeneSet <- reactive({
    list(gene_set = GeneSet()$raw_gene_set)
  })
  
  callModule(
    module = extendedTableDownload,
    id = "download",
    Table = TableData,
    Subject = GeneSetScoreTissue,
    Item = RawGeneSet,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "gene_set"
  )
}
