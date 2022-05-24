expressionTabUI <- function(id){
  ns <- NS(id)
  
  tbl <- tabsetPanel(
    id = ns("tabset"), 
    tabPanel("result"),
    tabPanel("volcano"),
    tabPanel("heatmap"),
    tabPanel("GSEA"),
    tabPanel("gene clustering")
  )
  
  makeCondition <- function(x, negate = FALSE){
    tabInput <- paste0("input['", ns("tabset"), "']")
    operator <- if (negate) " != " else " == "
    
    paste0(tabInput, " && ", tabInput, operator, "'", x, "'")
  }
  
  conditionalContent <- function(name, ...) {
    conditionalPanel(
      condition = makeCondition(name),
      ...
    )
  }
  
  tagList(
    tbl,
    fluidRow(
      column_2(expressionTabUI_sidebar(id, makeCondition)),
      column_10(
        conditionalContent("result", expressionResultsTabUI_main(ns("result"))),
        conditionalContent("volcano", expressionVolcanoTabUI_main(ns("volcano"))),
        conditionalContent("heatmap", expressionHeatmapTabUI_main(ns("heatmap"))),
        conditionalContent("GSEA", expressionGseaTabUI_main(ns("gsea"))),
        conditionalContent("gene clustering", XIFF::geneExpressionDimRedTabUI_main(ns("dimred")))
      )
    )
  )
  
}

expressionTabUI_sidebar <- function(id, makeCondition){
  ns <- NS(id)
  
  ret <- tagList(
    conditionalPanel(
      condition = makeCondition("GSEA", TRUE),
      checkboxInput(
        inputId = ns("result_filter"), 
        label = "filter for high expressors",
        value = FALSE
      ),
      conditionalPanel(
        condition = makeCondition("result", TRUE),
        checkboxInput(
          inputId = ns("useDT"), 
          label = "Use Filtering/Sorting of Results", 
          value = FALSE
        )
      ),
      selectInput(
        inputId = ns("geneset"), 
        label = "filter Hallmark genes", 
        choices = NULL,
        selectize = FALSE
      ),
      hr()
    ),
    conditionalPanel(
      condition = makeCondition("result"),
      expressionResultsTabUI_sidebar(ns("result"))
    ),
    conditionalPanel(
      condition = makeCondition("volcano"),
      expressionVolcanoTabUI_sidebar(ns("volcano"))
    ),
    conditionalPanel(
      condition = makeCondition("heatmap"),
      expressionHeatmapTabUI_sidebar(ns("heatmap"))
    ),
    conditionalPanel(
      condition = makeCondition("GSEA"),
      expressionGseaTabUI_sidebar(ns("gsea"))
    ),
    conditionalPanel(
      condition = makeCondition("gene clustering"),
      XIFF::geneExpressionDimRedTabUI_sidebar(ns("dimred"))
    )
  )
  
  div(
    uiOutput(ns("sidebar_run")),
    hr(),
    ret,
    uiOutput(ns("indicator"))
  )
}

expressionTab <- function(input, output, session, fm, classSelection, classLabel, gene_anno, 
                          TissueAnnotationFocus, gsea_data_hallmark, species, msigDBLink){
  # Sidebar -------------------------------------------------------------------
  ns <- session$ns
  
  output$indicator <- renderUI({
    hallmark_names <- names(gsea_data_hallmark())
    
    updateSelectInput(
      session = session, 
      inputId = "geneset",
      choices = c(
        "show all genes" = "all", 
        getHallmarkGeneSetChoices(sort(hallmark_names))
      ),
      selected = isolate(input$geneset) %||% "all"
    )
    
    NULL
  })
  
  output$sidebar_run <- renderUI({
    fmButtonBothClasses(ns("run"), fm, classSelection)
  })
  
  # FutureManager -------------------------------------------------------------
  Results <- reactiveVal()
  Args <- reactive({
    list(
      cs = reactiveValuesToList(classSelection),
      anno = gene_anno(),
      TissueAnnotationFocus = TissueAnnotationFocus(),
      gsea_data_hallmark = gsea_data_hallmark()
    )
  })
  
  fm$registerRunObserver(
    inputId = ns("run"),
    label = "Gene expression",
    statusVar = Results,
    longFun = expressionLongFun,
    Args = Args
  )
  
  resultClassLabel <- registerFreezedClassLabel(output, classLabel, Results, fm, ns("run"))
  
  TableData <- reactive({
    geneSet <- input$geneset
    filter_low <- input$result_filter
    validate(
      need(geneSet, "prefilter missing"),
      need(is.logical(filter_low), "prefilter missing")
    )
    res <- Results()
    fmValidate(res)
    r <- res$df
    
    if (filter_low){
      r <- r[res$correctTpms, ]
    }

    if (geneSet != "all") {
      gene_set_ensg <- res$gsea_data_hallmark[[geneSet]]
      r <- r %>% filter(ensg %in% gene_set_ensg)
    }
    
    cl <- reactiveValuesToList(resultClassLabel)
    
    r %>% 
      mutate_at(c("logFC", "P.Value", "adj.p.val"), signif, 3) %>% 
      mutate(higher = styleHigherCol(logFC < 0, cl)) %>%
      relocate(higher, .after = location)
  })
  
  FilteredData <- reactive({
    rows <- AllRows()
    tab <- TableData()
    
    if (!is.null(rows) && input$useDT){
      tab <- tab %>% slice(rows)
    }
    
    tab
  })
  
  # Result subtab -------------------------------------------------------------
  AllRows <- callModule(
    module = expressionResultsTab,
    id = "result",
    classSelection = classSelection,
    classLabel = resultClassLabel,
    Results = Results,
    TableData = TableData,
    species = species
  )
  
  # Volcano subtab ------------------------------------------------------------
  callModule(
    module = expressionVolcanoTab,
    id = "volcano",
    classLabel = resultClassLabel,
    TableData = FilteredData
  )
  
  # Heatmap subtab ------------------------------------------------------------
  callModule(
    module = expressionHeatmapTab,
    id = "heatmap",
    classLabel = resultClassLabel,
    Results = Results,
    TableData = FilteredData,
    AllRows = AllRows,
    TissueAnnotationFocus = TissueAnnotationFocus
  )
  
  # GSEA subtab ---------------------------------------------------------------
  callModule(
    module = expressionGseaTab,
    id = "gsea",
    fm = fm,
    classLabel = resultClassLabel,
    Results = Results,
    species = species,
    msigDBLink = msigDBLink
  )
  
  callModule(
    module = XIFF::geneExpressionDimRedTab,
    id = "dimred",
    fm = fm,
    Results = Results,
    TableData = FilteredData,
    classLabel = resultClassLabel,
    getDataGeneExpression = getTissueDataGeneExpression
  )
}

expressionLongFun <- function(task, cs, anno, TissueAnnotationFocus, gsea_data_hallmark){
  limmaData <- differentialGeneExpression_LimmaVoom(cs, anno, p = task)
  if (is.null(limmaData)) return() # handle the task cancel
  if (is.fmError(limmaData)) return(limmaData)
  df <- limmaData %>% 
    select(ensg, symbol, name, location, logFC, P.Value, adj.p.val)
  
  list(
    df = df,
    cs = cs,
    TissueAnnotationFocus = TissueAnnotationFocus,
    gsea_data_hallmark = gsea_data_hallmark,
    correctTpms = limmaData$correctTpms
  )
}
