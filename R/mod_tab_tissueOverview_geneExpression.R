tissueOverviewGeneExpressionTabUI_main <- function(id){
  ns <- NS(id)
  
  div(
    br(),
    dimensionReductionUI(ns("dimRed"))
  )
}

tissueOverviewGeneExpressionTabUI_sidebar <- function(id, mode = "gene_set"){
  ns <- NS(id)
  
  cluster_methods <- c(
    "pam k = automatic", "pam k = 2", "pam k = 3", "pam k = 4", 
    "affinity propagation",
    "k-means k = automatic", "k-means k = 2", "k-means k = 3", "k-means k = 4"
  )
  
  retVar <- if (mode == "gene_set"){
    selectInput(
      inputId = ns("gene_set"), 
      label = "Hallmark gene set:", 
      choices = NULL,
      selectize = FALSE
    )
  } else {
    sliderInput(
      inputId = ns("ngenes"), 
      label = "number of most varying genes", 
      min = 10, 
      max = 1000, 
      step = 10, 
      value = 300
    )
  }
  
  ret <- div(
    radioButtons(
      inputId = ns("dimRedMethod"), 
      label = "method:",
      choices = XIFF::dimRedAvailableMethods(), 
      selected = "tsne"
    ),
    hr(),
    retVar,
    hr(),
    selectInput(
      inputId = ns("labelDimRed"), 
      label = "label", 
      choices = NULL,
      selectize = FALSE
    ),
    selectInput(
      inputId = ns("colorDimRed"), 
      label = "colour", 
      choices = NULL,
      selectize = FALSE
    ),
    sliderInput(
      inputId = ns("fontsizeOverview"), 
      label = "font size", 
      min = 8, 
      max = 30, 
      step = 1, 
      value = 15
    ),
    hr(),
    selectInput(
      inputId = ns("cluster_method"), 
      label = "Select Cluster Method", 
      choices = cluster_methods,
      selected = cluster_methods[[1]]
    ),
    hr(),
    textInput(
      inputId = ns("filename"), 
      label = "Choose file name", 
      value = "filename"
    ),
    downloadButton(
      outputId = ns("download"), 
      label = "Download", 
      style = "padding:4px; margin:3px; font-size:90%"
    ),
    hr(),
    actionButton(
      inputId = ns("save_properties"), 
      label = "Save tissue properties",
      style = "padding:4px; margin:3px; font-size:90%"
    )
  )
  
  div(
    uiOutput(ns("sidebar_run")),
    hr(),
    ret,
    uiOutput(ns("indicator"))
  )
}

tissueOverviewGeneExpressionTab <- function(input, output, session, fm, classSelection, classLabel, classStack, mode, label,
                                            TissueAnnotation, TissueAnnotationFocus, gsea_data_hallmark){
  # Sidebar -------------------------------------------------------------------
  ns <- session$ns
  output$sidebar_run <- renderUI({
    fmButtonOneClass(ns("run"), fm, classSelection)
  })
  
  output$indicator <- renderUI({
    ti_anno <- TissueAnnotationFocus()
    gsea_data <- gsea_data_hallmark()
    
    if (mode == "gene_set"){
      updateSelectInput(
        session = session,
        inputId = "gene_set",
        choices = getHallmarkGeneSetChoices(sort(names(gsea_data))),
        selected = isolate(input$gene_set) %||% "HALLMARK_APOPTOSIS"
      )
    }
    
    updateSelectInput(
      session = session,
      inputId = "labelDimRed",
      choices = c("none", "cluster", colnames(ti_anno)),
      selected = isolate(input$labelDimRed) %||% "none"
    )
    
    updateSelectInput(
      session = session,
      inputId = "colorDimRed",
      choices = c("class selection", "cluster", setdiff(colnames(ti_anno), "tissuename")),
      selected = isolate(input$colorDimRed) %||% "class selection"
    )
    
    NULL
  })
  
  InputData <- reactiveVal()
  
  ArgsBasic <- reactive({
    list(
      cs = reactiveValuesToList(classSelection),
      ensg_gene_set = NA,
      TissueAnnotationFocus = TissueAnnotationFocus()
    )
  })
  
  ArgsSet <- reactive({
    validate(need(input$gene_set, "no gene set selected"))
    
    args <- ArgsBasic()
    args$ensg_gene_set <- input$gene_set
    args
  })
  
  Args <- if (mode == "gene_set"){
    ArgsSet
  } else {
    ArgsBasic
  }
  
  fm$registerRunObserver(
    inputId = ns("run"),
    label = label,
    statusVar = InputData,
    longFun = tissueOverviewGeneExpressionLongFun,
    Args = Args
  )
  
  resultClassLabel <- registerFreezedClassLabel(output, classLabel, Results, fm, ns("run"))
  
  AnalysisParams <- reactive({
    n <- if (mode == "varying_genes"){
      input$ngenes
    } else {
      50 # any number is ok
    }
    
    list(
      method = input$dimRedMethod,
      dataSource = mode,
      n = n,
      unit = "log2tpm"
    )
  })
  
  ClusterMethod <- reactive({
    input$cluster_method
  })
  
  PlotParams <- reactive({
    list(
      labelDimRed = input$labelDimRed,
      colorDimRed = input$colorDimRed,
      fontsizeOverview = input$fontsizeOverview
    )
  })
  
  Results <- callModule(
    module = dimensionReduction,
    id = "dimRed",
    InputData = InputData,
    AnalysisParams = AnalysisParams,
    ClusterMethod = ClusterMethod,
    PlotParams = PlotParams,
    classLabel = resultClassLabel
  )
  
  output$download <- downloadHandler (
    filename = function() { 
      paste0(input$filename, ".csv") 
    },
    content = function(file) {
      res <- Results()$d$data %>% 
        left_join(TissueAnnotation(), by = "tissuename") %>%
        select(tissuename, everything())
      
      readr::write_excel_csv(res, file, na = "")
    }
  )
  
  observeEvent(
    eventExpr = input$save_properties, 
    handlerExpr = {
      saveProperties(
        res = Results()$d$data, 
        classSelection = classSelection, 
        classLabel = resultClassLabel,
        classStack = classStack, 
        Annotation = TissueAnnotation, 
        msg = "dimension reduction data, cluster properties and other annotations"
      )
    }
  )
}

tissueOverviewGeneExpressionLongFun <- function(task, cs, ensg_gene_set, TissueAnnotationFocus){
  
  if(all(is.na(ensg_gene_set))) {
    ensgs <- ensg_gene_set
  } else {
    ensgs <- getGSEAdata("human", gene_set = ensg_gene_set)
  }
  df <- getExpressionDimRedData(cs, ensgs, p = task)
  if (is.null(df)) return() # handle the task cancel
  if (is.fmError(df)) return(df)
  
  list(
    df = df,
    cs = cs,
    annotationFocus = TissueAnnotationFocus
  )
}
