tissueOverviewHallmarkTabUI_main <- function(id){
  ns <- NS(id)
  plotWrapperUI(
    id = ns("heatmap"), 
    height = "800px"
  )
}

tissueOverviewHallmarkTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  sortId <- ns("sort")
  makeSortCondition <- function(x){
    conditionId <- paste0("input['", sortId, "']")
    paste0(conditionId, " && ", conditionId, " == '", x, "'")
  }
  
  div(
    checkboxGroupInput(
      inputId = ns("show_anno"),
      label = NULL, 
      choices = c(
        "show tumortype" = "tumortype", 
        "show organ" = "organ", 
        "show tumortype adjacent" = "tumortype_adj"
      )
    ),
    radioButtons(
      inputId = sortId,
      label = "sort by:",
      choices = c("none", "class", "gene set"),
      selected = "none"
    ),
    conditionalPanel(
      condition = makeSortCondition("gene set"),
      selectInput(
        inputId = ns("gene_set"),
        label = "Hallmark gene set:",
        choices = NULL,
        selectize = FALSE
      )
    ),
    hr(),
    radioButtons(
      inputId = ns("scale"), 
      label = "heatmap scaling:", 
      choices = c(
        "none",
        "spread weight over tissues (columns)" = "column", 
        "spread weight over gene sets (rows)" = "row"
      ),
      selected = "column"
    ),
    conditionalPanel(
      condition = makeSortCondition("none"),
      radioButtons(
        inputId = ns("clustering_column"), 
        label = "clustering dist. column:", 
        choices = c(
          "euclidean distance" = "euclidean", 
          "correlation" = "correlation"
        ), 
        selected = "euclidean"
      )
    ),
    radioButtons(
      inputId = ns("clustering_row"), 
      label = "clustering dist. row:", 
      choices = c(
        "euclidean distance" = "euclidean", 
        "correlation" = "correlation"
      ),
      selected = "euclidean"
    ),
    uiOutput(ns("indicator"))
  )
}

tissueOverviewHallmarkTab <- function(input, output, session, classSelection, classLabel, 
                                      TissueAnnotationFocus, gsea_data_hallmark){
  # Sidebar --------------------------------------------------------------------
  output$indicator <- renderUI({
    gsea_data <- gsea_data_hallmark()
    
    updateSelectInput(
      session = session,
      inputId = "gene_set",
      choices = getHallmarkGeneSetChoices(sort(names(gsea_data))),
      selected = isolate(input$gene_set) %||% "HALLMARK_APOPTOSIS"
    )
    
    NULL
  })
  
  HeatmapData <- reactive({
    anno <- TissueAnnotationFocus()
    cs <- reactiveValuesToList(classSelection)
    
    validate(need(
      length(cs$class1) > 0 || length(cs$class2) > 0,
      "goto input tab and select tissues for any class"
    ))
    
    data <- getHeatmapDataHallmark(cs)
    validate(need(data, "not enough data"))
    
    anno <- TissueAnnotationFocus() %>% 
      select(tissuename, tumortype, tumortype_adjacent, organ)
    
    pheno <- stackClasses(cs) %>% 
      left_join(anno, by = "tissuename") %>% 
      column_to_rownames("tissuename")
    pheno <- pheno[colnames(data), ]
    
    list(
      data = data,
      pheno = pheno
    )
  })
  
  HallmarkPlotExpr <- reactive({
    res <- HeatmapData()
    req(res)
    
    annoToShow <- input$show_anno
    
    cl <- reactiveValuesToList(classLabel)
    
    generateHallmarkHeatmap(
      mat = as.matrix(res$data),
      pheno = res$pheno,
      cl = cl,
      scale = input$scale,
      sortMethod = input$sort,
      clustering_distance_rows = input$clustering_row,
      clustering_distance_cols = input$clustering_column,
      gene_set = input$gene_set,
      showTumortype = "tumortype" %in% annoToShow,
      showOrgan = "organ" %in% annoToShow,
      showTumortypeAdj = "tumortype_adj" %in% annoToShow
    )
  })
  
  callModule(
    module = plotWrapper,
    id = "heatmap",
    PlotExpr = HallmarkPlotExpr
  )
}

getHeatmapAnnotation <- function(pheno, cl, showTumortype = FALSE, 
                                 showOrgan = FALSE, showTumortypeAdj = FALSE){
  cl_1 <- cl$class1_name
  cl_2 <- cl$class2_name
  
  pheno <- droplevels(pheno)
  tmpRownames <- rownames(pheno)
  
  classColors <- c(plotColors[1], plotColors[2])
  names(classColors) <- c(cl_1, cl_2)
  phenoColors <- list(class = classColors)
  
  pheno <- pheno %>%
    mutate(class = factor(
      x = ifelse(class == "class1", cl_1, cl_2),
      levels = c(cl_1, cl_2)
    ))
  
  annoToDisplay <- "class"
  if (showTumortype) annoToDisplay <- c(annoToDisplay, "tumortype")
  if (showOrgan) annoToDisplay <- c(annoToDisplay, "organ")
  if (showTumortypeAdj) annoToDisplay <- c(annoToDisplay, "tumortype_adjacent")
  
  pheno <- pheno[, annoToDisplay, drop = FALSE]
  
  onlyNa <- sapply(pheno, function(x) all(is.na(x)))
  onlyNa <- names(which(onlyNa))
  
  if (length(onlyNa) > 0){
    pheno <- pheno[, setdiff(annoToDisplay, onlyNa), drop = FALSE]
  }
  
  sortedLegend <- sort(sapply(pheno, nlevels), decreasing = TRUE)
  toBeRemoved <- names(sortedLegend)[sortedLegend > 10]
  
  pheno <- pheno %>% 
    select_at(names(sortedLegend)) %>% # ensure correct order
    mutate_at(toBeRemoved, minifyText) # ensure minimal required legend width
  tmp <- lapply(
    X = pheno %>% select(-class),
    FUN = getScaleColors
  )
  phenoColors <- c(phenoColors, tmp)
  
  rownames(pheno) <- tmpRownames
  
  list(
    pheno = pheno,
    colors = phenoColors,
    toBeRemoved = toBeRemoved
  )
}

generateHallmarkHeatmap <- function(mat, pheno, cl, scale, sortMethod, clustering_distance_rows, 
                                    clustering_distance_cols, gene_set, showTumortype = FALSE, 
                                    showOrgan = FALSE, showTumortypeAdj = FALSE){
  anno <- getHeatmapAnnotation(
    pheno = pheno, 
    cl = cl, 
    showTumortype = showTumortype, 
    showOrgan = showOrgan, 
    showTumortypeAdj = showTumortypeAdj
  )
  
  if (scale != "none"){
    mat <- if (scale == "row"){
      scaleRows(mat)
    } else if (scale == "column"){
      t(scaleRows(t(mat)))
    }
  }
  
  cluster_cols <- TRUE
  cutree_cols <- NA
  
  if (sortMethod == "gene set"){
    cluster_cols <- FALSE
    ix <- sort(unlist(mat[gene_set,]), index.return = TRUE)$ix
    mat <- mat[, ix]
    
  } else if (sortMethod == "class" && n_distinct(pheno$class) == 2){
    cutree_cols <- 2
    clustering_distance_cols <- getClassDistances(mat, pheno)
  }
  
  rownames(mat) <- gsub("^HALLMARK_", "", rownames(mat))
  p <- getPheatmap(
    mat = mat, 
    annotation_col = anno$pheno,
    annotation_colors = anno$colors,
    scale = "none", 
    cluster_cols = cluster_cols,
    cutree_cols = cutree_cols,
    clustering_distance_rows = clustering_distance_rows, 
    clustering_distance_cols = clustering_distance_cols
  )
  
  removeAnnotationLegends(p, anno$toBeRemoved)
}

getScaleColors <- function(x){
  items <- levels(x)
  structure(
    ggColorHue(length(items)),
    names = items
  )
}

removeAnnotationLegends <- function(p, toBeRemoved){
  if (length(toBeRemoved) == 0) return(p)
  
  layout <- p$gtable$layout
  idx <- which(layout$name == "annotation_legend")
  
  if (length(idx) == 0) return(p)
  
  toBeRemoved <- unlist(lapply(toBeRemoved, paste0, c("", " r", " t")))
  legendGrobOrder <- p$gtable$grobs[[idx]]$childrenOrder
  idxToRemove <- names(legendGrobOrder) %in% toBeRemoved
  
  p$gtable$grobs[[idx]]$childrenOrder <- legendGrobOrder[!idxToRemove]
  p$gtable$grobs[[idx]]$children <- p$gtable$grobs[[idx]]$children[!idxToRemove]
  p
}

minifyText <- function(x){
  x <- as.character(x)
  uq <- sort(unique(x))
  uqMinified <- as.character(seq_along(uq))
  
  factor(
    x = uqMinified[match(x, uq)],
    levels = uqMinified
  )
}
