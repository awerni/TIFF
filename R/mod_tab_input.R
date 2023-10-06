inputTabUI_main <- function(id){
  ns <- NS(id)
  
  makeCondition <- function(x){
    selectionId <- ns("selection")
    paste0("input['", selectionId, "'] == '", x, "'")
  }
  
  shinyjs::hidden(div(
    id = ns("panel"),
    conditionalPanel(
      condition = makeCondition("tissue_annotation"),
      tissueAnnotationInputModeUI(ns("tissue_anno"))
    ),
    conditionalPanel(
      condition = makeCondition("patient_annotation"),
      patientAnnotationInputModeUI(ns("patient_anno"))
    ),
    conditionalPanel(
      condition = makeCondition("mutational_burden"),
      mutationalBurdenInputModeUI(ns("mut_burden"))
    ),
    conditionalPanel(
      condition = makeCondition("mutation"),
      mutationInputModeUI(ns("mutation"))
    ),
    conditionalPanel(
      condition = makeCondition("copy_number"),
      copyNumberInputModeUI(ns("cn"))
    ),
    conditionalPanel(
      condition = makeCondition("gene_expression"),
      geneExpressionInputModeUI(ns("gene_expr"))
    ),
    conditionalPanel(
      condition = makeCondition("protein_expression"),
      proteinExpressionRppaInputModeUI(ns("protein_expr"))
    ),
    conditionalPanel(
      condition = makeCondition("immune_cells"),
      immuneCellsInputModeUI(ns("immune_cells"))
    ),
    conditionalPanel(
      condition = makeCondition("metabolics"),
      metabolicsInputModeUI(ns("metabolics"))
    ),
    conditionalPanel(
      condition = makeCondition("signaling"),
      signalingInputModeUI(ns("signaling"))
    ),
    conditionalPanel(
      condition = makeCondition("clones"),
      clonesInputModeUI(ns("clones"))
    ),
    conditionalPanel(
      condition = makeCondition("gene_signature"),
      geneSignatureInputModeUI(ns("gene_sig"))
    ),
    conditionalPanel(
      condition = makeCondition("hallmark_set"),
      hallmarkSetsInputModeUI(ns("hallmark"))
    ),
    conditionalPanel(
      condition = makeCondition("multiple_features"),
      multipleFeaturesInputModeUI(ns("multifeatures"))
    ),
    conditionalPanel(
      condition = makeCondition("file_upload"),
      uploadInputModeUI(ns("upload"), allowRds = TRUE)
    ),
    conditionalPanel(
      condition = makeCondition("restore_selection"),
      restoreSelectionInputModeUI(ns("restore"))
    ),
    conditionalPanel(
      condition = makeCondition("programmer_import"),
      programmerImportInputModeUI(ns("programmer"))
    )
  ))
}

inputTabUI_selector <- function(id) {
  
  ns <- NS(id)
  selectionId <- ns("selection")
  selector <- selectInput(
    inputId = ns("selection"), 
    label = "Select tissue by", 
    choices = c(
      "Tissue Annotation"      = "tissue_annotation",
      "File Upload"            = "file_upload",
      "Patient Annotation"     = "patient_annotation",
      "Mutational Burden"      = "mutational_burden",
      "Mutation Status"        = "mutation", 
      "Copy Number Alteration" = "copy_number", 
      "Gene Expression"        = "gene_expression",
      "Protein Expression"     = "protein_expression",
      "Immune Cells"           = "immune_cells",
      "Metabolic Pathways"     = "metabolics",
      "Signaling Pathways"     = "signaling",
      "Number of Clones"       = "clones",
      "Gene Signatures"        = "gene_signature",
      "Hallmark Sets"          = "hallmark_set",
      "Multiple Features"      = "multiple_features",
      "Restore Selection"      = "restore_selection",
      "Programmer's import"    = "programmer_import"
    ), 
    selected = "tissue_annotation", width = "80%"
  )
  
  tagList(div(selector, class = "input-method-selector"), hr())
}

#' @importFrom glue glue
#' @importFrom logger log_trace
inputTabUI_sidebar <- function(id){
  ns <- NS(id)
  selectionId <- ns("selection")
  
  makePrefilterCondition <- function(forbidden, allowed){
    forbidden <- paste0("input['", selectionId, "'] != '",
           forbidden, "'", collapse = " && ")
    allowed <- paste(glue::glue("output['{allowed}'] == 'show'"), collapse = " && ")
    cond <- glue::glue("({forbidden}) || ({allowed})")
    log_trace("Prefilter condition: {cond}")
    cond
  }
  
  list(
    conditionalPanel(
      condition = makePrefilterCondition(
        c("file_upload", "restore_selection", "programmer_import"),
        allowed = ns("show_prefilter_when_ml")),
      list(
        tissueSamplePrefilterUI(ns("prefilter")),
        hr()
      )
    ),
    classDetailsWrapperUI_main(
      ns("classDetailsWrapper"),
      c("class1", "class2")
    ),
    hr(),
    actionButton(
      inputId = ns("save"), 
      label = "Save selection", 
      style = "padding:4px; margin:3px; font-size:90%"
    ),
    actionButton(
      inputId = ns("discard_unbalanced_tt"), 
      label = "Discard unbalanced tumortypes", 
      style = 'padding:4px; margin:3px; font-size:90%'
    ),
    uiOutput(ns("indicator"))
  )
}

inputTab <- function(input, output, session, classSelection, classLabel, classStack, fm, species,
                     TissueAnnotation, TissueAnnotationFocus, PatientAnnotationFuller, 
                     geneSignatures, hallmarkGeneSets){
  output$indicator <- renderUI({
    shinyjs::show("panel")
    NULL
  })
  
  TissuePrefilter <- callModule(
    module = tissueSamplePrefilter, 
    id = "prefilter"
  )
  
  # Input modules -------------------------------------------------------------
  TissueAnnotation_ti <- callModule(
    module = tissueAnnotationInputMode,
    id = "tissue_anno",
    TissuePrefilter = TissuePrefilter, 
    TissueAnnotation = TissueAnnotation
  )
  
  PatientAnnotation_ti <- callModule(
    module = patientAnnotationInputMode,
    id = "patient_anno",
    TissuePrefilter = TissuePrefilter, 
    PatientAnnotationFuller = PatientAnnotationFuller
  )
  
  MutationalBurden_ti <- callModule(
    module = mutationalBurdenInputMode,
    id = "mut_burden",
    TissuePrefilter = TissuePrefilter
  )
  
  Mutation_ti <- callModule(
    module = mutationInputMode,
    id = "mutation",
    species = species, 
    TissuePrefilter = TissuePrefilter
  )
  
  CopyNumber_ti <- callModule(
    module = copyNumberInputMode,
    id = "cn",
    species = species, 
    TissuePrefilter = TissuePrefilter
  )
  
  Expression_ti <- callModule(
    module = geneExpressionInputMode,
    id = "gene_expr",
    species = species, 
    TissuePrefilter = TissuePrefilter
  )
  
  ProteinExpression_ti <- callModule(
    module = proteinExpressionRppaInputMode,
    id = "protein_expr",
    TissuePrefilter = TissuePrefilter
  )
  
  ImmuneCell_ti <- callModule(
    module = immuneCellsInputMode,
    id = "immune_cells",
    TissuePrefilter = TissuePrefilter
  )
  
  Metabolics_ti <- callModule(
    module = metabolicsInputMode,
    id = "metabolics",
    TissuePrefilter = TissuePrefilter
  )
  
  Signaling_ti <- callModule(
    module = signalingInputMode,
    id = "signaling",
    TissuePrefilter = TissuePrefilter
  )
  
  Clones_ti <- callModule(
    module = clonesInputMode,
    id = "clones",
    TissuePrefilter = TissuePrefilter
  )
  
  GeneSignature_ti <- callModule(
    module = geneSignatureInputMode,
    id = "gene_sig",
    TissuePrefilter = TissuePrefilter,
    geneSignatures = geneSignatures
  )
  
  HallmarkSet_ti <- callModule(
    module = hallmarkSetsInputMode,
    id = "hallmark",
    TissuePrefilter = TissuePrefilter,
    geneSets = hallmarkGeneSets
  )
  
  MultipleFeatures_ti <- callModule(
    module = multipleFeaturesInputMode,
    id = "multifeatures",
    species = species,
    TissuePrefilter = TissuePrefilter
  )
  
  AllTumortype <- reactive({
    req(TissuePrefilter())
    TissuePrefilter()$tissue
  })
  
  MlAnnotationFiltered <- reactive({
    req(TissueAnnotation(), TissuePrefilter())
    prefilter <- TissuePrefilter()
    anno <- TissueAnnotation()
    anno[anno[[prefilter$db_col]] %in% prefilter$tissue,]
  })
  
  Upload_ti <- callModule(
    module = uploadInputMode,
    id = "upload",
    AnnotationFull = TissueAnnotation,
    translationFun = getTissuenameTranslation,
    mlUseTumortypeFilter = FALSE,
    AnnotationFiltered = MlAnnotationFiltered
  )
  
  output$show_prefilter_when_ml <- renderText({
    input$selection
    mlUploadStatus <- input[["upload-ml-show"]]
    
    res <-if(is.null(mlUploadStatus) || input$selection != "file_upload") {
      "hide"
    } else {
      "show"
    }
    log_trace("show_prefilter_when_ml is present - {res}")    
    res
  })
  outputOptions(output, "show_prefilter_when_ml", suspendWhenHidden = FALSE)
  
  Restore_ti <- callModule(
    module = restoreSelectionInputMode,
    id = "restore",
    classStack = classStack
  )
  
  Programmer_ti <- callModule(
    module = programmerImportInputMode,
    id = "programmer",
    Annotation = TissueAnnotation
  )
  
  # Selected tissue -----------------------------------------------------------
  SelectedTissue <- reactive({
    ti <- switch(
      EXPR = input$selection,
      tissue_annotation = TissueAnnotation_ti(),
      patient_annotation = PatientAnnotation_ti(),
      mutational_burden = MutationalBurden_ti(),
      mutation = Mutation_ti(),
      copy_number = CopyNumber_ti(),
      gene_expression = Expression_ti(),
      protein_expression = ProteinExpression_ti(),
      immune_cells = ImmuneCell_ti(),
      metabolics = Metabolics_ti(),
      signaling = Signaling_ti(),
      clones = Clones_ti(),
      gene_signature = GeneSignature_ti(),
      hallmark_set = HallmarkSet_ti(),
      multiple_features = MultipleFeatures_ti(),
      file_upload = Upload_ti(),
      restore_selection = Restore_ti(),
      programmer_import = Programmer_ti()
    )
    
    req(ti)
    ti
  })
  
  # Class operations ----------------------------------------------------------
  checkTissueOverlap <- function(new_ti, orig_ti) {
    overlap_ti <- intersect(new_ti, orig_ti)
    if (length(overlap_ti) > 0) {
      new_ti <- setdiff(new_ti, overlap_ti)
      
      showModal(div(
        class = "danger",
        modalDialog(
          title = "Tissue overlap",
          tags$p("Already in the other class:"),
          tags$ul(lapply(overlap_ti, tags$li)),
          footer = modalButton("Ok")
        )
      ))
    }
    
    new_ti
  }
  
  onChange <- function(Selected, old, new, reset = TRUE, checkOverlap = FALSE, anotherSelection = NULL){
    if (reset){
      s <- SelectedTissue()$source
      if (!is.na(s)) session$resetBrush(s)
    }
    
    res <- if (checkOverlap){
      checkTissueOverlap(new, anotherSelection)
    } else {
      new
    }
    
    if (isDifferent(old, res)){
      fm$outdateRuns(TRUE)
    }
    
    res
  }
  
  callModule(
    module = classDetailsWrapper,
    id = "classDetailsWrapper",
    classLabel = classLabel,
    classSelection = classSelection,
    Selected = SelectedTissue,
    onChange = onChange
  )
  
  # Other actions -------------------------------------------------------------
  observeEvent(
    eventExpr = input$save,
    handlerExpr = {
      saveProperties(
        res = NULL, 
        classSelection = classSelection,
        classLabel = classLabel,
        classStack = classStack, 
        Annotation = TissueAnnotation, 
        msg = "tissue annotations"
      )
    }
  )
  
  observeEvent(
    eventExpr = input$discard_unbalanced_tt,
    handlerExpr = {
      shinyDropUnbalancedTumortypes(
        AnnotationFocus = TissueAnnotationFocus,
        classSelection = classSelection
      )
    }
  )
}
