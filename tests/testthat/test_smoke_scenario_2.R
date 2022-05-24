library(testthat)
library(shinytest)

source("commonFunctions.R")

# Parameters ------------------------------------------------------------------
checkInterval = 500 # 500ms
timeout <- 25000 # 25s

goToInputAndTestObject <- function(app, inputId, outputId, ..., isPlot = TRUE){
  setValue <- function(...){
    app$setInputs(..., wait_ = FALSE, values_ = FALSE)
  }
  
  setValue(`app-input-selection` = inputId)
  if (length(list(...)) > 0){
    setValue(...)
  }
  
  obj <- app$waitForValue(
    name = outputId,
    iotype = "output",
    timeout = timeout,
    checkInterval = checkInterval
  )
  
  if (is.list(obj)){
    type <- obj[["type"]]
    if (is.character(type) && "validation" %in% type){
      # wait until everything is ready
      obj <- app$waitForValue(
        name = outputId,
        iotype = "output",
        timeout = timeout,
        checkInterval = checkInterval,
        ignore = list(obj)
      )
    }
  }
  
  if (isPlot){
    testInputPlot(app, obj, outputId)
  } else {
    testInputTable(app, obj, outputId)
  }
  
  app$snapshot(filename = paste0(outputId, ".json"), screenshot = TRUE)
}

testInputPlot <- function(app, plotItem, plotId){
  itemNames <- c("src", "width", "height", "alt", "coordmap")
  expect_named(plotItem, itemNames)
  
  invisible(app$waitFor(
    expr = paste0("$('#", plotId, " > img').length;"),
    timeout = timeout,
    checkInterval = checkInterval
  ))
  
  expect_is(app$findElement(paste0("#", plotId, " > img")), "Element")
}

testInputTable <- function(app, tableItem, tableId){
  expect_is(tableItem, "json")
  
  invisible(app$waitFor(
    expr = paste0("$('#", tableId, " tbody').length;"),
    timeout = timeout,
    checkInterval = checkInterval
  ))
  
  rows <- app$findElements(paste0("#", tableId, " tbody tr"))
  expect_true(length(rows) > 0)
}

test_that(
  desc = "Scenario #2",
  code = {
    path <- system.file("shiny", package = "TIFF")
    app <- ShinyDriver$new(path, loadTimeout = 45000)
    
    initAppTestDir(
      app, 
      dirpath = Sys.getenv("TIFF_TESTS_ARTIFACTS_INPUTS", tempdir())
    )
    
    invisible(app$waitFor(
      expr = "$('#app-input-tissue_anno-handler-table-table').find('table').length;",
      timeout = timeout,
      checkInterval = checkInterval
    ))
    
    # Patient Annotation
    goToInputAndTestObject(
      app = app, 
      inputId = "patient_annotation", 
      outputId = "app-input-patient_anno-handler-table-table", 
      isPlot = FALSE
    )
    
    # Mutational Burden
    goToInputAndTestObject(
      app = app, 
      inputId = "mutational_burden", 
      outputId = "app-input-mut_burden-barplot-plot"
    )
    
    # Mutation Status
    goToInputAndTestObject(
      app = app, 
      inputId = "mutation",
      outputId = "app-input-mutation-table",
      `app-input-mutation-symbol` = "KRAS",
      isPlot = FALSE
    )
    
    # Copy Number Alteration 
    goToInputAndTestObject(
      app = app, 
      inputId = "copy_number",
      outputId = "app-input-cn-barplot-plot",
      `app-input-cn-symbol` = "KRAS"
    )
    
    goToInputAndTestObject(
      app = app, 
      inputId = "copy_number",
      outputId = "app-input-cn-barplot-plot",
      `app-input-cn-plot_type` = "types"
    )
    
    # Gene Expression
    goToInputAndTestObject(
      app = app, 
      inputId = "gene_expression",
      outputId = "app-input-gene_expr-barplot-plot",
      `app-input-gene_expr-symbol` = "KRAS"
    )
    
    # Protein Expr
    goToInputAndTestObject(
      app = app, 
      inputId = "protein_expression", 
      outputId = "app-input-protein_expr-barplot-plot"
    )

    # Immune Cells
    goToInputAndTestObject(
      app = app, 
      inputId = "immune_cells",
      outputId = "app-input-immune_cells-barplot-plot"
    )

    # Metabolics
    goToInputAndTestObject(
      app = app, 
      inputId = "metabolics", 
      outputId = "app-input-metabolics-barplot-plot"
    )
    
    # Signaling
    goToInputAndTestObject(
      app = app, 
      inputId = "signaling", 
      outputId = "app-input-signaling-barplot-plot"
    )
    
    # Clones
    goToInputAndTestObject(
      app = app, 
      inputId = "clones", 
      outputId = "app-input-clones-barplot-plot"
    )
    
    # Gene Signatures
    goToInputAndTestObject(
      app = app, 
      inputId = "gene_signature", 
      outputId = "app-input-gene_sig-barplot-plot"
    )
    
    # Hallmark Sets
    goToInputAndTestObject(
      app = app, 
      inputId = "hallmark_set", 
      outputId = "app-input-hallmark-barplot-plot"
    )
    
    # Multiple Features
    goToInputAndTestObject(
      app = app, 
      inputId = "multiple_features", 
      outputId = "app-input-multifeatures-plot-plot"
    )
    
    cleanupTestDir(app)
  }
)
