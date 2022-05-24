library(testthat)
library(shinytest)

source("commonFunctions.R")

# Parameters ------------------------------------------------------------------
checkInterval = 500 # 500ms
timeout <- 25000 # 25s

# Functions -------------------------------------------------------------------
runClick <- function(app, id, cancel = FALSE){
  app$executeScript(sprintf("$('#%s').click();", id))
  app$waitForValue(
    name = id, 
    iotype = "input",
    ignore = list(cancel, NULL),
    timeout = timeout,
    checkInterval = checkInterval
  )
}

goToTab <- function(app, tabsetId, tabId, mainTab = "Analysis"){
  # set tab on navbar page
  app$setValue("app-MainMenu", value = mainTab)
  app$waitForValue(name = "app-MainMenu")  
  
  app$setValue(tabsetId, value = tabId)
  app$waitForValue(name = tabsetId)
}

startProcess <- function(app, tabId, buttonId, navMenu = "Analysis",
                         mainTabsetId = "app-AnalysisMenu",
                         subTabsetId = NULL, subtabId = NULL){
  goToTab(app, mainTabsetId, tabId, navMenu)
  if (!is.null(subtabId)){
    goToTab(app, subTabsetId, subtabId, navMenu)
  }
  
  app$waitForValue(
    name = buttonId,
    iotype = "input",
    timeout = timeout,
    checkInterval = checkInterval
  )
  runClick(app, buttonId)
  app$snapshot(filename = paste0(buttonId, ".json"), screenshot = TRUE)
}

testResults <- function(app, tabId, prefix, leftPlot = TRUE, middlePlot = TRUE, rightPlot = TRUE,
                        navMenu = "Analysis", mainTabsetId = "app-AnalysisMenu",
                        subTabsetId = NULL, subtabId = NULL){
  goToTab(app, mainTabsetId, tabId, navMenu)
  if (!is.null(subtabId)){
    goToTab(app, subTabsetId, subtabId, navMenu)
  }
  
  containerId <- paste0(prefix, "-plot_left-container")
  
  invisible(app$waitFor(
    expr = paste0("$('#", containerId, "').find('img').length;"),
    timeout = 45000, # 45 s
    checkInterval = checkInterval
  ))
  
  appState <- app$getAllValues(input = FALSE)
  
  if (leftPlot) testPlot(app, appState, paste0(prefix, "-plot_left-plot.*"))
  if (middlePlot) testPlot(app, appState, paste0(prefix, "-plot_middle-plot.*"))
  if (rightPlot) testPlot(app, appState, paste0(prefix, "-plot_right-plot.*"))
  testTable(app, appState, paste0(prefix, "-table"))
  app$snapshot(filename = paste0(prefix, ".json"), screenshot = TRUE)
}

testPlot <- function(app, appState, plotIdPattern){
  outputNames <- names(appState$output)
  plotId <- grep(plotIdPattern, outputNames, value = TRUE)
  
  isTooltip <- grepl("plot_tooltip", plotId)
  plotItem <- appState$output[[plotId]]
  
  if (isTooltip){
    itemNames <- c("html", "deps")
    itemType <- "svg"
  } else {
    itemNames <- c("src", "width", "height", "alt", "coordmap")
    itemType <- "img"
  }
  
  expect_named(plotItem, itemNames)
  expect_is(app$findElement(paste0("#", plotId, " > ", itemType)), "Element")
}

testTable <- function(app, appState, tableId){
  tableItem <- appState$output[[tableId]]
  
  expect_is(tableItem, "json")
  rows <- app$findElements(paste0("#", tableId, " tbody tr"))
  expect_true(length(rows) > 0)
}

# Go Go Go! -------------------------------------------------------------------
# Start the app
test_that(
  desc = "Scenario #1",
  code = {
    path <- system.file("shiny", package = "TIFF")
    app <- ShinyDriver$new(path, loadTimeout = 45000)
    
    initAppTestDir(
      app, 
      dirpath = Sys.getenv("TIFF_TESTS_ARTIFACTS_TABS", tempdir())
    )
    
    invisible(app$waitFor(
      expr = "$('#app-input-tissue_anno-handler-table-table').find('table').length;",
      timeout = timeout,
      checkInterval = checkInterval
    ))
    
    # Go to Mutation Burden input selection
    app$setInputs(`app-input-selection` = "mutational_burden", wait_ = FALSE, values_ = FALSE)
    invisible(app$waitForValue(
      name = "app-input-mut_burden-barplot-plot",
      iotype = "output",
      timeout = timeout,
      checkInterval = checkInterval
    ))
    
    # Limit tumortypes
    app$setInputs(
      `app-input-prefilter-tumor_selector` = c(
        "breast invasive carcinoma", 
        "liver hepatocellular carcinoma"
      ), 
      values_ = FALSE
    )
    
    # Select tissues for class 1
    invisible(app$waitForValue(
      name = "app-input-mut_burden-barplot-brush_options",
      iotype = "input",
      timeout = timeout,
      checkInterval = checkInterval
    ))
    app$setInputs(`app-input-mut_burden-barplot-brush_options` = "click", values_ = FALSE)
    app$setInputs(`app-input-mut_burden-barplot-score_method_y` = "number", values_ = FALSE)
    invisible(app$waitForValue(
      name = "app-input-mut_burden-barplot-number_action_y",
      iotype = "input",
      timeout = timeout,
      checkInterval = checkInterval
    ))
    app$setInputs(`app-input-mut_burden-barplot-brush_button` = "click", values_ = FALSE)
    statInfo <- app$waitForValue(
      name = "app-input-mut_burden-barplot-selectionStat",
      iotype = "output",
      timeout = timeout,
      checkInterval = checkInterval,
      ignore = list(NULL, "-")
    )
    expect_match(statInfo, "^10 tissues selected")
    
    app$setInputs(`app-input-classDetailsWrapper-class1-add` = "click", values_ = FALSE)
    
    # Select tissues for class 2
    app$setInputs(`app-input-mut_burden-barplot-brush_options` = "click", values_ = FALSE)
    app$setInputs(`app-input-mut_burden-barplot-number_y` = 15, values_ = FALSE)
    app$setInputs(`app-input-mut_burden-barplot-number_action_y` = "highest", values_ = FALSE)
    
    app$setInputs(`app-input-mut_burden-barplot-brush_button` = "click", values_ = FALSE)
    statInfo <- app$waitForValue(
      name = "app-input-mut_burden-barplot-selectionStat",
      iotype = "output",
      timeout = timeout,
      checkInterval = checkInterval,
      ignore = list(NULL, "-")
    )
    expect_match(statInfo, "^15 tissues selected")
    
    app$setInputs(`app-input-classDetailsWrapper-class2-add` = "click", values_ = FALSE)
    
    # Run some processes
    startProcess(
      app = app, 
      tabId = "copy number", 
      buttonId = "app-copy_number-run"
    )
    startProcess(
      app = app, 
      tabId = "gene expression", 
      subTabsetId = "app-expression-tabset",
      subtabId = "result",
      buttonId = "app-expression-run"
    )
    
    # And verify results
    # 1) no background process
    testResults(
      app = app, 
      tabId = "mutation", 
      subTabsetId = "app-mutation-tabset",
      subtabId = "per_gene",
      prefix = "app-mutation-per_gene-tab", 
      middlePlot = FALSE
    )
    testResults(
      app = app, 
      tabId = "mutation", 
      subTabsetId = "app-mutation-tabset",
      subtabId = "mutational_fraction",
      prefix = "app-mutation-fraction-tab"
    )
    testResults(
      app = app, 
      tabId = "protein expression",
      prefix = "app-protein-tab"
    )
    testResults(
      app = app, 
      tabId = "derived data", 
      subTabsetId = "app-derived_data-tabset",
      subtabId = "immune cells",
      prefix = "app-derived_data-immune_cells-tab"
    )
    testResults(
      app = app,
      tabId = "derived data", 
      subTabsetId = "app-derived_data-tabset",
      subtabId = "metabolics",
      prefix = "app-derived_data-metabolics-tab"
    )
    testResults(
      app = app,
      tabId = "derived data", 
      subTabsetId = "app-derived_data-tabset",
      subtabId = "signaling",
      prefix = "app-derived_data-signaling-tab",
      middlePlot = FALSE
    )
    testResults(
      app = app,
      tabId = "derived data", 
      subTabsetId = "app-derived_data-tabset",
      subtabId = "clones",
      prefix = "app-derived_data-clones-tab"
    )
    testResults(
      app = app,
      tabId = "additional properties", 
      subTabsetId = "app-additional_props-tabset",
      subtabId = "gene signatures",
      prefix = "app-additional_props-gene_signatures-tab"
    )
    testResults(
      app = app,
      tabId = "additional properties", 
      subTabsetId = "app-additional_props-tabset",
      subtabId = "hallmark gene set comparison",
      prefix = "app-additional_props-hallmark-tab"
    )
    
    # 2) with background process
    testResults(
      app = app, 
      tabId = "copy number", 
      prefix = "app-copy_number-tab"
    )
    testResults(
      app = app, 
      tabId = "gene expression", 
      subTabsetId = "app-expression-tabset",
      subtabId = "result",
      prefix = "app-expression-result-tab"
    )
    
    cleanupTestDir(app)
  }
)
