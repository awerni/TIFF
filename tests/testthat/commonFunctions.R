initAppTestDir <- function(app, dirpath, .verbose = interactive()) {
  
  path <- file.path(dirpath, "tests/snapshot-current")
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # remove previous files from the test directory
  exFiles <- dir(path, full.names = TRUE)
  if(length(exFiles) > 0) {
    unlink(exFiles)
  }
  
  
  # The line below is a real insanity. However, when the shiny app
  # is not initialized using record tests, it does not sets the app
  # directory in a proper way. However it is required to make the 
  # snapshotDownload works correctly.
  app$.__enclos_env__$private$path <- dirpath
  
  if(.verbose) message("Output path: ", path)
  
  invisible(app)
}

cleanupTestDir <- function(app) {
  
  if(Sys.getenv("TIFF_KEEP_TESTS_ARTIFACTS", "FALSE") == "TRUE") {
    resultPath <- getAppFileFromTestDir(app, "")
    message("Path to files generated during tests: ", resultPath)
    return()
  }
  
  files <- dir(getAppFileFromTestDir(app, ""), full.names = TRUE)
  if(length(files) > 0) {
    unlink(files)
  }
  invisible()
}

getAppFileFromTestDir <- function(app, file) {
  
  file.path(app$.__enclos_env__$private$path,
            "tests/snapshot-current",
            file)
}
