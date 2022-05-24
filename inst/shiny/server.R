shinyServer(function(input, output, session) {
  status <- initialDbCheck() # check for each new session, will return a modal window with the error
  if (!status){
    errorScreen(attr(status, "reason"))
    stop("Initial check failed!")
  }
  
  fm <- FutureManager::FutureManager$new(
    input = input,
    session = session,
    opts = c(
      "dbname", "dbhost", "dbport", "dbuser", "dbpass",
      "xiff.column", "xiff.label", "xiff.name", "xiff.schema"
    ),
    keepPreviousResults = TRUE
  )
  
  callModule(
    module = app,
    id = "app",
    fm = fm,
    settings = settings,
    geneSignatures = geneSignatures
  )
})
