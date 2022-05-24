#' Run the TIFF application
#'
#' @param ... additional params passed into \code{shiny::runApp function}.
#'
#' @return nothing, called for side effects
#' @rdname run
#' @export
run <- function(...){
  appPath <- system.file("shiny", package = "TIFF")
  runApp(appPath, ...)
}
