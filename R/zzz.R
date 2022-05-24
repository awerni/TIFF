#' @rdname setTIFF
#' @export
TIFFOptions <- list(
  "xiff.column" = "tissuename",
  "xiff.label" = "tissue",
  "xiff.name" = "tissues",
  "xiff.schema" = "tissue",
  "xiff.tooltipCallbackFun" = XIFF:::tooltipCallbackFun
)

#' Set TIFF parameters.
#'
#' @param code to be evaluated.
#' 
#' @details 
#' 
#' \code{setTIFF} sets the options globally.
#' 
#' \code{withTIFF} evaluates a block of code within \code{TIFF} namespace
#' and with \code{TIFF} options.
#' 
#' \code{TIFFOptions} list with \code{TIFF} options.
#'
#' @return No value. Used for side-effects.
#' @export
#' 
#' @importFrom withr with_namespace
#'
#' @rdname setTIFF
#' @examples
#' 
#' # ops is used to revert options at the end of the example
#' ops <- options("xiff.column" = "test") 
#' 
#' getOption("xiff.column") # returns 'test'
#' 
#' 
#' withTIFF({
#'   x <- getOption("xiff.column") # "tissuename"
#'   print(x)
#' })
#' 
#' oldOpts <- setTIFF() # set TIFF options globally
#' 
#' # TIFF code to evaluate
#' getOption("xiff.column") # "tissuename"
#' 
#' options(oldOpts) # revert old options
#' 
#' getOption("xiff.column") # returns 'test'
#' options(ops) # revert options
#' 
#' setTIFF()
#' 
withTIFF <- function(code) {
  old <- options(TIFFOptions)
  on.exit(options(old))
  withr::with_namespace("TIFF", code)
}

#' @export
#' @rdname setTIFF
setTIFF <- function() {
  options(TIFFOptions)
}

#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import XIFF
#' @import FutureManager
#' @import shiny
.onLoad <- function(libname, pkgname){
  setTIFF()
  
  options("bitmapType" = "cairo")
  options("dplyr.summarise.inform" = FALSE)
}
