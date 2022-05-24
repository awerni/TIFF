library(testthat)
library(TIFF)

if (!shinytest::dependenciesInstalled()){
  shinytest::installDependencies()
}

test_check("TIFF")
