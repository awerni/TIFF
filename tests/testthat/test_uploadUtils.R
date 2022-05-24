library(testthat)

anno <- tibble::tibble(
  tissuename = c("ti_a", "ti.b", "ti c", "ti-d", "tie")
)

df <- tibble::tibble(
  tissuename = c("ti c", "ti-d", "tie", "tiz"),
  other = c("dummy", "dummy2", "tie", "dummy4"),
  other2 = c(1, 3, 1, 5),
  other3 = c(43.2, 12.53, 16.54, 37.4)
)

test_that(
  desc = "Tissuename column is detected properly",
  code = {
    df1 <- df
    expect_equal(
      object = TIFF:::findItemColumn(df1, anno),
      expected = "tissuename"
    )
    
    expect_null(TIFF:::findItemColumn(NULL, anno)) # NULL input
    expect_null(TIFF:::findItemColumn(list(), anno)) # no-data.frame object
    expect_null(TIFF:::findItemColumn(data.frame(), anno)) # data.frame with no rows
    expect_null(TIFF:::findItemColumn(df1[c("other2", "other3")], anno)) # only numeric columns
    
    expect_equal(
      object = TIFF:::findItemColumn(df1[c("other", "other2", "other3")], anno),
      expected = "other" # 1 common value (3rd row)
    ) 
    
    names(df1)[1] <- "dummy"
    expect_equal(
      object = TIFF:::findItemColumn(df1, anno),
      expected = "dummy"
    )
    
    df1 <- df1[c(2,4,1,3)]
    expect_equal(
      object = TIFF:::findItemColumn(df1, anno),
      expected = "dummy"
    )
  }
)

test_that(
  desc = "Tissue upload data is parsed properly",
  code = {
    df2 <- df
    
    res <- getTissuenameTranslation(df2, anno)
    expect_is(res, "list")
    expect_named(res, c("df", "missing"))
    expect_named(res[["df"]], c("tissuename", "other", "other2", "other3"))
    expect_equal(res[["missing"]], "tiz")
    
    df2 <- df[-4, ]
    res <- getTissuenameTranslation(df2, anno)
    expect_is(res, "list")
    expect_named(res, c("df", "missing"))
    expect_named(res[["df"]], c("tissuename", "other", "other2", "other3"))
    expect_length(res[["missing"]], 0)
    
    names(df2)[1] <- "dummy"
    res <- getTissuenameTranslation(df2, anno)
    expect_is(res, "list")
    expect_named(res, c("df", "missing"))
    expect_named(res[["df"]], c("tissuename", "other", "other2", "other3"))
    expect_length(res[["missing"]], 0)
    
    df2 <- df2[c(2,4,1,3)]
    res <- getTissuenameTranslation(df2, anno)
    expect_is(res, "list")
    expect_named(res, c("df", "missing"))
    expect_named(res[["df"]], c("tissuename", "other", "other3", "other2"))
    expect_length(res[["missing"]], 0)
    
    df2 <- df2[c("other", "other2", "other3")]
    res <- getTissuenameTranslation(df2, anno)
    expect_is(res, "list")
    expect_named(res, c("df", "missing"))
    expect_named(res[["df"]], c("tissuename", "other2", "other3"))
    expect_equal(res[["missing"]], c("dummy", "dummy2"))
  }
)
