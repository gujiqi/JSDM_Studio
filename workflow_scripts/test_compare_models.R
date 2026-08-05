source("app.R", local = TRUE)

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("shiny is required for Compare Models test", call. = FALSE)
}

example_dir <- normalizePath(
  file.path("output", "HmscHPC_real_01_fixed_poisson_categorical_20260530_141309"),
  winslash = "/",
  mustWork = TRUE
)

shiny::testServer(server, {
  session$setInputs(compare_engines = c("Hmsc-HPC"))
  session$setInputs(compare_hmschpc_dir = example_dir)
  session$setInputs(compare_run = 1)
  stopifnot(!is.null(rv$compare))
  stopifnot("Status_class" %in% names(rv$compare))
  stopifnot("Standard_diagnostics" %in% names(rv$compare))
  stopifnot(isTRUE(rv$compare$Fitted_for_metric_comparison[1]))
  invisible(output$compare_table)
  invisible(output$compare_outputs_matrix)
  invisible(output$compare_notes)
})

cat("COMPARE_TEST_OK\n")
