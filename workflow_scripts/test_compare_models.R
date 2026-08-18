source("app.R", local = TRUE)

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("shiny is required for Compare Models test", call. = FALSE)
}

example_dir <- make_engine_run_dir("CompareModels_test_Hmsc", "CompareModelsTest")
Y <- matrix(c(0, 1, 1, 0, 1, 1), nrow = 3, ncol = 2,
            dimnames = list(paste0("site_", 1:3), paste0("sp_", 1:2)))
X <- data.frame(pH = c(6.1, 6.8, 7.2), moisture = c(0.2, 0.5, 0.7),
                row.names = rownames(Y), check.names = FALSE)
status <- list(engine = "Hmsc", status = "fitted", runtime_seconds = 0.1,
               warnings = character(), errors = character())
if (requireNamespace("yaml", quietly = TRUE)) {
  yaml::write_yaml(list(engine = "Hmsc", test = "compare_models"), file.path(example_dir, "used_config.yml"))
} else {
  writeLines(c("engine: Hmsc", "test: compare_models"), file.path(example_dir, "used_config.yml"))
}
write_data_check_messages(example_dir, "Synthetic Compare Models fixture passed data checks.")
write_engine_status(example_dir, status)
write_standard_outputs(example_dir, "Hmsc", status, Y, X)
write.csv(data.frame(engine = "Hmsc", response_id = "sp_1", predictor = "pH",
                     direction = "positive", estimate = 0.42, lower = 0.05,
                     upper = 0.81, notes = "test fixture", stringsAsFactors = FALSE),
          file.path(example_dir, "standard", "effects_long.csv"), row.names = FALSE)
write.csv(data.frame(engine = "Hmsc", site_id = rownames(Y), response_id = "sp_1",
                     observed = Y[, "sp_1"], predicted_mean = c(0.15, 0.74, 0.63),
                     predicted_lower = c(0.04, 0.52, 0.41), predicted_upper = c(0.39, 0.91, 0.83),
                     stringsAsFactors = FALSE),
          file.path(example_dir, "standard", "predictions_long.csv"), row.names = FALSE)
write.csv(data.frame(engine = "Hmsc", response_1 = "sp_1", response_2 = "sp_2",
                     association_type = "Omega_random_level", estimate = 0.12,
                     comparable_level = "qualitative", stringsAsFactors = FALSE),
          file.path(example_dir, "standard", "associations_long.csv"), row.names = FALSE)
write.csv(data.frame(Step = "fixture", Expected_file = "standard/run_summary.csv",
                     Meaning = "Compare Models test fixture", stringsAsFactors = FALSE),
          file.path(example_dir, "tables", "Hmsc_result_workflow_map.csv"), row.names = FALSE)
make_html_report(example_dir, "Hmsc", list(engine = "Hmsc"), "fitted")
ensure_output_contract(example_dir, "Hmsc", status, Y, X)
zip_path <- make_zip(example_dir)
stopifnot(file.exists(zip_path), file.info(zip_path)$size > 0)

shiny::testServer(server, {
  session$setInputs(compare_engines = c("Hmsc"))
  session$setInputs(compare_hmsc_dir = normalizePath(example_dir, winslash = "/", mustWork = TRUE))
  session$setInputs(compare_run = 1)
  stopifnot(!is.null(rv$compare))
  stopifnot("Status_class" %in% names(rv$compare))
  stopifnot("Engine_status_raw" %in% names(rv$compare))
  stopifnot("Standard_diagnostics" %in% names(rv$compare))
  stopifnot("Standard_effects_species_environment" %in% names(rv$compare))
  stopifnot("Standard_predictions_site_species" %in% names(rv$compare))
  stopifnot("Standard_associations_species_species" %in% names(rv$compare))
  stopifnot("Effects_direction_comparison_allowed" %in% names(rv$compare))
  stopifnot("Prediction_comparison_allowed" %in% names(rv$compare))
  stopifnot("Association_numeric_comparison_allowed" %in% names(rv$compare))
  stopifnot("Comparable_outputs" %in% names(rv$compare))
  stopifnot("Non_comparable_outputs" %in% names(rv$compare))
  stopifnot("Primary_diagnostics" %in% names(rv$compare))
  stopifnot(identical(rv$compare$Engine_status[1], "fitted"))
  stopifnot(isTRUE(rv$compare$Fitted_for_metric_comparison[1]))
  stopifnot(identical(rv$compare$Association_numeric_comparison_allowed[1], FALSE))
  invisible(output$compare_table)
  invisible(output$compare_outputs_matrix)
  invisible(output$compare_notes)
  invisible(output$comparison_meaning_table)
})

cat("COMPARE_TEST_OK\n")
cat("COMPARE_TEST_FIXTURE=", normalizePath(example_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("COMPARE_TEST_ZIP=", normalizePath(zip_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
