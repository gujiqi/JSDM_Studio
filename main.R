# main.R
# Command-line HMSC workflow driven by config.yml.
# Most users can use the Shiny GUI in app.R instead.

library(yaml)
library(Hmsc)
library(coda)
library(ape)
source("R/helpers.R", encoding = "UTF-8")

cfg <- yaml::read_yaml("config.yml")
outdir <- make_output_dir("cli")
cat("Output folder:", outdir, "\n")

yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))

Y <- safe_read_csv(cfg$input$Y, "Y.csv")
XData <- safe_read_csv(cfg$input$XData, "XData.csv")
TrData <- NULL
if (!is.null(cfg$input$TrData) && cfg$input$TrData != "") TrData <- safe_read_csv(cfg$input$TrData, "TrData.csv")
studyDesign <- NULL
if (!is.null(cfg$input$studyDesign) && cfg$input$studyDesign != "") studyDesign <- safe_read_csv(cfg$input$studyDesign, "studyDesign.csv")
coordinates <- NULL
if (!is.null(cfg$input$coordinates) && cfg$input$coordinates != "") coordinates <- safe_read_csv(cfg$input$coordinates, "coordinates.csv")
phyloTree <- NULL
if (!is.null(cfg$input$phyloTree) && cfg$input$phyloTree != "") phyloTree <- ape::read.tree(cfg$input$phyloTree)

log_file <- file.path(outdir, "log.txt")
log_fun <- function(txt) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", txt)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

res <- fit_hmsc_workflow(Y, XData, TrData, phyloTree, studyDesign, coordinates, cfg, outdir, log_fun)
cat("Done. Results saved in:", outdir, "\n")
