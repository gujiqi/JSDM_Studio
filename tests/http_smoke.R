args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args) >= 1) as.integer(args[[1]]) else 7837L
if (!is.finite(port)) port <- 7837L
shiny::runApp(".", port = port, host = "127.0.0.1", launch.browser = FALSE)
