port <- as.integer(Sys.getenv("JSDM_STUDIO_SMOKE_PORT", "7799"))
shiny::runApp(appDir = ".", host = "127.0.0.1", port = port, launch.browser = FALSE)

