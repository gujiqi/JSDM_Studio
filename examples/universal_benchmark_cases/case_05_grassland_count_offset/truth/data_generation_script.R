# Rebuild the JSDM Studio universal benchmark data.
runner <- file.path('..', '..', 'workflow_scripts', 'universal_benchmark_runner.R')
system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'),
        c(runner, '--action=generate'))
