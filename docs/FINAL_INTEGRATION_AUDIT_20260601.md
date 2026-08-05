# JSDM Studio Final Integration Audit - 2026-06-01

## Static and smoke tests
- APP_PARSE_OK
- APP_SOURCE_OK
- STRICT_AUDIT_OK ui_ids=519 input_refs=519 output_handlers=70 download_handlers=8
- COMPARE_TEST_OK
- HTTP_SMOKE_STATUS=200 CONTENT_LEN=641594 HTTP_SMOKE_CONTENT_OK
- COMPARE_MODEL_DEFINED_STATUS_OK

## Workflow suite command results
- Hmsc: exit=0, seconds=56.9, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\Hmsc.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\Hmsc.err.log
- Hmsc-HPC: exit=0, seconds=159.1, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\Hmsc-HPC.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\Hmsc-HPC.err.log
- jSDM: exit=0, seconds=9, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\jSDM.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\jSDM.err.log
- GJAM: exit=0, seconds=10.8, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\GJAM.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\GJAM.err.log
- spOccupancy: exit=0, seconds=15.6, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\spOccupancy.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\spOccupancy.err.log
- sjSDM: exit=0, seconds=213.4, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\sjSDM.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\sjSDM.err.log
- boral: exit=0, seconds=11.7, stdout=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\boral.out.log, stderr=C:\Users\Google\Documents\CodeX\JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530\JSDMStudio\output\FINAL_VALIDATION_20260601_000903\boral.err.log

## Workflow case status counts
- boral: fitted = 6
- boral: model_defined = 1
- GJAM: fitted = 7
- Hmsc: fitted = 8
- Hmsc-HPC: fitted = 8
- Hmsc-HPC: model_defined = 1
- jSDM: fitted = 7
- sjSDM: fitted = 9
- spOccupancy: fitted = 7

## Included latest summaries
- output/Hmsc_real_example_suite_summary_20260601_001000.csv
- output/HmscHPC_real_example_suite_summary_20260601_001239.csv
- output/jSDM_real_example_suite_summary_20260601_001248.csv
- output/GJAM_real_example_suite_summary_20260601_001259.csv
- output/spOccupancy_real_example_suite_summary_20260601_001314.csv
- output/sjSDM_real_example_suite_summary_20260601_001648.csv
- examples/boral/last_real_example_suite_summary.csv and output/boral_example_suite_20260601_001648/boral_suite_summary.csv

## Installer
- Rebuilt installer: installer/output/JSDMStudio_Setup.exe
- Installer size bytes: 142835620
- Installer LastWriteTime: 06/01/2026 00:18:45

## Remaining scientific/runtime risks
- Full production interpretation still requires adequate MCMC/optimizer settings, convergence checks and domain-specific validation beyond smoke-suite sizes.
- Hmsc-HPC is packaged as the CPU pyhmsc/HDF5 workflow, not a GPU/Slurm front end and not a full replacement for every R Hmsc feature.
- sjSDM depends on reticulate/PyTorch; boral depends on system JAGS; spOccupancy requires proper detection-nondetection replicate structure.
- Compare Models can compare statuses and standardized summaries, but raw parameters and association matrices are engine-specific and not directly interchangeable.
