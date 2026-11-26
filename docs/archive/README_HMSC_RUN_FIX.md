# HMSC run fix

This build fixes the HMSC run-time error:

`Error in data.frame: arguments imply differing number of rows: 1, 0`

Cause:
Empty standard comparison tables were created with one constant column and several zero-length columns.

Fix:
All empty standard tables now use zero-length columns consistently.

Also, if HMSC data validation fails, the workflow now creates a diagnostic ZIP containing:
- used_config.yml
- tables/data_summary.csv
- diagnostics/data_check_messages.csv
- diagnostics/engine_status.json
- results/README_HMSC_CHECK_FAILED.txt
