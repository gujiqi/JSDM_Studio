# Hmsc-HPC workflow audit

Generated: 2026-06-01 00:12:39

The suite runs synthetic Hmsc-HPC cases through the JSDM Studio R adapter and Python runner.
Expected successful statuses are `fitted` for sampler cases and `model_defined` for the guarded random-slope compile-only case.

## Case summary

                                  case        status zip_size_kb file_count
          01_fixed_poisson_categorical        fitted       117.4         87
                02_fixed_probit_binary        fitted       102.3         87
            03_fixed_normal_hmc_toggle        fitted        46.3         84
                   04_traits_phylo_cov        fitted       110.7         91
               05_iid_random_intercept        fitted       114.8         98
      06_spatial_full_random_intercept        fitted       110.3        100
                 07_phylo_newick_fixed        fitted       104.3         90
          08_random_slope_compile_only model_defined        28.7         55
 09_gaussian_alias_tfd_predictions_off        fitted       127.2         82
 missing_required empty_dirs
                            
                            
                            
                            
                            
                            
                            
                            
                            

## Diagnostics

For any failed case, open `diagnostics/engine_status.json` first, then `diagnostics/HmscHPC_error.txt`, `diagnostics/HmscHPC_compile.log`, `diagnostics/HmscHPC_validate_init.log` and `diagnostics/HmscHPC_sample.log` in that output folder.

Machine summary CSV: C:/Users/Google/Documents/CodeX/JSDMStudio_BIG_WITH_HMSCHPC_WORK_20260530/JSDMStudio/output/HmscHPC_real_example_suite_summary_20260601_001239.csv
