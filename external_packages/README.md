# External Packages

This GitHub package only keeps the lightweight `hmsc-hpc-main` Python source needed by the Hmsc-HPC CPU workflow.

Large reference PDFs, archived source ZIPs and package documentation snapshots are intentionally excluded from this source package. They should be stored outside Git or attached to a GitHub Release if needed.

To prepare the Hmsc-HPC Python environment on Windows, run:

```bat
install_hmschpc_python_packages.bat
```

The other modelling engines are installed through R packages by:

```bat
install_packages.bat
```
