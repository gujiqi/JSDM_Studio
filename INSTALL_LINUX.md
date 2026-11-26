# Install JSDM Studio on Linux

The guided online installer targets **Ubuntu 24.04 LTS x86_64** with AVX.
It is not an offline or system-wide application installer. Extract the
ZIP to a writable directory and enter the extracted `JSDMStudio` folder.

```bash
bash Install_JSDM_Studio_Linux.sh --check
bash Install_JSDM_Studio_Linux.sh --install-system
bash Launch_in_browser.sh
```

`--check` reports missing Ubuntu packages without installing anything.
`--install-system` authorizes the apt step if needed; only that step
uses sudo. If all system packages are present, omit that option. The installer
downloads R and CPU Python packages into an isolated per-user runtime
under `${XDG_DATA_HOME:-$HOME/.local/share}/JSDMStudio`. It preserves
the prior runtime until the new dependency checks pass. The system
step includes fonts for all interface scripts, including CJK and Bengali.

If PyPI downloads are slow or blocked, you can set
`JSDM_STUDIO_PYPI_INDEX` to a trusted HTTPS PyPI mirror for this command
only. This does not change any global pip configuration or the PyTorch
CPU wheel source. For example:

```bash
JSDM_STUDIO_PYPI_INDEX=https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple \
  bash Install_JSDM_Studio_Linux.sh --install-system
```

Follow the printed localhost URL if the browser does not open. Leave the
terminal running; use Ctrl+C to stop Studio. Results are stored in the
project's `output/` directory. Installation logs are
`install-linux-*.log` in the runtime location above.

For other Linux distributions, install R, the required packages and JAGS
using instructions for your distribution, then use `Rscript install_packages.R`
and `bash Launch_in_browser.sh`. Engine-specific Python requirements are listed
in `installer/platforms/requirements-linux-cpu.txt`; use an isolated compatible
environment and set `RETICULATE_PYTHON` and `JSDM_HMSCHPC_PYTHON`.
