# Install JSDM Studio on macOS

Choose the kit for your Mac from **Apple menu > About This Mac**. Use the
`arm64` kit for an Apple chip and the `x86_64` kit for Intel. The CPU
runtime bundle targets macOS 14 or later because the bundled R package
fs uses that deployment target.

1. Extract the ZIP in Finder.
2. Open `Install_JSDM_Studio.command` as your regular user, without
   prefixing the entire command with `sudo`. The script uses the included R,
   Python and JAGS installers when needed. macOS may ask for an
   administrator password. The R and Python analysis packages install
   from the files in the kit. It verifies the payload before asking for
   elevation, and activates the new runtime only after dependency checks.
3. When Terminal reports **Installation complete**, open
   `Launch_JSDM_Studio.command` in the extracted kit or
   `~/Applications/JSDM Studio.command`. The app opens in your default browser.

To start from Terminal, run `bash Install_JSDM_Studio.command` from the
extracted folder, then `bash Launch_JSDM_Studio.command`.

The application, logs and analysis outputs live under
`~/Library/Application Support/JSDMStudio`. Existing `app/input`,
`app/output` and `app/config.yml` are retained. Keep the Terminal window
running during analysis. Close it or press Ctrl+C there to stop the app.
If installation stops, open the newest `install-macos-*.log` in that
folder. For a first installation, allow at least 10 GB of free space;
an upgrade retains the prior runtime. If Finder refuses to open a
downloaded script, verify that the ZIP matches its release checksum and
use the documented Terminal command, subject to local security policy.
These ZIPs provide a browser launcher, not a signed/notarized native app.
