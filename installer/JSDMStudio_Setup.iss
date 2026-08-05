; JSDMStudio_Setup.iss
; Build with Inno Setup Compiler (ISCC.exe)
; This installer packages the Shiny/JSDM Studio and creates desktop/start-menu shortcuts.

#define MyAppName "JSDM Studio"
#define MyAppVersion "0.1.0"
#define MyAppPublisher "Jiqi Gu"
#define MyAppURL "https://github.com/"
#define MyAppExeName "JSDMStudioLauncher.exe"

[Setup]
AppId={{D98C4A60-40A9-4C7F-8E3F-JSDMSTUDIO01}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={localappdata}\JSDMStudio
DefaultGroupName={#MyAppName}
AllowNoIcons=yes
OutputDir=output
OutputBaseFilename=JSDMStudio_Setup
Compression=lzma
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
SetupIconFile=..\assets\JSDMStudio.ico
UninstallDisplayName={#MyAppName}
UninstallDisplayIcon={app}\assets\JSDMStudio.ico

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create desktop icon: JSDM Studio Window Mode"; GroupDescription: "Additional icons:"
Name: "browsericon"; Description: "Create desktop icon: JSDM Studio Browser Mode"; GroupDescription: "Additional icons:"


[Files]
Source: "..\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs; Excludes: "installer\output\*,installer\JSDMStudio_Setup.iss,installer\build_installer.bat,output\*,JSDMStudioLauncher.exe.WebView2\*,launcher_build\*,webview2_launcher\bin\*,webview2_launcher\obj\*,diagnostics\*,input\*,models\*,tables\*,plots\*,predictions\*,report\*,results\*,chains\*,missing_data\*,ordination\*,sensitivity\*,standard\*,*.log,*.err,*.out,*.tmp,*.pyc,*.pyo,tmp_*"

[Icons]
Name: "{group}\JSDM Studio - Window Mode"; Filename: "{app}\JSDMStudioLauncher.exe"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Comment: "Launch JSDM Studio in independent WebView2 window"
Name: "{group}\JSDM Studio - Browser Mode"; Filename: "{app}\Launch_in_default_browser.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Comment: "Launch JSDM Studio in default browser"
Name: "{group}\Repair R packages"; Filename: "{app}\Repair_packages.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Comment: "Install or repair required R packages"
Name: "{group}\Open startup log"; Filename: "{app}\Open_startup_log.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Comment: "Open startup log for troubleshooting"
Name: "{group}\Troubleshoot R"; Filename: "{app}\Troubleshoot_R.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Comment: "Check R installation"
Name: "{group}\Uninstall JSDM Studio"; Filename: "{uninstallexe}"
Name: "{autodesktop}\JSDM Studio"; Filename: "{app}\JSDMStudioLauncher.exe"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Tasks: desktopicon; Comment: "Launch JSDM Studio in independent WebView2 window"
Name: "{autodesktop}\JSDM Studio Browser Mode"; Filename: "{app}\Launch_in_default_browser.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\JSDMStudio.ico"; Tasks: browsericon; Comment: "Launch JSDM Studio in default browser"


[Run]
Filename: "{app}\JSDMStudioLauncher.exe"; Description: "Launch JSDM Studio Window Mode now"; Flags: postinstall skipifsilent unchecked


[Code]
function IsRInstalled(): Boolean;
var
  ResultCode: Integer;
begin
  Result := Exec('cmd.exe', '/C where Rscript.exe', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) and (ResultCode = 0);
  if not Result then begin
    Result := DirExists(ExpandConstant('{pf}\R'));
  end;
end;

function InitializeSetup(): Boolean;
begin
  Result := True;
  if not IsRInstalled() then begin
    MsgBox(
      'R for Windows was not detected.' + #13#10 + #13#10 +
      'JSDM Studio needs R and the Hmsc R package to run.' + #13#10 +
      'The installer can continue, but the first launch will ask the user to install R.',
      mbInformation, MB_OK);
  end;
end;
