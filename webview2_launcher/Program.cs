using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Text.RegularExpressions;
using System.Runtime.InteropServices;
using System.Threading;
using System.Windows.Forms;
using Microsoft.Web.WebView2.Core;
using Microsoft.Web.WebView2.WinForms;

namespace JSDMStudioLauncher
{


internal static class NativeDpi
{
    // Fallback for Windows 8.1+ when Application.SetHighDpiMode is not enough under older runtimes.
    [DllImport("shcore.dll")]
    private static extern int SetProcessDpiAwareness(int value);

    [DllImport("user32.dll")]
    private static extern bool SetProcessDPIAware();

    public static void Enable()
    {
        try
        {
            // 2 = PROCESS_PER_MONITOR_DPI_AWARE
            SetProcessDpiAwareness(2);
        }
        catch
        {
            try { SetProcessDPIAware(); } catch { }
        }
    }
}

internal static class Program
{
    private static Process? rProcess;
    private static int Port = 0;
    private static readonly string Host = "127.0.0.1";
    private static string? globalLogFile;

    [STAThread]
    static void Main()
    {
        // High-DPI clarity fix: prevent Windows from bitmap-scaling the whole app.
        // This makes text and the embedded WebView2 much sharper on 125%/150%/175% displays.
        NativeDpi.Enable();
        try { Application.SetHighDpiMode(HighDpiMode.PerMonitorV2); } catch { }
        Application.EnableVisualStyles();
        Application.SetCompatibleTextRenderingDefault(false);

        Application.ThreadException += (sender, e) =>
        {
            WriteLog("UI exception: " + e.Exception.ToString());
            MessageBox.Show(e.Exception.ToString(), "JSDM Studio UI error", MessageBoxButtons.OK, MessageBoxIcon.Error);
        };

        AppDomain.CurrentDomain.UnhandledException += (sender, e) =>
        {
            WriteLog("Unhandled exception: " + (e.ExceptionObject?.ToString() ?? "unknown"));
        };

        try
        {
            string appDir = FindAppDirectory();
            Directory.SetCurrentDirectory(appDir);

            globalLogFile = Path.Combine(appDir, "startup_log.txt");
            WriteLog("Launcher started");
            WriteLog("App directory: " + appDir);

            Port = FindFreePort();
            WriteLog("Selected free port: " + Port);

            string? rscript = FindRscript();
            if (rscript == null)
            {
                MessageBox.Show(
                    "JSDM Studio needs R for Windows.\n\nPlease install R first:\nhttps://cran.r-project.org/bin/windows/base/\n\nAfter installing R, restart JSDM Studio.",
                    "R was not found",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Error);
                OpenUrl("https://cran.r-project.org/bin/windows/base/");
                return;
            }

            WriteLog("Rscript: " + rscript);

            string backendScript = Path.Combine(appDir, "run_shiny_backend.R");
            if (!File.Exists(backendScript))
            {
                MessageBox.Show($"Cannot find backend script:\n{backendScript}", "Missing file", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            StartRBackend(rscript, backendScript, appDir);

            bool ready = WaitForPort(Host, Port, TimeSpan.FromMinutes(15));
            if (!ready)
            {
                string msg =
                    "R/Shiny backend did not start in time.\n\nFirst launch can take several minutes if R packages are being installed.\n\n" +
                    "Please try these steps:\n" +
                    "1. Wait a few minutes and try again if this was the first launch\n2. Run install_packages.bat\n" +
                    "2. Check startup_log.txt in the installation folder\n" +
                    "3. Try Launch_in_default_browser.bat as fallback\n\n" +
                    $"Log file:\n{globalLogFile}";
                MessageBox.Show(msg, "Shiny backend did not start", MessageBoxButtons.OK, MessageBoxIcon.Error);
                TryKillR();
                return;
            }

            WriteLog("Shiny backend is ready. Opening WebView2 window.");
            Application.Run(new MainWindow($"http://{Host}:{Port}", globalLogFile));
        }
        catch (Exception ex)
        {
            WriteLog("Launcher fatal error: " + ex.ToString());
            MessageBox.Show(ex.ToString(), "JSDM Studio Launcher Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            TryKillR();
        }
    }

    private static int FindFreePort()
    {
        // Ask Windows for a free ephemeral port. This avoids fixed-port conflicts such as 3838 already in use.
        var listener = new TcpListener(IPAddress.Loopback, 0);
        listener.Start();
        int port = ((IPEndPoint)listener.LocalEndpoint).Port;
        listener.Stop();
        Thread.Sleep(100);
        return port;
    }

    public static void WriteLog(string text)
    {
        try
        {
            string log = globalLogFile ?? Path.Combine(AppContext.BaseDirectory, "startup_log.txt");
            File.AppendAllText(log, $"{DateTime.Now:yyyy-MM-dd HH:mm:ss} - {text}\n", Encoding.UTF8);
        }
        catch { }
    }

    private static string FindAppDirectory()
    {
        string baseDir = AppContext.BaseDirectory;
        if (File.Exists(Path.Combine(baseDir, "app.R"))) return baseDir;

        DirectoryInfo? dir = new DirectoryInfo(baseDir);
        for (int i = 0; i < 8 && dir != null; i++, dir = dir.Parent)
        {
            if (File.Exists(Path.Combine(dir.FullName, "app.R"))) return dir.FullName;
        }

        return baseDir;
    }

    private static string? FindRscript()
    {
        var candidates = new List<(string Path, Version Version, int Priority)>();

        foreach (var root in new[]
        {
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles), "R"),
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFilesX86), "R"),
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "Programs", "R")
        })
        {
            try
            {
                if (!Directory.Exists(root)) continue;
                foreach (var dir in Directory.GetDirectories(root, "R-*"))
                {
                    string candidate = Path.Combine(dir, "bin", "Rscript.exe");
                    if (File.Exists(candidate))
                    {
                        candidates.Add((candidate, ParseRVersion(Path.GetFileName(dir)), 1));
                    }
                }
            }
            catch { }
        }

        string? pathEnv = Environment.GetEnvironmentVariable("PATH");
        if (!string.IsNullOrWhiteSpace(pathEnv))
        {
            foreach (var p in pathEnv.Split(Path.PathSeparator))
            {
                try
                {
                    string candidate = Path.Combine(p.Trim(), "Rscript.exe");
                    if (File.Exists(candidate))
                    {
                        candidates.Add((candidate, new Version(0, 0, 0), 2));
                    }
                }
                catch { }
            }
        }

        foreach (var item in candidates.OrderBy(x => x.Priority).ThenByDescending(x => x.Version))
        {
            if (IsUsableRscript(item.Path)) return item.Path;
        }

        return null;
    }

    private static Version ParseRVersion(string folder)
    {
        var m = Regex.Match(folder, @"R-(\d+)\.(\d+)\.?(\d+)?");
        if (!m.Success) return new Version(0, 0, 0);
        int major = int.Parse(m.Groups[1].Value);
        int minor = int.Parse(m.Groups[2].Value);
        int patch = m.Groups[3].Success ? int.Parse(m.Groups[3].Value) : 0;
        return new Version(major, minor, patch);
    }

    private static bool IsUsableRscript(string path)
    {
        try
        {
            var psi = new ProcessStartInfo
            {
                FileName = path,
                Arguments = "--version",
                UseShellExecute = false,
                CreateNoWindow = true,
                RedirectStandardError = true,
                RedirectStandardOutput = true
            };
            using var p = Process.Start(psi);
            if (p == null) return false;
            p.WaitForExit(8000);
            string outText = p.StandardOutput.ReadToEnd();
            string errText = p.StandardError.ReadToEnd();
            return p.ExitCode == 0 || errText.Contains("R scripting front-end", StringComparison.OrdinalIgnoreCase) || outText.Contains("R scripting front-end", StringComparison.OrdinalIgnoreCase);
        }
        catch
        {
            return false;
        }
    }

    private static void StartRBackend(string rscript, string backendScript, string appDir)
    {
        var psi = new ProcessStartInfo
        {
            FileName = rscript,
            Arguments = $"\"{backendScript}\" {Port}",
            WorkingDirectory = appDir,
            UseShellExecute = false,
            CreateNoWindow = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true
        };

        rProcess = new Process { StartInfo = psi, EnableRaisingEvents = true };

        rProcess.OutputDataReceived += (_, e) =>
        {
            if (e.Data != null) WriteLog("R OUT: " + e.Data);
        };
        rProcess.ErrorDataReceived += (_, e) =>
        {
            if (e.Data != null) WriteLog("R ERR: " + e.Data);
        };
        rProcess.Exited += (_, e) =>
        {
            try
            {
                WriteLog("R backend exited. ExitCode=" + rProcess?.ExitCode);
            }
            catch
            {
                WriteLog("R backend exited.");
            }
        };

        WriteLog("Starting R backend on port " + Port);
        rProcess.Start();
        rProcess.BeginOutputReadLine();
        rProcess.BeginErrorReadLine();
    }

    private static bool WaitForPort(string host, int port, TimeSpan timeout)
    {
        DateTime end = DateTime.Now.Add(timeout);
        while (DateTime.Now < end)
        {
            if (!IsRRunning())
            {
                WriteLog("R process exited before port became ready.");
                return false;
            }

            try
            {
                using var client = new TcpClient();
                var result = client.BeginConnect(host, port, null, null);
                bool success = result.AsyncWaitHandle.WaitOne(TimeSpan.FromMilliseconds(700));
                if (success && client.Connected) return true;
            }
            catch { }

            Application.DoEvents();
            Thread.Sleep(500);
        }
        return false;
    }

    public static bool IsRRunning()
    {
        try
        {
            return rProcess != null && !rProcess.HasExited;
        }
        catch
        {
            return false;
        }
    }

    public static void TryKillR()
    {
        try
        {
            if (rProcess != null && !rProcess.HasExited) rProcess.Kill(true);
        }
        catch { }
    }

    private static void OpenUrl(string url)
    {
        try
        {
            Process.Start(new ProcessStartInfo { FileName = url, UseShellExecute = true });
        }
        catch { }
    }
}

public class MainWindow : Form
{
    private readonly WebView2 webView = new();
    private readonly string url;
    private readonly string logFile;
    private readonly System.Windows.Forms.Timer backendTimer = new();

    public MainWindow(string url, string logFile)
    {
        this.url = url;
        this.logFile = logFile;

        Text = "JSDM Studio";
        Width = 1500;
        Height = 950;
        MinimumSize = new Size(1100, 720);
        StartPosition = FormStartPosition.CenterScreen;
        AutoScaleMode = AutoScaleMode.Dpi;
        Font = new Font("Microsoft YaHei UI", 9F, FontStyle.Regular, GraphicsUnit.Point);

        var menu = new MenuStrip();
        var file = new ToolStripMenuItem("File");
        file.DropDownItems.Add(new ToolStripMenuItem("Reload", null, (_, _) => SafeReload()));
        file.DropDownItems.Add(new ToolStripMenuItem("Open startup log", null, (_, _) => OpenLog()));
        file.DropDownItems.Add(new ToolStripSeparator());
        file.DropDownItems.Add(new ToolStripMenuItem("Exit", null, (_, _) => Close()));

        var help = new ToolStripMenuItem("Help");
        help.DropDownItems.Add(new ToolStripMenuItem("About", null, (_, _) =>
            MessageBox.Show("JSDM Studio\nWebView2 window + R/Shiny/Hmsc backend", "About", MessageBoxButtons.OK, MessageBoxIcon.Information)));

        menu.Items.Add(file);
        menu.Items.Add(help);
        MainMenuStrip = menu;
        Controls.Add(menu);
        menu.Dock = DockStyle.Top;

        webView.Dock = DockStyle.Fill;
        Controls.Add(webView);
        webView.BringToFront();

        backendTimer.Interval = 3000;
        backendTimer.Tick += (_, _) =>
        {
            if (!Program.IsRRunning())
            {
                backendTimer.Stop();
                Program.WriteLog("Backend monitor detected that R is no longer running.");
                MessageBox.Show(
                    "The R/Shiny backend has stopped.\n\nPlease open startup_log.txt to see the error.\nYou can also try Launch_in_default_browser.bat as fallback.",
                    "R backend stopped",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Warning);
            }
        };
    }

    protected override async void OnLoad(EventArgs e)
    {
        base.OnLoad(e);

        try
        {
            await webView.EnsureCoreWebView2Async();
            webView.ZoomFactor = 1.0;

            webView.CoreWebView2.ProcessFailed += (_, ev) =>
            {
                Program.WriteLog("WebView2 process failed: " + ev.ProcessFailedKind.ToString());
                MessageBox.Show(
                    "WebView2 rendering process failed.\n\nPlease reopen JSDM Studio.\nIf it keeps happening, use Launch_in_default_browser.bat as fallback.",
                    "WebView2 process failed",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Error);
            };

            webView.CoreWebView2.NavigationCompleted += (_, ev) =>
            {
                if (!ev.IsSuccess)
                {
                    Program.WriteLog("Navigation failed. Error=" + ev.WebErrorStatus.ToString());
                }
            };

            Program.WriteLog("Navigating WebView2 to " + url);
            webView.CoreWebView2.Navigate(url);
            backendTimer.Start();
        }
        catch (Exception ex)
        {
            Program.WriteLog("WebView2 startup error: " + ex.ToString());
            MessageBox.Show(
                "WebView2 failed to start.\n\nPlease install Microsoft Edge WebView2 Runtime.\n\n" + ex.Message,
                "WebView2 Error",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error);
            try
            {
                Process.Start(new ProcessStartInfo
                {
                    FileName = "https://developer.microsoft.com/microsoft-edge/webview2/",
                    UseShellExecute = true
                });
            }
            catch { }
        }
    }

    private void SafeReload()
    {
        try
        {
            if (webView.CoreWebView2 != null) webView.Reload();
        }
        catch (Exception ex)
        {
            Program.WriteLog("Reload error: " + ex.ToString());
        }
    }

    protected override void OnFormClosed(FormClosedEventArgs e)
    {
        backendTimer.Stop();
        base.OnFormClosed(e);
        Program.TryKillR();
    }

    private void OpenLog()
    {
        try
        {
            if (File.Exists(logFile))
            {
                Process.Start(new ProcessStartInfo { FileName = logFile, UseShellExecute = true });
            }
            else
            {
                MessageBox.Show("startup_log.txt not found yet.", "Log", MessageBoxButtons.OK, MessageBoxIcon.Information);
            }
        }
        catch (Exception ex)
        {
            MessageBox.Show(ex.Message, "Could not open log", MessageBoxButtons.OK, MessageBoxIcon.Error);
        }
    }
}

}
