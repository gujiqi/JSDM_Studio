(function () {
  "use strict";

  const scriptUrl = new URL(document.currentScript.src);
  const roots = ".navbar, .tab-content > .tab-pane, #studio-settings, #studio-toolbar";
  const excluded = "script, style, svg, [translate=no], [data-studio-no-i18n], table, pre, code, textarea, select, .selectize-control";
  const textAttrs = ["title", "aria-label", "placeholder", "data-title"];
  const nativeKeys = [
    "File", "Project", "Open example files", "Open results folder", "Open separate browser session", "Exit",
    "Analysis", "Engine Guide", "Model workflows", "Universal Benchmark", "Compare Models",
    "View", "Home", "Reload workspace...", "Zoom in", "Zoom out", "Reset zoom", "Settings", "Help",
    "Manuals and examples", "Open startup log", "Check for updates...", "Release notes", "Send feedback...",
    "About JSDM Studio"
  ];
  const extraKeys = [
    "+2 detection data", "+6 objective", "direct objective", "alternative objective", "response match",
    "spatial input", "new-site prediction", "detection data", "traits", "phylogeny", "Choose", "pts",
    "Preferences could not be saved. Changes apply to this session only.",
    "The language file could not be loaded. The previous language is still active.",
    "Toggle navigation", "Open workflow", "Close", "Browse...", "No file selected"
  ];
  const defaults = { locale: "en", textSize: "standard", inlineHelp: "translated" };
  const cache = new Map([["en", Object.create(null)]]);
  const originals = new WeakMap();
  const originalAttrs = new WeakMap();
  const originalOptions = new WeakMap();
  let options, ui, prefs = { ...defaults }, messages = cache.get("en");
  let pending = prefs;
  let revision = 0, request = 0, userChanged = false, queued = false;
  const tidy = (text) => String(text).replace(/\s+/gu, " ").trim();
  const lookup = (text) => Object.hasOwn(messages, text) ? messages[text] : text;

  function formattedLabel(text) {
    // Keep the invariant portion of compound labels and translate only the UI phrase.
    if (!text) return text;
    if (/^\+\d+ /.test(text)) {
      return text.split("; ").map((part) => {
        const score = /^(\+\d+ )(.+)$/.exec(part);
        return score ? score[1] + lookup(score[2]) : part;
      }).join("; ");
    }
    const compact = /^(.*?\d+) pts(\)?)$/.exec(text);
    if (compact) return compact[1] + " " + lookup("pts") + compact[2];
    const model = /^Choose ((?:int|svc|lfMs|spMs|ms|sp)?PGOcc)$/.exec(text);
    if (model) return lookup("Choose") + " " + model[1];
    const limit = /^(Min|Max) ([-\d.,]+)$/.exec(text);
    if (limit) return lookup(limit[1]) + " " + limit[2];
    const heading = /^([1-9]\. )(.+)$/.exec(text);
    if (heading) {
      const label = lookup(heading[2]);
      return heading[1] + (label === heading[2] ? formattedLabel(heading[2]) : label);
    }
    const engineLabel = /^(Hmsc-HPC|HMSC|Hmsc|jSDM|GJAM|spOccupancy|sjSDM|boral) (.+)$/.exec(text);
    if (engineLabel) {
      const suffix = engineLabel[2];
      const label = messages[suffix] ?? messages[suffix[0].toUpperCase() + suffix.slice(1)];
      if (label) return engineLabel[1] + " · " + label;
    }
    const command = /^(Check|Run|Download) (HMSC|Hmsc|Hmsc-HPC|jSDM|GJAM|spOccupancy|sjSDM|boral) (data|workflow|ZIP)$/.exec(text);
    if (command) {
      const phrase = command[1] + " " + command[3];
      if (messages[phrase]) return messages[phrase] + " · " + command[2];
    }
    const file = /^([A-Za-z_0-9.]+\.(?:csv|nwk|yml))\s*-\s*(.+)$/.exec(text);
    if (file && messages[file[2]]) return file[1] + " · " + messages[file[2]];
    return text;
  }

  function translation(text, originalEnglish = false) {
    if (originalEnglish) return text;
    const key = tidy(text);
    const result = lookup(key);
    const rendered = key === result ? formattedLabel(key) : result;
    if (rendered === key) return text;
    return text.match(/^\s*/u)[0] + rendered + text.match(/\s*$/u)[0];
  }

  function translatableRoots() {
    // The toolbar lives in the navbar when JavaScript is enabled.
    return Array.from(document.querySelectorAll(roots))
      .filter((node) => !node.parentElement?.closest(roots));
  }

  function staticSelect(select) {
    return !select.closest("table, [translate=no], [data-studio-no-i18n]") &&
      !/(?:_file_preview_choice|^univ_builtin_case)$/.test(select.id);
  }

  function scan(root, textNode, attribute, choice) {
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode(node) {
        return !tidy(node.data) || node.parentElement?.closest(excluded) ?
          NodeFilter.FILTER_REJECT : NodeFilter.FILTER_ACCEPT;
      }
    });
    for (let node = walker.nextNode(); node; node = walker.nextNode()) textNode(node);

    const selector = textAttrs.map((name) => `[${name}]`).join(",");
    for (const el of root.querySelectorAll(selector)) {
      if (el.closest("svg, table, pre, [data-studio-no-i18n]")) continue;
      for (const name of textAttrs) if (el.hasAttribute(name)) attribute(el, name);
    }

    for (const select of root.querySelectorAll("select")) {
      if (!staticSelect(select)) continue;
      const selectize = select.selectize;
      if (selectize) {
        for (const [key, option] of Object.entries(selectize.options)) {
          choice(option, selectize.settings.labelField, selectize, key);
        }
      } else if (!select.classList.contains("selectized") &&
          !root.querySelector(`script[data-for="${CSS.escape(select.id)}"]`)) {
        for (const option of select.options) if (option.hasAttribute("value")) choice(option, "textContent");
      }
    }
  }

  function localize() {
    for (const root of translatableRoots()) scan(root, (node) => {
      let saved = originals.get(node);
      if (!saved || node.data !== saved.displayed) {
        saved = { source: node.data, displayed: node.data };
        originals.set(node, saved);
      }
      const englishHelp = prefs.inlineHelp === "english" && node.parentElement.closest(".param-help");
      const value = translation(saved.source, englishHelp);
      if (node.data !== value) node.data = value;
      saved.displayed = value;
    }, (el, name) => {
      let map = originalAttrs.get(el);
      if (!map) originalAttrs.set(el, map = new Map());
      const current = el.getAttribute(name);
      let saved = map.get(name);
      if (!saved || current !== saved.displayed) map.set(name, saved = { source: current, displayed: current });
      const value = translation(saved.source);
      if (value !== current) el.setAttribute(name, value);
      saved.displayed = value;
    }, (option, field, selectize, key) => {
      if (typeof option[field] !== "string") return;
      const source = originalOptions.get(option) ?? option[field];
      originalOptions.set(option, source);
      const label = translation(source);
      if (label === option[field]) return;
      if (!selectize) {
        option[field] = label;
        return;
      }
      const next = { ...option, [field]: label };
      originalOptions.set(next, source);
      selectize.updateOption(key, next);
      originalOptions.set(selectize.options[key], source);
    });
  }

  function queueLocalization() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      localize();
    });
  }

  function mirror() {
    ui.language.value = prefs.locale;
    ui.dialogLanguage.value = prefs.locale;
    ui.size.value = prefs.textSize;
    ui.help.value = prefs.inlineHelp;
  }

  function status(message) {
    ui.status.textContent = message || "";
  }

  function persist() {
    if (!window.Shiny?.shinyapp?.isConnected()) return;
    window.Shiny.setInputValue("studio_preference_change", { prefs, revision: ++revision }, { priority: "event" });
  }

  async function change(values, save) {
    const target = { ...pending, ...values };
    if (!options.has(target.locale)) return;
    pending = target;
    const currentRequest = ++request;
    try {
      if (!cache.has(target.locale)) {
        const url = new URL("studio-locale-" + target.locale + ".json", scriptUrl);
        url.search = scriptUrl.search;
        const response = await fetch(url, { cache: "default" });
        if (!response.ok) throw new Error("Missing locale pack");
        const pack = await response.json();
        if (!pack || Array.isArray(pack) || typeof pack.Settings !== "string") throw new Error("Invalid locale pack");
        cache.set(target.locale, pack);
      }
    } catch {
      if (currentRequest !== request) return;
      pending = prefs;
      mirror();
      status("The language file could not be loaded. The previous language is still active.");
      return;
    }
    if (currentRequest !== request) return;
    messages = cache.get(target.locale);
    prefs = target;
    document.documentElement.lang = prefs.locale;
    document.documentElement.dir = prefs.locale === "ar" ? "rtl" : "ltr";
    document.documentElement.dataset.studioTextSize = prefs.textSize;
    document.documentElement.dataset.studioInlineHelp = prefs.inlineHelp;
    mirror();
    status("");
    localize();
    if (window.chrome?.webview?.postMessage) {
      window.chrome.webview.postMessage({
        type: "studio.locale", locale: prefs.locale,
        labels: Object.fromEntries(nativeKeys.map((key) => [key, lookup(key)]))
      });
    }
    if (save) persist();
  }

  function initialize() {
    const get = (id) => document.getElementById("studio-" + id);
    ui = {
      language: get("language"), dialogLanguage: get("settings-language"), settings: get("settings"),
      size: get("settings-size"), help: get("settings-help"), status: get("settings-status")
    };
    if (!ui.language || !ui.settings) return;
    options = new Set(Array.from(ui.language.options, (option) => option.value));

    const nav = document.querySelector(".navbar > .container-fluid");
    if (nav) {
      const marker = nav.querySelector(":scope > .navbar-collapse");
      nav.insertBefore(get("toolbar"), marker);
      nav.classList.add("studio-has-tools");
    }
    const use = (key, value) => {
      userChanged = true;
      change({ [key]: value }, true);
    };
    ui.language.addEventListener("change", () => use("locale", ui.language.value));
    ui.dialogLanguage.addEventListener("change", () => use("locale", ui.dialogLanguage.value));
    ui.size.addEventListener("change", () => use("textSize", ui.size.value));
    ui.help.addEventListener("change", () => use("inlineHelp", ui.help.value));
    get("open-settings").addEventListener("click", () => ui.settings.showModal());
    ui.settings.querySelectorAll("[data-studio-dismiss]").forEach((button) => {
      button.addEventListener("click", () => ui.settings.close());
    });
    get("reset-settings").addEventListener("click", () => {
      userChanged = true;
      change(defaults, true);
    });
    window.studioOpenSettings = () => {
      if (!ui.settings.open) ui.settings.showModal();
    };

    new MutationObserver((records) => {
      if (records.some(({ target }) => {
        const el = target.nodeType === Node.ELEMENT_NODE ? target : target.parentElement;
        return el?.closest(roots) && !el.closest("table, pre, svg");
      })) queueLocalization();
    }).observe(document.body, { childList: true, characterData: true, subtree: true });

    const shiny = window.Shiny;
    shiny.addCustomMessageHandler("studio.preferences.load", (saved) => {
      if (userChanged) persist();
      else change(saved, false);
    });
    shiny.addCustomMessageHandler("studio.preferences.saved", (result) => {
      if (result.revision !== revision) return;
      status(result.ok ? "" : "Preferences could not be saved. Changes apply to this session only.");
    });
    const ask = () => shiny.setInputValue("studio_preference_request", true, { priority: "event" });
    window.jQuery(document).on("shiny:connected.studio-preferences", ask);
    if (shiny.shinyapp?.isConnected()) ask();
  }

  window.StudioLocalization = Object.freeze({
    inventory() {
      const strings = new Set([...nativeKeys, ...extraKeys]);
      const add = (text) => {
        const key = tidy(text);
        if (/[A-Za-z]{2}/.test(key)) strings.add(key);
      };
      for (const root of translatableRoots()) scan(root,
        (node) => add(originals.get(node)?.source ?? node.data),
        (el, name) => add(originalAttrs.get(el)?.get(name)?.source ?? el.getAttribute(name)),
        (option, field) => add(originalOptions.get(option) ?? option[field]));
      return Array.from(strings).sort();
    }
  });
  document.addEventListener("DOMContentLoaded", initialize);
})();
