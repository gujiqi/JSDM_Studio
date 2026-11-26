(function () {
  "use strict";

  const pages = new Set([
    "Hmsc Workflow", "Hmsc-HPC Workflow", "jSDM Workflow", "GJAM Workflow",
    "spOccupancy Workflow", "s-jSDM Workflow", "boral Workflow",
    "Universal Benchmark", "Project", "Compare Models"
  ]);
  const controlSelector = ".sync-pair,.shiny-input-container";
  const commandSelector = ".action-button, .shiny-download-link, .btn-file";
  const iconSelector = "[data-studio-icon]";
  const updateSelector = [controlSelector, commandSelector, iconSelector].join(",");
  const engineNav = {
    "Hmsc Workflow": "engine-hmsc",
    "Hmsc-HPC Workflow": "engine-hmschpc",
    "jSDM Workflow": "engine-jsdm",
    "GJAM Workflow": "engine-gjam",
    "spOccupancy Workflow": "engine-spoccupancy",
    "s-jSDM Workflow": "engine-sjsdm",
    "boral Workflow": "engine-boral"
  };
  const icons = window.StudioIcons;
  const alignIcon = typeof ResizeObserver === "function" ? new ResizeObserver(function (entries, observer) {
    entries.forEach(function ({ target }) {
      if (!target.isConnected) {
        observer.unobserve(target);
        return;
      }
      const row = target.closest(".studio-setting");
      const icon = row?.querySelector(".studio-setting-icon");
      const body = row?.querySelector(".studio-setting-body");
      if (!icon || !body) return;
      const title = target.getBoundingClientRect();
      if (!title.height) return;
      const shift = (title.top + title.bottom - icon.offsetHeight) / 2 - body.getBoundingClientRect().top;
      icon.style.marginBlockStart = Math.max(0, Math.round(shift)) + "px";
    });
  }) : null;

  function layout(root) {
    root.classList.add("studio-workbench");
    const isEngineWorkflow = Boolean(engineNav[root.dataset.value]);
    if (isEngineWorkflow) root.classList.add("studio-engine-workflow");
    root.querySelectorAll(controlSelector).forEach(function (control) {
      if (control.parentElement.closest(controlSelector) ||
          control.closest(".studio-setting,.studio-compact,.dataTables_wrapper,.dt-container")) return;
      const label = control.querySelector("label");
      if (!label || !label.textContent.trim()) return;
      const logicalId = (control.dataset.syncId || label.htmlFor).replace(/-selectized$/, "");
      const field = document.getElementById(logicalId) ||
        control.querySelector("input:not([type=range]),select,textarea") || control;
      const kind = icons.forInput(field.id);
      const row = document.createElement("div");
      row.className = "studio-setting";
      row.dataset.studioControl = field.id;
      if (isEngineWorkflow && field.type === "checkbox") row.classList.add("studio-setting-toggle");
      const icon = icons.create(kind, "art", "studio-setting-icon");
      const body = document.createElement("div");
      body.className = "studio-setting-body";
      const successor = control.nextElementSibling;
      const ownHelp = control.querySelector(".param-help");
      const help = ownHelp || (successor?.classList.contains("param-help") ? successor : null);
      control.before(row);
      row.append(icon, body);
      body.append(control);
      label.classList.add("studio-setting-label");
      if (isEngineWorkflow && control.classList.contains("sync-pair")) {
        control.querySelector(".sync-number label")?.classList.add("visually-hidden");
        field.setAttribute("aria-label", label.textContent.trim() + " value");
      }
      alignIcon?.observe(label);
      if (help) {
        if (field.id) {
          help.id = field.id + "-meaning";
          const focused = document.getElementById(label.htmlFor);
          for (const described of new Set([field, focused].filter(Boolean))) {
            described.setAttribute("aria-describedby",
              [described.getAttribute("aria-describedby"), help.id].filter(Boolean).join(" "));
          }
        }
        body.append(help);
      }
      for (let ancestor = row.parentElement; ancestor && ancestor !== root; ancestor = ancestor.parentElement) {
        if (ancestor.classList.contains("row")) ancestor.classList.add("studio-form-stack");
      }
    });
  }

  function decorateCommands() {
    document.querySelectorAll(commandSelector).forEach(function (button) {
      if (button.dataset.studioCommand || button.closest(".dataTables_wrapper,.dt-container")) return;
      const kind = icons.forCommand(button);
      button.dataset.studioCommand = kind;
      button.classList.add("studio-has-command-icon");
      button.prepend(icons.create(kind, "glyph"));
    });
  }

  function decorateEngineMenu() {
    document.querySelectorAll("#main_nav .dropdown-menu a[data-value]").forEach(function (link) {
      const kind = engineNav[link.dataset.value];
      if (kind && !link.querySelector(iconSelector)) link.prepend(icons.create(kind, "glyph"));
    });
    const toggle = document.querySelector(".navbar .navbar-toggler,.navbar .navbar-toggle");
    if (toggle && !toggle.querySelector(iconSelector)) {
      toggle.querySelectorAll(".navbar-toggler-icon,.icon-bar").forEach((symbol) => symbol.remove());
      toggle.append(icons.create("menu", "glyph"));
      toggle.classList.add("studio-nav-toggle");
    }
  }

  function navigate(value) {
    const link = Array.from(document.querySelectorAll("#main_nav a[data-value]"))
      .find(function (item) { return item.dataset.value === value; });
    if (!link) return;
    link.click();
    const menu = link.closest(".dropdown-menu");
    const menuToggle = menu?.parentElement.querySelector(".dropdown-toggle");
    if (menu?.classList.contains("show") && menuToggle) {
      if (window.bootstrap?.Dropdown) {
        window.bootstrap.Dropdown.getOrCreateInstance(menuToggle).hide();
      } else {
        menu.classList.remove("show");
        menuToggle.classList.remove("show");
        menuToggle.setAttribute("aria-expanded", "false");
      }
    }
    const toggle = document.querySelector(".navbar .navbar-toggler[aria-expanded=true],.navbar .navbar-toggle[aria-expanded=true]");
    if (toggle) toggle.click();
    window.scrollTo({ left: 0, top: 0, behavior: "instant" });
  }

  window.studioNavigate = navigate;
  document.addEventListener("click", function (event) {
    const stage = event.target.closest("[data-workflow-stage]");
    if (stage) {
      const pane = stage.closest(".studio-engine-workflow");
      const item = pane?.querySelector(
        `.studio-workflow-form > .accordion-item[data-value="${stage.dataset.workflowStage}"]`
      );
      const toggle = item?.querySelector(":scope > .accordion-header > button");
      if (toggle) {
        if (toggle.getAttribute("aria-expanded") === "false") toggle.click();
        toggle.focus({ preventScroll: true });
        item.scrollIntoView({ block: "start" });
      }
      return;
    }
    const link = event.target.closest("[data-studio-nav]");
    if (!link) return;
    event.preventDefault();
    navigate(link.dataset.studioNav);
  });

  document.addEventListener("DOMContentLoaded", function () {
    let queued = false;
    function refresh() {
      if (queued) return;
      queued = true;
      requestAnimationFrame(function () {
        queued = false;
        document.querySelectorAll(".tab-pane[data-value]").forEach(function (pane) {
          if (pages.has(pane.dataset.value)) layout(pane);
        });
        document.querySelectorAll(iconSelector).forEach(icons.paint);
        decorateCommands();
        decorateEngineMenu();
      });
    }
    refresh();
    new MutationObserver(function (changes) {
      if (changes.some(function (change) {
        return Array.from(change.addedNodes).some(function (node) {
          return node.nodeType === Node.ELEMENT_NODE &&
            (node.matches(updateSelector) || node.querySelector(updateSelector));
        });
      })) refresh();
    }).observe(document.body, { childList: true, subtree: true });
  });
})();
