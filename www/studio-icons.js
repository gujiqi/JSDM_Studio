(function () {
  "use strict";

  // Role-specific illustrations and compact Lucide command glyphs share these IDs.
  const roles = {
    response: ["table-2", "blue", "Community responses"],
    predictor: ["layers", "teal", "Environmental predictors"],
    trait: ["leaf", "leaf", "Species traits"],
    species: ["shapes", "leaf", "Species or response identities"],
    group: ["list-tree", "clay", "Grouping or sample design"],
    coordinate: ["map-pin", "cyan", "Spatial coordinates"],
    newcoords: ["map-pin-plus", "cyan", "Coordinates for new sites"],
    spatial: ["map", "cyan", "Spatial structure"],
    distance: ["ruler", "cyan", "Distances"],
    phylogeny: ["git-fork", "violet", "Phylogeny or taxonomy"],
    fold: ["columns-3", "amber", "Validation groups"],
    visits: ["calendar-days", "clay", "Surveys or observation effort"],
    occurrence: ["target", "teal", "Occurrence process"],
    detection: ["radar", "rose", "Detection process"],
    heldout: ["file-check-2", "blue", "Held-out observations"],
    prediction: ["chart-spline", "cyan", "Prediction"],
    newdata: ["file-plus-2", "cyan", "Predictors for new sites"],
    family: ["chart-no-axes-column-increasing", "amber", "Response distribution"],
    mixed: ["combine", "amber", "Mixed responses or observation limits"],
    formula: ["square-function", "indigo", "Model formula"],
    model: ["boxes", "indigo", "Model specification"],
    "reduced-rank": ["layers", "indigo", "Reduced-rank covariates"],
    "predictor-matrix": ["table-2", "teal", "Covariate matrix"],
    "rank-dimension": ["chart-no-axes-column-increasing", "violet", "Reduced-rank dimension"],
    latent: ["orbit", "violet", "Latent factors"],
    mcmc: ["repeat-2", "violet", "Sampler iterations"],
    chain: ["link-2", "violet", "MCMC chains"],
    parallel: ["cpu", "blue", "Compute resources"],
    seed: ["dice-5", "amber", "Random seed"],
    prior: ["chart-no-axes-combined", "violet", "Prior distribution"],
    start: ["flag", "clay", "Starting values"],
    tuning: ["sliders-vertical", "indigo", "Fitting or tuning controls"],
    effect: ["sigma", "rose", "Coefficients and effects"],
    scale: ["scaling", "slate", "Scaling or transformation"],
    offset: ["square-function", "clay", "Fixed offset on the linear predictor"],
    association: ["network", "rose", "Species associations"],
    partition: ["chart-column-stacked", "teal", "Variation partitioning"],
    order: ["columns-3", "slate", "Plot ordering"],
    ordination: ["chart-scatter", "indigo", "Ordination"],
    diagnostic: ["activity", "leaf", "Diagnostics"],
    accuracy: ["gauge", "leaf", "Model fit or predictive assessment"],
    plot: ["chart-line", "teal", "Plots"],
    save: ["save", "slate", "Saved model"],
    report: ["file-text", "slate", "Reports and tables"],
    script: ["file-code-2", "slate", "Reproducible scripts"],
    archive: ["package", "clay", "Results archive"],
    export: ["download", "clay", "Export"],
    folder: ["folder-open", "blue", "Project or run folder"],
    question: ["circle-help", "slate", "Research question"],
    setting: ["settings-2", "slate", "Workflow setting"],
    log: ["list", "slate", "Run messages"],

    "engine-hmsc": ["git-fork", "leaf", "Hmsc"],
    "engine-hmschpc": ["cpu", "blue", "Hmsc-HPC"],
    "engine-jsdm": ["orbit", "violet", "jSDM"],
    "engine-gjam": ["chart-column-stacked", "amber", "GJAM"],
    "engine-spoccupancy": ["radar", "rose", "spOccupancy"],
    "engine-sjsdm": ["grid-3x3", "teal", "sjSDM"],
    "engine-boral": ["chart-scatter", "indigo", "boral"],
    benchmark: ["flask-conical", "violet", "Universal Benchmark"],
    compare: ["columns-2", "teal", "Compare models"],
    guide: ["book-open", "indigo", "Guides"],
    update: ["refresh-cw", "teal", "Updates"],
    feedback: ["message-square", "rose", "Feedback"],
    upload: ["upload", "blue", "Upload"],
    check: ["shield-check", "leaf", "Validation"],
    run: ["play", "violet", "Fit model"],
    preview: ["scan-eye", "cyan", "Preview"],
    add: ["plus", "teal", "Add"],
    remove: ["minus", "slate", "Remove"],
    search: ["search", "teal", "Inspect"],
    arrow: ["arrow-right", "slate", "Open"],
    external: ["external-link", "slate", "Open external page"],
    menu: ["menu", "slate", "Toggle navigation"]
  };

  const illustrations = new Set(`
    response predictor trait species group coordinate newcoords spatial distance phylogeny
    fold visits occurrence detection heldout prediction newdata family mixed formula model
    reduced-rank predictor-matrix rank-dimension latent mcmc chain parallel seed prior start tuning effect scale offset
    association partition order ordination
    diagnostic accuracy plot save report script archive export folder question setting log
    benchmark compare guide feedback update
    upload check run preview engine-sjsdm
  `.trim().split(/\s+/));
  const engineArt = {
    "engine-hmsc": "phylogeny",
    "engine-hmschpc": "parallel",
    "engine-jsdm": "latent",
    "engine-gjam": "mixed",
    "engine-spoccupancy": "detection",
    "engine-boral": "ordination"
  };

  // Stable Shiny input IDs, not words in labels or engine names, define meaning.
  const fields = {
    folder: "project_name builtin_case",
    question: "project_question question_template",
    response: `y_file long_file long_response_col response_argument empirical_y
      project_many_zeros trimy trim_minobs missingy n_sites project_n_sites`,
    predictor: `x_file env_file xdata_file empirical_x occ_file env_list missingx
      x_ind project_n_predictors assembly_predictor`,
    "reduced-rank": "use_xrrr",
    "predictor-matrix": "xrrr_file",
    "rank-dimension": "ncrrr",
    trait: `tr_file trait_file traits_file trait_spec_file trait_types_file use_traits allow_traits
      which_traits trait_list do_traits project_has_traits traits_intercept_warning`,
    species: "species_file species_list long_species_col ind_spp n_species project_n_responses",
    group: `study_file units_column random_column random_effect_column random_mode random_name
      random_level_type random_n ranlevelsused random rowids_file ranefids_file row_eff
      use_ranef site_effect long_site_col group_file do_ranefsplot`,
    coordinate: "coord_file coords_file lon_col lat_col coord_x coord_y longlat project_has_spatial",
    newcoords: "newcoords_file",
    spatial: `spatial_file spatial_model s_knot_file sknot_file spatial_method smethod nneighbours
      nngp n_neighbors search_type cov_model svc_cols generate_spatial_ev spatial_ev_k
      spatial_ev_threshold spatial_new_file include_space_in_anova conv_alpha`,
    distance: "distmat_file",
    phylogeny: "phylo_file phylo_tree_file taxonomy_file use_phylogeny phylogeny_mode project_has_phylogeny plottree conv_rho",
    fold: "cv_file nfolds partition_column holdout_file holdoutn cv_k out_cv",
    visits: "trials_file trials trial_size effort_file use_effort weights_file do_weights n_visits n_reps project_has_integrated_sources",
    occurrence: "occ_formula",
    detection: "det_file det_formula empirical_visits",
    heldout: "empirical_test_y",
    newdata: "newdata_file empirical_new_x newstudy_file prediction_ids_file",
    prediction: `project_need_prediction do_predict do_conditional do_fitted predict_type predict_probs
      predict_max_sites id_sites id_species predictx inverse_prediction out_pred
      out_predicted out_prediction out_predictions out_gradients prediction_histograms table_predictions`,
    family: "distr distribution family family_text project_response_structure project_response_family type_single type_file typenames_text",
    mixed: "mixed censor_file use_censor censor_columns censor_values censor_intervals fcgroups ccgroups composition_reference",
    formula: `xformula trformula xrrrformula trait_formula site_formula random_slope_formula formula
      env_formula spatial_formula formula_x include_intercept`,
    model: `model_type model_mode model_name env_model run_sampler full dtype precision tnlib
      eager reduct reduct_n reduct_r use_ssvs ssvs_index ssvs_traitsindex bias real_fit do_fit`,
    latent: `n_latent latent_factors nf nfmin nfmax n_factors num_lv lv_type
      truncatenumberoffactors constrained_latent save_eta do_lvsplot dnn_hidden`,
    mcmc: `burnin mcmc samples transient thin iter sampling n_burn n_thin n_batch
      batch_length n_burnin n_iteration ng hmcleapfrog hmcthin sample_prior save_mcmc_samples anova_samples
      out_mcmc_rds save_samples out_chains`,
    chain: "chains nchains n_chains chains_to_run constrained_nchains pool_chains alignpost",
    parallel: "n_omp_threads nparallel parallel device",
    seed: "seed",
    prior: `prior_file prior_mode use_prior_template use_default_priors a1 b1 a2 b2 alphapw
      v_alpha shape_valpha rate_valpha shape_v rate_v mu_beta v_beta mu_gamma v_gamma
      mu_lambda v_lambda ssvs_g prior_type hypparams ematalpha`,
    start: "initpar beta_start gamma_start lambda_start w_start alpha_start v_start",
    tuning: `ropt accept_rate hmc_update step_size env_lambda env_alpha spatial_lambda
      spatial_alpha biotic_lambda biotic_alpha biotic_df on_diag reg_on_cov inverse
      tune_regularization tune_steps activation dropout optimizer learning_rate
      weight_decay scheduler lr_reduce_factor early_stopping`,
    effect: `update_beta_eta updater_gammaeta updater_beta updater_gamma table_beta table_gamma
      table_lambda table_alpha beta_support gamma_support plot_beta plot_gamma out_params
      out_parameters do_coef do_coefsplot`,
    scale: "xrrrscale xscale trscale yscale scale_site_data notstandard effect_scale",
    offset: "offset_file use_offset loff_file",
    association: `updater_omega conv_omega maxomega omega_support out_omega
      show_sp_names_omega cor_prob cor_type plot_residual_cor plot_associations
      out_resid out_env do_env_cor do_resid_cor do_cov_cor association_strength`,
    partition: "out_vp do_varpart do_anova do_internal internal_fractions do_iie",
    order: "varpart_order_explained varpart_order_raw omega_order",
    ordination: "do_ordination",
    diagnostic: `out_diag out_diagnostics conv_beta conv_gamma effective_size gelman_psrf
      diag_beta diag_lambda diag_w diag_alpha diag_valpha diag_v diag_deviance
      out_trace do_ppc do_residuals do_residual_plot out_missing save_hpd se`,
    accuracy: `out_fit out_waic do_waic calc_ics save_dic do_rsquared do_importance
      do_sensitivity compute_sair residual_sd observation_noise_sd`,
    plot: "do_plots out_plots out_figures show_sp_names_beta do_assembly",
    save: "save_model out_model",
    report: `out_report out_tables out_csv_tables do_summary do_tidy coda_summary
      table_model_spec create_standard out_config out_copy_inputs`,
    script: "export_scripts out_jags python python_source",
    archive: "out_zip data_zip integrated_file",
    log: "verbose n_report",
    setting: "profile preset empirical_profile",
    benchmark: "engines benchmark_engines empirical_engines truth_profile data_source generating_structure allow_large",
    compare: "compare_engines"
  };
  const inputRoles = new Map(Object.entries(fields).flatMap(([kind, ids]) =>
    ids.trim().split(/\s+/).map((id) => [id, kind])));
  const specificInputRoles = {
    spocc_y_file: "detection",
    spocc_occ_file: "occurrence",
    hmschpc_alpha: "spatial",
    hmschpc_random_slope_formula: "formula",
    hmsc_random_effect_column: "group",
    hmsc_smethod: "spatial"
  };
  const enginePrefix = /^(?:hmschpc|hmsc|jsdm|gjam|spocc|sjsdm|boral)_/;

  function forInput(id) {
    const key = id.toLowerCase();
    if (/^hmsc_group_col_\d+$/.test(key)) return "group";
    if (/^compare_[a-z]+_dir$/.test(key)) return "folder";
    if (/_file_preview(?:_|$)/.test(key)) return "preview";
    return specificInputRoles[key] ||
      inputRoles.get(key.replace(/^univ_/, "").replace(enginePrefix, "")) ||
      "setting";
  }

  const commands = {
    home_open_project: "folder",
    home_open_engine_guide: "guide",
    hmsc_add_group: "add",
    univ_generate: "benchmark",
    univ_use_builtin: "folder",
    univ_import_zip: "upload",
    univ_preflight: "check",
    univ_compare: "update",
    compare_run: "search"
  };
  function forCommand(button) {
    if (button.matches(".btn-file")) return "upload";
    if (button.matches(".shiny-download-link")) return "export";
    const id = button.id;
    if (commands[id]) return commands[id];
    if (/^hmsc_group_remove_\d+$/.test(id)) return "remove";
    if (/_check$/.test(id)) return "check";
    if (/_run(?:_all|_matrix)?$/.test(id)) return "run";
    return "setting";
  }

  const svgNs = "http://www.w3.org/2000/svg";
  const assetRevision = document.currentScript ? new URL(document.currentScript.src).search : "";

  function paint(icon) {
    if (icon.dataset.iconReady === "true") return;
    const kind = icon.dataset.studioIcon;
    const [symbol, tone, title] = roles[kind] || roles.setting;
    const artwork = engineArt[kind] || kind;
    icon.dataset.iconMode ||= "art";
    const illustrated = icon.dataset.iconMode === "art" && illustrations.has(artwork);
    icon.dataset.iconReady = "true";
    icon.dataset.iconTone = tone;
    if (icon.dataset.iconMode === "art") icon.title = title;
    if (illustrated) icon.dataset.iconIllustrated = "true";

    const svg = document.createElementNS(svgNs, "svg");
    svg.setAttribute("viewBox", illustrated ? "0 0 40 40" : "0 0 24 24");
    svg.setAttribute("focusable", "false");
    const use = document.createElementNS(svgNs, "use");
    const file = illustrated ? "studio-pictograms.svg" : "studio-glyphs.svg";
    use.setAttribute("href", file + assetRevision + "#" + (illustrated ? artwork : symbol));
    svg.append(use);
    icon.append(svg);
  }

  function create(kind, mode = "art", extraClass = "") {
    const icon = document.createElement("span");
    icon.className = ("studio-icon " + extraClass).trim();
    icon.dataset.studioIcon = kind;
    icon.dataset.iconMode = mode;
    icon.setAttribute("aria-hidden", "true");
    paint(icon);
    return icon;
  }

  window.StudioIcons = Object.freeze({ paint, create, forInput, forCommand });
})();
