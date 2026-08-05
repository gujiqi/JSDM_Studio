"""Run a CPU Hmsc-HPC / pyhmsc workflow inside a JSDM Studio output folder."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any


def _ensure_packages(source_dir: str | None) -> None:
    if source_dir:
        src = Path(source_dir)
        if src.exists():
            sys.path.insert(0, str(src))


def _read_yaml(path: Path) -> dict[str, Any]:
    import yaml

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return data


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _run_cmd(cmd: list[str], cwd: Path, log_path: Path, env: dict[str, str]) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write("COMMAND: " + " ".join(cmd) + "\n\n")
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return int(proc.returncode)


def _read_csv(path: Path):
    import pandas as pd

    return pd.read_csv(path, index_col=0)


def _safe_summary(fit, param: str, out_path: Path, warnings: list[str]) -> None:
    try:
        tab = fit.summary(param)
        tab.to_csv(out_path, index=param == "Beta")
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"{param} summary skipped: {exc}")
        out_path.write_text(f"summary skipped: {exc}\n", encoding="utf-8")


def _array_to_csv(value: Any, out_path: Path) -> None:
    import numpy as np
    import pandas as pd

    if hasattr(value, "to_dataframe"):
        value.to_dataframe().reset_index().to_csv(out_path, index=False)
        return
    if hasattr(value, "to_array"):
        value.to_array().to_dataframe(name="value").reset_index().to_csv(out_path, index=False)
        return
    arr = np.asarray(value)
    flat = arr.reshape((-1, arr.shape[-1])) if arr.ndim > 1 else arr.reshape((-1, 1))
    pd.DataFrame(flat).to_csv(out_path, index=False)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return data if isinstance(data, dict) else {}


def _write_session_info(path: Path, warnings: list[str]) -> None:
    import importlib.metadata
    import platform

    packages = ["numpy", "pandas", "h5py", "tensorflow", "tensorflow-probability", "pyhmsc"]
    lines = [
        "Hmsc-HPC session information",
        "============================",
        f"Python: {sys.version.replace(os.linesep, ' ')}",
        f"Executable: {sys.executable}",
        f"Platform: {platform.platform()}",
        "",
        "Packages:",
    ]
    for pkg in packages:
        try:
            version = importlib.metadata.version(pkg)
        except Exception as exc:  # noqa: BLE001
            version = f"not available ({exc})"
            if pkg in {"tensorflow", "pyhmsc"}:
                warnings.append(f"Session info could not read {pkg}: {exc}")
        lines.append(f"- {pkg}: {version}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_hdf5_summary(path: Path, out_path: Path, warnings: list[str]) -> None:
    import pandas as pd

    rows: list[dict[str, Any]] = []
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            def visit(name: str, obj: Any) -> None:
                if hasattr(obj, "shape"):
                    rows.append(
                        {
                            "object": name,
                            "kind": "dataset",
                            "shape": " x ".join(str(x) for x in obj.shape),
                            "dtype": str(obj.dtype),
                        }
                    )
                else:
                    rows.append({"object": name, "kind": "group", "shape": "", "dtype": ""})

            handle.visititems(visit)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"HDF5 posterior summary skipped: {exc}")
        rows.append({"object": str(path), "kind": "error", "shape": "", "dtype": str(exc)})
    pd.DataFrame(rows).to_csv(out_path, index=False)


def _long_matrix(df: Any, row_name: str, col_name: str, value_name: str):
    import pandas as pd

    out = df.copy()
    out.index.name = row_name
    return out.reset_index().melt(id_vars=row_name, var_name=col_name, value_name=value_name)


def _support_label(lower: Any, upper: Any, mean: Any) -> str:
    try:
        lo = float(lower)
        hi = float(upper)
        mu = float(mean)
    except Exception:  # noqa: BLE001
        return "unknown"
    if lo > 0:
        return "supported_positive"
    if hi < 0:
        return "supported_negative"
    if mu > 0:
        return "uncertain_positive"
    if mu < 0:
        return "uncertain_negative"
    return "uncertain_zero"


def _write_data_check_messages(workdir: Path, cfg: dict[str, Any], X: Any, Y: Any, warnings: list[str]) -> None:
    import numpy as np
    import pandas as pd

    key = str(cfg.get("model", {}).get("distribution", "")).lower()
    ym = Y.to_numpy(dtype=float)
    if key == "probit":
        family_ok = bool(np.isin(ym[~np.isnan(ym)], [0, 1]).all())
        family_msg = "probit response check: Y values are 0/1." if family_ok else "probit response check failed: Y must be 0/1."
    elif key == "poisson":
        vals = ym[~np.isnan(ym)]
        family_ok = bool(((vals >= 0) & (np.abs(vals - np.round(vals)) < 1e-8)).all())
        family_msg = "poisson response check: Y values are non-negative integers." if family_ok else "poisson response check failed: Y must be non-negative integers."
    else:
        family_ok = True
        family_msg = "normal/gaussian response check: continuous numeric Y is allowed."
    cat_cols = [name for name in X.columns if not pd.api.types.is_numeric_dtype(X[name])]
    rows = [
        {"level": "info", "check": "dimensions", "ok": True, "message": f"Y rows={Y.shape[0]}, Y responses={Y.shape[1]}, X rows={X.shape[0]}, X predictors={X.shape[1]}."},
        {"level": "info", "check": "distribution", "ok": family_ok, "message": family_msg},
        {"level": "info", "check": "formula", "ok": True, "message": f"X formula: {cfg.get('model', {}).get('XFormula', '~ .')}"},
        {"level": "info", "check": "categorical_predictors", "ok": True, "message": "Categorical predictors: " + (", ".join(cat_cols) if cat_cols else "none")},
        {"level": "info", "check": "random_effects", "ok": True, "message": f"Random mode: {cfg.get('random_effects', {}).get('mode', 'none')}"},
        {"level": "info", "check": "outputs", "ok": True, "message": f"predictions={cfg.get('outputs', {}).get('predictions', True)}, diagnostics={cfg.get('outputs', {}).get('diagnostics', True)}, plots={cfg.get('outputs', {}).get('plots', True)}"},
    ]
    if warnings:
        rows.extend({"level": "warning", "check": "runtime_warning", "ok": False, "message": msg} for msg in warnings)
    pd.DataFrame(rows).to_csv(workdir / "diagnostics" / "data_check_messages.csv", index=False)


def _write_hmschpc_step_scripts(workdir: Path, cfg: dict[str, Any]) -> None:
    workflow = workdir / "workflow_scripts"
    workflow.mkdir(parents=True, exist_ok=True)
    step_specs = [
        ("S1_define_models.py", "Define the pyhmsc model boundary from used_config.yml and workflow_scripts/hmschpc_model.yaml."),
        ("S2_fit_models.py", "Run pyhmsc compile, validate-init and the CPU TensorFlow sampler."),
        ("S3_evaluate_convergence.py", "Read samples/posterior.h5 and export Rhat/ESS diagnostics where ArviZ supports them."),
        ("S4_compute_model_fit.py", "Compute training fitted means and per-response RMSE/R2 summaries."),
        ("S5_show_model_fit.py", "Write model-fit tables and observed-versus-fitted plots."),
        ("S6_show_parameter_estimates.py", "Export Beta, Gamma, sigma, rho and random-level Eta/Lambda summaries."),
        ("S7_make_predictions.py", "Export training and optional newdata predictions."),
    ]
    for filename, description in step_specs:
        (workflow / filename).write_text(
            "\n".join(
                [
                    "#!/usr/bin/env python",
                    '"""Hmsc-HPC workflow step exported by JSDM Studio."""',
                    "from pathlib import Path",
                    "",
                    "root = Path(__file__).resolve().parents[1]",
                    f"description = {description!r}",
                    "print(description)",
                    "print('Output folder:', root)",
                    "print('This step is implemented by reproducible_script/run_this_HmscHPC_analysis.py.')",
                    "print('Re-run the master script to reproduce all Hmsc-HPC S1-S7 outputs.')",
                    "",
                ]
            ),
            encoding="utf-8",
        )
    (workflow / "S1_to_S7_step_map.csv").write_text(
        "step,file,meaning\n"
        "S1,workflow_scripts/hmschpc_model.yaml,pyhmsc model definition\n"
        "S2,samples/posterior.h5,CPU sampler posterior file\n"
        "S3,diagnostics/Beta_rhat.csv and diagnostics/Beta_ess.csv,convergence diagnostics\n"
        "S4,tables/fit_metrics.csv,model-fit metrics\n"
        "S5,plots/prediction_observed_vs_fitted.png,model-fit figure\n"
        "S6,tables/Beta_summary.csv and random-level tables,parameter estimates\n"
        "S7,predictions/predicted_mean.csv,predictions\n",
        encoding="utf-8",
    )


def _write_outputs(workdir: Path, cfg: dict[str, Any], posterior: Path, warnings: list[str]) -> None:
    import numpy as np
    import pandas as pd
    from pyhmsc.posterior import HmscFit

    tables = workdir / "tables"
    preds_dir = workdir / "predictions"
    plots = workdir / "plots"
    diagnostics = workdir / "diagnostics"
    standard = workdir / "standard"
    models = workdir / "models"
    results = workdir / "results"
    report_dir = workdir / "report"
    for directory in (tables, preds_dir, plots, diagnostics, standard, models, results, report_dir):
        directory.mkdir(parents=True, exist_ok=True)

    model_meta = _read_json(workdir / "models" / "compiled_model" / "init.json")
    _write_hmschpc_step_scripts(workdir, cfg)
    _write_hdf5_summary(posterior, diagnostics / "posterior_hdf5_summary.csv", warnings)

    fit = HmscFit.from_file(posterior)
    beta = fit.summary("Beta")
    beta.to_csv(tables / "Beta_summary.csv", index=False)
    beta.to_csv(results / "S6_parameter_estimates_Beta.csv", index=False)
    beta_mean = fit.beta_mean()
    beta_mean.to_csv(tables / "Beta_mean_matrix.csv")
    _long_matrix(beta_mean, "predictor", "response_id", "estimate").to_csv(
        tables / "Beta_mean_long.csv", index=False
    )
    beta_support = beta.copy()
    if all(name in beta_support.columns for name in ("lower", "upper", "mean")):
        beta_support["support"] = [
            _support_label(lo, hi, mu)
            for lo, hi, mu in zip(beta_support["lower"], beta_support["upper"], beta_support["mean"])
        ]
    beta_support.to_csv(tables / "Beta_support.csv", index=False)
    beta_support.to_csv(results / "S6_Beta_support.csv", index=False)
    _safe_summary(fit, "Gamma", tables / "Gamma_summary.csv", warnings)
    _safe_summary(fit, "sigma", tables / "sigma_summary.csv", warnings)
    try:
        gamma_mean = fit.gamma_mean()
        gamma_names = model_meta.get("names", {}).get("traits") or list(gamma_mean.columns)
        cov_names = model_meta.get("names", {}).get("covariates") or list(gamma_mean.index)
        if len(cov_names) == gamma_mean.shape[0]:
            gamma_mean.index = cov_names
        if len(gamma_names) == gamma_mean.shape[1]:
            gamma_mean.columns = gamma_names
        gamma_mean.to_csv(tables / "Gamma_mean_matrix.csv")
        _long_matrix(gamma_mean, "predictor", "trait", "estimate").to_csv(
            results / "S6_parameter_estimates_Gamma.csv", index=False
        )
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"Gamma mean export skipped: {exc}")
    try:
        sigma_mean = fit.sigma_mean()
        sigma_mean.rename("sigma_mean").reset_index(names="response_id").to_csv(
            tables / "sigma_mean.csv", index=False
        )
        sigma_mean.rename("sigma_mean").reset_index(names="response_id").to_csv(
            results / "S6_sigma_summary.csv", index=False
        )
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"sigma mean export skipped: {exc}")
    try:
        rho_mean = fit.rho_mean()
        rho_mean.rename("rho_mean").reset_index(names="rho_parameter").to_csv(
            tables / "rho_mean.csv", index=False
        )
        rho_mean.rename("rho_mean").reset_index(names="rho_parameter").to_csv(
            results / "S6_rho_summary.csv", index=False
        )
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"rho mean export skipped: {exc}")
    for param, fname in (("Beta", "Beta"), ("Gamma", "Gamma"), ("sigma", "sigma"), ("rhoInd", "rhoInd")):
        try:
            _array_to_csv(fit._samples(param), tables / f"{fname}_draws_flat.csv")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{param} draws export skipped: {exc}")

    x_path = workdir / "data" / "X.csv"
    y_path = workdir / "data" / "Y.csv"
    X = _read_csv(x_path)
    Y = _read_csv(y_path)
    _write_data_check_messages(workdir, cfg, X, Y, warnings)
    model_rows = [
        {"field": "engine", "value": "Hmsc-HPC"},
        {"field": "distribution", "value": cfg.get("model", {}).get("distribution")},
        {"field": "XFormula", "value": cfg.get("model", {}).get("XFormula")},
        {"field": "trait_formula", "value": cfg.get("model", {}).get("trait_formula")},
        {"field": "use_traits", "value": cfg.get("model", {}).get("use_traits")},
        {"field": "phylogeny_mode", "value": cfg.get("model", {}).get("phylogeny_mode")},
        {"field": "random_mode", "value": cfg.get("random_effects", {}).get("mode")},
        {"field": "random_name", "value": cfg.get("random_effects", {}).get("name")},
        {"field": "n_sites", "value": int(Y.shape[0])},
        {"field": "n_responses", "value": int(Y.shape[1])},
        {"field": "n_predictors", "value": int(X.shape[1])},
        {"field": "posterior_file", "value": str(posterior.relative_to(workdir)).replace("\\", "/")},
    ]
    pd.DataFrame(model_rows).to_csv(tables / "S1_model_definition.csv", index=False)
    pd.DataFrame(model_rows).to_csv(results / "S1_model_definition.csv", index=False)
    pd.DataFrame(
        [
            {
                "step": "S2_fit_models",
                "status": "fitted",
                "posterior_file": str(posterior.relative_to(workdir)).replace("\\", "/"),
                "samples": cfg.get("sampler", {}).get("samples"),
                "transient": cfg.get("sampler", {}).get("transient"),
                "thin": cfg.get("sampler", {}).get("thin"),
                "chains": cfg.get("sampler", {}).get("chains"),
            }
        ]
    ).to_csv(results / "S2_fit_models.csv", index=False)
    export_predictions = bool(cfg.get("outputs", {}).get("predictions", True))
    pred = fit.predict_mean(X, random_effects="none")
    pred_long = (
        pred.reset_index(names="site_id")
        .melt(id_vars="site_id", var_name="response_id", value_name="predicted_mean")
    )
    obs_long = (
        Y.reset_index(names="site_id")
        .melt(id_vars="site_id", var_name="response_id", value_name="observed")
    )
    pred_long = pred_long.merge(obs_long, on=["site_id", "response_id"], how="left")
    pred_long["engine"] = "Hmsc-HPC"
    pred_long["predicted_lower"] = np.nan
    pred_long["predicted_upper"] = np.nan
    pred_long["prediction_set"] = "training"
    prediction_cols = [
        "engine",
        "site_id",
        "response_id",
        "observed",
        "predicted_mean",
        "predicted_lower",
        "predicted_upper",
        "prediction_set",
    ]
    pred_long = pred_long[prediction_cols]
    prediction_tables = [pred_long]
    if export_predictions:
        pred.to_csv(preds_dir / "predicted_mean.csv")
        pred.to_csv(preds_dir / "training_predicted_mean.csv")
        pred_long.to_csv(preds_dir / "training_predictions_long.csv", index=False)
        pred_long.to_csv(results / "S7_predictions_training.csv", index=False)
        newdata_path = workdir / "data" / "newdata.csv"
        if newdata_path.exists():
            try:
                X_new = _read_csv(newdata_path)
                pred_new = fit.predict_mean(X_new, random_effects="none")
                pred_new.to_csv(preds_dir / "newdata_predicted_mean.csv")
                new_long = (
                    pred_new.reset_index(names="site_id")
                    .melt(id_vars="site_id", var_name="response_id", value_name="predicted_mean")
                )
                new_long["engine"] = "Hmsc-HPC"
                new_long["observed"] = np.nan
                new_long["predicted_lower"] = np.nan
                new_long["predicted_upper"] = np.nan
                new_long["prediction_set"] = "newdata"
                prediction_tables.append(new_long[prediction_cols])
                new_long[prediction_cols].to_csv(preds_dir / "newdata_predictions_long.csv", index=False)
                new_long[prediction_cols].to_csv(results / "S7_predictions_newdata.csv", index=False)
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"newdata prediction skipped: {exc}")
        pred_long_all = pd.concat(prediction_tables, ignore_index=True)
        pred_long_all.to_csv(standard / "predictions_long.csv", index=False)
        pred_long_all.to_csv(preds_dir / "predictions_long.csv", index=False)
    else:
        pd.DataFrame(columns=prediction_cols).to_csv(standard / "predictions_long.csv", index=False)
        (preds_dir / "predictions_disabled.txt").write_text(
            "Prediction export was disabled by the Hmsc-HPC output settings. "
            "Fit metrics were still computed internally from training fitted means.\n",
            encoding="utf-8",
        )
        pd.DataFrame(columns=prediction_cols).to_csv(results / "S7_predictions_training.csv", index=False)

    metrics = []
    for sp in Y.columns:
        observed = pd.to_numeric(Y[sp], errors="coerce").to_numpy(dtype=float)
        fitted = pd.to_numeric(pred[sp], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(observed) & np.isfinite(fitted)
        if ok.any():
            rmse = float(np.sqrt(np.mean((observed[ok] - fitted[ok]) ** 2)))
            ss_res = float(np.sum((observed[ok] - fitted[ok]) ** 2))
            ss_tot = float(np.sum((observed[ok] - np.mean(observed[ok])) ** 2))
            r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
        else:
            rmse = float("nan")
            r2 = float("nan")
        metrics.append({"engine": "Hmsc-HPC", "metric": "RMSE", "response_id": sp, "value": rmse})
        metrics.append({"engine": "Hmsc-HPC", "metric": "R2_training", "response_id": sp, "value": r2})
    pd.DataFrame(metrics).to_csv(standard / "fit_metrics.csv", index=False)
    pd.DataFrame(metrics).to_csv(tables / "fit_metrics.csv", index=False)
    metric_df = pd.DataFrame(metrics)
    metric_df.to_csv(results / "S4_model_fit_by_response.csv", index=False)
    fit_summary = (
        metric_df.groupby("metric", dropna=False)["value"]
        .agg(["mean", "median", "min", "max"])
        .reset_index()
        .rename(columns={"mean": "mean_value", "median": "median_value", "min": "min_value", "max": "max_value"})
    )
    fit_summary.insert(0, "engine", "Hmsc-HPC")
    fit_summary.to_csv(tables / "model_fit_summary.csv", index=False)
    fit_summary.to_csv(results / "S4_model_fit_summary.csv", index=False)
    prediction_summary = (
        pred_long.groupby("response_id", dropna=False)
        .agg(
            observed_mean=("observed", "mean"),
            predicted_mean=("predicted_mean", "mean"),
            residual_mean=("observed", lambda x: np.nan),
        )
        .reset_index()
    )
    residual_lookup = pred_long.assign(residual=pred_long["observed"] - pred_long["predicted_mean"])
    residual_summary = residual_lookup.groupby("response_id", dropna=False)["residual"].agg(["mean", "std"]).reset_index()
    prediction_summary = prediction_summary.drop(columns=["residual_mean"]).merge(
        residual_summary.rename(columns={"mean": "residual_mean", "std": "residual_sd"}),
        on="response_id",
        how="left",
    )
    prediction_summary.insert(0, "engine", "Hmsc-HPC")
    prediction_summary.to_csv(tables / "prediction_summary_by_species.csv", index=False)
    prediction_summary.to_csv(results / "S5_model_fit_prediction_summary.csv", index=False)

    effects = beta.rename(
        columns={"covariate": "predictor", "species": "response_id", "mean": "estimate"}
    )
    effects["engine"] = "Hmsc-HPC"
    effects["direction"] = np.where(effects["estimate"] >= 0, "positive", "negative")
    effects["notes"] = "CPU pyhmsc Beta posterior summary"
    effects[["engine", "response_id", "predictor", "direction", "estimate", "lower", "upper", "notes"]].to_csv(
        standard / "effects_long.csv", index=False
    )
    species_names = model_meta.get("names", {}).get("species") or list(Y.columns)
    association_rows: list[dict[str, Any]] = []
    for level_index, level in enumerate(model_meta.get("random_levels", []) or []):
        level_name = str(level.get("name", f"level_{level_index}"))
        try:
            eta_mean = fit.eta_mean(level_index)
            eta_mean.to_csv(tables / f"random_level_{level_index}_{level_name}_Eta_mean.csv")
            eta_mean.to_csv(results / f"S6_random_level_{level_index}_{level_name}_Eta_mean.csv")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"Eta mean export skipped for random level {level_name}: {exc}")
        try:
            lambda_mean = fit.lambda_mean(level_index)
            lambda_mean.to_csv(tables / f"random_level_{level_index}_{level_name}_Lambda_mean.csv")
            lambda_mean.to_csv(results / f"S6_random_level_{level_index}_{level_name}_Lambda_mean.csv")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"Lambda mean export skipped for random level {level_name}: {exc}")
        try:
            eta_samples = fit.eta_samples(level_index)
            lambda_samples = fit.lambda_samples(level_index)
            _array_to_csv(eta_samples, tables / f"random_level_{level_index}_{level_name}_Eta_draws_flat.csv")
            _array_to_csv(lambda_samples, tables / f"random_level_{level_index}_{level_name}_Lambda_draws_flat.csv")
            lam = np.asarray(lambda_samples, dtype=float)
            lam = lam.reshape((-1, lam.shape[-2], lam.shape[-1]))
            cov_stack = np.array([draw.T @ draw for draw in lam])
            cov_mean = np.nanmean(cov_stack, axis=0)
            diag = np.sqrt(np.diag(cov_mean))
            corr = cov_mean / np.outer(diag, diag)
            corr[~np.isfinite(corr)] = np.nan
            cov_df = pd.DataFrame(cov_mean, index=species_names, columns=species_names)
            corr_df = pd.DataFrame(corr, index=species_names, columns=species_names)
            cov_df.to_csv(tables / f"random_level_{level_index}_{level_name}_association_covariance.csv")
            corr_df.to_csv(tables / f"random_level_{level_index}_{level_name}_association_correlation.csv")
            corr_df.to_csv(results / f"S6_random_level_{level_index}_{level_name}_association_correlation.csv")
            for i, sp1 in enumerate(species_names):
                for j, sp2 in enumerate(species_names):
                    if j <= i:
                        continue
                    association_rows.append(
                        {
                            "engine": "Hmsc-HPC",
                            "response_1": sp1,
                            "response_2": sp2,
                            "association_type": f"lambda_correlation:{level_name}",
                            "estimate": float(corr_df.loc[sp1, sp2]),
                            "comparable_level": "random_level_latent_factor",
                        }
                    )
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"Random-level association export skipped for {level_name}: {exc}")
    if not association_rows:
        association_rows = [
            {
                "engine": "Hmsc-HPC",
                "response_1": None,
                "response_2": None,
                "association_type": "not_available_for_selected_hmschpc_model",
                "estimate": None,
                "comparable_level": "not_estimated",
            }
        ]
    pd.DataFrame(association_rows).to_csv(standard / "associations_long.csv", index=False)
    pd.DataFrame(association_rows).to_csv(tables / "associations_long.csv", index=False)
    pd.DataFrame(association_rows).to_csv(results / "S6_associations_long.csv", index=False)

    if bool(cfg.get("outputs", {}).get("diagnostics", True)):
        try:
            rhat = fit.rhat("Beta")
            ess = fit.ess("Beta")
            _array_to_csv(rhat, diagnostics / "Beta_rhat.csv")
            _array_to_csv(ess, diagnostics / "Beta_ess.csv")
            rhat_arr = np.asarray(rhat, dtype=float).ravel()
            ess_arr = np.asarray(ess, dtype=float).ravel()
            conv = pd.DataFrame(
                [
                    {
                        "engine": "Hmsc-HPC",
                        "parameter": "Beta",
                        "rhat_mean": float(np.nanmean(rhat_arr)) if rhat_arr.size else np.nan,
                        "rhat_max": float(np.nanmax(rhat_arr)) if rhat_arr.size else np.nan,
                        "ess_mean": float(np.nanmean(ess_arr)) if ess_arr.size else np.nan,
                        "ess_min": float(np.nanmin(ess_arr)) if ess_arr.size else np.nan,
                        "note": "Rhat/ESS from pyhmsc ArviZ conversion where available.",
                    }
                ]
            )
            conv.to_csv(diagnostics / "S3_convergence_summary.csv", index=False)
            conv.to_csv(results / "S3_convergence_summary.csv", index=False)
            conv.to_csv(standard / "diagnostics_long.csv", index=False)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"ArviZ diagnostics skipped: {exc}")
            (diagnostics / "HmscHPC_diagnostics_warning.txt").write_text(str(exc), encoding="utf-8")
            conv = pd.DataFrame(
                [
                    {
                        "engine": "Hmsc-HPC",
                        "parameter": "Beta",
                        "rhat_mean": np.nan,
                        "rhat_max": np.nan,
                        "ess_mean": np.nan,
                        "ess_min": np.nan,
                        "note": f"diagnostics skipped: {exc}",
                    }
                ]
            )
            conv.to_csv(diagnostics / "S3_convergence_summary.csv", index=False)
            conv.to_csv(results / "S3_convergence_summary.csv", index=False)
            conv.to_csv(standard / "diagnostics_long.csv", index=False)
    else:
        (diagnostics / "HmscHPC_diagnostics_disabled.txt").write_text(
            "Posterior diagnostics were disabled by the Hmsc-HPC output settings. "
            "engine_status.json is still written for run auditing.\n",
            encoding="utf-8",
        )
        conv = pd.DataFrame(
            [
                {
                    "engine": "Hmsc-HPC",
                    "parameter": "Beta",
                    "rhat_mean": np.nan,
                    "rhat_max": np.nan,
                    "ess_mean": np.nan,
                    "ess_min": np.nan,
                    "note": "diagnostics disabled by user",
                }
            ]
        )
        conv.to_csv(diagnostics / "S3_convergence_summary.csv", index=False)
        conv.to_csv(results / "S3_convergence_summary.csv", index=False)
        conv.to_csv(standard / "diagnostics_long.csv", index=False)
    if not (standard / "diagnostics_long.csv").exists() and (diagnostics / "S3_convergence_summary.csv").exists():
        pd.read_csv(diagnostics / "S3_convergence_summary.csv").to_csv(
            standard / "diagnostics_long.csv", index=False
        )

    if not bool(cfg.get("outputs", {}).get("plots", True)):
        (plots / "README.txt").write_text(
            "Plot export was disabled by the Hmsc-HPC output settings.\n",
            encoding="utf-8",
        )
    else:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            beta_matrix = fit.beta_mean()
            fig, ax = plt.subplots(figsize=(7.5, 4.8))
            im = ax.imshow(beta_matrix.to_numpy(dtype=float), aspect="auto", cmap="viridis")
            ax.set_xticks(range(beta_matrix.shape[1]), beta_matrix.columns, rotation=45, ha="right")
            ax.set_yticks(range(beta_matrix.shape[0]), beta_matrix.index)
            ax.set_title("Hmsc-HPC Beta posterior mean")
            fig.colorbar(im, ax=ax, shrink=0.8)
            fig.tight_layout()
            fig.savefig(plots / "Beta_heatmap.png", dpi=180)
            plt.close(fig)

            if export_predictions:
                fig, ax = plt.subplots(figsize=(5.8, 5.2))
                ax.scatter(pred_long["observed"], pred_long["predicted_mean"], s=28, alpha=0.75)
                ax.set_xlabel("Observed")
                ax.set_ylabel("Predicted mean")
                ax.set_title("Hmsc-HPC fitted prediction check")
                fig.tight_layout()
                fig.savefig(plots / "prediction_observed_vs_fitted.png", dpi=180)
                plt.close(fig)
            else:
                (plots / "prediction_plot_disabled.txt").write_text(
                    "Prediction plots were skipped because prediction export was disabled.\n",
                    encoding="utf-8",
                )
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"Plot export skipped: {exc}")

    run_summary = pd.DataFrame(
        [
            {
                "run_id": workdir.name,
                "engine": "Hmsc-HPC",
                "status": "fitted",
                "n_sites": int(Y.shape[0]),
                "n_responses": int(Y.shape[1]),
                "n_predictors": int(X.shape[1]),
                "distribution": cfg.get("model", {}).get("distribution"),
                "random_mode": cfg.get("random_effects", {}).get("mode"),
                "predictions_exported": bool(export_predictions),
                "diagnostics_exported": bool(cfg.get("outputs", {}).get("diagnostics", True)),
                "plots_exported": bool(cfg.get("outputs", {}).get("plots", True)),
                "hmschpc_workflow_steps": "S1-S7",
                "posterior_file": str(posterior.relative_to(workdir)).replace("\\", "/"),
            }
        ]
    )
    run_summary.to_csv(standard / "run_summary.csv", index=False)
    run_summary.to_csv(results / "S2_fit_run_summary.csv", index=False)
    pd.concat(
        [
            pd.DataFrame({"result_group": ["model_definition"], "file": ["results/S1_model_definition.csv"]}),
            pd.DataFrame({"result_group": ["fit"], "file": ["results/S2_fit_models.csv"]}),
            pd.DataFrame({"result_group": ["convergence"], "file": ["results/S3_convergence_summary.csv"]}),
            pd.DataFrame({"result_group": ["model_fit"], "file": ["results/S4_model_fit_summary.csv"]}),
            pd.DataFrame({"result_group": ["model_fit"], "file": ["results/S5_model_fit_prediction_summary.csv"]}),
            pd.DataFrame({"result_group": ["parameters"], "file": ["results/S6_parameter_estimates_Beta.csv"]}),
            pd.DataFrame({"result_group": ["predictions"], "file": ["results/S7_predictions_training.csv"]}),
        ],
        ignore_index=True,
    ).to_csv(tables / "HmscHPC_S1S7_result_index.csv", index=False)

    (workdir / "results" / "README_HmscHPC_results.txt").write_text(
        "Hmsc-HPC CPU workflow completed.\n"
        "The output follows the Hmsc-style S1-S7 result chain:\n"
        "S1 model definition: results/S1_model_definition.csv\n"
        "S2 fitted posterior: samples/posterior.h5 and results/S2_fit_models.csv\n"
        "S3 convergence: diagnostics/S3_convergence_summary.csv\n"
        "S4 model fit: results/S4_model_fit_summary.csv\n"
        "S5 model-fit display tables/plots: results/S5_model_fit_prediction_summary.csv and plots/\n"
        "S6 parameter estimates: Beta/Gamma/sigma/rho and random-level Eta/Lambda tables\n"
        "S7 predictions: predictions/ and standard/predictions_long.csv\n",
        encoding="utf-8",
    )
    _write_session_info(diagnostics / "session_info.txt", warnings)
    report = report_dir / "Hmsc-HPC_report.html"
    report.write_text(
        "<!DOCTYPE html><html><head><meta charset='UTF-8'><title>Hmsc-HPC report</title>"
        "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:40px;line-height:1.55}"
        "table{border-collapse:collapse}td,th{border:1px solid #dbe5f0;padding:6px 10px}</style></head>"
        "<body><h1>Hmsc-HPC CPU workflow report</h1>"
        "<p>Status: fitted.</p>"
        "<p>The workflow compiled raw CSV inputs to pyhmsc JSON/HDF5, sampled with TensorFlow on CPU, "
        "and exported reproducible summaries.</p>"
        "<h2>Hmsc-style result chain</h2>"
        "<table><tr><th>Step</th><th>Primary output</th></tr>"
        "<tr><td>S1 define models</td><td>results/S1_model_definition.csv</td></tr>"
        "<tr><td>S2 fit models</td><td>samples/posterior.h5, results/S2_fit_models.csv</td></tr>"
        "<tr><td>S3 convergence</td><td>diagnostics/S3_convergence_summary.csv</td></tr>"
        "<tr><td>S4 compute model fit</td><td>results/S4_model_fit_summary.csv</td></tr>"
        "<tr><td>S5 show model fit</td><td>plots/prediction_observed_vs_fitted.png</td></tr>"
        "<tr><td>S6 parameter estimates</td><td>tables/Beta_summary.csv, random-level Eta/Lambda tables</td></tr>"
        "<tr><td>S7 predictions</td><td>predictions/predictions_long.csv</td></tr>"
        "</table>"
        "</body></html>",
        encoding="utf-8",
    )


def _write_manifest(workdir: Path) -> None:
    import pandas as pd

    files = [str(path.relative_to(workdir)).replace("\\", "/") for path in workdir.rglob("*") if path.is_file()]
    pd.DataFrame({"file": sorted(files)}).to_csv(workdir / "standard" / "output_manifest.csv", index=False)


def _fill_empty_dirs(workdir: Path) -> None:
    for path in workdir.rglob("*"):
        if path.is_dir() and not any(path.iterdir()):
            (path / "README.txt").write_text(
                "This Hmsc-HPC output folder was created by the workflow. "
                "No file was produced here for the selected settings.\n",
                encoding="utf-8",
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--python-source")
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    _ensure_packages(args.python_source)
    diagnostics = workdir / "diagnostics"
    diagnostics.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    errors: list[str] = []
    status = "fit_failed"
    cfg: dict[str, Any] = {}

    env = os.environ.copy()
    if args.python_source:
        env["PYTHONPATH"] = str(Path(args.python_source).resolve()) + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "1")

    try:
        cfg = _read_yaml(workdir / "used_config.yml")
        model_yaml = workdir / "workflow_scripts" / "hmschpc_model.yaml"
        compile_dir = workdir / "models" / "compiled_model"
        posterior = workdir / "samples" / "posterior.h5"
        python = cfg.get("runtime", {}).get("python", "") or sys.executable

        code = _run_cmd(
            [python, "-m", "pyhmsc", "compile", str(model_yaml), "--output", str(compile_dir)],
            cwd=workdir,
            log_path=diagnostics / "HmscHPC_compile.log",
            env=env,
        )
        if code != 0:
            raise RuntimeError("pyhmsc compile failed; see diagnostics/HmscHPC_compile.log")

        sampler = cfg.get("sampler", {})
        validate_cmd = [python, "-m", "pyhmsc", "validate-init", str(compile_dir / "init.json")]
        if bool(sampler.get("run_sampler", True)):
            validate_cmd.append("--strict")
        code = _run_cmd(
            validate_cmd,
            cwd=workdir,
            log_path=diagnostics / "HmscHPC_validate_init.log",
            env=env,
        )
        if code != 0:
            raise RuntimeError("pyhmsc validate-init failed; see diagnostics/HmscHPC_validate_init.log")

        if not bool(cfg.get("sampler", {}).get("run_sampler", True)):
            status = "model_defined"
            warnings.append("Sampler was intentionally skipped; compiled Hmsc-HPC model artifacts were produced.")
        else:
            save_eta = bool(sampler.get("save_eta", True))
            if not save_eta:
                warnings.append(
                    "save Eta was requested as FALSE, but the current Hmsc-HPC HDF5 exporter "
                    "expects Eta entries. The runner forced --fse 1 to keep posterior.h5 valid."
                )
                save_eta = True
            cmd = [
                python,
                "-m",
                "hmsc.run_gibbs_sampler",
                "--input",
                str(compile_dir / "init.json"),
                "--output",
                str(posterior),
                "--samples",
                str(int(sampler.get("samples", 10))),
                "--transient",
                str(int(sampler.get("transient", 10))),
                "--thin",
                str(int(sampler.get("thin", 1))),
                "--verbose",
                str(max(1, int(sampler.get("verbose", 5)))),
                "--rngseed",
                str(int(sampler.get("seed", 1234))),
                "--hmcleapfrog",
                str(int(sampler.get("hmcleapfrog", 10))),
                "--hmcthin",
                str(int(sampler.get("hmcthin", 0))),
                "--updbe",
                "1" if bool(sampler.get("update_beta_eta", False)) else "0",
                "--tnlib",
                str(sampler.get("truncated_normal_library", "tf")),
                "--fse",
                "1" if save_eta else "0",
                "--profile",
                "1" if bool(sampler.get("profile", False)) else "0",
                "--fp",
                str(int(sampler.get("precision", 64))),
                "--eager",
                "1" if bool(sampler.get("eager", False)) else "0",
            ]
            chains = sampler.get("chains_to_run")
            if chains:
                cmd.extend(["--chains", *[str(int(x)) for x in chains]])
            code = _run_cmd(cmd, cwd=workdir, log_path=diagnostics / "HmscHPC_sample.log", env=env)
            if code != 0:
                raise RuntimeError("Hmsc-HPC sampler failed; see diagnostics/HmscHPC_sample.log")
            _write_outputs(workdir, cfg, posterior, warnings)
            status = "fitted"

    except Exception as exc:  # noqa: BLE001
        status = "fit_failed"
        errors.append(str(exc))
        (diagnostics / "HmscHPC_error.txt").write_text(
            "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
            encoding="utf-8",
        )
    finally:
        import pandas as pd

        for directory in [
            "inputs",
            "data",
            "models",
            "results",
            "tables",
            "plots",
            "diagnostics",
            "report",
            "predictions",
            "reproducible_script",
            "standard",
            "workflow_scripts",
            "samples",
        ]:
            (workdir / directory).mkdir(parents=True, exist_ok=True)
        _write_hmschpc_step_scripts(workdir, cfg)
        if status != "fitted":
            pd.DataFrame(
                [
                    {
                        "run_id": workdir.name,
                        "engine": "Hmsc-HPC",
                        "status": status,
                        "n_sites": None,
                        "n_responses": None,
                        "n_predictors": None,
                    }
                ]
            ).to_csv(workdir / "standard" / "run_summary.csv", index=False)
            pd.DataFrame(
                [
                    {"field": "engine", "value": "Hmsc-HPC"},
                    {"field": "status", "value": status},
                    {"field": "distribution", "value": cfg.get("model", {}).get("distribution")},
                    {"field": "random_mode", "value": cfg.get("random_effects", {}).get("mode")},
                ]
            ).to_csv(workdir / "results" / "S1_model_definition.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "step": "S2_fit_models",
                        "status": status,
                        "posterior_file": "samples/posterior.h5",
                        "samples": cfg.get("sampler", {}).get("samples"),
                        "transient": cfg.get("sampler", {}).get("transient"),
                        "thin": cfg.get("sampler", {}).get("thin"),
                        "chains": cfg.get("sampler", {}).get("chains"),
                    }
                ]
            ).to_csv(workdir / "results" / "S2_fit_models.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "engine": "Hmsc-HPC",
                        "parameter": "Beta",
                        "rhat_mean": None,
                        "rhat_max": None,
                        "ess_mean": None,
                        "ess_min": None,
                        "note": status,
                    }
                ]
            ).to_csv(workdir / "results" / "S3_convergence_summary.csv", index=False)
            pd.DataFrame(
                [{"engine": "Hmsc-HPC", "metric": "status", "response_id": None, "value": None, "notes": status}]
            ).to_csv(workdir / "results" / "S4_model_fit_summary.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "engine": "Hmsc-HPC",
                        "parameter": "Beta",
                        "rhat_mean": None,
                        "rhat_max": None,
                        "ess_mean": None,
                        "ess_min": None,
                        "note": status,
                    }
                ]
            ).to_csv(workdir / "standard" / "diagnostics_long.csv", index=False)
            pd.DataFrame(
                [{"engine": "Hmsc-HPC", "response_id": None, "observed_mean": None, "predicted_mean": None, "residual_mean": None, "notes": status}]
            ).to_csv(workdir / "results" / "S5_model_fit_prediction_summary.csv", index=False)
            pd.DataFrame(
                [{"engine": "Hmsc-HPC", "response_id": None, "predictor": None, "estimate": None, "notes": status}]
            ).to_csv(workdir / "results" / "S6_parameter_estimates_Beta.csv", index=False)
            pd.DataFrame(
                [{"engine": "Hmsc-HPC", "site_id": None, "response_id": None, "observed": None, "predicted_mean": None, "notes": status}]
            ).to_csv(workdir / "results" / "S7_predictions_training.csv", index=False)
            pd.DataFrame(
                [
                    {"result_group": "model_definition", "file": "results/S1_model_definition.csv"},
                    {"result_group": "fit", "file": "results/S2_fit_models.csv"},
                    {"result_group": "convergence", "file": "results/S3_convergence_summary.csv"},
                    {"result_group": "model_fit", "file": "results/S4_model_fit_summary.csv"},
                    {"result_group": "model_fit", "file": "results/S5_model_fit_prediction_summary.csv"},
                    {"result_group": "parameters", "file": "results/S6_parameter_estimates_Beta.csv"},
                    {"result_group": "predictions", "file": "results/S7_predictions_training.csv"},
                ]
            ).to_csv(workdir / "tables" / "HmscHPC_S1S7_result_index.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "level": "info",
                        "check": "workflow_status",
                        "ok": status == "model_defined",
                        "message": f"Hmsc-HPC workflow ended with status {status}.",
                    },
                    {
                        "level": "info",
                        "check": "distribution",
                        "ok": True,
                        "message": f"distribution={cfg.get('model', {}).get('distribution')}",
                    },
                    {
                        "level": "info",
                        "check": "random_effects",
                        "ok": True,
                        "message": f"random_mode={cfg.get('random_effects', {}).get('mode')}",
                    },
                ]
            ).to_csv(workdir / "diagnostics" / "data_check_messages.csv", index=False)
            pd.DataFrame(
                {"engine": ["Hmsc-HPC"], "response_id": [None], "predictor": [None], "direction": [None], "estimate": [None], "lower": [None], "upper": [None], "notes": [status]}
            ).to_csv(workdir / "standard" / "effects_long.csv", index=False)
            pd.DataFrame(
                {"engine": ["Hmsc-HPC"], "site_id": [None], "response_id": [None], "observed": [None], "predicted_mean": [None], "predicted_lower": [None], "predicted_upper": [None], "prediction_set": [status]}
            ).to_csv(workdir / "standard" / "predictions_long.csv", index=False)
            pd.DataFrame(
                {"engine": ["Hmsc-HPC"], "metric": [None], "response_id": [None], "value": [None], "notes": [status]}
            ).to_csv(workdir / "standard" / "fit_metrics.csv", index=False)
            pd.DataFrame(
                {"engine": ["Hmsc-HPC"], "response_1": [None], "response_2": [None], "association_type": [None], "estimate": [None], "comparable_level": [status]}
            ).to_csv(workdir / "standard" / "associations_long.csv", index=False)
            (workdir / "results" / "README_HmscHPC_results.txt").write_text(
                f"Hmsc-HPC workflow status: {status}\n"
                "Open diagnostics/engine_status.json and diagnostics/HmscHPC_error.txt if present.\n",
                encoding="utf-8",
            )
            (workdir / "report" / "Hmsc-HPC_report.html").write_text(
                f"<!DOCTYPE html><html><body><h1>Hmsc-HPC report</h1><p>Status: {status}</p></body></html>",
                encoding="utf-8",
            )
        try:
            _write_session_info(workdir / "diagnostics" / "session_info.txt", warnings)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"session_info export skipped: {exc}")
        payload = {
            "engine": "Hmsc-HPC",
            "status": status,
            "warnings": warnings,
            "errors": errors,
            "python": sys.executable,
        }
        _write_json(diagnostics / "engine_status.json", payload)
        import pandas as pd

        pd.DataFrame(
            [{"engine": "Hmsc-HPC", "status": status, "message": "; ".join(warnings + errors)}]
        ).to_csv(workdir / "tables" / "engine_status.csv", index=False)
        _fill_empty_dirs(workdir)
        _write_manifest(workdir)

    return 0 if status in {"fitted", "model_defined"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
