# -*- coding: utf-8 -*-
"""
Fixed-Effects Asymmetric Soil Temperature Response Analysis
-----------------------------------------------------------

This script performs a fixed-effects-based asymmetric response analysis of
soil temperature to snow-depth changes.

Workflow:
1. Load the pre-filtered panel dataset:
   - sampled_panel_filtered_with_landuse.parquet
   - or sampled_panel_filtered_with_landuse.csv.gz

2. Recalculate model coefficients for each ecosystem:
   - snow_pos_lag1
   - snow_neg_lag1
   - t2m
   - ssr
   - beta_pos, beta_neg, control-variable coefficients, p values, and rho_snow

3. Calculate binned soil temperature responses under equivalent snow-depth changes:
   - snow addition
   - snow reduction
   - mean, SE, 95% CI, p value, and significance stars

4. Generate final figures:
   - one combined 1 x 4 panel figure:
     Forest | Wetland | Grassland | Tundra
   - four separate asymmetry-ratio figures:
     forest_ratio_response.tif
     wetland_ratio_response.tif
     grassland_ratio_response.tif
     tundra_ratio_response.tif

5. Export summary tables:
   - model_coefficients_by_ecosystem_recomputed.csv
   - binned_response_4ecosystems.csv
   - ratio_response_4ecosystems.csv
   - rho_snow_4ecosystems.csv

Main model:
    delta_ts = beta_pos * snow_pos_lag1 + beta_neg * snow_neg_lag1
               + theta_t2m * t2m + theta_ssr * ssr
               + fixed effects + error

Asymmetry metric:
    rho_snow = |beta_neg| / |beta_pos|

Note:
- This script starts from an existing pre-filtered panel dataset.
- It does not resample the original TIFF files.
- The input panel must contain t2m and ssr columns.
- All paths and output names are written in English.
"""

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# 0) Paths and parameters
# =========================================================
# Put the input panel file in this directory, or modify this path.
RESULT_DIR = r"I:\ERA_TIFF\snow_depth_soil_temperature_ecosystem_asymmetry_final_outputs"

PANEL_PARQUET = os.path.join(RESULT_DIR, "sampled_panel_filtered_with_landuse.parquet")
PANEL_CSVGZ = os.path.join(RESULT_DIR, "sampled_panel_filtered_with_landuse.csv.gz")

# Final output folder
PLOT_OUT_DIR = os.path.join(RESULT_DIR, "fixed_effects_asymmetric_response_results")
os.makedirs(PLOT_OUT_DIR, exist_ok=True)

# Output file for recomputed model coefficients
COEF_OUT_CSV = os.path.join(PLOT_OUT_DIR, "model_coefficients_by_ecosystem_recomputed.csv")

# =========================================================
# 1) Unit and grouping settings
# =========================================================
# Use "m" if the original delta_snow_lag1 is in meters.
# Use "cm" if the original delta_snow_lag1 is already in centimeters.
SDE_UNIT = "m"

# Snow-depth amplitude bins up to 60 cm
AMP_CM = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60]

# Four ecosystems used in the final analysis
ECO_GROUPS = ["forest", "wetland", "grassland", "approx_tundra"]

GROUP_LABELS_EN = {
    "forest": "Forest",
    "wetland": "Wetland",
    "grassland": "Grassland",
    "approx_tundra": "Tundra"
}

OUTPUT_PREFIX = {
    "forest": "forest",
    "wetland": "wetland",
    "grassland": "grassland",
    "approx_tundra": "tundra"
}

# =========================================================
# 2) Model settings
# =========================================================
RECOMPUTE_MODEL = True

# Control variables included in the fixed-effects model
CONTROL_VARS = ["t2m", "ssr"]

# Fixed effects used in the model.
# Recommended setting: pixel_id + year + month.
# If the dataset is too large or slow to process, change this to ["year", "month"].
FIXED_EFFECTS = ["pixel_id", "year", "month"]

# Iteration settings for alternating demeaning
FE_MAX_ITER = 20
FE_TOL = 1e-7

# Minimum sample size required for model fitting
MIN_N_MODEL = 100

# Definition of snow_neg_lag1:
# True: snow_neg_lag1 = absolute magnitude of negative snow anomaly.
# This is recommended because beta_neg can then be interpreted as the effect of snow-reduction magnitude.
USE_NEGATIVE_MAGNITUDE = True

# =========================================================
# 3) Plot settings
# =========================================================
MAX_SCATTER_PER_SIGN = 25000
RANDOM_SEED = 42
DPI = 600

# Combined 1 x 4 figure
MAIN_FIG_W = 16.0
MAIN_FIG_H = 4.8

MAIN_YMIN = -10.0
MAIN_YMAX = 7.5

PLOT_XMIN = 0.0
PLOT_XMAX = 62.0

# Separate ratio figures
RATIO_FIG_W = 6.0
RATIO_FIG_H = 4.5

# Font settings
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "stix"

# ---------------------------------------------------------
# Main four-ecosystem figure settings
# ---------------------------------------------------------
# Ecosystem title font size
TITLE_SIZE = 14

# Common x/y axis title font size
LABEL_SIZE = 14

# x/y tick-label font size
TICK_SIZE = 12

# Significance star font size
STAR_SIZE = 12

# Axis frame and tick width for the main four-ecosystem figure
MAIN_SPINE_WIDTH = 0.6

# ---------------------------------------------------------
# Separate ratio figure settings
# ---------------------------------------------------------
RATIO_TITLE_SIZE = 20
RATIO_LABEL_SIZE = 18
RATIO_TICK_SIZE = 15
RATIO_LEGEND_SIZE = 13

# Keep ratio figure frame width unchanged
RATIO_SPINE_WIDTH = 1.2

# Colors
COLOR_INC = (255 / 255, 192 / 255, 127 / 255)
COLOR_DEC = (143 / 255, 196 / 255, 222 / 255)

# Ratio curve color
COLOR_RATIO_LINE = (70 / 255, 70 / 255, 70 / 255)

# rho_snow line color: RGB(70, 131, 180)
COLOR_RHO = (70 / 255, 131 / 255, 180 / 255)

COLOR_STAR = (70 / 255, 70 / 255, 70 / 255)

SHOW_STARS = True
SHOW_RATIO_LEGEND = True

# Fine-tune star positions to avoid clipping and overlap near 1 cm and 5 cm
STAR_X_SHIFT_FIRST = 1.15
STAR_X_SHIFT_SECOND = 0.35
STAR_Y_BASE_OFFSET = 0.030
STAR_Y_TIER_OFFSET = 0.050

# =========================================================
# 4) Basic utility functions
# =========================================================
def load_table_any(parquet_path=None, csv_gz_path=None):
    if parquet_path is not None and os.path.exists(parquet_path):
        return pd.read_parquet(parquet_path)

    if csv_gz_path is not None and os.path.exists(csv_gz_path):
        return pd.read_csv(csv_gz_path, compression="gzip")

    raise FileNotFoundError(
        f"Cannot find panel file:\n{parquet_path}\nor\n{csv_gz_path}"
    )


def load_panel():
    print(f"[INFO] Panel parquet: {PANEL_PARQUET}")
    print(f"[INFO] Panel csv.gz : {PANEL_CSVGZ}")

    panel = load_table_any(PANEL_PARQUET, PANEL_CSVGZ)

    return panel


def ensure_needed_columns(panel):
    needed_cols = [
        "pixel_id",
        "year",
        "month",
        "delta_snow_lag1",
        "delta_ts",
        "ecosystem"
    ] + CONTROL_VARS

    missing_cols = [c for c in needed_cols if c not in panel.columns]

    if missing_cols:
        raise RuntimeError(
            "The panel dataset is missing required columns: "
            f"{missing_cols}\n"
            "Please make sure the input panel contains t2m and ssr."
        )

    return panel


def snow_to_cm(x):
    x = pd.to_numeric(x, errors="coerce")

    if SDE_UNIT.lower() == "m":
        return x * 100.0

    return x


def star_from_p(p):
    if pd.isna(p):
        return ""

    if p < 0.001:
        return "***"

    if p < 0.01:
        return "**"

    if p < 0.05:
        return "*"

    return ""


def pvalue_from_t(t_val, df):
    if not np.isfinite(t_val) or df <= 0:
        return np.nan

    try:
        from scipy.stats import t as tdist
        p = 2 * (1 - tdist.cdf(abs(t_val), df=df))
        return float(p)
    except Exception:
        from math import erf, sqrt
        p = 2 * (1 - 0.5 * (1 + erf(abs(t_val) / sqrt(2))))
        return float(p)


def one_sample_pvalue(vals):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]

    n = len(vals)

    if n < 2:
        return np.nan

    mean = np.mean(vals)
    sd = np.std(vals, ddof=1)

    if (not np.isfinite(sd)) or sd == 0:
        return np.nan

    t_val = mean / (sd / math.sqrt(n))

    return pvalue_from_t(t_val, n - 1)


def make_bin_edges(centers):
    centers = np.array(sorted(centers), dtype=float)

    mids = (centers[:-1] + centers[1:]) / 2.0
    left = max(0.0, centers[0] - (centers[1] - centers[0]) / 2.0)
    right = centers[-1] + (centers[-1] - centers[-2]) / 2.0

    edges = np.concatenate([[left], mids, [right]])

    return edges


BIN_EDGES = make_bin_edges(AMP_CM)


def assign_amp_bin_cm(x_cm):
    idx = np.digitize(x_cm, BIN_EDGES, right=False) - 1
    idx = np.where((idx >= 0) & (idx < len(AMP_CM)), idx, -1)

    vals = np.array(AMP_CM, dtype=float)
    out = np.where(idx >= 0, vals[idx], np.nan)

    return out


# =========================================================
# 5) Prepare panel data for analysis
# =========================================================
def prepare_panel_for_analysis(panel):
    panel = panel.copy()
    panel = ensure_needed_columns(panel)

    panel = panel.replace([np.inf, -np.inf], np.nan)

    panel["delta_snow_lag1_cm"] = snow_to_cm(panel["delta_snow_lag1"])
    panel["amp_cm"] = np.abs(panel["delta_snow_lag1_cm"])

    panel["sign_group"] = np.where(
        panel["delta_snow_lag1_cm"] > 0,
        "increase",
        np.where(panel["delta_snow_lag1_cm"] < 0, "decrease", "zero")
    )

    panel["snow_pos_lag1"] = np.where(
        panel["delta_snow_lag1_cm"] > 0,
        panel["delta_snow_lag1_cm"],
        0.0
    )

    if USE_NEGATIVE_MAGNITUDE:
        panel["snow_neg_lag1"] = np.where(
            panel["delta_snow_lag1_cm"] < 0,
            np.abs(panel["delta_snow_lag1_cm"]),
            0.0
        )
    else:
        panel["snow_neg_lag1"] = np.where(
            panel["delta_snow_lag1_cm"] < 0,
            panel["delta_snow_lag1_cm"],
            0.0
        )

    for var in CONTROL_VARS:
        panel[var] = pd.to_numeric(panel[var], errors="coerce")

    panel["group_name"] = panel["ecosystem"].astype(str)

    panel = panel.dropna(
        subset=[
            "delta_snow_lag1_cm",
            "delta_ts",
            "amp_cm",
            "snow_pos_lag1",
            "snow_neg_lag1",
            "ecosystem"
        ] + CONTROL_VARS
    )

    panel = panel[panel["sign_group"].isin(["increase", "decrease"])].copy()
    panel = panel[panel["group_name"].isin(ECO_GROUPS)].copy()

    panel = panel[(panel["amp_cm"] >= 0) & (panel["amp_cm"] <= 60)].copy()

    return panel


# =========================================================
# 6) Fixed-effect residualization and OLS
# =========================================================
def residualize_by_fixed_effects(df, cols, fe_cols, max_iter=20, tol=1e-7):
    """
    Residualize y and X using alternating demeaning for multiple fixed effects.
    """
    if len(fe_cols) == 0:
        return df[cols].astype(float).copy()

    work = df[cols + fe_cols].copy()

    for c in cols:
        work[c] = pd.to_numeric(work[c], errors="coerce")

    work = work.dropna(subset=cols + fe_cols).copy()

    resid = work[cols].astype(float).copy()
    resid = resid - resid.mean(axis=0)

    last_ss = np.inf

    for _ in range(max_iter):
        for fe in fe_cols:
            means = resid.groupby(work[fe]).transform("mean")
            resid = resid - means

        ss = float(np.nansum(resid.values ** 2))

        if last_ss < np.inf:
            rel_change = abs(last_ss - ss) / max(last_ss, 1e-12)
            if rel_change < tol:
                break

        last_ss = ss

    resid.index = work.index

    return resid


def fit_ols_on_residuals(y, X):
    """
    Fit OLS on residualized y and X.
    No intercept is added because the variables have already been demeaned.
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)

    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[mask]
    X = X[mask, :]

    n = len(y)
    k = X.shape[1]

    if n <= k + 2:
        raise RuntimeError(f"Insufficient sample size for OLS: n={n}, k={k}")

    beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)

    y_hat = X @ beta
    e = y - y_hat

    df_resid = n - k
    sigma2 = float(np.sum(e ** 2) / df_resid)

    xtx_inv = np.linalg.pinv(X.T @ X)
    var_beta = sigma2 * xtx_inv
    se = np.sqrt(np.diag(var_beta))

    t_vals = beta / se
    p_vals = np.array([pvalue_from_t(t, df_resid) for t in t_vals], dtype=float)

    return {
        "coef": beta,
        "se": se,
        "t": t_vals,
        "p": p_vals,
        "n": n,
        "df_resid": df_resid,
        "rank": rank
    }


def fit_model_for_group(df_group, group_name):
    """
    Model:
    delta_ts = beta_pos * snow_pos_lag1
             + beta_neg * snow_neg_lag1
             + theta_t2m * t2m
             + theta_ssr * ssr
             + fixed effects
             + error

    beta_pos represents the marginal effect of snow addition after controlling for t2m, ssr,
    and multidimensional fixed effects.
    beta_neg represents the marginal effect of snow reduction after controlling for t2m, ssr,
    and multidimensional fixed effects.
    """
    df = df_group.copy()

    model_cols = ["delta_ts", "snow_pos_lag1", "snow_neg_lag1"] + CONTROL_VARS
    fe_cols = [c for c in FIXED_EFFECTS if c in df.columns]

    df = df.dropna(subset=model_cols + fe_cols).copy()

    if len(df) < MIN_N_MODEL:
        print(f"[SKIP MODEL] {group_name}: n={len(df)} < {MIN_N_MODEL}")
        return []

    print(f"[MODEL] {group_name}: n={len(df)}, controls={CONTROL_VARS}, fixed effects={fe_cols}")

    resid_df = residualize_by_fixed_effects(
        df=df,
        cols=model_cols,
        fe_cols=fe_cols,
        max_iter=FE_MAX_ITER,
        tol=FE_TOL
    )

    y = resid_df["delta_ts"].values
    x_vars = ["snow_pos_lag1", "snow_neg_lag1"] + CONTROL_VARS
    X = resid_df[x_vars].values

    fit = fit_ols_on_residuals(y, X)

    rows = []

    for j, var in enumerate(x_vars):
        rows.append({
            "group_name": group_name,
            "group_label": GROUP_LABELS_EN.get(group_name, group_name),
            "variable": var,
            "coef": float(fit["coef"][j]),
            "se": float(fit["se"][j]),
            "tvalue": float(fit["t"][j]),
            "pvalue": float(fit["p"][j]),
            "n": int(fit["n"]),
            "df_resid": int(fit["df_resid"]),
            "fixed_effects": "+".join(fe_cols),
            "controls": "+".join(CONTROL_VARS),
            "snow_unit": "cm",
            "negative_term": "magnitude" if USE_NEGATIVE_MAGNITUDE else "signed_negative"
        })

    return rows


def recompute_model_coefficients(plot_panel):
    all_rows = []

    for group_name in ECO_GROUPS:
        df_g = plot_panel[plot_panel["group_name"] == group_name].copy()

        if len(df_g) == 0:
            print(f"[SKIP MODEL] {group_name}: no rows")
            continue

        rows = fit_model_for_group(df_g, group_name)
        all_rows.extend(rows)

    coef_df = pd.DataFrame(all_rows)

    coef_df.to_csv(
        COEF_OUT_CSV,
        index=False,
        encoding="utf-8-sig"
    )

    print(f"[SAVE] Model coefficients -> {COEF_OUT_CSV}")

    return coef_df


# =========================================================
# 7) Binned response and ratio calculation
# =========================================================
def calc_binned_stats(df_group):
    df = df_group.copy()

    df = df[(df["amp_cm"] >= 0) & (df["amp_cm"] <= 60)].copy()
    df["amp_bin_cm"] = assign_amp_bin_cm(df["amp_cm"].values)
    df = df.dropna(subset=["amp_bin_cm"])

    rows = []

    for sign in ["increase", "decrease"]:
        d0 = df[df["sign_group"] == sign].copy()

        for a in AMP_CM:
            dd = d0[d0["amp_bin_cm"] == a].copy()

            y = dd["delta_ts"].values.astype(float)
            y = y[np.isfinite(y)]

            n = len(y)

            if n == 0:
                rows.append({
                    "sign_group": sign,
                    "amp_bin_cm": a,
                    "n": 0,
                    "mean": np.nan,
                    "se": np.nan,
                    "ci95": np.nan,
                    "pvalue": np.nan,
                    "stars": ""
                })
                continue

            mean = float(np.mean(y))
            se = float(np.std(y, ddof=1) / np.sqrt(n)) if n >= 2 else np.nan
            ci95 = float(1.96 * se) if np.isfinite(se) else np.nan
            pvalue = one_sample_pvalue(y)

            rows.append({
                "sign_group": sign,
                "amp_bin_cm": a,
                "n": int(n),
                "mean": mean,
                "se": se,
                "ci95": ci95,
                "pvalue": pvalue,
                "stars": star_from_p(pvalue)
            })

    return pd.DataFrame(rows)


def calc_observed_ratio_curve(bin_df):
    inc = bin_df[
        bin_df["sign_group"] == "increase"
    ][["amp_bin_cm", "mean", "ci95", "n"]].copy()

    dec = bin_df[
        bin_df["sign_group"] == "decrease"
    ][["amp_bin_cm", "mean", "ci95", "n"]].copy()

    inc = inc.rename(
        columns={
            "mean": "mean_inc",
            "ci95": "ci95_inc",
            "n": "n_inc"
        }
    )

    dec = dec.rename(
        columns={
            "mean": "mean_dec",
            "ci95": "ci95_dec",
            "n": "n_dec"
        }
    )

    rr = pd.merge(inc, dec, on="amp_bin_cm", how="outer")

    rr["ratio_obs"] = rr.apply(
        lambda r: np.nan if (
            pd.isna(r["mean_inc"]) or
            pd.isna(r["mean_dec"]) or
            abs(r["mean_inc"]) < 1e-12
        ) else abs(r["mean_dec"]) / abs(r["mean_inc"]),
        axis=1
    )

    rr = rr.sort_values("amp_bin_cm").reset_index(drop=True)

    return rr


def sample_scatter(df_group, max_points=25000, seed=42):
    rng = np.random.default_rng(seed)

    out_list = []

    for sign in ["increase", "decrease"]:
        dd = df_group[df_group["sign_group"] == sign].copy()

        dd = dd.dropna(subset=["amp_cm", "delta_ts"])
        dd = dd[(dd["amp_cm"] >= 0) & (dd["amp_cm"] <= 60)].copy()

        if len(dd) > max_points:
            idx = rng.choice(len(dd), size=max_points, replace=False)
            dd = dd.iloc[idx].copy()

        out_list.append(dd)

    return pd.concat(out_list, axis=0, ignore_index=True)


def get_coef_value(coef_df, group_name, variable, field):
    if coef_df is None or len(coef_df) == 0:
        return np.nan

    needed_cols = ["group_name", "variable", field]

    for c in needed_cols:
        if c not in coef_df.columns:
            return np.nan

    tmp = coef_df[
        (coef_df["group_name"] == group_name) &
        (coef_df["variable"] == variable)
    ]

    if len(tmp) == 0:
        return np.nan

    return float(tmp[field].iloc[0])


def calc_rho_snow_model(coef_df, group_name):
    beta_pos = get_coef_value(
        coef_df,
        group_name,
        "snow_pos_lag1",
        "coef"
    )

    beta_neg = get_coef_value(
        coef_df,
        group_name,
        "snow_neg_lag1",
        "coef"
    )

    p_pos = get_coef_value(
        coef_df,
        group_name,
        "snow_pos_lag1",
        "pvalue"
    )

    p_neg = get_coef_value(
        coef_df,
        group_name,
        "snow_neg_lag1",
        "pvalue"
    )

    if np.isfinite(beta_pos) and np.isfinite(beta_neg) and abs(beta_pos) > 1e-12:
        rho = abs(beta_neg) / abs(beta_pos)
    else:
        rho = np.nan

    return {
        "group_name": group_name,
        "beta_pos": beta_pos,
        "beta_neg": beta_neg,
        "p_pos": p_pos,
        "p_neg": p_neg,
        "rho_snow": rho
    }


# =========================================================
# 8) Plot style functions
# =========================================================
def set_axes_style(ax, tick_size=14, spine_width=1.2):
    for spine in ax.spines.values():
        spine.set_linewidth(spine_width)

    ax.tick_params(
        axis="both",
        labelsize=tick_size,
        width=spine_width,
        pad=5
    )

    for lbl in ax.get_xticklabels():
        lbl.set_fontweight("bold")
        lbl.set_fontfamily("Times New Roman")

    for lbl in ax.get_yticklabels():
        lbl.set_fontweight("bold")
        lbl.set_fontfamily("Times New Roman")


def set_legend_bold(leg):
    if leg is None:
        return

    for txt in leg.get_texts():
        txt.set_fontweight("bold")
        txt.set_fontfamily("Times New Roman")


def add_star_text(ax, x, y, star, above=True, tier=0):
    """
    Add significance stars with adjusted x/y offsets.
    This version reduces clipping at the first snow-depth bin
    and avoids overlap between the first and second bins.
    """
    if (not SHOW_STARS) or star is None or star == "" or (not np.isfinite(y)):
        return

    y0, y1 = ax.get_ylim()
    yr = y1 - y0

    offset = (STAR_Y_BASE_OFFSET + STAR_Y_TIER_OFFSET * tier) * yr

    if np.isclose(x, 1.0):
        x_text = x + STAR_X_SHIFT_FIRST
    elif np.isclose(x, 5.0):
        x_text = x + STAR_X_SHIFT_SECOND
    else:
        x_text = x

    yy = y + offset if above else y - offset

    margin = 0.040 * yr
    yy = min(max(yy, y0 + margin), y1 - margin)

    va = "bottom" if above else "top"

    ax.text(
        x_text,
        yy,
        star,
        ha="center",
        va=va,
        fontsize=STAR_SIZE,
        color=COLOR_STAR,
        fontweight="bold",
        fontfamily="Times New Roman",
        zorder=10,
        clip_on=False
    )


def add_stars_without_overlap(ax, bin_df):
    """
    Add significance stars and reduce overlap across adjacent bins.
    This is especially useful for the 1 cm and 5 cm bins.
    """
    if not SHOW_STARS:
        return

    star_df = bin_df.copy()
    star_df = star_df[star_df["stars"].astype(str) != ""].copy()
    star_df = star_df.dropna(subset=["mean", "amp_bin_cm"])

    if len(star_df) == 0:
        return

    star_df = star_df.sort_values(
        ["amp_bin_cm", "sign_group", "mean"]
    ).reset_index(drop=True)

    upper_tier_near_left = 0
    lower_tier_near_left = 0

    for x_val, dfx in star_df.groupby("amp_bin_cm"):
        dfx = dfx.copy().sort_values(["sign_group", "mean"]).reset_index(drop=True)

        above_count = 0
        below_count = 0

        for _, row in dfx.iterrows():
            y = row["mean"]
            sign = row["sign_group"]
            star = row["stars"]

            above = True if y >= 0 else False

            if (sign == "decrease") and np.isclose(x_val, 1.0):
                above = False

            if np.isclose(x_val, 1.0) or np.isclose(x_val, 5.0):
                if above:
                    tier = upper_tier_near_left
                    upper_tier_near_left += 1
                else:
                    tier = lower_tier_near_left
                    lower_tier_near_left += 1
            else:
                if above:
                    tier = above_count
                    above_count += 1
                else:
                    tier = below_count
                    below_count += 1

            add_star_text(
                ax=ax,
                x=x_val,
                y=y,
                star=star,
                above=above,
                tier=tier
            )


# =========================================================
# 9) Draw combined four-ecosystem response figure
# =========================================================
def draw_four_ecosystems_response(plot_panel, out_path):
    fig, axes = plt.subplots(
        nrows=1,
        ncols=4,
        figsize=(MAIN_FIG_W, MAIN_FIG_H),
        sharey=True
    )

    all_bin_rows = []

    for i, group_name in enumerate(ECO_GROUPS):
        ax = axes[i]
        group_label = GROUP_LABELS_EN[group_name]

        df_g = plot_panel[plot_panel["group_name"] == group_name].copy()

        if len(df_g) == 0:
            ax.set_title(
                group_label,
                fontsize=TITLE_SIZE,
                fontweight="bold",
                fontfamily="Times New Roman"
            )

            ax.text(
                0.5,
                0.5,
                "No data",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=12,
                fontweight="bold",
                fontfamily="Times New Roman"
            )

            ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
            ax.set_ylim(MAIN_YMIN, MAIN_YMAX)
            ax.set_yticks(np.arange(MAIN_YMIN, MAIN_YMAX + 0.001, 5))

            set_axes_style(ax, tick_size=TICK_SIZE, spine_width=MAIN_SPINE_WIDTH)

            if i > 0:
                ax.tick_params(axis="y", left=False, labelleft=False)

            continue

        scatter_df = sample_scatter(
            df_g,
            max_points=MAX_SCATTER_PER_SIGN,
            seed=RANDOM_SEED
        )

        bin_df = calc_binned_stats(df_g)
        bin_df["group_name"] = group_name
        bin_df["group_label"] = group_label
        all_bin_rows.append(bin_df)

        sc_inc = scatter_df[scatter_df["sign_group"] == "increase"]
        sc_dec = scatter_df[scatter_df["sign_group"] == "decrease"]

        ax.scatter(
            sc_inc["amp_cm"],
            sc_inc["delta_ts"],
            s=6,
            alpha=0.10,
            color=COLOR_INC,
            edgecolors="none",
            zorder=1
        )

        ax.scatter(
            sc_dec["amp_cm"],
            sc_dec["delta_ts"],
            s=6,
            alpha=0.10,
            color=COLOR_DEC,
            edgecolors="none",
            zorder=1
        )

        for sign, color in [
            ("increase", COLOR_INC),
            ("decrease", COLOR_DEC)
        ]:
            dd = bin_df[
                bin_df["sign_group"] == sign
            ].copy().sort_values("amp_bin_cm")

            x = dd["amp_bin_cm"].values
            y = dd["mean"].values
            e = dd["ci95"].values

            ax.plot(
                x,
                y,
                color=color,
                linewidth=2.2,
                marker="o",
                markersize=4.8,
                zorder=5
            )

            if np.isfinite(e).any():
                y_low = y - np.where(np.isfinite(e), e, 0)
                y_high = y + np.where(np.isfinite(e), e, 0)

                ax.fill_between(
                    x,
                    y_low,
                    y_high,
                    color=color,
                    alpha=0.20,
                    zorder=4
                )

        ax.axhline(
            0,
            linestyle="--",
            linewidth=0.8,
            color="#777777",
            zorder=2
        )

        ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
        ax.set_ylim(MAIN_YMIN, MAIN_YMAX)
        ax.set_xticks([0, 10, 20, 30, 40, 50, 60])
        ax.set_yticks(np.arange(MAIN_YMIN, MAIN_YMAX + 0.001, 5))

        ax.set_title(
            group_label,
            fontsize=TITLE_SIZE,
            fontweight="bold",
            fontfamily="Times New Roman",
            pad=8
        )

        set_axes_style(ax, tick_size=TICK_SIZE, spine_width=MAIN_SPINE_WIDTH)

        add_stars_without_overlap(ax, bin_df)

        if i > 0:
            ax.tick_params(axis="y", left=False, labelleft=False)

        leg = ax.get_legend()

        if leg is not None:
            leg.remove()

    fig.text(
        0.5,
        0.055,
        "Snow depth (cm)",
        ha="center",
        va="center",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        fontfamily="Times New Roman"
    )

    fig.text(
        0.032,
        0.52,
        "Soil temperature (°C)",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=LABEL_SIZE,
        fontweight="bold",
        fontfamily="Times New Roman"
    )

    fig.subplots_adjust(
        left=0.070,
        right=0.995,
        bottom=0.16,
        top=0.86,
        wspace=0.035
    )

    fig.savefig(
        out_path,
        dpi=DPI,
        facecolor="white",
        format="tif"
    )

    plt.close(fig)

    if len(all_bin_rows) > 0:
        all_bin_df = pd.concat(all_bin_rows, axis=0, ignore_index=True)
    else:
        all_bin_df = pd.DataFrame()

    return all_bin_df


# =========================================================
# 10) Draw separate ratio figures
# =========================================================
def draw_single_ratio_figure(ratio_df, rho_info, group_name, out_path):
    group_label = GROUP_LABELS_EN[group_name]

    fig, ax = plt.subplots(figsize=(RATIO_FIG_W, RATIO_FIG_H))

    rr = ratio_df.copy()

    ax.plot(
        rr["amp_bin_cm"],
        rr["ratio_obs"],
        color=COLOR_RATIO_LINE,
        linewidth=2.2,
        marker="o",
        markersize=5.2,
        label=r"$|Effect_{reduction}| / |Effect_{addition}|$",
        zorder=5
    )

    ax.axhline(
        1.0,
        linestyle="--",
        linewidth=1.2,
        color="#888888",
        label="ratio = 1",
        zorder=2
    )

    if np.isfinite(rho_info["rho_snow"]):
        ax.axhline(
            rho_info["rho_snow"],
            linestyle="-.",
            linewidth=1.5,
            color=COLOR_RHO,
            label=rf"$\rho_{{snow}}$ = {rho_info['rho_snow']:.2f}",
            zorder=3
        )

    ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
    ax.set_xticks(AMP_CM)
    ax.set_xticklabels([str(v) for v in AMP_CM])

    ratio_valid = rr["ratio_obs"].values
    ratio_valid = ratio_valid[np.isfinite(ratio_valid)]

    if len(ratio_valid) > 0:
        ymax = max(2.5, float(np.nanmax(ratio_valid)) * 1.35)

        if np.isfinite(rho_info["rho_snow"]):
            ymax = max(ymax, rho_info["rho_snow"] * 1.25)

        ymax = min(ymax, 12.0)
        ax.set_ylim(0, ymax)
    else:
        ax.set_ylim(0, 2.5)

    ax.set_title(
        group_label,
        fontsize=RATIO_TITLE_SIZE,
        fontweight="bold",
        fontfamily="Times New Roman",
        pad=8
    )

    ax.set_xlabel(
        "Snow depth (cm)",
        fontsize=RATIO_LABEL_SIZE,
        fontweight="bold",
        fontfamily="Times New Roman",
        labelpad=8
    )

    ax.set_ylabel(
        "Asymmetry ratio",
        fontsize=RATIO_LABEL_SIZE,
        fontweight="bold",
        fontfamily="Times New Roman",
        labelpad=8
    )

    set_axes_style(
        ax,
        tick_size=RATIO_TICK_SIZE,
        spine_width=RATIO_SPINE_WIDTH
    )

    if SHOW_RATIO_LEGEND:
        leg = ax.legend(
            loc="upper left",
            fontsize=RATIO_LEGEND_SIZE,
            frameon=False,
            handlelength=1.8,
            labelspacing=0.6
        )

        set_legend_bold(leg)

    fig.subplots_adjust(
        left=0.16,
        right=0.97,
        bottom=0.16,
        top=0.88
    )

    fig.savefig(
        out_path,
        dpi=DPI,
        facecolor="white",
        format="tif"
    )

    plt.close(fig)


def draw_all_ratio_figures(plot_panel, coef_df):
    all_ratio_rows = []
    rho_rows = []

    for group_name in ECO_GROUPS:
        group_label = GROUP_LABELS_EN[group_name]

        df_g = plot_panel[
            plot_panel["group_name"] == group_name
        ].copy()

        if len(df_g) == 0:
            print(f"[SKIP RATIO] {group_label}: no rows")
            continue

        bin_df = calc_binned_stats(df_g)
        ratio_df = calc_observed_ratio_curve(bin_df)
        rho_info = calc_rho_snow_model(coef_df, group_name)

        ratio_df["group_name"] = group_name
        ratio_df["group_label"] = group_label
        all_ratio_rows.append(ratio_df)

        rho_rows.append({
            "group_name": group_name,
            "group_label": group_label,
            "beta_pos": rho_info["beta_pos"],
            "beta_neg": rho_info["beta_neg"],
            "p_pos": rho_info["p_pos"],
            "p_neg": rho_info["p_neg"],
            "rho_snow_model": rho_info["rho_snow"]
        })

        out_name = f"{OUTPUT_PREFIX[group_name]}_ratio_response.tif"
        out_path = os.path.join(PLOT_OUT_DIR, out_name)

        print(f"[DRAW RATIO] {group_label} -> {out_path}")

        draw_single_ratio_figure(
            ratio_df=ratio_df,
            rho_info=rho_info,
            group_name=group_name,
            out_path=out_path
        )

    if len(all_ratio_rows) > 0:
        all_ratio_df = pd.concat(
            all_ratio_rows,
            axis=0,
            ignore_index=True
        )
    else:
        all_ratio_df = pd.DataFrame()

    if len(rho_rows) > 0:
        rho_df = pd.DataFrame(rho_rows)
    else:
        rho_df = pd.DataFrame()

    return all_ratio_df, rho_df


# =========================================================
# 11) Main workflow
# =========================================================
def main():
    print("=" * 90)
    print("Fixed-effects asymmetric soil temperature response analysis")
    print("=" * 90)

    panel = load_panel()
    plot_panel = prepare_panel_for_analysis(panel)

    print(f"[INFO] Valid rows for analysis: {len(plot_panel):,}")
    print(f"[INFO] Control variables used in the model: {CONTROL_VARS}")
    print(f"[INFO] Fixed effects used in the model: {FIXED_EFFECTS}")

    if RECOMPUTE_MODEL:
        coef_df = recompute_model_coefficients(plot_panel)
    else:
        if not os.path.exists(COEF_OUT_CSV):
            raise FileNotFoundError(
                f"RECOMPUTE_MODEL=False, but the coefficient file was not found:\n{COEF_OUT_CSV}"
            )

        coef_df = pd.read_csv(COEF_OUT_CSV, encoding="utf-8-sig")

    main_out = os.path.join(
        PLOT_OUT_DIR,
        "four_ecosystems_soil_temperature_response.tif"
    )

    print(f"[DRAW MAIN] Four-ecosystem response figure -> {main_out}")

    all_bin_df = draw_four_ecosystems_response(
        plot_panel=plot_panel,
        out_path=main_out
    )

    all_ratio_df, rho_df = draw_all_ratio_figures(
        plot_panel=plot_panel,
        coef_df=coef_df
    )

    if len(all_bin_df) > 0:
        out_csv = os.path.join(PLOT_OUT_DIR, "binned_response_4ecosystems.csv")
        all_bin_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
        print(f"[SAVE] {out_csv}")

    if len(all_ratio_df) > 0:
        out_csv = os.path.join(PLOT_OUT_DIR, "ratio_response_4ecosystems.csv")
        all_ratio_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
        print(f"[SAVE] {out_csv}")

    if len(rho_df) > 0:
        out_csv = os.path.join(PLOT_OUT_DIR, "rho_snow_4ecosystems.csv")
        rho_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
        print(f"[SAVE] {out_csv}")

    print("=" * 90)
    print("All tasks completed")
    print("=" * 90)
    print(f"Output directory: {PLOT_OUT_DIR}")
    print("Main output files:")
    print(" - model_coefficients_by_ecosystem_recomputed.csv")
    print(" - four_ecosystems_soil_temperature_response.tif")
    print(" - forest_ratio_response.tif")
    print(" - wetland_ratio_response.tif")
    print(" - grassland_ratio_response.tif")
    print(" - tundra_ratio_response.tif")
    print(" - binned_response_4ecosystems.csv")
    print(" - ratio_response_4ecosystems.csv")
    print(" - rho_snow_4ecosystems.csv")


if __name__ == "__main__":
    main()
