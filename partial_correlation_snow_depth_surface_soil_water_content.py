# -*- coding: utf-8 -*-
"""
Integrated full English script:
Partial correlation analysis between snow depth and surface soil water content
(pixel-wise workflow, including calculation, raster/table/scatter outputs,
and the final publication-style large-scale pixel correlation map).

Main features
-------------
1. Performs the original pixel-wise partial-correlation analysis.
2. Saves raster outputs, tables, scatter-plot outputs, and diagnostic outputs.
3. Generates the final publication-style large-scale pixel correlation map.
4. All comments, labels, and newly modified sections are provided in English.

Key map-plot updates
--------------------
- Remove background grid lines.
- Remove longitude/latitude coordinates, labels, and tick marks.
- Use Arial for the correlation legend.
- Use the requested diverging color scheme:
    negative correlation: RGB(32, 56, 136)
    positive correlation: RGB(225, 156, 102)
- Use the real world-boundary shapefile path:
    I:\\<world_map_folder>\\<world_map_folder>\\global_all_country.shp

Notes
-----
- This script preserves the original full analytical workflow as much as possible.
- The large-scale publication map is generated from the existing / newly created
  partial-correlation raster result.
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import math
import warnings
from multiprocessing import freeze_support
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.transform import Affine

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties

import geopandas as gpd

from scipy import stats
from scipy.spatial.distance import cdist
from tqdm import tqdm

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =========================================================
# 1. Basic parameters
# =========================================================
OUT_ROOT = r"I:\ERA_TIFF\!partial_corr_swvl1_core1_outputs_v1_scatter"
OUT_RASTER = os.path.join(OUT_ROOT, "01_rasters")
OUT_FIG = os.path.join(OUT_ROOT, "02_figures")
OUT_TABLE = os.path.join(OUT_ROOT, "03_tables")
OUT_TEXT = os.path.join(OUT_ROOT, "04_texts")
OUT_SCATTER_FIG = os.path.join(OUT_ROOT, "05_scatter_figures")
OUT_SCATTER_TABLE = os.path.join(OUT_ROOT, "06_scatter_tables")
OUT_PUB_MAP = os.path.join(OUT_ROOT, "07_publication_maps")

for d in [OUT_ROOT, OUT_RASTER, OUT_FIG, OUT_TABLE, OUT_TEXT, OUT_SCATTER_FIG, OUT_SCATTER_TABLE, OUT_PUB_MAP]:
    os.makedirs(d, exist_ok=True)

RUN_LOG = os.path.join(OUT_ROOT, "run_log.txt")
if os.path.exists(RUN_LOG):
    os.remove(RUN_LOG)

YEARS = list(range(2000, 2026))
MONTHS = list(range(1, 13))
TIME_STEPS = [(y, m) for y in YEARS for m in MONTHS]
MONTH_IDS = np.array([m for _, m in TIME_STEPS], dtype=np.int16)

TARGET_VAR = "swvl1"
TARGET_DESC = "Surface soil water content (volumetric soil water layer 1)"
X_VAR = "sde"

TARGET_DISPLAY_NAME = "surface soil water content"
X_DISPLAY_NAME = "snow depth"
ANALYSIS_TITLE = "Partial correlation analysis between snow depth and surface soil water content"
OUTPUT_PREFIX = "swvl1"

COVARS = [
    "t2m",
    "snowc",
    "sf",
    "smlt",
    "slhf",
    "v10",
    "lai_hv"
]

NODATA_OUT = -9999.0
COMPRESS = "LZW"
ALPHA = 0.05

REMOVE_MONTHLY_CLIM = True

# Domain definition: retain only pixels with snow-month frequency >= 0.30
SNOWC_MIN_FRACTION = 0.30
SNOWC_THRESHOLD = 0.0

N_JOBS = 4
BLOCK_ROWS = 8
PIXEL_CHUNK = 2000

SIG_POINT_STRIDE = 8
SIG_POINT_SIZE = 1.2

R_VMIN, R_VMAX = -1.0, 1.0
P_VMIN, P_VMAX = 0.0, 0.05

CHECK_ALL_FILES = True
CHECK_ALIGNMENT = True

# Number of displayed scatter points: only plotting is sampled,
# whereas statistics use all samples
PLOT_POINTS_PER_BLOCK_PER_GROUP = 1200
MAX_PLOT_POINTS_PER_GROUP = 30000
SCATTER_RANDOM_STATE = 42

# =========================================================
# 2. Plotting style parameters
# =========================================================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.unicode_minus"] = False

ARIAL_FONT = FontProperties(family="Arial")

# Scatter plot style
SCATTER_FIGSIZE = (10, 7.5)   # 4:3
SCATTER_SPINE_WIDTH = 1.0
SCATTER_LABEL_SIZE = 26
SCATTER_TICK_SIZE = 20
SCATTER_TITLE_SIZE = 20
SCATTER_TEXT_SIZE = 19

SCATTER_FACE = (145/255, 171/255, 210/255)
SCATTER_EDGE = (145/255, 171/255, 210/255)
REG_LINE = (175/255, 175/255, 175/255)

# Standard pixel-map style
MAP_FIGSIZE = (10, 6)
MAP_SPINE_WIDTH = 1.0
MAP_TITLE_SIZE = 15
MAP_CBAR_LABEL_SIZE = 12
MAP_CBAR_TICK_SIZE = 10

# Publication-style world map
FIG_DPI = 600
PUB_BG_RGB = (243 / 255, 243 / 255, 243 / 255)
PUB_EDGE_RGB = (175 / 255, 175 / 255, 175 / 255)

# Requested endpoint colors for the final correlation map
NEGATIVE_RGB = (32 / 255, 56 / 255, 136 / 255)
ZERO_RGB = (1.0, 1.0, 1.0)
POSITIVE_RGB = (225 / 255, 156 / 255, 102 / 255)

CORR_CMAP = LinearSegmentedColormap.from_list(
    "custom_partial_correlation",
    [NEGATIVE_RGB, ZERO_RGB, POSITIVE_RGB],
    N=256
)
CORR_CMAP.set_bad((1.0, 1.0, 1.0, 0.0))

# =========================================================
# 3. Path configuration
# =========================================================
ERA_BASE_DIR = r"I:\ERA_TIFF"

ERA_STYLE_VARS = {
    "swvl1": "swvl1",
    "sde": "sde",
    "snowc": "snowc",
    "sf": "sf",
    "smlt": "smlt",
    "lai_hv": "lai_hv",
}

SPECIAL_MONTHLY_VARS = {
    "t2m": {
        "candidates": [
            {
                "dir": r"I:\new_factor_tif1\de34709100aa97f6f68a02131b8fde57\t2m\t2m",
                "prefix": "de34709100aa97f6f68a02131b8fde57_t2m_"
            }
        ]
    },
    "v10": {
        "candidates": [
            {
                "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\v10\v10",
                "prefix": "9ce01cf43be7606676d84be71fe19678_v10_"
            }
        ]
    },
    "slhf": {
        "candidates": [
            {
                "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\slhf\slhf",
                "prefix": "9ce01cf43be7606676d84be71fe19678_slhf_"
            },
            {
                "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\slhf\sIhf",
                "prefix": "9ce01cf43be7606676d84be71fe19678_slhf_"
            }
        ]
    }
}

# World boundary shapefile for the final large-scale publication map
WORLD_SHP = "I:\\\u4e16\u754c\u5730\u56fe\\\u4e16\u754c\u5730\u56fe\\global_all_country.shp"

# =========================================================
# 4. Template global variables
# =========================================================
TEMPLATE_CRS = None
TEMPLATE_TRANSFORM = None
TEMPLATE_WIDTH = None
TEMPLATE_HEIGHT = None


def log_msg(msg):
    print(msg)
    with open(RUN_LOG, "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")


def ensure_file_exists(path, label="file"):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} does not exist: {path}")


def ensure_file_saved(path, label="file"):
    if not os.path.exists(path):
        raise IOError(f"{label} was not written successfully: {path}")
    if os.path.isdir(path):
        return
    if os.path.getsize(path) == 0:
        raise IOError(f"{label} is empty: {path}")


def save_df_csv(df, out_csv):
    df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    ensure_file_saved(out_csv, "CSV")


def save_txt(text, out_txt):
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text)
    ensure_file_saved(out_txt, "TXT")


def init_worker_template(crs, transform, width, height):
    global TEMPLATE_CRS, TEMPLATE_TRANSFORM, TEMPLATE_WIDTH, TEMPLATE_HEIGHT
    TEMPLATE_CRS = crs
    TEMPLATE_TRANSFORM = transform
    TEMPLATE_WIDTH = width
    TEMPLATE_HEIGHT = height


# =========================================================
# 5. Basic path and reading functions
# =========================================================
def _special_candidate_paths(var_name, year, month):
    spec = SPECIAL_MONTHLY_VARS[var_name]
    ymd = f"{year}{month:02d}01"
    paths = []
    for cand in spec["candidates"]:
        paths.append(os.path.join(cand["dir"], f"{cand['prefix']}{ymd}.tif"))
    return paths


def tif_path(var_name, year, month):
    if var_name in SPECIAL_MONTHLY_VARS:
        candidates = _special_candidate_paths(var_name, year, month)
        for p in candidates:
            if os.path.exists(p):
                return p
        return candidates[0]

    if var_name in ERA_STYLE_VARS:
        return os.path.join(ERA_BASE_DIR, ERA_STYLE_VARS[var_name], f"{year}_{month:02d}.tif")

    raise KeyError(f"Unknown variable: {var_name}")


def get_template_info():
    p = tif_path(TARGET_VAR, YEARS[0], MONTHS[0])
    if not os.path.exists(p):
        raise FileNotFoundError(f"Template file does not exist: {p}")

    with rasterio.open(p) as ds:
        info = {
            "path": p,
            "width": ds.width,
            "height": ds.height,
            "transform": ds.transform,
            "crs": ds.crs,
            "count": ds.count,
            "profile": ds.profile.copy()
        }
    return info


def check_inputs(template_info):
    all_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))

    if CHECK_ALL_FILES:
        missing = []
        for var in all_vars:
            for year, month in TIME_STEPS:
                p = tif_path(var, year, month)
                if not os.path.exists(p):
                    missing.append(p)

        if missing:
            print("The following files are missing (showing only the first 30):")
            for p in missing[:30]:
                print("  ", p)
            raise FileNotFoundError(f"A total of {len(missing)} missing files were found. Please complete them first.")

    if CHECK_ALIGNMENT:
        mismatched = []
        for var in all_vars:
            p = tif_path(var, YEARS[0], MONTHS[0])
            with rasterio.open(p) as ds:
                if ds.count != 1:
                    raise ValueError(f"{var} is not a single-band raster: {p}")

                same_grid = (
                    ds.width == template_info["width"] and
                    ds.height == template_info["height"] and
                    ds.transform == template_info["transform"] and
                    ds.crs == template_info["crs"]
                )
                if not same_grid:
                    mismatched.append((var, p))

        if mismatched:
            print(f"[Notice] {len(mismatched)} variables do not match the template grid and will be automatically aligned to the template during reading.")
            for var, p in mismatched[:10]:
                print(f"  Auto-aligning: {var} -> {p}")


def block_ranges(height, block_rows):
    ranges = []
    r0 = 0
    while r0 < height:
        r1 = min(r0 + block_rows, height)
        ranges.append((r0, r1))
        r0 = r1
    return ranges


def get_resampling_method(var_name):
    return Resampling.nearest if var_name == "snowc" else Resampling.bilinear


def read_one_block(var_name, path, row_start, row_end, width):
    h = row_end - row_start
    win = Window(col_off=0, row_off=row_start, width=width, height=h)
    resampling = get_resampling_method(var_name)

    with rasterio.open(path) as ds:
        same_grid = (
            ds.width == TEMPLATE_WIDTH and
            ds.height == TEMPLATE_HEIGHT and
            ds.transform == TEMPLATE_TRANSFORM and
            ds.crs == TEMPLATE_CRS
        )

        if same_grid:
            arr = ds.read(1, window=win).astype(np.float32)
            nd = ds.nodata
            if nd is not None:
                arr[arr == nd] = np.nan
            return arr

        with WarpedVRT(
            ds,
            crs=TEMPLATE_CRS,
            transform=TEMPLATE_TRANSFORM,
            width=TEMPLATE_WIDTH,
            height=TEMPLATE_HEIGHT,
            resampling=resampling
        ) as vrt:
            arr = vrt.read(1, window=win).astype(np.float32)
            nd = vrt.nodata if vrt.nodata is not None else ds.nodata
            if nd is not None:
                arr[arr == nd] = np.nan
            return arr


def load_var_cube(var_name, row_start, row_end, width):
    h = row_end - row_start
    cube = np.empty((len(TIME_STEPS), h, width), dtype=np.float32)
    for i, (year, month) in enumerate(TIME_STEPS):
        cube[i] = read_one_block(var_name, tif_path(var_name, year, month), row_start, row_end, width)
    return cube


# =========================================================
# 6. Time-series processing
# =========================================================
def remove_monthly_climatology_3d(cube):
    out = np.empty_like(cube, dtype=np.float32)
    for m in range(1, 13):
        idx = np.where(MONTH_IDS == m)[0]
        clim = np.nanmean(cube[idx], axis=0)
        out[idx] = cube[idx] - clim
    return out


def remove_monthly_climatology_1d(series):
    out = np.empty_like(series, dtype=np.float64)
    for m in range(1, 13):
        idx = np.where(MONTH_IDS == m)[0]
        clim = np.nanmean(series[idx])
        out[idx] = series[idx] - clim
    return out


def effective_sample_size(series):
    s = np.asarray(series, dtype=np.float64)
    s = s[np.isfinite(s)]

    n = len(s)
    if n < 10:
        return float(n)

    s_center = s - np.mean(s)
    x0 = s_center[:-1]
    x1 = s_center[1:]
    den = np.sum(x0 * x0)
    num = np.sum(x0 * x1)

    rho1 = 0.0 if den <= 0 else num / den
    rho1 = np.clip(rho1, -0.99, 0.99)

    n_eff = n * (1 - rho1) / (1 + rho1)
    return float(max(10, min(n, n_eff)))


# =========================================================
# 7. Core partial-correlation calculation
# =========================================================
def partial_corr_and_p_batch_with_eff_df(x, y, z, pixel_chunk=2000):
    T, N = x.shape
    K = z.shape[1]

    r_out = np.full(N, np.nan, dtype=np.float32)
    p_out = np.full(N, np.nan, dtype=np.float32)
    eff_n_out = np.full(N, np.nan, dtype=np.float32)

    valid = (
        np.all(np.isfinite(x), axis=0) &
        np.all(np.isfinite(y), axis=0) &
        np.all(np.isfinite(z), axis=(0, 1))
    )

    idx = np.where(valid)[0]
    if idx.size == 0:
        return r_out, p_out, eff_n_out

    n_chunks = math.ceil(idx.size / pixel_chunk)

    for i in range(n_chunks):
        sub = idx[i * pixel_chunk:(i + 1) * pixel_chunk]
        if sub.size == 0:
            continue

        x0 = x[:, sub].T.astype(np.float64)
        y0 = y[:, sub].T.astype(np.float64)
        z0 = np.transpose(z[:, :, sub], (2, 0, 1)).astype(np.float64)

        ones = np.ones((z0.shape[0], z0.shape[1], 1), dtype=np.float64)
        A = np.concatenate([ones, z0], axis=2)

        AtA = np.einsum("ntp,ntq->npq", A, A)
        AtA_inv = np.linalg.pinv(AtA)

        Atx = np.einsum("ntp,nt->np", A, x0)
        Aty = np.einsum("ntp,nt->np", A, y0)

        beta_x = np.einsum("npq,nq->np", AtA_inv, Atx)
        beta_y = np.einsum("npq,nq->np", AtA_inv, Aty)

        rx = x0 - np.einsum("ntp,np->nt", A, beta_x)
        ry = y0 - np.einsum("ntp,np->nt", A, beta_y)

        rx = rx - rx.mean(axis=1, keepdims=True)
        ry = ry - ry.mean(axis=1, keepdims=True)

        num = np.sum(rx * ry, axis=1)
        den = np.sqrt(np.sum(rx * rx, axis=1) * np.sum(ry * ry, axis=1))

        r = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
        r = np.clip(r, -0.999999, 0.999999)

        for j, pixel_idx in enumerate(sub):
            n_eff_x = effective_sample_size(x[:, pixel_idx])
            n_eff_y = effective_sample_size(y[:, pixel_idx])
            n_eff = min(n_eff_x, n_eff_y)
            eff_n_out[pixel_idx] = n_eff

            if np.isfinite(r[j]) and (n_eff - K - 2 > 0):
                df_eff = n_eff - K - 2
                t_stat = r[j] * np.sqrt(df_eff / (1.0 - r[j] * r[j]))
                p_out[pixel_idx] = 2.0 * stats.t.sf(np.abs(t_stat), df_eff)
            else:
                p_out[pixel_idx] = np.nan

        r_out[sub] = r.astype(np.float32)

    return r_out, p_out, eff_n_out


def residualize_batch(x, y, z):
    """
    x: (T, N)
    y: (T, N)
    z: (T, K, N)
    return rx, ry with shape (N, T)
    """
    x0 = x.T.astype(np.float64)
    y0 = y.T.astype(np.float64)
    z0 = np.transpose(z, (2, 0, 1)).astype(np.float64)

    ones = np.ones((z0.shape[0], z0.shape[1], 1), dtype=np.float64)
    A = np.concatenate([ones, z0], axis=2)

    AtA = np.einsum("ntp,ntq->npq", A, A)
    AtA_inv = np.linalg.pinv(AtA)

    Atx = np.einsum("ntp,nt->np", A, x0)
    Aty = np.einsum("ntp,nt->np", A, y0)

    beta_x = np.einsum("npq,nq->np", AtA_inv, Atx)
    beta_y = np.einsum("npq,nq->np", AtA_inv, Aty)

    rx = x0 - np.einsum("ntp,np->nt", A, beta_x)
    ry = y0 - np.einsum("ntp,np->nt", A, beta_y)

    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean(axis=1, keepdims=True)

    return rx, ry


# =========================================================
# 8. Diagnostic and statistical functions
# =========================================================
def check_multicollinearity(Z, var_names):
    Z = np.asarray(Z, dtype=np.float64)
    Z = Z[np.all(np.isfinite(Z), axis=1)]

    if Z.shape[0] == 0:
        return pd.DataFrame({"variable": var_names, "VIF": [np.nan] * len(var_names)})

    means = np.mean(Z, axis=0)
    stds = np.std(Z, axis=0, ddof=0)
    stds[stds == 0] = 1.0
    Z_scaled = (Z - means) / stds

    vif = []
    for i in range(Z_scaled.shape[1]):
        X_i = Z_scaled[:, i]
        X_others = np.delete(Z_scaled, i, axis=1)

        try:
            if X_others.shape[1] == 0:
                vif_i = 1.0
            else:
                beta = np.linalg.lstsq(X_others, X_i, rcond=None)[0]
                y_pred = X_others @ beta
                if np.std(X_i) == 0 or np.std(y_pred) == 0:
                    r2 = 0.0
                else:
                    r2 = np.corrcoef(X_i, y_pred)[0, 1] ** 2
                r2 = min(max(r2, 0.0), 0.999999)
                vif_i = 1.0 / (1.0 - r2)
        except Exception:
            vif_i = np.nan

        vif.append(vif_i)

    return pd.DataFrame({"variable": var_names, "VIF": vif})


def spatial_autocorr_test(r_map, max_points=2000):
    valid = np.isfinite(r_map)
    if np.sum(valid) < 100:
        return None, None

    coords = np.array(np.where(valid)).T
    values = r_map[valid].astype(np.float64)

    if len(values) > max_points:
        np.random.seed(42)
        idx = np.random.choice(len(values), size=max_points, replace=False)
        coords = coords[idx]
        values = values[idx]

    if len(values) < 50:
        return None, None

    dist_matrix = cdist(coords, coords)
    w_matrix = (dist_matrix < 1.5) & (dist_matrix > 0)

    if np.sum(w_matrix) == 0:
        return None, None

    n = len(values)
    z = values - np.mean(values)
    numerator = np.sum(w_matrix * np.outer(z, z))
    denominator = np.sum(z ** 2)
    s0 = np.sum(w_matrix)

    if denominator == 0 or s0 == 0:
        return None, None

    moran_i = (n / s0) * (numerator / denominator)

    try:
        ei = -1.0 / (n - 1)
        vari = (
            n ** 2 * np.sum((w_matrix + w_matrix.T) ** 2)
            - n * np.sum((np.sum(w_matrix, axis=1) + np.sum(w_matrix, axis=0)) ** 2)
            + s0 ** 2
        ) / (s0 ** 2 * (n ** 2 - 1))

        if vari <= 0 or not np.isfinite(vari):
            return moran_i, np.nan

        z_score = (moran_i - ei) / np.sqrt(vari)
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    except Exception:
        p_value = np.nan

    return moran_i, p_value


def summarize_pixel_distribution(r_vals, p_vals, n_eff_vals, label, alpha=0.05):
    r_vals = np.asarray(r_vals, dtype=np.float64)
    p_vals = np.asarray(p_vals, dtype=np.float64)
    n_eff_vals = np.asarray(n_eff_vals, dtype=np.float64)

    valid = np.isfinite(r_vals) & np.isfinite(p_vals)
    r = r_vals[valid]
    p = p_vals[valid]

    if len(r) == 0:
        return {
            "region": label,
            "valid_pixels": 0,
            "sig_pixels": 0,
            "sig_ratio": np.nan,
            "sig_pos_pixels": 0,
            "sig_pos_ratio": np.nan,
            "sig_neg_pixels": 0,
            "sig_neg_ratio": np.nan,
            "nonsig_pixels": 0,
            "nonsig_ratio": np.nan,
            "mean_r": np.nan,
            "median_r": np.nan,
            "std_r": np.nan,
            "min_r": np.nan,
            "max_r": np.nan,
            "p25_r": np.nan,
            "p75_r": np.nan,
            "mean_abs_r": np.nan,
            "mean_eff_n": np.nan,
            "median_eff_n": np.nan,
            "min_eff_n": np.nan,
            "max_eff_n": np.nan
        }

    sig = p < alpha
    sig_pos = sig & (r > 0)
    sig_neg = sig & (r < 0)
    eff = n_eff_vals[np.isfinite(n_eff_vals)]

    return {
        "region": label,
        "valid_pixels": int(len(r)),
        "sig_pixels": int(np.sum(sig)),
        "sig_ratio": float(np.mean(sig)),
        "sig_pos_pixels": int(np.sum(sig_pos)),
        "sig_pos_ratio": float(np.mean(sig_pos)),
        "sig_neg_pixels": int(np.sum(sig_neg)),
        "sig_neg_ratio": float(np.mean(sig_neg)),
        "nonsig_pixels": int(np.sum(~sig)),
        "nonsig_ratio": float(np.mean(~sig)),
        "mean_r": float(np.mean(r)),
        "median_r": float(np.median(r)),
        "std_r": float(np.std(r, ddof=0)),
        "min_r": float(np.min(r)),
        "max_r": float(np.max(r)),
        "p25_r": float(np.percentile(r, 25)),
        "p75_r": float(np.percentile(r, 75)),
        "mean_abs_r": float(np.mean(np.abs(r))),
        "mean_eff_n": float(np.mean(eff)) if len(eff) > 0 else np.nan,
        "median_eff_n": float(np.median(eff)) if len(eff) > 0 else np.nan,
        "min_eff_n": float(np.min(eff)) if len(eff) > 0 else np.nan,
        "max_eff_n": float(np.max(eff)) if len(eff) > 0 else np.nan
    }


def fisher_ci_from_r(r, n):
    if (not np.isfinite(r)) or (n is None) or (n <= 3) or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    z_low = z - 1.96 * se
    z_high = z + 1.96 * se
    return float(np.tanh(z_low)), float(np.tanh(z_high))


def make_integrated_conclusion(overall_stats, stratified_df):
    lines = [
        "Integrated conclusion for pixel-wise partial correlation between snow depth and surface soil water content (snow-month frequency >= 0.30; ssr covariate removed)",
        "=" * 50,
        "",
        "1. Overall result for the full domain",
        f"The full domain contains {overall_stats['valid_pixels']} valid pixels, including {overall_stats['sig_pixels']} significantly correlated pixels ({overall_stats['sig_ratio'] * 100:.2f}%).",
        f"Significantly positive pixels: {overall_stats['sig_pos_pixels']} ({overall_stats['sig_pos_ratio'] * 100:.2f}%); significantly negative pixels: {overall_stats['sig_neg_pixels']} ({overall_stats['sig_neg_ratio'] * 100:.2f}%).",
        f"The mean partial correlation coefficient for the full domain is {overall_stats['mean_r']:.4f}, the median is {overall_stats['median_r']:.4f}, and the IQR is [{overall_stats['p25_r']:.4f}, {overall_stats['p75_r']:.4f}].",
        ""
    ]

    if stratified_df.shape[0] > 0:
        lines.append("2. Results for different snow zones")
        for _, row in stratified_df.iterrows():
            region_name = row.get("region_name", row.get("region", "Unknown zone"))
            valid_pixels = int(row.get("valid_pixels", 0))
            if valid_pixels == 0:
                lines.append(f"{region_name}: no valid pixels.")
            else:
                lines.append(
                    f"{region_name}: {valid_pixels} valid pixels, "
                    f"significant-pixel ratio = {row['sig_ratio'] * 100:.2f}%, "
                    f"mean r = {row['mean_r']:.4f}, median r = {row['median_r']:.4f}."
                )

    return "\n".join(lines)


# =========================================================
# 9. Raster outputs and standard map plots
# =========================================================
def save_raster(out_tif, arr, profile, nodata=NODATA_OUT, dtype="float32"):
    profile_out = profile.copy()
    profile_out.update(dtype=dtype, count=1, nodata=nodata, compress=COMPRESS)

    if np.issubdtype(np.dtype(dtype), np.integer):
        out_arr = np.where(np.isfinite(arr), arr, nodata).astype(dtype)
    else:
        out_arr = np.where(np.isfinite(arr), arr.astype(np.float32), nodata).astype(dtype)

    with rasterio.open(out_tif, "w", **profile_out) as dst:
        dst.write(out_arr, 1)

    ensure_file_saved(out_tif, "TIF")


def _style_colorbar(cb, label_text):
    cb.set_label(label_text, fontsize=MAP_CBAR_LABEL_SIZE)
    cb.ax.tick_params(labelsize=MAP_CBAR_TICK_SIZE, width=1.0)


def _style_pixel_map_axes(ax, title):
    ax.set_title(title, fontsize=MAP_TITLE_SIZE)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(MAP_SPINE_WIDTH)


def save_corr_png(out_png, arr, title):
    fig, ax = plt.subplots(figsize=MAP_FIGSIZE, facecolor="white")
    masked = np.ma.masked_invalid(arr)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="white")

    im = ax.imshow(masked, cmap=cmap, vmin=R_VMIN, vmax=R_VMAX)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    _style_colorbar(cbar, "Partial correlation coefficient (r)")
    _style_pixel_map_axes(ax, title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


def save_corr_png_with_sig_points(out_png, r_arr, p_arr, title):
    fig, ax = plt.subplots(figsize=MAP_FIGSIZE, facecolor="white")
    masked = np.ma.masked_invalid(r_arr)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="white")

    im = ax.imshow(masked, cmap=cmap, vmin=R_VMIN, vmax=R_VMAX)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    _style_colorbar(cbar, "Partial correlation coefficient (r)")

    sig = np.isfinite(p_arr) & (p_arr < ALPHA)
    rr, cc = np.where(sig)
    if len(rr) > 0:
        keep = (rr % SIG_POINT_STRIDE == 0) & (cc % SIG_POINT_STRIDE == 0)
        rr2 = rr[keep]
        cc2 = cc[keep]
        ax.scatter(cc2, rr2, s=SIG_POINT_SIZE, c="k", alpha=0.7, linewidths=0)

    _style_pixel_map_axes(ax, title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


def save_p_png(out_png, p_arr, title):
    fig, ax = plt.subplots(figsize=MAP_FIGSIZE, facecolor="white")
    masked = np.ma.masked_invalid(p_arr)
    cmap = plt.cm.viridis_r.copy()
    cmap.set_bad(color="white")

    im = ax.imshow(masked, cmap=cmap, vmin=P_VMIN, vmax=P_VMAX)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    _style_colorbar(cbar, "p-value")
    _style_pixel_map_axes(ax, title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


def save_n_png(out_png, n_arr, title):
    fig, ax = plt.subplots(figsize=MAP_FIGSIZE, facecolor="white")
    masked = np.ma.masked_invalid(n_arr)
    cmap = plt.cm.plasma.copy()
    cmap.set_bad(color="white")

    if np.any(np.isfinite(n_arr)):
        vmin = np.nanpercentile(n_arr, 5)
        vmax = np.nanpercentile(n_arr, 95)
    else:
        vmin, vmax = 0, 1

    if not np.isfinite(vmin):
        vmin = 0
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = vmin + 1

    im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    _style_colorbar(cbar, "Effective sample size (n_eff)")
    _style_pixel_map_axes(ax, title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


def save_rpn_triptych(out_png, r_arr, p_arr, n_arr, main_title):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.8), facecolor="white")

    masked_r = np.ma.masked_invalid(r_arr)
    cmap_r = plt.cm.RdBu_r.copy()
    cmap_r.set_bad(color="white")
    im1 = axes[0].imshow(masked_r, cmap=cmap_r, vmin=R_VMIN, vmax=R_VMAX)
    cbar1 = fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.03)
    _style_colorbar(cbar1, "r")
    _style_pixel_map_axes(axes[0], "Partial correlation (r)")

    masked_p = np.ma.masked_invalid(p_arr)
    cmap_p = plt.cm.viridis_r.copy()
    cmap_p.set_bad(color="white")
    im2 = axes[1].imshow(masked_p, cmap=cmap_p, vmin=P_VMIN, vmax=P_VMAX)
    cbar2 = fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.03)
    _style_colorbar(cbar2, "p")
    _style_pixel_map_axes(axes[1], "P-value")

    masked_n = np.ma.masked_invalid(n_arr)
    cmap_n = plt.cm.plasma.copy()
    cmap_n.set_bad(color="white")

    if np.any(np.isfinite(n_arr)):
        vmin_n = np.nanpercentile(n_arr, 5)
        vmax_n = np.nanpercentile(n_arr, 95)
    else:
        vmin_n, vmax_n = 0, 1

    if not np.isfinite(vmin_n):
        vmin_n = 0
    if not np.isfinite(vmax_n) or vmax_n <= vmin_n:
        vmax_n = vmin_n + 1

    im3 = axes[2].imshow(masked_n, cmap=cmap_n, vmin=vmin_n, vmax=vmax_n)
    cbar3 = fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.03)
    _style_colorbar(cbar3, "n_eff")
    _style_pixel_map_axes(axes[2], "Effective sample size")

    fig.suptitle(main_title, fontsize=18)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


# =========================================================
# 10. Final publication-style large-scale correlation map
# =========================================================
def wrap_longitude_if_needed(arr, transform):
    """
    Convert a 0–360° raster arrangement to -180–180° when necessary.
    This preserves raster values and only reorders columns for display.
    """
    height, width = arr.shape
    right = transform.c + transform.a * width

    if right > 180:
        split_col = width // 2
        arr_fixed = np.hstack([arr[:, split_col:], arr[:, :split_col]])

        shift = transform.a * split_col
        new_transform = Affine(
            transform.a,
            transform.b,
            transform.c - shift,
            transform.d,
            transform.e,
            transform.f
        )
        return arr_fixed, new_transform

    return arr.copy(), transform


def save_lonfixed_raster(arr, profile, transform, out_tif):
    out_profile = profile.copy()
    out_profile.update(
        transform=transform,
        dtype="float32",
        nodata=NODATA_OUT,
        count=1,
        compress="LZW"
    )

    save_arr = np.where(
        np.isfinite(arr),
        arr.astype(np.float32),
        NODATA_OUT
    ).astype(np.float32)

    with rasterio.open(out_tif, "w", **out_profile) as dst:
        dst.write(save_arr, 1)

    ensure_file_saved(out_tif, "Longitude-adjusted GeoTIFF")


def save_publication_world_map(r_arr, profile, transform, crs, out_tif, out_png):
    """
    Draw the final publication-style large-scale pixel map
    with the required world background and horizontal correlation colorbar.
    """
    ensure_file_exists(WORLD_SHP, "World boundary shapefile")

    world = gpd.read_file(WORLD_SHP)

    if world.crs is None:
        raise ValueError(f"The world shapefile has no CRS information: {WORLD_SHP}")
    if crs is None:
        raise ValueError("The raster CRS is missing.")

    if world.crs != crs:
        world = world.to_crs(crs)

    height, width = r_arr.shape

    left = transform.c
    right = transform.c + transform.a * width
    top = transform.f
    bottom = transform.f + transform.e * height

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    ax = fig.add_axes([0.055, 0.12, 0.89, 0.78])
    ax.set_facecolor("white")

    # Global land background
    world.plot(
        ax=ax,
        facecolor=PUB_BG_RGB,
        edgecolor=PUB_EDGE_RGB,
        linewidth=0.5,
        zorder=1
    )

    # Partial-correlation raster
    masked_r = np.ma.masked_invalid(r_arr)
    image = ax.imshow(
        masked_r,
        extent=[left, right, bottom, top],
        origin="upper",
        cmap=CORR_CMAP,
        vmin=R_VMIN,
        vmax=R_VMAX,
        interpolation="none",
        zorder=2
    )

    # Country boundaries
    world.boundary.plot(
        ax=ax,
        color=PUB_EDGE_RGB,
        linewidth=0.5,
        zorder=3
    )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)

    # Remove background grid lines and all coordinates/ticks
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False
    )

    # Keep a regular black outer frame
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")

    ax.set_title("")

    # Horizontal correlation colorbar
    cax = ax.inset_axes([0.39, 0.045, 0.22, 0.028])
    colorbar = plt.colorbar(
        image,
        cax=cax,
        orientation="horizontal",
        extend="both"
    )
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])

    colorbar.ax.set_title(
        "Partial correlation coefficient (r)",
        fontsize=12,
        pad=6,
        fontproperties=ARIAL_FONT
    )
    colorbar.ax.tick_params(
        labelsize=11,
        direction="out",
        length=3,
        width=0.8,
        pad=2
    )
    for tick_label in colorbar.ax.get_xticklabels():
        tick_label.set_fontproperties(ARIAL_FONT)

    # Save TIFF
    try:
        fig.savefig(
            out_tif,
            dpi=FIG_DPI,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
            pil_kwargs={"compression": "tiff_lzw"}
        )
    except TypeError:
        fig.savefig(
            out_tif,
            dpi=FIG_DPI,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor()
        )

    # Save PNG
    fig.savefig(
        out_png,
        dpi=FIG_DPI,
        format="png",
        bbox_inches="tight",
        facecolor=fig.get_facecolor()
    )

    plt.close(fig)

    ensure_file_saved(out_tif, "Publication map TIFF")
    ensure_file_saved(out_png, "Publication map PNG")


# =========================================================
# 11. Stratified analysis within the snow-domain pixels
# =========================================================
def stratified_analysis_by_snow_freq(snow_freq_map, r_map, p_map, n_eff_map, width, height, profile):
    layers = {
        "permanent": {
            "name": "Permanent snow zone",
            "mask": snow_freq_map >= 0.8,
            "min_freq": 0.8,
            "max_freq": 1.0
        },
        "seasonal": {
            "name": "Seasonal snow zone",
            "mask": (snow_freq_map >= SNOWC_MIN_FRACTION) & (snow_freq_map < 0.8),
            "min_freq": SNOWC_MIN_FRACTION,
            "max_freq": 0.8
        }
    }

    results = {}

    for key, layer in layers.items():
        mask = layer["mask"]

        r_layer = r_map[mask]
        p_layer = p_map[mask]
        n_eff_layer = n_eff_map[mask]

        stats_dict = summarize_pixel_distribution(r_layer, p_layer, n_eff_layer, layer["name"], alpha=ALPHA)
        stats_dict["region"] = key
        stats_dict["region_name"] = layer["name"]
        stats_dict["min_freq"] = layer["min_freq"]
        stats_dict["max_freq"] = layer["max_freq"]

        results[key] = stats_dict

        r_layer_map = np.full((height, width), np.nan, dtype=np.float32)
        p_layer_map = np.full((height, width), np.nan, dtype=np.float32)
        n_layer_map = np.full((height, width), np.nan, dtype=np.float32)
        sig_layer_map = np.full((height, width), 255, dtype=np.uint8)

        r_layer_map[mask] = r_layer
        p_layer_map[mask] = p_layer
        n_layer_map[mask] = n_eff_layer

        sig_local = np.isfinite(p_layer) & (p_layer < ALPHA)
        sig_layer_map[mask] = sig_local.astype(np.uint8)

        out_r_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_r_{key}.tif")
        out_p_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_p_{key}.tif")
        out_n_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_effective_n_{key}.tif")
        out_sig_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_sigmask_{key}.tif")

        save_raster(out_r_tif, r_layer_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_p_tif, p_layer_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_n_tif, n_layer_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_sig_tif, sig_layer_map, profile, nodata=255, dtype="uint8")

        out_r_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_r_{key}.png")
        out_r_sig_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_r_{key}_sig_points.png")
        out_p_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_p_{key}.png")
        out_n_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_effective_n_{key}.png")
        out_rpn_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_rpn_triptych_{key}.png")

        save_corr_png(out_r_png, r_layer_map, f"Partial correlation: snow depth vs surface soil water content ({layer['name']})")
        save_corr_png_with_sig_points(
            out_r_sig_png,
            r_layer_map,
            p_layer_map,
            f"Partial correlation: snow depth vs surface soil water content ({layer['name']}, p<{ALPHA})"
        )
        save_p_png(out_p_png, p_layer_map, f"P-value: snow depth vs surface soil water content ({layer['name']})")
        save_n_png(out_n_png, n_layer_map, f"Effective sample size ({layer['name']})")
        save_rpn_triptych(
            out_rpn_png,
            r_layer_map,
            p_layer_map,
            n_layer_map,
            f"Partial correlation result: {layer['name']}"
        )

    return results


# =========================================================
# 12. Main raster-calculation worker
# =========================================================
def process_one_block(args):
    row_start, row_end, width = args
    h = row_end - row_start
    needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))

    snowc_raw = load_var_cube("snowc", row_start, row_end, width)
    snow_fraction = np.nanmean(snowc_raw > SNOWC_THRESHOLD, axis=0)
    snow_mask_block = np.isfinite(snow_fraction) & (snow_fraction >= SNOWC_MIN_FRACTION)

    if not np.any(snow_mask_block):
        return (
            row_start, row_end,
            np.full((h, width), np.nan, dtype=np.float32),
            np.full((h, width), np.nan, dtype=np.float32),
            np.full((h, width), np.nan, dtype=np.float32),
            np.full((h, width), np.nan, dtype=np.float32),
            snow_mask_block.astype(np.uint8)
        )

    data_proc = {}
    for var in needed_vars:
        cube_raw = snowc_raw if var == "snowc" else load_var_cube(var, row_start, row_end, width)
        cube_proc = remove_monthly_climatology_3d(cube_raw) if REMOVE_MONTHLY_CLIM else cube_raw
        data_proc[var] = cube_proc.reshape(len(TIME_STEPS), -1).astype(np.float32)

        if var != "snowc":
            del cube_raw
        del cube_proc

    del snowc_raw

    x = data_proc[X_VAR]
    y = data_proc[TARGET_VAR]
    z = np.stack([data_proc[c] for c in COVARS], axis=1)

    r_flat, p_flat, eff_n_flat = partial_corr_and_p_batch_with_eff_df(
        x=x, y=y, z=z, pixel_chunk=PIXEL_CHUNK
    )

    mask_flat = snow_mask_block.reshape(-1)
    r_flat[~mask_flat] = np.nan
    p_flat[~mask_flat] = np.nan
    eff_n_flat[~mask_flat] = np.nan

    return (
        row_start, row_end,
        r_flat.reshape(h, width).astype(np.float32),
        p_flat.reshape(h, width).astype(np.float32),
        eff_n_flat.reshape(h, width).astype(np.float32),
        snow_fraction.astype(np.float32),
        snow_mask_block.astype(np.uint8)
    )


# =========================================================
# 13. Scatter plot statistics and outputs
# =========================================================
def init_scatter_acc(title, out_png, out_csv):
    return {
        "title": title,
        "out_png": out_png,
        "out_csv": out_csv,
        "n": 0,
        "sum_x": 0.0,
        "sum_y": 0.0,
        "sum_x2": 0.0,
        "sum_y2": 0.0,
        "sum_xy": 0.0,
        "sample_x_chunks": [],
        "sample_y_chunks": []
    }


def remove_soil_water_content_outliers_iqr(x, y, k=1.5):
    """
    Remove outliers only based on surface soil water content (y).
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if y.size < 8:
        return x, y

    q1, q3 = np.percentile(y, [25, 75])
    iqr = q3 - q1

    if (not np.isfinite(iqr)) or (iqr <= 0):
        return x, y

    lower = q1 - k * iqr
    upper = q3 + k * iqr
    keep = (y >= lower) & (y <= upper)

    return x[keep], y[keep]


def update_scatter_acc(acc, x_vals, y_vals, rng):
    if x_vals.size == 0:
        return

    valid = np.isfinite(x_vals) & np.isfinite(y_vals)
    x = x_vals[valid].astype(np.float64, copy=False)
    y = y_vals[valid].astype(np.float64, copy=False)

    if x.size == 0:
        return

    x, y = remove_soil_water_content_outliers_iqr(x, y, k=1.5)
    if x.size == 0:
        return

    acc["n"] += int(x.size)
    acc["sum_x"] += float(np.sum(x))
    acc["sum_y"] += float(np.sum(y))
    acc["sum_x2"] += float(np.sum(x * x))
    acc["sum_y2"] += float(np.sum(y * y))
    acc["sum_xy"] += float(np.sum(x * y))

    keep = min(x.size, PLOT_POINTS_PER_BLOCK_PER_GROUP)
    if x.size > keep:
        idx = rng.choice(x.size, size=keep, replace=False)
        xs = x[idx].astype(np.float32)
        ys = y[idx].astype(np.float32)
    else:
        xs = x.astype(np.float32)
        ys = y.astype(np.float32)

    acc["sample_x_chunks"].append(xs)
    acc["sample_y_chunks"].append(ys)


def finalize_scatter_acc(acc):
    n = acc["n"]

    if n < 3:
        r = np.nan
        p = np.nan
        slope = np.nan
        intercept = np.nan
        ci_low = np.nan
        ci_high = np.nan
    else:
        mean_x = acc["sum_x"] / n
        mean_y = acc["sum_y"] / n
        sxx = acc["sum_x2"] - acc["sum_x"] * acc["sum_x"] / n
        syy = acc["sum_y2"] - acc["sum_y"] * acc["sum_y"] / n
        sxy = acc["sum_xy"] - acc["sum_x"] * acc["sum_y"] / n

        if sxx <= 0 or syy <= 0:
            r = np.nan
            p = np.nan
            slope = np.nan
            intercept = np.nan
            ci_low = np.nan
            ci_high = np.nan
        else:
            r = float(np.clip(sxy / np.sqrt(sxx * syy), -0.999999, 0.999999))
            slope = float(sxy / sxx)
            intercept = float(mean_y - slope * mean_x)

            df = n - 2
            if df > 0 and abs(r) < 1:
                t_stat = r * np.sqrt(df / (1.0 - r * r))
                p = float(2.0 * stats.t.sf(abs(t_stat), df))
            else:
                p = np.nan

            ci_low, ci_high = fisher_ci_from_r(r, n)

    if len(acc["sample_x_chunks"]) > 0:
        sx = np.concatenate(acc["sample_x_chunks"])
        sy = np.concatenate(acc["sample_y_chunks"])
        if sx.size > MAX_PLOT_POINTS_PER_GROUP:
            rng = np.random.default_rng(SCATTER_RANDOM_STATE)
            idx = rng.choice(sx.size, size=MAX_PLOT_POINTS_PER_GROUP, replace=False)
            sx = sx[idx]
            sy = sy[idx]
    else:
        sx = np.array([], dtype=np.float32)
        sy = np.array([], dtype=np.float32)

    return {
        "title": acc["title"],
        "out_png": acc["out_png"],
        "out_csv": acc["out_csv"],
        "n": int(n),
        "r": r,
        "p": p,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "slope": slope,
        "intercept": intercept,
        "sample_x": sx,
        "sample_y": sy,
    }


def format_p_value(p):
    if not np.isfinite(p):
        return "NaN"
    if p < 0.001:
        return "< 0.001"
    return f"{p:.3f}"


def save_partial_scatter_png(out_png, sx, sy, r, p, n, slope, intercept, title):
    sx = np.asarray(sx, dtype=np.float64).ravel()
    sy = np.asarray(sy, dtype=np.float64).ravel()

    valid = np.isfinite(sx) & np.isfinite(sy)
    sx = sx[valid]
    sy = sy[valid]

    sx, sy = remove_soil_water_content_outliers_iqr(sx, sy, k=1.5)

    fig, ax = plt.subplots(figsize=SCATTER_FIGSIZE, facecolor="white")

    if sx.size > 0 and sy.size > 0:
        ax.scatter(
            sx, sy,
            s=18,
            facecolor=SCATTER_FACE,
            edgecolor=SCATTER_EDGE,
            linewidth=0.35,
            alpha=0.60
        )

        if np.isfinite(slope) and np.isfinite(intercept):
            xs = np.linspace(np.nanmin(sx), np.nanmax(sx), 200)
            ys = slope * xs + intercept
            ax.plot(xs, ys, color=REG_LINE, linewidth=2.0)

    txt = (
        f"r = {r:.3f}\n"
        f"p = {format_p_value(p)}\n"
        f"n = {int(n):,}"
    )

    ax.text(
        0.04, 0.06, txt,
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=SCATTER_TEXT_SIZE,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.92)
    )

    ax.set_xlabel(f"Residualized {X_DISPLAY_NAME}", fontsize=SCATTER_LABEL_SIZE)
    ax.set_ylabel(f"Residualized {TARGET_DISPLAY_NAME}", fontsize=SCATTER_LABEL_SIZE)
    ax.set_title(title, fontsize=SCATTER_TITLE_SIZE)

    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(SCATTER_SPINE_WIDTH)

    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=5,
        width=1.0,
        labelsize=SCATTER_TICK_SIZE
    )

    ax.grid(False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    ensure_file_saved(out_png, "PNG")


def generate_scatter_outputs(template, r_map, p_map, snow_freq_map):
    valid_pixels = np.isfinite(r_map) & np.isfinite(p_map)
    sig_pos = valid_pixels & (p_map < ALPHA) & (r_map > 0)
    sig_neg = valid_pixels & (p_map < ALPHA) & (r_map < 0)
    nonsig = valid_pixels & ~(sig_pos | sig_neg)

    permanent = valid_pixels & (snow_freq_map >= 0.8)
    seasonal = valid_pixels & (snow_freq_map >= SNOWC_MIN_FRACTION) & (snow_freq_map < 0.8)

    groups = {
        "overall": {
            "mask": valid_pixels,
            "title": "Partial regression scatter: snow depth vs surface soil water content (snow-month frequency >= 0.30)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_overall_all_valid_pixels.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_overall_all_valid_pixels.csv"),
        },
        "significant_positive": {
            "mask": sig_pos,
            "title": "Partial regression scatter: snow depth vs surface soil water content (significant positive pixels)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_significant_positive_pixels.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_significant_positive_pixels.csv"),
        },
        "significant_negative": {
            "mask": sig_neg,
            "title": "Partial regression scatter: snow depth vs surface soil water content (significant negative pixels)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_significant_negative_pixels.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_significant_negative_pixels.csv"),
        },
        "non_significant": {
            "mask": nonsig,
            "title": "Partial regression scatter: snow depth vs surface soil water content (non-significant pixels)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_non_significant_pixels.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_non_significant_pixels.csv"),
        },
        "permanent_zone": {
            "mask": permanent,
            "title": "Partial regression scatter: snow depth vs surface soil water content (permanent snow zone)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_permanent_snow_zone.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_permanent_snow_zone.csv"),
        },
        "seasonal_zone": {
            "mask": seasonal,
            "title": "Partial regression scatter: snow depth vs surface soil water content (seasonal snow zone, snow-month frequency >= 0.30)",
            "png": os.path.join(OUT_SCATTER_FIG, "scatter_seasonal_snow_zone_snowfreq_ge_0p30.png"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_seasonal_snow_zone_snowfreq_ge_0p30.csv"),
        },
    }

    accs = {k: init_scatter_acc(v["title"], v["png"], v["csv"]) for k, v in groups.items()}
    rng = np.random.default_rng(SCATTER_RANDOM_STATE)

    height = template["height"]
    width = template["width"]
    ranges = block_ranges(height, BLOCK_ROWS)

    log_msg("Starting generation of partial regression scatter plots based on the original monthly raster series...")

    for row_start, row_end in tqdm(ranges, desc="scatter", ncols=100):
        block_valid = valid_pixels[row_start:row_end, :]
        if not np.any(block_valid):
            continue

        needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))
        cubes = {}
        for var in needed_vars:
            cube_raw = load_var_cube(var, row_start, row_end, width)
            cubes[var] = remove_monthly_climatology_3d(cube_raw) if REMOVE_MONTHLY_CLIM else cube_raw

        x_all = cubes[X_VAR].reshape(len(TIME_STEPS), -1)
        y_all = cubes[TARGET_VAR].reshape(len(TIME_STEPS), -1)
        z_all = np.stack([cubes[c].reshape(len(TIME_STEPS), -1) for c in COVARS], axis=1)

        block_valid_flat = block_valid.reshape(-1)
        idx0 = np.where(block_valid_flat)[0]
        if idx0.size == 0:
            del cubes, x_all, y_all, z_all
            continue

        strict_ok = (
            np.all(np.isfinite(x_all[:, idx0]), axis=0) &
            np.all(np.isfinite(y_all[:, idx0]), axis=0) &
            np.all(np.isfinite(z_all[:, :, idx0]), axis=(0, 1))
        )
        idx = idx0[strict_ok]
        if idx.size == 0:
            del cubes, x_all, y_all, z_all
            continue

        x_sel = x_all[:, idx]
        y_sel = y_all[:, idx]
        z_sel = z_all[:, :, idx]

        rx, ry = residualize_batch(x_sel, y_sel, z_sel)

        for gname, ginfo in groups.items():
            gmask_flat = ginfo["mask"][row_start:row_end, :].reshape(-1)
            membership = gmask_flat[idx]
            if not np.any(membership):
                continue

            xg = rx[membership].reshape(-1)
            yg = ry[membership].reshape(-1)
            update_scatter_acc(accs[gname], xg, yg, rng)

        del cubes, x_all, y_all, z_all, x_sel, y_sel, z_sel, rx, ry

    rows = []
    for gname, acc in accs.items():
        res = finalize_scatter_acc(acc)

        save_partial_scatter_png(
            res["out_png"],
            res["sample_x"],
            res["sample_y"],
            res["r"],
            res["p"],
            res["n"],
            res["slope"],
            res["intercept"],
            res["title"]
        )

        row = {
            "group": gname,
            "title": res["title"],
            "figure_file": os.path.basename(res["out_png"]),
            "n": res["n"],
            "partial_r": res["r"],
            "p_val": res["p"],
            "ci_low": res["ci_low"],
            "ci_high": res["ci_high"],
            "slope": res["slope"],
            "intercept": res["intercept"],
        }
        rows.append(row)
        save_df_csv(pd.DataFrame([row]), res["out_csv"])

    df_all = pd.DataFrame(rows)
    save_df_csv(df_all, os.path.join(OUT_SCATTER_TABLE, "scatter_groups_summary.csv"))

    readme = [
        "Partial regression scatter plot notes",
        "=" * 40,
        "These scatter plots are reconstructed from the original monthly raster series.",
        "Workflow:",
        "1. Group pixels using the existing pixel-wise r_map / p_map / snow_freq_map within the snow-domain pixels;",
        "   The snow-domain criterion is snow-month frequency >= 0.30, with snow month defined as snowc > 0.",
        "2. Return to the original monthly series of all pixels within each group;",
        "3. Remove covariate effects from snow depth and surface soil water content separately to obtain residualized x/y;",
        "4. Pool residual samples from all pixels and all months within each group for plotting;",
        "5. r / p / n / CI shown in each plot are recalculated using all samples in that group;",
        "6. Outliers in surface soil water content (y) are removed using the IQR method;",
        f"7. Displayed points are sampled for visualization only (maximum {MAX_PLOT_POINTS_PER_GROUP}), whereas statistics use all valid samples.",
        "8. The covariate ssr has been removed in this version."
    ]
    save_txt("\n".join(readme), os.path.join(OUT_SCATTER_TABLE, "scatter_readme.txt"))

    log_msg("Scatter plots and corresponding statistical tables have been saved.")


# =========================================================
# 14. Main program
# =========================================================
def main():
    log_msg(f"========== {ANALYSIS_TITLE} (ssr covariate removed) ==========")
    log_msg(f"Time range: {YEARS[0]}-{YEARS[-1]}, {len(TIME_STEPS)} months in total")
    log_msg(f"Dependent variable: {TARGET_VAR} ({TARGET_DESC})")
    log_msg(f"Independent variable: {X_VAR}")
    log_msg(f"Covariates: {', '.join(COVARS)}")
    log_msg(f"Remove monthly climatology: {REMOVE_MONTHLY_CLIM}")
    log_msg(f"Output root directory: {OUT_ROOT}")
    log_msg(f"Snow-domain criterion: snow-month frequency >= {SNOWC_MIN_FRACTION:.2f}; snow month is defined as snowc > {SNOWC_THRESHOLD}")

    template = get_template_info()
    init_worker_template(template["crs"], template["transform"], template["width"], template["height"])
    check_inputs(template)

    width = template["width"]
    height = template["height"]
    profile = template["profile"]
    ranges = block_ranges(height, BLOCK_ROWS)

    log_msg(f"Raster size: {width} x {height}")
    log_msg(f"Total row blocks: {len(ranges)}")

    r_map = np.full((height, width), np.nan, dtype=np.float32)
    p_map = np.full((height, width), np.nan, dtype=np.float32)
    eff_n_map = np.full((height, width), np.nan, dtype=np.float32)
    snow_freq_map = np.full((height, width), np.nan, dtype=np.float32)
    snow_mask_map = np.zeros((height, width), dtype=np.uint8)

    tasks = [(r0, r1, width) for (r0, r1) in ranges]

    with ProcessPoolExecutor(
        max_workers=N_JOBS,
        initializer=init_worker_template,
        initargs=(template["crs"], template["transform"], template["width"], template["height"])
    ) as ex:
        futures = [ex.submit(process_one_block, task) for task in tasks]

        for fut in tqdm(as_completed(futures), total=len(futures), desc="swvl1", ncols=100):
            row_start, row_end, r_block, p_block, eff_n_block, snow_freq_block, mask_block = fut.result()
            r_map[row_start:row_end, :] = r_block
            p_map[row_start:row_end, :] = p_block
            eff_n_map[row_start:row_end, :] = eff_n_block
            snow_freq_map[row_start:row_end, :] = snow_freq_block
            snow_mask_map[row_start:row_end, :] = mask_block

    sig_mask = (np.isfinite(p_map) & (p_map < ALPHA)).astype(np.uint8)

    # 1) Save main full-domain results
    out_r_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_r.tif")
    out_p_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_p.tif")
    out_eff_n_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_effective_n.tif")
    out_sig_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_partial_corr_sigmask.tif")
    out_mask_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_snowfreq_ge_0p30_domain_mask.tif")
    out_snow_freq_tif = os.path.join(OUT_RASTER, f"{OUTPUT_PREFIX}_snowc_frequency.tif")

    save_raster(out_r_tif, r_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_p_tif, p_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_eff_n_tif, eff_n_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_sig_tif, sig_mask, profile, nodata=255, dtype="uint8")
    save_raster(out_mask_tif, snow_mask_map, profile, nodata=255, dtype="uint8")
    save_raster(out_snow_freq_tif, snow_freq_map, profile, nodata=NODATA_OUT, dtype="float32")

    out_r_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_r.png")
    out_r_sig_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_r_sig_points.png")
    out_p_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_p.png")
    out_n_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_effective_n.png")
    out_rpn_png = os.path.join(OUT_FIG, f"{OUTPUT_PREFIX}_partial_corr_rpn_triptych.png")

    save_corr_png(out_r_png, r_map, "Pixel-wise partial correlation: snow depth vs surface soil water content")
    save_corr_png_with_sig_points(
        out_r_sig_png,
        r_map,
        p_map,
        f"Pixel-wise partial correlation: snow depth vs surface soil water content (significant points, p<{ALPHA})"
    )
    save_p_png(out_p_png, p_map, "Pixel-wise p-value: snow depth vs surface soil water content")
    save_n_png(out_n_png, eff_n_map, "Effective sample size: snow depth vs surface soil water content")
    save_rpn_triptych(
        out_rpn_png,
        r_map,
        p_map,
        eff_n_map,
        "Overall partial correlation result: snow depth vs surface soil water content"
    )

    log_msg("Full-domain TIF and PNG outputs have been saved.")

    # 2) Full-domain summary table
    overall_mask = np.isfinite(r_map) & np.isfinite(p_map)
    overall_stats = summarize_pixel_distribution(
        r_map[overall_mask],
        p_map[overall_mask],
        eff_n_map[overall_mask],
        "Overall",
        alpha=ALPHA
    )
    overall_stats["region_name"] = "Overall"

    df_overall = pd.DataFrame([overall_stats])
    out_overall_csv = os.path.join(OUT_TABLE, f"{OUTPUT_PREFIX}_overall_snowfreq_ge_0p30_summary.csv")
    save_df_csv(df_overall, out_overall_csv)
    log_msg("Full-domain summary table has been saved.")

    # 3) Stratified figures + tables
    stratified_results = stratified_analysis_by_snow_freq(
        snow_freq_map, r_map, p_map, eff_n_map, width, height, profile
    )

    df_stratified = pd.DataFrame(list(stratified_results.values()))
    out_stratified_csv = os.path.join(OUT_TABLE, f"{OUTPUT_PREFIX}_stratified_snowfreq_ge_0p30_summary.csv")
    save_df_csv(df_stratified, out_stratified_csv)
    log_msg("Figures and statistical tables for the snow-zone classes have been saved.")

    # 4) Scatter plots
    generate_scatter_outputs(template, r_map, p_map, snow_freq_map)

    # 5) Integrated summary table + integrated conclusion
    df_integrated = pd.concat(
        [df_overall.assign(region="overall"), df_stratified],
        ignore_index=True,
        sort=False
    )
    out_integrated_csv = os.path.join(OUT_TABLE, f"{OUTPUT_PREFIX}_integrated_partial_corr_snowfreq_ge_0p30_summary.csv")
    save_df_csv(df_integrated, out_integrated_csv)

    out_conclusion_txt = os.path.join(OUT_TEXT, f"{OUTPUT_PREFIX}_integrated_partial_corr_snowfreq_ge_0p30_conclusion.txt")
    save_txt(make_integrated_conclusion(overall_stats, df_stratified), out_conclusion_txt)
    log_msg("Integrated summary table and integrated conclusion text have been saved.")

    # 6) Diagnostics
    try:
        log_msg("Starting diagnostic analysis (VIF / Moran)...")

        needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))
        mean_series = {}

        for var in needed_vars:
            vals_sum = np.zeros(len(TIME_STEPS), dtype=np.float64)
            vals_cnt = np.zeros(len(TIME_STEPS), dtype=np.int64)

            for i, (year, month) in enumerate(TIME_STEPS):
                arr = read_one_block(var, tif_path(var, year, month), 0, height, width)
                valid_mask = np.isfinite(arr) & (snow_mask_map == 1)
                vals_sum[i] = np.nansum(arr[valid_mask])
                vals_cnt[i] = np.sum(valid_mask)

            s = np.divide(vals_sum, vals_cnt, out=np.full(len(TIME_STEPS), np.nan), where=vals_cnt > 0)
            if REMOVE_MONTHLY_CLIM:
                s = remove_monthly_climatology_1d(s)

            mean_series[var] = s

        Z = np.column_stack([mean_series[c] for c in COVARS])
        vif_df = check_multicollinearity(Z, COVARS)
        out_vif_csv = os.path.join(OUT_TABLE, "vif_diagnosis.csv")
        save_df_csv(vif_df, out_vif_csv)

        moran_i, moran_p = spatial_autocorr_test(r_map)
        df_diag = pd.DataFrame([{
            "moran_i": moran_i if moran_i is not None else np.nan,
            "moran_p": moran_p if moran_p is not None else np.nan
        }])
        out_diag_csv = os.path.join(OUT_TABLE, "spatial_diagnosis.csv")
        save_df_csv(df_diag, out_diag_csv)

        log_msg("VIF and spatial autocorrelation diagnostic tables have been saved.")

    except Exception as e:
        err_txt = os.path.join(OUT_TEXT, "diagnostic_error.txt")
        save_txt(str(e), err_txt)
        log_msg(f"An error occurred during diagnostics, but the main figures and statistical tables have been generated. See error information at: {err_txt}")

    # 7) Final publication-style large-scale pixel map
    try:
        log_msg("Starting the final publication-style large-scale pixel map...")

        r_map_pub, transform_pub = wrap_longitude_if_needed(r_map, template["transform"])

        out_lonfix_tif = os.path.join(OUT_PUB_MAP, f"{OUTPUT_PREFIX}_partial_corr_r_lonfix.tif")
        out_pub_tif = os.path.join(OUT_PUB_MAP, f"{OUTPUT_PREFIX}_partial_corr_world_map_pub.tif")
        out_pub_png = os.path.join(OUT_PUB_MAP, f"{OUTPUT_PREFIX}_partial_corr_world_map_pub.png")

        save_lonfixed_raster(
            r_map_pub,
            profile,
            transform_pub,
            out_lonfix_tif
        )

        save_publication_world_map(
            r_arr=r_map_pub,
            profile=profile,
            transform=transform_pub,
            crs=template["crs"],
            out_tif=out_pub_tif,
            out_png=out_pub_png
        )

        log_msg("The final publication-style large-scale pixel map has been saved.")

    except Exception as e:
        pub_err_txt = os.path.join(OUT_TEXT, "publication_map_error.txt")
        save_txt(str(e), pub_err_txt)
        log_msg(f"An error occurred while generating the publication-style map. See: {pub_err_txt}")

    log_msg("=" * 60)
    log_msg(f"Result root directory: {OUT_ROOT}")
    log_msg(f"Raster outputs: {OUT_RASTER}")
    log_msg(f"Figures: {OUT_FIG}")
    log_msg(f"Tables: {OUT_TABLE}")
    log_msg(f"Scatter figures: {OUT_SCATTER_FIG}")
    log_msg(f"Scatter tables: {OUT_SCATTER_TABLE}")
    log_msg(f"Publication maps: {OUT_PUB_MAP}")
    log_msg("=" * 60)


if __name__ == "__main__":
    freeze_support()
    main()