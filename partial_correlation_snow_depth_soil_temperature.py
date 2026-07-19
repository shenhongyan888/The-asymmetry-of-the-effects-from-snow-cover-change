# -*- coding: utf-8 -*-
"""
Pixel-wise partial correlation between snow depth and surface soil temperature
================================================================================

This script combines:
1. The calculation workflow:
   - Monthly raster time series from 2000-01 to 2025-12
   - Monthly climatology removal
   - AR(1)-based effective sample size correction for significance testing
   - Pixel-wise partial correlation between snow depth and surface soil temperature
   - Full covariate control
   - Snow-cover-domain masking based on snow-month frequency >= 0.30
   - Automatic grid alignment to the SDE template

2. The plotting and reporting workflow:
   - GeoTIFF outputs for r, p, significance mask, effective sample size, snow frequency, and snow-domain mask
   - Quicklook PNG figures
   - Publication-style world map in TIF format
   - Stratified statistics by snow-frequency zone
   - Publication-style partial-regression scatter plots
   - Domain-mean partial-regression scatter plot
   - VIF and spatial autocorrelation diagnostics

Notes
-----
- All text in comments, logs, figure labels, and output tables is written in English.
- The scatter plots are displayed using sampled points for readability, but their statistics
  are calculated from all valid samples in each group.
- If the CRS or grid differs among input rasters, rasters are read through WarpedVRT and
  resampled to the SDE template grid.
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


# =============================================================================
# 1. Main configuration
# =============================================================================
OUT_ROOT = r"I:\ERA_TIFF\!partial_corr_stl1_final_outputs_v2_english_full"

OUT_RASTER = os.path.join(OUT_ROOT, "01_rasters")
OUT_FIG = os.path.join(OUT_ROOT, "02_quicklook_figures")
OUT_TABLE = os.path.join(OUT_ROOT, "03_tables")
OUT_TEXT = os.path.join(OUT_ROOT, "04_texts")
OUT_SCATTER_FIG = os.path.join(OUT_ROOT, "05_scatter_figures")
OUT_SCATTER_TABLE = os.path.join(OUT_ROOT, "06_scatter_tables")
OUT_PUB_MAP = os.path.join(OUT_ROOT, "07_publication_maps")
OUT_DOMAIN_SCATTER = os.path.join(OUT_ROOT, "08_domain_mean_scatter")

for folder in [
    OUT_ROOT, OUT_RASTER, OUT_FIG, OUT_TABLE, OUT_TEXT,
    OUT_SCATTER_FIG, OUT_SCATTER_TABLE, OUT_PUB_MAP, OUT_DOMAIN_SCATTER
]:
    os.makedirs(folder, exist_ok=True)

RUN_LOG = os.path.join(OUT_ROOT, "run_log.txt")
if os.path.exists(RUN_LOG):
    os.remove(RUN_LOG)

YEARS = list(range(2000, 2026))
MONTHS = list(range(1, 13))
TIME_STEPS = [(year, month) for year in YEARS for month in MONTHS]
DATES = pd.to_datetime([f"{year}-{month:02d}-01" for year, month in TIME_STEPS])
MONTH_IDS = np.array([month for _, month in TIME_STEPS], dtype=np.int16)

TARGET_VAR = "stl1"
TARGET_DESC = "Surface soil temperature"
# Use TARGET_VAR = "skt" only if you intentionally want skin temperature instead of soil temperature.
X_VAR = "sde"
X_DESC = "Snow depth equivalent"

COVARS = [
    "t2m",
    "ssr",
    "snowc",
    "sf",
    "smlt",
    "v10",
    "lai_hv"
]

NODATA_OUT = -9999.0
COMPRESS = "LZW"
ALPHA = 0.05

REMOVE_MONTHLY_CLIM = True
USE_AR1_PREWHITEN = False

SNOWC_MASK_RULE = "frequency_ge_0_30"
SNOWC_THRESHOLD = 0.0
SNOWC_MIN_FRACTION = 0.30
USE_EFFECTIVE_SAMPLE_SIZE = True

# Snow-frequency classes used only for stratified summaries and scatter plots.
PERMANENT_SNOW_MIN = 0.80
SEASONAL_SNOW_MIN = 0.30
EPHEMERAL_SNOW_MIN = 0.05

N_JOBS = 4
BLOCK_ROWS = 8
PIXEL_CHUNK = 2000

SIG_POINT_STRIDE = 8
SIG_POINT_SIZE = 1.2

R_VMIN, R_VMAX = -1.0, 1.0
P_VMIN, P_VMAX = 0.0, 0.05

CHECK_ALL_FILES = True
CHECK_ALIGNMENT = True


# =============================================================================
# 2. Publication plotting configuration
# =============================================================================
FIG_DPI = 600

# Confirmed world boundary path.
WORLD_SHP = "I:\\u4e16\u754c\u5730\u56fe\\u4e16\u754c\u5730\u56fe\\global_all_country.shp"

# Keep the original global font setting for figures that were not requested
# to change. The correlation legend is explicitly formatted in Arial below.
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False

BG_RGB = (243 / 255, 243 / 255, 243 / 255)
EDGE_RGB = (175 / 255, 175 / 255, 175 / 255)

# Requested correlation colors:
# negative correlation = RGB (32, 56, 136)
# zero correlation     = white
# positive correlation = RGB (225, 156, 102)
NEGATIVE_CORR_RGB = (32 / 255, 56 / 255, 136 / 255)
ZERO_CORR_RGB = (1.0, 1.0, 1.0)
POSITIVE_CORR_RGB = (225 / 255, 156 / 255, 102 / 255)

CORR_CMAP = LinearSegmentedColormap.from_list(
    "custom_partial_correlation",
    [NEGATIVE_CORR_RGB, ZERO_CORR_RGB, POSITIVE_CORR_RGB],
    N=256
)
CORR_CMAP.set_bad(color="white")

ARIAL_FONT = FontProperties(family="Arial")

SCATTER_FACE = (238 / 255, 156 / 255, 74 / 255)
SCATTER_EDGE = (175 / 255, 175 / 255, 175 / 255)
REG_LINE = (175 / 255, 175 / 255, 175 / 255)

PLOT_POINTS_PER_BLOCK_PER_GROUP = 1200
MAX_PLOT_POINTS_PER_GROUP = 30000
SCATTER_RANDOM_STATE = 42


def apply_arial_to_colorbar(colorbar):
    """Apply Arial to a correlation colorbar title, label, and tick labels."""
    colorbar.ax.title.set_fontproperties(ARIAL_FONT)
    colorbar.ax.xaxis.label.set_fontproperties(ARIAL_FONT)
    colorbar.ax.yaxis.label.set_fontproperties(ARIAL_FONT)

    for tick_label in colorbar.ax.get_xticklabels():
        tick_label.set_fontproperties(ARIAL_FONT)

    for tick_label in colorbar.ax.get_yticklabels():
        tick_label.set_fontproperties(ARIAL_FONT)


# =============================================================================
# 3. Input path configuration
# =============================================================================
ERA_BASE_DIR = r"I:\ERA_TIFF"

ERA_STYLE_VARS = {
    "stl1": "stl1",
    "skt": "skt",
    "sde": "sde",
    "es": "es",
    "tp": "tp",
    "lai_hv": "lai_hv",
    "lai_lv": "lai_lv",
    "sf": "sf",
    "smlt": "smlt",
    "tsn": "tsn",
    "snowc": "snowc",
    "src": "src",
}

SPECIAL_MONTHLY_VARS = {
    "asn": {
        "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\asn\asn",
        "prefix": "a25738c980fe5e74f46d84d68f745a99_asn_"
    },
    "slhf": {
        "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\slhf\slhf",
        "prefix": "9ce01cf43be7606676d84be71fe19678_slhf_"
    },
    "ssr": {
        "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\ssr\ssr",
        "prefix": "9ce01cf43be7606676d84be71fe19678_ssr_"
    },
    "str": {
        "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\str\str",
        "prefix": "9ce01cf43be7606676d84be71fe19678_str_"
    },
    "v10": {
        "dir": r"I:\new_factor_tif\9ce01cf43be7606676d84be71fe19678\v10\v10",
        "prefix": "9ce01cf43be7606676d84be71fe19678_v10_"
    },
    "t2m": {
        "dir": r"I:\new_factor_tif1\de34709100aa97f6f68a02131b8fde57\t2m\t2m",
        "prefix": "de34709100aa97f6f68a02131b8fde57_t2m_"
    },
}


# =============================================================================
# 4. Template globals for multiprocessing workers
# =============================================================================
TEMPLATE_CRS = None
TEMPLATE_TRANSFORM = None
TEMPLATE_WIDTH = None
TEMPLATE_HEIGHT = None


def log_msg(message):
    print(message)
    with open(RUN_LOG, "a", encoding="utf-8") as file:
        file.write(str(message) + "\n")


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
    with open(out_txt, "w", encoding="utf-8") as file:
        file.write(text)
    ensure_file_saved(out_txt, "TXT")


def init_worker_template(crs, transform, width, height):
    global TEMPLATE_CRS, TEMPLATE_TRANSFORM, TEMPLATE_WIDTH, TEMPLATE_HEIGHT
    TEMPLATE_CRS = crs
    TEMPLATE_TRANSFORM = transform
    TEMPLATE_WIDTH = width
    TEMPLATE_HEIGHT = height


# =============================================================================
# 5. Paths and raster reading utilities
# =============================================================================
def tif_path(var_name, year, month):
    if var_name in ERA_STYLE_VARS:
        return os.path.join(ERA_BASE_DIR, ERA_STYLE_VARS[var_name], f"{year}_{month:02d}.tif")

    if var_name in SPECIAL_MONTHLY_VARS:
        spec = SPECIAL_MONTHLY_VARS[var_name]
        ymd = f"{year}{month:02d}01"
        return os.path.join(spec["dir"], f"{spec['prefix']}{ymd}.tif")

    raise KeyError(f"Unknown variable: {var_name}")


def get_template_info():
    path = tif_path(X_VAR, YEARS[0], MONTHS[0])
    if not os.path.exists(path):
        raise FileNotFoundError(f"Template file does not exist: {path}")

    with rasterio.open(path) as dataset:
        info = {
            "path": path,
            "width": dataset.width,
            "height": dataset.height,
            "transform": dataset.transform,
            "crs": dataset.crs,
            "count": dataset.count,
            "profile": dataset.profile.copy()
        }
    return info


def check_inputs(template_info):
    all_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))

    if CHECK_ALL_FILES:
        missing_files = []
        for var in all_vars:
            for year, month in TIME_STEPS:
                path = tif_path(var, year, month)
                if not os.path.exists(path):
                    missing_files.append(path)

        if missing_files:
            print("Missing files are listed below. Only the first 30 files are shown:")
            for path in missing_files[:30]:
                print("  ", path)
            raise FileNotFoundError(f"A total of {len(missing_files)} files are missing.")

    if CHECK_ALIGNMENT:
        mismatched = []
        for var in all_vars:
            path = tif_path(var, YEARS[0], MONTHS[0])
            with rasterio.open(path) as dataset:
                if dataset.count != 1:
                    raise ValueError(f"{var} is not a single-band raster: {path}")

                same_grid = (
                    dataset.width == template_info["width"] and
                    dataset.height == template_info["height"] and
                    dataset.transform == template_info["transform"] and
                    dataset.crs == template_info["crs"]
                )

                if not same_grid:
                    mismatched.append((var, path))

        if mismatched:
            print(
                f"{len(mismatched)} variables do not match the template grid. "
                "They will be aligned to the template during reading."
            )
            for var, path in mismatched[:10]:
                print(f"  Auto-aligning: {var} -> {path}")


def block_ranges(height, block_rows):
    ranges = []
    row_start = 0
    while row_start < height:
        row_end = min(row_start + block_rows, height)
        ranges.append((row_start, row_end))
        row_start = row_end
    return ranges


def get_resampling_method(var_name):
    return Resampling.nearest if var_name == "snowc" else Resampling.bilinear


def read_one_block(var_name, path, row_start, row_end, width):
    block_height = row_end - row_start
    window = Window(col_off=0, row_off=row_start, width=width, height=block_height)
    resampling = get_resampling_method(var_name)

    with rasterio.open(path) as dataset:
        same_grid = (
            dataset.width == TEMPLATE_WIDTH and
            dataset.height == TEMPLATE_HEIGHT and
            dataset.transform == TEMPLATE_TRANSFORM and
            dataset.crs == TEMPLATE_CRS
        )

        if same_grid:
            arr = dataset.read(1, window=window).astype(np.float32)
            nodata = dataset.nodata
            if nodata is not None:
                arr[arr == nodata] = np.nan
            return arr

        with WarpedVRT(
            dataset,
            crs=TEMPLATE_CRS,
            transform=TEMPLATE_TRANSFORM,
            width=TEMPLATE_WIDTH,
            height=TEMPLATE_HEIGHT,
            resampling=resampling
        ) as vrt:
            arr = vrt.read(1, window=window).astype(np.float32)
            nodata = vrt.nodata if vrt.nodata is not None else dataset.nodata
            if nodata is not None:
                arr[arr == nodata] = np.nan
            return arr


def load_var_cube(var_name, row_start, row_end, width):
    block_height = row_end - row_start
    cube = np.empty((len(TIME_STEPS), block_height, width), dtype=np.float32)

    for idx, (year, month) in enumerate(TIME_STEPS):
        cube[idx] = read_one_block(
            var_name,
            tif_path(var_name, year, month),
            row_start,
            row_end,
            width
        )

    return cube


# =============================================================================
# 6. Time-series preprocessing
# =============================================================================
def remove_monthly_climatology_3d(cube):
    output = np.empty_like(cube, dtype=np.float32)
    for month in range(1, 13):
        idx = np.where(MONTH_IDS == month)[0]
        climatology = np.nanmean(cube[idx], axis=0)
        output[idx] = cube[idx] - climatology
    return output


def remove_monthly_climatology_1d(series):
    series = np.asarray(series, dtype=np.float64)
    output = np.empty_like(series, dtype=np.float64)

    for month in range(1, 13):
        idx = np.where(MONTH_IDS == month)[0]
        climatology = np.nanmean(series[idx])
        output[idx] = series[idx] - climatology

    return output


def ar1_prewhiten_2d(arr_2d):
    valid_full = np.all(np.isfinite(arr_2d), axis=0)
    time_len, n_pixels = arr_2d.shape

    prewhitened = np.full((time_len - 1, n_pixels), np.nan, dtype=np.float32)
    phi = np.full(n_pixels, np.nan, dtype=np.float32)

    if not np.any(valid_full):
        return prewhitened, phi

    values = arr_2d[:, valid_full].astype(np.float64)
    centered = values - np.mean(values, axis=0, keepdims=True)

    x0 = centered[:-1]
    x1 = centered[1:]

    denominator = np.sum(x0 * x0, axis=0)
    numerator = np.sum(x0 * x1, axis=0)

    phi_sub = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > 0
    )
    phi_sub = np.clip(phi_sub, -0.99, 0.99)

    pw_sub = values[1:] - phi_sub[None, :] * values[:-1]

    prewhitened[:, valid_full] = pw_sub.astype(np.float32)
    phi[valid_full] = phi_sub.astype(np.float32)

    return prewhitened, phi


def ar1_prewhiten_1d(series):
    series = np.asarray(series, dtype=np.float64)
    valid = np.isfinite(series)

    if not np.all(valid):
        return np.full(len(series) - 1, np.nan, dtype=np.float64), np.nan

    centered = series - np.mean(series)
    x0 = centered[:-1]
    x1 = centered[1:]

    denominator = np.sum(x0 * x0)
    numerator = np.sum(x0 * x1)

    phi = 0.0 if denominator <= 0 else numerator / denominator
    phi = float(np.clip(phi, -0.99, 0.99))

    prewhitened = series[1:] - phi * series[:-1]
    return prewhitened, phi


# =============================================================================
# 7. Partial-correlation calculation
# =============================================================================
def effective_sample_size_from_ar1(series):
    """
    Estimate the effective sample size using the lag-1 autocorrelation coefficient.

    The estimate is bounded between 10 and the number of valid observations to avoid
    unrealistically small or large degrees of freedom in the pixel-wise t test.
    """
    values = np.asarray(series, dtype=np.float64)
    values = values[np.isfinite(values)]

    n = len(values)
    if n < 10:
        return float(n)

    centered = values - np.mean(values)
    x0 = centered[:-1]
    x1 = centered[1:]

    denominator = np.sum(x0 * x0)
    numerator = np.sum(x0 * x1)

    rho1 = 0.0 if denominator <= 0 else numerator / denominator
    rho1 = float(np.clip(rho1, -0.99, 0.99))

    n_eff = n * (1.0 - rho1) / (1.0 + rho1)
    return float(max(10.0, min(float(n), n_eff)))


def partial_corr_and_p_batch(x, y, z, pixel_chunk=2000):
    """
    Calculate pixel-wise partial correlation and p values.

    The original anomaly series are not AR(1)-prewhitened. Instead, the p value is
    computed using an AR(1)-based effective sample size for each pixel, and the
    effective sample size is exported as a spatial raster.
    """
    time_len, n_pixels = x.shape
    n_covariates = z.shape[1]

    r_out = np.full(n_pixels, np.nan, dtype=np.float32)
    p_out = np.full(n_pixels, np.nan, dtype=np.float32)
    n_out = np.full(n_pixels, np.nan, dtype=np.float32)

    valid = (
        np.all(np.isfinite(x), axis=0) &
        np.all(np.isfinite(y), axis=0) &
        np.all(np.isfinite(z), axis=(0, 1))
    )

    valid_idx = np.where(valid)[0]
    if valid_idx.size == 0:
        return r_out, p_out, n_out

    n_chunks = math.ceil(valid_idx.size / pixel_chunk)

    for chunk_idx in range(n_chunks):
        sub = valid_idx[chunk_idx * pixel_chunk:(chunk_idx + 1) * pixel_chunk]
        if sub.size == 0:
            continue

        x0 = x[:, sub].T.astype(np.float64)
        y0 = y[:, sub].T.astype(np.float64)
        z0 = np.transpose(z[:, :, sub], (2, 0, 1)).astype(np.float64)

        ones = np.ones((z0.shape[0], z0.shape[1], 1), dtype=np.float64)
        design = np.concatenate([ones, z0], axis=2)

        xtx = np.einsum("ntp,ntq->npq", design, design)
        xtx_inv = np.linalg.pinv(xtx)

        xtx_target_x = np.einsum("ntp,nt->np", design, x0)
        xtx_target_y = np.einsum("ntp,nt->np", design, y0)

        beta_x = np.einsum("npq,nq->np", xtx_inv, xtx_target_x)
        beta_y = np.einsum("npq,nq->np", xtx_inv, xtx_target_y)

        residual_x = x0 - np.einsum("ntp,np->nt", design, beta_x)
        residual_y = y0 - np.einsum("ntp,np->nt", design, beta_y)

        residual_x = residual_x - residual_x.mean(axis=1, keepdims=True)
        residual_y = residual_y - residual_y.mean(axis=1, keepdims=True)

        numerator = np.sum(residual_x * residual_y, axis=1)
        denominator = np.sqrt(np.sum(residual_x * residual_x, axis=1) * np.sum(residual_y * residual_y, axis=1))

        r = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan),
            where=denominator > 0
        )
        r = np.clip(r, -0.999999, 0.999999)

        for local_idx, pixel_idx in enumerate(sub):
            if not np.isfinite(r[local_idx]):
                continue

            if USE_EFFECTIVE_SAMPLE_SIZE:
                n_eff_x = effective_sample_size_from_ar1(residual_x[local_idx, :])
                n_eff_y = effective_sample_size_from_ar1(residual_y[local_idx, :])
                n_eff = min(n_eff_x, n_eff_y)
            else:
                n_eff = float(time_len)

            df_eff = n_eff - n_covariates - 2
            n_out[pixel_idx] = n_eff

            if df_eff > 0 and abs(r[local_idx]) < 1:
                t_stat = r[local_idx] * np.sqrt(df_eff / (1.0 - r[local_idx] * r[local_idx]))
                p_out[pixel_idx] = 2.0 * stats.t.sf(np.abs(t_stat), df_eff)

        r_out[sub] = r.astype(np.float32)

    return r_out, p_out, n_out

def residualize_batch(x, y, z):
    """
    Residualize x and y against z for multiple pixels.

    Parameters
    ----------
    x : ndarray, shape (T, N)
        Predictor time series.
    y : ndarray, shape (T, N)
        Response time series.
    z : ndarray, shape (T, K, N)
        Covariate time series.

    Returns
    -------
    residual_x : ndarray, shape (N, T)
    residual_y : ndarray, shape (N, T)
    """
    x0 = x.T.astype(np.float64)
    y0 = y.T.astype(np.float64)
    z0 = np.transpose(z, (2, 0, 1)).astype(np.float64)

    ones = np.ones((z0.shape[0], z0.shape[1], 1), dtype=np.float64)
    design = np.concatenate([ones, z0], axis=2)

    xtx = np.einsum("ntp,ntq->npq", design, design)
    xtx_inv = np.linalg.pinv(xtx)

    xtx_target_x = np.einsum("ntp,nt->np", design, x0)
    xtx_target_y = np.einsum("ntp,nt->np", design, y0)

    beta_x = np.einsum("npq,nq->np", xtx_inv, xtx_target_x)
    beta_y = np.einsum("npq,nq->np", xtx_inv, xtx_target_y)

    residual_x = x0 - np.einsum("ntp,np->nt", design, beta_x)
    residual_y = y0 - np.einsum("ntp,np->nt", design, beta_y)

    residual_x = residual_x - residual_x.mean(axis=1, keepdims=True)
    residual_y = residual_y - residual_y.mean(axis=1, keepdims=True)

    return residual_x, residual_y


def partial_corr_1d(x, y, z):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)

    valid = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(z), axis=1)
    x = x[valid]
    y = y[valid]
    z = z[valid]

    n = len(x)
    k = z.shape[1] if z.ndim == 2 else 0
    df = n - k - 2

    if n == 0 or df <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.array([]), np.array([]), np.nan, np.nan, 0

    design = np.column_stack([np.ones(n), z])

    beta_x, *_ = np.linalg.lstsq(design, x, rcond=None)
    beta_y, *_ = np.linalg.lstsq(design, y, rcond=None)

    residual_x = x - design @ beta_x
    residual_y = y - design @ beta_y

    residual_x = residual_x - np.mean(residual_x)
    residual_y = residual_y - np.mean(residual_y)

    numerator = np.sum(residual_x * residual_y)
    denominator = np.sqrt(np.sum(residual_x * residual_x) * np.sum(residual_y * residual_y))

    if denominator <= 0:
        return np.nan, np.nan, np.nan, np.nan, residual_x, residual_y, np.nan, np.nan, n

    r = float(np.clip(numerator / denominator, -0.999999, 0.999999))
    t_stat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2.0 * stats.t.sf(abs(t_stat), df))

    if n - k - 3 > 0:
        z_value = np.arctanh(r)
        se = 1.0 / np.sqrt(n - k - 3)
        z_low = z_value - 1.96 * se
        z_high = z_value + 1.96 * se
        ci_low = float(np.tanh(z_low))
        ci_high = float(np.tanh(z_high))
    else:
        ci_low, ci_high = np.nan, np.nan

    slope, intercept = np.polyfit(residual_x, residual_y, 1)

    return r, p, ci_low, ci_high, residual_x, residual_y, float(slope), float(intercept), n


# =============================================================================
# 8. Statistics and diagnostics
# =============================================================================
def summarize_pixel_distribution(r_vals, p_vals, n_vals, label, alpha=0.05):
    r_vals = np.asarray(r_vals, dtype=np.float64)
    p_vals = np.asarray(p_vals, dtype=np.float64)
    n_vals = np.asarray(n_vals, dtype=np.float64)

    valid = np.isfinite(r_vals) & np.isfinite(p_vals)
    r = r_vals[valid]
    p = p_vals[valid]
    n_valid = n_vals[np.isfinite(n_vals)]

    if len(r) == 0:
        return {
            "region": label,
            "valid_pixels": 0,
            "significant_pixels": 0,
            "significant_ratio": np.nan,
            "significant_positive_pixels": 0,
            "significant_positive_ratio": np.nan,
            "significant_negative_pixels": 0,
            "significant_negative_ratio": np.nan,
            "non_significant_pixels": 0,
            "non_significant_ratio": np.nan,
            "mean_r": np.nan,
            "median_r": np.nan,
            "std_r": np.nan,
            "min_r": np.nan,
            "max_r": np.nan,
            "p25_r": np.nan,
            "p75_r": np.nan,
            "mean_abs_r": np.nan,
            "mean_effective_n": np.nan,
            "median_effective_n": np.nan,
            "min_effective_n": np.nan,
            "max_effective_n": np.nan
        }

    significant = p < alpha
    significant_positive = significant & (r > 0)
    significant_negative = significant & (r < 0)

    return {
        "region": label,
        "valid_pixels": int(len(r)),
        "significant_pixels": int(np.sum(significant)),
        "significant_ratio": float(np.mean(significant)),
        "significant_positive_pixels": int(np.sum(significant_positive)),
        "significant_positive_ratio": float(np.mean(significant_positive)),
        "significant_negative_pixels": int(np.sum(significant_negative)),
        "significant_negative_ratio": float(np.mean(significant_negative)),
        "non_significant_pixels": int(np.sum(~significant)),
        "non_significant_ratio": float(np.mean(~significant)),
        "mean_r": float(np.mean(r)),
        "median_r": float(np.median(r)),
        "std_r": float(np.std(r, ddof=0)),
        "min_r": float(np.min(r)),
        "max_r": float(np.max(r)),
        "p25_r": float(np.percentile(r, 25)),
        "p75_r": float(np.percentile(r, 75)),
        "mean_abs_r": float(np.mean(np.abs(r))),
        "mean_effective_n": float(np.mean(n_valid)) if len(n_valid) > 0 else np.nan,
        "median_effective_n": float(np.median(n_valid)) if len(n_valid) > 0 else np.nan,
        "min_effective_n": float(np.min(n_valid)) if len(n_valid) > 0 else np.nan,
        "max_effective_n": float(np.max(n_valid)) if len(n_valid) > 0 else np.nan
    }


def fisher_ci_from_r(r, n):
    if (not np.isfinite(r)) or (n is None) or (n <= 3) or abs(r) >= 1:
        return np.nan, np.nan

    z_value = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    z_low = z_value - 1.96 * se
    z_high = z_value + 1.96 * se

    return float(np.tanh(z_low)), float(np.tanh(z_high))


def check_multicollinearity(z, var_names):
    z = np.asarray(z, dtype=np.float64)
    z = z[np.all(np.isfinite(z), axis=1)]

    if z.shape[0] == 0:
        return pd.DataFrame({"variable": var_names, "VIF": [np.nan] * len(var_names)})

    means = np.mean(z, axis=0)
    stds = np.std(z, axis=0, ddof=0)
    stds[stds == 0] = 1.0
    z_scaled = (z - means) / stds

    vif_values = []

    for idx in range(z_scaled.shape[1]):
        target = z_scaled[:, idx]
        others = np.delete(z_scaled, idx, axis=1)

        try:
            if others.shape[1] == 0:
                vif = 1.0
            else:
                beta = np.linalg.lstsq(others, target, rcond=None)[0]
                predicted = others @ beta

                if np.std(target) == 0 or np.std(predicted) == 0:
                    r2 = 0.0
                else:
                    r2 = np.corrcoef(target, predicted)[0, 1] ** 2

                r2 = min(max(r2, 0.0), 0.999999)
                vif = 1.0 / (1.0 - r2)

        except Exception:
            vif = np.nan

        vif_values.append(vif)

    return pd.DataFrame({"variable": var_names, "VIF": vif_values})


def spatial_autocorr_test(r_map, max_points=2000):
    valid = np.isfinite(r_map)

    if np.sum(valid) < 100:
        return None, None

    coords = np.array(np.where(valid)).T
    values = r_map[valid].astype(np.float64)

    if len(values) > max_points:
        rng = np.random.default_rng(42)
        selected = rng.choice(len(values), size=max_points, replace=False)
        coords = coords[selected]
        values = values[selected]

    if len(values) < 50:
        return None, None

    dist_matrix = cdist(coords, coords)
    weight_matrix = (dist_matrix < 1.5) & (dist_matrix > 0)

    if np.sum(weight_matrix) == 0:
        return None, None

    n = len(values)
    centered = values - np.mean(values)
    numerator = np.sum(weight_matrix * np.outer(centered, centered))
    denominator = np.sum(centered ** 2)
    s0 = np.sum(weight_matrix)

    if denominator == 0 or s0 == 0:
        return None, None

    moran_i = (n / s0) * (numerator / denominator)

    try:
        expected_i = -1.0 / (n - 1)
        variance_i = (
            n ** 2 * np.sum((weight_matrix + weight_matrix.T) ** 2)
            - n * np.sum((np.sum(weight_matrix, axis=1) + np.sum(weight_matrix, axis=0)) ** 2)
            + s0 ** 2
        ) / (s0 ** 2 * (n ** 2 - 1))

        if variance_i <= 0 or not np.isfinite(variance_i):
            return moran_i, np.nan

        z_score = (moran_i - expected_i) / np.sqrt(variance_i)
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    except Exception:
        p_value = np.nan

    return moran_i, p_value


def make_integrated_conclusion(overall_stats, stratified_df):
    lines = [
        "Integrated summary of pixel-wise partial correlation",
        "=" * 60,
        "",
        "1. Overall snow-domain result",
        (
            f"The snow-domain analysis included {overall_stats['valid_pixels']} valid pixels. "
            f"Among them, {overall_stats['significant_pixels']} pixels were significant "
            f"({overall_stats['significant_ratio'] * 100:.2f}%)."
        ),
        (
            f"Significant positive pixels accounted for "
            f"{overall_stats['significant_positive_ratio'] * 100:.2f}%, whereas significant "
            f"negative pixels accounted for {overall_stats['significant_negative_ratio'] * 100:.2f}%."
        ),
        (
            f"The mean partial correlation coefficient was {overall_stats['mean_r']:.4f}, "
            f"the median was {overall_stats['median_r']:.4f}, and the interquartile range was "
            f"[{overall_stats['p25_r']:.4f}, {overall_stats['p75_r']:.4f}]."
        ),
        ""
    ]

    if stratified_df.shape[0] > 0:
        lines.append("2. Stratified results by snow-frequency zone")
        for _, row in stratified_df.iterrows():
            region_name = row.get("region_name", row.get("region", "Unknown zone"))
            valid_pixels = int(row.get("valid_pixels", 0))

            if valid_pixels == 0:
                lines.append(f"{region_name}: no valid pixels.")
            else:
                lines.append(
                    f"{region_name}: {valid_pixels} valid pixels, "
                    f"{row['significant_ratio'] * 100:.2f}% significant pixels, "
                    f"mean r = {row['mean_r']:.4f}, median r = {row['median_r']:.4f}."
                )

    return "\n".join(lines)


# =============================================================================
# 9. Raster and quicklook output functions
# =============================================================================
def save_raster(out_tif, arr, profile, nodata=NODATA_OUT, dtype="float32"):
    profile_out = profile.copy()
    profile_out.update(dtype=dtype, count=1, nodata=nodata, compress=COMPRESS)

    if np.issubdtype(np.dtype(dtype), np.integer):
        out_arr = np.where(np.isfinite(arr), arr, nodata).astype(dtype)
    else:
        out_arr = np.where(np.isfinite(arr), arr.astype(np.float32), nodata).astype(dtype)

    with rasterio.open(out_tif, "w", **profile_out) as dst:
        dst.write(out_arr, 1)

    ensure_file_saved(out_tif, "GeoTIFF")


def save_corr_png(out_png, arr, title):
    plt.figure(figsize=(10, 6))
    masked = np.ma.masked_invalid(arr)
    image = plt.imshow(masked, cmap=CORR_CMAP, vmin=R_VMIN, vmax=R_VMAX)
    colorbar = plt.colorbar(image, fraction=0.035, pad=0.03)
    colorbar.set_label("Partial correlation coefficient (r)", fontsize=11)
    apply_arial_to_colorbar(colorbar)

    plt.title(title, fontsize=13)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    ensure_file_saved(out_png, "PNG")


def save_corr_png_with_sig_points(out_png, r_arr, p_arr, title):
    plt.figure(figsize=(10, 6))
    masked = np.ma.masked_invalid(r_arr)
    image = plt.imshow(masked, cmap=CORR_CMAP, vmin=R_VMIN, vmax=R_VMAX)
    colorbar = plt.colorbar(image, fraction=0.035, pad=0.03)
    colorbar.set_label("Partial correlation coefficient (r)", fontsize=11)
    apply_arial_to_colorbar(colorbar)

    significant = np.isfinite(p_arr) & (p_arr < ALPHA)
    rows, cols = np.where(significant)

    if len(rows) > 0:
        keep = (rows % SIG_POINT_STRIDE == 0) & (cols % SIG_POINT_STRIDE == 0)
        plt.scatter(
            cols[keep],
            rows[keep],
            s=SIG_POINT_SIZE,
            c="k",
            alpha=0.7,
            linewidths=0
        )

    plt.title(title, fontsize=13)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    ensure_file_saved(out_png, "PNG")


def save_p_png(out_png, p_arr, title):
    plt.figure(figsize=(10, 6))
    masked = np.ma.masked_invalid(p_arr)
    cmap = plt.cm.viridis_r.copy()
    cmap.set_bad(color="white")

    image = plt.imshow(masked, cmap=cmap, vmin=P_VMIN, vmax=P_VMAX)
    colorbar = plt.colorbar(image, fraction=0.035, pad=0.03)
    colorbar.set_label("p-value", fontsize=11)

    plt.title(title, fontsize=13)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    ensure_file_saved(out_png, "PNG")


def save_n_png(out_png, n_arr, title):
    plt.figure(figsize=(10, 6))
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

    image = plt.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax)
    colorbar = plt.colorbar(image, fraction=0.035, pad=0.03)
    colorbar.set_label("Effective sample size", fontsize=11)

    plt.title(title, fontsize=13)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    ensure_file_saved(out_png, "PNG")


def save_rpn_triptych(out_png, r_arr, p_arr, n_arr, main_title):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.8))

    masked_r = np.ma.masked_invalid(r_arr)
    image_r = axes[0].imshow(masked_r, cmap=CORR_CMAP, vmin=R_VMIN, vmax=R_VMAX)
    colorbar_r = fig.colorbar(image_r, ax=axes[0], fraction=0.046, pad=0.03)
    colorbar_r.set_label("r", fontsize=10)
    apply_arial_to_colorbar(colorbar_r)
    axes[0].set_title("Partial correlation", fontsize=11)
    axes[0].axis("off")

    masked_p = np.ma.masked_invalid(p_arr)
    cmap_p = plt.cm.viridis_r.copy()
    cmap_p.set_bad(color="white")
    image_p = axes[1].imshow(masked_p, cmap=cmap_p, vmin=P_VMIN, vmax=P_VMAX)
    colorbar_p = fig.colorbar(image_p, ax=axes[1], fraction=0.046, pad=0.03)
    colorbar_p.set_label("p", fontsize=10)
    axes[1].set_title("P-value", fontsize=11)
    axes[1].axis("off")

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

    image_n = axes[2].imshow(masked_n, cmap=cmap_n, vmin=vmin_n, vmax=vmax_n)
    colorbar_n = fig.colorbar(image_n, ax=axes[2], fraction=0.046, pad=0.03)
    colorbar_n.set_label("n", fontsize=10)
    axes[2].set_title("Effective sample size", fontsize=11)
    axes[2].axis("off")

    fig.suptitle(main_title, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    ensure_file_saved(out_png, "PNG")


# =============================================================================
# 10. Publication-style world map
# =============================================================================
def format_lon_dms(x):
    if abs(x) < 1e-9:
        return "0°0′0″"
    degree = int(abs(round(x)))
    hemisphere = "E" if x > 0 else "W"
    return f"{degree}°0′0″ {hemisphere}"


def format_lat_dms(y):
    if abs(y) < 1e-9:
        return "0°0′0″"
    degree = int(abs(round(y)))
    hemisphere = "N" if y > 0 else "S"
    return f"{degree}°0′0″ {hemisphere}"


def wrap_longitude_if_needed(arr, transform):
    height, width = arr.shape
    right = transform.c + transform.a * width

    if right > 180:
        split_col = width // 2
        arr_fixed = np.hstack([arr[:, split_col:], arr[:, :split_col]])
        shift = transform.a * split_col
        new_transform = Affine(
            transform.a, transform.b, transform.c - shift,
            transform.d, transform.e, transform.f
        )
        return arr_fixed, new_transform

    return arr.copy(), transform


def save_publication_world_map_tif(r_arr, profile):
    """
    Save the publication-style large-scale pixel partial-correlation map.

    Plot-specific requirements:
    1. No background grid lines.
    2. No longitude or latitude coordinates, labels, ticks, or tick marks.
    3. Arial font for the partial-correlation legend.
    4. Negative correlation color: RGB (32, 56, 136).
    5. Positive correlation color: RGB (225, 156, 102).
    """
    crs = profile["crs"]
    transform = profile["transform"]
    width = profile["width"]
    height = profile["height"]

    arr_fixed, new_transform = wrap_longitude_if_needed(r_arr, transform)

    out_r_corrected = os.path.join(
        OUT_PUB_MAP,
        f"{TARGET_VAR}_partial_corr_r_lonfix.tif"
    )
    out_fig = os.path.join(
        OUT_PUB_MAP,
        f"{TARGET_VAR}_partial_corr_world_map_pub.tif"
    )

    meta = profile.copy()
    meta.update({
        "transform": new_transform,
        "crs": crs,
        "dtype": "float32",
        "nodata": NODATA_OUT,
        "count": 1,
        "compress": COMPRESS
    })

    save_arr = np.where(
        np.isfinite(arr_fixed),
        arr_fixed.astype(np.float32),
        NODATA_OUT
    ).astype(np.float32)

    with rasterio.open(out_r_corrected, "w", **meta) as dst:
        dst.write(save_arr, 1)

    ensure_file_saved(out_r_corrected, "GeoTIFF")

    if not os.path.exists(WORLD_SHP):
        raise FileNotFoundError(
            f"World boundary shapefile does not exist: {WORLD_SHP}"
        )

    world = gpd.read_file(WORLD_SHP)
    if world.crs != crs:
        world = world.to_crs(crs)

    left = new_transform.c
    right = new_transform.c + new_transform.a * width
    top = new_transform.f
    bottom = new_transform.f + new_transform.e * height

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    ax = fig.add_axes([0.055, 0.12, 0.89, 0.78])
    ax.set_facecolor("white")

    world.plot(
        ax=ax,
        facecolor=BG_RGB,
        edgecolor=EDGE_RGB,
        linewidth=0.5,
        zorder=1
    )

    masked_r = np.ma.masked_invalid(arr_fixed)

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

    world.boundary.plot(
        ax=ax,
        color=EDGE_RGB,
        linewidth=0.5,
        zorder=3
    )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)

    # Remove background grid lines and all longitude/latitude information.
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

    # Retain the outer map frame.
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")

    ax.set_title("")

    cax = ax.inset_axes([0.39, 0.045, 0.22, 0.028])
    colorbar = plt.colorbar(
        image,
        cax=cax,
        orientation="horizontal",
        extend="both"
    )
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    colorbar.ax.tick_params(
        labelsize=11,
        direction="out",
        length=3,
        width=0.8,
        pad=2
    )
    colorbar.ax.set_title(
        "Partial correlation coefficient (r)",
        fontsize=12,
        pad=6,
        fontproperties=ARIAL_FONT
    )
    apply_arial_to_colorbar(colorbar)

    try:
        fig.savefig(
            out_fig,
            dpi=FIG_DPI,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
            pil_kwargs={"compression": "tiff_lzw"}
        )
    except TypeError:
        fig.savefig(
            out_fig,
            dpi=FIG_DPI,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor()
        )

    plt.close(fig)
    ensure_file_saved(out_fig, "TIF")

    return out_r_corrected, out_fig


# =============================================================================
# 11. Stratified analysis by snow-frequency zone
# =============================================================================
def stratified_analysis_by_snow_frequency(snow_freq_map, r_map, p_map, n_map, width, height, profile):
    zones = {
        "permanent": {
            "name": "Permanent snow zone",
            "mask": snow_freq_map >= PERMANENT_SNOW_MIN,
            "min_freq": PERMANENT_SNOW_MIN,
            "max_freq": 1.0
        },
        "seasonal": {
            "name": "Seasonal snow zone",
            "mask": (snow_freq_map >= SEASONAL_SNOW_MIN) & (snow_freq_map < PERMANENT_SNOW_MIN),
            "min_freq": SEASONAL_SNOW_MIN,
            "max_freq": PERMANENT_SNOW_MIN
        },
        "ephemeral": {
            "name": "Ephemeral snow zone",
            "mask": (snow_freq_map >= EPHEMERAL_SNOW_MIN) & (snow_freq_map < SEASONAL_SNOW_MIN),
            "min_freq": EPHEMERAL_SNOW_MIN,
            "max_freq": SEASONAL_SNOW_MIN
        }
    }

    results = {}

    for key, zone in zones.items():
        mask = zone["mask"]

        r_zone = r_map[mask]
        p_zone = p_map[mask]
        n_zone = n_map[mask]

        stats_dict = summarize_pixel_distribution(r_zone, p_zone, n_zone, zone["name"], alpha=ALPHA)
        stats_dict["region"] = key
        stats_dict["region_name"] = zone["name"]
        stats_dict["min_snow_frequency"] = zone["min_freq"]
        stats_dict["max_snow_frequency"] = zone["max_freq"]

        results[key] = stats_dict

        r_zone_map = np.full((height, width), np.nan, dtype=np.float32)
        p_zone_map = np.full((height, width), np.nan, dtype=np.float32)
        n_zone_map = np.full((height, width), np.nan, dtype=np.float32)
        sig_zone_map = np.full((height, width), 255, dtype=np.uint8)

        r_zone_map[mask] = r_zone
        p_zone_map[mask] = p_zone
        n_zone_map[mask] = n_zone

        sig_local = np.isfinite(p_zone) & (p_zone < ALPHA)
        sig_zone_map[mask] = sig_local.astype(np.uint8)

        out_r_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_r_{key}.tif")
        out_p_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_p_{key}.tif")
        out_n_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_effective_n_{key}.tif")
        out_sig_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_sigmask_{key}.tif")

        save_raster(out_r_tif, r_zone_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_p_tif, p_zone_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_n_tif, n_zone_map, profile, nodata=NODATA_OUT, dtype="float32")
        save_raster(out_sig_tif, sig_zone_map, profile, nodata=255, dtype="uint8")

        out_r_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_r_{key}.png")
        out_r_sig_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_r_{key}_sig_points.png")
        out_p_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_p_{key}.png")
        out_n_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_effective_n_{key}.png")
        out_rpn_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_rpn_triptych_{key}.png")

        save_corr_png(out_r_png, r_zone_map, f"Partial correlation: {TARGET_VAR} vs {X_VAR} ({zone['name']})")
        save_corr_png_with_sig_points(
            out_r_sig_png,
            r_zone_map,
            p_zone_map,
            f"Partial correlation: {TARGET_VAR} vs {X_VAR} ({zone['name']}, p < {ALPHA})"
        )
        save_p_png(out_p_png, p_zone_map, f"P-value: {TARGET_VAR} vs {X_VAR} ({zone['name']})")
        save_n_png(out_n_png, n_zone_map, f"Effective sample size: {TARGET_VAR} vs {X_VAR} ({zone['name']})")
        save_rpn_triptych(
            out_rpn_png,
            r_zone_map,
            p_zone_map,
            n_zone_map,
            f"Partial correlation result: {zone['name']}"
        )

    return results


# =============================================================================
# 12. Block-level worker
# =============================================================================
def process_one_block(args):
    row_start, row_end, width = args
    block_height = row_end - row_start

    needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))

    snowc_raw = load_var_cube("snowc", row_start, row_end, width)
    snow_present = np.isfinite(snowc_raw) & (snowc_raw > SNOWC_THRESHOLD)
    snow_freq_block = (np.sum(snow_present, axis=0) / float(len(TIME_STEPS))).astype(np.float32)

    if SNOWC_MASK_RULE == "frequency_ge_0_30":
        snow_mask_block = np.isfinite(snow_freq_block) & (snow_freq_block >= SNOWC_MIN_FRACTION)
    else:
        raise ValueError(f"Unsupported SNOWC_MASK_RULE: {SNOWC_MASK_RULE}")

    ts_sums = {var: np.zeros(len(TIME_STEPS), dtype=np.float64) for var in needed_vars}
    ts_counts = {var: np.zeros(len(TIME_STEPS), dtype=np.int64) for var in needed_vars}

    if not np.any(snow_mask_block):
        return (
            row_start, row_end,
            np.full((block_height, width), np.nan, dtype=np.float32),
            np.full((block_height, width), np.nan, dtype=np.float32),
            np.full((block_height, width), np.nan, dtype=np.float32),
            snow_freq_block,
            snow_mask_block.astype(np.uint8),
            ts_sums,
            ts_counts
        )

    data_proc = {}

    for var in needed_vars:
        cube_raw = snowc_raw if var == "snowc" else load_var_cube(var, row_start, row_end, width)

        values_in_domain = cube_raw[:, snow_mask_block]
        ts_sums[var] = np.nansum(values_in_domain, axis=1).astype(np.float64)
        ts_counts[var] = np.sum(np.isfinite(values_in_domain), axis=1).astype(np.int64)

        cube_proc = remove_monthly_climatology_3d(cube_raw) if REMOVE_MONTHLY_CLIM else cube_raw
        arr_2d = cube_proc.reshape(len(TIME_STEPS), -1)

        if USE_AR1_PREWHITEN:
            arr_proc, _phi = ar1_prewhiten_2d(arr_2d)
            del _phi
        else:
            arr_proc = arr_2d.astype(np.float32)

        data_proc[var] = arr_proc

        del values_in_domain
        if var != "snowc":
            del cube_raw
        del cube_proc
        del arr_2d

    del snowc_raw

    x = data_proc[X_VAR]
    y = data_proc[TARGET_VAR]
    z = np.stack([data_proc[covar] for covar in COVARS], axis=1)

    r_flat, p_flat, n_flat = partial_corr_and_p_batch(
        x=x,
        y=y,
        z=z,
        pixel_chunk=PIXEL_CHUNK
    )

    mask_flat = snow_mask_block.reshape(-1)
    r_flat[~mask_flat] = np.nan
    p_flat[~mask_flat] = np.nan
    n_flat[~mask_flat] = np.nan

    return (
        row_start, row_end,
        r_flat.reshape(block_height, width).astype(np.float32),
        p_flat.reshape(block_height, width).astype(np.float32),
        n_flat.reshape(block_height, width).astype(np.float32),
        snow_freq_block,
        snow_mask_block.astype(np.uint8),
        ts_sums,
        ts_counts
    )


# =============================================================================
# 13. Publication-style partial-regression scatter plots
# =============================================================================
def save_figure_tif(fig, out_tif, dpi=600):
    try:
        fig.savefig(
            out_tif,
            dpi=dpi,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
            pil_kwargs={"compression": "tiff_lzw"}
        )
    except TypeError:
        fig.savefig(
            out_tif,
            dpi=dpi,
            format="tif",
            bbox_inches="tight",
            facecolor=fig.get_facecolor()
        )

    ensure_file_saved(out_tif, "TIF")


def init_scatter_acc(title, out_tif, out_csv):
    return {
        "title": title,
        "out_tif": out_tif,
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


def update_scatter_acc(acc, x_vals, y_vals, rng):
    if x_vals.size == 0:
        return

    valid = np.isfinite(x_vals) & np.isfinite(y_vals)
    x = x_vals[valid].astype(np.float64, copy=False)
    y = y_vals[valid].astype(np.float64, copy=False)

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
        selected = rng.choice(x.size, size=keep, replace=False)
        xs = x[selected].astype(np.float32)
        ys = y[selected].astype(np.float32)
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
        sample_x = np.concatenate(acc["sample_x_chunks"])
        sample_y = np.concatenate(acc["sample_y_chunks"])

        if sample_x.size > MAX_PLOT_POINTS_PER_GROUP:
            rng = np.random.default_rng(SCATTER_RANDOM_STATE)
            selected = rng.choice(sample_x.size, size=MAX_PLOT_POINTS_PER_GROUP, replace=False)
            sample_x = sample_x[selected]
            sample_y = sample_y[selected]
    else:
        sample_x = np.array([], dtype=np.float32)
        sample_y = np.array([], dtype=np.float32)

    return {
        "title": acc["title"],
        "out_tif": acc["out_tif"],
        "out_csv": acc["out_csv"],
        "n": int(n),
        "r": r,
        "p": p,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "slope": slope,
        "intercept": intercept,
        "sample_x": sample_x,
        "sample_y": sample_y
    }


def format_p_value(p):
    if not np.isfinite(p):
        return "NaN"
    if p < 0.001:
        return "< 0.001"
    return f"{p:.3f}"


def residual_axis_suffix():
    if USE_AR1_PREWHITEN:
        return "after covariate adjustment"
    return "after covariate adjustment"


def save_partial_scatter_tif(out_tif, sample_x, sample_y, r, p, n, slope, intercept, title):
    sample_x = np.asarray(sample_x, dtype=np.float64).ravel()
    sample_y = np.asarray(sample_y, dtype=np.float64).ravel()

    valid = np.isfinite(sample_x) & np.isfinite(sample_y)
    sample_x = sample_x[valid]
    sample_y = sample_y[valid]

    fig, ax = plt.subplots(figsize=(11.2, 6.2), facecolor="white")

    if sample_x.size > 0 and sample_y.size > 0:
        ax.scatter(
            sample_x,
            sample_y,
            s=18,
            facecolor=SCATTER_FACE,
            edgecolor=SCATTER_EDGE,
            linewidth=0.35,
            alpha=0.60
        )

        if np.isfinite(slope) and np.isfinite(intercept):
            xs = np.linspace(np.nanmin(sample_x), np.nanmax(sample_x), 200)
            ys = slope * xs + intercept
            ax.plot(xs, ys, color=REG_LINE, linewidth=2.8)

    text = f"r = {r:.3f}\np = {format_p_value(p)}\nn = {int(n):,}"

    ax.text(
        0.05,
        0.95,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12
    )

    ax.set_xlabel(f"Residualized {X_VAR} ({residual_axis_suffix()})", fontsize=13)
    ax.set_ylabel(f"Residualized {TARGET_VAR} ({residual_axis_suffix()})", fontsize=13)
    ax.set_title(title, fontsize=14)

    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(1.0)

    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=4,
        width=1.0,
        labelsize=11
    )

    ax.grid(False)

    save_figure_tif(fig, out_tif, dpi=FIG_DPI)
    plt.close(fig)


def generate_scatter_outputs(template, r_map, p_map, snow_freq_map):
    valid_pixels = np.isfinite(r_map) & np.isfinite(p_map)
    sig_pos = valid_pixels & (p_map < ALPHA) & (r_map > 0)
    sig_neg = valid_pixels & (p_map < ALPHA) & (r_map < 0)
    nonsig = valid_pixels & ~(sig_pos | sig_neg)

    permanent = valid_pixels & (snow_freq_map >= PERMANENT_SNOW_MIN)
    seasonal = valid_pixels & (snow_freq_map >= SEASONAL_SNOW_MIN) & (snow_freq_map < PERMANENT_SNOW_MIN)
    ephemeral = valid_pixels & (snow_freq_map >= EPHEMERAL_SNOW_MIN) & (snow_freq_map < SEASONAL_SNOW_MIN)

    groups = {
        "overall": {
            "mask": valid_pixels,
            "title": "All valid pixels",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_overall_all_valid_pixels.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_overall_all_valid_pixels.csv"),
        },
        "significant_positive": {
            "mask": sig_pos,
            "title": "Significant positive pixels",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_significant_positive_pixels.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_significant_positive_pixels.csv"),
        },
        "significant_negative": {
            "mask": sig_neg,
            "title": "Significant negative pixels",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_significant_negative_pixels.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_significant_negative_pixels.csv"),
        },
        "non_significant": {
            "mask": nonsig,
            "title": "Non-significant pixels",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_non_significant_pixels.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_non_significant_pixels.csv"),
        },
        "permanent_zone": {
            "mask": permanent,
            "title": "Permanent snow zone",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_permanent_snow_zone.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_permanent_snow_zone.csv"),
        },
        "seasonal_zone": {
            "mask": seasonal,
            "title": "Seasonal snow zone",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_seasonal_snow_zone.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_seasonal_snow_zone.csv"),
        },
        "ephemeral_zone": {
            "mask": ephemeral,
            "title": "Ephemeral snow zone",
            "tif": os.path.join(OUT_SCATTER_FIG, "scatter_ephemeral_snow_zone.tif"),
            "csv": os.path.join(OUT_SCATTER_TABLE, "scatter_ephemeral_snow_zone.csv"),
        },
    }

    accs = {key: init_scatter_acc(value["title"], value["tif"], value["csv"]) for key, value in groups.items()}
    rng = np.random.default_rng(SCATTER_RANDOM_STATE)

    height = template["height"]
    width = template["width"]
    ranges = block_ranges(height, BLOCK_ROWS)

    log_msg("Generating publication-style partial-regression scatter plots...")

    for row_start, row_end in tqdm(ranges, desc="scatter", ncols=100):
        block_valid = valid_pixels[row_start:row_end, :]
        if not np.any(block_valid):
            continue

        needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))
        processed = {}

        for var in needed_vars:
            cube_raw = load_var_cube(var, row_start, row_end, width)
            cube_proc = remove_monthly_climatology_3d(cube_raw) if REMOVE_MONTHLY_CLIM else cube_raw
            arr_2d = cube_proc.reshape(len(TIME_STEPS), -1)

            if USE_AR1_PREWHITEN:
                arr_proc, _phi = ar1_prewhiten_2d(arr_2d)
                del _phi
            else:
                arr_proc = arr_2d.astype(np.float32)

            processed[var] = arr_proc

            del cube_raw
            del cube_proc
            del arr_2d

        x_all = processed[X_VAR]
        y_all = processed[TARGET_VAR]
        z_all = np.stack([processed[covar] for covar in COVARS], axis=1)

        block_valid_flat = block_valid.reshape(-1)
        idx0 = np.where(block_valid_flat)[0]
        if idx0.size == 0:
            del processed, x_all, y_all, z_all
            continue

        strict_ok = (
            np.all(np.isfinite(x_all[:, idx0]), axis=0) &
            np.all(np.isfinite(y_all[:, idx0]), axis=0) &
            np.all(np.isfinite(z_all[:, :, idx0]), axis=(0, 1))
        )
        idx = idx0[strict_ok]
        if idx.size == 0:
            del processed, x_all, y_all, z_all
            continue

        x_sel = x_all[:, idx]
        y_sel = y_all[:, idx]
        z_sel = z_all[:, :, idx]

        residual_x, residual_y = residualize_batch(x_sel, y_sel, z_sel)

        for group_name, group_info in groups.items():
            group_mask_flat = group_info["mask"][row_start:row_end, :].reshape(-1)
            membership = group_mask_flat[idx]

            if not np.any(membership):
                continue

            x_group = residual_x[membership].reshape(-1)
            y_group = residual_y[membership].reshape(-1)
            update_scatter_acc(accs[group_name], x_group, y_group, rng)

        del processed, x_all, y_all, z_all, x_sel, y_sel, z_sel, residual_x, residual_y

    rows = []

    for group_name, acc in accs.items():
        result = finalize_scatter_acc(acc)

        save_partial_scatter_tif(
            result["out_tif"],
            result["sample_x"],
            result["sample_y"],
            result["r"],
            result["p"],
            result["n"],
            result["slope"],
            result["intercept"],
            result["title"]
        )

        row = {
            "group": group_name,
            "title": result["title"],
            "figure_file": os.path.basename(result["out_tif"]),
            "n": result["n"],
            "partial_r": result["r"],
            "p_val": result["p"],
            "ci_low": result["ci_low"],
            "ci_high": result["ci_high"],
            "slope": result["slope"],
            "intercept": result["intercept"]
        }
        rows.append(row)
        save_df_csv(pd.DataFrame([row]), result["out_csv"])

    df_all = pd.DataFrame(rows)
    save_df_csv(df_all, os.path.join(OUT_SCATTER_TABLE, "scatter_groups_summary.csv"))

    readme = [
        "Partial-regression scatter plot notes",
        "=" * 45,
        "These scatter plots were reconstructed from the original monthly raster time series.",
        "1. Pixels were first grouped using the pixel-wise r map, p map, and snow-frequency map.",
        "2. For each group, the original monthly time series were reloaded.",
        "3. Monthly climatology removal was applied consistently with the main analysis; AR(1)-based effective sample size correction was used for pixel-wise significance testing.",
        "4. SDE and surface soil temperature were separately residualized against all covariates.",
        "5. r, p, n, and confidence intervals were recalculated using all valid residual samples in each group.",
        f"6. The displayed points were sampled for readability, with a maximum of {MAX_PLOT_POINTS_PER_GROUP} points per group."
    ]
    save_txt("\n".join(readme), os.path.join(OUT_SCATTER_TABLE, "scatter_readme.txt"))

    log_msg("Scatter figures and scatter summary tables have been saved.")


# =============================================================================
# 14. Domain-mean partial-regression scatter
# =============================================================================
def save_domain_mean_scatter(out_tif, residual_x, residual_y, r, p, ci_low, ci_high, slope, intercept, n):
    fig, ax = plt.subplots(figsize=(8.2, 6.2), facecolor="white")

    if len(residual_x) > 0 and len(residual_y) > 0:
        ax.scatter(
            residual_x,
            residual_y,
            s=28,
            facecolor=SCATTER_FACE,
            edgecolor=SCATTER_EDGE,
            linewidth=0.35,
            alpha=0.75
        )

        if np.isfinite(slope) and np.isfinite(intercept):
            xs = np.linspace(np.nanmin(residual_x), np.nanmax(residual_x), 200)
            ys = slope * xs + intercept
            ax.plot(xs, ys, color=REG_LINE, linewidth=2.6)

    text = (
        f"r = {r:.3f}\n"
        f"p = {format_p_value(p)}\n"
        f"95% CI = [{ci_low:.3f}, {ci_high:.3f}]\n"
        f"n = {int(n):,}"
    )

    ax.text(
        0.05,
        0.95,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12
    )

    ax.set_xlabel(f"Residualized {X_VAR} ({residual_axis_suffix()})", fontsize=13)
    ax.set_ylabel(f"Residualized {TARGET_VAR} ({residual_axis_suffix()})", fontsize=13)
    ax.set_title(f"Domain-mean partial regression: {TARGET_VAR} vs {X_VAR}", fontsize=14)

    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(1.0)

    ax.tick_params(axis="both", which="major", direction="out", length=4, width=1.0, labelsize=11)
    ax.grid(False)

    save_figure_tif(fig, out_tif, dpi=FIG_DPI)
    plt.close(fig)


def generate_domain_mean_partial_regression(total_sums, total_counts, needed_vars):
    mean_series = {}
    for var in needed_vars:
        mean_series[var] = np.divide(
            total_sums[var],
            total_counts[var],
            out=np.full(len(TIME_STEPS), np.nan, dtype=np.float64),
            where=total_counts[var] > 0
        )

    processed_series = {}
    phi_values = {}

    for var in needed_vars:
        series = mean_series[var].astype(np.float64)

        if REMOVE_MONTHLY_CLIM:
            series = remove_monthly_climatology_1d(series)

        if USE_AR1_PREWHITEN:
            series, phi = ar1_prewhiten_1d(series)
        else:
            phi = np.nan

        processed_series[var] = series
        phi_values[var] = phi

    x = processed_series[X_VAR]
    y = processed_series[TARGET_VAR]
    z = np.column_stack([processed_series[covar] for covar in COVARS])

    r, p, ci_low, ci_high, residual_x, residual_y, slope, intercept, n_used = partial_corr_1d(x, y, z)

    out_tif = os.path.join(OUT_DOMAIN_SCATTER, f"{TARGET_VAR}_domain_mean_partial_regression_scatter.tif")
    save_domain_mean_scatter(out_tif, residual_x, residual_y, r, p, ci_low, ci_high, slope, intercept, n_used)

    dates_used = DATES[1:] if USE_AR1_PREWHITEN else DATES
    out_series_csv = os.path.join(OUT_DOMAIN_SCATTER, f"{TARGET_VAR}_domain_mean_partial_regression_series.csv")
    df_series = pd.DataFrame({
        "date": dates_used[:len(residual_x)],
        f"residual_{X_VAR}": residual_x,
        f"residual_{TARGET_VAR}": residual_y
    })
    save_df_csv(df_series, out_series_csv)

    out_summary_csv = os.path.join(OUT_DOMAIN_SCATTER, f"{TARGET_VAR}_domain_mean_partial_corr_summary.csv")
    df_summary = pd.DataFrame([{
        "variable": TARGET_VAR,
        "partial_r": r,
        "p_val": p,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "n_samples": n_used,
        "n_times": len(TIME_STEPS),
        "x_var": X_VAR,
        "covariates": ",".join(COVARS),
        "remove_monthly_clim": REMOVE_MONTHLY_CLIM,
        "use_ar1_prewhiten": USE_AR1_PREWHITEN,
        "snowc_mask_rule": SNOWC_MASK_RULE,
        "snowc_threshold": SNOWC_THRESHOLD,
        "sig_alpha": ALPHA,
        "slope_residual_regression": slope,
        "intercept_residual_regression": intercept
    }])
    save_df_csv(df_summary, out_summary_csv)

    out_phi_csv = os.path.join(OUT_DOMAIN_SCATTER, "domain_mean_ar1_phi_values.csv")
    save_df_csv(pd.DataFrame([phi_values]), out_phi_csv)

    return mean_series, processed_series


# =============================================================================
# 15. Main workflow
# =============================================================================
def main():
    log_msg("========== Pixel-wise partial correlation analysis started ==========")
    log_msg(f"Time range: {YEARS[0]}-{YEARS[-1]}, {len(TIME_STEPS)} monthly steps")
    log_msg(f"Response variable: {TARGET_VAR} ({TARGET_DESC})")
    log_msg(f"Predictor variable: {X_VAR} ({X_DESC})")
    log_msg(f"Covariates: {', '.join(COVARS)}")
    log_msg(f"Remove monthly climatology: {REMOVE_MONTHLY_CLIM}")
    log_msg("AR(1) prewhitening: disabled")
    log_msg(f"Effective sample size correction: {USE_EFFECTIVE_SAMPLE_SIZE}")
    log_msg(f"SNOWC mask rule: {SNOWC_MASK_RULE}, threshold = {SNOWC_THRESHOLD}, minimum snow frequency = {SNOWC_MIN_FRACTION}")
    log_msg(f"N_JOBS = {N_JOBS}")
    log_msg(f"BLOCK_ROWS = {BLOCK_ROWS}")
    log_msg(f"PIXEL_CHUNK = {PIXEL_CHUNK}")
    log_msg(f"Output root: {OUT_ROOT}")
    log_msg("Publication map: no grid lines or longitude/latitude ticks")
    log_msg("Correlation colors: negative RGB (32, 56, 136); positive RGB (225, 156, 102)")
    log_msg("Correlation legend font: Arial")
    log_msg("")

    template = get_template_info()
    init_worker_template(template["crs"], template["transform"], template["width"], template["height"])
    check_inputs(template)

    width = template["width"]
    height = template["height"]
    profile = template["profile"]
    ranges = block_ranges(height, BLOCK_ROWS)

    log_msg(f"Raster size: {width} x {height}")
    log_msg(f"Number of row blocks: {len(ranges)}")

    r_map = np.full((height, width), np.nan, dtype=np.float32)
    p_map = np.full((height, width), np.nan, dtype=np.float32)
    n_map = np.full((height, width), np.nan, dtype=np.float32)
    snow_freq_map = np.full((height, width), np.nan, dtype=np.float32)
    snow_mask_map = np.zeros((height, width), dtype=np.uint8)

    needed_vars = list(dict.fromkeys([TARGET_VAR, X_VAR] + COVARS))
    total_sums = {var: np.zeros(len(TIME_STEPS), dtype=np.float64) for var in needed_vars}
    total_counts = {var: np.zeros(len(TIME_STEPS), dtype=np.int64) for var in needed_vars}

    tasks = [(row_start, row_end, width) for row_start, row_end in ranges]

    with ProcessPoolExecutor(
        max_workers=N_JOBS,
        initializer=init_worker_template,
        initargs=(template["crs"], template["transform"], template["width"], template["height"])
    ) as executor:
        futures = [executor.submit(process_one_block, task) for task in tasks]

        for future in tqdm(as_completed(futures), total=len(futures), desc=TARGET_VAR, ncols=100):
            (
                row_start, row_end,
                r_block, p_block, n_block, snow_freq_block, mask_block,
                ts_sums, ts_counts
            ) = future.result()

            r_map[row_start:row_end, :] = r_block
            p_map[row_start:row_end, :] = p_block
            n_map[row_start:row_end, :] = n_block
            snow_freq_map[row_start:row_end, :] = snow_freq_block
            snow_mask_map[row_start:row_end, :] = mask_block

            for var in needed_vars:
                total_sums[var] += ts_sums[var]
                total_counts[var] += ts_counts[var]

    sig_mask = (np.isfinite(p_map) & (p_map < ALPHA)).astype(np.uint8)

    # -------------------------------------------------------------------------
    # 1) Save full-domain raster outputs
    # -------------------------------------------------------------------------
    out_r_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_r.tif")
    out_p_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_p.tif")
    out_n_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_effective_n.tif")
    out_sig_tif = os.path.join(OUT_RASTER, f"{TARGET_VAR}_partial_corr_sigmask.tif")
    out_mask_tif = os.path.join(OUT_RASTER, "snowc_domain_mask.tif")
    out_snow_freq_tif = os.path.join(OUT_RASTER, "snowc_frequency.tif")

    save_raster(out_r_tif, r_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_p_tif, p_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_n_tif, n_map, profile, nodata=NODATA_OUT, dtype="float32")
    save_raster(out_sig_tif, sig_mask, profile, nodata=255, dtype="uint8")
    save_raster(out_mask_tif, snow_mask_map, profile, nodata=255, dtype="uint8")
    save_raster(out_snow_freq_tif, snow_freq_map, profile, nodata=NODATA_OUT, dtype="float32")

    # -------------------------------------------------------------------------
    # 2) Save quicklook figures
    # -------------------------------------------------------------------------
    out_r_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_r.png")
    out_r_sig_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_r_sig_points.png")
    out_p_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_p.png")
    out_n_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_effective_n.png")
    out_rpn_png = os.path.join(OUT_FIG, f"{TARGET_VAR}_partial_corr_rpn_triptych.png")

    save_corr_png(out_r_png, r_map, f"Pixel-wise partial correlation: {TARGET_VAR} vs {X_VAR}")
    save_corr_png_with_sig_points(
        out_r_sig_png,
        r_map,
        p_map,
        f"Pixel-wise partial correlation: {TARGET_VAR} vs {X_VAR} (significant points, p < {ALPHA})"
    )
    save_p_png(out_p_png, p_map, f"Pixel-wise p-value: {TARGET_VAR} vs {X_VAR}")
    save_n_png(out_n_png, n_map, f"Effective sample size: {TARGET_VAR} vs {X_VAR}")
    save_rpn_triptych(
        out_rpn_png,
        r_map,
        p_map,
        n_map,
        "Overall partial-correlation result"
    )

    log_msg("Full-domain GeoTIFFs and quicklook PNG figures have been saved.")

    # -------------------------------------------------------------------------
    # 3) Save publication-style world map
    # -------------------------------------------------------------------------
    profile_for_pub = profile.copy()
    profile_for_pub.update({
        "width": width,
        "height": height,
        "transform": profile["transform"],
        "crs": profile["crs"]
    })

    try:
        out_r_corrected, out_pub_map = save_publication_world_map_tif(r_map, profile_for_pub)
        log_msg(f"Publication-style world map saved: {out_pub_map}")
    except Exception as error:
        error_txt = os.path.join(OUT_TEXT, "publication_world_map_error.txt")
        save_txt(str(error), error_txt)
        log_msg(f"Publication-style world map failed. Error information saved to: {error_txt}")

    # -------------------------------------------------------------------------
    # 4) Save overall and stratified summaries
    # -------------------------------------------------------------------------
    overall_mask = np.isfinite(r_map) & np.isfinite(p_map)
    overall_stats = summarize_pixel_distribution(
        r_map[overall_mask],
        p_map[overall_mask],
        n_map[overall_mask],
        "Overall snow domain",
        alpha=ALPHA
    )
    overall_stats["region_name"] = "Overall snow domain"

    df_overall = pd.DataFrame([overall_stats])
    out_overall_csv = os.path.join(OUT_TABLE, "overall_valid_pixels_summary.csv")
    save_df_csv(df_overall, out_overall_csv)

    stratified_results = stratified_analysis_by_snow_frequency(
        snow_freq_map,
        r_map,
        p_map,
        n_map,
        width,
        height,
        profile
    )
    df_stratified = pd.DataFrame(list(stratified_results.values()))
    out_stratified_csv = os.path.join(OUT_TABLE, "stratified_valid_pixels_summary.csv")
    save_df_csv(df_stratified, out_stratified_csv)

    df_integrated = pd.concat(
        [df_overall.assign(region="overall"), df_stratified],
        ignore_index=True,
        sort=False
    )
    out_integrated_csv = os.path.join(OUT_TABLE, "integrated_partial_corr_summary.csv")
    save_df_csv(df_integrated, out_integrated_csv)

    out_conclusion_txt = os.path.join(OUT_TEXT, "integrated_partial_corr_conclusion.txt")
    save_txt(make_integrated_conclusion(overall_stats, df_stratified), out_conclusion_txt)

    log_msg("Overall and stratified summary tables have been saved.")

    # -------------------------------------------------------------------------
    # 5) Save publication-style scatter plots
    # -------------------------------------------------------------------------
    generate_scatter_outputs(template, r_map, p_map, snow_freq_map)

    # -------------------------------------------------------------------------
    # 6) Save domain-mean partial-regression scatter and series
    # -------------------------------------------------------------------------
    mean_series, processed_series = generate_domain_mean_partial_regression(
        total_sums,
        total_counts,
        needed_vars
    )
    log_msg("Domain-mean partial-regression results have been saved.")

    # -------------------------------------------------------------------------
    # 7) Diagnostics
    # -------------------------------------------------------------------------
    try:
        z_for_vif = np.column_stack([processed_series[covar] for covar in COVARS])
        vif_df = check_multicollinearity(z_for_vif, COVARS)
        out_vif_csv = os.path.join(OUT_TABLE, "vif_diagnosis.csv")
        save_df_csv(vif_df, out_vif_csv)

        moran_i, moran_p = spatial_autocorr_test(r_map)
        diag_df = pd.DataFrame([{
            "moran_i": moran_i if moran_i is not None else np.nan,
            "moran_p": moran_p if moran_p is not None else np.nan
        }])
        out_diag_csv = os.path.join(OUT_TABLE, "spatial_diagnosis.csv")
        save_df_csv(diag_df, out_diag_csv)

        log_msg("VIF and spatial-autocorrelation diagnostics have been saved.")

    except Exception as error:
        error_txt = os.path.join(OUT_TEXT, "diagnostic_error.txt")
        save_txt(str(error), error_txt)
        log_msg(f"Diagnostics failed, but main outputs were saved. Error information saved to: {error_txt}")

    # -------------------------------------------------------------------------
    # 8) Final log
    # -------------------------------------------------------------------------
    log_msg("=" * 70)
    log_msg(f"Output root: {OUT_ROOT}")
    log_msg(f"Raster outputs: {OUT_RASTER}")
    log_msg(f"Quicklook figures: {OUT_FIG}")
    log_msg(f"Publication-style world map: {OUT_PUB_MAP}")
    log_msg(f"Summary tables: {OUT_TABLE}")
    log_msg(f"Scatter figures: {OUT_SCATTER_FIG}")
    log_msg(f"Scatter summary tables: {OUT_SCATTER_TABLE}")
    log_msg(f"Domain-mean scatter outputs: {OUT_DOMAIN_SCATTER}")
    log_msg("Key publication outputs:")
    log_msg(f"1) 07_publication_maps/{TARGET_VAR}_partial_corr_world_map_pub.tif")
    log_msg("2) 05_scatter_figures/scatter_overall_all_valid_pixels.tif")
    log_msg("3) 05_scatter_figures/scatter_significant_positive_pixels.tif")
    log_msg("4) 05_scatter_figures/scatter_significant_negative_pixels.tif")
    log_msg("5) 05_scatter_figures/scatter_non_significant_pixels.tif")
    log_msg("6) 05_scatter_figures/scatter_permanent_snow_zone.tif")
    log_msg("7) 05_scatter_figures/scatter_seasonal_snow_zone.tif")
    log_msg("8) 05_scatter_figures/scatter_ephemeral_snow_zone.tif")
    log_msg("9) 06_scatter_tables/scatter_groups_summary.csv")
    log_msg(f"10) 08_domain_mean_scatter/{TARGET_VAR}_domain_mean_partial_regression_scatter.tif")
    log_msg("=" * 70)


if __name__ == "__main__":
    freeze_support()
    main()