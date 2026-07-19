# -*- coding: utf-8 -*-
"""
Plotting only script: 4 ecosystems response + separate ratio figures + overall ratio figure
------------------------------------------------------------------------------------------
只读取已有结果，不重新计算模型。

输入：
1) sampled_panel_filtered_with_landuse.parquet / csv.gz
2) model_coefficients_by_ecosystem.csv

说明：
- 不再读取 binned_observed_response_by_ecosystem.csv。
- 在原代码生成的 all_bin_df 中新增 overall 分组。
- Figure 2b 直接从 all_bin_df 中筛选 group_name == "overall" 后计算。

输出目录：
RESULT_DIR/6.9结果_全Arial_overall左对齐修正

输出：
1) four_ecosystems_soil_temperature_response.tif
   - Forest, Wetland, Grassland, Tundra 横向排列
   - 共用一个 y 轴
   - 只保留 Forest 的 y 轴刻度
   - y 轴范围：-10 ~ 10
   - y 轴刻度：-10、-5、0、5、10
   - 横坐标：Snow depth (cm)
   - 纵坐标：Soil temperature (℃)
   - 不嵌入 ratio 小图
   - 显著性星号上下错开，避免重叠
   - 四个子图之间间距缩短，且三个间距一致

2) Figure 2b：global_overall_asymmetry_ratio.tif
   - 直接从原代码生成的 all_bin_df 中提取 overall 结果
   - 计算：|Effect_Decreased| / |Effect_Increased|
   - 文字内容、图例内容与 4 张单独 ratio 图保持一致
   - 不显示 rho_snow 蓝色水平线，但在图例中保留 rho_snow 数值文字
   - 版式改为竖向长方形
   - x 轴：0、20、40、60
   - y 轴：仅正值
   - 字体：Arial

3) 4 张单独的不对称性 ratio 小图：
   - forest_ratio_response.tif
   - wetland_ratio_response.tif
   - grassland_ratio_response.tif
   - tundra_ratio_response.tif
   - 横坐标均为：Snow depth (cm)
   - 比值曲线保持灰黑色
   - 不显示 rho_snow 蓝色水平线，但在图例中保留 rho_snow 数值文字
   - ratio 公式统一为：|Effect_Decreased| / |Effect_Increased|

4) 对应统计表：
   - binned_response_4ecosystems.csv
   - overall_binned_response_from_all_bin_df.csv
   - ratio_response_4ecosystems.csv
   - rho_snow_4ecosystems.csv
   - overall_ratio_from_all_bin_df.csv
"""

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# 0) 路径与参数
# =========================================================
RESULT_DIR = r"I:\ERA_TIFF\最终版积雪厚度-土壤温度-生态系统ecosystem_stl1_stable_snow_lag1_same_amplitude_1_100cm_final_figs"

# 最终输出路径
PLOT_OUT_DIR = os.path.join(RESULT_DIR, "6.9结果_全Arial_overall左对齐修正")
os.makedirs(PLOT_OUT_DIR, exist_ok=True)

PANEL_PARQUET = os.path.join(RESULT_DIR, "sampled_panel_filtered_with_landuse.parquet")
PANEL_CSVGZ = os.path.join(RESULT_DIR, "sampled_panel_filtered_with_landuse.csv.gz")
COEF_CSV = os.path.join(RESULT_DIR, "model_coefficients_by_ecosystem.csv")

# 积雪深度单位
# 如果原始 delta_snow_lag1 是 m，保留 "m"
# 如果原始 delta_snow_lag1 已经是 cm，改成 "cm"
SDE_UNIT = "m"

# 同幅度分箱，只画到 60 cm
AMP_CM = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60]

# 散点抽样上限
MAX_SCATTER_PER_SIGN = 25000
RANDOM_SEED = 42

# 输出分辨率
DPI = 600

# 4 个生态系统主图参数
MAIN_FIG_W = 16.0
MAIN_FIG_H = 4.8
MAIN_YMIN = -20.0
MAIN_YMAX = 10.0

PLOT_XMIN = 0.0
PLOT_XMAX = 62.0

# 单独 ratio 小图参数
RATIO_FIG_W = 6.0
RATIO_FIG_H = 4.5

# overall ratio 图参数（Figure 2b）
OVERALL_RATIO_FIG_W = 4.8
OVERALL_RATIO_FIG_H = 7.0

# 字体
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.unicode_minus"] = False
# 强制所有数学公式（包括下标和希腊字母）使用 Arial 字体
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["mathtext.sf"] = "Arial"
plt.rcParams["mathtext.tt"] = "Arial"
plt.rcParams["mathtext.cal"] = "Arial"

# 主图字号
TITLE_SIZE = 18
LABEL_SIZE = 18
TICK_SIZE = 13
STAR_SIZE = 7

# ratio 图字号
RATIO_TITLE_SIZE = 18
RATIO_LABEL_SIZE = 16
RATIO_TICK_SIZE = 13
RATIO_LEGEND_SIZE = 13

SPINE_WIDTH = 0.8

# 颜色
COLOR_INC = (255 / 255, 192 / 255, 127 / 255)      # Increased snowpack thickness
COLOR_DEC = (143 / 255, 196 / 255, 222 / 255)      # Decreased snowpack thickness

# ratio 图中“比值曲线”保持灰黑色
COLOR_RATIO_LINE = (70 / 255, 70 / 255, 70 / 255)

# 星号颜色
COLOR_STAR = (70 / 255, 70 / 255, 70 / 255)

# ratio 图中 rho_snow 线颜色：RGB(70, 131, 180)
COLOR_RHO = (70 / 255, 131 / 255, 180 / 255)

# 是否显示显著性星号
SHOW_STARS = True

# ratio 小图是否显示图例
SHOW_RATIO_LEGEND = True

# =========================================================
# 1) 生态系统顺序
# =========================================================
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
# 2) 基础函数
# =========================================================
def load_table_any(parquet_path=None, csv_gz_path=None):
    if parquet_path is not None and os.path.exists(parquet_path):
        return pd.read_parquet(parquet_path)

    if csv_gz_path is not None and os.path.exists(csv_gz_path):
        return pd.read_csv(csv_gz_path, compression="gzip")

    raise FileNotFoundError(
        f"Cannot find file:\n{parquet_path}\nor\n{csv_gz_path}"
    )


def load_panel_and_coef():
    print(f"[INFO] Panel parquet: {PANEL_PARQUET}")
    print(f"[INFO] Panel csv.gz : {PANEL_CSVGZ}")
    print(f"[INFO] Coef csv     : {COEF_CSV}")

    panel = load_table_any(PANEL_PARQUET, PANEL_CSVGZ)

    if not os.path.exists(COEF_CSV):
        raise FileNotFoundError(f"Cannot find file:\n{COEF_CSV}")

    coef = pd.read_csv(COEF_CSV, encoding="utf-8-sig")

    return panel, coef


def ensure_needed_columns(panel):
    need_cols = [
        "pixel_id",
        "year",
        "month",
        "delta_snow_lag1",
        "delta_ts",
        "ecosystem"
    ]

    miss = [c for c in need_cols if c not in panel.columns]

    if miss:
        raise RuntimeError(f"panel 缺少必要字段：{miss}")

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

    try:
        from scipy.stats import t as tdist
        p = 2 * (1 - tdist.cdf(abs(t_val), df=n - 1))
        return float(p)
    except Exception:
        try:
            from math import erf, sqrt
            p = 2 * (1 - 0.5 * (1 + erf(abs(t_val) / sqrt(2))))
            return float(p)
        except Exception:
            return np.nan


# =========================================================
# 3) 数据准备与统计
# =========================================================
def prepare_plot_panel(panel):
    panel = panel.copy()
    panel = ensure_needed_columns(panel)

    panel["amp_cm"] = snow_to_cm(np.abs(panel["delta_snow_lag1"]))

    panel["sign_group"] = np.where(
        panel["delta_snow_lag1"] > 0,
        "increase",
        np.where(panel["delta_snow_lag1"] < 0, "decrease", "zero")
    )

    panel = panel.replace([np.inf, -np.inf], np.nan)

    panel = panel.dropna(
        subset=[
            "delta_snow_lag1",
            "delta_ts",
            "amp_cm",
            "ecosystem"
        ]
    )

    panel = panel[panel["sign_group"].isin(["increase", "decrease"])].copy()

    panel["group_name"] = panel["ecosystem"].astype(str)

    # 只保留森林、湿地、草地、苔原
    panel = panel[panel["group_name"].isin(ECO_GROUPS)].copy()

    return panel


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
    inc = bin_df[bin_df["sign_group"] == "increase"][["amp_bin_cm", "mean", "ci95", "n"]].copy()
    dec = bin_df[bin_df["sign_group"] == "decrease"][["amp_bin_cm", "mean", "ci95", "n"]].copy()

    inc = inc.rename(columns={"mean": "mean_inc", "ci95": "ci95_inc", "n": "n_inc"})
    dec = dec.rename(columns={"mean": "mean_dec", "ci95": "ci95_dec", "n": "n_dec"})

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

    need_cols = ["group_name", "variable", field]
    for c in need_cols:
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
    beta_pos = get_coef_value(coef_df, group_name, "snow_pos_lag1", "coef")
    beta_neg = get_coef_value(coef_df, group_name, "snow_neg_lag1", "coef")
    p_pos = get_coef_value(coef_df, group_name, "snow_pos_lag1", "pvalue")
    p_neg = get_coef_value(coef_df, group_name, "snow_neg_lag1", "pvalue")

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


def calc_overall_rho_snow_model(coef_df):
    """尽量从已有系数表中提取 overall 的 rho_snow，不重新计算模型。"""
    if coef_df is None or len(coef_df) == 0 or "group_name" not in coef_df.columns:
        return {
            "group_name": "overall",
            "beta_pos": np.nan,
            "beta_neg": np.nan,
            "p_pos": np.nan,
            "p_neg": np.nan,
            "rho_snow": np.nan
        }

    group_candidates = ["overall", "global_overall", "global", "all", "pooled"]
    group_series = coef_df["group_name"].astype(str).str.strip().str.lower()

    for candidate in group_candidates:
        if group_series.eq(candidate).any():
            matched_name = coef_df.loc[group_series.eq(candidate), "group_name"].iloc[0]
            return calc_rho_snow_model(coef_df, matched_name)

    return {
        "group_name": "overall",
        "beta_pos": np.nan,
        "beta_neg": np.nan,
        "p_pos": np.nan,
        "p_neg": np.nan,
        "rho_snow": np.nan
    }


# =========================================================
# 4) 绘图样式
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
        lbl.set_fontfamily("Arial")

    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily("Arial")


def set_legend_font(leg):
    if leg is None:
        return
    for txt in leg.get_texts():
        txt.set_fontfamily("Arial")


def add_star_text(ax, x, y, star, above=True, tier=0, group_name=None):
    """
    tier is used to vertically stagger multiple significance stars at the same snow-depth bin.
    Wetland uses a larger tier spacing to avoid overlap between the first and second ***.
    """
    if (not SHOW_STARS) or star is None or star == "" or (not np.isfinite(y)):
        return

    y0, y1 = ax.get_ylim()
    yr = y1 - y0

    if group_name == "wetland":
        offset = (0.050 + 0.110 * tier) * yr
    else:
        offset = (0.040 + 0.060 * tier) * yr

    x_offsets = {
        1: 1.5,
        5: -0.35,
        10: 0.35,
        15: -0.35,
        20: 0.35,
        25: -0.25,
        30: 0.25,
        40: -0.2,
        50: 0.2,
        60: 0.0
    }
    x_text = x + x_offsets.get(int(round(x)), 0.0)

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
        fontfamily="Arial",
        zorder=100,
        clip_on=False
    )


def add_stars_without_overlap(ax, bin_df, group_name=None):
    if not SHOW_STARS:
        return

    star_df = bin_df.copy()
    star_df = star_df[star_df["stars"].astype(str) != ""].copy()
    star_df = star_df.dropna(subset=["mean", "amp_bin_cm"])

    if len(star_df) == 0:
        return

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
                tier=tier,
                group_name=group_name
            )


# =========================================================
# 5) 主响应图：4个生态系统横向排列
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
                fontfamily="Arial"
            )
            ax.text(
                0.5,
                0.5,
                "No data",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=16,
                fontfamily="Arial"
            )

            ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
            ax.set_ylim(MAIN_YMIN, MAIN_YMAX)
            ax.set_yticks([10, 0, -10, -20])

            set_axes_style(ax, tick_size=TICK_SIZE, spine_width=SPINE_WIDTH)
            for lbl in ax.get_yticklabels():
                lbl.set_fontfamily("Arial")

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

        # 原始散点
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

        # 分箱均值曲线 + 95% CI
        for sign, color in [
            ("increase", COLOR_INC),
            ("decrease", COLOR_DEC)
        ]:
            dd = bin_df[bin_df["sign_group"] == sign].copy().sort_values("amp_bin_cm")
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
                markeredgecolor=(0.75, 0.75, 0.75),
                markeredgewidth=0.5,
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
            linewidth=1.0,
            color="#777777",
            zorder=2
        )

        ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
        ax.set_ylim(MAIN_YMIN, MAIN_YMAX)
        ax.set_xticks([0, 20, 40, 60])
        ax.set_yticks([10, 0, -10, -20])

        ax.set_title(
            group_label,
            fontsize=TITLE_SIZE,
            fontfamily="Arial",
            pad=8
        )

        set_axes_style(ax, tick_size=TICK_SIZE, spine_width=SPINE_WIDTH)
        for lbl in ax.get_yticklabels():
            lbl.set_fontfamily("Arial")

        # 组合图不显示显著性星号
        # add_stars_without_overlap(ax, bin_df, group_name=group_name)

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
        fontfamily="Arial"
    )

    fig.text(
        0.032,
        0.52,
        "Soil temperature (℃)",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=LABEL_SIZE,
        fontfamily="Arial"
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

    # ---------------------------------------------------------
    # 直接用完整 plot_panel 计算 overall 分箱结果，并加入 all_bin_df。
    # 这样 Figure 2b 可直接从 all_bin_df 中筛选 overall，
    # 不再依赖任何外部 binned_observed_response_by_ecosystem.csv。
    # ---------------------------------------------------------
    if len(plot_panel) > 0:
        overall_bin_df = calc_binned_stats(plot_panel)
        overall_bin_df["group_name"] = "overall"
        overall_bin_df["group_label"] = "Global overall"
        all_bin_rows.append(overall_bin_df)

    if len(all_bin_rows) > 0:
        all_bin_df = pd.concat(all_bin_rows, axis=0, ignore_index=True)
        all_bin_df = all_bin_df.sort_values(
            ["group_name", "sign_group", "amp_bin_cm"]
        ).reset_index(drop=True)
    else:
        all_bin_df = pd.DataFrame()

    return all_bin_df


# =========================================================
# 6) 单独 ratio 小图：4个生态系统分别输出
# =========================================================
def draw_single_ratio_figure(ratio_df, rho_info, group_name, out_path):
    group_label = GROUP_LABELS_EN[group_name]

    fig, ax = plt.subplots(figsize=(RATIO_FIG_W, RATIO_FIG_H))
    rr = ratio_df.copy()

    # 比值结果曲线：保持灰黑色
    ax.plot(
        rr["amp_bin_cm"],
        rr["ratio_obs"],
        color=COLOR_RATIO_LINE,
        linewidth=2.2,
        marker="o",
        markersize=5.2,
        label=r"$|Effect_{Decreased}| / |Effect_{Increased}|$",
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

    # 不绘制 rho_snow 蓝色水平线，也不将其加入左上角图例。
    # 仅生成文字，并在后面单独放到右上角。
    if np.isfinite(rho_info["rho_snow"]):
        rho_text = rf"$\rho_{{snow}}$ = {rho_info['rho_snow']:.2f}"
    else:
        rho_text = r"$\rho_{snow}$ = NA"

    ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xticklabels(["0", "20", "40", "60"])

    ratio_valid = rr["ratio_obs"].values
    ratio_valid = ratio_valid[np.isfinite(ratio_valid)]
    values = ratio_valid.tolist()

    if len(values) > 0:
        # 在曲线最大值上方至少保留一个完整整数刻度，
        # 为左上角图例和右上角 rho_snow 文字预留空间，避免重叠。
        data_max = float(np.nanmax(values))
        ymax = max(1, int(np.ceil(data_max)) + 1)
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(0, ymax + 1, 1))
        ax.set_yticklabels([str(i) for i in range(0, ymax + 1)])
    else:
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["0", "1"])

    ax.set_title(
        group_label,
        fontsize=RATIO_TITLE_SIZE,
        fontfamily="Arial",
        pad=8
    )

    ax.set_xlabel(
        "Snow depth (cm)",
        fontsize=RATIO_LABEL_SIZE,
        fontfamily="Arial",
        labelpad=8
    )

    ax.set_ylabel(
        "Asymmetry ratio",
        fontsize=RATIO_LABEL_SIZE,
        fontfamily="Arial",
        labelpad=8
    )

    set_axes_style(ax, tick_size=RATIO_TICK_SIZE, spine_width=SPINE_WIDTH)

    if SHOW_RATIO_LEGEND:
        leg = ax.legend(
            loc="upper left",
            fontsize=RATIO_LEGEND_SIZE,
            frameon=False,
            handlelength=1.8,
            handletextpad=0.6,
            labelspacing=1.5
        )
        set_legend_font(leg)

    ax.text(
        0.985,
        0.98,
        rho_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=RATIO_LEGEND_SIZE,
        fontfamily="Arial"
    )

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

        df_g = plot_panel[plot_panel["group_name"] == group_name].copy()

        if len(df_g) == 0:
            print(f"[SKIP] {group_label}: no rows for ratio figure")
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
        all_ratio_df = pd.concat(all_ratio_rows, axis=0, ignore_index=True)
    else:
        all_ratio_df = pd.DataFrame()

    if len(rho_rows) > 0:
        rho_df = pd.DataFrame(rho_rows)
    else:
        rho_df = pd.DataFrame()

    return all_ratio_df, rho_df


# =========================================================
# 7) Figure 2b：直接从 all_bin_df 中提取 overall 并绘图
# =========================================================
def extract_overall_binned_from_all_bin_df(all_bin_df):
    """从 draw_four_ecosystems_response() 返回的 all_bin_df 中提取 overall。"""
    required_cols = {
        "group_name",
        "sign_group",
        "amp_bin_cm",
        "mean",
        "ci95",
        "n"
    }

    missing = sorted(required_cols.difference(all_bin_df.columns))
    if missing:
        raise RuntimeError(
            f"all_bin_df 缺少必要字段：{missing}。"
        )

    overall_bin_df = all_bin_df[
        all_bin_df["group_name"].astype(str).str.lower().eq("overall")
    ].copy()

    if len(overall_bin_df) == 0:
        raise RuntimeError(
            "all_bin_df 中未找到 group_name == 'overall' 的结果。"
        )

    overall_bin_df = overall_bin_df.sort_values(
        ["sign_group", "amp_bin_cm"]
    ).reset_index(drop=True)

    return overall_bin_df


def draw_global_overall_ratio_figure(ratio_df, rho_info, out_path):
    fig, ax = plt.subplots(figsize=(OVERALL_RATIO_FIG_W, OVERALL_RATIO_FIG_H))
    rr = ratio_df.copy().sort_values("amp_bin_cm")

    ax.plot(
        rr["amp_bin_cm"],
        rr["ratio_obs"],
        color=COLOR_RATIO_LINE,
        linewidth=2.2,
        marker="o",
        markersize=5.6,
        label=r"$|Effect_{Decreased}| / |Effect_{Increased}|$",
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

    # 不绘制 rho_snow 蓝色水平线，也不将其加入左上角图例。
    # 仅生成文字，并在后面单独放到右上角。
    if np.isfinite(rho_info["rho_snow"]):
        rho_text = rf"$\rho_{{snow}}$ = {rho_info['rho_snow']:.2f}"
    else:
        rho_text = r"$\rho_{snow}$ = NA"

    ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xticklabels(["0", "20", "40", "60"])

    ratio_valid = rr["ratio_obs"].values
    ratio_valid = ratio_valid[np.isfinite(ratio_valid)]
    values = ratio_valid.tolist()

    if len(values) > 0:
        # 在曲线最大值上方至少保留一个完整整数刻度，
        # 为左上角图例和右上角 rho_snow 文字预留空间，避免重叠。
        data_max = float(np.nanmax(values))
        ymax = max(1, int(np.ceil(data_max)) + 1)
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(0, ymax + 1, 1))
        ax.set_yticklabels([str(i) for i in range(0, ymax + 1)])
    else:
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["0", "1"])

    ax.set_title(
        "Global overall",
        fontsize=RATIO_TITLE_SIZE,
        fontfamily="Arial",
        pad=8
    )

    ax.set_xlabel(
        "Snow depth (cm)",
        fontsize=RATIO_LABEL_SIZE,
        fontfamily="Arial",
        labelpad=8
    )

    ax.set_ylabel(
        "Asymmetry ratio",
        fontsize=RATIO_LABEL_SIZE,
        fontfamily="Arial",
        labelpad=8
    )

    set_axes_style(ax, tick_size=RATIO_TICK_SIZE, spine_width=SPINE_WIDTH)

    leg = ax.legend(
        loc="upper left",
        fontsize=RATIO_LEGEND_SIZE,
        frameon=False,
        handlelength=1.8,
        handletextpad=0.6,
        labelspacing=1.5
    )
    set_legend_font(leg)

    ax.text(
        0.045,
        0.82,
        rho_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=RATIO_LEGEND_SIZE,
        fontfamily="Arial"
    )

    fig.subplots_adjust(
        left=0.20,
        right=0.96,
        bottom=0.11,
        top=0.93
    )

    fig.savefig(
        out_path,
        dpi=DPI,
        facecolor="white",
        format="tif"
    )
    plt.close(fig)



# =========================================================
# 8) 主程序
# =========================================================
def main():
    panel, coef_df = load_panel_and_coef()
    plot_panel = prepare_plot_panel(panel)

    print("=" * 90)
    print("开始绘图：overall 直接从 all_bin_df 中提取，不读取外部分箱文件")
    print("=" * 90)

    # 1) 四生态系统主图
    main_out = os.path.join(PLOT_OUT_DIR, "four_ecosystems_soil_temperature_response.tif")
    print(f"[DRAW MAIN] Four ecosystems response -> {main_out}")

    all_bin_df = draw_four_ecosystems_response(
        plot_panel=plot_panel,
        out_path=main_out
    )

    # 2) 四张单独 ratio 图
    all_ratio_df, rho_df = draw_all_ratio_figures(
        plot_panel=plot_panel,
        coef_df=coef_df
    )

    # 3) Figure 2b：直接从 all_bin_df 中提取 overall 并绘图
    overall_bin_df = extract_overall_binned_from_all_bin_df(all_bin_df)
    overall_ratio_df = calc_observed_ratio_curve(overall_bin_df)
    overall_rho_info = calc_overall_rho_snow_model(coef_df)

    overall_out = os.path.join(PLOT_OUT_DIR, "global_overall_asymmetry_ratio.tif")
    print(f"[DRAW OVERALL] Global overall asymmetry ratio -> {overall_out}")

    draw_global_overall_ratio_figure(
        ratio_df=overall_ratio_df,
        rho_info=overall_rho_info,
        out_path=overall_out
    )

    # 4) 保存统计表
    if len(all_bin_df) > 0:
        # 四类生态系统结果：保持原输出文件名不变
        four_ecosystem_bin_df = all_bin_df[
            all_bin_df["group_name"].isin(ECO_GROUPS)
        ].copy()

        four_ecosystem_bin_df.to_csv(
            os.path.join(PLOT_OUT_DIR, "binned_response_4ecosystems.csv"),
            index=False,
            encoding="utf-8-sig"
        )

        # overall 结果明确保存，数据来源就是 all_bin_df
        overall_bin_df.to_csv(
            os.path.join(PLOT_OUT_DIR, "overall_binned_response_from_all_bin_df.csv"),
            index=False,
            encoding="utf-8-sig"
        )

    if len(all_ratio_df) > 0:
        all_ratio_df.to_csv(
            os.path.join(PLOT_OUT_DIR, "ratio_response_4ecosystems.csv"),
            index=False,
            encoding="utf-8-sig"
        )

    if len(rho_df) > 0:
        rho_df.to_csv(
            os.path.join(PLOT_OUT_DIR, "rho_snow_4ecosystems.csv"),
            index=False,
            encoding="utf-8-sig"
        )

    if len(overall_ratio_df) > 0:
        overall_ratio_df.to_csv(
            os.path.join(PLOT_OUT_DIR, "overall_ratio_from_all_bin_df.csv"),
            index=False,
            encoding="utf-8-sig"
        )

    pd.DataFrame([overall_rho_info]).to_csv(
        os.path.join(PLOT_OUT_DIR, "overall_rho_snow_from_coef.csv"),
        index=False,
        encoding="utf-8-sig"
    )

    print("=" * 90)
    print("绘图完成")
    print("=" * 90)
    print(f"输出目录：{PLOT_OUT_DIR}")
    print("输出文件：")
    print(" - four_ecosystems_soil_temperature_response.tif")
    print(" - global_overall_asymmetry_ratio.tif")
    print(" - forest_ratio_response.tif")
    print(" - wetland_ratio_response.tif")
    print(" - grassland_ratio_response.tif")
    print(" - tundra_ratio_response.tif")
    print(" - binned_response_4ecosystems.csv")
    print(" - ratio_response_4ecosystems.csv")
    print(" - rho_snow_4ecosystems.csv")
    print(" - overall_binned_response_from_all_bin_df.csv")
    print(" - overall_ratio_from_all_bin_df.csv")
    print(" - overall_rho_snow_from_coef.csv")


if __name__ == "__main__":
    main()