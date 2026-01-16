#!/usr/bin/env python3
"""
TSR (Terrain Solar Radiation) 效果对比分析脚本

比较 TSR=OFF（基准）与 TSR=ON 的模拟结果差异，
验证地形辐射修正对水文变量的影响。
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# 添加 shud_reader 路径
sys.path.insert(0, str(Path(__file__).parent.parent / "validation" / "tsr" / "py"))
from shud_reader import DatReader, read_dat


def load_dat(path: Path) -> pd.DataFrame:
    """读取 .dat 文件为 DataFrame"""
    return read_dat(path, time_mode="index")


def compute_stats(df: pd.DataFrame) -> dict:
    """计算统计量"""
    values = df.values.flatten()
    valid = values[np.isfinite(values)]
    return {
        "mean": np.mean(valid),
        "std": np.std(valid),
        "min": np.min(valid),
        "max": np.max(valid),
        "median": np.median(valid),
        "q25": np.percentile(valid, 25),
        "q75": np.percentile(valid, 75),
    }


def compare_variable(base_path: Path, tsr_path: Path, var_name: str) -> dict:
    """比较单个变量的差异"""
    df_base = load_dat(base_path)
    df_tsr = load_dat(tsr_path)

    # 确保维度一致
    assert df_base.shape == df_tsr.shape, f"Shape mismatch: {df_base.shape} vs {df_tsr.shape}"

    # 计算差异
    diff = df_tsr.values - df_base.values
    rel_diff = np.where(np.abs(df_base.values) > 1e-10,
                        diff / df_base.values * 100,
                        np.nan)

    valid_diff = diff[np.isfinite(diff)]
    valid_rel = rel_diff[np.isfinite(rel_diff)]

    return {
        "variable": var_name,
        "base_stats": compute_stats(df_base),
        "tsr_stats": compute_stats(df_tsr),
        "diff_mean": np.mean(valid_diff),
        "diff_std": np.std(valid_diff),
        "diff_abs_max": np.max(np.abs(valid_diff)),
        "rel_diff_mean": np.nanmean(valid_rel),
        "rel_diff_abs_max": np.nanmax(np.abs(valid_rel)),
        "n_changed": np.sum(np.abs(valid_diff) > 1e-10),
        "n_total": len(valid_diff),
    }


def analyze_tsr_diagnostics(tsr_dir: Path) -> dict:
    """分析 TSR 诊断输出"""
    results = {}

    # TSR factor
    rn_factor_path = tsr_dir / "ccw.rn_factor.dat"
    if rn_factor_path.exists():
        df = load_dat(rn_factor_path)
        values = df.values.flatten()
        valid = values[np.isfinite(values)]
        results["rn_factor"] = {
            "mean": np.mean(valid),
            "std": np.std(valid),
            "min": np.min(valid),
            "max": np.max(valid),
            "pct_gt_1": np.sum(valid > 1.0) / len(valid) * 100,
            "pct_lt_1": np.sum(valid < 1.0) / len(valid) * 100,
            "pct_eq_1": np.sum(np.abs(valid - 1.0) < 0.01) / len(valid) * 100,
        }

    # 水平面辐射 vs 地形修正辐射
    rn_h_path = tsr_dir / "ccw.rn_h.dat"
    rn_t_path = tsr_dir / "ccw.rn_t.dat"
    if rn_h_path.exists() and rn_t_path.exists():
        df_h = load_dat(rn_h_path)
        df_t = load_dat(rn_t_path)

        h_vals = df_h.values.flatten()
        t_vals = df_t.values.flatten()

        # 只在有辐射时计算比值
        mask = (h_vals > 1.0) & np.isfinite(h_vals) & np.isfinite(t_vals)
        ratio = t_vals[mask] / h_vals[mask]

        results["radiation"] = {
            "rn_h_mean": np.mean(h_vals[np.isfinite(h_vals)]),
            "rn_t_mean": np.mean(t_vals[np.isfinite(t_vals)]),
            "ratio_mean": np.mean(ratio),
            "ratio_std": np.std(ratio),
            "ratio_min": np.min(ratio),
            "ratio_max": np.max(ratio),
        }

    return results


def main():
    # 目录配置
    base_dir = Path("output/ccw_base.out")
    tsr_dir = Path("output/ccw_tsr")

    print("=" * 70)
    print("TSR (Terrain Solar Radiation) 效果对比分析")
    print("=" * 70)
    print(f"\n基准目录 (TSR=OFF): {base_dir}")
    print(f"TSR目录 (TSR=ON):   {tsr_dir}")

    # 1. 分析 TSR 诊断输出
    print("\n" + "=" * 70)
    print("1. TSR 诊断变量分析")
    print("=" * 70)

    tsr_diag = analyze_tsr_diagnostics(tsr_dir)

    if "rn_factor" in tsr_diag:
        f = tsr_diag["rn_factor"]
        print(f"\n[TSR Factor (地形辐射因子)]")
        print(f"  均值: {f['mean']:.4f}")
        print(f"  标准差: {f['std']:.4f}")
        print(f"  范围: [{f['min']:.4f}, {f['max']:.4f}]")
        print(f"  因子>1 (增强): {f['pct_gt_1']:.1f}%")
        print(f"  因子<1 (削弱): {f['pct_lt_1']:.1f}%")
        print(f"  因子≈1 (无变化): {f['pct_eq_1']:.1f}%")

    if "radiation" in tsr_diag:
        r = tsr_diag["radiation"]
        print(f"\n[辐射对比]")
        print(f"  水平面辐射均值 (rn_h): {r['rn_h_mean']:.2f} W/m²")
        print(f"  地形修正辐射均值 (rn_t): {r['rn_t_mean']:.2f} W/m²")
        print(f"  修正比值: {r['ratio_mean']:.4f} ± {r['ratio_std']:.4f}")
        print(f"  比值范围: [{r['ratio_min']:.4f}, {r['ratio_max']:.4f}]")

    # 2. 关键水文变量对比
    print("\n" + "=" * 70)
    print("2. 关键水文变量对比 (TSR=ON vs TSR=OFF)")
    print("=" * 70)

    # 需要比较的变量
    variables = [
        ("eleveta", "实际蒸散发 (ET)"),
        ("elevetp", "潜在蒸散发 (PET)"),
        ("elevettr", "蒸腾 (Transpiration)"),
        ("elevetev", "蒸发 (Evaporation)"),
        ("elevetic", "冠层截留蒸发 (Interception)"),
        ("eleygw", "地下水位 (GW)"),
        ("eleyunsat", "非饱和带水量 (Unsat)"),
        ("eleysurf", "地表水深 (Surf)"),
        ("rivqdown", "河道流量 (Q)"),
    ]

    comparison_results = []

    for var_code, var_name in variables:
        base_path = base_dir / f"ccw.{var_code}.dat"
        tsr_path = tsr_dir / f"ccw.{var_code}.dat"

        if not base_path.exists() or not tsr_path.exists():
            print(f"\n[{var_name}] - 文件缺失，跳过")
            continue

        result = compare_variable(base_path, tsr_path, var_name)
        comparison_results.append(result)

        print(f"\n[{var_name}] ({var_code})")
        print(f"  基准均值: {result['base_stats']['mean']:.6e}")
        print(f"  TSR均值:  {result['tsr_stats']['mean']:.6e}")
        print(f"  差异均值: {result['diff_mean']:.6e}")
        print(f"  差异最大绝对值: {result['diff_abs_max']:.6e}")
        print(f"  相对差异均值: {result['rel_diff_mean']:.2f}%")
        print(f"  变化单元数: {result['n_changed']:,} / {result['n_total']:,}")

    # 3. 汇总表格
    print("\n" + "=" * 70)
    print("3. 差异汇总")
    print("=" * 70)

    print(f"\n{'变量':<25} {'基准均值':>12} {'TSR均值':>12} {'相对差异%':>10} {'变化率%':>10}")
    print("-" * 70)

    for r in comparison_results:
        change_pct = r['n_changed'] / r['n_total'] * 100 if r['n_total'] > 0 else 0
        print(f"{r['variable']:<25} {r['base_stats']['mean']:>12.4e} {r['tsr_stats']['mean']:>12.4e} "
              f"{r['rel_diff_mean']:>10.2f} {change_pct:>10.1f}")

    # 4. 结论
    print("\n" + "=" * 70)
    print("4. 分析结论")
    print("=" * 70)

    # 检查是否有显著差异
    significant_vars = [r for r in comparison_results if abs(r['rel_diff_mean']) > 0.1]

    if significant_vars:
        print("\n[有显著影响的变量] (相对差异 > 0.1%)")
        for r in significant_vars:
            print(f"  - {r['variable']}: {r['rel_diff_mean']:.2f}%")
    else:
        print("\n[注意] 所有变量的相对差异均 < 0.1%")
        print("  可能原因：")
        print("  1. 研究区地形较平坦，TSR修正效果有限")
        print("  2. 输出为日均值，瞬时差异被平滑")
        print("  3. 需检查 TSR 是否正确启用")

    # 保存详细结果
    output_file = Path("post_analysis/tsr_comparison_results.csv")
    df_results = pd.DataFrame(comparison_results)
    df_results.to_csv(output_file, index=False)
    print(f"\n详细结果已保存至: {output_file}")


if __name__ == "__main__":
    main()
