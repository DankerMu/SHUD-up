# Issue 12 TSR 调试输出：端到端验证报告（`./shud ccw`）

## 结论摘要
- ✅ 新增输出文件已生成：`ccw.rn_h.dat`、`ccw.rn_t.dat`、`ccw.rn_factor.dat`
- ✅ 二进制 header（前 1024 字节）包含 `TSR/radiation/lon/lat` 配置信息
- ✅ 数学一致性通过：`rn_t = rn_h * factor`（最大绝对误差 `8.53e-14`，小时测试样例）
- ⚠️ **TSR=ON 时 `rn_h` 与 `factor` 的昼夜相位明显不一致**：`rn_h>0` 与 `factor>0` 的重叠仅 `2/48` 小时；在 `rn_h` 最大的小时（`441.233 W/m²`）域内 `factor` 全为 `0`，导致 `rn_t` 全域为 `0`。这更像是**时间基准/时区不一致（forcing 时间轴 vs solarPosition 时间轴）**导致的系统性偏移，而非地形遮蔽本身。

## 测试说明
本次验证包含两次运行/两套输出，用于同时满足“文件存在/头信息”和“昼夜合理性”的检查需求：

1) **默认 `ccw` 配置（TSR=OFF，日尺度输出）**
- 命令：`./shud ccw`
- 输出目录：`output/ccw.out`
- 说明：`DT_QE_ET=1440`（日平均），无法直接从输出中判断“夜间=0”，但可验证文件存在、header、以及 `TSR=OFF` 下 `factor≈1` 与 `rn_t=rn_h`。

2) **短窗口小时输出（TSR=ON，用于昼夜检查）**
- 命令：`./shud ccw`
- 运行时临时将 `input/ccw/ccw.cfg.para` 调整为：`END=2`、`DT_QE_ET=60`、`TERRAIN_RADIATION=1`
- 输出目录：`output/ccw.out.tsr_on_hourly_2d`
- 说明：运行后已将 `input/ccw/ccw.cfg.para` 恢复；上述临时配置保存在 `input/ccw/ccw.cfg.para.bak_tsr_test` 便于复现实验。

## 输出文件存在性
- `output/ccw.out/ccw.rn_h.dat`、`output/ccw.out/ccw.rn_t.dat`、`output/ccw.out/ccw.rn_factor.dat`
- `output/ccw.out.tsr_on_hourly_2d/ccw.rn_h.dat`、`output/ccw.out.tsr_on_hourly_2d/ccw.rn_t.dat`、`output/ccw.out.tsr_on_hourly_2d/ccw.rn_factor.dat`

## Header 校验（前 1024 字节）
两套输出的 header 均包含以下字段（示例为 `rn_h.dat`）：
- `# Radiation input mode: SWDOWN`
- `# Terrain radiation (TSR): OFF/ON`
- `# Solar lon/lat mode: FORCING_FIRST`
- `# Solar lon/lat (deg): lon=-122.710000, lat=39.195000`

## 数值合理性与一致性校验

### A. 默认运行（`output/ccw.out`，TSR=OFF，日平均）
- 记录数：`1827`（与 `END=1827` 天一致），变量数：`1147`（NumEle）
- `rn_h`（W/m²，日平均）整体统计：`min=30.120`，`max=384.389`，`mean=220.906`
- `factor`：全时空恒为 `1.0`（符合 TSR=OFF 预期：不做地形修正）
- `rn_t`：与 `rn_h` 完全一致
- `rn_t = rn_h * factor`：最大绝对误差 `0`

### B. 小时验证运行（`output/ccw.out.tsr_on_hourly_2d`，TSR=ON，小时平均）
- 记录数：`48`（2 天 * 24 小时），变量数：`1147`
- `rn_h`（W/m²，小时平均，取第 1 个单元代表）：`>0` 的小时数 `20`，`=0` 的小时数 `28`；`max=441.233`
- `factor`：
  - “全域为 0”的小时数：`30`
  - “存在非 0”的小时数：`18`
  - 全时空统计：`min=0`，`max=5`（命中 `RAD_FACTOR_CAP=5`），`mean≈0.410`
  - “水平面≈1”检验（以 `t=2100 min` 的一小时为代表）：`factor_mean≈0.998`，且 `factor∈[0.95,1.05]` 的单元占比约 `11.77%`
- `rn_t = rn_h * factor`：最大绝对误差 `8.53e-14`（数值一致性通过）
- **异常现象（强烈建议关注）**：
  - `rn_h>0` 且 `factor>0` 的小时仅 `2/48`（`t=960`、`2400`）
  - 在 `rn_h` 达到最大值的小时（`t=1200 min`，`rn_h=441.233 W/m²`），`factor` 在 **1147 个单元上全部为 `0`**，导致 `rn_t` 全域为 `0` —— 这通常不符合物理直觉（至少近平坦面应有 `factor≈1`）。

## 建议/下一步
- 若 `ccw` forcing 时间轴为 UTC（从 `forcing.csv` 首日峰值出现在 `20:00` 这一特征看很像 UTC），而 TSR 的 `solarPosition()` 默认使用经度推算时区（更像“本地时间”假设），则会出现上述相位错配。建议明确 forcing 的时间基准（UTC vs 本地时间），并在 TSR 太阳几何计算处增加可配置的时区/时间基准选项，或在 forcing 预处理阶段完成时区转换。

