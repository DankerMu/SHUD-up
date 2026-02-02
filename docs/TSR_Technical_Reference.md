# TSR (Terrain Solar Radiation) Technical Reference

本文档详细说明 SHUD v2.1 新增的地形太阳辐射修正模块。

## 1. 物理原理

TSR 模块根据地形坡度、坡向和太阳位置，修正每个单元接收的短波辐射：

```
Rn_terrain = Rn_horizontal × TSR_factor
```

其中 TSR_factor 由入射角余弦与太阳高度角余弦的比值决定。

## 2. 配置参数

在 `*.cfg.para` 文件中配置：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `TERRAIN_RADIATION` | 0 | 0=关闭, 1=启用 |
| `SOLAR_UPDATE_INTERVAL` | 60 | 太阳位置更新间隔 (分钟) |
| `RAD_FACTOR_CAP` | 5.0 | TSR因子上限，防止极端几何导致放大过大 |
| `RAD_COSZ_MIN` | 0.05 | cosZ下限截断，避免日出/日落附近数值发散 |
| `TSR_FACTOR_MODE` | INSTANT | TSR 因子语义：INSTANT=瞬时几何因子，FORCING_INTERVAL=forcing 区间等效因子 |
| `TSR_INTEGRATION_STEP_MIN` | 60 | forcing 区间等效因子积分步长 (分钟)，仅在 `TSR_FACTOR_MODE=FORCING_INTERVAL` 时生效 |
| `RADIATION_INPUT_MODE` | SWDOWN | SWDOWN=下行短波，SWNET=净短波 |
| `SOLAR_LONLAT_MODE` | FORCING_FIRST | 经纬度来源选择 |

### SOLAR_LONLAT_MODE 选项

- `FORCING_FIRST`: 使用 forcing 列表第一个站点的经纬度
- `FORCING_MEAN`: 对 forcing 列表中有效经纬度取均值
- `FIXED`: 使用 `SOLAR_LON_DEG` / `SOLAR_LAT_DEG` 指定值

### TSR_FACTOR_MODE 选项（重要）

- `INSTANT`（默认）：在一个时间点（由 `SOLAR_UPDATE_INTERVAL` 的 bucket 对齐）计算太阳位置与 `TSR_factor`，并在该 bucket 内复用。
  - 适用：forcing 较高时间分辨率（小时/子小时），且更关注运行速度时。
- `FORCING_INTERVAL`：对每个 forcing 记录区间 `[t0,t1)` 计算一个等效 `TSR_factor`，并在该 forcing 区间内保持常数。
  - 计算形式为 cosZ 加权平均：`F_eff = (∫ max(cosZ,0) * f(t) dt) / (∫ max(cosZ,0) dt)`，其中 `f(t)` 仍遵循 `terrainFactor()`（含 `RAD_COSZ_MIN`/`RAD_FACTOR_CAP`）。
  - 适用：coarse forcing（尤其逐日短波 forcing），可避免“夜间 factor=0”把区间等效因子显著拉低的问题。

## 3. 诊断输出

当 `TERRAIN_RADIATION=1` 且 `DT_QE_ET > 0` 时，额外输出：

| 文件 | 变量 | 单位 | 说明 |
|------|------|------|------|
| `*.rn_h.dat` | Rn_horizontal | W/m² | 水平面短波辐射 (forcing 原值) |
| `*.rn_t.dat` | Rn_terrain | W/m² | 地形修正后短波辐射 |
| `*.rn_factor.dat` | TSR_factor | - | 几何因子 (cosθ/cosZ) |

**时间戳语义**：输出时间为区间左端点，值为该输出区间内的均值。

## 4. 使用示例

### 4.1 启用 TSR

在 `input/<project>/<project>.cfg.para` 中添加：

```
TERRAIN_RADIATION	1
```

### 4.2 运行模型

```bash
./shud -o output/project_tsr project
```

### 4.3 对比分析

```bash
python3 post_analysis/compare_tsr.py
```

## 5. 验证结果

CCW 流域测试 (1827天模拟)：

### 5.1 辐射修正效果

| 指标 | 值 |
|------|-----|
| 水平面辐射均值 | 220.91 W/m² |
| 地形修正辐射均值 | 216.62 W/m² |
| 辐射增强比例 (ratio > 1) | 42.2% |
| 辐射削弱比例 (ratio < 1) | 57.8% |
| 最大增强倍数 | 2.09x |

### 5.2 水文变量响应

| 变量 | 相对变化 |
|------|----------|
| 潜在蒸散发 (PET) | -1.2% |
| 实际蒸散发 (ET) | -0.5% |
| 蒸发 | -1.1% |
| 非饱和带水量 | +0.2% |
| 河道流量 | -0.8% |

## 6. 源码结构

| 文件 | 功能 |
|------|------|
| `src/ModelData/MD_ET.cpp` | 辐射读取与 TSR 因子应用 |
| `src/Equations/SolarRadiation.cpp` | 太阳位置与地形因子计算 |
| `src/classes/TimeContext.cpp` | 儒略日与时间基准 |
| `src/classes/Model_Control.cpp` | 参数解析与输出头信息 |
| `src/classes/Element.cpp` | 单元地形法向量/坡度/坡向 |

## 7. 限制与注意事项

1. **单一太阳位置近似**：当前使用全局统一的太阳位置，不支持每个单元独立计算
2. **无地形遮蔽**：不考虑地平线遮挡/阴影效果
3. **时区处理**：太阳位置按 UTC 计算，避免经度推断时区的相位偏移

## 8. 相关文档

- [Model Upgrade As-Is vs To-Be](model_upgrade.md)
- [TSR Physical Assumptions](../TSR_Physical_Assumptions.md)
- [Validation Framework](../validation/tsr/README.md)
