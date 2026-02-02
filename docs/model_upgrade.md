# Model Upgrade：As-Is（当前实现） vs To-Be（计划实现）

本文档用于把“规划/设想”与“仓库当前真实实现”明确区分，避免新读者把计划误读为已完成。

- **As-Is（当前实现）**：以 `src/` 为准（本文中的 As-Is 结论均可在源码中定位）。
- **To-Be（计划实现）**：以 `docs/spec.md` 与 `epic/`（Epic backlog）为准。

## 文档入口（必读）

- 计划与范围（Spec）：[`docs/spec.md`](spec.md)
- Epic backlog（需求拆分与追踪）：[`epic/README.md`](../epic/README.md)
- TSR 验证（已完成工作）：[`validation/tsr/README.md`](../validation/tsr/README.md)
- TSR 物理假设与尺度限制（背景）：[`TSR_Physical_Assumptions.md`](../TSR_Physical_Assumptions.md)

## Current Status (As-Is)

> 本节描述“现在仓库里已经有什么”，以 `src/` 与 `validation/` 为准。

### 1) Terrain Solar Radiation (TSR) 开关与计算链路

**现状：已实现（可用）**

- **开关**：`TERRAIN_RADIATION`（0/1），在参数文件中配置后生效。
- **forcing 区间等效因子**：TSR 因子按 forcing 记录区间 `[t0,t1)` 计算，并在该 forcing 区间内复用（避免 coarse forcing 下“夜间 factor=0”稀释区间效应）。
- **区间内采样步长**：`TSR_INTEGRATION_STEP_MIN`（分钟）；控制 `[t0,t1)` 内太阳几何采样的步长（逐日短波 forcing 默认 60 分钟即 24 个采样点）。
- **数值边界参数**：
  - `RAD_FACTOR_CAP`：TSR 因子上限（防止极端几何导致放大过大）
  - `RAD_COSZ_MIN`：分母 `cosZ` 的下限截断（避免日出/日落附近数值发散）
- **短波辐射输入语义**：`RADIATION_INPUT_MODE`
  - `SWDOWN`（默认）：forcing 第 6 列为下行短波，模型内部会再乘 `(1 - Albedo)` 得到净短波
  - `SWNET`：forcing 第 6 列已是净短波，不再乘 `(1 - Albedo)`
- **全局太阳经纬度选择**：`SOLAR_LONLAT_MODE`
  - `FORCING_FIRST`（默认）：使用 forcing 列表第一个站点的经纬度
  - `FORCING_MEAN`：对 forcing 列表中有效经纬度取均值
  - `FIXED`：使用 `SOLAR_LON_DEG` / `SOLAR_LAT_DEG`
- **时间基准与时区约定**：
  - 基准日期来自 forcing 列表头部的 `ForcStartTime`（`YYYYMMDD`），用于计算儒略日/太阳几何
  - TSR 的太阳位置计算按 UTC 对齐（显式传入 `timezone_hours=0.0`），避免由经度推断时区造成相位偏移

**源码入口（可定位）**

- 辐射读取与 TSR 因子应用：`src/ModelData/MD_ET.cpp`
- 太阳位置与地形因子：`src/Equations/SolarRadiation.cpp`、`src/Equations/SolarRadiation.hpp`
- 基准日期/儒略日（供太阳几何使用）：`src/classes/TimeContext.cpp`、`src/classes/TimeContext.hpp`
- 单元地形法向量/坡度/坡向：`src/classes/Element.cpp`、`src/classes/Element.hpp`
- 参数解析与输出头信息：`src/classes/Model_Control.cpp`、`src/classes/Model_Control.hpp`

### 2) TSR 相关输出（用于诊断/验证）

**现状：已实现（可用）**

- 当 `DT_QE_ET > 0`（启用 element ET 输出）时，会额外输出三类 TSR 诊断量：
  - `*.rn_h`：水平面短波（forcing 原值，W/m²）
  - `*.rn_t`：地形修正后短波（W/m²）
  - `*.rn_factor`：TSR 因子（无量纲）
- 输出文件既支持二进制 `.dat`（默认）也支持 ASCII `.csv`（取决于 `ASCII/BINARY` 配置）。
- **时间戳语义**：输出时间采用区间左端点（`t - Interval`），值为该输出区间内的**均值**（见 `Print_Ctrl::PrintData`）。

### 3) TSR 验证现状（已完成工作）

**现状：已完成一套可复现实证/回归流程**

- `validation/tsr/` 提供 baseline（TSR=OFF）与 TSR=ON 的端到端运行脚本与深度校验：
  - 目录说明与使用方法：[`validation/tsr/README.md`](../validation/tsr/README.md)
  - Python 参考实现（`.dat` 读取 + TSR 计算对照）与单元测试：[`validation/tsr/py/`](../validation/tsr/py/)

### 4) 明确“当前没有”的内容（避免误读）

以下内容**在当前仓库中并未实现**（若在其他文档中出现，应视为 To-Be/提案，或需要补充实现后再宣称“已完成”）：

- **分布式太阳位置**（例如“每个 element 使用自身中心点经纬度计算太阳位置”）：当前 TSR 使用单个全局太阳位置近似。
- **地形遮蔽/地平线遮挡（horizon/shadowing）**：当前 `terrainFactor()` 为坡面入射几何修正，不包含视域/遮挡建模。
- **C++ 单元测试框架/`tests/` 目录**：当前主要采用 `validation/` 下的脚本与 Python `unittest` 做验证。

## Planned Features (To-Be)

> 本节描述“计划做什么”，不代表已经实现；请以 `docs/spec.md` 与 `epic/README.md` 为准。

- 计划与需求基线：[`docs/spec.md`](spec.md)
- 任务拆分与追踪（Epic backlog）：[`epic/README.md`](../epic/README.md)

为避免 To-Be 与 As-Is 混淆，本文档仅保留摘要；所有 To-Be 的细节、边界条件与验收标准请在上述两处维护。
