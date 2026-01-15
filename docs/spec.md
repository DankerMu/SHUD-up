# SHUD-up Spec（To-Be / Draft）

本文件描述 **To-Be（计划实现）**：我们希望 SHUD-up 在后续迭代中达成的目标与约束。

- **当前实现（As-Is）**：请以 [`docs/model_upgrade.md`](model_upgrade.md) 的 As-Is 部分与 `src/` 为准。
- **任务拆分与追踪（Epic backlog）**：请以 [`epic/README.md`](../epic/README.md) 为准。

## 目标（Goals）

- 让 TSR（Terrain Shortwave Radiation）能力“可解释、可配置、可验证”，并且在不引入过高计算开销的前提下提升适用范围。
- 让关键输出具备可追溯性（运行时能知道关键开关/经纬度选择等元数据），便于复现实验与回归比较。
- 形成持续可用的验证基线（validation workflow），将“研究原型”沉淀为可复现流程。

## 非目标（Non-goals）

- 不在短期内引入完整辐射传输/遮蔽/地平线视域等重物理模型（如需应另立 Epic）。
- 不强制引入新的 C++ 测试框架；当前以 `validation/` 工作流为主（如需 CI 化应另立 Epic）。

## To-Be：计划功能与改进方向（摘要）

> 具体拆分到 issue 的粒度与优先级由 Epic backlog 管理；此处仅给出方向与边界。

### 1) TSR：太阳位置空间分布能力

现状使用“单个全局太阳位置”近似（见 [`TSR_Physical_Assumptions.md`](../TSR_Physical_Assumptions.md)），计划提供可选的更高精度模式：

- **分区/分带太阳位置**：对流域划分少量区域，每区计算独立太阳位置与/或 TSR 因子，用插值或分配方式作用到单元，兼顾精度与性能。
- **ELEMENT_CENTER（提案）**：每个 element 使用自身中心点经纬度计算太阳位置，并保留时间桶更新策略控制开销。
  - 前置条件：需要明确 element 的经纬度来源（例如输入数据提供，或定义投影与转换方法）。

### 2) TSR：参数与边界行为明确化

- 继续保持关键参数可配置（例如 `SOLAR_UPDATE_INTERVAL`、`RAD_FACTOR_CAP`、`RAD_COSZ_MIN`）。
- 对低太阳高度角附近的边界行为（例如 `RAD_COSZ_MIN` 截断导致的数值特性）在文档中给出明确解释与建议默认值范围。

### 3) 验证与回归（Validation / Regression）

- 扩展 [`validation/tsr/`](../validation/tsr/) 覆盖更多示例流域/工况，形成“可重复跑”的回归集合。
- 为关键数值关系定义验收规则（例如“TSR=ON 输出应与 OFF 有显著差异，但形状/维度一致”），并在脚本中固化。
- 如引入 CI（可选）：优先做“轻量级”校验（小步长/短时间窗），避免拉长构建时间。

### 4) 文档一致性（Docs）

- `docs/model_upgrade.md` 始终保持 As-Is/To-Be 对齐：As-Is 以 `src/` 为准；To-Be 以本 spec 与 Epic backlog 为准。
- 对“尚未实现/未来提案/待定”的内容必须显式标注，避免误导新读者。
