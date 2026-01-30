# TSR validation (baseline vs TSR=ON)

此目录用于验证 Terrain Solar Radiation (TSR) 开关对 `ccw` 示例流域输出的影响，保证：
- **不修改** `input/ccw` 原文件（在 `validation/tsr/tmp/` 生成工作副本并在副本上改配置）
- 两组输出互不覆盖：`output/ccw.base/` vs `output/ccw.tsr/`

## 目录结构

- `validation/tsr/run_baseline.sh`：TSR=OFF（`TERRAIN_RADIATION=0`）
- `validation/tsr/run_tsr.sh`：TSR=ON（`TERRAIN_RADIATION=1`）
- `validation/tsr/integration_test.sh`：端到端运行 + 深度校验（读取二进制 `.dat` 输出并检查数值关系）
- `validation/tsr/tmp/`：运行时创建的输入工作副本（已加入 `validation/tsr/.gitignore`）

## 前置条件

- 已编译二进制：仓库根目录存在可执行文件 `./shud`
  - 如果没有：`make shud`
- `python3`（用于 `integration_test.sh` 深度校验）

## 使用方法

从仓库根目录执行：

```bash
bash validation/tsr/run_baseline.sh
bash validation/tsr/run_tsr.sh
```

运行端到端集成验证：

```bash
bash validation/tsr/integration_test.sh
```

## 可调参数（环境变量）

默认行为：
- 复制 `input/ccw` 到 `validation/tsr/tmp/...` 后，**保持** `ccw.cfg.para` 里的 `END` / `DT_QE_ET` 原值（即 `input/ccw/ccw.cfg.para` 的值）
- 仅修改 `TERRAIN_RADIATION`（baseline=0 / tsr=1）

可选：为了加速验证，可通过环境变量覆盖（输出频率绑定在 `DT_QE_ET`）：

```bash
SHUD_VALIDATION_END_DAYS=2 SHUD_VALIDATION_DT_QE_ET_MIN=60 bash validation/tsr/run_tsr.sh
```

`integration_test.sh` 也会继承同样的环境变量；需要加速时可在运行前 `export` 这两个变量。

## 输出位置

- Baseline：`output/ccw.base/`
- TSR ON：`output/ccw.tsr/`

重复运行说明：
- 若目标输出目录已存在且非空，脚本会自动将其移动为 `*.bak.YYYYmmdd-HHMMSS`，再生成新的输出目录
- 模型标准输出/错误输出写入对应目录的 `run.log`

关键 TSR 输出文件（两组目录都会生成）：
- `ccw.rn_h.dat`：水平面太阳辐射（W/m²）
- `ccw.rn_t.dat`：地形修正后辐射（W/m²）
- `ccw.rn_factor.dat`：地形因子（无量纲）

## 故障排查

- 提示 `required file missing: .../shud`：先在仓库根目录编译 `make shud`，并确认 `./shud` 可执行。
- 提示 `required directory missing: .../input/ccw`：确认示例输入目录存在且完整（应包含 `ccw.cfg.para`、`ccw.tsd.forc` 等文件）。
- 提示 `directory not writable` 或无法创建目录：检查 `output/` 和 `validation/tsr/tmp/` 是否有写权限；必要时清理旧目录或修复权限。
- `integration_test.sh` 报 `.dat` 解析/截断：优先查看对应输出目录下的 `run.log`，通常是模型中途退出导致文件不完整。
- `integration_test.sh` 报 “baseline vs tsr ... identical”：查看 `output/ccw.tsr/run.log` 确认 `TERRAIN_RADIATION=1` 生效；确认运行的是最新编译的 `./shud`。
