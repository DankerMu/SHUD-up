# 与远程仓库差异：更改说明（工作区 vs 远程分支）

> 生成时间：2026-01-30  
> 说明：本文档用于汇总**当前工作区（未提交修改）**相对远程仓库分支的差异，并解释关键改动的目的、影响与复现方式。

## 1. 对比基准（Remote refs）

本仓库配置了两个远程：

- `origin`：`git@github.com:DankerMu/SHUD-up.git`
- `upstream`：`git@github.com:SHUD-System/SHUD.git`

本次差异说明以 **`origin/master`** 为主基准（因为当前 `HEAD == origin/master`，差异全部来自未提交的工作区改动）。

- `origin/master`：`bfd9673c29bb0130276a4ac493d68096fd4efa66`
- `upstream/master`：`9b55b0c88fb89eabebdf91b64b0daf9176ad660d`

与 `upstream/master` 的关系（用于背景说明）：
- `HEAD` 相对 `upstream/master` **ahead 22 commits**（`upstream/master` 是当前 `HEAD` 的祖先，非分叉对立）。

## 2. Diff 摘要（相对 `origin/master`）

### 2.1 已修改文件（tracked, modified）

`git diff --numstat origin/master` 统计：

- `README.md`（+2/-0）
- `src/Model/f.cpp`（+1/-0）
- `src/Model/shud.cpp`（+32/-7）
- `src/ModelData/MD_ET.cpp`（+1/-1）
- `src/ModelData/Model_Data.hpp`（+4/-0）
- `validation/tsr/README.md`（+6/-5）
- `validation/tsr/_common.sh`（+31/-5）
- `validation/tsr/integration_test.sh`（+18/-8）

核心含义（按功能分组）：

1) **水量平衡诊断接入求解循环**  
   - `src/Model/shud.cpp`：在 `SHUD()` 主循环中增加水量平衡诊断模块的启用、采样与写盘钩子（受 `SHUD_WB_DIAG=1` 控制）。
   - `src/ModelData/Model_Data.hpp`：新增 `Model_Data::wbdiag` 指针用于在 `ET()` 之后 snapshot `qEleE_IC`。
   - `src/Model/f.cpp`：补充 include（编译依赖）。

2) **拦截库容（IS）守恒性 bug 修复**  
   - `src/ModelData/MD_ET.cpp`：修复拦截库容 `yEleIS` 的“写入乘 `VegFrac`、读取未除回”问题，避免系统性非守恒漂移。

3) **TSR 验证脚本默认行为调整 + 集成测试稳定化**  
   - `validation/tsr/_common.sh`：默认不再强行覆盖 `END/DT_QE_ET`；仅在显式设置环境变量时才覆盖（便于“默认全时段验证”与“短跑加速”共存）。
   - `validation/tsr/integration_test.sh`：默认强制用 `END=2 day, DT_QE_ET=60 min` 做快速集成验证，并使用隔离输出目录 `output/ccw.{base,tsr}.itest/` 防止覆盖你的长时段结果。
   - `validation/tsr/README.md`、`README.md`：补充/更新上述行为说明。

### 2.2 新增文件（untracked, new）

当前工作区新增但尚未提交的文件：

- `docs/water_balance_verification.md`：水量守恒验证报告（全域 + 逐 element），包含默认 1827 天与 2 天短跑的量化结果，以及 TSR 阴阳坡（南北坡）Unsat/GW 差异分析。
- `src/Model/WaterBalanceDiag.hpp`、`src/Model/WaterBalanceDiag.cpp`：水量平衡诊断实现（逐 element 残差 + 全域 basinwbfull）。
- `post_analysis/analyze_aspect_response.py`：基于坡向/坡度分组（南北坡窗口）计算 `rn_t / ETa / yUnsat / yGW / UnsatRatio` 的面积加权时间均值与季节差异。
- `validation/tsr/py/check_water_balance.py`：辅助扫描 `*.elewb*_resid.dat` 的最大残差/异常值（用于排查，不等同于“严格为 0”）。

> 注意：`output/` 下的运行产物目录（例如 `output/ccw.base/`、`output/ccw.tsr/`、`output/ccw.*.itest/`）为运行生成，不属于 git diff 内容。

## 3. 关键行为变化（What changes for users）

### 3.1 新增水量平衡诊断输出（需显式开启）

启用方式：

```bash
SHUD_WB_DIAG=1 ./shud -p <your_project.SHUD>
```

新增输出（位于项目输出目录，如 `output/ccw.tsr/`）：

- `ccw.elewb3_resid.dat`：逐 element 的 A 类残差（`ΔS3 − ∫(dS3/dt)dt`，导数来自 `CVodeGetDky(k=1)`）
- `ccw.elewbfull_resid.dat`：逐 element 的 A 类残差（包含 Snow/IS 的 `Sfull`）
- `ccw.elewb3_budget_resid.dat`：逐 element 的 B 类残差（`ΔS3 − ∫(RHS_budget)dt`，RHS 由通量拼装）
- `ccw.elewbfull_budget_resid.dat`：逐 element 的 B 类残差（包含拦截蒸发、降水与净雨差）
- `ccw.basinwbfull.dat`：全域水量平衡（m³/区间）：`ΔS_total` 与 `P + BC + SS − ET − Qout − Qedge` 的闭合残差

### 3.2 修复拦截库容读写口径不一致导致的非守恒

在 `ET()` 中：

- 写入：`yEleIS = icStg * VegFrac`（表示“按单元面积折算”的平均拦截水深）
- 修复前读取：`icStg = yEleIS`（错误：把“已按 VegFrac 缩放”的量当成 canopy 部分库容）
- 修复后读取：`icStg = yEleIS / VegFrac`（正确：还原到 canopy 覆盖部分的库容变量）

效果：避免降雨/拦截过程引入系统性水量丢失（全域累计残差漂移显著降低）。

### 3.3 TSR 验证脚本行为更贴合“默认全期 + 可选短跑”

- `validation/tsr/run_{baseline,tsr}.sh` 默认**保持** `input/ccw/ccw.cfg.para` 的 `END/DT_QE_ET`，只切换 `TERRAIN_RADIATION`。
- 需要短跑时再显式指定：

```bash
SHUD_VALIDATION_END_DAYS=2 SHUD_VALIDATION_DT_QE_ET_MIN=60 bash validation/tsr/run_tsr.sh
```

- `validation/tsr/integration_test.sh` 默认始终跑 2 天（用于快速 CI / 本地 sanity），且输出到 `output/ccw.*.itest/`，避免污染你的长时段结果目录。

## 4. 复现与验证建议（How to reproduce）

1) 编译：

```bash
make -j4
```

2) 默认全期（ccw）+ 水量平衡诊断：

```bash
SHUD_WB_DIAG=1 bash validation/tsr/run_baseline.sh
SHUD_WB_DIAG=1 bash validation/tsr/run_tsr.sh
```

3) 快速集成验证（默认 2 天，且不覆盖长时段输出）：

```bash
bash validation/tsr/integration_test.sh
```

4) TSR 阴阳坡响应分析：

```bash
python3 post_analysis/analyze_aspect_response.py
```

5) 详细结果与口径说明：

- `docs/water_balance_verification.md`
