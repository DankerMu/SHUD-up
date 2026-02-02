# 与远程仓库差异：更改说明（工作区 vs 远程分支）

> 生成时间：2026-01-30  
> 说明：本文档用于汇总**当前工作区（未提交修改）**相对远程仓库分支的差异，并解释关键改动的目的、影响与复现方式。

## 1. 对比基准（Remote refs）

本仓库配置了两个远程：

- `origin`：`git@github.com:DankerMu/SHUD-up.git`
- `upstream`：`git@github.com:SHUD-System/SHUD.git`

本次差异说明以 **`origin/master`** 为主基准（因为当前 `HEAD == origin/master`，差异全部来自未提交的工作区改动）。

- `origin/master`：`a705d2f967b982bf3228f42481da81ab3e3f43d3`
- `upstream/master`：`9b55b0c88fb89eabebdf91b64b0daf9176ad660d`

与 `upstream/master` 的关系（用于背景说明）：
- `HEAD` 相对 `upstream/master` **ahead 26 commits**（`upstream/master` 是当前 `HEAD` 的祖先，非分叉对立）。

## 2. Diff 摘要（相对 `origin/master`）

### 2.1 已修改文件（tracked, modified）

`git diff --numstat origin/master` 统计：

- `.gitignore`（+6/-0）
- `docs/remote_diff_notes.md`（+38/-30）
- `docs/water_balance_verification.md`（+126/-2）
- `post_analysis/plot_water_balance_verification.py`（+307/-0）
- `src/Model/WaterBalanceDiag.cpp`（+203/-12）
- `src/Model/WaterBalanceDiag.hpp`（+41/-0）
- `src/Model/shud.cpp`（+16/-4）
- `src/ModelData/MD_f.cpp`（+61/-0）
- `validation/tsr/_common.sh`（+4/-0）

核心含义（按功能分组）：

1) **水量平衡诊断：积分方法增强（默认关闭）**  
   - 新增两种“更稳健”的区间积分方式并纳入开关控制：
     - `SHUD_WB_DIAG_TRAPZ=1`：对外步采样的 rate 做**梯形积分**（低开销）
     - `SHUD_WB_DIAG_QUAD=1`：用 `CV_ONE_STEP` 驱动 CVODE，并在**内部自适应步长**上累计 basin 通量；额外输出 `*.basinwbfull_quad.dat`（更慢）
   - `src/Model/WaterBalanceDiag.*`：实现 TRAPZ 与 basin 内部步积分累计，并写盘输出。
   - `src/ModelData/MD_f.cpp`：在 RHS 计算中同步计算并缓存 basin 通量 rate（供 QUAD 模式使用）。
   - `src/Model/shud.cpp`：在 QUAD 模式下改用 `CV_ONE_STEP` 循环并在每个内部步刷新通量场，确保积分 rate 与 accepted solution 对齐。

2) **可视化与报告增强**  
   - `post_analysis/plot_water_balance_verification.py`：新增
     - `MAX_SOLVER_STEP` 收敛对比图（CDF / zoom / cumsum）
     - 诊断积分方法对比图（BE vs TRAPZ vs QUAD；以及 QUAD-run 内 sampled vs internal-step）
   - `docs/water_balance_verification.md`：新增 §6.5（诊断积分方法对比）并补充图表、CDF 含义解释与结论。

3) **仓库卫生（忽略本地/生成产物）**  
   - `.gitignore`：忽略本地编辑器配置（`.obsidian/`）与本地生成的 `docs/water_balance_verification.pdf`（以 `.md` 为源文件）。

### 2.2 新增文件（untracked, new）

当前工作区新增但尚未提交的文件（主要为图与脚本）：

- `validation/tsr/run_tsr_maxstep_experiments.sh`：生成 `MAX_SOLVER_STEP` 收敛性实验输出（`ms5/ms2`）。
- `validation/tsr/run_tsr_wb_integrators.sh`：生成 BE/TRAPZ/QUAD 三组“诊断积分方法对比”输出目录。
- `docs/figures/water_balance/` 下新增图：
  - `basin_residual_convergence_*.png`
  - `basin_residual_integrators_*.png`
  - `basin_residual_quadrun_*.png`

> 注意：`output/` 下的运行产物目录为运行生成（本仓库默认忽略），不属于 git diff 内容。

## 3. 关键行为变化（What changes for users）

### 3.1 水量平衡诊断输出开关与新增输出（需显式开启）

启用方式（默认关闭，避免影响运行速度）：

```bash
SHUD_WB_DIAG=1 ./shud -p <your_project.SHUD>
```

新增/增强开关：

- `SHUD_WB_DIAG_TRAPZ=1`：对采样积分使用梯形积分（低开销）。
- `SHUD_WB_DIAG_QUAD=1`：在 CVODE 内部步长上累计 basin 通量，并输出 `*.basinwbfull_quad.dat`（更慢）。

主要输出（位于项目输出目录，如 `output/ccw.tsr/`）：

- `ccw.elewb3_resid.dat`：逐 element 的 A 类残差（`ΔS3 − ∫(dS3/dt)dt`，导数来自 `CVodeGetDky(k=1)`）
- `ccw.elewbfull_resid.dat`：逐 element 的 A 类残差（包含 Snow/IS 的 `Sfull`）
- `ccw.elewb3_budget_resid.dat`：逐 element 的 B 类残差（`ΔS3 − ∫(RHS_budget)dt`，RHS 由通量拼装）
- `ccw.elewbfull_budget_resid.dat`：逐 element 的 B 类残差（包含拦截蒸发、降水与净雨差）
- `ccw.basinwbfull.dat`：全域水量平衡（m³/区间）：`ΔS_total` 与 `P + BC + SS − ET − Qout − Qedge` 的闭合残差
- `ccw.basinwbfull_quad.dat`：全域水量平衡（m³/区间）：与上一条同口径，但通量积分使用**内部步长**（仅在 `SHUD_WB_DIAG_QUAD=1` 时输出）

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
