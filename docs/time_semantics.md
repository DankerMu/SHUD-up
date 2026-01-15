# SHUD 时间语义（输入 / 运行 / 输出）

本文面向 **SHUD-up 当前实现（以 `src/` 为准）** 的用户，说明“时间轴”在输入、运行与输出中的真实语义与关键约束，避免把相对时间误读为固定日历日期。

## 0. 统一约定：t(min) 是“从基准日开始的分钟数”

- **模型内部时间轴**使用 `t(min)`：单位为 **分钟**，零点为 `t=0`。
- **基准日**由 `ForcStartTime(YYYYMMDD)` 指定（见下文），它提供：
  - 时间轴的“日期标签/基准日”（`t=0` 对应基准日 `00:00`）
  - 太阳几何/儒略日等需要日历信息的计算基准（如 TSR 相关）
- `t(min)` 与相对天数的换算：
  - `Time_min = Time_Day * 1440`
  - `Time_Day = Time_min / 1440`

> 重要：**START/END 不是日历日期**，它们是“相对基准日的天数偏移”，最终会被换算成 `t(min)`。

## 1. Input Time：forcing/TSD 的时间列语义与约束

### 1.1 `ForcStartTime(YYYYMMDD)` 的来源与用途

`ForcStartTime` 从 forcing 列表文件（`<project>.tsd.forc`）首行读取：

```text
# 以 ccw 为例：input/ccw/ccw.tsd.forc
1 20000101
```

- 第 2 个值 `20000101` 即 `ForcStartTime`（基准日）。
- 该值会被要求与各类 TSD 文件的 `StartTime`（头行第 3 列，`YYYYMMDD`）一致，否则直接报错退出（forcing/LAI/MF/各类 BC 均同理）。

相关实现入口：`src/ModelData/MD_readin.cpp`、`src/classes/TimeContext.cpp`

### 1.2 TSD 文件：第一列必须是 `Time_Day`（会乘 1440），不能直接写累计分钟

SHUD 读取各类 TSD（forcing/LAI/MF/BC 等）时，**数据区第一列按“天”读取并乘 1440 转为分钟**（`Time_Day -> Time_min`），因此：

- 第一列应写 **相对天数**（可为小数，表示日内时刻）
- 不要把“累计分钟”直接写进第一列（例如 `60` 会被当成 **60 天**，而不是 60 分钟）

相关实现入口：`src/classes/TimeSeriesData.cpp`（`timeDay * 1440.0`）

### 1.3 `sample-and-hold`（零阶保持）对齐语义

模型在任意时刻 `t(min)` 使用 **最近一条 `Time_min <= t` 的记录值**，并保持到下一条时间戳到来：

- 记录 `Time_Day = d_i` 的值生效区间为：`[d_i, d_{i+1})`（换算到分钟同理）
- 若存在相同 `Time_Day` 的重复时间戳，后读到的会覆盖先前值（时间列允许“非递减”）

相关实现入口：`src/ModelData/MD_update.cpp`（`movePointer()` 调用链）、`src/classes/TimeSeriesData.cpp`（`movePointer()`）

### 1.4 覆盖不足策略（强约束）

- forcing（逐站点 forcing csv）必须覆盖整个模拟区间 `[CS.StartTime, CS.EndTime]`（单位分钟），否则直接报错退出，并提示“延长 forcing 或调整 START/END”。
- 时间列必须 **单调非递减**（允许相等），否则直接报错退出。

相关实现入口：`src/ModelData/MD_readin.cpp`（`validateTimeStamps()`）、`src/classes/TimeSeriesData.cpp`（单调性校验）

### 1.5 ccw forcing.csv：2–3 行示例（Time_Day -> t(min)）

`input/ccw/forcing.csv` 是单站点 forcing TSD 文件，其头两行与前三条数据如下：

```text
87667 6 20000101 20100101
Time_Day    APCP  TMP   SPFH   UGRD  DSWRF
0           0     9.56  0.4758 2.03  88.6
0.0416666666666667 0 7.66 0.5266 2.04 0
0.0833333333333333 0 5.75 0.5841 2.06 0
```

解释（以基准日 `ForcStartTime=20000101` 为例）：

- `Time_Day=0` → `t=0` min → 2000-01-01 00:00
- `Time_Day=1/24=0.041666...` → `t=60` min → 2000-01-01 01:00
- `Time_Day=2/24=0.083333...` → `t=120` min → 2000-01-01 02:00

## 2. Runtime：START/END 的真实语义（相对天数偏移 → t(min)）

模拟起止由 `<project>.cfg.para` 中的 `START` / `END` 控制（单位：**天**，可为小数）：

```text
# 以 ccw 为例：input/ccw/ccw.cfg.para
START 0
END   1827
```

运行时会做以下换算（单位换算发生在控制参数读取后）：

- `CS.StartTime = START * 1440`（分钟）
- `CS.EndTime   = END * 1440`（分钟）

因此：

- `START=0` 表示从基准日 `t=0` 开始（不是“某年某月某日”）
- 若你想从基准日之后第 `k` 天开始跑，就设 `START=k`（或设小数表示日内时刻）

相关实现入口：`src/classes/Model_Control.cpp`（`updateSimPeriod()`）

## 3. Output Time：`Time_min` / `Time_Day` 与 `t-Interval`（左端点）语义

### 3.1 时间列：`Time_min` 与 `Time_Day`

ASCII 输出（`.csv`）的第一列固定为：

- `Time_min`：输出时间戳（分钟），与输入/运行共用同一 `t(min)` 时间轴
- 需要“天”时，可自行换算：`Time_Day = Time_min / 1440`

相关实现入口：`src/classes/Model_Control.cpp`（`Print_Ctrl::open_file()` / `fun_printASCII()`）

### 3.2 `t-Interval`：左端点时间戳 + 输出触发策略

每个输出文件都有一个输出间隔 `Interval`（单位分钟，来自 `DT_*` 参数，例如 `DT_YE_SURF=1440` 表示逐日输出）。

输出采用 **区间左端点时间戳**：

- 当模型时间到达某个输出边界（`t` 落在 `Interval` 的整数倍附近）时触发写出
- 写出的时间戳不是“当前 t”，而是：`Time_min = t - Interval`（左端点）
- 因此第一条记录常见为 `Time_min=0.0`，它代表区间 `[0, Interval)` 的统计量

实现细节（便于对齐排查）：触发条件约为 `floor(t + 0.001) % Interval == 0`，写出时间为 `floor(t + 0.001) - Interval`。

相关实现入口：`src/classes/Model_Control.cpp`（`Print_Ctrl::PrintData()`）

### 3.3 输出头里的基准日

输出 `.csv` 的第二行包含基准日（`YYYYMMDD`），用于把 `Time_min` 映射回日历时间轴：

```text
0    <NumVar>   <ForcStartTime>
```

该 `<ForcStartTime>` 与输入侧保持一致（来自 forcing 列表文件），用于“时间轴基准日”标签。

## 4. 源码索引（高级用户）

- 基准日与日历换算：`src/classes/TimeContext.hpp`、`src/classes/TimeContext.cpp`
- forcing 与时间戳一致性/覆盖校验：`src/ModelData/MD_readin.cpp`
- TSD 时间列解析（`Time_Day * 1440`）与 sample-and-hold：`src/classes/TimeSeriesData.cpp`
- START/END（天→分钟）与输出时间戳（左端点 `t-Interval`）：`src/classes/Model_Control.cpp`

