# TSR (Terrain Solar Radiation) Technical Reference

本文档详细说明 SHUD 的地形太阳辐射修正模块（TSR, Terrain Solar Radiation）的**当前实现**：

- forcing 短波被视为“水平面短波通量（W/m²）”
- TSR 将其几何修正为“坡面有效入射短波通量（W/m²）”
- 修正系数按 forcing 的时间间隔（$[t_0,t_1)$）计算**区间等效因子**，避免“固定 60min bucket”在逐日/粗时间步 forcing 下引入系统误差

## 1. 物理原理

TSR 模块根据地形坡度、坡向和太阳位置，修正每个单元接收的短波辐射。其核心是把“水平面下行短波”换算为“坡面有效入射短波”：

$$
R_{n,t} = R_{n,h}\,F_{\mathrm{eff}}
$$

其中 $F_{\mathrm{eff}}$ 本质上来自几何关系（同一束太阳光在坡面与水平面上的投影比值；见 §6）。

以“仅直射束”为例，若直射法向辐照度为 $I_n$，则水平面接收 $I_n\cos Z$、坡面接收 $I_n\cos i$，因此瞬时几何因子为：

$$
f(t) \approx \frac{\cos i(t)}{\cos Z(t)}
$$

TSR 的作用是把 forcing 中的“水平面短波通量”按该几何关系做**相对修正**（注意是 $\cos i/\cos Z$，不是再乘一次 $\cos Z$）。

### 1.1 记号与坐标系

- 时间 $t$：以分钟计，$t=0$ 对应 forcing 列表文件中的 `ForcStartTime`（yyyymmdd），详见 `src/ModelData/MD_readin.cpp` 与 `src/classes/TimeContext.*`。
- 坐标系：$x$=East（东）、$y$=North（北）、$z$=Up（天顶）。
- $\mathbf{n}=(n_x,n_y,n_z)$：单元地形面单位法向量（$n_z\ge 0$）。
- 太阳位置（`SolarPosition`）：
  - $\cos Z$：太阳天顶角余弦（水平面的入射余弦）。
  - $\mathrm{azimuth}$：太阳方位角（弧度），定义为 $\operatorname{atan2}(east,\,north)$，范围 $[0,2\pi)$，即 North=0、East=$\pi/2$。
- 太阳方向单位向量 $\mathbf{s}=(s_x,s_y,s_z)$：由 $\cos Z$ 与 $\mathrm{azimuth}$ 转换得到（见 §4.3）。
- 坡面入射余弦：$\cos i = \mathbf{n}\cdot \mathbf{s}$。

## 2. 配置参数

在 `*.cfg.para` 文件中配置：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `TERRAIN_RADIATION` | 0 | 0=关闭, 1=启用 |
| `RAD_FACTOR_CAP` | 5.0 | TSR因子上限，防止极端几何导致放大过大 |
| `RAD_COSZ_MIN` | 0.05 | $\cos Z$ 下限截断，避免日出/日落附近数值发散 |
| `TSR_INTEGRATION_STEP_MIN` | 60 | forcing 区间等效因子积分步长 (分钟) |
| `RADIATION_INPUT_MODE` | SWDOWN | SWDOWN=下行短波，SWNET=净短波 |
| `SOLAR_LONLAT_MODE` | FORCING_FIRST | 经纬度来源选择 |

### SOLAR_LONLAT_MODE 选项

- `FORCING_FIRST`: 使用 forcing 列表第一个站点的经纬度
- `FORCING_MEAN`: 对 forcing 列表中有效经纬度取均值
- `FIXED`: 使用 `SOLAR_LON_DEG` / `SOLAR_LAT_DEG` 指定值

### 2.1 TSR 因子语义（当前实现：forcing-interval 等效）

TSR 因子按 forcing 的“当前记录区间” $[t_0,t_1)$ 计算一个等效因子 $F_{\mathrm{eff}}$，并在该 forcing 区间内保持常数。

区间等效形式为 $\cos Z$ 加权平均：

$$
F_{\mathrm{eff}}
= \frac{\int_{t_0}^{t_1}\max(\cos Z(t),0)\,f(t)\,dt}{\int_{t_0}^{t_1}\max(\cos Z(t),0)\,dt}
$$

其中 $f(t)$ 为瞬时几何因子（遵循 `terrainFactor()`，含 `RAD_COSZ_MIN`/`RAD_FACTOR_CAP`；见 §5）。

`TSR_INTEGRATION_STEP_MIN` 控制 $[t_0,t_1)$ 内采样步长；对逐日短波 forcing，默认 60 分钟对应每天 24 个采样点。

> 兼容性说明：历史参数 `SOLAR_UPDATE_INTERVAL` / `TSR_FACTOR_MODE` 已废弃；若仍出现在配置中会被解析，但不再改变 TSR 的计算方法。

## 3. 地形几何（$\mathbf{n}$ 向量）

每个 element 的地形面单位法向量 $\mathbf{n}=(n_x,n_y,n_z)$ 来自三角形单元的三节点坐标（使用表面高程 `zmax`）：

给定三点 $\mathbf{p}_1=(x_1,y_1,z_1)$、$\mathbf{p}_2=(x_2,y_2,z_2)$、$\mathbf{p}_3=(x_3,y_3,z_3)$：

1. 边向量：$\mathbf{v}_1=\mathbf{p}_2-\mathbf{p}_1$，$\mathbf{v}_2=\mathbf{p}_3-\mathbf{p}_1$
2. 法向量：$\mathbf{n}_{raw} = \mathbf{v}_1 \times \mathbf{v}_2$
3. 单位化：$\mathbf{n} = \mathbf{n}_{raw} / \lVert \mathbf{n}_{raw} \rVert$
4. 保证 $n_z\ge 0$：若 $n_z<0$ 则 $\mathbf{n}=-\mathbf{n}$（使法向量朝上）

对应源码：`src/classes/Element.cpp`。

## 4. 太阳位置（solarPosition）与太阳方向向量

### 4.1 太阳位置输入：时间与经纬度

- 时间基准：`ForcStartTime` 来自 forcing 列表文件（`*.tsd.forc` 第一行），模型把它设置为 `TimeContext` 的基准日期（yyyymmdd）。
- 经纬度：由 `SOLAR_LONLAT_MODE` 从 forcing 列表（Lon/Lat 列）或固定值选择，并存入 `CS.solar_lon_deg/CS.solar_lat_deg`。
- **时区**：TSR 明确把 forcing 时间当作 **UTC**，因此在计算 solarPosition 时显式传入 `timezone_hours = 0.0`，避免按经度推断时区（`round(lon/15)`）造成相位偏移。

### 4.2 solarPosition 使用的公式（与 `src/Equations/SolarRadiation.cpp` 一致）

设：

- 纬度 $\varphi$（rad），经度 $\lambda$（deg，范围 wrap 到 $[-180,180]$）
- $\mathrm{doy}$：年内第几天（1–366）
- $\mathrm{hour}$：当天小时（$\mathrm{minuteOfDay}(t)/60$）
- $\gamma$：fractional year（rad）
- $E$：equation of time（min）
- $\delta$：declination（rad）
- $\mathrm{tz}$：时区（hours；TSR 固定为 0）
- $\mathrm{TST}$：true solar time（min）
- $\omega$：hour angle（rad）

计算：

$$
\begin{aligned}
\gamma &= \frac{2\pi}{365}\left((\mathrm{doy}-1) + \frac{\mathrm{hour}-12}{24}\right) \\
E &= 229.18\Bigl(0.000075 + 0.001868\cos\gamma - 0.032077\sin\gamma \\
&\quad\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\; - 0.014615\cos 2\gamma - 0.040849\sin 2\gamma\Bigr) \\
\delta &= 0.006918 - 0.399912\cos\gamma + 0.070257\sin\gamma \\
&\quad - 0.006758\cos 2\gamma + 0.000907\sin 2\gamma \\
&\quad - 0.002697\cos 3\gamma + 0.00148\sin 3\gamma \\
\mathrm{time\_offset} &= E + 4\lambda - 60\,\mathrm{tz} \\
\mathrm{TST} &= \operatorname{wrap}_{1440}\!\left(\mathrm{minuteOfDay}(t) + \mathrm{time\_offset}\right) \\
\omega &= \operatorname{deg2rad}\!\left(\mathrm{TST}/4 - 180\right) \\
\cos Z &= \sin\varphi\,\sin\delta + \cos\varphi\,\cos\delta\,\cos\omega \\
Z &= \arccos(\cos Z) \\
east &= -\cos\delta\,\sin\omega \\
north &= \cos\varphi\,\sin\delta - \sin\varphi\,\cos\delta\,\cos\omega \\
\mathrm{azimuth} &= \operatorname{wrap}_{2\pi}\!\left(\operatorname{atan2}(east,\,north)\right)
\end{aligned}
$$

### 4.3 由 ($\cos Z$, $\mathrm{azimuth}$) 得到太阳方向向量 $\mathbf{s}$

令 $\sin Z = \sqrt{\max\left(0,\,1-\cos^2 Z\right)}$，则：

$$
\begin{aligned}
s_x &= \sin Z\,\sin(\mathrm{azimuth}) \\
s_y &= \sin Z\,\cos(\mathrm{azimuth}) \\
s_z &= \cos Z
\end{aligned}
$$

## 5. 瞬时 TSR 因子 $f(t)$

对任意时刻 $t$（且太阳在地平线上方），瞬时几何因子 $f(t)$ 与 `terrainFactor()` 一致：

1. 坡面入射余弦：$\cos i = \mathbf{n}\cdot\mathbf{s}$
2. 条件与截断：

$$
f(t)=
\begin{cases}
0, & \cos Z(t)\le 0 \\
0, & \cos i(t)\le 0 \\
\min\!\left(\dfrac{\cos i(t)}{\max(\cos Z(t),\,\mathrm{RAD\_COSZ\_MIN})},\,\mathrm{RAD\_FACTOR\_CAP}\right), & \text{otherwise}
\end{cases}
$$

说明：

- `RAD_FACTOR_CAP`：防止在日出/日落附近 $\cos Z$ 很小导致 $\cos i/\cos Z$ 过大。
- `RAD_COSZ_MIN`：分母下限（默认 0.05，对应 $Z\approx 87.13^\circ$），用于数值稳定；这会导致**极低太阳高度**时即使水平面也可能出现 $f(t)<1$，但此时 $\cos Z$ 权重很小，对日尺度均值影响通常可忽略。

> 水平面一致性：当 $\mathbf{n}=(0,0,1)$ 且 $\cos Z\ge \mathrm{RAD\_COSZ\_MIN}$ 时有 $f(t)=1$；只有在 $\cos Z<\mathrm{RAD\_COSZ\_MIN}$（极低太阳高度角）时会因分母截断产生 $f(t)<1$。

## 6. forcing 区间等效 TSR 因子 $F_{\mathrm{eff}}$（当前实现）

### 6.1 forcing 区间定义

weather forcing CSV 的每条记录带有时间戳 $t_k$（以 day 为单位存储，内部转为分钟），模型把该记录视为在区间 $[t_k, t_{k+1})$ 内保持常数（通过 `TimeSeriesData::movePointer()` 切换到下一条记录）。

因此，设当前 forcing 区间为 $[t_0,t_1)$（分钟），forcing 中的水平面短波值为常数 $R_{n,h}$（W/m²）。

### 6.2 在 $[t_0,t_1)$ 内积分近似

令：

$$
\begin{aligned}
\Delta t_{\mathrm{forc}} &= t_1 - t_0 \\
\Delta t_{\mathrm{int}}  &= \min(\mathrm{TSR\_INTEGRATION\_STEP\_MIN},\,\Delta t_{\mathrm{forc}}) \\
n &= \left\lceil \frac{\Delta t_{\mathrm{forc}}}{\Delta t_{\mathrm{int}}} \right\rceil \quad (n\ge 1) \\
\Delta t_{\mathrm{seg}} &= \frac{\Delta t_{\mathrm{forc}}}{n} \\
t_k &= t_0 + (k+0.5)\,\Delta t_{\mathrm{seg}},\quad k=0,\dots,n-1
\end{aligned}
$$

对每个 $t_k$ 计算 $\cos Z_k$、$\mathrm{azimuth}_k$、$f_k=f(t_k)$。

### 6.3 $\cos Z$ 加权的区间等效因子

当前实现使用 $w(t)=\max(\cos Z(t),0)$ 作为权重（近似把短波日变化形状视为与 $\cos Z$ 同相）：

$$
\begin{aligned}
w_k &= \max(\cos Z_k,\,0)\,\Delta t_{\mathrm{seg}} \\
F_{\mathrm{eff}} &= \frac{\sum_k w_k f_k}{\sum_k w_k}
\end{aligned}
$$

若 $\sum_k w_k = 0$（整段为夜间或极端情况）则 $F_{\mathrm{eff}}=0$。

**最终在该 forcing 区间内使用常数因子：**

$$
R_{n,t} = R_{n,h}\,F_{\mathrm{eff}}
$$

对应源码：`src/ModelData/MD_ET.cpp`（区间采样/加权/缓存）与 `src/Equations/SolarRadiation.cpp`（`solarPosition()`/`terrainFactor()`）。

## 7. 与 RADIATION_INPUT_MODE 的关系（SWDOWN vs SWNET）

在 `tReadForcing()` 中，TSR 总是在短波输入层面先应用到 $R_{n,h}$：

$$
R_{n,t} = R_{n,h}\,F_{\mathrm{eff}}
$$

随后根据 `RADIATION_INPUT_MODE` 决定是否乘以 $(1-\mathrm{Albedo})$：

### 7.1 `SWDOWN`（默认）

forcing 第 6 列为下行短波 $SW_{\downarrow}$（W/m²），模型内部净短波为：

$$
R_{n,\mathrm{sw,net}} = R_{n,t}\,(1-\mathrm{Albedo})
$$

### 7.2 `SWNET`

forcing 第 6 列已是净短波（已包含地表反照率等效），模型内部净短波为：

$$
R_{n,\mathrm{sw,net}} = R_{n,t}
$$

> 注意：本模块只做“坡面/水平面”的几何修正，不区分直射/散射分量，也不做天空视域/地形遮蔽修正。

## 8. 诊断输出与时间戳语义

当 `TERRAIN_RADIATION=1` 且 `DT_QE_ET > 0` 时，额外输出：

| 文件 | 变量 | 单位 | 说明 |
|------|------|------|------|
| `*.rn_h.dat` | `Rn_horizontal` | W/m² | 水平面短波辐射（forcing 原值） |
| `*.rn_t.dat` | `Rn_terrain` | W/m² | 地形修正后短波辐射 |
| `*.rn_factor.dat` | `TSR_factor` | - | 区间等效几何因子 $F_{\mathrm{eff}}$ |

**时间戳语义**：输出时间为区间左端点，值为该输出区间内的均值（见 `src/classes/Model_Control.cpp::Print_Ctrl::PrintData()`）。

设输出间隔为 $\Delta t_{\mathrm{out}}$，则输出是时间平均：

$$
\begin{aligned}
rn_{h,\mathrm{out}}(t)      &= \frac{1}{\Delta t_{\mathrm{out}}}\int_{t}^{t+\Delta t_{\mathrm{out}}} rn_h(\tau)\,d\tau \\
rn_{\mathrm{factor,out}}(t) &= \frac{1}{\Delta t_{\mathrm{out}}}\int_{t}^{t+\Delta t_{\mathrm{out}}} F_{\mathrm{eff}}(\tau)\,d\tau \\
rn_{t,\mathrm{out}}(t)      &= \frac{1}{\Delta t_{\mathrm{out}}}\int_{t}^{t+\Delta t_{\mathrm{out}}} rn_h(\tau)\,F_{\mathrm{eff}}(\tau)\,d\tau
\end{aligned}
$$

因此一般有 $rn_{t,\mathrm{out}} \ne rn_{h,\mathrm{out}}\;rn_{\mathrm{factor,out}}$（两者存在协方差；尤其当 forcing 为子日尺度而输出为日尺度时更明显）。

## 9. 使用示例

### 9.1 启用 TSR

在 `input/<project>/<project>.cfg.para` 中添加：

```
TERRAIN_RADIATION	1
```

### 9.2 运行模型

```bash
./shud -o output/project_tsr project
```

### 9.3 对比分析

```bash
python3 post_analysis/compare_tsr.py
```

## 10. 验证结果（示例）

CCW 流域测试 (1827天模拟)：

### 10.1 辐射修正效果

| 指标 | 值 |
|------|-----|
| 水平面辐射均值 | 220.91 W/m² |
| 地形修正辐射均值 | 216.62 W/m² |
| 辐射增强比例 (ratio > 1) | 42.2% |
| 辐射削弱比例 (ratio < 1) | 57.8% |
| 最大增强倍数 | 2.09x |

### 10.2 水文变量响应

| 变量 | 相对变化 |
|------|----------|
| 潜在蒸散发 (PET) | -1.2% |
| 实际蒸散发 (ET) | -0.5% |
| 蒸发 | -1.1% |
| 非饱和带水量 | +0.2% |
| 河道流量 | -0.8% |

## 11. 源码结构（关键实现点）

| 文件 | 功能 |
|------|------|
| `src/ModelData/MD_ET.cpp` | 辐射读取与 TSR 因子应用 |
| `src/Equations/SolarRadiation.cpp` | 太阳位置与地形因子计算 |
| `src/classes/TimeContext.cpp` | 儒略日与时间基准 |
| `src/classes/Model_Control.cpp` | 参数解析与输出头信息 |
| `src/classes/Element.cpp` | 单元地形法向量/坡度/坡向 |

## 12. 限制与注意事项

1. **单一太阳位置近似**：当前使用全局统一的太阳位置，不支持每个单元独立计算
2. **无地形遮蔽**：不考虑地平线遮挡/阴影效果
3. **直射/散射未分离**：对总下行短波统一乘以几何因子，未对 diffuse 采用天空视域等修正
4. **时区处理**：太阳位置按 UTC 计算（`timezone_hours=0.0`），避免经度推断时区的相位偏移

## 13. 相关文档

- [Model Upgrade As-Is vs To-Be](model_upgrade.md)
- [TSR Physical Assumptions](../TSR_Physical_Assumptions.md)
- [Validation Framework](../validation/tsr/README.md)
