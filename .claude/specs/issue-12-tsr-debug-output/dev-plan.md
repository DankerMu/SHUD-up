# TSR 调试/验证输出通道 - Development Plan

## Overview
新增 TSR（Terrain Solar Radiation）调试输出功能，输出水平面太阳辐射（rn_h）、地形校正后辐射（rn_t）和地形因子（factor），用于验证 TSR 计算正确性和区分 TSR 开关影响。

## Task Breakdown

### Task 1: 新增输出文件路径字段
- **ID**: task-1
- **type**: default
- **Description**: 在 FileOut 类中新增 3 个 TSR 输出文件路径字段（ele_rn_h、ele_rn_t、ele_tsr_factor），并在构造函数中初始化路径
- **File Scope**:
  - src/classes/IO.hpp（类定义，约第 75-161 行）
  - src/classes/IO.cpp（构造函数，约第 105 行附近）
- **Dependencies**: None
- **Test Command**:
  ```bash
  make clean && make shud && ./shud ccw
  ```
- **Test Focus**:
  - 编译通过无警告
  - 新增字段在 FileOut 对象中正确初始化
  - 路径字符串格式与现有输出文件一致

### Task 2: 新增数组并在 tReadForcing 中赋值
- **ID**: task-2
- **type**: default
- **Description**: 在 Model_Data 类中新增 rn_h、rn_t、tsr_factor 公有数组指针，在 malloc_Y() 中分配内存，在 FreeData() 中释放，在 tReadForcing() 中计算并赋值（复用已有 TSR 逻辑）
- **File Scope**:
  - src/ModelData/Model_Data.hpp（约第 140 行，公有成员声明区域）
  - src/ModelData/Model_Data.cpp（malloc_Y 约第 85 行，FreeData 相应位置）
  - src/ModelData/MD_readin.cpp（tReadForcing 约第 772 行，TSR 计算逻辑附近）
  - src/ModelData/MD_ET.cpp（ET 函数约第 17 行，访问强制数据区域）
- **Dependencies**: None
- **Test Command**:
  ```bash
  make clean && make shud && ./shud ccw
  ```
- **Test Focus**:
  - 内存分配/释放无泄漏（可用 valgrind 或 Address Sanitizer）
  - tReadForcing 正确计算 rn_h、rn_t、factor（与现有 t_rn 逻辑一致）
  - TSR ON 时 factor 在白天 <1（陡坡）或 >1（平地），夜间 =0
  - TSR OFF 时 factor ≈1

### Task 3: 接入输出系统并记录 TSR 配置
- **ID**: task-3
- **type**: default
- **Description**: 在 initialize_output() 中调用 CS.PCtrl[ip++].Init() 注册 3 个新输出通道，在 Print_Ctrl::open_file() 中扩展 header 记录 TSR on/off 状态（基于 radiation_input_mode 参数）
- **File Scope**:
  - src/ModelData/MD_initialize.cpp（initialize_output 约第 242 行，输出注册区域）
  - src/classes/Model_Control.hpp（Print_Ctrl 类约第 35 行，open_file 签名）
  - src/classes/Model_Control.cpp（open_file 实现约第 375 行，header 构建逻辑）
- **Dependencies**: task-1, task-2
- **Test Command**:
  ```bash
  ./shud ccw
  ```
- **Test Focus**:
  - 输出目录中出现 `*.rn_h.dat`、`*.rn_t.dat`、`*.tsr_factor.dat` 文件
  - 二进制 header（前 1024 字节）包含 TSR/radiation/lon/lat 配置信息
  - CSV header（如果启用 ASCII 输出）同样包含配置注释
  - 数据行数 = 模拟步数 * 输出频率

### Task 4: 集成验证
- **ID**: task-4
- **type**: default
- **Description**: 端到端验证输出数据物理合理性：检查 rn_h 昼夜模式、factor 夜间归零、rn_t = rn_h * factor 关系、TSR ON/OFF 差异
- **File Scope**:
  - 无源码修改，仅验证输出文件
- **Dependencies**: task-3
- **Test Command**:
  ```bash
  ./shud ccw
  # 建议补充 Python/R 脚本验证输出数据
  ```
- **Test Focus**:
  - rn_h: 白天 >0 (典型 200-1000 W/m²)，夜间 =0
  - factor: 夜间 =0，白天 TSR ON 时 0.5-1.5 范围（取决于地形），TSR OFF 时 ≈1.0
  - rn_t = rn_h * factor（数值误差 <1e-6）
  - TSR ON 与 OFF 对比：rn_t 在坡度大的网格差异显著
  - Header 正确反映 TSR 配置（SWNET/SWDOWN、lon/lat mode、坐标值）

## Acceptance Criteria
- [ ] 运行 `./shud ccw` 后 output 目录下存在 rn_h、rn_t、tsr_factor 三个输出文件
- [ ] rn_h 符合太阳辐射日循环模式（白天 >0，夜间 =0）
- [ ] factor 在夜间为 0，白天 TSR ON 时体现地形影响（非恒定 1），TSR OFF 时约为 1
- [ ] 满足物理关系 rn_t = rn_h * factor（相对误差 <0.1%）
- [ ] 二进制文件 header 包含完整的 TSR、radiation mode、lon/lat 配置信息
- [ ] 通过 TSR ON/OFF 对比实验，验证输出能区分两种模式的辐射差异
- [ ] 编译通过无警告，运行无内存泄漏，代码覆盖率 ≥90%

## Technical Notes
- **输出系统集成**: 复用现有 Print_Ctrl 框架，保持与 ele_q_prcp、ele_q_ETP 等输出通道一致的调用模式
- **数组生命周期**: rn_h/rn_t/factor 数组遵循 Model_Data 其他时变数组的管理模式（malloc_Y 分配，FreeData 释放）
- **TSR 逻辑复用**: tReadForcing() 中已有 TSR 计算逻辑（调用 SolarRadiation 模块），新任务仅需保存中间结果到新增数组
- **Header 扩展性**: 当前 header 固定 1024 字节，新增 TSR 状态字段需在现有空间内追加（已有 radiation mode、lon/lat 信息，补充 TSR on/off 标志）
- **性能影响**: 每个输出通道增加约 NumEle * sizeof(double) 字节内存，对于千元素级别网格（如 ccw 示例）影响可忽略
- **向后兼容**: 新增输出为可选（通过 cfgout 文件控制），不影响现有输出通道的行为
