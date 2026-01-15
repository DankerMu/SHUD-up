# TSR Validation Scripts - Development Plan

## Overview
创建自动化验证脚本，对比 TSR 开关前后的模型输出，确保 TSR 特性可复现且输出隔离。

## Task Breakdown

### Task 1: 搭建验证目录结构
- **ID**: task-1
- **type**: quick-fix
- **Description**: 创建 validation/tsr/ 目录结构，配置 .gitignore（排除 tmp/ 和 output/），编写 README.md 说明脚本用途和使用方法
- **File Scope**: validation/tsr/README.md, validation/tsr/.gitignore
- **Dependencies**: None
- **Test Command**: `test -d validation/tsr && test -f validation/tsr/README.md && test -f validation/tsr/.gitignore && echo "Directory structure validated"`
- **Test Focus**:
  - 验证目录创建成功
  - README.md 包含脚本使用说明
  - .gitignore 正确排除 tmp/ 和 output/ 目录

### Task 2: 实现 TSR 基准和对比运行脚本
- **ID**: task-2
- **type**: default
- **Description**: 实现三个脚本：
  - `_common.sh`：公共函数（错误检查、输入隔离复制、配置修改、输出备份）
  - `run_baseline.sh`：TSR=OFF，输出到 output/ccw.base/
  - `run_tsr.sh`：TSR=ON，输出到 output/ccw.tsr/
  每个脚本需包含：二进制存在性检查、输入目录检查、工作目录设置、cfg.para 的 TERRAIN_RADIATION 配置 upsert、已有输出目录自动备份
- **File Scope**: validation/tsr/_common.sh, validation/tsr/run_baseline.sh, validation/tsr/run_tsr.sh
- **Dependencies**: Depends on task-1
- **Test Command**: `bash validation/tsr/run_baseline.sh && test -d validation/tsr/output/ccw.base && bash validation/tsr/run_tsr.sh && test -d validation/tsr/output/ccw.tsr && echo "Scripts executed successfully"`
- **Test Focus**:
  - 缺少二进制文件时报错退出
  - 缺少输入目录时报错退出
  - 成功复制输入到 validation/tsr/tmp/
  - cfg.para 中 TERRAIN_RADIATION 正确设置为 0（baseline）或 1（tsr）
  - 输出目录已存在时自动备份为 .bak{timestamp}
  - 两次运行输出到不同目录且互不覆盖

### Task 3: 实现集成测试脚本
- **ID**: task-3
- **type**: default
- **Description**: 实现 integration_test.sh，执行以下深度验证：
  - 调用 run_baseline.sh 和 run_tsr.sh
  - 检查两个输出目录存在且非空
  - 解析关键 .dat 文件（如 ccw.surf、ccw.rad）的 header（行数、列数、时间步）
  - 采样验证数值范围合理性（非 NaN、非全零）
  - 对比两组输出的差异（TSR=ON 时辐射应有差异）
  - 生成验证报告到 validation/tsr/test_report.txt
- **File Scope**: validation/tsr/integration_test.sh
- **Dependencies**: Depends on task-2
- **Test Command**: `bash validation/tsr/integration_test.sh && grep -q "PASS" validation/tsr/test_report.txt && echo "Integration test passed"`
- **Test Focus**:
  - 两个脚本顺利执行无错误
  - 输出文件存在且文件大小 > 0
  - .dat 文件 header 解析正确（行数、列数、时间步匹配预期）
  - 数值采样验证通过（无 NaN、非全零、范围合理）
  - TSR=ON 与 TSR=OFF 的辐射输出存在数值差异
  - 测试报告清晰标注 PASS/FAIL 和失败原因

## Acceptance Criteria
- [ ] 用户执行一条命令可复现 TSR 开关对比实验
- [ ] TSR=OFF 输出到 output/ccw.base/，TSR=ON 输出到 output/ccw.tsr/，互不覆盖
- [ ] 脚本在以下场景给出清晰错误信息：缺少二进制、缺少输入目录、缺少工作目录权限
- [ ] 集成测试深度验证输出文件内容（header 解析、数值采样、差异对比）
- [ ] 所有脚本通过集成测试且测试报告显示 PASS
- [ ] README.md 提供完整的使用说明和故障排查指南

## Technical Notes
- **输入隔离策略**：每次运行前复制 input/ccw/ 到 validation/tsr/tmp/ccw，避免污染原始输入
- **配置修改策略**：使用 sed/awk upsert cfg.para 中的 TERRAIN_RADIATION 行，保留其他配置不变
- **输出备份策略**：检测到已有输出目录时，重命名为 `.bak$(date +%Y%m%d_%H%M%S)` 后再运行
- **深度验证策略**：不仅检查文件存在，还需解析文件 header（行数、列数、时间步）和采样验证数值范围，确保输出有效性
- **错误处理原则**：所有错误场景必须打印可操作的错误信息（如"缺少二进制，请运行 make 编译"），并使用非零退出码
- **可移植性**：脚本应使用 POSIX 兼容的 shell 语法，避免依赖 bash 特定特性（除非明确声明 #!/bin/bash）
