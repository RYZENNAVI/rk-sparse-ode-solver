# 基于 Runge-Kutta 方法的稀疏 ODE 系统高效求解器

本项目配套论文 **《Efficiently Solving Sparse Linear ODE Systems with Explicit and Embedded Runge-Kutta Methods: Implementation, Parallelization and Optimization》**，提供了文中算法与实验的完整实现。

---

## 📖 项目概览

本代码库提供了一种结构化方法，用于求解**稀疏线性常微分方程组（ODE 系统）**，聚焦以下几个方面：

- **显式 Runge-Kutta 方法**（如 RK4）
- **嵌入式 Runge-Kutta 方法**（如 RK45 类型）
- 使用 **OpenMP** 进行并行化
- 各类**性能优化技术**：缓存调优、稀疏存储格式、数值精度控制等

这些方法被应用于稀疏结构的 **TASEP（全不对称简单排除过程）模型**，并进行了系统评估。

---

## 📐 核心模块

### 1. 数学基础

项目包括：

- 经典 Runge-Kutta 方法（RK4）
- 自适应步长的嵌入式 Runge-Kutta 方法
- 误差控制与步长调整策略

### 2. 实现方式

实现了三个版本：

- 两个用于 **显式 RK 方法**
- 一个用于 **嵌入式 RK 方法**

全部使用 **C 语言** 编写，注重 **模块化与数值稳定性**。

### 3. 并行化策略

所有方法均基于 **OpenMP** 实现了并行化：

- 针对行的并行计算方式，避免依赖冲突
- 通过调度策略精确控制线程行为

### 4. 性能优化

优化内容包括：

- **缓存优化**（如循环分块、数据局部性提升）
- **稀疏矩阵存储格式**（CRS、对角存储及改进版）
- **算法级优化**，权衡计算精度与性能

---

## 🚀 快速开始

### 🔧 环境要求

- 支持 OpenMP 的 GCC 或 Clang 编译器
- 可选：CMake 构建工具
- 推荐平台：Linux / macOS / Windows + WSL

### 📦 编译与运行

```bash
make run
# 或运行性能基准测试
make benchmark
```

---

## 📫 联系方式

如有问题或建议，欢迎提交 Issue 或发送邮件至：

```
liyixianggermany@163.com
```

---

## 📝 许可证

MIT License © Yixiang Li
