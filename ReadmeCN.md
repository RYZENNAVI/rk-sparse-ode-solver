# rk-sparse-ode-solver

该项目提供了一组 C 语言工具以及基于 Runge–Kutta 的求解器，用于处理来自 TASEP 模型的稀疏线性常微分方程组。工具可以生成矩阵和向量，在不同稀疏格式间转换，并执行显式或嵌入式 Runge–Kutta 方法。求解器采用 OpenMP 并行，并包含多种优化特性。

## 工具

以下辅助程序位于 `src` 目录中（详细用法见各源文件开头的注释）：

- **genmat** – 生成 COO 文本格式的稀疏矩阵，例如：
  ```bash
  ./genmat -n 4 -i zahlen.txt -o test_n_4.mat
  ```
- **gencrs** – 将 COO 格式的矩阵转换为 CRS：
  ```bash
  ./gencrs -i test_n_4.mat -o test_n_4.crs
  ```
- **gencrs_new** – 将 COO 格式转换为改进的 CRS：
  ```bash
  ./gencrs_new -i test_n_4.mat -o test_n_4_new.crs
  ```
- **gendia** – 将 COO 格式转换为对角存储：
  ```bash
  ./gendia -i test_n_4.mat -o test_n_4_dia.mat
  ```
- **geninit** – 生成具有正态分布的初始向量文件：
  ```bash
  ./geninit -d 4 -n test_n_4.init
  ```
- **vis** – 生成矩阵的文本可视化：
  ```bash
  ./vis -i test_n_4.mat -o test_n_4_vis_new.mat
  ```

## `rk_new` 求解器

`rk_new` 使用多种 Runge–Kutta 变体求解稀疏线性 ODE 系统，主要特性如下（完整用法见 `rk_new.c`）：

- 内置表格的标准 RK4
- 使用 MAT 或 CRS Butcher 表的显式方法（如 RK4 或 E5）
- 显式方法的 Richardson 外推
- 带容差控制的嵌入式方法（RKF45、BS23）
- 支持 CRS、改进 CRS 与对角矩阵格式
- 可配置 CPU 亲和性的多线程
- 两种版本的循环分块优化，可设定块大小

RK4 加载 MAT 表的一般调用方式：
```bash
./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
```

更多示例请参见 `rk_new.c` 中的注释。

### 构建

进入 `src` 目录运行 `make` 即可编译全部工具，生成的二进制文件会在同一目录下。

```bash
cd src
make
```

## 许可

本项目仅用于科研与教学目的。版权信息请参阅各源文件顶部的声明。
