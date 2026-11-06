# MaoyeLoc

**V1.0 · 局地地震绝对定位、双差重定位与震级计算**

MaoyeLoc 用 P、S 波到时估计震源位置和发震时刻，利用近邻地震的同站同相到时差改善相对位置。程序采用 Julia 编写，面向局地地震与微震研究，提供命令行工具 `myloc` 和 Julia 接口。任务使用 TOML 配置，观测与结果使用 CSV。

[下载可执行程序](https://github.com/fumeng6/maoyeloc/releases/latest) · [用户手册](docs/用户手册.md) · [数据约定](docs/数据约定.md) · [数学物理原理](docs/MaoyeLoc原理.md)

## 功能

| 研究任务 | 计算内容 |
| --- | --- |
| 绝对定位 | 根据 P、S 到时估计震源位置与发震时刻，支持固定深度、固定位置和已知事件参考 |
| 双差重定位 | 使用目录拾取差分 CT、外部互相关差分 CC 或二者联合，支持多事件群和联合绝对到时 |
| 定位误差 | 按明确的观测误差和估计条件计算形式误差、条件协方差、置信椭球及完整算法的观测误差传播 |
| 震级计算 | 根据外部测得的振幅、周期或尾波持续时间，计算 ML/MX、Md/Mt，应用仪器校准、站改正并汇总事件震级 |
| 结果检查 | 保存事件状态和资料计数，可开启逐观测残差、权重、射线及迭代诊断，支持批处理与逐读数复核 |

首波走时支持均匀介质、水平常速分层、连续及分段线性梯度模型和固定接收面的走时表。可配置独立 P/S 模型、已知台站延迟，以及按站或按区域的多模型组合。观测由已有的拾取、关联、互相关和振幅测量流程提供。

绝对定位采用逐次线性化、经验加权和 SVD 改正。双差定位采用分阶段经验权与列缩放阻尼，小系统使用密集 SVD，大系统使用 LSQR。公式、权重、参考条件和误差定义见原理文档。

## 下载与安装

从 [Releases](https://github.com/fumeng6/maoyeloc/releases/latest) 下载对应平台的程序包，完整解压到固定目录。可执行包自带 Julia 运行库。

| 平台 | 程序包 | 安装与运行 |
| --- | --- | --- |
| Windows 10/11 x64 | `maoyeloc-v1.0.0-windows-x64.zip` | 双击 `install.cmd`，重新打开终端后运行 `myloc --help` |
| Linux x86_64，glibc 2.34+ | `maoyeloc-v1.0.0-linux-x86_64.tar.gz` | 解压后执行 `sh install.sh`，将 `~/.local/bin` 加入 PATH |

也可直接运行 Windows 包中的 `.\myloc.cmd`、`.\bin\myloc.exe`，或 Linux 包中的 `./bin/myloc`。运行时保留完整的 `bin`、`lib`、`share` 等目录，移动程序目录后重新执行安装脚本。

## 准备输入并运行

在自己的研究目录建立配置文件。

```sh
mkdir my_case
myloc --template my_case/case.toml
```

填写模板中的速度模型、文件路径、坐标、垂直基准与共同 UTC 参考，然后检查并计算。

```sh
myloc --check my_case/case.toml
myloc my_case/case.toml
```

配置中的相对路径按 **配置文件所在目录** 解析。

| `task` | 输入资料 |
| --- | --- |
| `absolute` | 台站、P/S 到时和初始事件位置，可显式指定初始化深度 |
| `relative` | 台站、初始目录和同站同相到时差，可联合绝对到时 |
| `absolute_then_relative` | 先绝对定位，再根据输入到时差或显式 CT 配对条件进行重定位 |

输入和计算统一使用 **km、km/s、s**，局部坐标采用 **x 向东、y 向北、z 向下**。源、台站与速度层界共用垂直基准。`arrival_s` 和 `origin_s` 相对于同一 UTC 参考，`dt_s` 表示完整到时差 **i−j**。事件、台站与读数身份保留为字符串。CSV 字段及单位见数据约定。

常用命令包括显示生效配置、批处理和开启详细诊断。

```sh
myloc --version
myloc --show-config my_case/case.toml
myloc --batch batch.toml
myloc --diagnostics my_case/case.toml
```

## 阅读结果

默认保存定位目录、已请求的误差与运行信息。

| 文件 | 内容 |
| --- | --- |
| `calculation.toml` | 实际生效的模型、单位、时空参考与求解设置 |
| `run_summary.toml` | 各阶段完成情况、事件状态、观测计数与耗时 |
| `absolute/events.csv` 或 `relative/events.csv` | 震源位置、发震时刻、事件状态及所请求的误差字段 |
| 阶段目录中的 `geographic.csv` | 配置地理参考时输出经纬度、深度和 UTC 起时 |
| 已请求的误差与震级表 | 所选误差估计、协方差、震级结果及相关条件 |

使用 `--diagnostics` 或 `[output] diagnostics = true` 开启逐观测与迭代诊断。诊断模式下，`detail = "core"` 保存主要诊断，`detail = "full"` 增加完整历史和谱方向。各结果字段及事件状态见用户手册。

定位误差应结合估计条件阅读。条件误差取决于给定模型、台站与观测误差假设，野外准确度还受速度结构与资料质量影响。双差减小共同传播误差的效果取决于事件距离、台站几何和共同路径条件。

## Julia 运行与绘图

源码运行推荐 **Julia 1.13.2**。定位与震级计算使用 Julia 标准库。在源码目录运行同一套命令入口。

```sh
julia --project=. bin/myloc.jl --check /path/to/case.toml
julia --project=. bin/myloc.jl /path/to/case.toml
```

Julia 脚本可通过 `prepare_case` 和 `run_case` 执行配置，也可向 `travel_time`、`locate`、`locate_catalog` 和 `relocate` 传入内存数据。接口及返回值见用户手册。

可选绘图扩展使用 **CairoMakie**，绘制事件分布、深度剖面和累计事件数。绘图在安装并加载 CairoMakie 的 Julia 源码环境中运行，操作见用户手册的绘图章节。

## 项目内容

| 路径 | 内容 |
| --- | --- |
| `src/` | 走时、定位、震级、误差及数据读写源码 |
| `src/data/` | 通用地区名称与地理范围数据 |
| `ext/` | CairoMakie 绘图扩展 |
| `bin/` | 命令入口、安装、批处理、绘图、结果复核、走时表生成及编译脚本 |
| `docs/用户手册.md` | 安装、配置、运行、结果与 Julia 接口 |
| `docs/数据约定.md` | 输入输出字段、单位、时间、坐标与身份 |
| `docs/MaoyeLoc原理.md` | 方程、数值方法、物理假设与误差含义 |
| `Project.toml` | Julia 包信息与依赖 |

## 开源许可与引用

MaoyeLoc 自身实现与文档采用 **MIT License**，允许使用、修改、再分发和集成，分发时保留版权和许可声明。原作者作品与运行库的权利及条款见来源说明。项目提供 CITATION.cff 供学术引用。

[许可全文](LICENSE) · [来源与第三方说明](NOTICE.md) · [引用信息](CITATION.cff)

## 方法来源与引用

绝对定位采用逐次线性化与加权最小二乘，双差重定位采用 Waldhauser 与 Ellsworth 的到时差方程，稀疏求解采用 Paige 与 Saunders 的 LSQR。具体公式、经验关系及数值选择见原理文档，源码参考和数据表出处集中记于来源说明。研究中请引用实际采用的方法。

1. Klein, F. W. (2002). *User’s Guide to HYPOINVERSE-2000, a Fortran Program to Solve for Earthquake Locations and Magnitudes.* USGS Open-File Report 02-171. [用户指南](https://pubs.usgs.gov/of/2002/0171/)
2. Waldhauser, F., & Ellsworth, W. L. (2000). *A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California.* BSSA, 90(6), 1353–1368. [论文](https://doi.org/10.1785/0120000006)
3. Waldhauser, F. (2001). *hypoDD—A Program to Compute Double-Difference Hypocenter Locations.* USGS Open-File Report 01-113. [用户指南](https://pubs.usgs.gov/of/2001/0113/)
4. Paige, C. C., & Saunders, M. A. (1982). *LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares.* ACM TOMS, 8(1), 43–71. [论文](https://doi.org/10.1145/355984.355989)

项目由 [fumeng6](https://github.com/fumeng6) 维护。使用时请记录 MaoyeLoc 版本和计算配置。问题反馈可附命令、报错和可复现输入，提交至 [Issues](https://github.com/fumeng6/maoyeloc/issues)。
