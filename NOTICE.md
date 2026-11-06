# 来源、致谢与第三方许可

## MaoyeLoc 的许可范围

MaoyeLoc 维护者有权许可的 Julia 实现、脚本和文档采用 MIT License，版权及完整条款见 LICENSE。发行时保留本文件、版权声明和许可文本。

开发过程参考了下列原作者的论文、手册与程序源码，并采用其中的数学方法、经验关系和部分程序规则。各原作者保有其作品的权利，原作品及第三方组件按各自条款使用。MaoyeLoc 的 MIT 声明仅授予本项目有权授予的权利。

## 方法与实现来源

| 来源 | MaoyeLoc 中采用的内容 | 对应实现 |
| --- | --- | --- |
| Fred W. Klein，HYPOINVERSE-2000 及 1.40 扩充手册与源码 | 绝对定位的经验质量、相位、台站、距离和残差行因子，部分 SVD 谱规则、步长和停止控制，形式误差及震级经验关系 | `absolute*.jl`、`magnitude*.jl`、`amplitude_magnitudes.jl`、`coda_magnitudes.jl` |
| Felix Waldhauser 与 William L. Ellsworth，双差方法论文及 HypoDD 1.3 源码与手册 | 双差到时方程、CT/CC 观测定义、分阶段类别与距离/残差权、事件分群、列缩放阻尼和负深度退出重算规则 | `relative*.jl`、`pairing.jl`、`neighbor_pairing.jl` |
| Christopher C. Paige 与 Michael A. Saunders，LSQR | Golub–Kahan 双对角化与迭代最小二乘方法 | `relative_lsqr.jl` |
| USGS HYPOINVERSE 1.40 Geographic NET 表，Fred W. Klein | 地区名称、代码、地理范围及有序地区查询 | `src/data/geographic_regions.toml`、`geographic_regions.jl` |
| Klein 手册和源码所列的震级关系与表 | Eaton、Bakun/Joyner、Richter、Nordquist、CISN/Uhrhammer 距离改正，仪器响应与增益关系，尾波持续时间关系 | `amplitude_magnitudes.jl`、`magnitude_calibration.jl`、`coda_magnitudes.jl` |

表中实现路径均位于 `src/`。震级经验式、离散表及校准规则按所采用的具体关系保留原作者名称。公式适用范围、手册与源码存在的具体差异，以及本程序的数值选择，在数学物理原理和相应函数说明中逐项解释。

源码参考涵盖原程序的算法步骤和经验常数。研究使用时请同时引用所采用的方法与程序来源。MaoyeLoc 的接口、输出精度、误差含义及停止状态由其用户文档定义。

## 上游分发说明

HYPOINVERSE 的软件与手册来源为 USGS 官方发布。USGS 对其工作人员履行公务产生的软件说明了美国境内公共领域地位，并要求保留第三方作品的原有权利。参考数据的来源在 TOML 文件中另有标注。[HYPOINVERSE 软件页](https://www.usgs.gov/software/hypoinverse-earthquake-location) · [USGS 软件许可政策](https://www.usgs.gov/office-of-science-quality-and-integrity/fundamental-science-practices/faq/137-software-and-license)

HypoDD 1.3 的原包 README 以及作者仓库 README 包含以下请求。

> Users of the program are free to make modifications to the programs to meet their particular needs, but are asked not to distribute modified code to others.

原包还请求将改进与错误修正反馈给作者。该请求在此保留，MaoyeLoc 的 MIT 许可不对 HypoDD 原程序作重新许可。MaoyeLoc 的源码包和可执行包分发本项目 Julia 实现，原版 Fortran 源码、程序和观测资料由其作者渠道获取。[HypoDD 作者仓库](https://github.com/fwaldhauser/HypoDD) · [作者软件页](https://www.ldeo.columbia.edu/~felixw/hypoDD.html)

## 可执行包的运行库

可执行包包含 Julia 1.13.2 及其运行库。Julia 自身采用 MIT，链接库和部分标准库使用各自许可。包内 `third_party/LICENSE.md` 和 `third_party/THIRDPARTY.md` 保留 Julia 官方声明，其他可用的标准库及依赖许可保存在 `third_party/`。这些文件与对应组件一同分发。具体运行库名称和版本记录在 `third_party/runtime_components.json`，许可原文的来源记录在 `third_party/license_sources.json`。Julia 及其依赖的完整上游源码可通过 [Julia 1.13.2 完整源码包](https://github.com/JuliaLang/julia/releases/download/v1.13.2/julia-1.13.2-full.tar.gz) 获取。[Julia 官方许可](https://github.com/JuliaLang/julia/blob/v1.13.2/LICENSE.md) · [Julia 第三方说明](https://github.com/JuliaLang/julia/blob/v1.13.2/THIRDPARTY.md)

CairoMakie 为可选绘图依赖，按其自身 MIT 许可使用。绘图扩展在用户安装 CairoMakie 的 Julia 环境中加载。[CairoMakie 所属项目许可](https://github.com/MakieOrg/Makie.jl/blob/master/LICENSE.md)

## 学术引用

引用本程序可使用仓库中的 CITATION.cff，并记录实际使用的版本及计算设置。使用具体定位或震级方法时，再引用相应原始文献。引用是一项学术要求，与 MIT 的版权和许可保留条件分别适用。

1. Klein, F. W. (2002). *User’s Guide to HYPOINVERSE-2000, a Fortran Program to Solve for Earthquake Locations and Magnitudes.* U.S. Geological Survey Open-File Report 02-171. 同时参考原作者随 1.40 发行的扩充指南及源码。[官方出版页](https://pubs.usgs.gov/of/2002/0171/)
2. Waldhauser, F., & Ellsworth, W. L. (2000). *A Double-Difference Earthquake Location Algorithm: Method and Application to the Northern Hayward Fault, California.* Bulletin of the Seismological Society of America, 90(6), 1353–1368.[论文](https://doi.org/10.1785/0120000006)
3. Waldhauser, F. (2001). *hypoDD—A Program to Compute Double-Difference Hypocenter Locations.* U.S. Geological Survey Open-File Report 01-113. 源码规则核对使用 HypoDD 1.3。[官方出版页](https://pubs.usgs.gov/of/2001/0113/)
4. Paige, C. C., & Saunders, M. A. (1982). *LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares.* ACM Transactions on Mathematical Software, 8(1), 43–71.[论文](https://doi.org/10.1145/355984.355989)
