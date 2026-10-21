# MATLAB绘图入口

## 当前论文

- `plot_araim_fde_comparison`：读取2026-10-15三种方案全部启用ARAIM FDE的完整重跑结果，与保留的抗差方案FDE关闭结果比较；默认显示图窗，传入 `true` 同时导出。对比可用率、水平包络率和水平RMS，虚线为纯ARAIM（FDE开启）。详细水平/竖直指标及相同有效历元对照见实验目录的report.md与comparison.csv。

- `run_satellite_figures`：调用原有JC01/JC02/JC04脚本中的卫星数和天空图函数，输出到 `figures/satellites/<测站>/`。每站包含PNG、矢量PDF、可编辑FIG及绘图数据MAT。卫星数沿用原算法：按频点统计流动站RINEX有观测记录的卫星；不是15°截止后或FDE后使用的卫星数。天空图沿用原算法，只绘制GPS/BDS且高度角不低于0°的观测卫星轨迹。
- `plot_satellite_overview`：读取上一入口保存的MAT，输出三站合并的图2卫星数及图3天空图。天空图在相邻记录间隔超过该卫星中位采样间隔两倍时断线，保留所有原始轨迹点；不绘制密集PRN文字，避免遮挡。原脚本及其原图保留不变。
- `plot_current_pl_pe_ksigma`：读取现有九组完整运行结果CSV及相应Ksigma CSV，以红/蓝/绿散点绘制三种方法、三个测站的对数轴包络图，输出到 `figures/current_pl_pe_ksigma/`。不重新计算PL，不重复施加验后尺度，保留不可用历元和全部异常值。

在MATLAB中将当前目录设为本目录，运行上述函数即可；数据目录和输出目录集中在函数开头。

默认只显示并保留图窗，不导出PNG/PDF/FIG。需要同时导出时可传入 `true`，例如 `plot_satellite_overview(true)`。`run_satellite_figures`仍保存供合并图读取的MAT数据及调用用的临时脚本，但不会默认保存图片。

## 原脚本归档

- `result/<原场景名>/`：原result目录内的全部MATLAB脚本，原名及目录分组不变。
- `experiments/<原实验名>/`：原tmp目录内的MATLAB绘图脚本。
- `script_locations.csv`：28个移动脚本的原位置与现位置。

移动后的原脚本只调整路径，仍从原数据目录读取并向原目录输出；统计、配色与筛选逻辑未改变。`plot_paper.m`等旧入口包含早期的Slope/LS-SS对比、告警限值及假固定判据，应作为历史实验入口，不能与当前三方案统计混用。

当前PL/PE/Ksigma原图使用Python生成；新增的MATLAB入口直接读取同一份保存结果，便于后续在MATLAB中调整论文图。

本次已在MATLAB R2023b执行上述三个入口。图2/图3统计覆盖全部2864个输入历元。当前包络图采用Pfa=1e-4、Pmd=1e-3对应结果，所有PE/PL/Ksigma均为散点，对数纵轴；三种方法在同站同方向使用共同坐标范围。PDF矢量导出对密集散点给出性能提示，但运行完成；最终投稿字号和排版需按版面调整。
