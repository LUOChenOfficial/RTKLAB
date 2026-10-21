# RTK 完好性代码导航

本目录保留 RTKLIB 的 C 接口、列主序矩阵和现有输出格式。新增实验代码按职责组织：

| 文件 | 职责 |
| --- | --- |
| `src/rtk_experiment_config.h` | 实验主开关、方法选择及编译参数默认值；可用编译器定义覆盖 |
| `src/rtk_int_config.c` | 将所选实验配置应用到 `prcopt_t`；入口为 `rtkimcfg()` |
| `src/rtk_int.c` | 监测生命周期、FILTER_SS/LS_SS、FDE 判定、偏差表和监测结果输出 |
| `src/rtk_gls.c` | 当前历元 GLS 解分离及故障模式评估 |
| `src/rtk_slope.c` | 斜率法统计阈值、先验捕获及保护水平计算 |
| `src/ransac.c` | 载波单差一致集选择内核 |
| `src/rtk_robust.c` | 固定先验的 IGGIII 权值迭代与等效观测协方差 |
| `src/rtk_debug.c` | 调试文件、观测与模糊度诊断、在线蒙特卡洛 |
| `src/rtk_int_common.c` | 共用概率阈值与坐标工具 |
| `src/rtk_int_internal.h` | 监测模块内部接口；不增加到公共 `rtklib.h` |

`main.c` 在设置普通定位参数后调用 `rtkimcfg()`，再设置后续参数。配置函数会根据当前宏覆盖 AR、迭代次数、FDE 和风险参数；调用顺序与原代码一致。

`rtkpos.c` 保留必要的数据捕获与监测调用：保存 DD 模型，在 KF 更新前捕获先验，定位后检测并评估。已有 FDE 重算流程仍在定位入口中，保持原有状态生命周期。

后续新增方法应将数值计算放入对应方法文件，通过 `rtkim_t` 传递缓存与结果；调试文件句柄由 `rtk_debug.c` 管理，内部诊断使用 `rtkimdbg()`。保持矩阵维数、单位、调用阶段和失效返回语义明确，避免在方法文件中加入数据集路径或文件格式逻辑。

新增代码的命名和排版约定见 [STYLE.md](STYLE.md)。

工程使用 `rtklibdemo5.vcxproj` 的 `Debug|Win32` 配置。结构整理的回归记录位于仓库 `tmp/style_refactor/`，包括整理前后 GLS、SLOPE 和关闭组合开关三种配置的短时真实 Static 数据输出。短时一致性检查不等同于全数据集完好性验证。
