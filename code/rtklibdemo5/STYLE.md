# 新增 RTK 模块编码约定

以 RTKLIB 的 `lambda.c`、`rtkcmn.c` 等原有模块为参照。当前仓库的历史代码本身存在排版差异；新增完好性模块采用以下统一约定，不重排整个上游代码库。

- 函数采用简短小写名称：公共监测接口为 `rtkiminit/rtkimfree/rtkimdetect/rtkimassess`，方法接口为 `rtkglsupdate/rtkslopethres/rtkransac`，诊断接口使用 `rtkdbg` 前缀。文件内工具使用 `calcpl/lsqdd/chkgeom` 等含义明确的短名称。
- 类型保留 `_t` 后缀，如 `rtkim_t`、`rtkim_pl_t`。宏保留大写及下划线，并保留当前实验主开关名称。
- 变量采用 RTKLIB 常见的 `rtk/opt/mon/sat/ix/nx/nv/np`。矩阵用 `H/R/P/Q/F/K` 或带方法前缀的 `slK/glsQ`，状态与参数采用 `intmon/intpfa/glsready` 等短名称。
- 缩进四空格；函数体左花括号另起一行，控制语句左花括号同行；指针星号靠变量名。表达式采用紧凑写法，函数参数逗号后保留空格。长接口按参数换行并对齐，复杂循环保持层次清楚。
- 文件头说明文件职责；公共数值接口采用 `args / return / notes` 注释，注明输入输出、列主序、单位和维数。缩短变量名时用注释保留科学含义。
- 日志记录名、文件扩展名和配置文本键保持原值；代码符号重命名不改变输出协议。

| 旧代码名称 | 当前名称 | 含义 |
| --- | --- | --- |
| `rtkim_init` | `rtkiminit` | 初始化监测状态 |
| `rtkim_configure` | `rtkimcfg` | 应用实验处理参数 |
| `rtk_gls_position_update` | `rtkglsupdate` | 投影完整状态 KF 的位置更新 |
| `rtk_slope_thresholds` | `rtkslopethres` | 计算斜率检验阈值与非中心参数 |
| `enable_rtk_integrity_monitor` | `intmon` | 监测开关 |
| `rtk_integrity_false_alarm_prob` | `intpfa` | 虚警概率 |
| `rtk_integrity_miss_detect_prob` | `intpmd` | 漏检概率 |
| `slope_K / slope_Qinv` | `slK / slQi` | 位置增益与新息精度矩阵 |

旧的自定义函数和字段已更名，仓库内调用与测试同步更新；这些名称属于本项目新增接口。原有 RTKLIB 接口保持原名。
