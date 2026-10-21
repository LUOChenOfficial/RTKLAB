# 单滤波器 RANSAC-ARAIM 试验

## 2026-10-15 经重算验证的FDE恢复

当显式启用`RTK_GLS_FDE=1`时，三个新开关默认均为1：
`RTK_GLS_FDE_FLOAT_FIRST`、`RTK_GLS_FDE_TRANSACTIONAL`、`RTK_GLS_FDE_SEARCH`。
前者在浮点全集/子集检验通过、固定相关检验报警时，优先用同一历元前状态
关闭AR重算；只有通过完整监测，才接受浮点输出。该操作不是错误固定的确诊。
联合检验仍参与故障监测和PL，风险预算、验后尺度及PL公式不变。

卫星排除候选按单星优先、再双星排序，同阶按原始浮点解分离比值降序。
初始几何及数值模型有效的候选均可试算，统计量只用于排序，不直接声明故障。
每次试算均恢复`relpos`前完整状态，不重复SPP、不重复施加历元时间更新。
剔星重算后若仅固定相关检验失败，再试该观测集的浮点恢复。
只有位置状态有效、模型就绪、全部假设有效、无报警且HPL/VPL有限才提交后验。
监测仍可枚举最多8星，但自动恢复搜索仅支持单星和双星。

所有候选失败时，丢弃试算造成的状态及跟踪计数变化，恢复原始全观测定位后验
继续递推；当前历元明确保持PL=inf。该原始后验仍未通过完好性验证，不能视为
安全恢复，也没有消除历史故障风险。恢复只拷贝监测标量，保留当前动态内存所有权；
被试算覆盖的模型缓存标记为未就绪，当前诊断已完成，不再使用该缓存重新评估。

日志`GLS_FDE_TRY_SAT`记录每次尝试的卫星，`GLS_FDE_TRIAL`记录完整验证结果，
`GLS_FDE_RESULT`记录历元尝试总数、接受/回退和初始动作。`.pint`原FDEact字段
仍表示尝试动作，不能作为已验证剔除数；接受的最终候选须结合上述日志读取。
多次试算的`.glsmodes`可含同历元多组记录，不能当作仅有最终评估的一组假设。

将三个开关全部设为0可复现旧单次重算流程。分阶段全时段试验、统计和测试入口：
`tmp/gls_fde_verified_20261015/{run.ps1,analyze.py,test.ps1}`。
`gls_fde_test.c`直接验证深拷贝状态回退、监测缓存所有权、失败/缺失模式不可接受
以及候选排序；`gls_test.c`同时覆盖固定相关报警的浮点恢复选择。
可用率提高不等价于固定解性能或统计完好性保证提高，需同时报告解状态、PE及包络率。

## 2026-10-13 高阶纯 GLS / FDE 试验

`RTK_GLS_MAX_FAULTS` 控制枚举阶数，允许 1..8，默认 2。枚举所有 1..K 星删除组合，
按整颗卫星删除所有观测；参考星删除仍使用完整双差变换。`RTK_GLS_FDE=1`
打开 GLS 试验性恢复，默认 0；纯 GLS 编译同时设置 `Open_Robust=0` 和 RANSAC 关闭。
固定输出继续包含错误固定及观测故障与错误固定联合模式。

FDE 候选采用最小报警删除阶数，同阶选择最大解分离比值。这是探索性候选规则，
不是已证明的故障隔离检验。只重算一次：观测报警删除候选集合，恢复本历元
`relpos` 前保存的完整估计状态、协方差、卫星/模糊度历史，再用保留观测重算；
GLS 分支不再调用 `rtkfree/rtkinit`，不重复 SPP，也不从已更新的后验再次更新。
仅错误固定报警则关闭 AR 重算本历元。重算后重新完整评估，模式缺失/报警保持 PL=inf。
删除集合仅用于本历元，后续历元恢复原 AR 配置；共同先验污染未由此得到完整处理。
预测、周跳/失锁计数及 AR 从保存历史重放一次；监测器动态内存保持有效并以新观测重建。
失败重算仍不可用，不回退到被拒绝的全观测后验。保存状态分配失败时禁用本历元重算。
`.rtk_debug` 的 `GLS_FDE_PRIOR` 记录不重置恢复的状态维数和历元时间差。

高阶输出模式 6 为三颗及以上观测故障，7 为其与错误固定联合模式。
`.glsdiag` 的 Hsat1/Hsat2 仅记录集合前两星；完整控制包络的集合由
`.rtk_debug` 的 `GLS_HIGH_SOURCE` 记录，剔除集合由 `GLS_FDE_EXCLUDE` 记录。
离线风险加权脚本不支持模式 6/7，不能直接用于高阶输出。
短段真实数据、结果和复现脚本位于 `tmp/pure_gls_highorder_20261013/`；不作为全天验证。

`ENABLE_RTK_RANSAC_ARAIM=1` 默认选择 `RTK_INT_METHOD_GLS_SS`（3）。
主滤波器筛选后更新的完整浮点状态和协方差传到下一历元；连续AR只生成
当前固定候选输出，不执行fix-and-hold。各子集不创建私有`rtk_t`。
显式设置`RTK_INT_METHOD=RTK_INT_METHOD_FILTER_SS`（0）仍可运行原并行滤波分支。
纯GLS对照：组合和RANSAC关闭，显式设置`RTK_INT_METHOD=3`，监测保持开启。

## 当前历元子集

在实际`filter()`前，按相同有效状态索引构建完整创新协方差
`Q=H' P_prior H+R`及位置块`F=(P_prior H)[0:3,:]`。
子集只对创新模型施加稀疏变换`T`，使用`T Q T'`、`F T'`及`T v`更新位置块。
这是完整状态KF位置后验的等价计算，保留位置—模糊度相关项，不是三维位置LS。
所有子集共用更新前先验；主浮点解与重算全观测GLS的位置差超过0.1mm时不可用。

单星和双星故障都按卫星删除全部观测。参考星被删除时，选择同双差组的剩余
目标星作锚点：`(ref-target)-(ref-anchor)=anchor-target`；观测和完整协方差
同时变换。该过程不重新调用星历、零差、周跳及整个RTK解算。
使用缩放Cholesky求解，避免显式求逆导致近低秩分离协方差失去半正定性。

若有n颗实际参与双差的卫星，当前解算次数为`1+n+n*(n-1)/2`。
不维护子滤波器，不受原`MAXSAT+1`子滤波器数组上限截断。
每个子集要求当前伪距位置几何满秩；不能只凭先验产生有限协方差就判作可监测。

## 错误固定与保护水平

固定输出增加全观测浮点保护子解；所有删星保护子解也保持浮点，不继承主整数。
同一次删星浮点解复用于观测故障和观测故障+错误固定联合模式，分别计入风险分配。
固定输出的故障假设数为`2*(n+n*(n-1)/2)+1`；浮点输出不增加错误固定模式。
模式编码：0名义、1单星、2错误固定、3双星、4单星+错误固定、5双星+错误固定。

2026-10-13 修正：每轴采用解分离检验阈值加子解误差界：

`PL_axis_i=K_FA*sigma_SS_axis+K_MD*sigma_subset_axis`。

`sigma_SS=sqrt(diag(P_subset-P_output))`（转至 ENU 后）；浮点输出使用主浮点
协方差，固定输出使用当前整数候选的条件固定协方差。实际位置分离量只用于检验
和日志，不再参与 PL 第一项。差协方差不满足数值半正定检查时模式不可用。
固定输出增加固定—删星浮点子解联合检验，与其 PL 阈值一致；虚警预算按包括联合
模式的全部假设数分配。纯浮点名义模式解分离标准差为零。
`.pldiag` 的 Hsss/Vsss 保存主导模式解分离标准差，Hsep/Vsep 仍是实际分离量。
K_FA/K_MD 沿用总预算等分方式，不等于已实现 ISM 故障先验风险求和模型。

每模式HPL由其E/N轴界的二范数得到，VPL取U轴界；最终取模式最大值。
子解高斯双侧尾部预算按三轴及全部假设分配。
观测故障检测用主浮点/子浮点分离和`P_subset-P_main_float`；正确整数名义模式下
的固定/浮点检测用`P_main_float-P_fixed`。这些差协方差需通过数值半正定检查。
检测与自动FDE分离：检测超限输出inf，不触发主滤波器重新解算。
缺失必需子集、模型异常、未启用固定保护而输出固定解时也输出inf。

没有从PE调整PL，也未加载原经验rbias偏差项。本轮不是完整标准ARAIM风险模型；
RANSAC选择分布、整数选择效应、相关/高阶故障和历史先验偏差仍待验证。
解除当前整数约束不清除历史先验污染，因此本方法禁止fix-and-hold。

## 首轮复现

`tmp/ransac_gls_20261010/build.ps1`构建RANSAC开/关Release版本并运行`gls_test.c`。
`run.ps1`运行Open/Static/Wall/Tree/Build首300秒和JC01/JC02/JC04首3000秒，AR开启且非hold；每种设置
再关闭监测运行一次，检查监测没有改变位置输出。
`analyze.py`生成逐历元数据、统计和报告；Open/Wall/Tree/Build用已有ref.txt，JC组用已有
GEOLIB相对坐标，未当作独立绝对真值。Static只匹配质量3的KSXT参考轨迹，时差不超过0.05秒。
覆盖超限仅在有参考且有限PL的历元统计；同时报告计划历元、定位输出和有限PL数。
`.glsdiag`记录假设总数/有效数、缺失数、告警、主解一致性、耗时和PL来源两颗卫星。
`.pint`的subA在GLS模式下表示有效假设数；原并行分支格式保持不变。

测试覆盖完整/单双星（含参考星）子集与实际KF的位置/协方差比较、错误固定告警、
禁用固定保护以及缺失模型时不可用。真实数据结果不预设包络或分米级目标成立。
# Deformation-monitoring risk profile (2026-10-13)

Epoch covariance scale (2026-10-14): `RTK_GLS_POST_SCALE=1` estimates the
common augmented-GLS variance scale as `v_pre' S^-1 v_pre / nv`, using the
full active-state prior and current effective measurement covariance.
All output position/velocity and GLS main/subset/fixed covariances use
the same scale; nominal internal P/Pa remain unscaled. The scale may be
below one. This is a reporting/common-scale approximation; Gaussian
K multipliers do not account for estimated-scale uncertainty or selection.
Use `RTK_GLS_POST_SCALE=0` to restore nominal reporting. POST_SCALE logs
contain time, redundancy, total quadratic, variance factor and sigma0.

The default total GLS experimental budgets are now `RTK_GLS_PFA=1e-4`
and `RTK_GLS_PMD=1e-3`, before equal mode/axis allocation. Compiler
definitions can restore the previous budgets (`1e-6`, `5e-8`).
This is a changed application risk setting, not posterior risk calibration.
Increasing PFA lowers solution-separation thresholds and can reduce PL
availability even while increasing PMD tightens candidate protection levels.
The slope profile retains its own existing budgets.
