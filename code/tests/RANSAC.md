# RANSAC 迁移说明

来源：`F:\完好性\RANSAC_RAIM-RTK\RANSAC-RTK\RANSACRTK.cpp`，
以及其 `rtkpos.cpp` 中的 `sdres()`、`relpos()` 调用链。
目标实现为纯 C：`src/ransac.c`，适配入口为 `rtkpos.c` 的 `ransac_sd()`。

在 `src/rtklib.h` 用两个互斥的总开关控制组合方式，修改后需要重新编译：

|方式|ENABLE_RTK_RANSAC_ARAIM|ENABLE_RTK_RANSAC_SLOPE|
|---|---:|---:|
|RANSAC + 当前历元GLS解分离（当前）|1|0|
|RANSAC + 创新量斜率保护级|0|1|
|原有 main 配置：整数固定开启、监测关闭|0|0|

组合开启时自动启用 RANSAC、选择对应完好性算法，在 main 中启用监测，
使用一次 KF 更新。斜率组合关闭整数固定；ARAIM组合默认使用GLS_SS，
允许连续AR但不hold，传递主浮点KF状态，监测最多两星故障及错误固定。
GLS子集不维护跨历元滤波器，详见`GLS.md`。
显式选择FILTER_SS仍可运行下面记录的原并行单星浮点试验。
两者关闭自动FDE恢复，PFA为 `1e-6`、PMD预算为 `5e-8`。
两种组合同时启用、组合开启时显式关闭RANSAC或选择不匹配算法，均报编译错误。
两种组合关闭时，还可独立将 `ENABLE_RTK_RANSAC` 设为 `1` 使用观测筛选。
其他参数为 `RTK_RANSAC_MAX_ITER=1000`、
`RTK_RANSAC_SIGMA=12.0`、`RTK_RANSAC_PROBABILITY=0.99`，对应原程序
`RANSACRTK::multiThreshold=12.0`、`probability_=0.99`。
`RTK_RANSAC_SIGMA` 是阈值倍数；观测标准差由 `sqrt(varerr(...))` 取得，单位为米。
原程序另有一个未使用的固定阈值成员 `threshold=30`；RTK默认 `thresholdMethod=1`，
实际使用的是 `12*sqrt(1/P[i])`，不是固定30米，也不是3sigma。

流程：预测状态和零差残差 → 载波单差模型 → RANSAC 一致集 →
双差参考星选择与双差构建 → KF → 原有 PAR/固定解验证。
掩码每历元重置，同一历元滤波迭代及浮点/固定后验残差复用同一掩码。
完整双差协方差仍由原有 `ddcov()` 构建，共享参考星相关项不变。

单差模型为 `v_s=(y_rover-y_base)-lambda_s*N_float_s`，设计列为
`[-e_x,-e_y,-e_z,lambda_s*I(system)]`。
随机取模型维数个观测做加权最小二乘，按
`|v_s-A_s*x| < RTK_RANSAC_SIGMA*sqrt(var_s)` 判断一致性。
原程序将所有载波频率放在一次 `calRANSAC()` 中，使用位置列和各系统公共偏移列；
当前 C 版已改为联合多频调用。此前按频率分别调用的结果不代表本次修正后的表现。
优先选择一致观测最多的模型，同数时选择归一化残差均方最小的模型。
模型中的位置与公共偏移仅用于观测筛选，不直接写入滤波状态。

迁移中的明确适配：

- 源工程实际只构建载波单差，目标也只筛选载波，不剔除伪距。
- 所有频率联合拟合，系统偏移列按实际出现的系统分配，各行使用本频率波长。
- 按卫星和频率索引掩码，在参考星选择前排除粗差。
- 修正未初始化控制量、跨调用迭代计数、随机不可复现、无界尝试及源权重映射越界问题。
- 源默认硬阈值和最大一致集保留；未迁移可选 IGGIII、相关性筛选、局部优化、DOP 和 PL 实验分支。
- 不迁移 `RANSACStatic` 的静态跨历元坐标反算流程。
- 无冗余、秩亏、非法方差或无有效一致集时保留原观测并记录未应用。
- 当前支持 GPS/BDS/Galileo/QZSS 原始频率载波；GLONASS/SBAS 不参与筛选。
  DGPS、IFLC、估计电离层或估计对流层模式跳过 RANSAC。

启用后 `.rtk_debug` 中新增：

`RANSAC_SD time sat freq model_dim observation_count consensus_count keep sd_m variance_m2`

`keep=0` 表示此卫星此频率的载波观测被排除。
`consensus_count=0` 表示未找到有冗余的一致集，保持原观测。
此筛选不提供独立完好性保证，也不能据此判断整数固定一定正确。

验证入口：在仓库根目录运行 `& code/tests/run_ransac_tests.ps1`。
脚本以 `e16b799` 为基线，分别构建原版、关闭、独立 RANSAC、组合斜率四种版本，
运行 Static 数据 2025/04/25 07:47:00–07:48:00（1 Hz），
比较关闭时逐行解算结果，并执行合成粗差/秩亏/无冗余/可复现性测试。
测试关闭在线 Monte Carlo 以缩短耗时；不改动正常配置。
输出与临时基线副本仅写入 `tmp/ransac_validation/`。
真实数据短时运行用于接口验证，不构成定位精度或错误固定率改进的结论。

本次验证：关闭时与基线 61 行解算结果完全一致；开启时输出 61 个历元，
1,125 条单差记录中剔除 355 条；24 个固定双差端点没有使用被剔除的载波。
该短时段关闭/开启时固定历元分别为 5/4，浮点历元为 56/57，
因此不能把这次迁移解释为固定率提升。正式精度评估还需完整数据与真值对照。

## 组合斜率实验

`rtkim_slope_prior()` 在实际 KF 更新之前按 `filter()` 的相同有效状态规则
保存位置增益 `G` 和创新精度矩阵 `Qinv`。完整模型为
`Q=H' P_prior H+R`，`G=(P_prior H Qinv)` 的前三行。
保留完整位置—模糊度相关项与双差观测相关项，不使用滤波后的位置协方差重算增益。

创新检测量为 `T=v'Qinv v`。在当前高斯随机先验模型下使用 `nv` 个自由度，
不能套用三维无先验最小二乘的 `nv-3`。
阈值由 `P(T>threshold|no fault)=pfa` 求解；非中心参数由
`P(T<=threshold|lambda)=pmd_mode` 求解。
单故障模式和名义噪声按 `pmd/(mode_count+1)` 分配预算。

每个保留的卫星、频率、观测类型定义一个单差故障方向 `c`：
对应参考星的双差行为 `+1`、对应目标星的行为 `-1`。
参考星故障因此同时影响整个双差组。
水平/垂直斜率分别为 `norm((E*G*c)[E,N])/sqrt(c'Qinv*c)` 和
`abs((E*G*c)[U])/sqrt(c'Qinv*c)`。
故障项取各自最大斜率乘 `sqrt(lambda)`，再叠加分别计算的名义噪声项。
最终 HPL/VPL 写入原有 `.pint` / `.pld`；`.rtk_debug` 新增 `SLOPE_KF` 诊断行。
检测超限、无有效模型、非浮点输出或多次迭代时输出 `inf`，表示不可用。
组合模式不触发 FDE，也不计算错误整数固定的保护水平。

该模型只覆盖当前历元单个单差观测故障，在固定筛选模型、高斯先验假设下计算。
RANSAC 的数据选择效应、跨历元故障传播、多故障和错误固定风险尚未纳入，
输出属于有待完整故障注入验证的实验性保护水平，不能视为已证明的端到端完好性保证。

组合版本验证：61 个历元全部为浮点解，61/61 输出有效保护水平。
合成测试对照二维正态数值积分的非中心卡方参考值，核对完整状态相关项、
有效状态筛选规则及检测超限/模型缺失时不可用标志；
在实际 `filter()` 上分别注入参考星及两个目标星的正负单差故障，
对应 ENU 位置偏移均落在计算的 HPL/VPL 内。
该有限合成测试不代表真实数据全部故障模式的覆盖率验证。

## RANSAC + FILTER_SS 单星故障实验（2026-10-09）

统一由 `ENABLE_RTK_RANSAC_ARAIM=1` 控制；当前 `ENABLE_RTK_RANSAC_SLOPE=0`。
原先独立的 `ENABLE_RTK_RANSAC_SHARED_MASK` 已合并到该总开关。
主解完成当前历元 RANSAC 后，将载波保留掩码复制给删星子滤波器；
子滤波器不重置掩码、不重新运行 RANSAC，再剔除假设故障卫星的全部观测。
因此每个子解建立在同一个主解一致集上。

可复现实验入口：`tmp/ransac_araim_20261009/run.ps1`（五组数据首300秒，
浮点解、FDE关闭），`run_probe.ps1`（实际RANSAC内核已知故障注入），
`analyze.py`（保护级与参考位置比较）。结果见同目录 `report.md`、
`summary.json`、`epochs.json`、`fault_probe.jsonl`。
这是现有 ARAIM式滤波解分离实现的串联试验，不是完整标准ARAIM风险模型。
筛选只覆盖载波；默认12sigma不能保证剩余最多一个故障。
仍需处理选择效应、历史故障、相关故障及子解缺失时的完整性不可用标志。
