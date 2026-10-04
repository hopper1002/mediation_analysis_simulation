# 模拟设计与原论文的对应关系

本文件保留字面正文与历史源码配置的说明。**当前主notebook使用paper_comparison/paper_adjusted：非零α均值改为0.35τ，其他主要分布沿用正文。** 当前设计和方法实现以[PAPER_MODEL_COMPARISON.md](PAPER_MODEL_COMPARISON.md)为准；下文“默认”指此前弱信号演示。

## 基础设计

来源为用户附件 MLFDR.pdf 的发表版，Roy & Zhang (2026)，第 2.2 节（第 5–8 页），不是只按照 README 中的简化示例。对每个条件生成 n 个受试者、m=1000 个中介-结局路径。

X 在所有路径间共享，服从 Bernoulli(0.1)；每条路径有自己的 M_i 和 Y_i。四个潜在状态按 H00、H10、H01、H11 顺序抽样，稀疏概率 (0.88,0.05,0.05,0.02)，密集概率 (0.4,0.2,0.2,0.2)。H11 是唯一真实备择，不用拟合结果确定真值。

默认 `paper2026`：非零 α=0.05τ+N(0,1/n)，非零 β=−0.5τ+N(0,4/n)；相应零状态设为精确的 0。直接效应 γ=N(1,0.5)，将第二参数解释为方差，即 R 的 `sd=sqrt(0.5)`。所有误差独立 N(0,1)。τ 均值不再除以 √n，这是第 2.2 节明确写出的形式。第 5.1 节另有 √n 参数化；本模拟明确采用 simulation 章节的定义，避免混用。

三种模型：

1. 线性：M_i=α_i X+e_i；Y_i=β_i M_i+γ_i X+ε_i。
2. 已测混杂：Z=N(0,1)，M 加 θ_i Z，Y 加 δ_i Z，θ_i、δ_i 独立 U(0,0.5)。两个回归正确加入 Z。
3. 二元结局：M 同线性，Y_i~Bernoulli(plogis(β_i M_i+γ_i X))；α 用 OLS，β 用 logistic GLM。

回归包含截距，与作者 `lm(M ~ X)` 的做法一致。二元结局 logistic 模型无额外 DGP 截距，拟合时有截距。没有将所有中介同时加入一个公共结局模型；这种扩展需要另写 DGP 和估计模块。

## 源码差异与明确选择

| 项目 | 2026 论文正文 | 作者 linear_model.R / binY.R | 工程选择 |
|:--|:--|:--|:--|
| α 非零噪声 | 方差 1/n | `rnorm(...,sd=kap)`，kap=1 | paper2026 用 sqrt(1/n)；source_code 用 sd=1 |
| β 非零噪声 | 方差 4/n | `rnorm(...,sd=psi)`，psi=4 | paper2026 用 sqrt(4/n)；source_code 用 sd=4（方差16） |
| 直接效应 | N(1,0.5) | `rnorm(1,0.5)`，均值0.5、sd1 | paper2026: mean1、sd sqrt(.5)；source_demo: mean.5、sd1 |
| n | 100 与 300 | 各脚本部分条件 n 不完全相同 | 两个 n 全交叉 |
| 重复数 | 正文第 2.2 节未明确说明 | n.sim=250 | 演示5，paper_grid 250 |
| MLFDR 输入/返回 | 排序 local FDR step-up | 旧脚本返回向量；包 localFDR 返回 list | 使用已安装包的 `$lfdr` |
| MLFDR 拒绝边界 | 最大合格累计均值 | package while(k<m) 存在零/全拒绝边界 | 使用独立修正版；并列分数按整组处理 |
| 原始矩阵方向 | 统计模型约定 | m×n | n×m，明确数据字典 |

工程不将这些配置称为论文逐位复现。改变 DGP 分布、样本量、方法实现或随机种子都会改变结果。

## 方法实现

HDMT：null_estimation 和 fdr_est 使用命名空间调用，不依赖 `library()` 顺序。默认 `exact=0`，与作者线性/二元脚本同样使用默认渐近估计。若配置 exact=1，记录成另一组推断配置。空合格集明确返回零拒绝。

MDACT：使用作者配套 `balancing_DACT_control_DR_adjust` 和 `DR_DACT_thr_adjust`。保留其 1e-3 混合权重下界、F00 数值积分、F01/F10 经验计算、基于统计量的二分阈值搜索和 `<`/`<=` 约定。二分搜索隐含单调性假设，工程没有声称验证其有限样本单调性或替换成全阈值扫描。只将积分失败后的无限 rounding 重试改成 k=10,…,0 的有限重试并明确报错。原副本未执行。

MLFDR：localFDR 的一步 EM 参数默认 eps=.01，其余采用已安装包默认值；输入 var_alpha/var_beta 是 SE²。`lfdr_step_up` 按累计均值 ≤ q 选择最大合格并列组的结尾，零拒绝和全拒绝均可返回长度 m 的合法逻辑向量。它修复了包阈值函数的边界，不改变包内 EM。

浮点误差可令 package local FDR 等于 1+2.22e-16。本工程容忍1e-12内的越界并截到[0,1]，记录 roundoff_clipped；超过该幅度仍报错，不能把数值失败伪装成有效概率。

算法接口只提供观测、系数统计量及 n/m/场景/case_id。case_id 用于识别文件，不能用它推导真值或参数。评估在算法运行后单独计算。

## 指标和局限

重复 b：R_b 为拒绝数，V_b 为错误拒绝，S_b 为正确拒绝，A_b 为 H11 数。FDP_b=V_b/max(R_b,1)，power_b=S_b/A_b。A_b=0 时 power=NA。经验 FDR 是 mean(FDP_b)，并非 sum(V_b)/sum(R_b)（后者是不同指标）。

每个条件独立按成功记录汇总；n_total、n_success、n_error、n_warning、每指标有效 n 同时报告。失败值是 NA，不算作零。logistic 警告和 EM 警告不全局屏蔽，记录到 diagnostics；有警告的有效输出仍计入成功统计，因此应查看原文警告内容。

SD 是重复之间的标准差，MCSE=SD/sqrt(B_valid)。图中区间是 t(B_valid−1) 的近似均值区间，截到 [0,1]；B=1 或无有效重复时区间为 NA。边界数据或极少重复会让该近似很不稳定。默认5次只是验证和展示流程；结果差异不支持稳定优劣结论。

不包括论文 SVA、交互、潜在混杂、多正态复合备择、TCGA 和 global-null 实验。未来可通过新 DGP/估计模块扩展；新比较算法只需插件接口。
