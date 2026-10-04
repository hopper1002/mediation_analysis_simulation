# ASHMED_Adaptive 实际开发与验证结果

日期：2026-10-05。原版插件与第三 notebook 保留；新增第四 notebook。以下数字全部来自实际运行的保存结果，没有重新生成数据。

## 选中了什么

选定 **null10**：学习α/β非零中心、两个尺度的矩形字典、非零中心点质量、去掉零中心替代字典，加三个复合零状态的总伪计数λ=10。典型64成分；任一学习中心与零距离不足0.5倍中位标准误时退回原版。q始终0.05。推导见 [ASHMED_ADAPTIVE_METHOD.md](ASHMED_ADAPTIVE_METHOD.md)。

1–3次比较结构，4–8次选择λ，9–12次确认，13–40次固定参数最终比较。最后每条件28次、两套共1008份数据、5040个五方法结果；只新拟合改进版，原四方法按相同case_id提取。该划分是历史保存数据中的内部保留重复，原对照此前已在这些数据上运行，不能称为独立外部验证。

## 开发尝试与失败记录

仅放开零中心尺度导致两套整体FDR约0.38/0.35；增加学习中心但保留零中心替代成分也出现明显膨胀。近零替代字典与精确零状态难以识别，灵活性增加并不自动使FDR改善。结构开发记录保存在 `structure_scores.csv`，所有候选和源代码快照保留。去掉零中心替代字典后，再比较正则强度：

| 候选 | 数据 | 平均FDR | 最差条件FDR | 平均Power | 相比原版Power差 | 开发门槛 |
| --- | --- | --- | --- | --- | --- | --- |
| null0 | paper_adjusted | 0.0671 | 0.1773 | 0.7355 | 0.3038 | FALSE |
| null0 | teaching_fixed | 0.0523 | 0.1146 | 0.7298 | 0.3217 | TRUE |
| null10 | paper_adjusted | 0.0545 | 0.1065 | 0.7165 | 0.2847 | TRUE |
| null10 | teaching_fixed | 0.0477 | 0.1059 | 0.7154 | 0.3073 | TRUE |
| null30 | paper_adjusted | 0.0451 | 0.0984 | 0.6929 | 0.2612 | TRUE |
| null30 | teaching_fixed | 0.0452 | 0.0879 | 0.6963 | 0.2882 | TRUE |
| null90 | paper_adjusted | 0.0403 | 0.1014 | 0.6484 | 0.2167 | TRUE |
| null90 | teaching_fixed | 0.0343 | 0.0582 | 0.6563 | 0.2482 | TRUE |

门槛为每套整体FDR≤0.06、最大条件FDR≤0.15、Power高于原版且无错误；合格者按平均Power选。门槛用于小样本开发，不是FDR控制保证。两套确认结果如下，确认通过后不再按13–40次结果改参数：

| 数据 | 确认FDR | 最差条件FDR | 确认Power | 原版Power | 确认门槛 |
| --- | --- | --- | --- | --- | --- |
| paper_adjusted | 0.0441 | 0.0817 | 0.7174 | 0.4495 | TRUE |
| teaching_fixed | 0.0438 | 0.0721 | 0.7037 | 0.3974 | TRUE |

## 固定方法后的结果

每次FDP=V/max(R,1)，Power=S/A；FDR为每条件28次FDP平均。以下整体表为18个条件等权平均，MCSE按重复编号块计算；必须结合逐条件FDR判断，整体均值不代表每个条件都控制在0.05。

### 第一套：论文模型增强信号

| 方法 | 整体FDR | FDR MCSE | 整体Power | Power MCSE |
| --- | --- | --- | --- | --- |
| HDMT | 0.0359 | 0.0017 | 0.5255 | 0.0078 |
| MDACT | 0.0394 | 0.0023 | 0.5628 | 0.0077 |
| MLFDR | 0.0469 | 0.0022 | 0.7402 | 0.0067 |
| ASHMED | 0.0230 | 0.0018 | 0.4610 | 0.0076 |
| ASHMED_Adaptive | 0.0455 | 0.0021 | 0.7347 | 0.0066 |

第一套：论文模型增强信号：原版平均Power 0.461 → 改进版 0.735，提高 **27.4个百分点**；改进版整体FDR 0.0455，逐条件均值范围 0.0274–0.0827。

最高FDR出现在 confounded / sparse / τ=1.3：均值0.0827，MCSE=0.0261，近似95%区间[0.0291, 0.1363]。该条件需保留，不能据整体平均删除或忽略。

- 相比HDMT整体Power差+20.9个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。
- 相比MDACT整体Power差+17.2个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。
- 相比MLFDR整体Power差-0.6个百分点；18条件中1个均值更高，0个配对近似95%区间下界>0。
- 相比ASHMED整体Power差+27.4个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。

这些区间未做多重比较校正，只用于描述蒙特卡洛差异；不能把18个条件的区间当作研究性显著性筛选。

预先指定示例条件：线性、密集、τ=1.5（28次均值）：

| 方法 | FDR | FDR下界 | FDR上界 | Power |
| --- | --- | --- | --- | --- |
| ASHMED | 0.0150 | 0.0106 | 0.0194 | 0.4043 |
| ASHMED_Adaptive | 0.0559 | 0.0490 | 0.0627 | 0.8112 |
| HDMT | 0.0358 | 0.0297 | 0.0419 | 0.5581 |
| MDACT | 0.0450 | 0.0394 | 0.0505 | 0.6011 |
| MLFDR | 0.0563 | 0.0487 | 0.0638 | 0.8142 |

改进版没有条件的近似FDR区间下界超过0.05；部分点估计高于q。这不构成普遍控制证明，区间包含q也不等于已经证明控制。

504次新拟合全部权重及中心迭代收敛；最大KKT gap=9.98783e-05；近零中心回退1次。完整逐条件表、配对表、诊断、示例后验及图形见 `results/ashmed_optimization_20261005/validation_paper_adjusted/`。

### 第二套：固定效应教学数据

| 方法 | 整体FDR | FDR MCSE | 整体Power | Power MCSE |
| --- | --- | --- | --- | --- |
| HDMT | 0.0353 | 0.0026 | 0.4818 | 0.0049 |
| MDACT | 0.0438 | 0.0035 | 0.5187 | 0.0057 |
| MLFDR | 0.0486 | 0.0019 | 0.7253 | 0.0050 |
| ASHMED | 0.0204 | 0.0020 | 0.3970 | 0.0042 |
| ASHMED_Adaptive | 0.0466 | 0.0022 | 0.7148 | 0.0052 |

第二套：固定效应教学数据：原版平均Power 0.397 → 改进版 0.715，提高 **31.8个百分点**；改进版整体FDR 0.0466，逐条件均值范围 0.0348–0.0591。

最高FDR出现在 confounded / sparse / τ=1.5：均值0.0591，MCSE=0.0144，近似95%区间[0.0297, 0.0886]。该条件需保留，不能据整体平均删除或忽略。

- 相比HDMT整体Power差+23.3个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。
- 相比MDACT整体Power差+19.6个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。
- 相比MLFDR整体Power差-1.0个百分点；18条件中0个均值更高，0个配对近似95%区间下界>0。
- 相比ASHMED整体Power差+31.8个百分点；18条件中18个均值更高，18个配对近似95%区间下界>0。

这些区间未做多重比较校正，只用于描述蒙特卡洛差异；不能把18个条件的区间当作研究性显著性筛选。

预先指定示例条件：线性、密集、τ=1.5（28次均值）：

| 方法 | FDR | FDR下界 | FDR上界 | Power |
| --- | --- | --- | --- | --- |
| ASHMED | 0.0105 | 0.0050 | 0.0159 | 0.3427 |
| ASHMED_Adaptive | 0.0520 | 0.0461 | 0.0578 | 0.7984 |
| HDMT | 0.0362 | 0.0286 | 0.0439 | 0.5289 |
| MDACT | 0.0416 | 0.0339 | 0.0493 | 0.5780 |
| MLFDR | 0.0517 | 0.0458 | 0.0577 | 0.8020 |

改进版没有条件的近似FDR区间下界超过0.05；部分点估计高于q。这不构成普遍控制证明，区间包含q也不等于已经证明控制。

504次新拟合全部权重及中心迭代收敛；最大KKT gap=9.99334e-05；近零中心回退0次。完整逐条件表、配对表、诊断、示例后验及图形见 `results/ashmed_optimization_20261005/validation_teaching_fixed/`。

## 怎样解释这次改进

原版将替代效应放在零中心，主要通过增大先验方差表达强信号；当前数据的效应有固定方向且较集中，这会浪费先验质量。学习非零中心后，小方差及非零点质量能表示更集中的效应，独立尺度也更适合两段信号不同的情况，因而识别效率上升。去掉近零替代字典与零状态正则化则抑制误把零路径归入H11。增加Power的同时是否控制FDR，仍必须由保存真值评估。

改进版不是全面替代MLFDR；MLFDR本来就含可学习的有方向成分，当前网格中可有更高Power。HDMT/MDACT使用p值校准，而本方法借助跨路径效应分布，优势依赖此分布建模质量。当前方法只有一个非零中心，对混合正负方向、多峰、强依赖、极弱或更稀疏信号仍需另行测试。二元结局采用GLM正态摘要近似；经验贝叶斯后验没有计入学习中心、字典和权重的估计不确定性。

## 核查与复现

```text
ASHMED_Adaptive final evaluation validation passed.
paper_adjusted: 504 held-out cases / 2520 five-method results; all adaptive fits converged; max gap 9.98783e-05; fallbacks 1; raw MD5 720/720.
teaching_fixed: 504 held-out cases / 2520 five-method results; all adaptive fits converged; max gap 9.99334e-05; fallbacks 0; raw MD5 720/720.
Total: 1008 saved cases and 5040 method results recomputed from truth and decision caches.
All summaries, MCSE/interval and paired tables recomputed and matched CSV.
Frozen method fingerprint and options match development selection; repeats 13:40 disjoint from development/confirmation.
Protected original files: 18/18 unchanged.
No original simulation data regeneration; baseline decisions and metrics exactly preserved.
```

数学单元检查 `tests/check_adaptive_ashmed.R` 已验证一般非零中心密度及共轭产品矩、已知惩罚最优解、目标单调性、原版等价开关与回退。notebook/HTML检查另外验证全部R单元已执行、无错误、10张图一致且无生成数据调用。

最终重跑：`Rscript --vanilla scripts/run_adaptive_ashmed.R`；完整验证：`Rscript --vanilla scripts/validate_adaptive_ashmed.R`、`python scripts/validate_adaptive_ashmed_notebook.py`。源码/选项变化要求重新锁定，不能反复使用本轮验证结果选最优仍声称验证独立。
