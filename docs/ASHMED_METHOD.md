# ASHMED：附件数学方法的独立 R 实现

来源为用户提供的 `D:/A_Mediation_analysis/论文和笔记/ASHMED.pdf`，题名 *ASHMED: Adaptive Shrinkage for High-Dimensional Mediation Analysis*，作者 Huaying Fang。实现依据第2.1–2.2节，尤其PDF第9–15页的似然、先验、EM和local FDR公式；未用该文实验部分设计数据、设置信号或选择排名。附件的软件地址仍为占位地址，因此本项目实现不称为作者官方软件，也不依赖ASHMED包。

代码在 `R/methods/plugins/ashmed.R`；比较调度在 `R/ashmed_comparison.R`。基础R即可运行新增方法。所有方法继续共用 `R/estimation.R` 的系数和标准误，接口不包含模拟真值。

## 1. 输入和采样模型

对每条候选路径i，输入 `alpha_hat`, `beta_hat`, `var_alpha`, `var_beta`，以及可选 `pathway_id`, `cov_alpha_beta`。这里测试数记为m，对应文中p；本实验每份m=1000。

\[
\widehat\theta_i=(\widehat\alpha_i,\widehat\beta_i)^T,
\qquad \widehat\theta_i\mid\theta_i\sim N_2(\theta_i,S_i).
\]

已有两套回归摘要未报告交叉采样协方差，因此按第9页使用工作对角近似

\[
S_i=\operatorname{diag}(v_{\alpha i},v_{\beta i}).
\]

若输入提供 `cov_alpha_beta`，插件使用完整2×2矩阵，并要求正定；不可把真实效应先验相关ρ代替采样误差相关。输入包含非有限估计、非正方差或不正定矩阵时明确失败，避免无效结果进入汇总。

四类状态：H00的α=β=0，H10的α有效/β=0，H01的α=0/β有效，H11两者均有效。连续线性模型下路径间接效应为αβ。

## 2. 固定先验网格

使用第10页的观测标准误中位数构造网格：

\[
r_j=r_{min}(r_{max}/r_{min})^{(j-1)/(J-1)},\quad
u_j=(r_j\widetilde s_\alpha)^2,\quad
v_j=(r_j\widetilde s_\beta)^2.
\]

主比较预先使用J=12、r_min=0.25、r_max=8、R={−0.5,0,0.5}。所有先验成分均以0为中心。

| 状态 | 先验分布 | 成分数 |
|:--|:--|--:|
| H00 | 原点点质量，U=0 | 1 |
| H10 | N2(0,diag(u_j,0))，β恰为0 | 12 |
| H01 | N2(0,diag(0,v_j))，α恰为0 | 12 |
| H11 | N2(0,Σ_jρ) | 36 |

\[
\Sigma_{j\rho}=
\begin{pmatrix}u_j&\rho\sqrt{u_jv_j}\\
\rho\sqrt{u_jv_j}&v_j\end{pmatrix}.
\]

H11严格使用same-scale，即同一个j决定两个先验方差。总数61，不展开为全部α网格×β网格组合，不添加非零先验均值，不用真值选择尺度。网格只通过观测标准误确定，不随EM更新。

## 3. 展开权重与EM的等价性

每个成分k有状态T(k)、固定先验U_k、状态概率π_T及状态内概率ω_k。令w_k=π_Tω_k，得到

\[
L_{ik}=\phi_2(\widehat\theta_i;0,S_i+U_k),\quad
f_i=\sum_kw_kL_{ik},\quad
\ell(w)=\sum_i\log f_i.
\]

约束w_k≥0、Σw_k=1。E步为r_ik=w_kL_ik/f_i，M步为w_new,k=Σ_i r_ik/m。

文中分层M步π_new,T=Σ_iγ_iT/m，ω_new,k=Σ_i r_ik/Σ_iγ_iT，两者相乘即w_new,k；因此展开权重估计的是同一个模型和同一个无惩罚边际似然。拟合后π_T=Σ_{k:T(k)=T}w_k，ω_k=w_k/π_T。

### 需要补全或修正的数值细节

- **初始π不满足归一化**：第13页的(0.85,0.10,0.10,0.05)合计1.10。插件除以1.10，保持相对比例；ω初始化均匀。初始化所有状态要求正权重，以便EM探索全部成分。
- **零状态质量**：若拟合π_T=0，状态内ω任意且不能由除法识别。输出用均匀ω占位；w及该状态后验仍为0。
- **密度下溢**：先计算2×2正态log密度，再逐行减最大log密度后指数化；每行减去的常数在最终loglik中加回。
- **退化先验**：H00及两轴成分原样保留。似然的S_i+U_k正定，不用给点质量添加小方差。
- **收敛精度**：文中以相对变化1e−6、最大200轮为示例。这里使用相对变化1e−8、最多5000轮，并额外检查KKT gap≤1e−4。它们是求解精度设置，不改变统计模型；初次检查发现少量高度重叠成分下严格权重收敛很慢，因此保留似然收敛和明确的优化误差界。
- **数值加速**：EM后作SQUAREM形式的外推提议，投影到单纯形后再EM。只接受不低于普通EM两步似然的提议，否则退回普通EM。每20轮或EM停滞时允许向最大梯度成分作顶点线搜索，可恢复被投影到0的有用成分。全部步骤只优化同一个无惩罚似然，每轮记录loglik、gap、相对变化和更新类型。`accelerate=FALSE`关闭外推；顶点收敛补救仍保留。

固定网格下，ℓ是非负权重单纯形上的凹函数；文中“local maximum”表述可进一步澄清为全局目标，成分重复或高度相似时最优权重可能不唯一。定义

\[
g_k=\frac1m\sum_i\frac{L_{ik}}{f_i},\qquad
\delta=\max(0,\max_kg_k-1).
\]

KKT条件为正权重成分g_k=1，零权重成分g_k≤1；又有Σw_kg_k=1。由凹函数切线界，任何全局最优解w*满足ℓ(w*)−ℓ(w)≤mδ。`dual_gap`记录δ；满足阈值后每条路径平均对数似然的最优值差≤1e−4。该界约束目标值，不保证不唯一的每个ω被精确识别。优化未收敛时明确产生警告且保存状态，不伪装为成功收敛。

主实验全量拟合1000条路径。可设`fit_subset`对观测路径随机取子集拟合，再在全部路径上分块预测；网格仍取全部观测标准误的中位数。随机性使用实验引擎固定的方法种子。主比较`fit_subset=NULL`，没有子集抽样。

## 4. 后验状态与筛选

\[
\gamma_{i,T}=\sum_{k:T(k)=T}r_{ik},\qquad
\operatorname{lfdr}_i=\gamma_{i,00}+\gamma_{i,10}+\gamma_{i,01}=1-\gamma_{i,11}.
\]

每条路径四状态后验之和为1。使用原项目 `lfdr_step_up()`：排序后找最大的k使前k项平均lfdr≤q。分数并列时只考虑并列组末尾；零发现和全发现均可。筛选向量保持原路径顺序。

选择集合平均local FDR描述插件经验贝叶斯模型中的预期FDP。真值定义的频率学经验FDR是重复间FDP平均，二者不能混称；先验/似然错设及权重估计误差可能使实际FDR偏离目标。

## 5. 正态共轭补全效应收缩与产品方差

附件引言提到收缩效应及不确定性，第2.2节输出主要明确了状态后验与local FDR。这里按同一模型推导并实现后验效应矩：

\[
b_{ik}=U_k(S_i+U_k)^{-1}\widehat\theta_i,\quad
C_{ik}=U_k-U_k(S_i+U_k)^{-1}U_k.
\]

无需U_k可逆，因此也覆盖两轴及点质量。混合后的E(θ|data)=Σr_ik b_ik，Cov(θ|data)=Σr_ik(C_ik+b_ikb_ikᵀ)−μμᵀ。

对于每个成分，记后验均值a,b，方差A,B，协方差C：

\[
E(\alpha\beta\mid k)=ab+C,
\]

\[
E(\alpha^2\beta^2\mid k)=a^2b^2+a^2B+b^2A+4abC+AB+2C^2.
\]

先按r_ik混合这两个矩，再用Var(αβ)=E(α²β²)−E(αβ)²。轴状态与H00的αβ恰为0。输出`alpha_mean/sd`、`beta_mean/sd`、`alpha_beta_cov`、`mediation_mean/sd`。其中`mediation_mean`不能用两个边际均值相乘代替。

这些是把拟合权重当作已知的经验贝叶斯后验矩，不包含权重估计误差。混合后验含原点/轴质量，因此不把正态近似“均值±1.96SD”标成精确95%可信区间。文中讨论的其他方向错误率等未给出完整操作定义，本比较不冒充实现额外的lfsr阈值。

## 6. 原有数据上的适用范围

连续线性和已测混杂数据使用原OLS摘要；混杂Z在两段回归均调整。二分类数据沿用原Logistic摘要，按大样本正态似然作**GLM摘要扩展**。PDF方法和后文局限将线性模型作为主要推导，因此二分类比较明确标注扩展；αβ在这里是条件系数乘积，不称为结局概率尺度的自然间接效应。

两套已有DGP均含有明确非零中心的有效效应，教学数据还有固定非零点质量。ASHMED的零中心、same-scale先验可能错设，使筛选较保守；MLFDR的非零均值成分可能较适合。此解释是模型与数据设计的比较，不引用ASHMED论文实验的排名。保持数据、q和原方法配置，不按结果加入非零均值、放大信号、改变网格或删重复。

## 7. 接口和输出

`register_plugin(registry)`自动注册ASHMED；`run_ashmed(input,q,options)`返回统一的`reject`、`score=lfdr`、`score_type="local_FDR"`，外加：

- `posterior`：逐路径四状态概率、lfdr、决定、效应后验矩；
- `components`：61个固定先验及拟合联合权重w、状态内权重ω；
- `diagnostics`：π、选中集合平均lfdr、优化轨迹、收敛、KKT gap、设置、全量/子集大小及线性/GLM标签。

插件源码、参数、R/包版本、数据MD5参与缓存/运行记录。调整ASHMED选项只影响其缓存；原方法可继续复用。入口参数集中在第三notebook开头，加入其他方法时追加插件和具名选项即可。

数值外推的背景参见[Varadhan与Roland（2008）的SQUAREM原始论文](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9469.2007.00585.x)。零中心混合先验的适应性及限制可参考[Stephens（2017）的原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5379932/)；这里仍严格拟合ASHMED附件规定的无惩罚61成分模型，没有替换为ashr包或该文的惩罚估计。

## 8. 数学与交付核查

`tests/check_ashmed.R`用独立矩阵运算验证正态密度及共轭矩，验证已知内点/边界最优解、目标单调性、数值加速与非外推求解同目标、保存数据上的概率归一化、行顺序及分块不变性、单行分块、不合法协方差和接口真值隔离。该测试不调用模拟数据生成器。

`scripts/validate_ashmed.R`核对两套全部保存数据、每次拒绝和评价指标、原方法基线不变、ASHMED后验/权重/收敛/平均lfdr阈值、重算汇总及配对结果。`scripts/validate_ashmed_notebook.py`检查R内核、所有代码单元执行、无错误输出及八幅嵌入图。实际运行统计见 `ASHMED_RUN_REPORT.md`。
