# 简单中介路径数据设计

入口为 `notebooks/02_teaching_mediation.ipynb`；完整模型函数为 `R/simulation.R` 中的 `simulate_teaching_dataset()`，参数预设为 `default_config("teaching")`。

## 固定机制、随机观测

n=300名受试者随机分为150名对照和150名暴露者。m=1000条候选路径共享X，每列分别有自己的M和Y。固定四类路径数量：稀疏880/50/50/20，密集400/200/200/200，按H00/H10/H01/H11排列；随后随机打乱路径列顺序。若m或比例改变，用最大余数法取整，保持总数为m。

非零α=0.2τ，非零β=0.5τ，直接效应γ=0.2；没有随机效应系数，只有受试者噪声与抽样随机性。H10只有α非零，H01只有β非零，H11两者均非零。两个系数scale须非零，以保持真值定义与数据一致。

- 线性：M=αX+e，Y=0.2X+βM+ε；e、ε独立N(0,1)。
- 已测混杂：Z~N(0,1)，M和Y各增加0.3Z，两个回归都正确调整Z。
- 二元：M不变，Y~Bernoulli(plogis(0.2X+βM))；β是条件log odds系数。

固定效应的`alpha_noise_variance`、`beta_noise_variance`、`direct_sd`为0。`confounder_max`在此profile中解释为固定θ/δ的大小，在旧profile中仍为Uniform上界。`exposure_probability`控制固定暴露人数，人数为round(n×p)，再随机分配受试者。

## 与01 notebook的区别

改变DGP：固定分组代替Bernoulli抽样，固定状态数量代替随机状态，固定同向效应代替随机高斯效应，固定直接效应与混杂系数代替随机效应。n/m/τ/重复数/q/稀疏密集比例/场景和方法选项一致；同一套回归、适配器、指标、区间、配对汇总和主图直接复用。

正式seed=20261010，独立于开发检查20261008。预先选α=0.2τ，是因为平衡分组时SE(α)约0.115，三个τ对应约2.3/2.6/3.3个SE，使功效可以逐步变化。这里不是原论文参数复现。

方法不接收状态或真实系数。MLFDR仍采用允许零先验方差的`paper_em`，不按真值初始化；该DGP与其可表达的固定效应模型比较匹配，应说明这一点。跨notebook的数字不能单独证明方法排名普遍成立。

## 结果和复现

各条件40次，共720份数据、2160次方法调用。FDR=mean(V/max(R,1))，Power=mean(S/A)，继续输出MCSE、近似95%区间和配对差异。先生成RDS、索引CSV和示例长表，再读取分析。生成文件只保存在本地，Git忽略data与tmp。

`docs/current_teaching.rds`与`.txt`保存新入口；不会改写`current_comparison`。实际记录见 `TEACHING_RUN_REPORT.md`。生成模块的代码指纹变化后，旧结果仍保留其运行时指纹，未来重跑旧参数会进入新的生成指纹目录，这是缓存隔离行为。

`scripts/execute_notebook.py`支持传入notebook路径，默认仍运行01；`scripts/validate_teaching.R`逐份验证新实验。主notebook五张图包含两张路径示例和三张比较图；新增方法接口沿用旧版。
