# MLFDR Simulation Lab

使用 R / Jupyter 建立可扩展的高维中介效应模拟实验。当前主入口是**论文模型的信号增强比较**：保留正文三种模型与主要分布，只增强α均值，展示HDMT / MDACT / MLFDR的真实FDR、功效与配对差异。参数控制、先生成文件再读取、拟合诊断和新方法接口集中在同一个notebook。此版不能称为原图的严格复现。

## 从这里开始

打开 **`notebooks/01_simulation_workbench.ipynb`**，选择 **R** 内核，按顺序运行。交付的 notebook 已执行，保留结果；同名 `.html` 可直接在浏览器中阅读，不需要 R 或 Jupyter。

新增 **`notebooks/02_teaching_mediation.ipynb`**：用各半随机分组、固定四类路径数量和固定效应生成更容易理解的数据，展示X→M与调整X后的M→Y。HDMT/MDACT/MLFDR、FDR/Power、样本量、信号网格、40次重复及主图与01一致。设计见[TEACHING_SIMULATION.md](docs/TEACHING_SIMULATION.md)，实际运行见[TEACHING_RUN_REPORT.md](docs/TEACHING_RUN_REPORT.md)。命令行重跑：`Rscript --vanilla scripts/run_experiment.R teaching`。

新增 **[第三个notebook：ASHMED比较](notebooks/03_ashmed_comparison.ipynb)**（同名HTML可直接阅读）：依据用户附件ASHMED的数学方法独立实现四状态零中心自适应收缩，在前两套已经保存的正式数据上比较ASHMED/HDMT/MDACT/MLFDR。**第三个notebook只读取数据，不生成数据**；沿用q、FDR/Power、40次重复、MCSE、配对比较和主图，并验证原三方法结果与数据文件不变。数学及补全见[ASHMED_METHOD.md](docs/ASHMED_METHOD.md)，实际结果见[ASHMED_RUN_REPORT.md](docs/ASHMED_RUN_REPORT.md)。二元结局明确标注GLM正态摘要扩展。新插件只依赖基础R。

第三个notebook的开头集中两个已有`data_id`及对应`baseline_run`、新增方法设置和worker数。缺少原始数据时先在01/02生成并保存数据，再把03入口改为实际保存的manifest与对应结果目录；03不会自动生成替代数据。命令行比较：`Rscript --vanilla scripts/run_ashmed.R`；检查：`Rscript --vanilla tests/check_ashmed.R`、`Rscript --vanilla scripts/validate_ashmed.R`、`python scripts/validate_ashmed_notebook.py`。

GitHub 仓库：[hopper1002/mediation_analysis_simulation](https://github.com/hopper1002/mediation_analysis_simulation)。仓库包含完整生成代码、已执行 notebook、文档和全部历史／当前实验结果；**不上传原始模拟数据、估计缓存或临时文件**。克隆后按顺序运行 notebook，即可先生成并保存数据，再从文件读取进行比较；无需下载额外数据附件。

在新机器的仓库根目录，也可通过 `Rscript --vanilla scripts/run_experiment.R paper_comparison` 重建当前主实验。生成参数和逐次种子由配置确定；位级复现仍需要与报告相同的 R／包版本。历史报告中的 `D:/...` 路径记录原执行位置，重跑时将使用当前工程路径。

在工程根目录启动（PowerShell）：

```powershell
Set-Location 'D:\A_Mediation_analysis\code\mlfdr_simulation_lab'
jupyter notebook notebooks/01_simulation_workbench.ipynb
```

本机 R 4.6.1、IRkernel 1.3.2、HDMT 1.0.5、MLFDR 0.1.0、ggplot2 4.0.3 已发现。项目使用已安装包；MDACT 没有独立安装包，使用作者源码中配套实现。依赖不会自动安装。移到另一台机器时先安装 `HDMT`, `MLFDR`, `ggplot2`, `IRkernel`，再运行 `IRkernel::installspec()` 注册 R 内核。MLFDR 的依赖 NMOF 会由常规包安装处理。

本机主 `python.exe` 缺少 nbformat，但 `D:\miniconda3\envs\svg310\python.exe` 已具备 notebook 工具；**R 代码全部在 R 内核运行，Python 只负责执行/导出 notebook**。需要批量执行并重新导出 HTML 时：

```powershell
& 'D:\miniconda3\envs\svg310\python.exe' scripts/execute_notebook.py
```

普通 Jupyter 菜单运行不需要这个 Python 脚本。`scripts/build_comparison_notebook.py` 是当前notebook的作者脚本，**会清除已有输出和手工改动**；日常调参数直接改notebook，不要再运行作者脚本。旧`build_notebook.py`仅用于历史弱信号演示。

新增 **[第四个 notebook：改进 ASHMED](notebooks/04_ashmed_adaptive.ipynb)**（[HTML](notebooks/04_ashmed_adaptive.html)）：测试学习效应中心、独立尺度和复合零状态正则化，保留失败尝试与原版参照。第1–8次用于开发、第9–12次确认，固定参数后仅读取第13–40次，按原FDR/Power指标与HDMT/MDACT/MLFDR比较。它是针对当前两套有方向效应DGP的探索方法，不代表原版论文或普遍FDR保证。推导见[ASHMED_ADAPTIVE_METHOD.md](docs/ASHMED_ADAPTIVE_METHOD.md)，实际结果见[ASHMED_ADAPTIVE_RUN_REPORT.md](docs/ASHMED_ADAPTIVE_RUN_REPORT.md)。

第四个 notebook 从 `results/ashmed_optimization_20261005/locked_method.rds` 读取选定配置；源码或依赖指纹变更时拒绝继续称为固定方法的验证。最终比较入口 `Rscript --vanilla scripts/run_adaptive_ashmed.R`；检查入口 `tests/check_adaptive_ashmed.R`、`scripts/validate_adaptive_ashmed.R`、`scripts/validate_adaptive_ashmed_notebook.py`。开发入口 `explore_adaptive_ashmed.R`、`tune_adaptive_ashmed.R`。在新机器重建需先跑01/02生成数据和原方法结果，再跑03生成原版参照，最后开发/锁定新方法。历史锁含执行环境指纹，软件版本改变后应重新锁定并明确标注复现结果。

## 默认实验

依据用户提供的 Roy & Zhang (2026) 论文第2.2节，覆盖线性、已测混杂、二元结局和稀疏/密集备择。`paper_comparison`预设：m=1000、n=300、τ=1.3/1.5/1.9，每条件40次，共 **720份原始数据 / 2160次方法调用**，原始矩阵未压缩约3.219 GiB。最终验证种子20261007，独立于开发试运行20261005和较弱信号评估20261006。

保持正文的X~Ber(0.1)、Var(h)=1/n、Var(g)=4/n、β均值−0.5τ、γ~N(1,0.5)、连续模型独立N(0,1)误差。**唯一增强的DGP系数是非零α均值：0.05τ→0.35τ**，目的是形成有区分度的中等到强信号。模型、参数来源与局限详见[当前比较设计](docs/PAPER_MODEL_COMPARISON.md)。

MLFDR默认使用论文四成分模型的向量化EM后端`paper_em`，允许0先验方差，使用√n缩放、log密度、3个固定起点并按观测似然选择。原包后端`package`仍可配置，不能将二者混称为原包默认结果。MDACT的F00积分改为数学等价解析CDF以避免数值失败；HDMT使用已安装包接口。方法均不接收模拟真值。

当前实际运行与核查见[COMPARISON_RUN_REPORT.md](docs/COMPARISON_RUN_REPORT.md)，入口在`docs/current_comparison.rds`。此前字面正文弱信号演示完整保存在`notebooks/archive/01_weak_signal_demo.ipynb`及HTML，旧结果与[复核报告](docs/RESULTS_AUDIT.md)保留。notebook另显示τ=1.1的历史压力测试及偏高FDR；开发试运行、编码修复前轮次和数值失败记录均独立保留，不与最终验证混合汇总。

`default_config("paper_grid")` 是**字面正文参数**的大网格，α均值仍为0.05τ，可能再次出现极低功效。其τ=0.1:0.2:1.9、每条件250次、n=100/300，共30,000份数据、未压缩约89.4 GiB，不建议直接作为日常演示。当前预设不含SVA、交互、复合备择和真实数据分析。

首个代码单元是所有参数的控制区；新增方法、q、EM选项只影响推断，不重新生成数据。`demo`和`source_demo`保留为历史配置，当前主比较不使用源码的未缩放系数噪声。

## 工程结构

```text
mlfdr_simulation_lab/
├── notebooks/        # 主 R notebook + 已执行的 HTML
├── R/
│   ├── config.R      # 预设、参数验证、网格与规模估计
│   ├── simulation.R  # 数据生成、RDS/CSV 导出与读取
│   ├── estimation.R  # 共享 OLS / logistic GLM，构造无真值输入
│   ├── experiment.R  # 调度、PSOCK、缓存、逐次指标、诊断
│   ├── reporting.R   # 汇总、MCSE、区间、图形
│   └── methods/      # 统一方法接口；plugins/ 自动发现插件
├── vendor/           # MDACT 作者原始文件、节选实现、来源记录
├── examples/         # 新方法模板（默认不执行）
├── scripts/          # 命令行入口与 notebook 执行/导出
├── tests/            # 统计计算和工程流程的检查
├── docs/             # 设计、数据字典、方法扩展和运行记录
├── data/             # 按生成指纹隔离的数据与估计缓存
└── results/          # 每次实验的 CSV、PNG、配置、环境、方法缓存
```

`data/` 和 `tmp/` 被 `.gitignore` 忽略，生成文件保留在本地；`results/` 中已生成的汇总、图形、逐次结果、配置和诊断随仓库上传。原 `../MLFDR-main` 没有改动。项目运行不依赖该目录存在。

## 命令行与检查

```powershell
Rscript.exe --vanilla scripts/run_experiment.R paper_comparison
Rscript.exe --vanilla tests/run_tests.R
Rscript.exe --vanilla tests/check_pipeline.R
Rscript.exe --vanilla tests/test_paper_em.R
Rscript.exe --vanilla tests/test_mdact_cdf.R
Rscript.exe --vanilla tests/test_hash_locale.R
Rscript.exe --vanilla scripts/validate_comparison.R
& 'D:\miniconda3\envs\svg310\python.exe' scripts/validate_notebook.py
```

集成测试的临时数据和故意失败记录写到 `tmp/pipeline_fixture`，不会混入交付主结果。

命令行和 notebook 调用同一批模块。默认 2 个 PSOCK worker，Windows 可用；调试时设为 1。工作区不执行 `rm(list=ls())`，不修改当前工作目录，不使用作者脚本中的 macOS 绝对路径或 `mclapply(mc.cores>1)`。

每份原始数据带种子、配置、MD5；每次结果保存环境和源码指纹。这是“记录与缓存验证”，不是跨 R/包版本保证位级一致的环境锁定。变化的包版本会让方法缓存失效；正式比较前应固定软件版本。

新增方法参见[ADDING_METHODS.md](docs/ADDING_METHODS.md)，文件内容参见[DATA_SCHEMA.md](docs/DATA_SCHEMA.md)。`RUN_REPORT.md`与`validate_results.R`是历史弱信号演示的记录/检查入口，不代表当前代码版本的验证；当前使用`validate_comparison.R`。

## 来源与许可

- Roy A, Zhang X. 2026. *Powerful large scale inference in high dimensional mediation analysis*. PLOS Computational Biology 22(1):e1013880. [DOI](https://doi.org/10.1371/journal.pcbi.1013880)。用户附件：`D:\A_Mediation_analysis\论文和笔记\MLFDR.pdf`。
- 作者仓库：[asmita112358/MLFDR](https://github.com/asmita112358/MLFDR)。本地参考：`D:\A_Mediation_analysis\code\MLFDR-main`。
- MDACT 节选代码的原始副本、变更说明、SHA256 见 [vendor/PROVENANCE.md](vendor/PROVENANCE.md)。本工程采用 GPL-2，与作者包 DESCRIPTION 声明一致；LICENSE 包含许可全文。
