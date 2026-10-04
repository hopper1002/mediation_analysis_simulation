# 运行与交付检查记录

执行日期：2026-10-04（Asia/Shanghai）。

主实验：`results/run_45a379a50352`。
数据入口：`data/data_6fdb451c6d0d/manifest.rds`。
源码设定对照：`results/run_5580e7f5e7ad`。

## 实际执行

- 主实验180份数据、540次方法调用，全部成功；3场景×2混合比例×2样本量×3个τ×5重复。
- 额外源码设定对照6份数据、18次方法调用，全部成功，每条件1次。
- 主实验有警告的调用：142；均保存在 diagnostics.csv。
- MLFDR 浮点截界次数（逐路径累计）：63；只接受1e-12内的舍入误差。
- 已核查全部主数据MD5、每方法拒绝向量与真值计算的FDP/power、CSV汇总与6张PNG。
- notebook逐单元执行；HTML导出和图形布局另由交付检查核查。

## 主实验结果的解释

各条件经验FDR均值范围：0.0000–0.6000；power均值范围：0.0000–0.0222。
发表版正文的α效应很弱，本次演示的power低。源码未缩放的系数噪声与正文不同，额外对照的检出明显不同。
5次重复的波动较大；某个条件的均值超出目标q不等于稳定违背或验证渐近保证。单次源码对照更不能作FDR结论。

| 条件 | 方法 | 单次FDP | 单次power | 发现数 |
|:--|:--|--:|--:|--:|
| linear/sparse, n=300, τ=1.9 | HDMT | 0.1000 | 0.6000 | 10 |
| linear/sparse, n=300, τ=1.9 | MDACT | 0.1111 | 0.5333 | 9 |
| linear/sparse, n=300, τ=1.9 | MLFDR | 0.1000 | 0.6000 | 10 |
| confounded/sparse, n=300, τ=1.9 | HDMT | 0.0000 | 0.5625 | 9 |
| confounded/sparse, n=300, τ=1.9 | MDACT | 0.0909 | 0.6250 | 11 |
| confounded/sparse, n=300, τ=1.9 | MLFDR | 0.0909 | 0.6250 | 11 |
| binary/sparse, n=300, τ=1.9 | HDMT | 0.0000 | 0.6957 | 16 |
| binary/sparse, n=300, τ=1.9 | MDACT | 0.0000 | 0.6957 | 16 |
| binary/sparse, n=300, τ=1.9 | MLFDR | 0.0000 | 0.7391 | 17 |
| linear/dense, n=300, τ=1.9 | HDMT | 0.0325 | 0.5920 | 123 |
| linear/dense, n=300, τ=1.9 | MDACT | 0.0400 | 0.5970 | 125 |
| linear/dense, n=300, τ=1.9 | MLFDR | 0.0496 | 0.6667 | 141 |
| confounded/dense, n=300, τ=1.9 | HDMT | 0.0234 | 0.6410 | 128 |
| confounded/dense, n=300, τ=1.9 | MDACT | 0.0234 | 0.6410 | 128 |
| confounded/dense, n=300, τ=1.9 | MLFDR | 0.0625 | 0.6923 | 144 |
| binary/dense, n=300, τ=1.9 | HDMT | 0.0500 | 0.6584 | 140 |
| binary/dense, n=300, τ=1.9 | MDACT | 0.0559 | 0.6683 | 143 |
| binary/dense, n=300, τ=1.9 | MLFDR | 0.0347 | 0.6881 | 144 |

## 正确性与工程检查

- run_tests.R：确定性种子、扩展网格不扰动原场景种子、OLS与lm等价、状态/方差生成、step-up零/全/并列/浮点边界、所有适配器。
- check_pipeline.R：文件往返、串行/Windows PSOCK、缓存命中、增加插件复用、q变化缓存隔离、checksum拒绝、方法失败为NA。
- 使用已安装HDMT和MLFDR；MDACT使用有来源和许可证记录的配套代码。
- 原MLFDR-main未修改；项目运行无需该原仓库路径。

环境版本见本次结果目录的 package_versions.csv 和 sessionInfo.txt。此记录是本机运行验证，不承诺其他包版本位级复现。
