# 数据与结果格式

## 原始 RDS

`readRDS("<case_id>.rds")` 返回列表：

| 名称 | 格式与含义 |
|:--|:--|
| schema_version | 1 |
| meta | case_id、scenario、mixture、profile、n、m、τ、replicate、seed 等 |
| X | 长度n，所有路径共享的暴露 |
| Z | 已测混杂变量，或 numeric(0) |
| M / Y | n×m，行=受试者，列=pathway_id |
| truth | m行，pathway_id、state、真实alpha/beta/direct/theta/delta、is_nonnull |
| dgp | 数据生成配置，用于读取验证 |

每条路径有独立结局，非所有中介解释公共结局。RDS 保存精确 R 结构；不只保存 z/p 值。真值只用于评估，不直接作为方法输入。

profile区分paper2026（字面正文）、paper_adjusted（论文模型、α均值增强）、source_code（历史源码分布）、teaching_fixed（固定分组/状态数量/效应的简单路径数据）。参数的实际数值以dgp为准，不依赖标签推断。逐次种子保存在meta与manifest；当前种子/缓存哈希用v2序列化，避免Windows命令行和IRkernel的native_encoding头部差异。原始RDS本身仍用v3保存，读取不受影响。

## 索引与跨语言 CSV

manifest.csv 每行一份 RDS：case_id、scenario、mixture、n、m、tau_index、tau、replicate、seed、file、md5。`file` 相对于该 manifest 所在目录，因此搬迁整个数据目录后仍可用 `read_manifest()` 读取；不要只移动 RDS 而丢失 manifest。

manifest.rds 保存 design、dgp 和生成代码指纹。重新读取使用 `read_manifest(path)` + `read_dataset(manifest, case_id)`；后者校验 MD5 和维度。

example_long.csv 列：subject_id、pathway_id、X、Z、M、Y。每行是一条路径的一个样本，X/Z 因共享而重复。无Z时对应列为空。example_long_truth.csv 是单独真值表，不与算法观测合并。

## 推断与结果

每个估计缓存 RDS 含 estimates 和 warnings，估计列在新增方法指南中说明。按数据文件 MD5、估计代码和 R 版本确定缓存文件名。

每方法缓存 RDS：status、error、warnings、elapsed、result；result 包含 reject、可选score、score_type、diagnostics。仅成功调用缓存。耗时是第一次成功调用本方法的耗时，复用时不重复计时；不包含共享回归、读盘或数据生成。

replicate_metrics.csv 每行一个“数据集×方法”，包含实验条件、status、cached、warning/error、原方法耗时、缓存key、R/V/S/A、FDP、power。summary.csv 每行“条件×方法”，提供有效重复数、均值、SD、MCSE、近似95%区间。错误重复明确为 NA 并报告次数。

paired_differences.csv在同一case_id配对方法的FDP/power差异，提供差值、MCSE和区间。fit_diagnostics.csv记录paper_em的MLFDR权重、原系数尺度的均值/先验方差、收敛、迭代及所选起点。当前notebook在docs/current_comparison.rds记录结果目录、manifest路径和配置；校验脚本另写current_comparison.txt供报告生成器读取。

example_decisions.csv 包含一次实验的 pathway_id、state、is_nonnull 和每方法的 reject，供核查。完整每次决策与评分在 method_cache 中，而不是全部展开成巨大CSV。

provenance.rds 含原始manifest、方法指纹、项目R代码MD5；配置、包版本和sessionInfo另存文本/CSV。缓存文件名不替代审计记录。工程内部路径记录可按根目录重建；不要只拷贝单个结果缓存。
