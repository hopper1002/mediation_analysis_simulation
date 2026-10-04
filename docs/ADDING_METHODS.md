# 新增方法与统一接口

## 文件插件（推荐）

复制 `examples/new_method_template.R` 到 `R/methods/plugins/my_method.R`，实现 `register_plugin(registry)`。

```r
register_plugin <- function(registry) {
  run <- function(input, q, options = list()) {
    e <- input$estimates
    adjusted <- p.adjust(pmax(e$p_alpha, e$p_beta), "BH")
    list(reject = adjusted <= q, score = adjusted,
         score_type = "BH_adjusted_p", diagnostics = list())
  }
  register_method(registry, "MyMethod", run, packages = "stats",
                  version = "0.1.0", description = "方法及引用")
}
```

然后修改 notebook 参数区：

```r
config$methods <- c("HDMT", "MDACT", "MLFDR", "MyMethod")
config$method_options$MyMethod <- list()
registry <- load_method_registry(PROJECT_ROOT)
experiment <- run_experiment(manifest, config, PROJECT_ROOT, registry)
save_figures(experiment)
```

已有可运行示范 `MaxP_BH`。可以直接加入完整网格验证扩展；默认不参与主结果。

## 输入

`input$estimates`：data.frame，每行一个 pathway；顺序与原 M/Y 列顺序相同。列为 pathway_id、alpha_hat、beta_hat、var_alpha、var_beta、se_alpha、se_beta、p_alpha、p_beta。p 值处于 [1e-17,1]，方差有限且正。

`input$raw`：X (长度n)、Z（长度n或numeric(0)）、M (n×m)、Y (n×m)。允许需要原始样本的方法使用这些数据，但不可改动返回顺序。

`input$meta`：case_id、n、m、scenario，仅用于识别与选择模型。没有 truth、state、α/β 真值、π、τ、种子。不要在插件内部通过外部文件或全局变量绕过接口取真值，否则比较失去意义。

`q`：目标 FDR，0<q<1。`options`：本方法参数（例如 bandwidth、迭代次数），从 config$method_options[[方法名]] 获取。

## 输出

必需 `reject`：logical 长度 m，无 NA，按输入原顺序。可选 score：有限 numeric 长度 m；score_type 说明含义（p 值、local FDR、其他统计量）。可选 diagnostics：list 保存收敛信息、权重、模型参数等。不要求所有方法都产生 p 值或 local FDR。

无发现返回 `rep(FALSE,m)`。计算失败应 `stop()`，不要返回全 FALSE 或伪造零 FDR；调度器会隔离方法错误并写诊断。函数中产生的 warning 会被捕获和保存。

新方法自动参与逐次指标、汇总、图例和结果导出。新增时只需要修改插件与方法列表。插件文件哈希、函数体、适配器版本、依赖包版本、数据/估计指纹、q 和 options 控制缓存。辅助文件或外部模型未自动扫描，变化时需更新 version；不要使用不可记录的隐式全局选项。

## 并行与随机性

每方法每数据集单独设置确定性种子，执行后恢复原 RNG 状态。随机方法也能在 worker 数变化时获得同一设置。插件的函数与 registry 会传给 PSOCK worker，基础模块在各 worker 初始化。尽量用命名空间调用 `yourPackage::function`，并在 packages 声明包名；不要依赖 notebook 里未声明的临时变量、打开的连接或外部指针。调试设 workers=1。

## 修改后检查

先对一份已读取数据调用方法，再执行完整网格：

```r
registry <- load_method_registry(PROJECT_ROOT)
input <- make_method_input(example_data, estimate_coefficients(example_data)$estimates)
result <- registry$MyMethod$run(input, config$q, config$method_options$MyMethod)
invisible(validate_method_result(result, ncol(example_data$M)))
```

评估真值只在测试/评估阶段使用。增加至少一个与本方法数值语义相关的检查，避免只测“输出长度”。最终查看 diagnostics.csv 和每条件 n_success/n_error，不只看均值曲线。
