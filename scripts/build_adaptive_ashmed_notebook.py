"""Author notebook 04. Rebuilding replaces its edits and executed outputs."""
from pathlib import Path
import nbformat as nbf

root = Path(__file__).resolve().parents[1]
cells = []
def md(s): cells.append(nbf.v4.new_markdown_cell(s.strip()))
def r(s): cells.append(nbf.v4.new_code_cell(s.strip()))

md(r"""
# 改进 ASHMED：方法开发、固定参数与两套数据的验证

原版 ASHMED 的零中心先验及共享尺度在当前两套数据上过于保守。本 notebook 保留它的**四状态贝叶斯收缩、混合先验、local FDR 筛选**，尝试学习非零效应中心、解除两段尺度绑定，并为复合零状态增加正则化。独立插件名为 **ASHMED_Adaptive**，不冒称论文原版或作者官方实现。完整推导见 `docs/ASHMED_ADAPTIVE_METHOD.md`。

**只读取前两个 notebook 已保存的数据，不重新生成。** 两套均有18个条件（线性/已测混杂/二元 × 稀疏/密集 × τ=1.3/1.5/1.9），n=300、m=1000、q=0.05。第一套非零效应有方向及随机波动，第二套有方向且效应固定；这也是本改进的适用范围。

本轮按重复编号分开：1–3测试结构，4–8选择正则强度，9–12确认候选，13–40用于最终比较（**每条件28次**）。最后一组未用于本轮选参数；它们此前已用于原版及对照方法，所以这是保存数据中的内部验证，不能称为独立外部验证。原版 ASHMED、HDMT、MDACT、MLFDR 与前三个 notebook 保留。
""")
r(r"""
find_project_root <- function(start=getwd()) {
  path <- normalizePath(start,winslash="/",mustWork=TRUE)
  for (i in 1:8) {
    if (file.exists(file.path(path,"R/bootstrap.R"))) return(path)
    candidate <- file.path(path,"mlfdr_simulation_lab")
    if (file.exists(file.path(candidate,"R/bootstrap.R"))) return(candidate)
    path <- dirname(path)
  }
  stop("请从工程目录或 notebooks 目录启动。")
}
PROJECT_ROOT <- find_project_root()
source(file.path(PROJECT_ROOT,"R/bootstrap.R"),encoding="UTF-8")
source(file.path(PROJECT_ROOT,"R/ashmed_optimization.R"),encoding="UTF-8")
options(repr.plot.width=10,repr.plot.height=7,repr.plot.res=145,
        repr.matrix.max.rows=90,repr.matrix.max.cols=18)
WORKERS <- 4L
campaign <- optimization_directory(PROJECT_ROOT)
lock <- readRDS(file.path(campaign,"locked_method.rds"))
SELECTED_OPTIONS <- lock$selection$options
registry <- load_method_registry(PROJECT_ROOT)
stopifnot(identical(method_fingerprint(registry$ASHMED_Adaptive,PROJECT_ROOT),lock$selection$fingerprint))
canonical <- environment(registry$ASHMED_Adaptive$run)$adaptive_options(SELECTED_OPTIONS)
canonical$base <- NULL
cat("选定方案：",lock$selection$arm,"；q=0.05；最终重复13:40\n")
print(canonical)
inputs <- optimization_inputs(PROJECT_ROOT)
IRdisplay::display(do.call(rbind,lapply(names(inputs),function(name) {
  x <- inputs[[name]]
  data.frame(dataset=name,data_id=x$manifest$data_id,saved_datasets=nrow(x$manifest$design),
    final_datasets=18*28,n=x$config$n,m=x$config$m,q=x$config$q)
})))
""")
md(r"""
## 方法改动与原因

设 $\hat\theta_i=(\hat\alpha_i,\hat\beta_i)^T\mid\theta_i\sim N_2(\theta_i,S_i)$，继续使用原 OLS/GLM 摘要与标准误。四状态 H00/H10/H01/H11 的定义不变，只有 H11 是真中介。

1. **学习方向与位置。** 在每个边际拟合“精确零点质量 + 自由均值正态”经验贝叶斯模型。用四个确定起点，按观测似然选结果，从摘要估计非零均值 $\mu_\alpha,\mu_\beta$；不读取真值。
2. **允许独立尺度及固定效应。** 围绕学习的均值放置6个正方差尺度和1个非零点质量，H11用两个尺度的笛卡尔积。典型模型共有 $1+7+7+49=64$ 成分，允许两段效应的离散程度不同。
3. **避免近零替代成分混淆。** 开发结果显示，保留大量零中心替代成分会抢走精确零状态的权重，导致 FDR 膨胀。选定版本只用有方向的非零中心替代字典。如果任一学习中心距零不足0.5倍中位标准误，则退回未修改的原版 ASHMED；不把近零点质量称作非零。
4. **对复合零状态正则化。** 最大化

\[
\sum_i\log\sum_kw_k f_{ik}+\frac{\lambda}{3}(\log\pi_{00}+\log\pi_{10}+\log\pi_{01}),
\quad \pi_s=\sum_{k:s(k)=s}w_k.
\]

这是状态权重的 Dirichlet 型 MAP 目标；EM 更新增加三个零状态的伪计数，保留非负、和为1的权重。固定字典的目标凹，优化同时检查目标单调性及 KKT gap。中心估计本身不保证全局最优。$\lambda$ 只在开发重复中选择，**q 始终0.05**，没有通过降低阈值或修改 DGP 调结果。

该思路参考 [Stephens (2017) 方法补充材料](https://stephenslab.uchicago.edu/assets/papers/Stephens2017-supplement.pdf) 中学习非零中心、通过混合比例惩罚稳定零成分的讨论；这里将其改为中介四状态模型，具体推导及实现差异见方法文档。
""")
md(r"""
## 完整开发记录：失败方案也显示

结构测试用每条件3次；随后只对更可识别的非零中心字典测试 λ=0/10/30/90，每条件5次。下表“平均”是18个固定条件等权平均；最大FDR是最差条件均值，不能用整体平均掩盖坏条件。

事先记录的**开发筛选门槛**为两套数据分别满足：整体平均 FDR≤0.06、最差条件 FDR≤0.15、Power高于原版、无方法错误。门槛仅用于少量开发重复筛选，不是0.05 FDR控制定理。合格方案中选平均Power最高者，再用9–12确认；确认失败时脚本会停止。最终28次比较不再据此改参数。
""")
r(r"""
structure_scores <- read.csv(file.path(campaign,"structure_scores.csv"),stringsAsFactors=FALSE)
development_scores <- read.csv(file.path(campaign,"development_scores.csv"),stringsAsFactors=FALSE)
confirmation_scores <- read.csv(file.path(campaign,"confirmation_scores.csv"),stringsAsFactors=FALSE)
cat("结构测试：\n"); IRdisplay::display(structure_scores)
cat("正则强度开发：\n"); IRdisplay::display(development_scores)
cat("独立于开发重复的候选确认：\n"); IRdisplay::display(confirmation_scores)
stopifnot(all(confirmation_scores$passes))
dev <- rbind(transform(structure_scores,phase="structure: B=3"),
             transform(development_scores,phase="penalty: B=5"))
dev$arm <- factor(dev$arm,levels=unique(dev$arm))
for (metric in c("mean_FDR","mean_power")) {
  g <- ggplot2::ggplot(dev,ggplot2::aes(arm,.data[[metric]],color=phase))+
    ggplot2::geom_point(size=3)+ggplot2::facet_wrap(~dataset,ncol=1)+
    ggplot2::theme_minimal(base_size=12)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="bottom")+
    ggplot2::labs(x="Candidate",y=metric,title=paste("Development |",metric))
  if (metric=="mean_FDR") g <- g+ggplot2::geom_hline(yintercept=.05,linetype=2)
  print(g)
  ggplot2::ggsave(file.path(campaign,paste0("development_",metric,".png")),g,
                   width=10,height=7,dpi=160,bg="white")
}
""")
md(r"""
## 固定参数后的最终四方法比较

读取13–40次保存数据，仅运行 ASHMED_Adaptive；HDMT、MDACT、MLFDR 及原版 ASHMED 从原实验提取**完全相同 case_id 的结果**。逐项核对原指标没有改变，检查720份原始文件/套的 MD5，另外检查原插件、三个 notebook、manifest 和原结果的保护指纹。

指标沿用原工程：每次 FDP=V/max(R,1)，Power=S/A；经验FDR是28次FDP平均，MCSE=SD/√28，阴影是同样的近似95%蒙特卡洛均值区间。并未改为合并V/合并R。整体表是18条件等权均值，区间按重复编号的平均块计算；条件图和表才是判断校准的主要依据。
""")
r(r"""
suite <- run_final_adaptive_comparison(PROJECT_ROOT,WORKERS)
for (name in names(suite)) {
  x <- suite[[name]]
  cat("\n",name,"：",x$directory,"\n")
  IRdisplay::display(x$macro)
  IRdisplay::display(x$data_audit)
  IRdisplay::display(x$reference_audit)
}
""")
for dataset, title in [("paper_adjusted","第一套：论文模型增强信号"),("teaching_fixed","第二套：固定效应教学数据")]:
    md(f"## {title}\n\n紫色为改进版，蓝/绿/红为HDMT/MDACT/MLFDR。左列FDR、右列Power，上排稀疏、下排密集；虚线q=0.05。")
    for scenario, label in [("linear","连续线性结局"),("confounded","已测混杂并调整Z"),("binary","二元结局：GLM正态摘要扩展")]:
        md(f"### {label}")
        r(f'''
x <- suite[["{dataset}"]]
print(plot_adaptive_comparison(x$main,"{scenario}"))
tab <- x$main[x$main$scenario=="{scenario}",c("mixture","tau","method","n_success","n_error",
  "FDR_mean","FDR_mcse","FDR_lo","FDR_hi","power_mean","power_mcse")]
numeric <- vapply(tab,is.numeric,logical(1)); tab[numeric] <- lapply(tab[numeric],round,4)
IRdisplay::display(tab)
''')
md(r"""
## 相对原版的改进，以及配对比较

原版与改进版在**同一28次重复**上比较；原版此前40次的均值不能直接拿来作这轮配对差。下面展示连续结局的原版/改进版曲线，然后列出 τ=1.5 的所有场景配对差。“改进版−对照”Power为正表示功效提升，FDR为正表示误报增加；应一起解释。单位为百分点，完整配对表保留所有τ。
""")
r(r"""
for (name in names(suite)) {
  x <- suite[[name]]
  original <- x$summary[x$summary$method %in% c("ASHMED","ASHMED_Adaptive"),]
  print(plot_adaptive_comparison(original,"linear"))
  tab <- x$paired[x$paired$tau==1.5,c("scenario","mixture","comparison","n_pairs",
    "delta_FDR_mean","delta_power_mean","delta_power_lo","delta_power_hi")]
  for (key in c("delta_FDR_mean","delta_power_mean","delta_power_lo","delta_power_hi"))
    tab[[key]] <- round(100*tab[[key]],2)
  cat("\n",name,"\n"); IRdisplay::display(tab)
}
""")
md(r"""
## 检查优化和FDR校准

固定字典的 KKT gap≤1e−4且相对目标变化≤1e−8才标记权重收敛；同时检查边际中心迭代。local FDR的选择均值≤q是算法行为，**不等于实际FDR必然≤q**。下面明确列出FDR均值超过0.05的条件，以及近似区间下界仍超过0.05的更明显偏高条件；不删除这些结果。
""")
r(r"""
for (name in names(suite)) {
  x <- suite[[name]]; f <- x$fits
  cat("\n",name,"拟合诊断\n")
  IRdisplay::display(data.frame(fits=nrow(f),weight_converged=sum(f$converged),
    location_converged=sum(f$location_converged),fallbacks=sum(f$fallback),
    max_gap=max(f$dual_gap),max_selected_lfdr=max(f$selected_mean_lfdr,na.rm=TRUE)))
  s <- x$main; s$clear_excess <- s$FDR_lo>.05
  IRdisplay::display(s[s$FDR_mean>.05,c("scenario","mixture","tau","method",
    "FDR_mean","FDR_lo","FDR_hi","clear_excess","power_mean")])
}
""")
md(r"""
## 一次拟合长什么样

预先取两套中“线性、密集、τ=1.5、第13次”作示例。均值和尺度只从估计系数学习；真值在拟合之后才加入此表。第一套β多为负、第二套β为正，改进版应由数据识别方向，而无需配置真实方向。总计64成分（其中H11为49成分）的字典完整权重、优化轨迹和1000条后验均有CSV。

后验采用非零中心正态共轭公式，产品均值包括分量内αβ协方差，产品SD包括混合不确定性。它们是经验贝叶斯插件后验矩，未包含超参数估计不确定性；不把均值±1.96SD当作精确可信区间。
""")
r(r"""
for (name in names(suite)) {
  x <- suite[[name]]; p <- x$example$posterior; d <- x$example$result$diagnostics
  cat("\n",name,"：",x$example$case_id,"\n")
  IRdisplay::display(data.frame(link=c("alpha","beta"),
    learned_mean=c(d$locations$alpha$mean,d$locations$beta$mean),
    fitted_alt_fraction=c(d$locations$alpha$pi1,d$locations$beta$pi1)))
  IRdisplay::display(data.frame(state=names(d$pi),fitted_pi=unname(d$pi)))
  IRdisplay::display(head(p[order(p$lfdr),c("pathway_id","truth_state","H00","H10","H01","H11",
    "lfdr","reject","alpha_hat","alpha_mean","beta_hat","beta_mean","mediation_mean","mediation_sd")],10))
}
""")
md(r"""
## 结论与扩展边界

本改进针对这两套**有方向、效应较集中**的数据。一个学习中心无法充分表示多峰或正负混合效应；稀疏弱信号时中心也可能不稳定。二元模型仍使用GLM正态摘要近似。有限条件的经验FDR、内部保留重复和KKT数值收敛都不构成一般FDR控制证明。应结合全部条件及不确定性判断，不能仅凭开发平均Power宣布优势。

最新逐条件结果和解读见 `docs/ASHMED_ADAPTIVE_RUN_REPORT.md`，完整核查见 `docs/ASHMED_ADAPTIVE_VALIDATION.txt`。后续方法仍是 `R/methods/plugins/` 插件；新方法只接收无真值摘要。原版 `ashmed.R` 保持不变，新增实现是 `ashmed_adaptive.R`，调度是 `R/ashmed_optimization.R`。

重跑最终比较：`Rscript --vanilla scripts/run_adaptive_ashmed.R`。重新执行/导出：`python scripts/execute_notebook.py notebooks/04_ashmed_adaptive.ipynb`。开发脚本分别为 `explore_adaptive_ashmed.R` 和 `tune_adaptive_ashmed.R`，会写开发记录及锁定配置；不要在看过最终结果后反复选参数仍把同一13–40次称为保留验证。若改动科学代码，指纹检查会要求重新开发/锁定；应安排新的验证数据或清楚标注探索结果。
""")
nb = nbf.v4.new_notebook(cells=cells,metadata={"kernelspec":{"display_name":"R","language":"R","name":"ir"},
  "language_info":{"name":"R","file_extension":".r","mimetype":"text/x-r-source"}})
path = root/"notebooks/04_ashmed_adaptive.ipynb"
nbf.write(nb,path)
print(path)
