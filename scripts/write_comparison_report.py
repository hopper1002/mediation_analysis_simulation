"""Write a Chinese run report from validated result files; no result editing."""
from pathlib import Path
import csv
from collections import defaultdict

root = Path(__file__).resolve().parents[1]
pointer = dict(line.split("=", 1) for line in (root / "docs" / "current_comparison.txt").read_text().splitlines())
run = Path(pointer["run_dir"])
def read(name):
    with (run / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))
summary = read("summary.csv")
metrics = read("replicate_metrics.csv")
paired = read("paired_differences.csv")
fits = read("fit_diagnostics.csv")
fdr = [float(x["FDR_mean"]) for x in summary]
power = [float(x["power_mean"]) for x in summary]
flagged = [x for x in summary if float(x["FDR_lo"]) > float(x["q"])]
by_method = defaultdict(list)
for x in summary:
    by_method[x["method"]].append(x)

text = f"""# 论文模型比较：执行和核查记录

执行日期：2026-10-04。当前结果：`{run.name}`。
数据入口：`{pointer['manifest_path']}`。

## 设计和可复现范围

依据论文第2.2节式(3)–(5)、图2–4，保持三种模型、m=1000、n=300、稀疏/密集状态概率、X~Ber(0.1)、hi方差1/n、gi方差4/n、β均值−0.5τ、γ~N(1,0.5)和两个连续模型误差N(0,1)。唯一增强的DGP系数为α均值：0.05τ→0.35τ。τ取原网格1.3、1.5、1.9，每条件40次。它是**论文模型的信号增强比较**，不是原图严格复现。

试运行seed=20261005，n=100/300、τ=0.7/1.3/1.9、每条件4次，run_23dcc0251117完整保留。较弱信号评估seed20261006、统一哈希后的run_88f15ca0d1e2包含τ1.1/1.5/1.9，每条件40次；其τ1.1二元稀疏MLFDR经验FDR约0.1002，因此作为压力测试列在主notebook，全部数据/重复保留，notebook/HTML另归档03_weaker_signal_stress。最终主比较选τ≥1.3，独立验证基础seed=20261007预先固定，没有按排名搜索种子或删除重复。

最终使用v2序列化哈希，修复v3头部native_encoding使Windows命令行和IRkernel生成不同逐次种子的问题。修复前run_7c638fcb10f3和run_b4ee92372f64均保留，后者归档02_locale_variant_comparison。历史字面正文弱信号演示也独立保留。

## 实际结果

- {len({x['case_id'] for x in metrics})}份数据，{len(metrics)}次方法调用；错误{sum(x['status'] != 'ok' for x in metrics)}次。
- 有警告的调用{sum(bool(x['warning']) for x in metrics)}次，完整文本在diagnostics.csv；主要为强β导致p值舍入并列后的KS警告。
- MLFDR收敛{sum(x['converged'] == 'TRUE' for x in fits)}/{len(fits)}，迭代范围{min(int(x['iterations']) for x in fits)}–{max(int(x['iterations']) for x in fits)}。
- 各条件经验FDR均值范围{min(fdr):.4f}–{max(fdr):.4f}；功效均值范围{min(power):.4f}–{max(power):.4f}。
- FDR近似95%区间下限高于目标0.05的条件数：{len(flagged)}。这是校准提示，不能作为多个条件的联合保证。
- MLFDR减HDMT/MDACT的{len(paired)}项配对功效比较中，有{sum(float(x['delta_power_lo']) > 0 for x in paired)}项的近似95%区间下限为正；所有比较都输出，不隐藏其余项。

| 方法 | FDR均值范围 | 功效均值范围 |
|:--|--:|--:|
"""
for method in ("HDMT", "MDACT", "MLFDR"):
    rows = by_method[method]
    fs = [float(x["FDR_mean"]) for x in rows]
    ps = [float(x["power_mean"]) for x in rows]
    text += f"| {method} | {min(fs):.4f}–{max(fs):.4f} | {min(ps):.4f}–{max(ps):.4f} |\n"
text += "\n## 中间信号τ=1.5的实际对比\n\n| 场景 | 比例 | 方法 | FDR | 功效 | FDR MCSE | 功效 MCSE |\n|:--|:--|:--|--:|--:|--:|--:|\n"
for x in summary:
    if float(x["tau"]) == 1.5:
        text += f"| {x['scenario']} | {x['mixture']} | {x['method']} | {float(x['FDR_mean']):.4f} | {float(x['power_mean']):.4f} | {float(x['FDR_mcse']):.4f} | {float(x['power_mcse']):.4f} |\n"
text += """
## 方法实现和数值修复

HDMT使用已安装包原接口；MDACT保留作者的统计公式、权重下界和阈值搜索，但F00数值积分改为同一分布的解析CDF；200个阈值/权重组合与独立数值积分一致到1e−9。第一轮run_7c638fcb10f3有一次MDACT积分失败，该数据已按原始种子单独原样复算成功，记录在docs/audit/mdact_repair.csv，没有删除失败数据。最终完整比较还包含编码哈希修复，基础种子和参数均不因观测排名调整。

MLFDR默认后端paper_em按论文四成分模型拟合，使用√n缩放、log密度和3个固定初始权重，允许先验方差为0，仅按观测似然选起点。该后端是论文模型的工程实现，不能冒充MLFDR包原默认接口；package后端仍可通过参数选择。它没有使用模拟状态、真实混合比例或真实系数来拟合。拟合诊断完整输出。

## 验证

已逐份核对原始文件MD5、每个方法拒绝向量重新计算的FDP和功效、汇总、配对差异、EM收敛及当前源代码指纹。另有OLS/lm等价、确定性种子、不同Windows字符环境下哈希一致、串行/PSOCK一致、缓存隔离、异常处理、step-up边界、EM混合恢复/似然递增/尺度不变、MDACT解析CDF检查。notebook逐单元执行并导出HTML；三张主比较图直接嵌入notebook。

FDR为逐次FDP的均值；MCSE为重复间SD/√40；区间是近似蒙特卡洛均值区间，边界全0/全1采用保守端点界。不能保证每个有限样本条件精确等于0.05，也不能将功效增加单独当作优势而忽视误差率。论文图的效应强度和本实验不同，数值不要求相同。

## 重跑

参数集中于notebook首个代码单元。R函数在R/，新增方法在R/methods/plugins/。只改q或新增算法可复用原始数据；更改DGP参数生成独立数据目录。执行scripts/validate_comparison.R后再执行scripts/write_comparison_report.py可重新生成本记录。源代码、软件版本和配置保存在本次结果目录。
"""
(root / "docs" / "COMPARISON_RUN_REPORT.md").write_text(text, encoding="utf-8")
print(root / "docs" / "COMPARISON_RUN_REPORT.md")
