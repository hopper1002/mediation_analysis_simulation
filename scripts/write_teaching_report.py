"""Write a Chinese report from the actual validated teaching experiment."""
from pathlib import Path
import csv

root = Path(__file__).resolve().parents[1]
pointer = dict(line.split("=",1) for line in (root/"docs/current_teaching.txt").read_text().splitlines())
run = Path(pointer["run_dir"])
def rows(name):
    with (run/name).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
metrics, summary = rows("replicate_metrics.csv"), rows("summary.csv")
paired, fits = rows("paired_differences.csv"), rows("fit_diagnostics.csv")
cfg = {x["parameter"]:x["value"] for x in rows("teaching_design.csv")}
taus = sorted({float(x["tau"]) for x in summary})
middle = taus[(len(taus)-1)//2]
flags = [x for x in summary if float(x["FDR_lo"])>float(x["q"])]
worst = max(summary,key=lambda x:float(x["FDR_mean"]))
lines = ["# 简单中介路径实验：实际运行记录", "",f"当前结果：`{run.name}`。正式基础种子：{cfg['seed']}。", "",
         "## 生成机制与比较设置", "",
         f"固定比例分组，暴露比例{cfg['exposure_probability']}；固定四类路径数量并随机打乱顺序。非零α={cfg['alpha_mean_scale']}τ，非零β={cfg['beta_mean_scale']}τ，直接效应={cfg['direct_mean']}；连续模型M/Y误差标准差为{cfg['error_sd_m']}/{cfg['error_sd_y']}。混杂场景Z~N(0,1)，固定θ=δ={cfg['confounder_max']}并正确调整Z；二元场景使用logistic结局。",
         "",f"n={cfg['n']}、m={cfg['m']}、τ={cfg['tau']}、每条件{cfg['repetitions']}次、q={cfg['q']}。三种场景、稀疏/密集比例、方法选项、回归、FDR/Power/MCSE/区间、配对差异和主图沿用01 notebook；只改变DGP及其基础种子。开发小试使用20261008，结果完整保留，不与正式重复合并。", "",
         "## 执行和校准", "",
         f"- 数据份数：{len({x['case_id'] for x in metrics})}；方法调用：{len(metrics)}；失败：{sum(x['status']!='ok' for x in metrics)}。",
         f"- 有警告的调用：{sum(bool(x['warning']) for x in metrics)}；完整文本保存在diagnostics.csv。主要是HDMT/MDACT零比例估计中的KS并列p值提示，强β的p值在下界1e-17出现并列；记录为警告并结合实际FDR检查。",
         f"- EM收敛：{sum(x['converged']=='TRUE' for x in fits)}/{len(fits)}。",
         f"- FDR均值范围：{min(float(x['FDR_mean']) for x in summary):.4f}–{max(float(x['FDR_mean']) for x in summary):.4f}；Power均值范围：{min(float(x['power_mean']) for x in summary):.4f}–{max(float(x['power_mean']) for x in summary):.4f}。",
         f"- FDR区间下限高于目标的条件：{len(flags)}；配对功效区间下限为正：{sum(float(x['delta_power_lo'])>0 for x in paired)}/{len(paired)}。",
         "", "区间为相同的近似蒙特卡洛均值区间；不能将未发现超标当作所有条件严格控制的证明。实际结果未截断到q，未按方法排名筛选重复。", "",
         f"最高FDR出现在{worst['scenario']}、{worst['mixture']}、τ={worst['tau']}、{worst['method']}：均值{float(worst['FDR_mean']):.4f}，近似95%区间{float(worst['FDR_lo']):.4f}–{float(worst['FDR_hi']):.4f}。稀疏条件只有20条真中介，发现数少时一两条误报就会明显影响FDP；当前40次重复的精度不足以对弱条件宣称严格5%控制。", "",
         f"## 中间τ={middle:g}的结果", "", "| 场景 | 比例 | 方法 | FDR | Power | FDR MCSE | Power MCSE |", "|:--|:--|:--|--:|--:|--:|--:|"]
for x in summary:
    if float(x['tau'])==middle:
        lines.append(f"| {x['scenario']} | {x['mixture']} | {x['method']} | {float(x['FDR_mean']):.4f} | {float(x['power_mean']):.4f} | {float(x['FDR_mcse']):.4f} | {float(x['power_mcse']):.4f} |")
lines += ["", "## 如何理解", "",
          "四类路径中只有H11是真中介。另一段无效时，即使单段显著也属于复合零假设。新DGP的α、β固定且同向，MLFDR允许零先验方差的混合模型能利用效应方向及两段的联合密度；HDMT/MDACT使用p值及零假设混合。功效差异需结合FDR判断；结果只适用于所展示的分布。新旧数据的效应分布和分组比例不同，跨notebook的Power变化不是方法实现变好的证据。", "",
          "## 文件和复核", "",
          "`notebooks/02_teaching_mediation.ipynb`和同名HTML包含五张图：α/β的直观路径示例，以及三个场景的FDR/Power主图。结果目录含完整逐次指标、汇总、配对差异、拟合诊断、每次拒绝向量及配置。数据仍在本地data目录，生成代码可重建；原notebook入口独立保留。", "",
          "已逐份核对MD5、分组和路径数量、固定效应、真实拒绝向量重算的发现数/FDP/Power、汇总、配对差异、EM收敛、源码指纹、图形及Rscript/IRkernel的种子和缓存一致性。另检查固定路径数据的OLS/GLM与标准R回归一致、RNG隔离、无真值输入，以及旧版基础检查。", "",
          "重跑：`Rscript --vanilla scripts/run_experiment.R teaching`；核查：`Rscript --vanilla scripts/validate_teaching.R`；随后运行`python scripts/write_teaching_report.py`更新本记录。"]
if flags:
    lines += ["", "## 实际校准提示", "", "| 场景 | 比例 | τ | 方法 | FDR | 区间下限 | 区间上限 |", "|:--|:--|--:|:--|--:|--:|--:|"]
    for x in flags:
        lines.append(f"| {x['scenario']} | {x['mixture']} | {x['tau']} | {x['method']} | {float(x['FDR_mean']):.4f} | {float(x['FDR_lo']):.4f} | {float(x['FDR_hi']):.4f} |")
path = root/"docs/TEACHING_RUN_REPORT.md"
path.write_text("\n".join(lines)+"\n",encoding="utf-8")
print(path)
