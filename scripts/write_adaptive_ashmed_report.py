"""Write a factual report from completed evaluation CSV; attach it to notebook 04."""
from pathlib import Path
import csv
import nbformat
from nbconvert import HTMLExporter

root = Path(__file__).resolve().parents[1]
campaign = root / "results/ashmed_optimization_20261005"
def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))
def num(row, key): return float(row[key])
def table(rows, columns):
    lines = ["| " + " | ".join(label for _, label in columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        values = []
        for key, _ in columns:
            value = row[key]
            try: value = f"{float(value):.4f}"
            except (ValueError, TypeError): pass
            values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)

validation = root / "docs/ASHMED_ADAPTIVE_VALIDATION.txt"
assert validation.exists(), "Validate all final results before writing the report."
assert "validation passed" in validation.read_text(encoding="utf-8")
datasets = [("paper_adjusted", "第一套：论文模型增强信号"), ("teaching_fixed", "第二套：固定效应教学数据")]
dev = read(campaign/"development_scores.csv")
confirm = read(campaign/"confirmation_scores.csv")
selected = confirm[0]["arm"]
assert all(x["passes"] == "TRUE" and x["arm"] == selected for x in confirm)
report = ["# ASHMED_Adaptive 实际开发与验证结果", "",
  "日期：2026-10-05。原版插件与第三 notebook 保留；新增第四 notebook。以下数字全部来自实际运行的保存结果，没有重新生成数据。",
  "", "## 选中了什么", "",
  f"选定 **{selected}**：学习α/β非零中心、两个尺度的矩形字典、非零中心点质量、去掉零中心替代字典，加三个复合零状态的总伪计数λ={selected.removeprefix('null')}。典型64成分；任一学习中心与零距离不足0.5倍中位标准误时退回原版。q始终0.05。推导见 [ASHMED_ADAPTIVE_METHOD.md](ASHMED_ADAPTIVE_METHOD.md)。",
  "", "1–3次比较结构，4–8次选择λ，9–12次确认，13–40次固定参数最终比较。最后每条件28次、两套共1008份数据、5040个五方法结果；只新拟合改进版，原四方法按相同case_id提取。该划分是历史保存数据中的内部保留重复，原对照此前已在这些数据上运行，不能称为独立外部验证。",
  "", "## 开发尝试与失败记录", "",
  "仅放开零中心尺度导致两套整体FDR约0.38/0.35；增加学习中心但保留零中心替代成分也出现明显膨胀。近零替代字典与精确零状态难以识别，灵活性增加并不自动使FDR改善。结构开发记录保存在 `structure_scores.csv`，所有候选和源代码快照保留。去掉零中心替代字典后，再比较正则强度：", "",
  table(dev, [("arm","候选"),("dataset","数据"),("mean_FDR","平均FDR"),("max_FDR","最差条件FDR"),("mean_power","平均Power"),("gain","相比原版Power差"),("passes","开发门槛")]),
  "", "门槛为每套整体FDR≤0.06、最大条件FDR≤0.15、Power高于原版且无错误；合格者按平均Power选。门槛用于小样本开发，不是FDR控制保证。两套确认结果如下，确认通过后不再按13–40次结果改参数：", "",
  table(confirm, [("dataset","数据"),("mean_FDR","确认FDR"),("max_FDR","最差条件FDR"),("mean_power","确认Power"),("original_power","原版Power"),("passes","确认门槛")]),
  "", "## 固定方法后的结果", "",
  "每次FDP=V/max(R,1)，Power=S/A；FDR为每条件28次FDP平均。以下整体表为18个条件等权平均，MCSE按重复编号块计算；必须结合逐条件FDR判断，整体均值不代表每个条件都控制在0.05。"]
short = ["## 本轮实际结果解读", "",
  f"固定方案为 **{selected}（λ={selected.removeprefix('null')}）**。以下均为第13–40次，每条件28次的结果；开发与确认未混入最终均值。"]
for label, title in datasets:
    directory = campaign / f"validation_{label}"
    macro = read(directory/"macro_summary.csv")
    summaries = read(directory/"summary.csv")
    pairs = read(directory/"paired_differences.csv")
    fits = read(directory/"adaptive_fit_diagnostics.csv")
    by_method = {r["method"]:r for r in macro}
    adaptive = by_method["ASHMED_Adaptive"]; original = by_method["ASHMED"]
    new = [r for r in summaries if r["method"] == "ASHMED_Adaptive"]
    flags = [r for r in new if num(r,"FDR_lo") > .05]
    cols = [("method","方法"),("FDR_mean","整体FDR"),("FDR_mcse","FDR MCSE"),
      ("power_mean","整体Power"),("power_mcse","Power MCSE")]
    gain = 100*(num(adaptive,"power_mean")-num(original,"power_mean"))
    sentence = (f"{title}：原版平均Power {num(original,'power_mean'):.3f} → 改进版 {num(adaptive,'power_mean'):.3f}，"
      f"提高 **{gain:.1f}个百分点**；改进版整体FDR {num(adaptive,'FDR_mean'):.4f}，"
      f"逐条件均值范围 {min(num(r,'FDR_mean') for r in new):.4f}–{max(num(r,'FDR_mean') for r in new):.4f}。")
    report += ["", f"### {title}", "", table(macro, cols), "", sentence]
    short += ["", sentence, "", table(macro,[("method","方法"),("FDR_mean","整体FDR"),("power_mean","整体Power")])]
    worst = max(new,key=lambda row:num(row,"FDR_mean"))
    worst_sentence = (f"最高FDR出现在 {worst['scenario']} / {worst['mixture']} / τ={worst['tau']}："
      f"均值{num(worst,'FDR_mean'):.4f}，MCSE={num(worst,'FDR_mcse'):.4f}，"
      f"近似95%区间[{num(worst,'FDR_lo'):.4f}, {num(worst,'FDR_hi'):.4f}]。"
      "该条件需保留，不能据整体平均删除或忽略。")
    report += ["",worst_sentence]
    short += ["",worst_sentence]
    win = []
    for other in ("HDMT", "MDACT", "MLFDR", "ASHMED"):
        comp = f"ASHMED_Adaptive - {other}"
        rows = [r for r in pairs if r["comparison"] == comp]
        positive = sum(num(r,"delta_power_mean") > 0 for r in rows)
        clear = sum(num(r,"delta_power_lo") > 0 for r in rows)
        delta = 100*(num(adaptive,"power_mean")-num(by_method[other],"power_mean"))
        win.append(f"相比{other}整体Power差{delta:+.1f}个百分点；18条件中{positive}个均值更高，{clear}个配对近似95%区间下界>0。")
    report += ["", *["- "+s for s in win], "",
      "这些区间未做多重比较校正，只用于描述蒙特卡洛差异；不能把18个条件的区间当作研究性显著性筛选。",
      "", "预先指定示例条件：线性、密集、τ=1.5（28次均值）：", "",
      table([r for r in summaries if r["scenario"] == "linear" and r["mixture"] == "dense" and num(r,"tau") == 1.5],
        [("method","方法"),("FDR_mean","FDR"),("FDR_lo","FDR下界"),("FDR_hi","FDR上界"),("power_mean","Power")])]
    if flags:
        report += ["", f"**需关注：改进版有{len(flags)}个条件的FDR近似区间下界超过0.05。**", "",
          table(flags,[("scenario","场景"),("mixture","混合"),("tau","τ"),("FDR_mean","FDR"),("FDR_lo","下界"),("FDR_hi","上界"),("power_mean","Power")])]
        short += ["", f"**校准限制：这套改进版有{len(flags)}个条件的FDR近似区间下界>0.05，见上面的偏高条件表。**"]
    else:
        report += ["", "改进版没有条件的近似FDR区间下界超过0.05；部分点估计高于q。这不构成普遍控制证明，区间包含q也不等于已经证明控制。"]
        short += ["", "改进版没有条件的近似FDR区间下界超过0.05；部分点估计高于q，应保留不确定性。区间包含q不等于已证明控制。"]
    report += ["", f"504次新拟合全部权重及中心迭代收敛；最大KKT gap={max(num(r,'dual_gap') for r in fits):.6g}；近零中心回退{sum(r['fallback']=='TRUE' for r in fits)}次。完整逐条件表、配对表、诊断、示例后验及图形见 `results/ashmed_optimization_20261005/validation_{label}/`。"]
report += ["", "## 怎样解释这次改进", "",
  "原版将替代效应放在零中心，主要通过增大先验方差表达强信号；当前数据的效应有固定方向且较集中，这会浪费先验质量。学习非零中心后，小方差及非零点质量能表示更集中的效应，独立尺度也更适合两段信号不同的情况，因而识别效率上升。去掉近零替代字典与零状态正则化则抑制误把零路径归入H11。增加Power的同时是否控制FDR，仍必须由保存真值评估。",
  "", "改进版不是全面替代MLFDR；MLFDR本来就含可学习的有方向成分，当前网格中可有更高Power。HDMT/MDACT使用p值校准，而本方法借助跨路径效应分布，优势依赖此分布建模质量。当前方法只有一个非零中心，对混合正负方向、多峰、强依赖、极弱或更稀疏信号仍需另行测试。二元结局采用GLM正态摘要近似；经验贝叶斯后验没有计入学习中心、字典和权重的估计不确定性。",
  "", "## 核查与复现", "", "```text", validation.read_text(encoding="utf-8").strip(), "```",
  "", "数学单元检查 `tests/check_adaptive_ashmed.R` 已验证一般非零中心密度及共轭产品矩、已知惩罚最优解、目标单调性、原版等价开关与回退。notebook/HTML检查另外验证全部R单元已执行、无错误、10张图一致且无生成数据调用。",
  "", "最终重跑：`Rscript --vanilla scripts/run_adaptive_ashmed.R`；完整验证：`Rscript --vanilla scripts/validate_adaptive_ashmed.R`、`python scripts/validate_adaptive_ashmed_notebook.py`。源码/选项变化要求重新锁定，不能反复使用本轮验证结果选最优仍声称验证独立。"]
(root/"docs/ASHMED_ADAPTIVE_RUN_REPORT.md").write_text("\n".join(report)+"\n",encoding="utf-8")
short += ["", "改进版明显提高当前两套数据的识别效率，但没有在所有条件超过MLFDR。优势与有方向、集中效应的建模匹配有关，混合方向等场景需要新测试；全部失败尝试和校准限制保留。更详细解读、逐条件表及核查见 `docs/ASHMED_ADAPTIVE_RUN_REPORT.md`。"]
path = root/"notebooks/04_ashmed_adaptive.ipynb"
nb = nbformat.read(path,as_version=4)
for existing in nb.cells:
    if existing.cell_type == "markdown":
        existing.source = existing.source.replace("H11字典64成分的完整权重", "总计64成分（其中H11为49成分）的字典完整权重")
nb.cells = [c for c in nb.cells if c.get("metadata",{}).get("adaptive_result_interpretation") is not True]
cell = nbformat.v4.new_markdown_cell("\n".join(short))
cell.metadata["adaptive_result_interpretation"] = True
nb.cells.append(cell)
nbformat.write(nb,path)
body,_ = HTMLExporter().from_notebook_node(nb)
body = body.replace("</head>","<style>.jp-RenderedMarkdown td{overflow-wrap:anywhere}.jp-RenderedHTMLCommon table{font-size:13px}</style></head>")
path.with_suffix(".html").write_text(body,encoding="utf-8")
print("Wrote actual run report and attached interpretation to executed notebook/HTML.")
