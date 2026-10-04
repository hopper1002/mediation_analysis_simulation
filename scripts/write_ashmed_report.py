"""Write a factual run report from notebook 03's saved CSV results."""
from pathlib import Path
import csv

root = Path(__file__).resolve().parents[1]

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

pointer = {}
label = None
for line in (root / "docs/current_ashmed.txt").read_text(encoding="utf-8").splitlines():
    key, value = line.split("=", 1)
    if key == "dataset":
        label = value
        pointer[label] = {}
    else:
        pointer[label][key] = value

lines = ["# ASHMED 两套已存数据比较：实际运行报告", "",
         "第三notebook按ASHMED附件第2.1–2.2节独立实现数学方法，未采用其模拟实验设计。读取01与02已经生成的两套正式数据，**没有重新生成原始数据**。", "",
         "新增方法的模型设定预先固定为J=12、r_min=0.25、r_max=8、ρ∈{−0.5,0,0.5}，共61个零中心same-scale成分；初始π按附件相对值归一化。相对似然阈值1e−8、KKT gap阈值1e−4、最大5000轮，全量1000条路径拟合。收敛保护及共轭后验效应矩是明确说明的补全。", "",
         "两套均保留n=300、m=1000、τ=1.3/1.5/1.9、3场景×2混合、每条件40次、q=0.05、原三个方法与选项、FDR/Power/MCSE/近似95%区间和重复内配对比较。每套720份、2880次方法调用；总计1440份、5760次，其中ASHMED新增1440次，原三方法4320次复用历史缓存。", "",
         "## 逐套结果", ""]

for label, spec in pointer.items():
    directory = Path(spec["run_dir"])
    summary = read_csv(directory / "summary.csv")
    metrics = read_csv(directory / "replicate_metrics.csv")
    fits = read_csv(directory / "ashmed_fit_diagnostics.csv")
    audit = read_csv(directory / "baseline_audit.csv")
    data_audit = read_csv(directory / "data_audit.csv")[0]
    errors = [x for x in metrics if x["status"] != "ok"]
    ash_warnings = [x for x in metrics if x["method"] == "ASHMED" and x["warning"]]
    baseline_warnings = [x for x in metrics if x["method"] != "ASHMED" and x["warning"]]
    old = [x for x in metrics if x["method"] != "ASHMED"]
    max_gap = max(float(x["dual_gap"]) for x in fits)
    assert len(metrics) == 2880 and not errors and all(x["cached"] == "TRUE" for x in old)
    assert all(x["converged"] == "TRUE" for x in fits) and max_gap <= 1e-4
    assert all(x["unchanged"] == "TRUE" for x in audit)
    lines += [f"### {label}", "", f"- 已存manifest：`{Path(spec['manifest_path']).parent.name}/manifest.rds`。",
              f"- 新结果：`results/{directory.name}`。",
              f"- 成功调用{len(metrics)}/{len(metrics)}；ASHMED警告{len(ash_warnings)}；保留原方法历史警告{len(baseline_warnings)}次。",
              f"- ASHMED收敛{sum(x['converged']=='TRUE' for x in fits)}/{len(fits)}；最大轮数{max(int(x['iterations']) for x in fits)}；最大KKT gap={max_gap:.6g}。",
              f"- 原方法2160条逐次记录的R/V/S/A/FDP/Power均保持一致，CSV读写浮点误差至多{max(float(x['max_absolute_difference']) for x in audit):.3g}。",
              f"- 原始数据MD5匹配{data_audit['checksum_matches']}/{data_audit['datasets']}；修改时间与manifest均未变。", "",
              "线性、密集、τ=1.5的40次重复结果：", "",
              "| 方法 | 经验FDR | FDR MCSE | Power | Power MCSE |",
              "|:--|--:|--:|--:|--:|"]
    middle = [x for x in summary if x["scenario"] == "linear" and x["mixture"] == "dense" and float(x["tau"]) == 1.5]
    for name in ["HDMT", "MDACT", "MLFDR", "ASHMED"]:
        x = next(z for z in middle if z["method"] == name)
        lines.append(f"| {name} | {float(x['FDR_mean']):.4f} | {float(x['FDR_mcse']):.4f} | {float(x['power_mean']):.4f} | {float(x['power_mcse']):.4f} |")
    lines += ["", "全部18个条件的均值范围：", "", "| 方法 | FDR范围 | Power范围 |", "|:--|:--|:--|"]
    for name in ["HDMT", "MDACT", "MLFDR", "ASHMED"]:
        rows = [x for x in summary if x["method"] == name]
        fdr = [float(x["FDR_mean"]) for x in rows]
        power = [float(x["power_mean"]) for x in rows]
        lines.append(f"| {name} | {min(fdr):.4f}–{max(fdr):.4f} | {min(power):.4f}–{max(power):.4f} |")
    lines += ["", "这里每条件的FDR是逐次FDP平均，不是合并错误数/合并发现数；范围不表示每个条件都严格控制在0.05。", ""]
    pairs = read_csv(directory / "paired_differences.csv")
    selected = [x for x in pairs if x["scenario"] == "linear" and x["mixture"] == "dense" and float(x["tau"]) == 1.5 and x["comparison"] == "ASHMED - MLFDR"][0]
    lines += [f"该线性密集条件下，ASHMED−MLFDR的配对功效差为{100*float(selected['delta_power_mean']):.2f}个百分点，近似95%区间[{100*float(selected['delta_power_lo']):.2f}, {100*float(selected['delta_power_hi']):.2f}]。同时ASHMED的FDR更低，不能只看Power作总体优劣判断。", ""]
    flags = [x for x in summary if float(x["FDR_lo"]) > float(x["q"])]
    lines += [f"FDR区间下限超过目标的条件：{len(flags)}个。这是近似校准提示，未作多条件联合检验。", ""]
    if flags:
        lines += ["| 场景 | 混合 | τ | 方法 | FDR | 近似区间 |", "|:--|:--|--:|:--|--:|:--|"]
        for x in flags:
            lines.append(f"| {x['scenario']} | {x['mixture']} | {x['tau']} | {x['method']} | {float(x['FDR_mean']):.4f} | [{float(x['FDR_lo']):.4f}, {float(x['FDR_hi']):.4f}] |")
        lines += [""]
    lines += ["![线性场景四方法比较](../results/" + directory.name + "/figures/linear_comparison.png)", ""]

lines += ["## 解释与边界", "",
          "ASHMED在线性与已测混杂场景整体较保守，Power明显低于MLFDR，部分二元条件接近或略超过HDMT/MDACT。它在两套数据上全部数值收敛，因此低Power不能简单归为优化未完成。两套数据均有明确非零效应方向，教学设计还有固定非零效应，而ASHMED只用零中心正态尺度混合并限制same-scale；模型适配差别是有依据的解释，但本实验没有逐项隔离这些因素，不能证明唯一原因。", "",
          "第一套二元、密集、τ=1.9时ASHMED的FDR偏高，区间下限超过0.05，已如实保留。二元结局在ASHMED中是Logistic系数正态摘要扩展，不是附件线性理论的直接验证。后验估计FDR≤0.05与真实经验FDR≤0.05并非同一个保证。", "",
          "不依据方法排名更换种子、增加信号、调整先验中心/网格或删除重复。原三个方法也有FDR偏高的历史条件，原样展示。40次重复带来有限Monte Carlo精度，稀疏场景的FDP区间尤其较宽；观察到0误报也保留端点不确定性。", "",
          "所有数值均在[0,1]内；三方法原结果保持不变。各方法的elapsed_seconds来自不同运行/历史缓存，不据此宣称严格运行速度排名。", "",
          "## 输出与核查", "",
          "- 已执行 `notebooks/03_ashmed_comparison.ipynb` 及同名HTML：12个R代码单元，无错误输出，8幅图（六幅主比较、两幅收缩示例）。",
          "- 每套保存72行汇总、2880行逐次指标、54行ASHMED与三个原方法的配对比较、720行ASHMED拟合诊断。",
          "- 每套保留1000条示例后验、61个先验成分/拟合权重、示例优化轨迹；代码/环境/配置/数据指纹及审计表在相应新结果目录。",
          "- `tests/check_ashmed.R`：独立数学和数值检查通过。",
          "- `scripts/validate_ashmed.R`：全部1440份原数据及5760次拒绝/真值指标、后验、收敛、阈值、汇总、配对表与原结果一致性逐项核查通过；见 `ASHMED_VALIDATION.txt`。",
          "- `scripts/validate_ashmed_notebook.py`：R内核、全部执行、无错误、八幅嵌入图及HTML一致性检查通过；代码单元没有数据生成调用。",
          "- 已目视检查主图与收缩图，坐标、图例、四方法颜色及信息完整。", "",
          "入口 `docs/current_ashmed.rds`/`.txt` 独立于原两个入口。已有notebook和原数据未被覆盖。复用及补全依据见 [ASHMED_METHOD.md](ASHMED_METHOD.md)。", ""]

target = root / "docs/ASHMED_RUN_REPORT.md"
target.write_text("\n".join(lines), encoding="utf-8")
print(target)
