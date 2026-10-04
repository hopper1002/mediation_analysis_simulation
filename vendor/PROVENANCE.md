# MDACT 来源与修改

作者仓库：[asmita112358/MLFDR](https://github.com/asmita112358/MLFDR)。本次取自用户提供本地副本：

`D:\A_Mediation_analysis\code\MLFDR-main\results_in_paper\MDACT_code\funcs.R`

原始文件 SHA256：

`CC5931D8827A596E2F38C23931824A0B06D7FDA6A1D46283A4F25725CC6D0060`

`mdact_original.R` 是该文件的逐字副本，供审计，**不由 bootstrap 执行**。`mdact_core.R` 节选 `balancing_DACT_control_DR_adjust` 和 `DR_DACT_thr_adjust`，另包含等价解析CDF辅助函数，在独立环境中加载，不执行原文件的 `library(HDMT)` 或其余方法。

保护性修改：将两处 `while(inherits(x,"error"))` 改成 `while(inherits(x,"error") && k >= 0)`，失败后明确报错，原rounding顺序保留。当前比较另将F00的数值积分改为相同零假设事件的解析CDF，解决第一轮一次`the integral is probably divergent`数值错误；原失败记录保留在run_7c638fcb10f3。FDR统计量、权重下界和二分阈值搜索不变。未用BH或其他DACT方法替代。FWER和size分支保留，但只宣称验证FDR分支。

## F00解析计算

令a=π01、b=π10、c=π00，U/V独立Uniform(0,1)。原F00计算的是P(aU+bV+c max(U,V)²≤t)。对V=v条件化，分界r满足cr²+(a+b)r=t。

- v≤r时，允许的U上界为(-a+√(a²+4c(t−bv)))/(2c)，截在[0,1]；当v≤(t−c−a)/b时上界为1。
- v>r时，上界为(t−cv²−bv)/a，直到cv²+bv=t后变为0。

分别积分平方根和多项式得到解析面积；代码对二次方程使用稳定根，并将两个3/2次幂的差因式分解，减少尾部消去误差。t≤0返回0，t≥a+b+c返回1。权重的原1e−3下界保证分母非零。

tests/test_mdact_cdf.R验证100组随机权重、每组2个内部阈值（含小尾部），与按分界点拆分的独立数值积分一致到1e−9，并检查端点、有界和单调性。本修复改变数值计算方式，不改变F00的统计定义。

本地作者包 DESCRIPTION 声明 GPL-2，工程 LICENSE 保存 GPL-2 全文。作者/方法学归属仍属于原作者；新增工程模块也采用 GPL-2。使用/发表结果时应同时引用 MLFDR 论文和各比较方法原文。MDACT 对应论文参考文献30；HDMT 对应参考文献3，使用论文给出的正式引用信息。
