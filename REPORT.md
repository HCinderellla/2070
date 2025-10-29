# Portfolio Optimisation Summary / 投资组合优化总结

## Question 1 / 问题1
- **English:** Using daily simple returns from 1 July 2022 to 30 June 2023, the mean return vector and covariance matrix (rounded to six decimals) match the reference values in the brief, confirming our data preparation pipeline. 【F:analysis.py†L31-L68】【F:output.txt†L1-L28】
- **中文：** 基于2022年7月1日至2023年6月30日的日简单收益率计算得到的平均收益向量与协方差矩阵（保留六位小数）与任务描述中的校验值一致，证明数据处理流程正确。 【F:analysis.py†L31-L68】【F:output.txt†L1-L28】

## Question 2 / 问题2
- **English:** Solving the Markowitz problem with only the budget constraint yields linear weight formulas of the form \(w_i(t)=a_i t + b_i\); BHP and TLS have positive intercepts and slopes, so risk-averse investors (\(t>0\)) never short them, while CSL, NAB, and TCL turn negative once \(t\) exceeds their respective break-even values. 【F:analysis.py†L115-L134】【F:analysis.py†L531-L546】【F:output.txt†L30-L41】
- **中文：** 在仅含资金约束的马科维茨模型中，最优权重呈线性形式 \(w_i(t)=a_i t + b_i\)；其中 BHP 与 TLS 的截距和斜率皆为正，因此对于任何风险厌恶程度 \(t>0\) 的投资者都不会被做空，而 CSL、NAB、TCL 则在 \(t\) 超过各自阈值时会变为负权重。 【F:analysis.py†L115-L134】【F:analysis.py†L531-L546】【F:output.txt†L30-L41】

## Question 3 / 问题3
- **English:** For \(t=0.08\) and a \$1,000,000 portfolio, the optimal unconstrained allocation invests heavily in TLS and modestly in BHP, while taking short positions in NAB and TCL; the expected daily profit is \$767 with an \$8,229 daily standard deviation. 【F:analysis.py†L548-L554】【F:output.txt†L43-L53】
- **中文：** 当 \(t=0.08\) 且总资金为100万美元时，最优无限制组合大量配置 TLS 并适度持有 BHP，同时在 NAB 与 TCL 上建立空头头寸；日均收益约为 767 美元，日标准差约为 8,229 美元。 【F:analysis.py†L548-L554】【F:output.txt†L43-L53】

## Question 4 / 问题4
- **English:** The global minimum-variance portfolio (GMV) has a 0.000455 mean return with 0.006538 risk, and the generated SVG illustrates assets, the efficient frontier, and the Q3 optimum as required. 【F:analysis.py†L556-L562】【F:analysis.py†L345-L415】【F:figures/question4.svg†L1-L69】【F:output.txt†L55-L59】
- **中文：** 全局最小方差组合的平均收益为0.000455，风险为0.006538，生成的SVG图展示了各资产、有效前沿以及问题3的最优点，符合题目要求。 【F:analysis.py†L556-L562】【F:analysis.py†L345-L415】【F:figures/question4.svg†L1-L69】【F:output.txt†L55-L59】

## Question 5 / 问题5
- **English:** Enumerating active sets provides the no-short-sale optimum (BHP 28.7%, TLS 71.3%), and the computed KKT multipliers verify optimality; relative to Question 3, this feasible portfolio sacrifices return for lower risk. 【F:analysis.py†L171-L234】【F:analysis.py†L564-L581】【F:output.txt†L61-L77】
- **中文：** 通过枚举可能的活跃约束得到禁止卖空时的最优组合（BHP 占28.7%，TLS 占71.3%），计算出的 KKT 乘子验证了其最优性；与问题3相比，该组合降低了风险但也减少了期望收益。 【F:analysis.py†L171-L234】【F:analysis.py†L564-L581】【F:output.txt†L61-L77】

## Question 7 / 问题7
- **English:** Incorporating the 0.00015 risk-free rate leads the \(t=0.08\) investor to lend 42.9% at the risk-free rate and lever the tangency portfolio for the remainder; the SVG shows the capital market line and inferred market portfolio. 【F:analysis.py†L237-L264】【F:analysis.py†L583-L606】【F:figures/question7.svg†L1-L79】【F:output.txt†L79-L87】
- **中文：** 将无风险利率0.00015纳入模型后，\(t=0.08\) 的投资者会将42.9%的资金配置于无风险资产，并把余下资金投入切线组合；SVG 图给出了资本市场线与由协方差推导出的市场组合。 【F:analysis.py†L237-L264】【F:analysis.py†L583-L606】【F:figures/question7.svg†L1-L79】【F:output.txt†L79-L87】
