# Note for Factor Investing

## 1.Basis

根据CAPM模型，资产的预期超额收益率由如下一元线性模型决定
$$
E[R_i]-R_f=\beta_i(E[R_M]-R_f)
$$
其中，$\beta_i=\frac{cov(R_i,R_M)}{var(R_M)}$ 

APT理论，在CAPM基础上构建了线性多因子定价模型，假设资产 $i$ 的预期超额收益由以下多元线性模型决定
$$
E[R^e_i]=\beta_i'\lambda
$$
其中，$E[R^e_i]$ 代表预期超额收益 $E[R_i]-R_f$ ，本书简称其为 **预期收益率**

异象
$$
E[R^e_i]=\alpha_i+\beta_i'\lambda
$$
如果 $\alpha_i$ 显著不为0，则资产 $i$ 是一个异象

## 2.Methodology

### 2.1 Portfolio Sorting

#### 2.1.1 Factor Mimicking Portfolio

投资组合需要满足以下两个条件：

1. 该投资组合仅在目标因子上有大于0的暴露，在其他因子上暴露为 0 
2. 在所有满足条件 1 的投资组合中，该投资组合的特质性风险（idiosyncratic risk）最小



##### Calculate idiosyncratic risk

在股票和因子暴露及特异风险的计算中，特质性风险指的是股票的个体特征或不可分解的风险，也称为特异风险（idiosyncratic risk）或非系统性风险（unsystematic risk）。特质性风险是由于股票自身特征或其他公司特定因素导致的风险，与市场整体风险和宏观经济因素无关。

计算特质性风险的方法主要包括两种，一种是使用资产定价模型（Asset Pricing Models）进行回归分析，另一种是使用协方差矩阵分解的方法。

1. 使用资产定价模型进行回归分析时，通常采用CAPM、Fama-French三因子模型或Carhart四因子模型等资产定价模型，将股票收益率与市场因子和其他因子进行回归分析。回归分析的残差项即为特质性风险，它代表了股票自身特征所导致的风险。
2. 另一种方法是使用协方差矩阵分解的方法，将整个投资组合的协方差矩阵分解为市场风险和特质性风险两部分。其中，特质性风险的计算方法为，将整个投资组合的协方差矩阵减去市场因子的协方差矩阵和其他因子的协方差矩阵，得到的即为特质性风险协方差矩阵。

以上两种方法都可以用来计算投资组合中的特质性风险，但使用资产定价模型进行回归分析需要选择合适的模型和数据，而使用协方差矩阵分解的方法需要估算协方差矩阵

#### 2.1.2 Single-sort

具体操作方法略。

投资组合排序检验最重要的是**检验因子预期收益率**，原假设通常为因子预期收益为 $0$ ，令 $\{\lambda_t \} (t=1,2,\dots,T)$ 代表因子收益率时间序列，则因子预期收益率的估计 $\hat{\lambda}$ 及标准误 $s.e.(\hat{\lambda})$ 为
$$
\begin{align*}
\hat{\lambda} =& \frac{1}{T}\sum_{t=1}^{T}\lambda_t \\
s.e.(\hat{\lambda})=&\frac{std(\lambda_t)}{\sqrt{T}}
\end{align*}
$$
可在原假设下计算 *t* 值，满足自由度为T-1的 t 分布
$$
t-value = \frac{\hat{\lambda}}{s.e.(\hat{\lambda})}
$$
排序法还关注依照排序变量高低得到的L个投资组合的收益率是否有较好的**单调性**，可以通过计算**收益率**和**排序变量分组**的**秩相关系数**来检验，最流行的为 $Spearman$ 秩相关系数
$$
\rho_s=\frac{cov(X_r,X_g)}{\sigma_{X_r}\sigma_{X_g}}
$$
其中，将L个投资组合的收益率高低排位记为 $X_r$ ，将他们依排序变量分组的高低排位记为 $X_g$

#### 2.2.3 Double-sort

##### Independent Double Sorting

操作方法略。

假设两个变量为 $X_1 、 X_2$ ，$L_1、L_2$分别为分组数，令 $R_{ij,t}$ 代表投资组合 $P_{ij}$ 第 $t$ 期的收益率，则因子（围绕$X_1$构建的因子投资组合）第 t 期的收益率 $\lambda_{X_1t}$ 为
$$
\lambda_{X_1t}=\frac{1}{L_2}\sum_{i=1}^{L_2}R_{L_1i,t}-\frac{1}{L_2}\sum_{i=1}^{L_2}R_{1i,t}
$$
同样，因子（围绕$X_2$构建的因子投资组合）第 t 期的收益率 $\lambda_{X_2t}$ 为
$$
\lambda_{X_2t}=\frac{1}{L_1}\sum_{i=1}^{L_1}R_{iL_2,t}-\frac{1}{L_1}\sum_{i=1}^{L_1}R_{i1,t}
$$
将上述两式种的投资组合收益率重新排列，可以得到如下等价形式
$$
\lambda_{X_1t}=\frac{1}{L_2}\sum_{i=1}^{L_2}(R_{L_1i,t}-R_{1i,t}) \\
\lambda_{X_2t}=\frac{1}{L_1}\sum_{i=1}^{L_1}(R_{iL_2,t}-R_{i1,t})
$$

##### Conditional Double Sorting

以$X_1$和$X_2$分别为第一、第二排序变量，在此方法中
$$
\begin{align*}
P_{L_2}^{top}=&P_{1L_2}\cup P_{2L_2} \cup \dots \cup P_{L_1L_2} \\
P_{L_2}^{bottom}=&P_{11}\cup P_{21} \cup \dots \cup P_{L_11}
\end{align*}
$$
得到此种方法下围绕变量 $X_2$ 构建的因子收益率
$$
\lambda_{X_2t}=R_{L_2}^{top} - R_{L_2}^{bottom}
$$

### 2.2 Multifactor Explanations of Asset Pricing Anomalies

资产预期（超额）收益和因子预期收益率之间满足如下关系
$$
E[R_i^e]=\alpha_i+\beta_i'\lambda
$$
多因子模型检验的三个部分：

1. 估计值： $\hat{\alpha}_i,\hat{\beta}_i,{\hat{{\lambda}}}$
2. 标准误：$\sigma(\hat{\alpha}_i),\sigma(\hat{\beta}_i),\sigma({\hat{{\lambda}}})$
3. 检验：
   - 联合检验所有N个资产的定价误差
   - 检验每个因子的预期收益率

多因子模型回归检验的三步：

1. 计算每个资产在所有因子上的暴露 $\beta_i$ ；
2. 通过回归分析对多因子模型进行估计；
3. 联合检验资产定价误差 $\alpha_i$ 及每个因子的预期收益率 $\lambda_k$ 。

#### 2.2.1 Time-Series Regression

##### Overview

需已知因子收益率，适合分析风格因子构成的多因子模型。

以因子收益率为自变量，资产超额收益率为因变量。

在时序上：
$$
R_{it}^{e}=\alpha_i+\beta'_i\lambda_t+\epsilon_{it} , \ t=1,2,\dots,T
$$
其中，$\lambda_t$ 表示t期因子收益率向量，$R_{it}^{e}$ 表示资产i 在t 期的超额收益

对每个资产 $i=1,2,\dots,N$ 用OLS估计模型，将$R_{it}^{e}$ 和$\lambda_t$在时序上取均值，可得
$$
E_T[R_i^e] = \hat{\alpha_i}+\hat{\beta_i'}\hat{\lambda} , \ i=1,2,\dots,N
$$
时间序列回归是对每个资产分别独立用多因子模型进行时序回归，不是最小化 $\hat{\alpha_i}$的平方和为目的求出的，区别于截面回归。

##### `GRS` test

假设 $\varepsilon_{it} \ I.I.D$，Gibbons et al.(1989)，检验 $\alpha_i$ 是否联合为 0

$H_0$ ：所有$\alpha_i$ 均为0

定义 $\hat{\alpha}=[\hat{\alpha_1}，\hat{\alpha_2}，\dots,\hat{\alpha_N}]'$ ，$\hat{\varepsilon_t}=[\hat{\varepsilon}_{1t},\hat{\varepsilon}_{2t},\dots,\hat{\varepsilon}_{Nt}]'$ 

$GRS$ 检验统计量：
$$
\frac{T-N-K}{N}(1+E[\lambda_t]'\hat{\Sigma}_\lambda^{-1}E[\lambda_t])^{-1}\hat{\alpha}'\hat{\Sigma}^{-1}\hat{\alpha} \sim F_{N,T-N-K}
$$
其中，
$$
\begin{align*}
\hat{\Sigma}_\lambda^{-1} =& \frac{1}{T}\sum_{t=1}^{T}[\lambda_t-E[\lambda_t]][\lambda_t-E[\lambda_t]]' \\
\hat{\Sigma}^{-1} =& \frac{1}{T}\sum_{t=1}^{T}\hat{\varepsilon}_t\hat{\varepsilon}_t'
\end{align*}
$$
PS，若误差项存在相关性或者异方差，ols出问题，此统计量也出问题，可用GMM来检验。

##### Reference

- Gibbons, M. R., Ross, S. A., & Shanken, J. (1989). A test of the efficiency of a given portfolio. *Econometrica: Journal of the Econometric Society*, 1121-1152.

#### 2.2.2 Cross-section Regression

##### Overview

适合GDP、CPI、利率等宏观经济因子，便于处理因子收益率未知的情况。截面回归考察 $R_{it}^{e}$ 和 $\beta_i$ 在截面上的关系。

- 1.先利用时序回归确定资产的因子暴露

假设$t$期一组因子的取值为 $f_t=[f_{1t},f_{2t},\dots,f_{Kt}]'$ ，K表示因子（特征）的数量，通过如下时间序列线性回归模型确定因子暴露
$$
R_{it}^{e}=a_i+\beta_i'f_t+\varepsilon_{it} , \ \ t=1,2,\dots,T,\forall i 
$$
注意，此处是$a_i$而非$\alpha_i$，因为解释变量不是因子收益率，那么截距项便不是定价误差。

使用OLS估计模型，得到资产的因子暴露**$\hat{\beta}_i$** 和残差$\hat{\varepsilon}_{it}$ 后，进入第二步

- 2.截面回归

使用第一步得到的因子暴露的估计$\hat{\beta}_i$作为解释变量，截面上满足线性回归：
$$
E_T[R_i^e] = \hat{\beta_i'}\hat{\lambda}+\alpha_i , \ i=1,2,\dots,N
$$
使用OLS估计，得到因子预期收益率的估计 $\hat{\lambda}$ 和每个资产的定价误差的估计 $\hat{\alpha}_i$ 

**$\star$** 此处，多因子模型假定不存在模型设定偏误时，资产的预期收益率 $E_T[R_i^e]$ 仅由因子暴露 $\hat{\beta_i'}$ 和因子预期收益率 $\hat{\lambda}$ 决定，因此无截距项，使用回归的残差直接作为定价误差 $\hat{\alpha}_i$ 

PS，Cochrane(2005)，指出也可以加截距项 $\gamma$ 

##### Reference

- Fama, E. F., & MacBeth, J. D. (1973). Risk, return, and equilibrium: Empirical tests. *Journal of political economy*, *81*(3), 607-636.
- Cochrane, J. (2009). *Asset pricing: Revised edition*. Princeton university press.



#### 2.2.3时间序列回归 vs 截面回归

#### 2.2.4Fama-MacBeth 回归

#### 2.2.5不同回归方法比较























