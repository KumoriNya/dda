<h1 class="title" style="text-align:center; font-weight:bold">Statistical Analysis of Dependant Data</h1>

1. [Multivariate Normal Distribution](#mnd)
2. [Time Series](#ts)
3. [Spatial Statistics](#ss)

# Questions for the exam
in case we need to apply any algorithm, would their definition be attached somewhere? (Like DL, innovations)

if model selection shows up in tut -> make (written) notes
Or check whether model selection would show up in the exam

# 2. Multivariate Normal Distribution (Lecture 2) <a name="mnd"></a>

## 2.2. Random Vectors (Lecture 2A)

### 2.2.2. Independence

The random variables $X_1, \dots, X_n$ are independent if their joint density function $F_X(x_1, ..., x_n) = F_{x_1}(x_1)···F_{x_n}(x_n)$ is equal to the product of their individual density function.

### 2.2.5. Covariance

$Cov(X, Y)_{ij}:=Cov(X_i, X_j)$  
$=E[(X_i-E[X_i])(Y_j-E[Y_j])]$

## 2.3. Multivariate Normal Model

### 2.3.2. Multivariate Normal

Definition 2.6: A d-dimensional random vector $X$ is called a **multivariate normal random vector** if $X = \mu + A^{T}Z$, where $Z$ is the n-dimensional standard normal random vector.

### 2.3.3. Multivariate Moment Generating Functions

Definition 2.7: The **mgf** of a d-dimensional random vector $X$ is $m_X(u)=E(e^{u^{T}X})$.

### 2.3.4. Degenerate and Nondegenerate.

Definition 2.9: If a d-dimensional normal random vector takes values only in a subspace of dimension smaller than d, it is **degenerate**. If it takes values in all of $\R^d$, it is called non-degenerate.

Theorem 2.7: A multivariate normal random vector is non-degenerate iff its covariance matrix is poeitive definite.

* A non-negative definite matrix has an inverse (i.e. the matrix is non-singular) iff it is positive definite, which has a positve determinant.

### 2.3.6. Conditional Distributions

Theorem 2.9: (In Lecture 2A) Suppose we have: 

$(X_1, X_2)^{T} \sim N_{d_1+d_2}( (\mu_1, \mu_2)^{T}, ( (\Sigma_{11}, \Sigma_{21}) ^{T} , (\Sigma_{12}, \Sigma_{22}) ^{T} ) )$

where $\Sigma_{ij} = Cov(X_i, X_j)$ and $i, j \in {1, 2}$, and where $|\Sigma_{22}| > 0$. Then the conditional distribution of $X_1$ given $X_2 = x_2$ is 

$N_{d_1}(\mu_1+\Sigma_{12}\Sigma_{22}^{-1}(x_2-\mu_2), \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21})$.

## 2.4. The Multivariate Random Sample (Lecture 2B)

Definition 2.11: the **sample mean** $\bar{X} = \frac{1}{n}\sum_{i=1}^{n}X_i$, and the **sample covariance matrix** $S = \frac{1}{n-1}\sum_{i=1}^{n}(X_i-\bar{X})(X_i-\bar{X})^{T}$.

Theorem 2.10: The sample mean vector and the sample covariance matrix are both **unbiased** estimators for $\mu$ and $\Sigma$. Proof in the lecture slides.

### 2.4.4. Statistical Distance

Definition 2.13: The **statistical distance** from its mean is the length of the vector $\Sigma^{-1/2}(X-\mu)$.

Corollary 2.13: Let $X_1, \dots, X_n$ be a sample for any $d$-dimensional multivariate distribution with mean vector $\mu$ and covariance matrix $\Sigma$. Then as the sample size $n$ tends to $\inf$, the random variable $n(\bar{X}-\mu)^{T}\Sigma^{-1}(\bar{X}-\mu)$ converges in distribution to chi-squared distribution with d dof.

* Chi-squared with dof = d is $\sum_{i=1}^{d}Z_i^2$, where $Z_i$'s are standard normal.

## 2.5. Multivariate Normal Samples

### 2.5.1. Confirming Normality

Definition 2.14: Let $X_1, X_2, \dots, X_n$ be a univariate sample drawn from the distribution with distribution function F. Write $X_{(1)} < X_{(2)} < \dots < X_{(n)}$ for the same observations in increasing order.

The **empirical distribution** has distribution function $F_n(x) =$ #{$i : X_i \leq x$}$ / n$. The $\frac{i}{n}$ **sample quantile** is the number $X_{(i)}$.

### 2.5.2. Maximum Likelihood Estimation

Definition 2.15: The **likelihood function** L for the multivariate normal distribution is the function $L(\mu, \Sigma; x_1, \dots, x_n) = \prod_{i=1}^{n}f_X(x_i)$. If the estimator ($\hat{\mu}, \hat{\Sigma}$) satisfies $L(\hat{\mu}, \hat{\Sigma}; x_1, \dots, x_n) \leq L(\mu, \Sigma; x_1, \dots, x_n)$ for all $\mu \in \R^d$ and all positive definite $d \times d$ matrices $\Sigma$, then it is the **Maximum Likelihood Estimator** for the multivariate normal.

Theorem 2.15: If $B$ is a $d \times d$ positive definite matrix and $b \gt 0$, then $\frac{1}{|\Sigma|^b}e^{-\frac{1}{2}tr(\Sigma^{-1}B)} \leq \frac{1}{B^b}(2b)^{db}e^{-db}$ for all positive definite $d \times d$ matrices $\Sigma$ where equality holds only for $\Sigma = \frac{1}{2b}B$.

Theorem 2.16: For a multivariate normal distribution, $\hat{\mu} = \bar{X}, \hat{\Sigma} = \frac{n-1}{n}S$ are the unique maximum likelihood estimators.

### 2.5.3. Sampling distribution of $\bar{X}$ and $S$

Definition 2.16: Let $Y_1, \dots, Y_n$ be a random sample of length $n$ taken from the distribution $N_d(0, \Sigma)$, then $A = \sum_{i=1}^{n}Y_iY_i^{T}$ is **Wishart** distributed with covariance matrix $\Sigma$ and $n$ degrees of freedom. When $d = 1$ and $\Sigma = 1$, $A \sim \chi_n^2$.

Theorem 2.17: For multivariate normal samples with $(\mu, \Sigma)$, 
1. $\bar{X} \sim N_d(\mu, n^{-1}\Sigma)$,
2. $(n-1)S \sim W_{n-1}(\Sigma)$, 
3. $\bar{X}$ and $S$ are independent.

Proof: Require Definition 2.17 on **data matrix** and Theorem 2.18.

## 2.6. Inference for a Multivariate Normal Mean

### 2.6.1. Univariate data, $\sigma$ known. Test $\mu = \mu_0$:

By CLT, $\frac{\bar{X}-\mu}{\sigma / \sqrt{n}} \sim N(0,1)$, $n(\bar{X}-\mu)(\sigma^2)^{-1}(\bar{X}-\mu) \sim \chi_1^2$.

Test statistic: $U^2 = (\bar{X}-\mu_0)(\sigma^2)^{-1}(\bar{X}-\mu_0)$, so that $U^2 \sim \chi_1^2$ when $H_0$ is true.

For a size $\alpha$ test of $H_0: \mu = \mu_0 \text{ vs. } H_1: \mu \neq \mu_0$, reject $H_0$ if $U^2 > \chi_{1,\alpha}^2$.

### 2.6.2. Multivariate data, $\Sigma$ known.

$n(\bar{X}-\mu)(\Sigma)^{-1}(\bar{X}-\mu) \sim \chi_d^2$

### 2.6.3. Univariate data, $\sigma$ unknown.

$\frac{\bar{X}-\mu}{s / \sqrt{n}} \sim t_{n-1}$, 

Test statistic: $T^2 = n(\bar{X}-\mu_0)(S^2)^{-1}(\bar{X}-\mu_0) \sim F_{1, n-1, \alpha}$.

### 2.6.4. Multivariate, $\Sigma$ unknown => Hotelling's $T^2\text{- statistic}$

If $H = n\mathbf{Z}^TA^{-1}\mathbf{Z}$, where $\mathbf{Z}$ is standard $d$-variate normal vector, and $A \sim W_n(I_d)$ is independently Wishard distributed with covariance matrix unity $I_d$ and $n$ degrees of freedom where $n \leq d$, $H$ is the **Hotelling's $T^2$-distribution** with $n$ degrees of freedom and dimensionality parameter $d$.

Remark 2.2: $\frac{n-d+1}{nd}H \sim F_{d, n-d+1}$.

Recall that a random variable $V$ is $F_{d1, d2}$ distributed if $V = \frac{U_1/d_1}{U_2/d_2}$, where $U_1 \sim \chi_{d1}^2, U_2 \sim \chi_{d2}^2$, and $U_1$ and $U_2$ are independent. A proof is in *T>W> Anderson, "An Introduction to Multivariate Statistical Analysis"*.

Theorem 2.19: $T^2=n(\bar{X}-\mu)^{T}(S)^{-1}(\bar{X}-\mu)$ follows Hotelling's T-squared with $n-1$ dof

-> $\frac{n-d}{(n-1)d}H \sim F_{d, n-d}$

### 2.6.5. The Multivariate $t$-test

Result: Reject $H_0$ when $P(T^2 \geq \frac{(n-1)d}{nd}F_{d, n-d} | \mathbf{\mu} = \mu_0)$

### 2.6.6. Confidence Regions

Since $\frac{\bar{X}-\mu}{s / \sqrt{n}} \sim t_{n-1}$, its rejection region is $\frac{\bar{X}-\mu}{s / \sqrt{n}} \gt t_{n-1, \alpha/2}$, the CI is then $\bar{X}-\frac{S}{\sqrt{n}}t_{n-1,\alpha/2} \lt \mu_0 \lt \bar{X}+\frac{S}{\sqrt{n}}t_{n-1,\alpha/2}$ .

For Multivariate Normal Mean, $R = \{\mu : n(\bar{X}-\mu)^{T}S^{-1}(\bar{X}-\mu) \leq \kappa^2 = T^2\}$

### 2.6.7. Simultaneous Confidence Statements

- Uniformity over all linear combinations

Result on page 21 of Lecture 2C.

### 2.6.8. The Bonferroni Method

### 2.6.9. Large samples from non-normal distributions

# 3. Time Series (Lecture 3) <a name="ts"></a>

## 3.1. Introduction

### 3.1.1. Basic Definitions

Definition 3.19: Univariate **time series model**s.

Definition 3.20: The **Mean Function** of a time series $X$ is $\mu_X(t) = \mathbb{E}(X_t), t \in \Z$.

Definition 3.21: The **Covariance Function** of $X$ is $\gamma_X(r, s) = Cov(X_r, X_s) = \mathbb{E}[X_r-\mu_X(r)][X_s-\mu_X(s)], r,s \in \Z$.

Definition 3.22: A time series $X$ is **weakly stationary** if the mean and covariance function satisfy $\mu_X(t) = \mu_X \in \R, \gamma_X(t, t+h) = \gamma_X(h) \forall t, h \in \Z$.

The $\gamma_X(h)$ is called the **autocovariance function (ACVF)** of $X$, and the **autocorrelation function (ACF)** of $X$ is $\rho_X(h) = \frac{\gamma_X(h)}{\gamma_X(0)}, h \in \Z$.

Definition 3.23: **Strict stationarity** if $(X_{t_1}, X_{t_2}, \dots, X_{t_n}) =^d (X_{t_1 + h}, \dots, X_{t_n+h})$ holds for all $n \in \N$ and $h \in \Z$.

### 3.1.2. Differencing

The **backshift operator** $B$ takes as input a time series and produces as output the series $BX$ which is shifted backwards in time by one time unit. $BX_t = X_{t-1} \forall t \in \Z$.

The **difference operator** $\nabla$ is defined as $\nabla = (1 - B)$.

The **seasonal difference operator** $\nabla_S$: $\nabla_S = (1 - B^S)$.

### 3.1.3. Linear filters

Definition 3.24: Let $\{\psi_j\}_{j \in \Z}$ be an absolutely summable sequence of constants: $\sum_{j\in\Z}|\psi_j| \lt \inf$. The operator $\psi(B) = \sum_{j \in \Z}\psi_jB^j$ is called a **linear filter**. For a time series $X$, the **filtered time series** is $\psi(B)X$, given by $\psi(B)X_t = \sum_{j\in\Z}\psi_jX_{t-j}, t \in \Z$.

Definition 3.25: A series $X$ is a **general linear process** if it can be represented as $X = \psi(B)Z$ (10), where $Z \sim WN(0, \sigma^2)$.

In the case of ARMA(p, q), 

## 3.2. Trends and Seasonal Components

### 3.2.1. Eliminiation of trend in the **absence** of seasonality

Method 1: **Least squares estimation of** $m_t$. (Page 37 of Lecture 3A)

Method 2: Smoothing via a moving average filter. (Page 41)

Method 3: Differencing to obtain stationary data (Page 53)

### 3.2.2. Elimination of both trend and seasonality

Method 1: The small trend method (Page 59)

Method 2: Moving average estimation (Page 71)

Method 3: Differencing at lag $d$.

### 3.2.3. Design of filters **EXAMPLES**

## 3.3. ARMA(p, q) processes

#### Definition 3.26: White Noise Process

White noise process with variance $\sigma^2$ is a weakly stationary time series and acvf evaluated to variance at lag zero, else zero.

### 3.3.1. AR Process

#### Definition 3.27: AR($p$) Process

$X_t = \phi_1X_{t-1} + \phi_2X_{t-2} + \dots + \phi_pX_{p} + Z_t$

=> $Z_t = \phi(B)X_t$, 

where the **AR operator** is defined to be $\phi(B) = 1 - \phi_1B - \phi_2B^2 - \dots \phi_pB^p$.

### 3.3.2. MA Models

#### Definition 3.28: MA($q$) Model

$X_t = Z_t + \theta_1Z_{t-1} + \theta_2Z_{t-2} + \dots + \theta_qZ_{t-q}$.

=> $X_t = \theta(B)Z_t$,

where the **MA operator** is $\theta(B) = 1 + \theta_1B + \theta_2B^2 + \dots + \theta_qB^q$.

### 3.3.3. The ARMA(p, q) Process

#### Definition 3.29: ARMA($p, q$)

$\phi(B)X = \theta(B)Z$.

#### Theorem 3.20: 

### 3.3.4. Causal and Invertible

#### Definition 3.30: Causality

$X = \psi(B)Z$ for an absolutely summable sequence $\psi$ satisfying $\psi_j = 0$ for $j < 0$.

#### Theorem 3.21: AR Causality

ARMA($p, q$) is **causal** iff the AR polynomial $\phi(z) \neq 0 \forall z \in \mathbb{C}$ with $|z| \leq 1$.

i.e. Roots of the AR polynomial all have absolute value greater than 1.

The sequence $\psi$ is determined by:

$\psi(z) = \sum_{j = 0}^{\inf} \psi_j z^j = \frac{\theta(z)}{\psi(z)}$,

or equivalently

$(1 - \phi_1z - \dots - \phi_pz^p) ( \psi_0 + \psi_1z + \dots) = 1 + \theta_1z + \dots + \theta_qz^q$.

#### Definition 3.31: Invertibility

**Invertible** if $Z = \pi(B)X$ for an absolutely summable sequence $\pi$ satisfying $\pi_j = 0$ for $j < 0$.

#### Theorem 3.22: MA Invertibility

Similar to 3.21. Roots to $\theta(z)$ lie outside of the unit circle.

### 3.3.4. ACF of ARMA($p,q$) Processes

If causal, then $X_t = \sum_{j=0}^{\inf}psi_jZ_{t-j}$. 

Since white noises are uncorrelated, 

$\gamma_X(j) = E(X_{t+h}X_t) = \sigma^2\sum_{j=0}^{\inf}\psi_j\psi_{j+|h|}$

## 3.4. Estimation of Means and Covariances

A stationary process $X$ is characterised (from a second-order point of view) by its mean $\mu_X$ and ACVF $\gamma_X$, the estimation of these parameters and the ACF $\rho_X(h) = \gamma(h) / \gamma(0)$ is thus crucial for inference.

### 3.4.1. Sample mean of a stationary time series

The sample mean 

$\bar{X_n} = \frac{1}{n}(X_1 + \dots + X_n)$ 

is an unbiased estimator of $\mu$. Its variance is 

$Cov(\bar{X_n}) = \frac{1}{n}\sum_{h=-n}^{n}(1-\frac{|h|}{n})\gamma(h)$

If the infinite sum for $h$ is absolutely summable as n approaches infinity, then the sum approaches sum of just ACVF for all $h$, and the covariance of the sample mean approaches the sum of all covariance divided by $n$.

### 3.4.2. Distribution of the sample mean

If $\{X_t\}$ is Gaussian, the sample mean then is exactly normally distributed

$\sqrt{n}(\bar{X_n} - \mu) \sim N( 0 , \sum_{h=-n}^{n} (1- \frac{|h|}{n})\gamma(h) )$.

### 3.4.3. Sample ACV and AC

The usual estimator of the ACVF is 

$\hat{\gamma}(-h) = \hat{\gamma}(h) = \frac{1}{n} \sum_{t=1}^{n-h}(x_t - \bar{x})(x_{t+h}-\bar{x}), h \geq 0$

and the estimator of the ACF is

$\hat{\rho}(h) = \hat{\gamma}(h) / \hat{\gamma}(0)$.

### 3.4.4. Distribution of the sample ACF

**Bartlett's formula** on Page 71 of lecture 3B

## 3.5. Forecasting Stationary Time Series (Lecture 3C)

### 3.5.1. The Best Linear Predictor

#### Definition 3.32: BLP

$P_n(X_{n+h}) := a_0 + a_1X_n + \dots + a_nX_1$ is the **best linear predictor (BLP)** of $X_{n+h}$ if the square error is minimal for all $i$ from $0$ to $n$.

We eventually obtain the **Yule-Walker** equations. (Via multiplying $X$ to $X_{p}$ on both sides of the equation then taking expectation)

### 3.5.2. General setting for Yule-Walker equations

$\Gamma \mathbf{a} = \gamma$, where $\gamma = Cov(\mathbf{W}, \mathbf{W})$.

### 3.5.3. The Durbin-Levinson algorithm

The D-L recursion uses the one-step predictor to simplify the calculation of the upcoming predictor in the next step.

#### Theorem 3.23: Durbin-Levinson

### 3.5.4. The innovations algorithm

Is also a recursive algorithm that only requires finite second moments and not necessarily stationarity.

#### Theorem 3.24: Innovations algorithm

### 3.5.6. Forecasting ARMA($p, q$) Proceses

## 3.6. Parameter Estimation (Lecture 3D)

We need to estimate the mean and the white noise variance;

Estimate the parameter values for $\psi$ and $\theta$;

Choose the **correct** model.

### 3.6.1. Preliniary estimation

### 3.6.2. Method 1: Yule-Walker estimation

Results displayed with equation (30) and (32)

### 3.6.3. Method 2: The Partial ACF

#### Theorem 3.25: PACF

### 3.6.5. Innovations estimation

### 3.6.7. Maximum Likelihood Estimation

## 3.7. Model Selection

### 3.7.1. Large sample distribution of MLE's

Results on equation 36

### 3.7.2. Final Prediction Error

Results on page 58 of Lecture 3D

### 3.7.3. Akaike Information Criterion

Obtained via logging the FPE, and apply approximation of non-sigma terms , then multiplying by $n$.

# 4. Spatial Statistics (Lecture 4) <a name="ss"></a>

## 4.1. Random Fields (Lecture 4A)

### 4.1.1. Definitions

#### Definition 4.35: Random field

**Covariance function** of a spatial process: $C(s, h)$

### 4.1.2. Stationarity

#### Weak stationarity

#### Strict stationarity

## 4.2. Semivariogram

#### Intrinsic Stationarity

#### Definition 4.37: The **Semivariogram**

is defined as $\gamma(h) := \frac{1}{2}Var(Z(s+h)-Z(s))$ for any **intrinsically stationary** random field $Z$.

(Because $Var(Z(s+h) - Z(s)) = 2 C(0) - 2 C(h)$ )

### The Nugget Effect

Those semivariograms that do **not** pass through the origin.

## 4.3. Useful Semivariograms (Lecture 4B)

### 4.3.1. Isotropy and Anisotropy

Let the distance between two points in a weakly stationary random field with CVF $C(h)$ be:

$|h| = \sqrt{h_1^2 + h_2^2 + \dots + h_p^2}$.

If $C(h) = C(|h|)$, it depdens only on the distance and not its orientation, hence **isotropic**.

If the orientation is also a dependent then the field is **anisotropic**.

### 4.3.2. The Matern class of covariance functions

### 4.3.3. Spherical family of Covariance functions

### 4.3.4. Isotropic models allowing negative correlations

### 4.3.5. Basic models not second-order stationary

The **power** and the **linear** model are **not** second order stationary isotropic semivariograms.

#### The Power semivariogram

$\gamma(h) = \theta h^\alpha$, $\theta \geq 0, 0 \leq \alpha < 2$.

- If alpha greater than or equal to 2 the model is no longer intrinsic. 
- Unbounded (not have a sill)
- Pure nugget effect as alpha -> 0
- alpha = 1 is the linear semivariogram.

### 4.3.6. Nuggest effects and nested models

### 4.3.7. Accomodating anisotropy

We consider two forms of anisotropic: geometric and zonal.


## 4.4. Sample Semivariograms (Lecture 4C)

## 4.5. Parametric Modelling of Semivariograms

Methods: Manual fitting; Least Squares; Likelihood-based.

### Manual fitting

### 4.5.1. Least Squares

#### Ordinary Least Squares (OLS)

minimises the distance between the estimate (observation?) and the true?

#### Generalised Least Squares (GLS)

#### Weighted Least Squares (WLS)

### 4.5.2. Likelihood-based methods


## 4.6. Kriging (Lecture 4D)

Kriging is a multi-dimensional extension of time series forecasting.

### Types of Kriging:
- Simple Kriging: when the mean is a known function
- Ordinary Kriging: The mean is a constant and needs to be estimated
- Universal Kriging: Mean unknown and may be estimated as a linear combination of a set of basis functions.

### Objectives of Kriging:

Given the observations, we would like to estimate the value of the random field at one or more **unsampled** locations from the sample data.

### 4.6.1. Simple Kriging

By applying Yule Walker with

$W = Z(s), Y = Z(s_0)$

$\Gamma = Cov(W) = Cov(Z(s)) = [Cov(Z(s_i), Z(s_j))]_{i, j = 1}^{n}$

$\gamma = Cov(Y, W) = Cov(Z(s_0), Z(s)) = [Cov(Z(s_i), Z(s_0))]_{i=1}^{n}$

We get the **simple kriging variance** results on page 9.

Note: If the nugget effect consists of micro-scale variability only, then we would honor the observed data, since there is **not** a structured portion of the spatial variability; If the nuggest effect contains a (significant) measurement error component, then we do not want the predictor to interpolate the observed data.

### 4.6.2. Ordinary Kriging

Ordinary kriging estimates a value for a point, $s_0$, in the region $\mathbb{D}$ for which the semivariogram is **known** using observed data in the neighbourhood of the estimation location.

The method of **Lagrange Multipliers** is used to find a minimum of a multivariate function subject to a side condition.

Note that OK is an **exact interpolator**.