# Canonical forms: linear and state-dependent local projections

## Common elements

**Observables.** $y_t \in \mathbb{R}^n$, $n = 5$, ordered as

$$
y_t = \big(\text{GDP}_t,\ \text{EM}_t,\ \pi_t,\ \text{INT}_t,\ \text{EN}_t\big)'
$$

| Symbol | Variable | Units |
|---|---|---|
| $\text{GDP}_t$ | real GDP | $100 \log(\cdot)$, levels |
| $\text{EM}_t$ | employment rate | levels |
| $\pi_t$ | inflation | levels |
| $\text{INT}_t$ | policy interest rate | levels |
| $\text{EN}_t$ | energy demand | $100 \log(\cdot)$, levels |

**Sample.** $t = 1, \dots, T$, where $t = 1$ is the first date at which the 4-period moving average $m_t$ (defined below) exists. The first 3 raw observations are dropped from **all** variables.

**Settings.** Lag length $p = 4$. Horizons $h = 0, 1, \dots, H$ with $H = 14$. Quarterly data.

**Estimation.** Every coefficient block is computed by OLS, separately for each horizon $h$. The regressors are always $y_{t-1}, \dots, y_{t-p}$ (dated before $t$), and $t$ runs over $p+1, \dots, T-h$, which gives $T-h-p$ observations.

**Dimensions.** $\alpha(h)$ is $n \times 1$. Each $B(i,h)$, $B_1(i,h)$, $B_2(i,h)$ is $n \times n$. $IR_h$ is $n \times n$.

**Reading $IR_h$.** Element $(k,j)$ of $IR_h$ is the response of variable $k$ at $t+h$ to structural shock $j$ at $t$. Shock $j$ is normalised to a **unit impact on its own variable**, so the diagonal of $IR_0$ is all ones.

---

## 1. Linear model (canonical form)

$$
\big[\alpha(h),\, B(1,h),\, \dots,\, B(p,h)\big] = Y_h X_h' \big(X_h X_h'\big)^{-1}
$$

$$
Y_h = \big[\, y_{p+1+h},\, \dots,\, y_T \,\big], \qquad
X_h = \begin{bmatrix} 1 & \cdots & 1 \\ y_p & \cdots & y_{T-h-1} \\ y_{p-1} & \cdots & y_{T-h-2} \\ \vdots & & \vdots \\ y_1 & \cdots & y_{T-h-p} \end{bmatrix}
$$

$$
y_{t+h} = \alpha(h) + \sum_{i=1}^{p} B(i,h)\, y_{t-i} + u_{t+h}(h)
$$

$$
IR_h = B(1,h-1)\, B_0^{-1}, \qquad B(1,-1) = I
$$

$$
B_0^{-1} = P \, \operatorname{diag}(P)^{-1}, \qquad P = \operatorname{chol}(\Sigma), \qquad \Sigma = \frac{1}{T-p} \sum_{t=p+1}^{T} u_t(0)\, u_t(0)'
$$

**Elements not written explicitly:**

- $Y_h$ is $n \times (T-h-p)$. $X_h$ is $(1+np) \times (T-h-p)$; each $y$ entry is an $n \times 1$ block, so column $t$ of $X_h$ is $(1, y_{t-1}', \dots, y_{t-p}')'$.
- $u_{t+h}(h)$ is the OLS residual at horizon $h$. $u_t(0)$ is the horizon-0 residual, which is the reduced-form residual of a VAR($p$) with intercept.
- **Error distribution.** The VAR innovations are assumed $u_t(0) \overset{\text{iid}}{\sim} \mathcal{N}(0, \Sigma)$. For $h \ge 1$, $u_{t+h}(h)$ is a moving average of $u_{t}(0), \dots, u_{t+h}(0)$. It is therefore serially correlated, with no separate distributional assumption, and inference on $B(i,h)$ uses Newey–West (HAC) standard errors.
- $\operatorname{chol}(\Sigma)$ is the **lower-triangular** Cholesky factor, $\Sigma = PP'$. The recursive ordering is the order of $y_t$: GDP responds on impact to all shocks, and EN only to its own.
- $\operatorname{diag}(P)$ is the diagonal matrix holding the diagonal of $P$. Because of this normalisation, $B_0^{-1}$ does not depend on any scalar rescaling of $\Sigma$, such as the degrees-of-freedom correction.
- $IR_0 = B_0^{-1}$.

---

## 2. State-dependent (nonlinear) model (canonical form)

States are indexed by $s \in \{1, 2\}$: $s = 1$ is expansion and $s = 2$ is recession.

$$
\big[\alpha(h),\, B_1(1,h), \dots, B_1(p,h),\, B_2(1,h), \dots, B_2(p,h)\big] = Y_h X_h' \big(X_h X_h'\big)^{-1}
$$

$$
Y_h = \big[\, y_{p+1+h},\, \dots,\, y_T \,\big], \qquad
X_h = \begin{bmatrix} 1 \;\cdots\; 1 \\[2pt] X^{(1)}_h \\[2pt] X^{(2)}_h \end{bmatrix}
$$

$$
\underbrace{X^{(1)}_h = \begin{bmatrix}
(1-F_p)\, y_p & \cdots & (1-F_{T-h-1})\, y_{T-h-1} \\
\vdots & & \vdots \\
(1-F_p)\, y_1 & \cdots & (1-F_{T-h-1})\, y_{T-h-p}
\end{bmatrix}}_{\text{multiplies } B_1(1,h), \dots, B_1(p,h)}
\qquad
\underbrace{X^{(2)}_h = \begin{bmatrix}
F_p\, y_p & \cdots & F_{T-h-1}\, y_{T-h-1} \\
\vdots & & \vdots \\
F_p\, y_1 & \cdots & F_{T-h-1}\, y_{T-h-p}
\end{bmatrix}}_{\text{multiplies } B_2(1,h), \dots, B_2(p,h)}
$$

$$
y_{t+h} = \alpha(h) + \sum_{i=1}^{p} \Big[(1-F_{t-1})\, B_1(i,h) + F_{t-1}\, B_2(i,h)\Big] y_{t-i} + u_{t+h}(h)
$$

$$
IR^{(s)}_h = B_s(1,h-1)\, B_0^{-1}, \qquad B_s(1,-1) = I, \qquad s \in \{1, 2\}
$$

$$
B_0^{-1} = P \, \operatorname{diag}(P)^{-1}, \qquad P = \operatorname{chol}(\Sigma), \qquad \Sigma = \frac{1}{T-p} \sum_{t=p+1}^{T} u_t(0)\, u_t(0)' \quad \text{(from the linear model)}
$$

$$
F_t = \frac{e^{-\gamma z_t}}{1 + e^{-\gamma z_t}}, \qquad \gamma = 2
$$

$$
z_t = \frac{c_t - \bar c}{\operatorname{sd}(c)}, \qquad c = \operatorname{HP}_{\lambda = 1600}(m), \qquad m_t = \tfrac14 \textstyle\sum_{j=0}^{3} \text{GDP}_{t-j}
$$

**Elements not written explicitly:**

- $X_h$ is $(1+2np) \times (T-h-p)$, and column $t$ is
  $$
  \big(1,\ (1-F_{t-1})y_{t-1}', \dots, (1-F_{t-1})y_{t-p}',\ F_{t-1}y_{t-1}', \dots, F_{t-1}y_{t-p}'\big)'.
  $$
  There is a single intercept $\alpha(h)$, common to both states.
- **Error distribution.** No distributional assumption is made on $u_{t+h}(h)$. It is serially correlated for $h \ge 1$, and inference uses Newey–West (HAC) standard errors. The only normality assumption enters through $\Sigma$, which is inherited from the linear model.
- The transition weight enters with **one lag**: the regression at date $t$ uses $F_{t-1}$, and the same $F_{t-1}$ multiplies all $p$ lags.
- **$\Sigma$, $u_t(0)$ and $B_0^{-1}$ are taken from the linear model** (Section 1, same $p$). Both states share the same impact matrix: $IR^{(1)}_0 = IR^{(2)}_0 = B_0^{-1}$. State dependence enters only through $B_s(1,h-1)$ for $h \ge 1$.
- **HP filter.** Given the series $m = (m_1, \dots, m_T)$, the trend $\tau$ solves
  $$
  \tau = \arg\min_{\tau} \sum_{t=1}^{T} (m_t - \tau_t)^2 + \lambda \sum_{t=2}^{T-1} \big[(\tau_{t+1} - \tau_t) - (\tau_t - \tau_{t-1})\big]^2,
  $$
  with $c = m - \tau$. The filter is estimated once on the **full sample** $1, \dots, T$, so it is two-sided: $c_t$ depends on $m_1, \dots, m_T$.
- $\bar c$ and $\operatorname{sd}(c)$ are the full-sample mean and standard deviation of $c$.
- **Interpretation of $F$.** $F_t \in (0,1)$ is decreasing in $z_t$, with $F_t = 0.5$ at $z_t = 0$. A high $z_t$ (GDP above trend) gives $F_t \to 0$, which weights state 1 (expansion). A low $z_t$ gives $F_t \to 1$, which weights state 2 (recession).
