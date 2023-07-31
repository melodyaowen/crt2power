# hybrid2power

## Overview

`hybrid2power` is an R package that allows users to calculate the statistical power or sample size of their hybrid type 2 cluster randomized trials (CRTs), given a set of input parameters. Hybrid 2 studies are studies where there are two co-primary outcomes, namely an implementation outcome (such as fidelity or reach) and a health outcome (such as infection rates, or change from baseline health scores). When powering these studies, cluster correlations and the inflation of the Type I error rate must be accounted for.

The five key study design approaches are included in this package that can be used to power hybrid 2 CRTs. 
1. P-Value Adjustments for Multiple Testing
2. Combined Outcomes Approach
3. Single 1-Degree of Freedom (DF) Combined Test for Two Outcomes
4. Disjunctive 2-DF Test for Two Outcomes
5. Conjunctive Intersection-Union Test for Two outcomes

 For details on the methods listed above, please refer to the publication that discusses these methods, available here. (Add link)

## Installation

```
install.packages("devtools")
require(devtools)
install_github("hybrid2power")
```

## Study Design Methods

_Table of Input Parameter Notation:_
| Parameter | Variable Name | Description | Required for... | 
| ---               | ---           | --- | --- |
| $K$               | `K_input`     | X | X |
| $m$               | `m_input`     | X | X |
| $\pi$             | `power_input` | X | X |
| $\alpha$          | `alpha_input` | X | X |
| $\beta_1^*$       | `beta1_input` | X | X |
| $\beta_2^*$       | `beta2_input` | X | X |
| $\beta_c^*$       | `betaC_input` | X | X |
| $\rho_0^{(1)}$    | `rho01_input` | X | X |
| $\rho_0^{(2)}$    | `rho02_input` | X | X |
| $\rho_0^{(c)}$    | `rho0C_input` | X | X |
| $\rho_1^{(1,2)}$  | `rho1_input`  | X | X |
| $\rho_2^{(1,2)}$  | `rho2_input`  | X | X |
| $\sigma_1^2$      | `var1_input`  | X | X |
| $\sigma_2^2$      | `var2_input`  | X | X |
| $\sigma_c^2$      | `varC_input`  | X | X |

### 1. P-Value Adjustments for Multiple Testing

**Hypothesis Test Framework**:  &nbsp;  &nbsp; $H_0$: $\beta_1^* = 0$ and $\beta_2^* = 0$  &nbsp;  &nbsp; vs.  &nbsp;  &nbsp;  $H_A$: $\beta_1^* \neq 0$ or $\beta_2^* \neq 0$

### 2. Combined Outcomes Approach

**Hypothesis Test Framework**:  &nbsp;  &nbsp; $H_0$: $\beta_c^* = 0$  &nbsp;  &nbsp; vs.  &nbsp;  &nbsp;  $H_A$: $\beta_c^* \neq 0$

### 3. Single 1-Degree of Freedom (DF) Combined Test for Two Outcomes

**Hypothesis Test Framework**:  &nbsp;  &nbsp; $H_0$: $\beta_1^* = \beta_2^* = 0$  &nbsp;  &nbsp; vs.  &nbsp;  &nbsp;  $H_A$: $\beta_1^* \neq 0$ or $\beta_2^* \neq 0$

### 4. Disjunctive 2-DF Test for Two Outcomes

**Hypothesis Test Framework**:  &nbsp;  &nbsp; $H_0$: $\boldsymbol{L} \boldsymbol{\beta}^* = 0$ &nbsp;  &nbsp; vs.  &nbsp;  &nbsp;  $\boldsymbol{L} \boldsymbol{\beta}^* \neq 0$

### 5. Conjunctive Intersection-Union Test for Two outcomes

**Hypothesis Test Framework**:  &nbsp;  &nbsp; $H_0$: $\beta_1^* = 0$ and $\beta_2^* = 0$  &nbsp;  &nbsp; vs.  &nbsp;  &nbsp;  $H_A$: $\beta_1^* \neq 0$ or $\beta_2^* \neq 0$


## Usage 

## Contact
