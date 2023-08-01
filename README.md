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

## Required Input Parameters

_Table of Required Input Parameters:_
| Parameter | Statistical Notation | Variable Name | Description |
| ---                             | ---           | --- | --- |
| Number of clusters              | $K$              | `K_input`     | Number of clusters in each treatment arm |
| Cluster size                    | $m$              | `m_input`     | Number of individuals in each cluster |
| Statistical Power               | $\pi$            | `power_input` | Probability of detecting a true effect under $H_A$ |
| Family-wise false positive rate | $\alpha$         | `alpha_input` | Probability of one or more Type I error(s) |
| Effect for $Y_1$                | $\beta_1^*$      | `beta1_input` | Estimated intervention effect on the first outcome ($Y_1$) |
| Effect for $Y_2$                | $\beta_2^*$      | `beta2_input` | Estimated intervention effect on the second outcome ($Y_2$) |
| Effect for $Y_c$                | $\beta_c^*$      | `betaC_input` | Estimated intervention effect on the combined outcome ($Y_c$) |
| Endpoint-specific ICC for $Y_1$ | $\rho_0^{(1)}$   | `rho01_input` | Correlation for $Y_1$ for two different individuals in the same cluster |
| Endpoint-specific ICC for $Y_2$ | $\rho_0^{(2)}$   | `rho02_input` | Correlation for $Y_2$ for two different individuals in the same cluster |
| Endpoint-specific ICC for $Y_c$ | $\rho_0^{(c)}$   | `rho0C_input` | Correlation for $Y_c$ for two different individuals in the same cluster |
| Inter-subject between-endpoint ICC | $\rho_1^{(1,2)}$ | `rho1_input`  | Correlation between $Y_1$ and $Y_2$ for two different individuals in the same cluster |
| Intra-subject between-endpoint ICC | $\rho_2^{(1,2)}$ | `rho2_input`  | Correlation between $Y_1$ and $Y_2$ for the same individual |
| Total variance of $Y_1$ | $\sigma_1^2$     | `varY1_input` | Total variance of the first outcome, $Y_1$ |
| Total variance of $Y_2$ | $\sigma_2^2$     | `varY2_input` | Total variance of the second outcome, $Y_2$ |
| Total variance of $Y_c$ | $\sigma_c^2$     | `varYC_input` | Total variance of the combined outcome, $Y_c$ |

## Function Description

Each method has a set of functions for calculating the statistical power ($\pi$), required number of clusters per treatment group ($K$), or cluster size ($m$) given a set of input parameters. The names of all functions offered in this package are listed below, organized by study design method.

### 1. P-Value Adjustments for Multiple Testing

- `calc_pwr_pval_adj()` calculates power for this method
- `calc_K_pval_adj()` calculates number of clusters per treatment group for this method
- `calc_m_pval_adj()` calculates cluster size for this method

### 2. Combined Outcomes Approach

- `calc_pwr_comb_outcome()` calculates power for this method
- `calc_K_comb_outcome()` calculates number of clusters per treatment group for this method
- `calc_m_comb_outcome()` calculates cluster size for this method

### 3. Single 1-Degree of Freedom (DF) Combined Test for Two Outcomes

- `calc_pwr_single_1dftest()` calculates power for this method
- `calc_K_single_1dftest()` calculates number of clusters per treatment group for this method
- `calc_m_single_1dftest()` calculates cluster size for this method

### 4. Disjunctive 2-DF Test for Two Outcomes

- `calc_pwr_disj_2dftest()` calculates power for this method
- `calc_K_disj_2dftest()` calculates number of clusters per treatment group for this method
- `calc_m_disj_2dftest()` calculates cluster size for this method

### 5. Conjunctive Intersection-Union Test for Two outcomes

- `calc_pwr_conj_test()` calculates power for this method
- `calc_K_conj_test()` calculates number of clusters per treatment group for this method
- `calc_m_conj_test()` calculates cluster size for this method

## Usage 

## Contact

For questions or comments, please email Melody Owen at melody.owen@yale.edu, or submit an issue to this repository. 
