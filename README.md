# garch_VaR
Simulate and estimate volatility by GARCH with/without leverage, riskmetriks. VaR compute and test on VaR Violation.
These codes use the package `rugarch` for Volatitly models

**Code files**
1. `Financial_Econometrics.R`: clean, summarize, and plot the daily return data
2. `garch_simulation.R`: simulate and compare the performance of different GARCH models under different distribution 
3. `interval_graph.R`: compare the estimated intervals of coefficient from different set-up of simulation
4. `garch_estimate.R`: GARCH and Value-ar-Risk estimate on daily return data + diagnostics in-sample/out-sample



## 1. GARCH

Several models would be applied: GARCH and GARCH-with leverage (eGARCH, TGARCH, APARCH). 
By AIC and BIC, GARCH(2,1) and eGARCH(2,1) outperforms. 

There are six models: (GARCH, eGARCH) x ("ged", "norm", "std"). 
They would be assessed in both in-sample and out-sample set.

## 2. VaR

We compute the VaR in-sample and out-sample (given theorerical distibutrion). 
Then create the hit sequence $H_t$ and conduct 3 tests in the behavior of VaR violations.
