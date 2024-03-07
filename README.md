We compare the relative performance of the ME method developed in 
"A mixed model approach to estimate the survivor average causal effect in cluster-randomized trials"
to other statistical estimands and estimators for assessing patient-centered outcomes missing due to death that researchers 
can or have used in practice.
The various methods are the fixed effects only method, the inverse probability weighting approach, 
the worst ranked value approach, and the win/ratio approach. 
We conduct the comparison through simulation studies focusing on the most common parallel cluster design.

R functions

fecode: SACE estimation with fixed effects only.

mecode: SACE estimation with random intercepts in the outcome models.

simfun: Functions used in the simulations.

simdata: Data simulation.

simfit: analysis of the simulated data.

simsummary: Performance comparison in terms of the power, the proportion of coverage, and the computation time.
