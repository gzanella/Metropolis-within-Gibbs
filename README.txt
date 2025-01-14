This folder contains the code to replicate Figure 1 in Scalability of Metropolis-within-Gibbs schemes for high- dimensional Bayesian models (Ascolani, Roberts and Zanella).

It includes the following files:

- "Main.R", which contains the code to replicate the simulations.

- "Functions_Gibbs.R", which contains the functions to implement an exact Gibbs.

- "Functions_MH.R", which contains the functions to implement a Metropolis-within-Gibbs, where the local parameters are updated according to independent Barker proposals.

- "Functions_MALA.R", which contains the functions to implement a Metropolis-Adjusted-Lanegvin (MALA) algorithm.

- "hierarchical_logit_stan.stan", which contains the definition of the hierarchical model suitable for rStan.


See Section 6.5.1 for additional details on the simulation.

NOTE: in the "Main.R" file the exact Gibbs is emulated by a Metropolis-within-Gibbs with multiple steps per iteration. However, the code for an exact Gibbs is also available and commented.