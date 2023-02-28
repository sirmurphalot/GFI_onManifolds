Generalized Fiducial Inference on Differentiable Manifolds
========

**author**: `Alexander Murph`

For questions, issues or clarifications please reach out to Murph:
acmurph@unc.edu.

Overview
----

Code implementations of the models discussed in *Generalized Fiducial Inference on Differentiable Manifolds* by A. C. Murph, J. Hannig, and J. Williams.  Manuscript on [arXiv](https://arxiv.org/abs/2209.15473).

This repository is organized according to the three difference experiments that use three different MCMC algorithms.  The first algorithm considered (HMC) is [Brubaker et al's Constrained Hamiltonian Monte Carlo](http://www.cs.toronto.edu/~mbrubake/projects/cmcmc/) method and the second algorithm considered (MH) is [Zappa et al's Constrained Metropolis-Hastings](https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21783) algorithm.  The three problems considered are: inference for data from a multivariate normal density with the mean parameters on a sphere, a linear logspline density estimation problem, and a reimagined approach to the AR(1) mode

Reference
----
<a id="1">[1]</a> 
Murph, A. C., Hannig, J., Williams, J. P.  Generalized Fiducial Inference on Differentiable Manifolds

Citations
----
<a id="1">[2]</a> 
Brubaker, M. A., Salzmann, M. and Urtasun, R. (2012) A Family of MCMC Methods on Implicitly Defined
Manifolds. In Proceedings of the Fifteenth International Conference on Artificial Intelligence and Statistics, vol. 22, 161–172. <br>
<a id="1">[3]</a> 
Zappa, E., Holmes-Cerfon, M. and Goodman, J. (2018) Monte Carlo on Manifolds: Sampling Densities and Integrating Functions. Communications on Pure and Applied Mathematics, 71, 2609–2647.


`R` Packages Required
----
`ggplot2`, `reshape`
