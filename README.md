# Latent-Factor Half-Trek Criterion

This repository contains the code used to compute the numerical experiments in Section 6 of the paper

Barber, R. F., Drton, M., Sturma, N., and Weihs L. (2022).
Half-Trek Criterion for Identifiability of Latent Variable Models.
[arXiv](https://arxiv.org/abs/2201.04457).

The code is written in `R` and `SageMath` and the repository is structured as follows:

* `R`: Code to generate all isomorphic classes of DAGs of experimental setups 1 and 2 and to check the latent-factor half-trek criterion on all of these graphs using the package `SEMID`.
* `SageMath`: Code to check identifiability by algeraic techniques as described in Appendix C of the paper.
* `experiments`: The results of the experiments.
