# BayTensor
R code for simulation studies in the paper:

Wang, Kunbo, and Yanxun Xu. "Bayesian Tensor-on-Tensor Regression with Efficient Computation." arXiv preprint arXiv:2210.11363 (2022).

The paper is available at https://arxiv.org/abs/2210.11363.

# Code
## File: Main.R
Simulation analyses under the proposed BayTensor MCMC algorithm (Algorithm1) and BayTensor Fast algorithm (Algorithm2) 

The Main.R script sources the following functions/scripts: Functions.R, Generate_Data.R, BayTensor_MCMC.R, BayTensor_Fast.R with the details explained below.

Full simulation setups can be found in Section 6 of the [paper](https://arxiv.org/abs/2210.11363) 

## File: Functions.R 
1. array_to_matrix: 
	A function that matricize a tensor to a matrix according to Equation (2) in the [paper](https://arxiv.org/abs/2210.11363)

2. tensor_tensor:
	Tensor contracted product defined in Section 2.1

3. tensor_matrix:
	Tensor matrix product, a special case of contracted product where the second input is a matrix

4. matrix_outer:
	Kronecker product of two matrices, defined in Section 2.1

5. Find_candidate:
	A function that returns core tensor dimensions in the neighbour of current core tensor dimension

## File: Generate_Data.R
1. Generate_data: Generate simulation data according to Section 6.2 in the [paper](https://arxiv.org/abs/2210.11363)

2. Generate_data_cor: Generate correlated simulation data according to Section 6.3 in the [paper](https://arxiv.org/abs/2210.11363)

## File: BayTensor_MCMC.R
1. Estimate: 
	A function that generate posterior samples of core tensor and factor matrices given the dimension of core tensor

2. BayTensor_MCMC:
	Full BayTensor MCMC algorithm (Algorithm 1) in the [paper](https://arxiv.org/abs/2210.11363)

## File: BayTensor_Fast.R
1. MAP:
	A function that returns MAP estimators of parameters given the dimension of core tensor

2. BayTensor_Fast:
	Full BayTensor Fast algorithm (Algorithm 2) in the [paper](https://arxiv.org/abs/2210.11363)


