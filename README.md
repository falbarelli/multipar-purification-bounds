## Purification-based bounds for multiparameter quantum channel estimation

[![arXiv](https://img.shields.io/badge/arXiv-2104.11264-b31b1b.svg)](https://arxiv.org/abs/2104.11264)

This repository provides MATLAB code to evaluate the bounds derived in the paper [F. Albarelli and R. Demkowicz-Dobrzański, "Probe incompatibility in multiparameter noisy quantum channel estimation"](http://arxiv.org/abs/2104.11264)

In the paper the general scenario of bounding the sum of the single-parameter QFI of different channels is considered (i.e. the random sensing scenario described in the paper).
However, the functions provided here work only in the multiparameter channel estimation setting, where there is only one channel.

### Usage

We provide the function
```
[trFIbound,hvec ] = totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,variant,dimIn,dimOut,wInv)
```
that evaluates the optimal total QFI `trFIbound`, either for a single use of the channel `variant='single'` or it gives the asymptotic standard quantum limit bound with `variant='asymptotic'`.
<!-- 
```
	[] = optistate_SDP()
```
gives an optimal state attaining the optimal total QFI and the corresponding QFI matrix.
 -->
Some possibly nontrivial conventions for using the function are that the Kraus operators `K1,K2,...` are stacked in column as `KrausOps=[ K1 ; K2 ; K3; ... ]` and their partial derivatives are also stacked in column, one parameter after the other, as `KrausOpsDeriv=[ d1K1 ; d1K2 ; d1K3; ... d2K1 ; d2K2; d2K3 ...   ]`.
The purification matrices are returned in a single tensor of dimension `[ nKraus, nKraus, npar ]` and are individually accessed as `hvec(:,:,i)`. 
The other input arguments should be intuitive.
The functions do not have sanity checks, so make sure to provide consistent input data.
A detailed description of the input and output arguments can be obtained by typing `help totalQFI_SDP` in the MATLAB Command Window.
<!-- and `help optistate_SDP` in the MATLAB Command Window. -->

### Examples

The files `hamiltonian_erasure.m`, `multiphase_loss.m`, `phase_loss.m`, `phase_dephasing.m`, `GADchannel.m` and `grover_dephasing.m` are the examples considered in the paper and can be inspected for further instruction on how to use the main function and to check numerically the analytical results reported in the paper.

### Requirements

* [MATLAB](https://www.mathworks.com/products/matlab.html). 
* [QETLAB](https://github.com/nathanieljohnston/QETLAB).
* [YALMIP](https://yalmip.github.io/) with an SDP solver: [SDPT3](https://github.com/SQLP/SDPT3) and [SCS](https://github.com/cvxgrp/scs) are free, while [MOSEK](http://cvxr.com/cvx/doc/mosek.html) is commercial, but free for academic use.
