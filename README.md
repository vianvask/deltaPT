# deltaPT
C++ code that computes the distribution of the density contrast generated in strongly supercooled phase transitions.

Please cite [2402.04158](https://arxiv.org/abs/2402.04158) if you use this code. The main ingredients of the computation are described in the supplemental material in that paper.

Running takes two input parameters beta/H0 and gamma/beta that parametrize the bubble nucleation rate. 

Output:
  - tkR...: largest k value ($k_{\rm max}$) that exits horizon during the thermal inflation: $\left( t_{k_{\rm max}}, k_{\rm max}, H(t_{k_{\rm max}}), R_*(t_{k_{\rm max}}) \right)$
  - klist...: list of k values and the times when they reenter horizon: $\left( k, t_k, H(t_k) \right)$
  - FkW...: the time evolution of $\bar{F}$, $\bar{N}_k$, $\partial_t \bar{N}_k$ and $F_k^{j>j_c}$ for $k = 0.9 k_{\rm max}$: $\left( t,\bar{F},F_k^{j>j_c},\bar{N}_k,\partial_t \bar{N}_k \right)$
