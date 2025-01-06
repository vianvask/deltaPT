# deltaPT
C++ code that computes the distribution of the density contrast generated in strongly supercooled phase transitions.

Please cite [2402.04158](https://arxiv.org/abs/2402.04158) if you use this code. The main ingredients of the computation are described in the supplemental material in that paper.

Running takes two input parameters beta/H0 and gamma/beta that parametrize the bubble nucleation rate. The code computes the distribution $P_k(\delta)$ of the density contrast $\delta$ for 18 wavenumbers $k$ by generating $10^6$ realizations of the evolution of the false vacuum fraction $F_k$. For each realization, the nucleation times and distances of the first $j_c=50$ bubbles are generated. The contribution from the $j>j_c$ bubbles is integrated. 

Output:
  - tkR...: largest wavenumber ($k_{\rm max}$) that exits horizon during the thermal inflation: $\left( t_{k_{\rm max}}, k_{\rm max}, H(t_{k_{\rm max}}), R_*(t_{k_{\rm max}}) \right)$
  - klist...: list of wavenumbers for which the distribution of $\delta$ is computed and the times when they reenter horizon: $\left( k, t_k, H(t_k) \right)$
  - FkW...: the time evolution of the average false vacuum fraction, the contribution to the false vacuum fraction from $j>j_c$ bubbles, the expected number of bubbles and its time derivative for $k = 0.9 k_{\rm max}$: $\left( t,\bar{F},F_k^{j>j_c},\bar{N}_k,\partial_t \bar{N}_k \right)$
  - Fk...: the time evolution of $F_k$ for $k = 0.9 k_{\rm max}$ in the first 10 simulations: $\left( F_k(t_1), F_k(t_2), \dots \right)$
  - deltak...: the time evolution of $\delta$ for $k = 0.9 k_{\rm max}$ in the first 10 simulations: $\left( \delta(t_1), \delta(t_2), \dots \right)$
  - tj...: the nucleation times time of $j<j_c$ bubbles for $k = 0.9 k_{\rm max}$: $\left( t_1, t_j, \dots, t_{j_c} \right)$
  - deltabinsk...: distribution of $\delta$: $\left(\delta, P_{k_1}(\delta), P_{k_2}(\delta), \dots \right)$
