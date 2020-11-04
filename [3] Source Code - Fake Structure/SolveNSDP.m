%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper:
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
%}

function S_star = SolveNSDP(Sigma_k, x_dim, gamma_1, gamma_2)
    S_star = gamma_2 * Sigma_k;
end

