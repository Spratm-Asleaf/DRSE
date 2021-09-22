%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE
%}

function [xhat, V, G, S] = MKF(sys, Y_T, x_0, V_0, gamma_1, gamma_2)
    T = size(Y_T, 2);
    n = size(sys.F, 1);
    m = size(sys.H, 1);
    d = n + m;

    x_prev = x_0;
    V_prev = V_0;
    xhat = zeros(n, T);
    V = zeros(n, n, T);
    G = zeros(n, m, T);
    S = zeros(d, d, T);
    for t = 1 : T
        [mu_t, Sigma_t] = predict(x_prev, V_prev, sys);
        [xhat(:,t), V(:,:,t)] = update(mu_t, Sigma_t, Y_T(:,t), length(x_0), gamma_1, gamma_2);

        x_prev = xhat(:,t);
        V_prev = V(:,:,t);
    end
end

function [mu_t, Sigma_t] = predict(x_prev, V_prev, sys)
    A_aug = [sys.F; sys.H * sys.F];
    B_aug = [
        sys.G*(sys.Q)*sys.G'                sys.G*(sys.Q)*sys.G'*sys.H'
        sys.H*sys.G*(sys.Q)*sys.G'          sys.H*sys.G*(sys.Q)*sys.G'*sys.H' + sys.R
    ];
    mu_t = A_aug * x_prev;
    Sigma_t = A_aug * V_prev * A_aug' + B_aug * B_aug';
end

function [xhat_t, V_t] = update(mu_t, Sigma_t, y_t, x_dim, gamma_1, gamma_2)
    S_star = SolveNSDP(Sigma_t, x_dim, gamma_1, gamma_2);

    S_xx = S_star(1:x_dim, 1:x_dim);
    S_xy = S_star(1:x_dim, end);
    S_yx = S_xy';
    S_yy = S_star(x_dim+1:end,x_dim+1:end);
    
    xhat_t = mu_t(1:x_dim,1) + S_xy * S_yy^-1 * (y_t - mu_t(x_dim+1:end,1));
    V_t =  S_xx - S_xy * S_yy^-1 * S_yx;
end