%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE
%}

function [x,y,y0,true_F] = GenerateData(sys, x_0, T, coeff, is_TV)
    [m,n] = size(sys.H);

    y0 = sys.H * x_0 + sys.D * chol(sys.R) * randn(m,1);
    x_prev = x_0;
    x = zeros(n,T);
    y = zeros(m,T);
    true_F = cell(1,T);
    Delta = 2 * rand - 1;   % for changable
%     Delta = 1;              % for fixed 
    A_purt = sys.F + [0, coeff * Delta; 0, 0];
    for t = 1 : T
        true_F{t} = A_purt;
        %A_purt
        x(:,t) = A_purt * x_prev + sys.G * chol(sys.Q) * randn(n,1);
        y(:,t) = sys.H * x(:,t) + sys.D * chol(sys.R) * randn(m,1);
        x_prev = x(:,t);
        if is_TV
            Delta = 2 * rand - 1;
            A_purt = sys.F + [0, coeff * Delta; 0, 0];
        end
    end
end