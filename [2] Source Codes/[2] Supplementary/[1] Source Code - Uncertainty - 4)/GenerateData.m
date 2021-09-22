%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE
%}

function [x,y,y0,true_d] = GenerateData(sys, x_0, T)
    [m,n] = size(sys.H);

    y0 = sys.H * x_0 + sys.D * chol(sys.R) * randn(m,1);
    x_prev = x_0;
    x = zeros(n,T);
    y = zeros(m,T);
    
    true_d = [];

    for t = 1 : T
        now_d = rand;
        true_d = [true_d now_d];
        x(:,t) = sys.F * x_prev + sys.Gamma * now_d + sys.G * chol(sys.Q) * randn(n,1);
        y(:,t) = sys.H * x(:,t) + sys.D * chol(sys.R) * randn(m,1);
        x_prev = x(:,t);
    end
end