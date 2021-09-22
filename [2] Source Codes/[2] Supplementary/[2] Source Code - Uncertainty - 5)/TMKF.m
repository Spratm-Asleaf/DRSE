%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE
%}

function [hat_X,P] = TMKF(Y,true_q,true_r,F,G,H,Q,R,Pi_0,hat_x_0)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
n = length(hat_x_0);
hat_X = [];

for i=1:N
    X = F*X;
    Z_ = H*X;

    P = F*P*F' + G*(Q + true_q(i))*G';

    K = P*H'*(H*P*H' + R + true_r(i))^-1;
    X = X + K*(Y(i) - Z_);
    
    hat_X = [hat_X X];

    P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*(R+true_r(i))*K';
end
