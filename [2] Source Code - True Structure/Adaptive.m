%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [14] and [13]
%}

function [hat_X,P] = Adaptive(Y,F,G,H,Q,R,Pi_0,hat_x_0,isEstimateQ)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
n = length(hat_x_0);
hat_X = [];

if isEstimateQ
    Q_hat = eye(size(G*Q*G'));
    X_last = hat_x_0;
    P_last = Pi_0;
    for i=1:N
        X = F*X;
        Z_ = H*X;

        P = F*P*F' + Q_hat;
        K = P*H'*(H*P*H' + R)^-1;
        X = X + K*(Y(i) - Z_);

        hat_X = [hat_X X];

        P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R*K';

        % Eq. (23) and (24)
        delta_x = X - F*X_last;
        alpha = 0.01;
        if 1/i > alpha
            Q_hat = ((i-1)/i)*Q_hat + (1/i)*(delta_x*delta_x' + P - F*P_last*F');
        else
            Q_hat = (1-alpha)*Q_hat + alpha*(delta_x*delta_x' + P - F*P_last*F');
        end

        % Record
        X_last = X;
        P_last = P;
    end
else
    R_hat = eye(size(R));
    for i=1:N
        X = F*X;
        Z_ = H*X;

        P = F*P*F' + G*Q*G';
        K = P*H'*(H*P*H' + R_hat)^-1;
        X = X + K*(Y(i) - Z_);

        hat_X = [hat_X X];

        P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R_hat*K';

        % Eq. (18) and (19)
        v = Y(i) - H*X;
        alpha = 0.01;
        if 1/i > alpha
            R_hat = ((i-1)/i)*R_hat + (1/i)*(v*v' + H*P*H');
        else
            R_hat = (1-alpha)*R_hat + alpha*(v*v' + H*P*H');
        end
    end
end