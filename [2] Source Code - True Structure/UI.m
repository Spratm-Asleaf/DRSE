%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [31]
%}

function [ hat_X,P ] = UI( Y,F,g,G,H,Q,R,Pi_0,hat_x_0 )
    x_hat = hat_x_0;
    P = Pi_0;
    
    n = length(hat_x_0);
    m = length(Y(:,1));
    hat_X = [];

    len = length(Y(1,:));

    if rank(H*g) ~= rank(g)
        warning('Unknown Input Filter :: Assumption 1 violated :: first =. Solution Unreliable !');
    elseif rank(g) ~= length(g(1,:))
        warning('Unknown Input Filter :: Assumption 1 violated :: second =. Solution Unreliable !');
    end
    for k = 1:len
        % Eq. (12)
        P_p = F*P*F' + G*Q*G';
        R_tilde = H*P_p*H' + R;
        
        % Eq. (13)
        f = H*g;
        M = (f'*R_tilde^-1*f)^-1 * f' * R_tilde^-1;
        
        % Eq. (3)
        x_hat_p = F*x_hat;
        
        % Eq. (4)
        d = M*(Y(:,k) - H*x_hat_p);
        
        % Eq. (5)
        x_hat_star = x_hat_p + g*d;
        
        % Eq. (20)
        K = P_p*H'*R_tilde^-1;
        
        % Eq. (6)
        x_hat = x_hat_star + K*(Y(:,k) - H*x_hat_star);
        
        % Eq. (25)
        P_star = (eye(n)-g*M*H)*P_p*(eye(n)-g*M*H)' + g*M*R*M'*g';
        
        % Eq. (28), (29)
        S = -g*M*R;
        V = P_star*H' + S;
        R_tilde_star = (eye(m)-H*g*M)*R_tilde*(eye(m)-H*g*M)';
        P = K*R_tilde_star*K' - V*K' - K*V' + P_star;
        
        % Save 
        hat_X = [hat_X x_hat];
    end
end

