%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [24]
%}

function [ hat_X,P ] = UB( Y,F,G,H,Q,R,Pi_0,hat_x_0 )

    if rank(H) ~= length(hat_x_0)
        warning('UB :: Assumption Voilated. Solution Unreliable !');
    end
    
    x_hat = hat_x_0;
    P = Pi_0;
    
    n = length(hat_x_0);
    hat_X = [];

    len = length(Y(1,:));
    for k = 1:len
        % Eq. (3)
        x_hat_p = F*x_hat;
        % Eq. (4)
        gamma = Y(:,k) - H*x_hat_p;
        
        % Eq. (26)
        %{
            alpha = sdpvar;
            obj = alpha;
            Constrj = [alpha >= 1.0, alpha*H*F*P*F'*H' + H*G*Q*G'*H' + R >= gamma*gamma'];
            solvesdp(Constrj,obj); % or use the function "optimize" to replace "solvesdp"
            alpha = value(alpha);
        %}
        % NOTE:
        % It is really time-consuming to invoke the above SDP solver, ...
        % so we use an approximated method below
        alpha = 1.0;
        while true
            if min(eig(alpha*H*F*P*F'*H' + H*G*Q*G'*H' + R - gamma*gamma')) >= 0
                break;
            end
            alpha = alpha + 0.1;
        end

        % Eq. (15)
        P_p = alpha*F*P*F' + G*Q*G';
        
        % Eq. (16)
        V = H*P_p*H' + R;
        
        % Eq. (25)
        K = P_p*H'*V^-1;
        
        % Eq. (5)
        x_hat = x_hat_p + K*gamma;
        
        % Eq. (17)
        P = (eye(n)-K*H)*P_p*(eye(n)-K*H)' + K*R*K';
        
        % Save 
        hat_X = [hat_X x_hat];
    end
end

