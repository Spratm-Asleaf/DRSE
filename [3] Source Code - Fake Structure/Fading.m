%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [11]
%}

function [hat_X,P] = Fading(Y,F,G,H,Q,R,Pi_0,hat_x_0)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
n = length(hat_x_0);
hat_X = [];

lamda = 1;
C0 = 1;
G1 = 0;
G2 = 0;
lastP = [];
for i=1:N
    % Eq. (6)
	X = F*X;
    
    % Eq. (7)
    Z_ = H*X;

    % Eq. (30), Eq. (31), Eq. (29)
    G1 = G1/lamda + (Y(:,i) - Z_)^2;
    G2 = G2/lamda + 1;
    
    if abs(G2) < 1e-6
        C0 = 1;
    else
        C0 = G1/G2;
    end
    
    % Eq. (44), Eq. (45)
    if i == 1
        MM = H*F*P*F'*H';
    else
        MM = H*F*lastP*F'*H';
    end
    
    NN = C0 - H*G*Q*G'*H' - R;
    
    % Eq. (52)
    lamda = max(1, trace(NN)/trace(MM));

    % Eq. (9)+(11)
    P = lamda*F*P*F' + G*Q*G';
    %P = F*P*F' + G*Q*G';
    
    % Eq. (8)
    K = P*H'*(H*P*H' + R)^-1;
    
    % Eq. (7)
    X = X + K*(Y(:,i) - Z_);
    hat_X = [hat_X X];

    % Eq. (10)
    P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R*K';

    lastP = P;
end
