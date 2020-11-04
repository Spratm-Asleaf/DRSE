%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [20]
%}

function [hat_X,P] = HINFTY(Y,F,G,H,Q,R,Pi_0,hat_x_0)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
[m,n] = size(H);

hat_X = [];

L = eye(n);

% in this case, we cannot know a gamma in advance such that it guarantees Eq. (8)
gamma = 25;     
% show the gamma
disp(['{{ H-infinity :: gamma: ' num2str(gamma) '}}']);

for i=1:N
    X = F*X;
    Z_ = H*X;

    % Eq. (10)
    R_e = [eye(m), zeros(m,n); zeros(n,m), -(gamma^2)*eye(n)] + [H; L]*P*[H' L'];
    % Eq. (9)
    P = F*P*F' + G*G' - F*P*[H' L']*R_e^-1*[H; L]*P*F';  % note that P' = P

    % Eq. (12)
    K = P*H'*(H*P*H' + eye(m))^-1;
    
    % Eq. (11)
    X = X + K*(Y(i) - Z_);
    
    hat_X = [hat_X X];
end
