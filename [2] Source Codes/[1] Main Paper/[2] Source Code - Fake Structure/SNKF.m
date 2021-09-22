% Adapted from: "https://asl.epfl.ch/software/"

%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Site: https://github.com/Spratm-Asleaf/DRSE
%}


function [hat_X,P] = SNKF(Y,y0,F,G,H,Q,R,M,E_f,P_0,x_0)
[p,N]=size(Y);
[q,n]=size(E_f);
[n,m]=size(G);

hat_X = zeros(n, N);

% subscript "f" is for "filtered"
% Table I: initial conditions
P_f = ((P_0)^-1+H'*(R)^-1*H)^-1;
hat_x_f = P_f * H' * (R)^-1 * y0;

% Eq. (45)
A = [H*F H*G];
E_g = zeros(q,m);

for i=1:N
    %% Step 1
    b = Y(:,i)-H*F*hat_x_f;             % Eq. (45)
    Qb = zeros(2*n,2*n);                % Eq. (48)
    Qb(1:n,1:n) = (P_f)^-1; 
    Qb(n+1:2*n,n+1:2*n) = (Q)^-1; 
    W = R^-1;                           % Eq. (49) 
    H_ = H*M;                           % Eq. (50)
    E_a = [E_f E_g];                    % Eq. (51)
    E_b = -E_f*hat_x_f;                 % Eq. (52)

    % Sovle Eq. (11)
    lambda = bdu_rls(A,b,Qb,W,H_,E_a,E_b); 

    %% Step 2
    % See also Subsection IV-D (pp. 1004) for E_g = 0
    hat_R = R - lambda^-1 * H*(M*M')*H';
    hat_P_f = ((P_f)^-1 + lambda*(E_f'*E_f))^-1;
    hat_F = F*(eye(n,n) - lambda*hat_P_f*(E_f'*E_f));
    % hat_Q = Q;
    % hat_G = G;

    %% Step 3
    P = F*hat_P_f*F' + G*Q*G';
    R_e = hat_R + H*P*H';
    P_f = P - P*H'*(R_e)^-1*H*P;
    predicted_hat_x = hat_F*hat_x_f;
    hat_x_f = predicted_hat_x + P_f*H'*(hat_R)^-1*(Y(:,i) - H*predicted_hat_x); 

    hat_X(:,i) = hat_x_f;
end

function [lambda_hat] = bdu_rls(A,b,Q,W,H,E_a,E_b)
[N,n] = size(A);  % dimensions;
[N,m] = size(H);

% Evaluation of lower bound on lambda
l_bound = norm(H'*W*H);

% Upper bound is arbitrarily set  (can be changed).
u_bound = 5*l_bound;     % Extensive simulations have shown 
                      % that the optimal lambda is usually
                      % very close to l_bound or can be taken
                      % as l_bound!

epsilon=1e-6;
options=[];            % Default options

lambda_hat =  fminbnd(@cost_rls,l_bound-epsilon,u_bound,options,A,b,Q,W,H,E_a,E_b);

%  min_cost = cost_rls(lambda_hat,A,b,Q,W,H,E_a,E_b);
% 
%  Q_hat  = Q + lambda_hat * E_a'*E_a;
%  X1_hat = lambda_hat*eye(m,m) - H'*W*H;   
%  W_hat  = W + W*H*pinv(X1_hat)*H'*W; 
% 
%  X2_hat = Q_hat + A'*W_hat*A;
%  x_hat  = inv(X2_hat)*(A'*W_hat*b + lambda_hat*E_a'*E_b);

 
function G = cost_rls(lambda,A,b,Q,W,H,E_a,E_b)
 [N,m] = size(H);      % dimensions

 % Evaluation of W(lambda)
 X1 = lambda*eye(m,m) - H'*W*H;   
 W_lambda = W + W*H*pinv(X1)*H'*W; 

 % Evaluation of x(lambda)
 X2 = Q + lambda*(E_a'*E_a) + A'*W_lambda*A;
 x_lambda = (X2)^-1*(A'*W_lambda*b + lambda*E_a'*E_b);

 % Evaluation of G(lambda)
 X3 = x_lambda'*Q*x_lambda;
 X4 = (A*x_lambda - b)'*W_lambda*(A*x_lambda - b);
 X5 = lambda*(E_a*x_lambda - E_b)'*(E_a*x_lambda - E_b);
 G = X3 + X4 + X5;