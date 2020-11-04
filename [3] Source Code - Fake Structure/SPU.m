%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper: [36]
%}

function [ hat_X,P ] = SPU( Measure,F,G,H,D,Q,R,As,Bs,Cs,Ds,Pi_0,hat_x_0)
    x_f = hat_x_0;  % Used in Eq. (7a)
    P = Pi_0;
    
    n = length(hat_x_0);
    p = 2^n;
    
    %{{ ----- X_i(0), i = 1,2,...,p
        X = cell(1,p);
        for i = 1:p
            X{i} = diag(eig(Pi_0));
        end
        radius = 0.5;
        X{1}(1,1) = X{1}(1,1) - radius;
        X{1}(2,2) = X{1}(2,2) - radius;

        X{2}(1,1) = X{2}(1,1) - radius;
        X{2}(2,2) = X{2}(2,2) + radius;

        X{3}(1,1) = X{3}(1,1) + radius;
        X{3}(2,2) = X{3}(2,2) - radius;

        X{4}(1,1) = X{4}(1,1) + radius;
        X{4}(2,2) = X{4}(2,2) + radius;
    %}}
    
    hat_X = [];

    [a,b] = size(H);
    A = F;
    B = [G*Q*G' zeros(b,a)];
    C = H;
    D = [zeros(a,b) R*D*R'];
    
    %% Save Data
    B_hat_last = cell(1,p);
    D_hat_last = cell(1,p);
    X_last = cell(1,p);


    %% Step 1: Initialization. k = 1.
    % when k = 1 using Eq. (22)
    for i = 1:p
        X_last{i} = X{i};       % save X_i(0)
        B_hat = [As*chol(X{i}) B Bs];   B_hat_last{i} = B_hat;  % B_hat_i(0)
        D_hat = [Cs*chol(X{i}) D Ds];   D_hat_last{i} = D_hat;  % D_hat_i(0)
        X{i} = A*X{i}*A' + B*B' + As*X{i}*As' + Bs*Bs';         % X_i(1)
    end
    
    % all indexed by k. Now, Q := Q_k = Q_1
    Q = sdpvar(2,2);    % Q = P^-1, where P is 2 by 2
    M = sdpvar(2,2);    % upper bound of MSE. See M in Eq. (14)
    Y = sdpvar(2,1);    % Y = P^-1*K, where P is 2 by 2 and K is 2 by 1. See K in Eq. (7b)
    F = sdpvar(2,1);    % F is 2 by 1. See F in Eq. (7b)
    
    obj = trace(M);     % Eq. (22)
    
    T12_last = Q*A + Y*C;
    Constrj = [];
    for i = 1:p
        T13_last = Q*B_hat_last{i} + Y*D_hat_last{i};
        
        B_hat = [As*chol(X{i}) B Bs];   B_hat_last{i} = B_hat;
        D_hat = [Cs*chol(X{i}) D Ds];   D_hat_last{i} = D_hat;
        
        [~,col] = size(T13_last);
        LeftMat1 = [
            Q,                          T12_last*X_last{i},                    	T13_last; 
            X_last{i}*T12_last',        X_last{i},                             	zeros(length(X{i}(:,1)), col); 
            T13_last',                  zeros(col,length(X{i}(:,1))),           eye(col)
        ];
    
        Constrj = [Constrj, LeftMat1>=0];    % Eq. (22)
        
        F_D_hat = F*D_hat;
        [~,col] = size(F_D_hat);
        LeftMat2 = [
            M,                                  eye(length(F(:,1)))+F*C,            F_D_hat; 
            eye(length(F(:,1)))+C'*F',          Q,                                 	zeros(length(Q(:,1)), col); 
            F_D_hat',                           zeros(col,length(Q(:,1))),        	eye(col)
        ];
    
        Constrj = [Constrj, LeftMat2>=0];    % Eq. (22)
    end

    diagnostics = optimize(Constrj,obj);
    if diagnostics.problem ~= 0
        if diagnostics.problem == 1
            error('SPU :: Infeasible !')
        else
            error('SPU :: Unkown Error Happened !')
        end
    end
    
    
    Q_last = value(Q);
    
    % Eq. (7)
    K_last = pinv(value(Q))*value(Y);
    x_f = A*x_f;
    x_hat_p = x_f;
    x_hat = x_hat_p + value(F)*(C*x_hat_p - Measure(:,1));  % k = 1
    hat_X = [hat_X x_hat];
    
    %% Step 2: Eq. (23)
    % implemented in Step 1 and Step 3, respectively

    %% Step 3: Filtering
    len = length(Measure(1,:));
    for k = 2:len
        %% Solve Eq. (16)
        % when k >= 2 using Eq. (16)
        for i = 1:p
            X_last{i} = X{i};       % save X_i(k-1)
            X{i} = A*X{i}*A' + B*B' + As*X{i}*As' + Bs*Bs';     % X_i(k)
        end

        % all indexed by k. Now, Q := Q_k
        Q = sdpvar(2,2);    % Q = P^-1, where P is 2 by 2
        M = sdpvar(2,2);    % upper bound of MSE. See M in Eq. (14)
        Y = sdpvar(2,1);    % Y = P^-1*K, where P is 2 by 2 and K is 2 by 1. See K in Eq. (7b)
        F = sdpvar(2,1);    % F is 2 by 1. See F in Eq. (7b)

        obj = trace(M);     % Eq. (16a)

        T12_last = Q*A + Y*C;     % Eq. (16d)
        Constrj = [];
        for i = 1:p
            T13_last = Q*B_hat_last{i} + Y*D_hat_last{i};     % Eq. (16d)

            B_hat = [As*chol(X{i}) B Bs];   B_hat_last{i} = B_hat;     % Eq. (16d)
            D_hat = [Cs*chol(X{i}) D Ds];   D_hat_last{i} = D_hat;     % Eq. (16d)

            [~,col] = size(T13_last);
            LeftMat1 = [
                Q,                  T12_last,                    	T13_last; 
                T12_last',          Q_last,                       	zeros(length(Q_last(:,1)), col); 
                T13_last',       	zeros(col,length(Q_last(:,1))),	eye(col)
            ];

            Constrj = [Constrj, LeftMat1>=0];    % Eq. (16b)

            F_D_hat = F*D_hat;
            [~,col] = size(F_D_hat);
            LeftMat2 = [
                M,                                  eye(length(F(:,1)))+F*C,            F_D_hat; 
                eye(length(F(:,1)))+C'*F',          Q,                                 	zeros(length(Q(:,1)), col); 
                F_D_hat',                           zeros(col,length(Q(:,1))),        	eye(col)
            ];

            Constrj = [Constrj, LeftMat2>=0];    % Eq. (16c)
        end

        diagnostics = optimize(Constrj,obj);
        if diagnostics.problem ~= 0
            if diagnostics.problem == 1
                error('SPU :: Infeasible !')
            else
                error('SPU :: Unkown Error Happened !')
            end
        end
        
        Q_last = value(Q);
        
        % Eq. (7)
        K_last = pinv(value(Q))*value(Y);
        x_f = A*x_f + K_last*(C*x_f - Measure(:,k-1));
        x_hat_p = x_f;
        x_hat = x_hat_p + value(F)*(C*x_hat_p - Measure(:,k));  % k = 1
        hat_X = [hat_X x_hat];
    end
end

