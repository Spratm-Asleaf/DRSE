%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim
From the Department of Industrial Systems Engineering and Management, 
National University of Singapore (S. Wang and A. Lim); and 
the School of Management Science and Engineering, 
Nanjing University of Information Science and Technology (Z. Wu).

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Acknowledgements:
1. "SNKF.m" is adapted from sources at 
    "https://asl.epfl.ch/software/".
2. "WKF.m", "FrankWolfe.m" and "tauKF.m" are adapted from sources at
    "https://github.com/sorooshafiee/WKF".
3. For "tauKF.m", we also acknowledge
    "https://mathworks.com/matlabcentral/fileexchange/54308-robust-kalman-filtering-package?focused=5751770&tab=function"
%}

clc;
clear all;
close all;

% Problem setting
% n = 2; % dimension of state vector
% m = 1; % dimension of measurement
delta_t = 1;
sys.F = [1, delta_t; 0, 1];
%sys.F = [sys.F,sys.F*0;sys.F*0,sys.F];
sys.G = [delta_t^2/2; delta_t];
%sys.G = [sys.G, sys.G*0;sys.G*0, sys.G];
sys.H = [1, 0];
%sys.H = [1, 0, 0, 0; 0, 0, 1, 0];
sys.D = 1;
sys.Q = 0.1;
sys.R = 400;

% Simulate 1000 steps
T = 1000;

x_0 = randn(2,1);
V_0 = eye(2);

% Parameters
[x, y, y0,true_q,true_r] = GenerateData(sys, x_0, T);

%% Optimal Kalman filter (i.e., KF with the true model)
% suffix "TM" is for "True Model"
x_0 = [0;0];
tic
xhat_TMKF = TMKF(y,true_q,true_r,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
err_TMKF = sum((x - xhat_TMKF).^2,1);


%% Standard Kalman filter
% The immediate line below is to test the validity of WKF.
% Note that when the radius of uncertainty set equals to 0, ...
% WKF degrades into the standard KF.
% xhat = WKF(sys, 0, y, x_0, V_0); err_KF_0 = sum((x-xhat).^2,1);
x_0 = [0;0];
tic
xhat_KF = KF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
err_KF = sum((x - xhat_KF).^2,1);

%% Moment filter
gamma = 0.005;   % 0.025
gamma_1 = 1 - gamma;
gamma_2 = 1 + gamma;
xhat_MKF_1 = MKF(sys, y, x_0, V_0, gamma_1, gamma_2);
err_MKF_1 = sum((x - xhat_MKF_1).^2,1);

gamma = 0.02;   % 0.025
gamma_1 = 1 - gamma;
gamma_2 = 1 + gamma;
xhat_MKF_2 = MKF(sys, y, x_0, V_0, gamma_1, gamma_2);
err_MKF_2 = sum((x - xhat_MKF_2).^2,1);

gamma = 0.05;   % 0.025
gamma_1 = 1 - gamma;
gamma_2 = 1 + gamma;
xhat_MKF_3 = MKF(sys, y, x_0, V_0, gamma_1, gamma_2);
err_MKF_3 = sum((x - xhat_MKF_3).^2,1);


%% Plot
figure;
smt = 50; %50
% can also try "semilogx" instead of "plot"
plot(smooth(10*log10(err_TMKF), smt), 'k', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_KF), smt), 'r', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_1), smt), 'm', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_2), smt), 'b', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_3), smt), 'g', 'LineWidth', 2);
sqrt(mean(err_TMKF))
sqrt(mean(err_KF))
sqrt(mean(err_MKF_1))
sqrt(mean(err_MKF_2))
sqrt(mean(err_MKF_3))
% axis([0 1000 -5 50]);

% Plot Format
set(gca, 'FontSize', 18);
ylabel('Estimation Error (dB)','FontSize', 20, 'Interpreter', 'latex');
xlabel('Time Step','FontSize', 20, 'Interpreter', 'latex');
leg1 = legend({'TMKF', 'KF', 'MKF (1.005)', 'MKF (1.02)', 'MKF (1.05)'}, 'Interpreter', 'latex', 'Location', 'northwest');

set(leg1,'FontSize',15);



