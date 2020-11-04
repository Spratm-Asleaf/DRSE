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
    "https://asl.epfl.ch/software/". See [4];
2. "WKF.m", "FrankWolfe.m" and "tauKF.m" are adapted from sources at
    "https://github.com/sorooshafiee/WKF". See [57].
3. For "tauKF.m", we also acknowledge M. Zorz [23].
    "https://mathworks.com/matlabcentral/fileexchange/54308-robust-kalman-filtering-package?focused=5751770&tab=function"
%}


clc;
clear all;
close all;

if ~exist('YALMIP','dir')
	msgbox('Make sure you have unzipped the YALMIP.zip in this folder !','Information','Warn');
    return;
end

% Problem setting
% n = 2; % dimension of state vector
% m = 1; % dimension of measurement
sys.F = [0.9802, 0.0196; 0, 0.9802];
sys.G = [1, 0; 0, 1];
sys.H = [1, -1];
sys.D = 1;
sys.Q = [1.9608, 0.0195; 0.0195, 1.9605];
sys.R = 1;

% Simulate 1000 steps
T = 1000;

x_0 = randn(2,1);
V_0 = eye(2);

% Parameters
isTimeVariant = true;
coeff_alpha = 5;
[x, y, y0, real_F] = GenerateData(sys, x_0, T, coeff_alpha, isTimeVariant);

%% Optimal Kalman filter (i.e., KF with the true model)
% suffix "TM" is for "True Model"
x_0 = [0;0];
tic
xhat_TMKF = TMKF(y,real_F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
time_TMKF = toc/T;
disp(['Avg Time of TMKF: ' num2str(time_TMKF)]);
err_TMKF = sum((x - xhat_TMKF).^2,1);
RMSE_TMKF = sqrt(mean(err_TMKF));
disp(['-------RMSE of TMKF: ' num2str(RMSE_TMKF)]);


%% Standard Kalman filter
% The immediate line below is to test the validity of WKF.
% Note that when the radius of uncertainty set equals to 0, ...
% WKF degrades into the standard KF.
% xhat = WKF(sys, 0, y, x_0, V_0); err_KF_0 = sum((x-xhat).^2,1);
x_0 = [0;0];
tic
xhat_KF = KF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
time_KF = toc/T;
disp(['Avg Time of KF: ' num2str(time_KF)]);
err_KF = sum((x - xhat_KF).^2,1);
RSME_KF = sqrt(mean(err_KF));
disp(['-------RMSE of KF: ' num2str(RSME_KF)]);


%% Adaptive Kalman filter
x_0 = [0;0];
isEstimateQ = true;     % we can choose to estimate Q (or R) given R (or Q)
tic
xhat_Adapt = Adaptive(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,isEstimateQ);
time_Adapt = toc/T;
disp(['Avg Time of Adaptive: ' num2str(time_Adapt)]);
err_Adapt = sum((x - xhat_Adapt).^2,1);
RMSE_Adapt = sqrt(mean(err_Adapt));
disp(['-------RMSE of Adaptive: ' num2str(RMSE_Adapt)]);


%% Fading Kalman filter
x_0 = [0;0];
tic
xhat_Fading = Fading(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
time_Fading = toc/T;
disp(['Avg Time of Fading: ' num2str(time_Fading)]);
err_Fading = sum((x - xhat_Fading).^2,1);
RMSE_Fading = sqrt(mean(err_Fading));
disp(['-------RMSE of Fading: ' num2str(RMSE_Fading)]);


%% H-infinity filter
x_0 = [0;0];
tic
xhat_HINFTY = HINFTY(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
time_HINFTY = toc/T;
disp(['Avg Time of HINFTY: ' num2str(time_HINFTY)]);
err_HINFTY = sum((x - xhat_HINFTY).^2,1);
RMSE_HINFTY = sqrt(mean(err_HINFTY));
disp(['-------RMSE of HINFTY: ' num2str(RMSE_HINFTY)]);


%% UB filter
addpath(genpath([pwd '\YALMIP']));      % Solve SDP
x_0 = [0;0];
tic
xhat_UB = UB(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0 );
time_UB = toc/T;
disp(['Avg Time of UB: ' num2str(time_UB)]);
err_UB = sum((x - xhat_UB).^2,1);
RMSE_UB = sqrt(mean(err_UB));
disp(['-------RMSE of UB: ' num2str(RMSE_UB)]);


%% UI filter
x_0 = [0;0];
g = [0; 1];          % G in Eq. (1)
tic
xhat_UI = UI(y,sys.F,g,sys.G,sys.H,sys.Q,sys.R,V_0,x_0 );
time_UI = toc/T;
disp(['Avg Time of UI: ' num2str(time_UI)]);
err_UI = sum((x - xhat_UI).^2,1);
RMSE_UI = sqrt(mean(err_UI));
disp(['-------RMSE of UI: ' num2str(RMSE_UI)]);


%% SPU filter
try
    addpath(genpath([pwd '\YALMIP']));      % Solve SDP
    As = [0 0; 3 0];
    Bs = zeros(size(sys.G));
    Cs = zeros(size(sys.H));
    Ds = zeros(size(sys.H));
    x_0 = [0;0];
    tic
    % use evalc to discard the temporary output of the SU filter 
    SPU_log = evalc('xhat_SPU = SPU(y,sys.F,sys.G,sys.H,sys.D,sys.Q,sys.R,As,Bs,Cs,Ds,V_0,x_0)');
    time_SPU = toc/T;
    disp(['Avg Time of SPU: ' num2str(time_SPU)]);
    err_SPU = sum((x - xhat_SPU).^2,1);
    RMSE_SPU = sqrt(mean(err_SPU));
    disp(['-------RMSE of SPU: ' num2str(RMSE_SPU)]);
catch
    warning('main :: SPU no longer works !');
end

%% Sayed's Norm-constrained Kalman filter
M = [0.0198; 0];
E_f = [5 0];       % [0 10]   [0 50]
x_0 = [0;0];
tic
xhat_SNKF = SNKF(y,y0,sys.F,sys.G,sys.H,sys.Q,sys.R,M,E_f,V_0,x_0);
time_SNKF = toc/T;
disp(['Avg Time of SNKF: ' num2str(time_SNKF)]);
err_SNKF = sum((x - xhat_SNKF).^2,1);
RMSE_SNKF = sqrt(mean(err_SNKF));
disp(['-------RMSE of SNKF: ' num2str(RMSE_SNKF)]);


%% tau filter (Note: when tau-divergence filter takes tau = 0, it gives the Kullback-Leibler divergence filter)
tau = 0;
c = 1.5e-4;
y_delay = [y0, y(:,1:end-1)]';
x_0 = [0;0];
tic
xhat_tauKF = tauKF(sys, c, tau, y_delay, x_0, V_0);
time_tauKF = toc/T;
disp(['Avg Time of tauKF: ' num2str(time_tauKF)]);
err_tauKF = sum((x - xhat_tauKF').^2,1);
RMSE_tauKF = sqrt(mean(err_tauKF));
disp(['-------RMSE of tauKF: ' num2str(RMSE_tauKF)]);


%% Wasserstein filter
x_0 = [0;0];
tic
rho = 0.1;
[xhat_WKF, V_WKF, G_WKF, S_WKF] = WKF(sys, rho, y, x_0, V_0);
time_WKF = toc/T;
disp(['Avg Time of WKF: ' num2str(time_WKF)]);
err_WKF = sum((x - xhat_WKF).^2,1);
RMSE_WKF = sqrt(mean(err_WKF));
disp(['-------RMSE of WKF: ' num2str(RMSE_WKF)]);


%% Moment filter
tic
gamma = 0.02;   % 0.025
gamma_1 = 1 - gamma;
gamma_2 = 1 + gamma;
xhat_MKF = MKF(sys, y, x_0, V_0, gamma_1, gamma_2);
time_MKF = toc/T;
disp(['Avg Time of Moment: ' num2str(time_MKF)]);
err_MKF = sum((x-xhat_MKF).^2,1);
RMSE_MKF = sqrt(mean(err_MKF));
disp(['-------RMSE of MKF: ' num2str(RMSE_MKF)]);


%% Save Data
name = [];
if isTimeVariant
    name = 'RANDOM';
else
    name = 'FIXED';
end
name = ['FAKE_' name '_' num2str(coeff_alpha)];
if exist(name,'var')
    delete([name '.mat']);
end
save(name);

%% Plot
figure;
smt = 50; %50
% can also try "semilogx" instead of "plot"
plot(smooth(10*log10(err_TMKF), smt), 'k', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_KF), smt), 'r', 'LineWidth', 2);
plot(smooth(10*log10(err_Adapt), smt), 'c', 'LineWidth', 2);
plot(smooth(10*log10(err_Fading), smt), 'LineWidth', 2);
plot(smooth(10*log10(err_HINFTY), smt), 'LineWidth', 2);
plot(smooth(10*log10(err_UB), smt), 'LineWidth', 2);
plot(smooth(10*log10(err_UI), smt), 'LineWidth', 2);
% plot(smooth(10*log10(err_SPU), smt), 'LineWidth', 2);
plot(smooth(10*log10(err_SNKF), smt), 'g', 'LineWidth', 2);
plot(smooth(10*log10(err_tauKF), smt), 'm', 'LineWidth', 2);
plot(smooth(10*log10(err_WKF), smt), 'LineWidth', 2);
plot(smooth(10*log10(err_MKF), smt), 'b', 'LineWidth', 2);

% Plot Format
set(gca, 'FontSize', 18);
ylabel('Estimation Error (dB)','FontSize', 20, 'Interpreter', 'latex');
xlabel('Time Step','FontSize', 20, 'Interpreter', 'latex');
leg1 = legend({'TMKF', 'KF', 'Adaptive', 'Fading', '$H_{\infty}$', 'UB', 'UI', 'SNKF', '$\tau$-KF', 'WKF', 'MKF'}, 'Interpreter', 'latex', 'Location', 'northwest');

set(leg1,'FontSize',15);


%% Show Statictics
disp(['Avg Time of TMKF: ' num2str(time_TMKF)]);
disp(['-------RMSE of TMKF: ' num2str(RMSE_TMKF)]);
disp(['Avg Time of KF: ' num2str(time_KF)]);
disp(['-------RMSE of KF: ' num2str(RSME_KF)]);
disp(['Avg Time of Adaptive: ' num2str(time_Adapt)]);
disp(['-------RMSE of Adaptive: ' num2str(RMSE_Adapt)]);
disp(['Avg Time of Fading: ' num2str(time_Fading)]);
disp(['-------RMSE of Fading: ' num2str(RMSE_Fading)]);
disp(['Avg Time of HINFTY: ' num2str(time_HINFTY)]);
disp(['-------RMSE of HINFTY: ' num2str(RMSE_HINFTY)]);
disp(['Avg Time of UB: ' num2str(time_UB)]);
disp(['-------RMSE of UB: ' num2str(RMSE_UB)]);
disp(['Avg Time of UI: ' num2str(time_UI)]);
disp(['-------RMSE of UI: ' num2str(RMSE_UI)]);
% disp(['Avg Time of SPU: ' num2str(time_SPU)]);
% disp(['-------RMSE of SPU: ' num2str(RMSE_SPU)]);
disp(['Avg Time of SNKF: ' num2str(time_SNKF)]);
disp(['-------RMSE of SNKF: ' num2str(RMSE_SNKF)]);
disp(['Avg Time of tauKF: ' num2str(time_tauKF)]);
disp(['-------RMSE of tauKF: ' num2str(RMSE_tauKF)]);
disp(['Avg Time of WKF: ' num2str(time_WKF)]);
disp(['-------RMSE of WKF: ' num2str(RMSE_WKF)]);
disp(['Avg Time of Moment: ' num2str(time_MKF)]);
disp(['-------RMSE of MKF: ' num2str(RMSE_MKF)]);