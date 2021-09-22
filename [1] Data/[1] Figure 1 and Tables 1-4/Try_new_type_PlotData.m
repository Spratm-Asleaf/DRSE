%{
Online supplementary materials of the paper titled 
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
Authored By: Shixiong Wang, Zhongming Wu, and Andrew Lim

@Author: Shixiong Wang
@Date: 3 Nov 2020
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE

Related Paper:
"Robust State Estimation for Linear Systems under Distributional Uncertainty"
%}

clear all;
close all;
clc;

%% Load Data
%TRUE: TRUE Structure
%FIXED: FIXED \Delta_k; RANDOM: RANDOM \Delta_k
%5: \alpha = 5

table = 1;
switch table
    case 1
        load TRUE_FIXED_5;   % Figure 1(a). Table I.
    case 2
        load TRUE_RANDOM_1;  % Figure 1(b). Table II.
    case 3
        load TRUE_RANDOM_5;  % Figure 1(c). Table III.
    case 4
        load TRUE_RANDOM_0;  % Figure 1(d). Table IV.
end

%% Plot
figure;
smt = 50; %50
% can also try "semilogx" instead of "plot"
markersize = 6;
plot(smooth(10*log10(err_TMKF), smt), 'k-.', 'LineWidth', 2,'MarkerSize',markersize); hold on;
plot(smooth(10*log10(err_KF), smt), 'r-o', 'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_Adapt), smt), 'c-x', 'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_Fading), smt), '-+','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_HINFTY), smt), '-*','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_UB), smt), '-s','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_UI), smt), '-d','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_SPU), smt), '-v','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_SNKF), smt), 'g-','Color',[rand rand rand], 'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_tauKF), smt), 'm-^', 'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_WKF), smt), '-p','Color',[rand rand rand],'LineWidth', 2,'MarkerSize',markersize);
plot(smooth(10*log10(err_MKF), smt), 'b-h', 'LineWidth', 2,'MarkerSize',markersize);

% Plot Format
set(gca, 'FontSize', 18);
ylabel('Estimation Error (dB)','FontSize', 20, 'Interpreter', 'latex');
xlabel('Time Step','FontSize', 20, 'Interpreter', 'latex');
leg1 = legend({'TMKF', 'KF', 'Adaptive', 'Fading', '$H_{\infty}$', 'UB', 'UI', 'SPU', 'SNKF', '$\tau$-KF', 'WKF', 'MKF'}, 'Interpreter', 'latex', 'Location', 'northwest');

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
disp(['Avg Time of SPU: ' num2str(time_SPU)]);
disp(['-------RMSE of SPU: ' num2str(RMSE_SPU)]);
disp(['Avg Time of SNKF: ' num2str(time_SNKF)]);
disp(['-------RMSE of SNKF: ' num2str(RMSE_SNKF)]);
disp(['Avg Time of tauKF: ' num2str(time_tauKF)]);
disp(['-------RMSE of tauKF: ' num2str(RMSE_tauKF)]);
disp(['Avg Time of WKF: ' num2str(time_WKF)]);
disp(['-------RMSE of WKF: ' num2str(RMSE_WKF)]);
disp(['Avg Time of Moment: ' num2str(time_MKF)]);
disp(['-------RMSE of MKF: ' num2str(RMSE_MKF)]);
