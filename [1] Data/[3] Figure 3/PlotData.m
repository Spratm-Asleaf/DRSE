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
load DIFFERENT_GAMMA;     % Figure 3.

%% Plot
figure;
smt = 50; %50
% can also try "semilogx" instead of "plot"
plot(smooth(10*log10(err_TMKF), smt), 'k', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_KF), smt), 'r', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_1), smt), 'm', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_2), smt), 'b', 'LineWidth', 2);
plot(smooth(10*log10(err_MKF_3), smt), 'g', 'LineWidth', 2);


% Plot Format
set(gca, 'FontSize', 18);
ylabel('Estimation Error (dB)','FontSize', 20, 'Interpreter', 'latex');
xlabel('Time Step','FontSize', 20, 'Interpreter', 'latex');
leg1 = legend({'TMKF', 'KF', 'MKF (1.005)', 'MKF (1.02)', 'MKF (1.05)'}, 'Interpreter', 'latex', 'Location', 'northwest');

set(leg1,'FontSize',15);
