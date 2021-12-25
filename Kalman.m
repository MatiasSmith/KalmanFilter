                                                             
clc;
clear;
close all;
formatlatex
format longG
warning('off',  'all')
mkdir('./Figures')
warning('on',  'all')
%% Problem 1. KF
close all
% Load workspaces
load OCV_table.mat
load OCV_slope_table.mat
load IV_data_linear.mat

T = 0.1;
RO = 0.01;
Rc = 0.015;
Ccap = 2400;
Cbat = 18000;
VocO = 3.435;
alp = 0.65;
n =  length(t);
                                    % System information
Qk = [2.5E-7,0;0,0];
% Qk = (0.0005, 0; 0,0];
Rk = 1E-4;

% Continuous-time SS
Ac = [0, 0; ... 
    0, -1 / (Ccap * Rc)];
Bc = [-1 / Cbat, 1 / Ccap].';
Cc = [alp, -1];
Dc = -R0;
sysc = ss(Ac, Bc, Cc, Dc);

% Discrete-time SS
A = expm(Ac*T);
B = [T * Bc(1, 1), Rc * exp(-T / (Rc * Ccap)) * (exp(T / (Rc * Ccap)) -1)].';
C = Cc;
D = Dc;

%Initial condition
xhatKf = zeros(2,n); %2xn array of zeros
xhatKf(1) = 1;
SOCdr = zeros(1,n);
SOCdr(1) = 1;
Pkf = zeros(2,2,n);
%Pkf(1:2, 1:2, 1) = Rk * eye(2);

%%% b. Kalman Filter and open loop estimate
for ii = 1:n -1
    [xhatKF(:, ii + 1), Pkf(:, :, ii+1)] = kf(A, B, C, D, Qk, Rk, xhatKF(:, ii), V(ii)-Voc0)  %NOT SURE
    SOCdr(ii+1) = SOCdr(ii) + B(1) * I(ii); %dead reckoning
end

SOChat = xhatKF(1, :);

%Plotting
figure
plot(t,  SOC_act,  'DisplayName',  "Actual SOC")
hold on
plot(t,  SOChat,  'DisplayName', "Kalman Filter")
plot(t,  SOCdr,  'DisplayName',  "Open Loop")
xlim([0, 9E3])
legend
title(["SOC Estimate Using Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])
xlabel("time (s)")
ylabel("SOC")
saveas(gcf,  "./Figures/lb.jpg")

%%% c.  Steady State Kalman Filter
[Pinf,  Kinf,  ~] = idare(A.',  C.',  Qk, Rk);
% Initial  condition
xhatKFss = zeros(2, n);
xhatKFss(1) = 1;
for ii = 1:n - 1
    xhatKFss(:, ii+1) = sskf(A, B, C, D, xhatKF(:, ii), V(ii)-Voc0, I(ii), Kinf.');
end

SOChatss = xhatKFss(1, :);

%Plotting
figure
plot(t,  SOC_act,  'DisplayName',  "Actual SOC")
hold on
plot(t,  SOChatss,  'DisplayName',  "Steady-State KF")
plot(t,  SOCdr,  'DisplayName',  "Open Loop")
xlim([0, 9E3])
legend
title["SOC Estimate Using Steady-State Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"]) 
xlabel("time (s)")
ylabel("SOC")
Saveas(gcf, "./Figures/1c.jpg")

figure('Position', [0.01, 0.01, 0.3, 0.3])
stopind = 700;

plot(t(1:stopind), ones(stopind, 1)*Pinf(1,1), 'DisplayName', "Steady-State Covariance")
hold on
plot(t(1:stopind), squeeze(Pkf(1, 1, 1:stopind)), 'DisplayName', "Time-Variant Covariance")
legend

saveas(gcf, "./Figures/1c_cov.jpg")

%%% d. SOC estimation errors
e = SOC_act - SOChat.';
mu = 0;
sig = sqrt(Pkf(1, 1, end));

histcomparison(e, t, sig, mu)
saveas(gcf, "./Figures/1hist.jpg")

%%Problem 2. EKF
close all
clear xhatKF SOChat
load IV_data_nonlinear.mat

dt = T;

%Nonlinear Voc(SOC) function
Voc0 = interp1(soc_intpts_OCV, OCV_intpts, 0);
Vocfn = @(SOC) interp1(soc_intpts_OCV_slope, OCV_slope_intpts, SOC, 'pchip') - Voc0;

%Nonlinear Vocdot(SOC) function
Vocdot = @(SOC) interp1(soc_intpts_OCV_slope, OCV_slope_intpts, SOC, 'pchip')

%Input function
ufn = @(t) interp1((1:n)*T, I, t);

% Continuous-time nonlinear battery model function
fc = @(t, x, u) [-u / Cbat; ...
    1 / Ccap * ufn(t) - 1 / (Ccap * Rc) * x(2)];

%Continuous-time output equation
hc = @(t, x, u) Vocfn(x(1)) - x(2) - R0 * u;

%Function to update states by 1 time step
fd = @(t, x, u) rk4s(fc, t, x, u, dt);

%Putting fd in a more general form:
ff = @(~, t, x, u) fd(t, x, u);

%Jacobian of state equations
Ap = [0, 0; 0, -1 / (Ccap * Rc)];

% Jacobian of output function
Cp = @(x) [Vocdot(x(1)), 0];

%Discretize the Jacobians
Ffn = @(varargin) expm(Ap * dt);
Hfn = @(~, x, ~) Cp(x);

% Initial conditions of estimate and covariance matrix
xhatEKF = zeros(2, n);
xhatEKF(1) = 1;
Pekf = zeros(2, 2, n);
Pekf(1:2, 1:2, 1) = Rk * eye(2);

tic
for ii = 2:n
    [xhatEKF(:, ii), Pekf(:, :, ii)] = ekf(t(ii-1), t(ii), xhatEKF(:, ii-1), ...
        Pekf(:, :, ii-1), ufn, V(ii)-Voc0, ff, hc, Ffn, Hfn, Qk, Rk);
    if mod(ii, 20) == 0
        clc
        fprintf("Running Time Step")
        fprintf("Current Run Time")
    end
end

SOChat = xhatEKF(1,:);

%Plotting
figure
plot(t, SOC_act, 'DisplayName', "Actual SOC")
hold on
plot(t, SOChat, 'DisplayName', "Extended Kalman Filter")
plot(t, SOCdr, 'DisplayName', "Open Loop")
xlim([0, 9E3])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])

% Page 14
xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/2ekf.jpg")

e = SOC_act - SOChat.';
mu = 0;
sig = sqrt(Pekf(1, 1, end-1));

histcomparison(e, t, sig, mu)
saveas(gcf, "./Figures/2hist.jpg")

%% Problem 3. KF on Nonlin Data
close all
clear xhatSOChat
load IV_data_nonlinear.mat
Voc0 = 3.435;

%%% Run the Kalman filter on the nonlinear data
%initial condition
xhatKF = zeros(2,n);
xhatKF(1) = 1;
Pkf = zeros(2, 2, n);
Pkf(1:2, 1:2, 1) = RK * eye(2);

for ii = 1;n - 1
    [xhatKF(:, ii + 1), Pkf(:,:, ii + 1)] = kf(A, B, C, D, Qk, Rk, xhatKF(:, ii), V(ii)-Voc0) %NOT SURE
end

SOChat = xhatKF(1,:);

figure
plot(t, SOC_act, 'DisplayName', "Actual SOC")
hold on
plot(t, SOChat, 'DisplayName', "Kalman Filter")
plot(t, SOCdr, 'DisplayName', "Open Loop")
xlim([0, 9E3])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs time"])
xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/3kf.jpg")

e = SOC_act - SOChat.';
mu = 0;
sig = sqrt(Pkf(1, 1, end));

histcomparison(e, t, sig, mu)
saveas(gcf, "./Figures/3hist.jpg")
%% kf FUNCTION DEFINITION
% Runs a single step of the Kalman Filter
function [xhat, P] = kf(A, B, C, D, Q, R, xhat, y, u, P)

%Propagation
xhat = A * xhat + B * u;
P = A * P * A.' + Q;

%Kalman Gain
Lk = P * C.' * (C * P * C.' + R)^-1;

%Correction
xhat = xhat + Lk * (y - C * xhat - D * u);
P = P - Lk * C * P;
end

%%sskf FUNCTION DEFINITION
% Runs a single step of the steady state Kalman Filter
function xhat = sskf(A, B, C, D, xhat, y, u, Kinf)

% Propagation
xhat = A * xhat + B * u;

%Correction
xhat = xhat + Kinf * (y - C * xhat - D * u);

end

%% euls FUNCTION DEFINITION
% Runs a single step of Euler integration
function xk = euls(~, ~, x, ~, dt)
% Inputs:
% f: function xdot = f(x)
% t: current time
% x: current state
% u: function u(t) that governs the input
% dt: time step
%
% Outputs:
% xk: states at next time step

xk = x + dt*x;
end

%% rk2s FUNCTION DEFINITION
% Runs a single step of the Runge_Kutta second order ode solver
function xk = rk2s(f, t, x, u, dt)
% Inputs:
% f: function xdot = f(x)
% t: current time
% x: current state
% u: function u(t) that governs the input
% dt: time step
%
% Outputs:
% xk: states at next time step

d1 = f(t, x, u(t));
d2 - f(t+0.5*dt, x+0.5*dt*d1, u(t + 0.5 * dt));
xk = x+dt / 2 * (d1 + d2);

end

 %%rk4s FUNCTION DEFINITION
 % Runs a single step of the Runge-Kutta fourth order ode solver
 function xk = rk4s(f, t, x, u, dt)
% Inputs:
% f: function xdot = f(x)
% t: current time
% x: current state
% u: function u(t) that governs the input
% dt: time step
%
% Outputs:
% xk: states at next time step

d1 = f(t, x, u(t));
d2 = f(t+0.5*dt, x+0.5*dt*d1, u(t + 0.5 * dt));
d3 = f(t+0.5*dt, x+0.5*dt*d2, u(t + 0.5 * dt));
d4 = f(t+dt, x+dt*d3, u(t + dt));
xk = x + dt / 6 * (d1 + 2 * d2 + 2 * d3 + d4);

end

%% ekf FUNCTION DEFINITION
function [xhat, P] = ekf(tkm1, tk, xhat, P, ufn, y, f, h, Ffn, Hfn, Q, R)

u = ufn(tk);

F = Ffn(tkm1, tk, xhat, u);

xhat = f(tkm1, tk, xhat, ufn);
P = F * P * F.' + Q;

H = Hfn(tk, xhat, u);

Lk = P * H.' * (H * P * H.' + R)^-1;

xhat = xhat + Lk*(y - h(tk, xhat, u));

P = P - Lk * H * P;
end

%% histcomparison FUNCTION DEFINITION
function histcomparison(e, t, sig, mu)
figure
histogram(e, 100, 'Normalization', 'pdf')
hold on

% s = trapz(t(1:end-1), abs(e(1:end-1)));

x = linspace(-5*sig, 5*sig, 1001);
% f = 1/ (sqrt(2 * pi) * sig) * exp(-(x - mu).^2./(2*sig^2));
f = normpdf(x, mu, sig);

plot(x, f)
xlim([-max(abs(e))*1.25, max(abs(e))*1.25])

title("SOC Estimation Error")
ylabel("Probability Density")
xlabel("Estimation Error")
end

%% formatlatex FUNCTION DEFINITION
function formatlatex
%latex figure formatting function
%Color plots version

reset(groot)
% FINISH THIS LINE set(groot, 'defaultaxeslinestyleorder', {'-', ':', '-.', '--'}, 'defaultaxesnextplot',   adc
set(groot, 'defaulttextinterpreter',  'latex')
set(groot, 'defaultcolorbarticklabelinterpreter', 'latex')
set(groot, 'defaultfigureunits',  'normalized')
set(groot, 'defaultfigureposition', [0.01, 0.01, 0.3, 0.4])
set(groot, 'defaultaxesticklabelinterpreter', 'latex')
set(groot, 'defaultlegendinterpreter',  'latex')
set(groot, 'defaultaxesfontsize', 24)
set(groot, 'defaultaxeslinewidth', 1)

set(groot, 'defaultscattermarker',  'o')
set(groot, 'defaultscatterlinewidth', 2)
set(groot, 'defaultlinemarkersize', 15)
set(groot, 'defaultlinelinewidth', 2.5)
set(groot, 'defaultAxesXgrid',  'on',  'defaultAxesYgrid', 'on', idefaultAxesZgridl,'on')
set(groot, 'defaultAxesGridLineStyle',  '-.')
set(groot, 'defaultAxesXlim', [0, 2 * pi]);
set(groot, 'defaultAxesYlim', [-0.6, 0.6]);
end
