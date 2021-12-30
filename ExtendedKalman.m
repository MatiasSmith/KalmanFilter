clc;
clear;
close all;
%formatlatex
format longG
warning('off',  'all')
mkdir('./Figures')
warning('on',  'all')

V = [10, 11, 12, 13, 15, 14, 17, 18]
t = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
I = [5, 5.5, 6, 6.5, 7.5, 7, 8.5, 9]

T = 0.1;   %Sampling Period
R0 = 0.01;
Rc = 0.015;
Ccap = 2400;
Cbat = 18000;
n = length(t);

% System information
Qk = [2.5E-7,0;0,0];
% Qk = (0.0005, 0; 0,0];
Rk = 1E-4;

soc_intpts_OCV = t;
OCV_intpts = V;
soc_intpts_OCV_slope = t;
OCV_slope_intpts = V;

dt = T;

%Nonlinear Voc(SOC) function
Voc0 = interp1(soc_intpts_OCV, OCV_intpts, 0);  %Using graph, predict value at t=0
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
SOC_act = V(1, :);
SOCdr = zeros(1,n);
SOCdr(1) = 1;

tic
for ii = 2:n        %Declaration on line 319
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
xlim([0, 1])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])

% Page 14
xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/2ekf.jpg")

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

%% ekf FUNCTION DEFINITION     %Called on Line 164
function [xhat, P] = ekf(tkm1, tk, xhat, P, ufn, y, f, h, Ffn, Hfn, Q, R)

u = ufn(tk);

F = Ffn(tkm1, tk, xhat, u);

xhat = f(tkm1, tk, xhat, ufn);   %Prediction
P = F * P * F.' + Q;             %Prediction

H = Hfn(tk, xhat, u);            

Lk = P * H.' * (H * P * H.' + R)^-1;    %Kalman Gain

xhat = xhat + Lk*(y - h(tk, xhat, u));  %Correction
P = P - Lk * H * P;                     %Correction
end