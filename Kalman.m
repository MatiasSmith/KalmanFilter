                                                             
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
Ac = [0, 0; ... 
    0, -1 / (Ccap * Rc)];
Bc = [-1 / Cbat, 1 / Ccap].;
Cc = [alp, -1];
Dc = -R0;
Sysc = ss(Ac, Bc, Cc, Dc);

% Discrete-time SS
A = expm(Ac*T);
B = [T * Bc(1, 1), Rc * exp(-T / (Rc * Ccap)) * (exp(T / (Rc * Ccap)) -1)].;
C = Cc;
D = Dc;

%Inidtial condition
xhatKf = zeros(2,n);
xhatKf(1) = 1;
SOCdr = zeros(1,n);
SOCdr(1) = 1;
Pkf = zeros(2,2,n);
%Pkf(1:2, 1:2, 1) = Rk * eye(2);

%%% b. Kalman Filter and open loop estimate
for ii = 1:n -1
    [xhatKF(:, ii + 1), Pkf(:, :, ii+1)] = kf(A, B, C, D, Qk, Rk, xhatKF(:, ii), V(ii)-VocO)  %NOT SURE
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
title(["SOC Estimate Using Kalman Filter, Open Loop Estimate and", "Actual SOC vs.  Time"])
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

e = SOC_act - SOChat.';
mu = 0;
sig = sqrt(Pkf(1, 1, end));

histcomparison(e, t, sig, mu)
saveas(gcf, "./Figures/1hist.jpg")

close all
clear xhatKF SOChat
load IV_data_nonlinear.mat

dt = T;

Voc0 = interp1(soc_intpts_OCV, OCV_intpts, 0);
Vocfn = @(SOC) interp1(soc_intpts_OCV_slope, OCV_slope_intpts, SOC, 'pchip');

ufn = @(t) interp1((1:n)*T, I, t);

fc = @(t, x, u) [-u / Cbat; ...
    1 / Ccap * ufn(t) - 1 / (Ccap * Rc) * x(2)];

hc = @(t, x, u) Vocfn(x(1)) - x(2) - R0 * u;

fd = @(t, x, u) rk4s(fc, t, x, u, dt);

ff = @(~, t, x, u) fd(t, x, u);

Ap = [0, 0; 0, -1 / (Ccap * Rc)];

Cp = @(x) [Vocdot(x(1)), 0];

Ffn = @(varargin) expm(Ap * dt);
Hfn = @(~, x, ~) Cp(x);

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

figure
plot(t, SOC_act, 'DisplayName', "Actual SOC")
hold on
plot(t, SOChat, 'DisplayName', "Extended Kalman Filter")
plot(t, SOChat, 'DisplayName', "Open Loop")
xlim([0, 9E3])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])


14

xlabel("time                                         (8)")
ylabel("SOC")
saveas(gcf,  "./Figures/2ekf.jpg")
e                                                    = SOC_act                                            - SOChat.'                 ;
mu = 0;
sig                                                  =  sqrt(Pekf(1,                                      1,  end-1));
histcomparison(e,  t,  sig, mu)
saveas(gcf,  "./Figures/2hist.jpg")
%% Problem 3.  KF on Nonlin Data
close all
clear xhat SOChat
load  IV_data_nonlinear.mat
VocO                                                 = 3.435;
                                                                                                                                                                     XXX Run  the Kalman filter on the nonlinear             data
                                                                                                          X Initial  condition
                                                                                                          xhatKF = zeros(2,  n);
xhatKF(1)                                                                                                 =                          1;
Pkf                                                                                                       = zeros(2,                                       2,  n);
Pkf(1:2,                                                                                                  1:2,                                        1)   = Rk      * eye(2);
for  ii                                              =                                                    1:n                        -                1
                                                                                                          [xhatKF(:,  ii                                   +         1),  Pkf(:,                                   +   1)]   = kf(A,   B,   C,   D,   Qk,   Rk,   xhatKF(:,   ii),   V(ii)-Vc
end
SOChat                                               = xhatKF(1,:);
figure
plot(t,  SOC_act,   'DisplayName',  "Actual SOC")
hold  on
plot(t,  SOChat,  'DisplayName',  "Kalman Filter")
plot(t,  SOCdr,  'DisplayName',  "Open Loop")
xlim([0,                                             9E3])
legend
                                                     title(["SOC Estimate Using Extended Kalman Filter,   Open Loop Estimate and",   "Actual SOC vs
xlabel("time                                         (s)")
ylabel("SOC")
saveas(gcf,  "                                       /Figures/3kf.jpg")
histcomparison(e,  t,  sig, mu)
15

saveas(gcf,  "./Figures/3hist.jpg")
XX kf FUNCTION DEFINITIW
X Runs  a single  step  of  the Kalman
function  Lshat,  P]                                        = kf(A,  B,  C,  D,  Q,  R,  xhat,  y,  u,  P)
X Kalman Gain
Lk = P * C.'                                                *                                                (C                            * P                                             * C.'   + R)--1;
XX  sskf FUNCTION DEFINITION
X Runs  a single step  of  the steady state Kalman filter
function xhat                                               = sskf(A,  B,  C,                                0,  xhat,  y,  u,  Kinf)
X Propagation
'that                                                       = A *  xhat + B                                  u;
X Correction
                                                            xhat = xhat                                      + Kinf                        *                                               (y      - C                                       xhat   -- D   -  u)   ;
end
                                                                                                             XX etas FUNCTION DEFINITION
                                                                                                                                           X Runs  a  single  step  of Euler integration
                                                            function xk                                      = euls(-,                                                                             x,                                        dt)
                                                            X Inputs:
                                                            f:                                               function xdot                                                                         = f(x)
%                                                           t:                                               current  time
X                                                           z:                                               current  state
X                                                           u:                                                                                                                                     function u(t)  that governs  the  input
dt:  time  step
X
Outputs:
X                                                           xk:  states  at next  time step
xk = x + dt*x;
end
XX rk2s FUNCTION DEFINITION

                                                                                                                                                                                            X Runs  a single  step  of  the Runge-Kutta second                                                                                    order   ode   solver
                                 function xk                                                                                                                                                = rk2s(f,  t,  x,  u,  dt)
                                 X Inputs:
X                                f:                                                 function xdot                                                                                                                                                = f(x)
X                                t:                                                 current  time
X                                x:                                                 current  state
X                                u:                                                                                                                                                                                                              function u(t)  that governs the input
X                                dt:            time step
X
                                 X Outputs:
                                                                                                                                                        xk:  states  at next  time step
dl                                              = f(t,  x,  u(t));
                                                                                                                  d2 = f(t+0.5*dt,  x+0 b*dt*dl, u(t                                                                                                                                                       +   0.5     *  dt));
                                 xk = x + dt    /                                   2                                                                   (d1                                                                                      + d2);
end
                                                                                    XX rk4s FUNCTION DEFINITION
                                                                                                                                                                                            X Runs  a single  step  of the Runge -Kutta fourth                                                                                    order   ode   solver
                                 function xk                                                                                                                                                = rk4s(f,  t,  x,  u,  dt)
                                 % Inputs:
                                 f:                                                 function xdot                                                                                                                                                = f(x)
X                                t:                                                 current  time
                                 m:                                                 current  state
X                                u:                                                                                                                                                                                                              function u(t)  that  governs  the  input
X                                               dt:  time  step
                                 X Outputs:
                                                                                                                                                        xk:  states  at  next  time  step
dl                                              = f(t,  x,  u(t));
                                                                                                                  d2 = f(t+0.5*dt,  x+0.5*dt*dl,  u(t                                                                                                                                                      +   0.5     * dt));
                                                                                                                  d3 = f(t+0.5*dt,  x+0.5*dt*d2,  u(t                                                                                                                                                      +   0.5     *  dt));
                                                d4 = f(t+dt,  x+dt*d3,  u(t                                                                                                                                                                                                                 + dt));
                                 xk = x +  dt   /                                   6                             *                                     (dl                                 +                                                    2                                          * d2      +    2   *  d3   +  d4);
end
X% ekf FUNCTION DEFINITION
function                         [xhat,  P]     = ekf(tkml,  tk,  xhat,  P,  ufn,   y,                            f,                                    h,                                  Ffn,                                                 Hfn,                                       Q,        R)
u = ufn(tk);
F = Ffn(tkml,  tk,  xhat,  u);
17

xhat                                                                                                                    = f(tkml,  tk,  xhat,  ufn);
P = F                                          * P                                            *- F.'                                                              + Q;
                                               H = Hfn(tk,  xhat,  u);
Lk = P                                         * H.'                                          *                         (H                                        * P                              * H.'                 + R)--1;
xhat                                           = xhat                                                                   + Lk*(y                                                                    - h(tk,  xhat, u));
P = P                                          - Lk                                           * H                       * P;
end
XX histcomparison FUNCTION DEFINITION
function histcomparison(e, t, sig, mu)
figure
                                               histogram(e,                                                                                                       100,  'Normalization',  'pdf')
hold on
% s                                                                                                                     =  trapz(t(1:end-1),  abs(e(1:end-1)));
                                                                                              x = linspace(-5*sig,      5*sig,                                    1001);
X f =                                          1                                              /  (sort(2                * pi)                                     * sig)                           * exp(-(x             -                        mu).   "2./(2   *   sig-2));
f                                                                                             = normpdf(x, mu,  sig);
plot(x,  f)
xlim([-max(abs(e))*1.25,  max(abs(e))*1.25])
title("SOC Estimation Error")
ylabel("Probability Density")
xlabel("Estimation Error")
end
formattater FUNCTION DEFINITION
function formatlatex
/ Latex figure formatting function
X                                              '  Color plots version
reset(groot)
set(groot,                                     'defaultaxeslinestyleorder',                                                                                                                                              'defaultaxesnextplot',   'adc
set(groot,                                     'defaulttextinterpreter',  'latex')
set(groot,                                     'defaultcolorbarticklabelinterpreter'                                                                              , 'latex')
set(groot,                                     'defaultfigureunits',  'normalized')
set(groot,                                     idefaultfigureposition',                       [0.01,                    0.01,                                     0.3,                             0.4])
set(groot,                                     'defaultaxesticklabelinterpreter',  'latex')
set(groot,                                     idefaultlegendinterpreter',  'latex')
set(groot,                                     'defaultaxesfontsize',                         24)
set(groot,                                     'defaultaxeslinewidth',                        1)
18

set(groot,   'defaultscattermarker',  'o')
set(groot,   'defaultscatterlinewidth',                                                                                             2)
set(groot,   'defaultlinemarkersize',                                                                                               15)
set(groot,   'defaultlinelinewidth',                                                                                                2.5)
set(groot,                                                                        'defaultAxesXgrid',  'on',  'defaultAxesYgrid',             on   , idefaultAxesZgridl,  'on')
set(groot,                                   idefaultAxesGridLineStyle',  '-.')
set(groot,   idefaultAxesXlimi,              [0,                                  2                                                 * pi]);
set(groot,   IdefaultAxesYlim',                                                   [-0.6,                                            0.6]);
end
