clc;
clear;
close all;
%formatlatex
format longG
warning('off',  'all')
mkdir('./Figures')
warning('on',  'all')

T = 0.1;   %Sampling Period
R0 = 0.01;
Rc = 0.015;  
Ccap = 2400;
Cbat = 18000;


% System information
Qk = [2.5E-7,0;0,0];
% Qk = (0.0005, 0; 0,0];
Rk = 1E-4;

length = 100;
TrueVal = 77;
r1 = (rand(length, 1) * 4) + TrueVal-2;
r5 = [1:-0.01:0.01];
r1 = r1.*r5';
r2 = zeros(length+1,1);
r3 = r2;
r2(1) = 78;   %initial est
r3(1) = 2;    %initial est error

Est = [74, 0, 0, 0, 0];
errEst = [4, 0, 0, 0, 0];
Mea = r1;
%Mea(49) = Mea(49)-50;
%Mea(50) = Mea(50)-50;
%Mea(51) = Mea(51)-50;
errMea = 4.* ones(length, 1);

V = [0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4]
t = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
t = [0:length-1];

%for i = 2:n        % for (int ii = 2; ii < n; i++)  Declaration on line 319
%    [xhatEKF(:, ii), Pekf(:, :, ii)] = ekf(t(ii-1), t(ii), xhatEKF(:, ii-1), ...
%        Pekf(:, :, ii-1), ufn, V(ii)-Voc0, ff, hc, Ffn, Hfn, Qk, Rk);
%n = length(Mea);
n = length;
for i = 2:n
     [Est(i), errEst(i)] = ekf(errEst(i-1), errMea(i), Est(i-1), Mea(i));
end


constant = TrueVal.*ones(length, 1).*r5'
%constant(49) = constant(49)-50;
%constant(50) = constant(50)-50;
%constant(51) = constant(51)-50;
%Plotting
figure
plot(t, constant, 'DisplayName', "Actual SOC")
hold on
plot(t, Est, 'DisplayName', "Extended Kalman Filter")
%plot(t, SOCdr, 'DisplayName', "Open Loop")
ylim = ([72-5, 72+5]);
xlim([0, length+2])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])

% Page 14
xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/2ekf.jpg")

function [nextEst, nextErrEst] = ekf(errEst, errMea, prevEst, Mea)

    KG = errEst / (errEst + errMea);

    nextEst = prevEst + KG*(Mea - prevEst);

 %   ErrEst = (errMea * errEst) / (errMea + errEst);
    nextErrEst = (1-KG)*errEst;

end
%function x = f(

%end
