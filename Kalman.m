clc;
clear;
close all;
%formatlatex
format longG
warning('off',  'all')
mkdir('./Figures')
warning('on',  'all')

%Parameters
T = 0.1;
R0 = 0.01;
Rc = 0.015;  
Ccap = 2400;
Cbat = 18000;


% System information
Qk = [2.5E-7,0;0,0];
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
errMea = 4.* ones(length, 1);

V = [0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4]
t = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
t = [0:length-1];

%Driving loop
n = length;
for i = 2:n
     [Est(i), errEst(i)] = ekf(errEst(i-1), errMea(i), Est(i-1), Mea(i));
end


constant = TrueVal.*ones(length, 1).*r5'

%Plotting
figure
plot(t, constant, 'DisplayName', "Actual SOC")
hold on
plot(t, Est, 'DisplayName', "Extended Kalman Filter")

ylim = ([72-5, 72+5]);
xlim([0, length+2])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])

xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/2ekf.jpg")

%Calculation
function [nextEst, nextErrEst] = ekf(errEst, errMea, prevEst, Mea)

    %Kalman Gain
    KG = errEst / (errEst + errMea);
    
    nextEst = prevEst + KG*(Mea - prevEst);

    ErrEst = (errMea * errEst) / (errMea + errEst);
    
    nextErrEst = (1-KG)*errEst;

end

