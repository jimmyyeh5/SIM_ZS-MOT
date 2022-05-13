%=========================================================
% Name: MOT_SIM.m
% What: 1. Simulate the trajectory of atoms with different initial velocity
%          with given MOT-B-profile and laser beam properties.
%       2. find the caputure velocity for MOT.
% Required: B_coil.csv / Z_coil.csv / BeamPower_MOT.m / trajectory_MOT.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
%=========================================================
clc;
clear;
close all;
%% Parameter setup
Beam_radius = 10E-3;% MOT beam in m
Detuning = -5; % Gamma
s0 = 1.2;%1.2

v_test_i = 15 ; % m/s initial value for test initial velocity 
v_test_f = 55; % m/s final value for test initial velocity
v_test_gap = 10; % m/s step size for test initail velocity


 Beam_radius_cm = Beam_radius*100;
 alpha = 0.37;
 Power = BeamPower_MOT(s0,Beam_radius_cm,alpha)

%% calculate parameters
v_test = v_test_i:v_test_gap:v_test_f;
n = numel(v_test);
Vfilename = strings([1:n]);
Zfilename = strings([1:n]);
V_legend = strings([1:n]);


%% calculate trajectory
i=1;
for v_i = v_test
    trajectory_MOT(v_i,Detuning,s0,Beam_radius);
    chr = int2str(v_i);
    Vfilename(i) = strcat('V_actual_',chr,'.csv');
    Zfilename(i) = strcat('Z_actual_',chr,'.csv');
    V_legend(i) = strcat('v = ',chr,'m/s');
    i = i + 1;
end

%% make plot
Trajactory = figure;
for j=1:n
    V = csvread(Vfilename(j));
    Z = csvread(Zfilename(j));
    Zcm = Z*100;
    plot(Zcm,V,'.');
    hold on
end

legend(V_legend)

%subject = strcat('Exact ',', Detuning = ',num2str(Detuning),', s0 =', num2str(s0),', Power = ',num2str(Power),'mW', ', r = ',num2str(Beam_radius_cm),'cm');    
subject='Trajectory (MOT region)';
title(subject)
xlabel('Position(cm)');
ylabel('Velocity (m/s)');


