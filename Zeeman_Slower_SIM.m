%=========================================================
% Name: Zeeman_Slower_SIM.m
% What: 1. Simulate the trajectory of atoms with different initial velocity
%          with given B-profile and laser beam properties.
%       2. Calculate final velocity of ZS, and find the caputure 
%          velocity for MOT.
% Required: B_coil.csv / Z_coil.csv / BeamPower.m / B_ideal.m /
%           trajectory.m
% Note: Simulate using ideal field(flag = 0) or coil field(flag = 1)
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
%=========================================================

clear;
close all;
clc;

%% Parameter setup
flag = 1; % choose to use ideal field(flag = 0) or coil field(flag = 1)
Detuning = -9; % Gamma
eta = 0.4;% Safety parameter
v_c = 1; % m/s Capture velocity
nI = 10; % intensity parameter (ZS laser beam power)
MOT_beam_radius = 10E-3; %in m 

v_test_i = 100; % m/s initial value for test initial velocity 
v_test_f = 1100; % m/s final value for test initial velocity
v_test_gap = 100; % m/s step size for test initail velocity

v_0 = 1000; % m/s velocity to generate ideal profile

%% ZS laser beam properties
r_min = 0.62;% cm
alpha = 0.37;
%Power = BeamPower(nI,eta,r_min,alpha)
Power = BeamPower(nI,eta,r_min,alpha)% nI = s0

%% calculated parameters
v_test = v_test_i:v_test_gap:v_test_f;
n = numel(v_test);
Vfilename = strings([1:n]);
Zfilename = strings([1:n]);
V_legend = strings([1:n]);


%% calculate B ideal field
B_ideal(Detuning,eta,v_c,v_0);

%% calculate trajectory
i=1;
for v_i = v_test
    trajectory(v_i,v_0,Detuning,eta,v_c,flag,nI,MOT_beam_radius);
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
    Zcm = 100*(Z-Z(end));
    plot(Zcm,V);
    hold on
    legend(V_legend)
end
title('Trajectory (ZS region)')
xlabel('Position(cm)');
ylabel('Velocity (m/s)');
limits = [-67.6 3 0 v_test_f+100];
axis(limits)

