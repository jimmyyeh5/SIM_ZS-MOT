%=========================================================
% Name: MOT_Coil_field.m
% What: Construct MOT coil field
% Required: MOT_coil.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/17
%=========================================================

%%Construct magnetic field and compare to the ideal one
close all
clc;
clear;

MOT_Beam_radius = 10E-3;
Loop = 1000;

z = linspace(-MOT_Beam_radius,MOT_Beam_radius,Loop);
%% Wire parameters
MOT_dL = 3E-3; % in m
MOT_dr = 3E-3; % in m

%% Fixed parameter
MOT_Rmin = 101E-3; %in m
MOT_Xp = 45E-3; %in m
MOT_Zp = 0; %in m

%% MOT Free parameters
%geometric parameters

MOT_n = 10; % in turns max 11
% free parameters
MOT_noL = 20; % in turns 20
MOT_I = 40; % in A Threshold 120 A36
%%

MOT_Bz = MOT_coil(MOT_noL,MOT_dL,MOT_n,MOT_dr,z,MOT_Rmin,MOT_I,MOT_Xp,MOT_Zp);

figure
plot(z,MOT_Bz)
csvwrite('B_MOT_coil.csv',MOT_Bz);
csvwrite('Z_MOT_coil.csv',z);
xlabel('Position(m)');
ylabel('B field (T)');
