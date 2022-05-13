%=========================================================
% Name: B_coil_field.m
% What: Construct magnetic field and compare to the ideal one
% Required: B_ideal.m / section.m / MOT_coil.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
%=========================================================

close all
clc;
clear;

%Ideal curve setting
Detuning = -15; % Gamma
eta = 0.4;% Safety parameter
v_c = 1; % m/s Capture velocity of MOT
v_0 = 1000; % m/s Initial velocity (Capure velocity of ZS)

%% constant
    u_Bohr_eff = 9.274E-24; %J/T effective Bohr magneton
    Gamma = 2*pi*5872400; % rad/s Natural linewidth
    h_bar = 6.626E-34/(2*pi); % redunced Planck constant
    a_max = 1.82E6; % m/s^2 maximum decceleration
    Lambda = 670.977338E-9; %m Resonance wavelength
%% calculated parameters
    k = 2*pi/Lambda; % rad/m wave number
    Delta_0 = Detuning*Gamma; % rad/s laser detuning
%% Calculate the ideal field B0 for reference and comparison
l0 = (v_0^2-v_c^2)/(2*eta*a_max); %length of Slower
B_ideal (Detuning,eta,v_c,v_0);
B0 = csvread('B_ideal.csv');
Z0 = csvread('Z_ideal.csv');
Loop = 1000;
z = linspace(0,l0,Loop);
%% Wire parameters
MOT_dL = 3E-3; % in m
MOT_dr = 3E-3; % in m
ZS_dL = 3E-3; % in m
ZS_dr = 3E-3; % in m

%% Fixed parameter
ZS_Rmin = 15E-3; % in m
MOT_Rmin = 101E-3; %in m
MOT_Xp = 45E-3; %in m
MOT_Zp = l0; %in m

%% ZS Free Parameters

% current
ZSI1 = 18.9; % in A
ZSI2 = 13.0; % in A

% number of turns
% in turns
ZS1_n = 16; ZS2_n = 12; ZS3_n = 11; ZS4_n = 10; 
ZS5_n = 9; ZS6_n = 8; ZS7_n = 7; ZS8_n = 6; 
% Z_noL total = 184
ZS1_noL = 23; ZS2_noL = 23; ZS3_noL = 23; ZS4_noL = 23;
ZS5_noL = 23; ZS6_noL = 23; ZS7_noL = 23; ZS8_noL = 20;%35
% space for spacer
d = 1E-3; 

ZS1_I = ZSI1;ZS2_I = ZSI1; ZS3_I = ZSI1; ZS4_I = ZSI1;
ZS5_I = ZSI1; ZS6_I = ZSI1; ZS7_I = ZSI1; 
ZS8_I = ZSI2;

%geometric parameters
% in m
ZS1_Zp = 0 ;
ZS2_Zp = d+ZS1_noL*ZS_dL; ZS3_Zp = d+ ZS2_Zp+ZS2_noL*ZS_dL; ZS4_Zp = d+ZS3_Zp+ZS3_noL*ZS_dL; 
ZS5_Zp = d+ZS4_Zp+ZS4_noL*ZS_dL; ZS6_Zp = d+ZS5_Zp+ZS5_noL*ZS_dL; ZS7_Zp = d+ZS6_Zp+ZS6_noL*ZS_dL; 
ZS8_Zp = d+ZS7_Zp+ZS7_noL*ZS_dL; 

%% MOT Free parameters
%geometric parameters
MOT_n = 10; % in turns max 11
MOT_noL = 20; % in turns 20
MOT_I = 40; % in A Threshold 120 A36

%% Calculate B field for each component
BZS1 = section(ZS1_noL,ZS_dL,ZS1_n,ZS1_Zp,ZS_dr,z,ZS_Rmin,ZS1_I);
BZS2 = section(ZS2_noL,ZS_dL,ZS2_n,ZS2_Zp,ZS_dr,z,ZS_Rmin,ZS2_I);
BZS3 = section(ZS3_noL,ZS_dL,ZS3_n,ZS3_Zp,ZS_dr,z,ZS_Rmin,ZS3_I);
BZS4 = section(ZS4_noL,ZS_dL,ZS4_n,ZS4_Zp,ZS_dr,z,ZS_Rmin,ZS4_I);
BZS5 = section(ZS5_noL,ZS_dL,ZS5_n,ZS5_Zp,ZS_dr,z,ZS_Rmin,ZS5_I);
BZS6 = section(ZS6_noL,ZS_dL,ZS6_n,ZS6_Zp,ZS_dr,z,ZS_Rmin,ZS6_I);
BZS7 = section(ZS7_noL,ZS_dL,ZS7_n,ZS7_Zp,ZS_dr,z,ZS_Rmin,ZS7_I);
BZS8 = section(ZS8_noL,ZS_dL,ZS8_n,ZS8_Zp,ZS_dr,z,ZS_Rmin,ZS8_I);
MOT_Bz = MOT_coil(MOT_noL,MOT_dL,MOT_n,MOT_dr,z,MOT_Rmin,MOT_I,MOT_Xp,MOT_Zp);

B_tot = BZS1 + BZS2 + BZS3 + BZS4 + BZS5 + BZS6 + BZS7 + BZS8 + MOT_Bz;

%% from T,m to G,cm

Z0cm = 100*(Z0 - l0); % cm
zcm = 100*(z - l0);
G = 1E4;
figure
Bideal = plot(Z0cm,B0*G,'--');
Bideal.Color = [0 0.4470 0.7410];
hold on
Btot = plot(zcm,B_tot*G);
Btot.Color = [0.6350 0.0780 0.1840];
hold on
BMOT = plot(zcm,MOT_Bz*G);
BMOT.Color = [0.9290 0.6940 0.1250];
hold on

ZScolor = [0.4660 0.6740 0.1880];
Bzs1 = plot(zcm,BZS1*G);
Bzs1.Color = ZScolor;
hold on;
Bzs2 = plot(zcm,BZS2*G);
Bzs2.Color = ZScolor;
hold on;
Bzs3 = plot(zcm,BZS3*G);
Bzs3.Color = ZScolor;
hold on;
Bzs4 = plot(zcm,BZS4*G);
Bzs4.Color = ZScolor;
hold on;
Bzs5 = plot(zcm,BZS5*G);
Bzs5.Color = ZScolor;
hold on;
Bzs6 = plot(zcm,BZS6*G);
Bzs6.Color = ZScolor;
hold on;
Bzs7 = plot(zcm,BZS7*G);
Bzs7.Color = ZScolor;
hold on;
Bzs8 = plot(zcm,BZS8*G);
Bzs8.Color = ZScolor;
hold on;


legend('B_{Ideal}','B_{ZS}','B_{MOT}','Components')
title('Magnetic Profile')
xlabel('Position(cm)');
ylabel('Magnetic field (G)');

limits = [-68.68 1 0 1100];
axis(limits)

csvwrite('B_coil.csv',B_tot);
csvwrite('Z_coil.csv',z);
