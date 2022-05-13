%=========================================================
% Name: R_scatt_MOT.m
% What: Compute the scattering rate and return the scattering rate and the
% detuing
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: input unit should be SI unit
%=========================================================
function R = R_scatt_MOT(B,v,s0,Delta_0)

u_Bohr_eff = 9.274E-24; %J/T effective Bohr magneton
Gamma = 2*pi*5872400; % rad/s Natural linewidth
h_bar = 6.626E-34/(2*pi); % redunced Planck constant
Lambda = 670.977338E-9; %m Resonance wavelength
%% calculate parameters
k = 2*pi/Lambda;

Del1 = Delta_0 + k*v - u_Bohr_eff*B/h_bar; %  detuning 1
Del2 = Delta_0 - k*v + u_Bohr_eff*B/h_bar; %  detuning 2
%% Return 
R = (Gamma/2)*(s0/(1+s0+4*(Del1^2/Gamma^2)))-(Gamma/2)*(s0/(1+s0+4*(Del2^2/Gamma^2))); % 1/s total scattering rate
end