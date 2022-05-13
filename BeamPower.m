%=========================================================
% Name: BeamPower.m
% What: Compute laser beam properties
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit for input / output [mW,cm] = BeamPower((none),(none),cm,(none))
%=========================================================

function [Power,waist] = BeamPower(n,eta,r,alpha)
    I_sat = 2.54; %mW/cm^2
    I_min = n*eta/(1-eta)*I_sat; %mW/cm^2
    waist = r/sqrt(-0.5*log(alpha)); % cm
    Power = pi*I_min*waist^2/(2*alpha); % mW
end