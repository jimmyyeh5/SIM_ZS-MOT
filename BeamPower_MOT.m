%=========================================================
% Name: BeamPower.m
% What: Compute laser beam properties
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/13
% Note: The unit for input / output [mW,cm] = BeamPower((none),cm,(none))
%=========================================================
% [mW,cm] = BeamPower((none),cm,(none))
function [Power,waist] = BeamPower_MOT(s0,r,alpha)
    I_sat = 2.54; %mW/cm^2
    I_min = s0*I_sat; %mW/cm^2
    waist = r/sqrt(-0.5*log(alpha));
    Power = pi*I_min*waist^2/(2*alpha);
end