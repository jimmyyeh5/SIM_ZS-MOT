%=========================================================
% Name: Vrec2.m
% What: generate random direction of given velocity
% Required: randvec.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
%=========================================================

function [vx,vy,vz] = Vrec2(v)
    [x,y,z] = randvec();
    vx = v*x;
    vy = v*y;
    vz = v*z;
end