%=========================================================
% Name: MOT_coil.m
% What: compute B field for MOT
% Required: MOT_multilayer_loop.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = section((none),m,(none),m,m,m,A,m,m)
%=========================================================
function [Bz,Bx,By] = MOT_coil(noL,dL,n,dr,z,R,I0,Xp,Zp)
    Bz = zeros(1,numel(z));
    Bx = zeros(1,numel(z));
    By = zeros(1,numel(z));
    
    [bbz,bbx,bby] =  MOT_section(noL,dL,n,dr,z,R,I0,Xp,Zp);
    Bz = Bz + bbz;
    Bx = Bx + bbx;
    By = By + bby;
    
    [bbz,bbx,bby] =  MOT_section(noL,dL,n,dr,z,R,-I0,-Xp,Zp);
    Bz = Bz + bbz;
    Bx = Bx + bbx;
    By = By + bby;
end