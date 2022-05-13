%=========================================================
% Name: MOT_multilayer_loop.m
% What: compute B field for MOT one turn
% Required: MOT_single_loop.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = MOT_multilayer_loop(m,m,m,(none),m,m,A)
%=========================================================
function [Bz,Bx,By] =  MOT_multilayer_loop(z,Zp,d,n,dr,R,I0)
    Bz = zeros(1,numel(z));
    Bx = zeros(1,numel(z));
    By = zeros(1,numel(z));
    for i=1:n
        Rn = R + (i-1)*dr;
        [bbz,bbx,bby] = MOT_single_loop(z,Zp,d,Rn,I0);
        Bz = Bz + bbz;
        Bx = Bx + bbx;
        By = By + bby;
    end
end