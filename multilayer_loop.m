%=========================================================
% Name: multilayer_loop.m
% What: compute B field for ZS one turn
% Required: single_loop.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = multilayer_loop(m,m,m,m,(none),A)
%=========================================================
function [Bz,Bx,By] =  multilayer_loop(z,Zp,R,dr,n,I0)
    Bz = zeros(1,numel(z));
    Bx = zeros(1,numel(z));
    By = zeros(1,numel(z));
    for i=1:n
        Rn = R + (i-1)*dr;
        [bbz,bbx,bby] = single_loop(z,Zp,Rn,I0);
        Bz = Bz + bbz;
        Bx = Bx + bbx;
        By = By + bby;
    end
end
