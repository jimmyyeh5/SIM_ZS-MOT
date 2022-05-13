%=========================================================
% Name: section.m
% What: compute B field for ZS section
% Required: multilayer_loop.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = section((none),m,(none),m,m,m,m,A)
%=========================================================
function [Bz,Bx,By] =  section(noL,dL,n,Zp,dr,z,R,I0)

    Bz = zeros(1,numel(z));
    Bx = zeros(1,numel(z));
    By = zeros(1,numel(z));
    for i = 1:noL
        Zn = Zp + (i-1)*dL;
        [bbz,bbx,bby] = multilayer_loop(z,Zn,R,dr,n,I0);
        Bz = Bz + bbz;
        Bx = Bx + bbx;
        By = By + bby;
    end
end