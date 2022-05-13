%=========================================================
% Name: MOT_section.m
% What: compute B field for MOT section
% Required: MOT_multilayer_loop.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = MOT_section((none),m,(none),m,m,m,A,m,m)
%=========================================================
function [Bz,Bx,By] =  MOT_section(noL,dL,n,dr,z,R,I0,Xp,Zp)

    Bz = zeros(1,numel(z));
    Bx = zeros(1,numel(z));
    By = zeros(1,numel(z));
    
    if Xp > 0%(God!! It's all about construction!!)
        for i = 1:noL
        dn = Xp + (i-1)*dL;
            [bbz,bbx,bby] = MOT_multilayer_loop(z,Zp,dn,n,dr,R,I0);
            Bz = Bz + bbz;
            Bx = Bx + bbx;
            By = By + bby;
        end
    else
        for i = 1:noL
            dn = Xp - (i-1)*dL;
            [bbz,bbx,bby] = MOT_multilayer_loop(z,Zp,dn,n,dr,R,I0);
            Bz = Bz + bbz;
            Bx = Bx + bbx;
            By = By + bby;
        end
    end
end