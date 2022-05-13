%=========================================================
% Name: single_loop.m
% What: compute B field for ZS single loop
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The unit of input / output [T,T,T] = single_loop(m,m,m,A)
%=========================================================
function [Bz,Bx,By] = single_loop(z,z_p,R,I0)

% Pts of interst is on z axis, with coil center at Z axis
x_p = 0;
y_p = 0;
x = 0;
y = 0;

global u0 %permeability of free space is a global variable

rc=((x-x_p).^2+(y-y_p).^2).^.5; %Radial component is required for cylindrical coordinate system.
m=(4.*R.*rc).*(((rc+R).^2)+((z-z_p).^2)).^(-1); %This is a parameter for calculating the Elliptical integrals
[kofkc,eofkc] = ellipke(m);
%Note for improved accuracy, Matlab has built in elliptical integral
%calculation but these expressions here are still very accurate when rc < R

Brc=(u0.*I0./(2.*pi.*rc)).*(z-z_p).*((((rc+R).^2)+((z-z_p).^2)).^(-.5)).*(-kofkc+eofkc.*((rc.^2+R.^2+(z-z_p).^2)./(((rc-R).^2)+((z-z_p).^2)))); %radial component of B%
Bz=(u0.*I0./(2.*pi)).*((((rc+R).^2)+((z-z_p).^2)).^(-.5)).*(kofkc-eofkc.*((rc.^2-R.^2+(z-z_p).^2)./(((rc-R).^2)+((z-z_p).^2)))); %axial component of B
Bx=Brc.*(x-x_p)./rc; %This converts the polar component into cartesian form.
By=Brc.*(y-y_p)./rc;

%The following sets any terms that result in Inf to zero, this occurs at
%the points near the coil itself.
Bx(isnan(Bx)) = 0 ;
By(isnan(By)) = 0 ;
Bz(isnan(Bz)) = 0 ;






   


   
