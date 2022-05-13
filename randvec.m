% This function is to generate random unit vector

function [x,y,z] = randvec()
    u=rand;
    v=rand;
    theta = acos(2*u-1);
    phi = 2*pi*v;
    x = sin(theta)*cos(phi);
    y = sin(theta)*sin(phi);
    z = cos(theta);
end