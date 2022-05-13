%=========================================================
% Name: B_ideal.m
% What: Calculate the ideal field profile 
% Required: B_ideal.m / section.m / MOT_coil.m
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
% Note: The function outputs are 'B_ideal.csv' / 'Z_ideal.csv'
%=========================================================


function B_ideal (Detuning,eta,v_c,v_0)
    %% constant
    u_Bohr_eff = 9.274E-24; %J/T effective Bohr magneton
    Gamma = 2*pi*5872400; % rad/s Natural linewidth
    h_bar = 6.626E-34/(2*pi); % redunced Planck constant
    a_max = 1.82E6; % m/s^2 maximum decceleration
    Lambda = 670.977338E-9; %m Resonance wavelength
    %% calculate parameters
    k = 2*pi/Lambda; % rad/m wave number
    Delta_0 = Detuning*Gamma; % rad/s laser detuning

    %% Calculate the ideal field B0 for reference and comparison
    l0 = (v_0^2-v_c^2)/(2*eta*a_max); %length of Slower
    b0 = h_bar*k*v_0/u_Bohr_eff;
    B_bias = h_bar*Delta_0/u_Bohr_eff;
    
    Loop = 1000; % Loop for calculating ideal field
    Z0 = linspace(0,l0,Loop); % Simulate the magnatic field from Z=0 m to Z = 0.4 m
    %Z0 = linspace(0,l0+0.02,Loop); % Simulate the magnatic field from Z=0 m to Z = 0.4 m
    B0 = zeros(1,Loop); % T ideal field
    for i = 1:Loop
        B0(i) = B_bias+b0*sqrt(1-Z0(i)/l0);% T magnetic field
    end
     csvwrite('B_ideal.csv',B0)
     csvwrite('Z_ideal.csv',Z0)

end