%=========================================================
% Name: trajectory.m
% What: 1. Simulate the trajectory of atoms with given initial velocity
%          with given B-profile and laser beam properties.
%       2. Calculate final velocity of ZS, and find the caputure 
%          velocity for MOT.
% Required: B_coil.csv / Z_coil.csv / R_scatt_s0.m / Vrec2.m
% Note: 1. Simulate using ideal field(flag = 0) or coil field(flag = 1)
%       2. This function is to know the trajectory given B field and initial
%          velocity, and output the V-Z data 
%       3. SI unit for input 
% Author: Yu-Hao Yeh (yuhyeh@iu.edu)
% Date: 2022/05/12
%=========================================================
% The input data is (test_velocity, velocity to generate ideal profile,
% laser detuning(Gamma), safety parameter, capture_velocity of MOT beam, 
% flag to determine use ideal field(0), or Coil field(1))

function trajectory(v_i,v_0,Detuning,eta,v_c,flag,nI,MOT_Beam_radius)
    if flag == 0
        Z0 = csvread('Z_ideal.csv');
        B0 = csvread('B_ideal.csv');
    else
        Z0 = csvread('Z_coil.csv');
        B0 = csvread('B_coil.csv');
    end

    %% constant
    u_Bohr_eff = 9.274E-24; %J/T effective Bohr magneton
    Gamma = 2*pi*5872400; % rad/s Natural linewidth
    h_bar = 6.626E-34/(2*pi); % redunced Planck constant
    a_max = 1.82E6; % m/s^2 maximum decceleration
    m_Li = 9.988E-27; % kg Lithium mass
    I_sat = 2.54; % mW/cm^2 for D2 line
    Lambda = 670.977338E-9; %m Resonance wavelength
    %% calculate parameters
    chr = int2str(v_i); % convert v_i to string for csv file
    %v_0 = sqrt(k_B*T_oven/m_Li); % m/s initial velocity
    %v_i = v_0; % test initital velocity
    k = 2*pi/Lambda; % rad/m wave number
    Delta_0 = Detuning*Gamma; % rad/s laser detuning
    %I = 2*P/(pi*Beam_radius^2); %mW/cm^2 laser intensity ( defined by safety parameter)
    I = nI*eta/(1-eta)*I_sat;% 
    s0 = I/I_sat;
    %% loop number
    t_total = (v_i-v_c)/(a_max*eta); % s
    dt = 1E-9; % s constant time step

    n = 100*round(t_total/dt);% number of loops ("70" is just a number to 
    % increase the loop and to compensate the inefficiency of slowering for 
    % atoms with initial velocity below v_0. In short, just for declaring a 
    % large enough memory for the loop to run) 

    N = 1:n; % Label scattering events from 1 to n
    
    
    l0 = (v_0^2-v_c^2)/(2*eta*a_max); %length of Slower
    
    %% Declare variables for iteration
    a = zeros(1,n); %deceleration of each event
    R = zeros(1,n); %scattering rate for each event
    Del = zeros(1,n); %scattering rate for each event
    B = zeros(1,n); %B field by interpolation
    
    VVz = zeros(1,n+1);% m/s velocity in z direction 
    VVx = zeros(1,n+1);% m/s velocity in x direction 
    VVy = zeros(1,n+1);% m/s velocity in y direction 
    SSz = zeros(1,n+1);% m position in z direction
    SSx = zeros(1,n+1);% m position in x direction
    SSy = zeros(1,n+1);% m position in y direction
    t_tot = zeros(1,n+1); % s accumulative time for scattering
   
    %% Initial condition
    VVz(1) = v_i; % m/s velocity in Z direction %
    i=1;
    %% Calculate all parameters2
    while (i<=n && SSz(i)<=l0-MOT_Beam_radius)
            B(i) = interp1(Z0,B0,SSz(i),'linear'); %1) T determine B field by interpolation
            [R(i),Del(i)] = R_scatt_s0(B(i),VVz(i),s0,Delta_0); %(1/s) call R_scatt to calculate the scattering rate of each event
            a(i) = (h_bar*k/m_Li)*R(i); % m/s^2 decceleration #####
            %% Absorption
            VVz(i+1) = VVz(i)-a(i)*dt; %Absorption 
            %% Distance and position btw absorbtion and emission.
            S_mid_absorb = VVz(i+1)*dt/2;
            SSz(i+1) = SSz(i)+S_mid_absorb;
    
            %% Spontaneous emission
            [vx,vy,vz] = Vrec2(a(i)*dt); % turn on the emission %sqrt(vx^2+vy^2+vz^2);% = 0.0989;
            VVz(i+1) = VVz(i+1)+vz; % m/s ~~(= v-v_rec), velocity in z direction #####v
    
            %% Distance and postion btw emission and next absoption
            S_final_emit = VVz(i+1)*dt/2;
            SSz(i+1) = SSz(i+1)+S_final_emit;% position for next absorption event
    
            VVx(i+1) = VVx(i)+vx; % m/s  velocity in z direction #####v
            SSx(i+1) = SSx(i)+VVx(i)*dt;% trajectory in x direction #####v
    
            VVy(i+1) = VVy(i)+vy; % m/s  velocity in y direction #####v
            SSy(i+1) = SSy(i)+VVy(i)*dt; % trajectory in y direction #####v
    
            %% End test 
            t_tot(i+1) = t_tot(i)+dt; % s accumulative time for scattering
            %% Calculate the test
            i=i+1;
    end
    
    end_num = nnz(SSz)+1;
    VVz = VVz(1:end_num);
    SSz = SSz(1:end_num);
    VVx = VVx(1:end_num);
    SSx = SSx(1:end_num);
    VVy = VVy(1:end_num);
    SSy = SSy(1:end_num);
    t_tot = t_tot(1:end_num);
    a = a(1:end_num);
    N = N(1:end_num);
    R = R(1:end_num);
    Del = Del(1:end_num);
    B = B(1:end_num);
    
    Vfilename = strcat('V_actual_',chr,'.csv');
    Zfilename = strcat('Z_actual_',chr,'.csv');
    
    csvwrite(Vfilename,VVz)
    csvwrite(Zfilename,SSz)
end