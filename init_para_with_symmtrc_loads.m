%---------------------------------------------------------------------
% init_para_with_symmtrc_loads.m
%---------------------------------------------------------------------
% 13 Feedertest Grid with balanced loads 
%---------------------------------------------------------------------
clear all


%% Simulation parameter
Ts = 50e-6;  % 50  mycro sec
Sn = 15e6;
H = 8;
Rdroop = 0.06;
Kdroop = Sn/(Rdroop*60);

%% Parameter for event based grid analysis and later for control 
%TestCase = 1;  % Power flow from External Grid into the 13 Feeder Grid 
                % P = +3MW
%TestCase = 2;   % Power flow from the 13 Feeder Grid into External Grid 
                % P = -6MW
TestCase = 3;   % Balanced 60Hertz

        

if TestCase == 1
    % Loads/Demands
    disp('Set Load Configuration for TestCase1')
    para_grid.n634YgPQ.ts = 0.1;    % switching on at ts
    para_grid.n672YgPQ.ts = 0.15;   % switching on at ts
    para_grid.n675YgPQ.ts = 0.2;    % switching on at ts
    para_grid.n680YgPQ.ts = 0.25;   % switching on at ts

    disp('Set Renewable Resources for TestCase1')
    % Renewable Resources (for all the same)
    para_PVgeno.S_ini1   = 0;
    para_PVgeno.S_ini2   = 0;
    para_PVgeno.S_ini3   = 0;
    para_PVgeno.S_ini4   = 0;
    para_PVgeno.S_step  = 0;
    para_PVgeno.t_step  = 0.2;
end

if TestCase == 2
    % Loads/Demands
    disp('Set Load Configuration for TestCase2')
    para_grid.n634YgPQ.ts = 100.0;   % switching on at ts
    para_grid.n672YgPQ.ts = 100.0;   % switching on at ts
    para_grid.n675YgPQ.ts = 100.0;   % switching on at ts
    para_grid.n680YgPQ.ts = 100.0;   % switching on at ts

    disp('Set Renewable Resources for TestCase2')
    % Renewable Resources (for all the same)
    para_PVgeno.S_ini1   = 400;
    para_PVgeno.S_ini2   = 400;
    para_PVgeno.S_ini3   = 400;
    para_PVgeno.S_ini4   = 400;
    para_PVgeno.S_step  = 400;
    para_PVgeno.t_step  = 0.2;
end

if TestCase == 3
    % Loads/Demands
    disp('Set Load Configuration for TestCase3')
    para_grid.n634YgPQ.ts = 0;   % switching on at ts
    para_grid.n672YgPQ.ts = 0;   % switching on at ts
    para_grid.n675YgPQ.ts = 0;   % switching on at ts
    para_grid.n680YgPQ.ts = 0;   % switching on at ts

    disp('Set Renewable Resources for TestCase3')
    % Renewable Resources (for all the same)
    para_PVgeno.S_ini1   = 186;
    para_PVgeno.S_ini2   = 186;
    para_PVgeno.S_ini3   = 186;
    para_PVgeno.S_ini4   = 186;
    para_PVgeno.S_step  = 0;
    para_PVgeno.t_step  = 0;

    para_PVgeno.S_step11  = 0;
    para_PVgeno.t_step11  = 0;
    para_PVgeno.S_step12  = 0;
    para_PVgeno.t_step12  = 0;

    para_PVgeno.S_step21  = 0;
    para_PVgeno.t_step21  = 0;
    para_PVgeno.S_step22  = 0;
    para_PVgeno.t_step22  = 0;

    para_PVgeno.S_step31  = 0;
    para_PVgeno.t_step31  = 0;
    para_PVgeno.S_step32  = 0;
    para_PVgeno.t_step32  = 0;

    para_PVgeno.S_step41  = 0;
    para_PVgeno.t_step41  = 0;
    para_PVgeno.S_step42  = 0;
    para_PVgeno.t_step42  = 0;

end



%% Parameter of the Dynamic Equivalent External Grid (DEG) connected to the 13 Feeder Grid
disp('Set parameter of the Dynamic Equivalent External Grid') 
fg          = 60;             % nominal frequency
U_rms       = 4160;           % nominal voltage of the grid in RMS ph-ph
kv_Vp       = [ 1.021 1.042 1.0174]; % asymmetric phase amplitde factor 
S_sc        = 100e6;                 % short circuit apparent power    
phi_vec_deg = [ -2.49 -121.72 117.83];  % phi values of the infinite bus in the 13 feeder grid
phi_vec     = phi_vec_deg/180*pi; 

% Grid Dynamic Equivalent Line parameters 
XoverR = 8;              % X/R ratio
I_rms  = S_sc / (sqrt(3) * U_rms);
Rg = (U_rms^2 / S_sc) / (sqrt(1 + XoverR^2));
Xg = Rg * XoverR;  
Lg = Xg /(2*pi*fg);



%% Power grid parameter
disp('Set parameter of the 13NodeTestFeeder Grid') 
%para_grid_13Node_Feeder_with_Asymmetric_Load_Para

set_line_parameter
set_asymmetric_load_para
set_symmetric_load_para    % transtion to symmetric load parameter



% have to be defined 
S_base_all = 1;


%% Renewable source integration to grid
% PV array one stage DC/AC VSC interface
% 
% Created: 15.11.2024
% Version: VSC AC Side connected to ideal grid
% Control structure: Inner Current Control Loop
%                    Phase-Locked Loop
%                    I/V abc to qd0 transformation
%                    id (reactive power reference)
%                    iq (active power reference)
%
% Mod: 27.11.2024
% Version: PV Array + DC/AC VSC connected to ideal grid
% Control structure: MPPT (open circuit method)
%                    Outer DC Voltage Loop
%                    Outer Reactive Power Loop
%                    Inner Current Control Loop
%                    Phase-Locked Loop
%                    I/V abc to qd0 transformation
%                    id (reactive power reference)
%                    iq (active power reference)
%
% Mod: 03.12.2024
% Included P-V curves of PV array
%
% Author: Ullón Huber, Mauricio

%% PV module (Model: KC200GT Solar)

disp('Set parameter of PV Generator System') 

Imp_n=7.61;                 % Nominal max. power current
Vmp_n=26.3;                 % Nominal max. power voltage
Pmax_e=200.143;             % Nominal max. power
Isc_n=8.21;                 % Nominal short-circuit current
Voc_n=32.9;                 % Nominal open-circuit voltage
KV=-0.123;                  % Temp. coeff. of Voc
KVmp=Vmp_n/Voc_n;           % Ratio max. power to open-circuit (MPPT)
KI=0.0032;                  % Temp. coeff. of Isc
Ns_cell=54;                 % Number of cells connected in series
a=1.3;                      % Diode ideality constant
Rp=425.405;                 % Eq. parallel resistance
Rs=0.221;                   % Eq. series resistance

% Size of the PV Array
Ns=45;                      % Number of modules in series
Np=345;                     % Number of modules in parallel

% Standard test conditions (STC)
Gn=1000;                    % Irradiance
Tn=25;                      % Temperature

% Aditional values for calculations
k=1.3806503e-23;            % Bolzmann's constant
q=1.6021764e-19;            % Charge of electron

% Nominal Ipv considering Rp and Rs
Ipv_n=Isc_n*(Rp+Rs)/Rp;

% Plot P-V curves
%get_PV_curves;

disp('Set parameter of the VSC controllers') 

%% Grid parameters
Ug=480;                     % VSC voltage (L2L-rms) , Terminal voltage
Vg=Ug/sqrt(3);              % L2L-rms to L2N-rms
%fg=60;                      % Nominal frequency -> already defined by grid parameter
Vp=Ug*sqrt(2)/sqrt(3);      % Peak voltage (L2G)

%% VSC parameters
Scn=3e6;                    % Nominal converter power
Ucn=480;                    % Nominal converter voltage (L2L-rms)
Vdc=1200;                   % Nominal DC voltage
Cdc=278e-3;                 % Capacitance of VSC DC bus
Rc=0.768e-3;                % Resistance of coupling filter
Lc=20.372e-6;               % Inductance of coupling filter
ma=1;                       % Modulation index

% Converter limits
Ic_max=Scn/(3*Vg);              % Current
Vc_max=(1/sqrt(2))*ma*(Vdc/2);  % Voltage (L2N-rms)
Uc_max=Vc_max*sqrt(3);          % Voltage (L2L-rms)

fprintf('VSC output limits\n');
fprintf('Current: Ic_max = %.2f A\n',Ic_max);
fprintf('Voltage: Uc_max = %.2f V\n',Uc_max);

%% Grid impedance
% XLratio=10;                 % Ratio between reactance and resistance X/R
% SCR=5;                      % Short-circuit ratio
% Sscg=SCR*Scn;               % Short-circuit power of the grid
% Zsc=Ug^2/Sscg;              % Short-circuit impedance
% Rg=Zsc/sqrt(1+XLratio^2);   % Resistance of the grid
% Xg=XLratio*Rg;              % Reactance of the grid
% Lg=Xg/(2*pi*fg);            % Inductance of the grid
              
%% VSC Control parameters



% PLL: All PLLs also the PLL for neg. and pos. seq. calculation=
xi_PLL=sqrt(2)/2;
omega_PLL=2*pi*fg;
Kp_PLL=xi_PLL*2*omega_PLL/Vp;
tau_PLL=2*xi_PLL/omega_PLL;
Ki_PLL=Kp_PLL/tau_PLL;

% Current control loop
tau_C=1e-3;
Kp_C=Lc/tau_C;
Ki_C=Rc/tau_C;

% Power control loop
Ki_pq=3*Vp/2;
tau_P=10*tau_C;
Kp_P=tau_C/(Ki_pq*tau_P);
Ki_P=1/(Ki_pq*tau_P);

% DC voltage Control
xi_DC=sqrt(2)/2;
omega_DC=4/(0.1*xi_DC);
Kp_DC=2*xi_DC*omega_DC*Cdc;
Ki_DC=omega_DC^2*Cdc;


%% DVPP Frequency Control
load estimated_tfs_pt1.mat

% Testing Parameters, delta_P Step
p_jump_time = 5;
p_jump_value = 2*-1e6;

ppv_nom = 4.9687e+05; % gemessen

% 1:1 aus Paper
b0 = -1.48911;
a0 = 0.01778;
a1 = 0.2667;
Tagg_pf = tf([b0],[1 a1 a0]);
m1paper = tf([1],[1 2]);
m2paper = tf([1],[1 3]);
m3paper = tf([1 3 1],[1 5 6]);

% % Eigene m Faktoren
% m1 = tf([1],[1 2]);
% m2 = tf([1 ],[1 3]);
% m3 = tf([1],[1 4]);
% m4 = tf([1 6 8 -2],[1 9 26 24]);

% Andere Variante m Faktoren
m1 = tf([1],[0.5 1]);
m2 = tf([0.5 0],[0.5 1]);
m3 = m1*m2;
m4 = 1-m1-m2-m3;

Gplant_real = tf_ges_real;
Gplant_simp = tf_ges_simp;


% PI Regler Entwurf
Kdcgain = dcgain(Gplant_real);
kp_dvpp = -(2*Kdcgain)^-1;
Ti_dvpp = 2;
ki_dvpp = kp_dvpp/Ti_dvpp;
s = tf('s');
C_agg = (kp_dvpp*s+ki_dvpp)/s;

T_pf_agg = C_agg;

