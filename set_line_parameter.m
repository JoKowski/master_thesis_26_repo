%-----------------------------------------------------------------------
% set_line_parameter.m
%-----------------------------------------------------------------------
%   Set grid line parameter for the IEEE 13 Node Test Feeder grid
%
%   1.) Definition of distributed Line Parameters defined in the orginal 
%       13-Feeder model power_13NodeTestFeeder.slx  
%   2.) Based on the definition in 1, a common R and L value (for each line 
%       the same) is calculated for all 9 line segments 
%-----------------------------------------------------------------------

disp('-> set line parameter')

% miles/km
mi2km = 1.60934;

% feet to km
ft2km = 0.0003048;

% microsiemens to Farads
ms2F = 1/2/pi/60*1e-6;


%%---------------------------------------------------------------------
%  Distributed Line Parameters defined in the orginal 13-Feeder model
%  power_13NodeTestFeeder.slx
%---------------------------------------------------------------------

% Configuration 601 - series reactance - ohm/mile
R_601 = [0.3465 0.1560 0.1580;0.1560 0.3375 0.1535;0.1580 0.1535 0.3414];
X_601 = [1.0179 0.5017 0.4236;0.5017 1.0478 0.3849;0.4236 0.3849 1.0348];

% charging susceptance - microsiemens/mile

B_601 = [6.2998 -1.9958 -1.2595;-1.9958 5.9597 -0.7417;-1.2595 -0.7417 5.6386];

% convert for SPS

R_601 = R_601/mi2km;

L_601 = X_601/mi2km/2/pi/60;

C_601 = B_601/mi2km*ms2F;


% Configuration 602 - series reactance - ohm/mile

R_602 = [0.7526 0.1580 0.1560;0.1580 0.7475 0.1535;0.1560 0.1535 0.7436];
X_602 = [1.1814 0.4236 0.5017;0.4236 1.1983 0.3849;0.5017 0.3849 1.2112];

% charging susceptance - microsiemens/mile

B_602 = [5.6990 -1.0817 -1.6905;-1.0817 5.1795 -0.6588;-1.6905 -0.6588 5.4246];

% convert for SPS

R_602 = R_602/mi2km;

L_602 = X_602/mi2km/2/pi/60;

C_602 = B_602/mi2km*ms2F;


% Configuration 603 - series reactance - ohm/mile

R_603 = [0 0 0;0 1.3294 0.2066;0 0.2066 1.3238];
X_603 = [0 0 0;0 1.3471 0.4591;0 0.4591 1.3569];

% charging susceptance - microsiemens/mile

B_603 = [0 0 0;0 4.7097 -0.8999;0 -0.8999 4.6658];

% convert for SPS

R_603 = R_603/mi2km;

L_603 = X_603/mi2km/2/pi/60;

C_603 = B_603/mi2km*ms2F;


% Configuration 604 - series reactance - ohm/mile

R_604 = [1.3238 0 0.2066;0 0 0;0.2066 0 1.3294];
X_604 = [1.3569 0 0.4591;0 0 0;0.4591 0 1.3471];

% charging susceptance - microsiemens/mile

B_604 = [4.6658 0 -0.8999;0 0 0;-0.8999 0 4.7097];


% convert for SPS

R_604 = R_604/mi2km;

L_604 = X_604/mi2km/2/pi/60;

C_604 = B_604/mi2km*ms2F;


% Configuration 605 - series reactance - ohm/mile

R_605 = [0 0 0;0 0 0;0 0 1.3292];
X_605 = [0 0 0;0 0 0;0 0 1.3475];

% charging susceptance - microsiemens/mile

B_605 = [0 0 0;0 0 0;0 0 4.5193];


% convert for SPS

R_605 = R_605/mi2km;

L_605 = X_605/mi2km/2/pi/60;

C_605 = B_605/mi2km*ms2F;

% Configuration 606 - series reactance - ohm/mile

R_606 = [0.7982 0.3192 0.2849;0.3192 0.7891 0.3192;0.2849 0.3192 0.7982];
X_606 = [0.4463 0.0328 0.0143;0.0328 0.4041 0.0328;0.0143 0.0328 0.4463];
%X_606 = [0.4463 0.0328 -0.0143;0.0328 0.4041 0.0328;-0.0143 0.0328 0.4463];

% charging susceptance - microsiemens/mile

B_606 = [96.8897 -1e-6 -1e-6;-1e-6 96.8897 -1e-6;-1e-6 -1e-6 96.8897];

% convert for SPS

R_606 = R_606/mi2km;

L_606 = X_606/mi2km/2/pi/60;

C_606 = B_606/mi2km*ms2F;


% Configuration 607 - series reactance - ohm/mile

R_607 = [1.3425 0 0;0 0 0;0 0 0];
X_607 = [0.5124 0 0;0 0 0;0 0 0];

% charging susceptance - microsiemens/mile

B_607 = [88.9912 0 0;0 0 0;0 0 0];

% convert for SPS

R_607 = R_607/mi2km;

L_607 = X_607/mi2km/2/pi/60;

C_607 = B_607/mi2km*ms2F;


%%---------------------------------------------------------------------
%  Based on previous definition in 1, a common R and L value (for 
%  each line the same) is calculated for all 9 line segments 
%----------------------------------------------------------------------
set_XtoR_ON    = 0; % {0,1}  
set_XtoR_value = 10;
if set_XtoR_ON == 1
    fprintf('Set X/R to %f\n',set_XtoR_value)
end

% 1.) #603  
para_grid.l645_646.length_feet = 300;
para_grid.l645_646.length_km   = para_grid.l645_646.length_feet * ft2km;
para_grid.l645_646.R           = R_603(2,2) * para_grid.l645_646.length_km;
if set_XtoR_ON == 1
   para_grid.l645_646.L        = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l645_646.R);
else
   para_grid.l645_646.L        = L_603(2,2) * para_grid.l645_646.length_km;
end
para_grid.l645_646.XtoR        = para_grid.l645_646.L*2*pi*fg/para_grid.l645_646.R;
%para_grid.l645_646


% 2.) #603  
para_grid.l644_645.length_feet = 500;
para_grid.l644_645.length_km   = 500 * ft2km;
para_grid.l644_645.R           = R_603(2,2) * para_grid.l644_645.length_km;
if set_XtoR_ON == 1
   para_grid.l644_645.L        = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l644_645.R);
else
   para_grid.l644_645.L        = L_603(2,2) * para_grid.l644_645.length_km;
end
para_grid.l644_645.XtoR        = para_grid.l644_645.L*2*pi*fg/para_grid.l644_645.R;
%para_grid.l644_645


% 3.) #602  
para_grid.l633_644.length_feet = 500;
para_grid.l633_644.length_km   = 500 * ft2km;
para_grid.l633_644.R           = R_602(1,1) * para_grid.l633_644.length_km;
if set_XtoR_ON == 1
    para_grid.l633_644.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l633_644.R);
else
    para_grid.l633_644.L           = L_602(1,1) * para_grid.l633_644.length_km;
end
para_grid.l633_644.XtoR        = para_grid.l633_644.L*2*pi*fg/para_grid.l633_644.R;
%para_grid.l633_644


% 4.) #601
para_grid.l632_671.length_feet = 2000;
para_grid.l632_671.length_km   = para_grid.l632_671.length_feet * ft2km;
para_grid.l632_671.R           = R_601(1,1) * para_grid.l632_671.length_km;
if set_XtoR_ON == 1
    para_grid.l632_671.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l632_671.R);
else
    para_grid.l632_671.L           = L_601(1,1) * para_grid.l632_671.length_km;
end
para_grid.l632_671.XtoR        = para_grid.l632_671.L*2*pi*fg/para_grid.l632_671.R;
%para_grid.l632_671


% 5.) #607
para_grid.l652_684.length_feet = 800;
para_grid.l652_684.length_km   = 800 * ft2km;
para_grid.l652_684.R           = R_607(1,1) * para_grid.l652_684.length_km;
if set_XtoR_ON == 1
  para_grid.l652_684.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l652_684.R);
else
  para_grid.l652_684.L           = L_607(1,1) * para_grid.l652_684.length_km;
end
para_grid.l652_684.XtoR        = para_grid.l652_684.L*2*pi*fg/para_grid.l652_684.R;
%para_grid.l652_684


% 6.) #604
para_grid.l671_684.length_feet = 300;
para_grid.l671_684.length_km   = 300 * ft2km;
para_grid.l671_684.R           = R_604(1,1) * para_grid.l671_684.length_km;
if set_XtoR_ON == 1
    para_grid.l671_684.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l671_684.R);  
else
    para_grid.l671_684.L           = L_604(1,1) * para_grid.l671_684.length_km;
end
para_grid.l671_684.XtoR        = para_grid.l671_684.L*2*pi*fg/para_grid.l671_684.R;
%para_grid.l671_684


% 7.) #606
para_grid.l675_692.length_feet = 500;
para_grid.l675_692.length_km   = 500 * ft2km;
para_grid.l675_692.R           = R_606(1,1) * para_grid.l675_692.length_km;
if set_XtoR_ON == 1
    para_grid.l675_692.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l675_692.R);  
else
    para_grid.l675_692.L           = L_606(1,1) * para_grid.l675_692.length_km;
end
para_grid.l675_692.XtoR        = para_grid.l675_692.L*2*pi*fg/para_grid.l675_692.R;
%para_grid.l675_692


% 8.) #605
para_grid.l611_684.length_feet = 300;
para_grid.l611_684.length_km   = 300 * ft2km;
para_grid.l611_684.R           = R_605(3,3) * para_grid.l611_684.length_km;
if set_XtoR_ON == 1
    para_grid.l611_684.L           = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l611_684.R);  
else
    para_grid.l611_684.L           = L_605(3,3) * para_grid.l611_684.length_km;
end
para_grid.l611_684.XtoR        = para_grid.l611_684.L*2*pi*fg/para_grid.l611_684.R;
%para_grid.l611_684


% 9.) #607  
para_grid.l671_680.length_feet = 1000;
para_grid.l671_680.length_km   = 1000 * ft2km;
para_grid.l671_680.R           = R_601(3,3) * para_grid.l671_680.length_km;
if set_XtoR_ON == 1
    para_grid.l671_680.L        = 1/(2*pi*fg) * (set_XtoR_value * para_grid.l671_680.R);  
else
    para_grid.l671_680.L       = L_601(3,3) * para_grid.l671_680.length_km;
end
para_grid.l671_680.XtoR        = para_grid.l671_680.L*2*pi*fg/para_grid.l671_680.R;
%para_grid.l671_680




