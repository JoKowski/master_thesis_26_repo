%-----------------------------------------------------------------------
% set_asymmetric_load_para.m
%-----------------------------------------------------------------------
% Set the asymmetric load parameter
%-----------------------------------------------------------------------


%% Load Parameters
disp('-> set loads as asymmetric loads (15 loads)')
         
% Load at node 611
para_grid.n611Z.V_LN_rms   = 4160/sqrt(3);
para_grid.n611Z.P          = 0;
para_grid.n611Z.QL         = 0;
para_grid.n611Z.QC         = 100e3;

% Load at node 611
para_grid.n611L.V_LN_rms   =  0.986941*4160/sqrt(3);
para_grid.n611L.P          =  170e3;
para_grid.n611L.QL         =  80e3;
para_grid.n611L.QC         =  0;

% Load at node 632
para_grid.n632YgPQ.V_LN_rms =  [ 1.021 1.042 1.0174]*2401.78;  %[V]
para_grid.n632YgPQ.P        =  [17 66  117]*1e3/2;             %[W]
para_grid.n632YgPQ.QL       =  [10 38  68]*1e3/2;              %[W]
para_grid.n632YgPQ.QC       =  [0 0 0];                        %[W]

% Switchable Load at low voltage (LV) node 634
para_grid.n634YgPQ.V_LN_rms  = [ 0.993805 1.0218 0.995862]*277.128;  %[V] 
para_grid.n634YgPQ.P         = [160 120 120]*1e3;  %[W]
para_grid.n634YgPQ.QL        = [110 90 90]*10e3;   %[W]
para_grid.n634YgPQ.QC        = [0 0 0];            %[W]


% Load at node 645
para_grid.n645PQ.V_LN_rms   = 1.03289*4160/sqrt(3);
para_grid.n645PQ.P          = 170e3; 
para_grid.n645PQ.QL         = 125e3; 
para_grid.n645PQ.QC         = 0;


% Load at node 646
para_grid.n646Z.V_LL_rms   =  4160;
para_grid.n646Z.P          =  230e3;
para_grid.n646Z.QL         =  132e3;
para_grid.n646Z.QC         =  0;
 

% Load at node 652
para_grid.n652Z.V_LN_rms =  4160/sqrt(3);
para_grid.n652Z.P        =  128e3;
para_grid.n652Z.QL       =  86e3;
para_grid.n652Z.QC       =  0;


% Load at node 671
para_grid.n671YgPQ.V_LN_rms =  [ 0.987076 1.0519 0.978043]*4160/sqrt(3);
para_grid.n671YgPQ.P        =  [17 66  117] * 1e3/2;
para_grid.n671YgPQ.QL       =  [10 38  68]  * 1e3/2;
para_grid.n671YgPQ.QC       =  [0 0 0];
% 

% Load at node 671
para_grid.n671DPQ.V_LN_rms =   [ 1.73773 1.77224 1.71484]*4160/sqrt(3);
para_grid.n671DPQ.P        =   [385 385 385 ]*1e3;
para_grid.n671DPQ.QL       =   [220 220 220]*1e3;
para_grid.n671DPQ.QC       =   [0  0  0];


% Switchable Load at low voltage (LV) node 672
para_grid.n672YgPQ.V_LN_rms =  [ 0.993805 1.0218 0.995862]*480/sqrt(3);
para_grid.n672YgPQ.P        =  [160 120 120]*1e3;
para_grid.n672YgPQ.QL       =  [110 90 90]*10e3;
para_grid.n672YgPQ.QC       =  [0 0 0];
 

% Load at node 675
para_grid.n675YgZ.V_LN_rms =  [1  1  1]*4160/sqrt(3);
para_grid.n675YgZ.P        =  [0 0 0];
para_grid.n675YgZ.QL       =  [0 0 0];
para_grid.n675YgZ.QC       =  [200 200 200]*1e3;


% Load at node 675
para_grid.n675YgPQa.V_LN_rms =  [ 0.98063 1.05431 0.975983]*4160/sqrt(3);
para_grid.n675YgPQa.P        =  [485  68  290]*1e3;
para_grid.n675YgPQa.QL       =  [190  60  212]*1e3;
para_grid.n675YgPQa.QC       =  [0 0 0];


% % Switchable Load at low voltage (LV) node 675
para_grid.n675YgPQb.V_LN_rms = [ 0.993805 1.0218 0.995862]*480/sqrt(3);
para_grid.n675YgPQb.P        = [160 120 120]*1e3;
para_grid.n675YgPQb.QL       = [110 90 90]*10e3;
para_grid.n675YgPQb.QC       = [0 0 0];


% Switchable Load at low voltage (LV) node 680
para_grid.n680YgPQ.V_LN_rms = [ 0.993805 1.0218 0.995862]*480/sqrt(3);
para_grid.n680YgPQ.P        = [120 120 120]*5e2;
para_grid.n680YgPQ.QL       = [0 0 0];
para_grid.n680YgPQ.QC       = [0 0 0];
% 

% Load at node 692
para_grid.n692DL.V_LN_rms =  [ 1.73488 1.75203 1.72341]*4160/sqrt(3);
para_grid.n692DL.P        =  [1e-6 1e-6 170 ]*1e3;
para_grid.n692DL.QL       =  [0 0 151]*1e3;
para_grid.n692DL.QC       =  [0  0  0];










