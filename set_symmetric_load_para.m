%-----------------------------------------------------------------------
% set_symmetric_load_para.m
%-----------------------------------------------------------------------
% Conversion of asymmetric to symmetric load parameter 
%-----------------------------------------------------------------------


%% Load Parameters
disp('-> convert asymmetric to symmetric loads (15 loads)')

% Load at node 611Z: set load and related voltage symmetrically PLUS changing line configuration
para_grid.n611YgZ.V_LN_rms  = [ 1 1 1 ] *  para_grid.n611Z.V_LN_rms;              % [V]
para_grid.n611YgZ.P         = [ 1/3 1/3 1/3 ] * para_grid.n611Z.P;                % [W]
para_grid.n611YgZ.QL        = [ 1/3 1/3 1/3 ] * para_grid.n611Z.QL;               % [W]
para_grid.n611YgZ.QC        = [ 1/3 1/3 1/3 ] * para_grid.n611Z.QC;               % [W]


% Load at node 611L: set load and related voltage symmetrically PLUS changing line configuration
para_grid.n611YgL.V_LN_rms  = [ 1 1 1 ] *  para_grid.n611L.V_LN_rms;              % [V]
para_grid.n611YgL.P         = [ 1/3 1/3 1/3 ] * para_grid.n611L.P;                % [W]
para_grid.n611YgL.QL        = [ 1/3 1/3 1/3 ] * para_grid.n611L.QL;               % [W]
para_grid.n611YgL.QC        = [ 1/3 1/3 1/3 ] * para_grid.n611L.QC;               % [W]


% Load at node 632 -> set load and related voltage symmetrically
para_grid.n632YgPQ.V_LN_rms =  [ 1/3 1/3 1/3 ] * sum(para_grid.n632YgPQ.V_LN_rms);  % [V]
para_grid.n632YgPQ.P        =  [ 1/3 1/3 1/3 ] * sum(para_grid.n632YgPQ.P);         % [W]
para_grid.n632YgPQ.QL       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n632YgPQ.QL);        % [W]
para_grid.n632YgPQ.QC       =  [0 0 0];                                             % [W]


% Switchable Load at low voltage (LV) node 634 -> set load and related voltage symmetrically
para_grid.n634YgPQ.V_LN_rms  = [ 1/3 1/3 1/3 ] * sum(para_grid.n634YgPQ.V_LN_rms);  % [V] 
para_grid.n634YgPQ.P         = [ 1/3 1/3 1/3 ] * sum(para_grid.n634YgPQ.P);         % [W]
para_grid.n634YgPQ.QL        = [ 1/3 1/3 1/3 ] * sum(para_grid.n634YgPQ.QL);        % [W]
para_grid.n634YgPQ.QC        = [0 0 0];                                             % [W]

% Load at node 645: set load and related voltage symmetrically PLUS changing line configuration
para_grid.n645YgPQ.V_LN_rms   = [ 1 1 1 ] * para_grid.n645PQ.V_LN_rms;               % [V]
para_grid.n645YgPQ.P          = [ 1/3 1/3 1/3 ] * para_grid.n645PQ.P;                % [W]
para_grid.n645YgPQ.QL         = [ 1/3 1/3 1/3 ] * para_grid.n645PQ.QL;               % [W]
para_grid.n645YgPQ.QC         = [0 0 0];                                             % [W]


% Load at node 646: set load and related voltage symmetrically PLUS changing line configuration
para_grid.n646YgPQ.V_LN_rms   =   [ 1 1 1 ] * para_grid.n646Z.V_LL_rms/sqrt(3);        % [V]
para_grid.n646YgPQ.P          =   [ 1/3 1/3 1/3 ] * para_grid.n646Z.P;                 % [W]
para_grid.n646YgPQ.QL         =   [ 1/3 1/3 1/3 ] * para_grid.n646Z.QL;                % [W] 
para_grid.n646YgPQ.QC         =   [0 0 0];                                             % [W]

  
% Load at node 652: set load and related voltage symmetrically PLUS changing line configuration
para_grid.n652YgZ.V_LN_rms    = [ 1 1 1 ]       *  para_grid.n652Z.V_LN_rms;      % [V]
para_grid.n652YgZ.P           = [ 1/3 1/3 1/3 ] *  para_grid.n652Z.P;             % [W]
para_grid.n652YgZ.QL          = [ 1/3 1/3 1/3 ] *  para_grid.n652Z.QL;            % [W]
para_grid.n652YgZ.QC          = [ 1/3 1/3 1/3 ] *  para_grid.n652Z.QC;            % [W]


% Load at node 671 -> set load and related voltage symmetrically 
para_grid.n671YgPQ.V_LN_rms =  [ 1/3 1/3 1/3 ] * sum(para_grid.n671YgPQ.V_LN_rms);
para_grid.n671YgPQ.P        =  [ 1/3 1/3 1/3 ] * sum(para_grid.n671YgPQ.P); 
para_grid.n671YgPQ.QL       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n671YgPQ.QL);
para_grid.n671YgPQ.QC       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n671YgPQ.QC);


% Load at node 671 -> set load and related voltage symmetrically 
para_grid.n671DPQ.V_LN_rms =   [ 1/3 1/3 1/3 ] * sum(para_grid.n671DPQ.V_LN_rms);
para_grid.n671DPQ.P        =   [ 1/3 1/3 1/3 ] * sum(para_grid.n671DPQ.P);
para_grid.n671DPQ.QL       =   [ 1/3 1/3 1/3 ] * sum(para_grid.n671DPQ.QL);
para_grid.n671DPQ.QC       =   [ 1/3 1/3 1/3 ] * sum(para_grid.n671DPQ.QC);


% Switchable Load at low voltage (LV) node 672 -> set load and related voltage symmetrically 
para_grid.n672YgPQ.V_LN_rms =  [ 1/3 1/3 1/3 ] * sum(para_grid.n672YgPQ.V_LN_rms);
para_grid.n672YgPQ.P        =  [ 1/3 1/3 1/3 ] * sum(para_grid.n672YgPQ.P);
para_grid.n672YgPQ.QL       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n672YgPQ.QL);
para_grid.n672YgPQ.QC       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n672YgPQ.QC);


% Load at node 675 -> set load and related voltage symmetrically 
para_grid.n675YgPQa.V_LN_rms =  [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQa.V_LN_rms);
para_grid.n675YgPQa.P        =  [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQa.P);
para_grid.n675YgPQa.QL       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQa.QL);
para_grid.n675YgPQa.QC       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQa.QC);



% % Switchable Load at low voltage (LV) node 675 is already symmetric 
para_grid.n675YgPQb.V_LN_rms = [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQb.V_LN_rms);
para_grid.n675YgPQb.P        = [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQb.P);
para_grid.n675YgPQb.QL       = [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQb.QL);
para_grid.n675YgPQb.QC       = [ 1/3 1/3 1/3 ] * sum(para_grid.n675YgPQb.QC);


% Switchable Load at low voltage (LV) node 680 is already symmetric 
para_grid.n680YgPQ.V_LN_rms = [ 1/3 1/3 1/3 ] * sum(para_grid.n680YgPQ.V_LN_rms);
para_grid.n680YgPQ.P        = [ 1/3 1/3 1/3 ] * sum(para_grid.n680YgPQ.P);
para_grid.n680YgPQ.QL       = [ 1/3 1/3 1/3 ] * sum(para_grid.n680YgPQ.QL);
para_grid.n680YgPQ.QC       = [ 1/3 1/3 1/3 ] * sum(para_grid.n680YgPQ.QC);


% Load at node 692: set load and related voltage symmetrically 
para_grid.n692DL.V_LN_rms =  [ 1/3 1/3 1/3 ] * sum(para_grid.n692DL.V_LN_rms); 
para_grid.n692DL.P        =  [ 1/3 1/3 1/3 ] * sum(para_grid.n692DL.P);
para_grid.n692DL.QL       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n692DL.QL);
para_grid.n692DL.QC       =  [ 1/3 1/3 1/3 ] * sum(para_grid.n692DL.QC);



