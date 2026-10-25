%% ========================================================================
%  build_Yb_13feeder.m
%  ------------------------------------------------------------------------
%  Knotenadmittanzmatrix Y_b(s) und Maschinenimpedanzen Z_m(s) des Modells
%       power_13NodeTestFeeder_4VSC_PV_DEG_symmtrc_loads_
%  nach Gu et al., Gl. (14)-(17), in der D'Q'-Basis (komplexe +/- -Basis).
%
%  Gl. (16):  Y_bkl(s) = -[ R_kl+(s+j w0)L_kl ,        0            ]^-1
%                         [        0          , R_kl+(s-j w0)L_kl   ]
%
%  Gl. (17):  Y_bkk(s) = - sum_{l~=k} Y_bkl(s)
%                        + [ 1/R_kk+(s+j w0)C_kk ,          0            ]
%                          [          0          , 1/R_kk+(s-j w0)C_kk   ]
%
%  ------------------------------------------------------------------------
%  DARSTELLUNG UND DIMENSION
%    Y_b ist nach Gl.(15) eine N x N BLOCKmatrix; jeder Eintrag Y_bkl ist
%    die 2x2-Matrix aus Gl.(16). Das ist hier die primaere Darstellung:
%       Ybkl{k,l}  -  N x N cell, jede Zelle 2x2   <- so wird es notiert
%    Ausgeschrieben wird nur fuer die Numerik, weil inv/det/eig/lu auf der
%    2N x 2N-Matrix arbeiten (blockweise ginge das nur, wenn die Bloecke
%    untereinander kommutieren - sie tun es nicht):
%       Yb = blk2full(Ybkl)  -  2N x 2N           <- dieselbe Matrix
%    N ist die Zahl der KNOTEN, nicht die der Zweige und nicht die der
%    Lasten. Standardeinstellung: N = 12, also 12 x 12 Bloecke.
%
%  KNOTEN UND MASCHINEN (Standardeinstellung)
%    12 Knoten, davon 11 mit Maschine. 684 ist eine reine Verzweigung vom
%    Grad 3 und traegt keine Maschine -> Y_m684 = 0, Z_m existiert nicht.
%    Es wird NICHT reduziert. assemble_Ysys.m benutzt stattdessen die duale
%    Form Gl.(20)/(21), die Z_b statt Z_m braucht. Y_b ist regulaer, weil
%    die beiden Kondensatorbaenke als C_kk darin stehen (M1).
%    11 Zweige bei 12 Knoten = N-1, also ein Baum - der Feeder ist radial.
%
%  ------------------------------------------------------------------------
%  TOPOLOGIE (vom Anwender bestaetigt)
%   T1  644 hat keine eigene physikalische Repraesentanz -> 644 == 632.
%       l644_645 ist die Leitung 632-645, l633_644 die Leitung 632-633.
%   T2  692_DI ist ein Messblock, kein RL-Zweig -> 692 == 671.
%       l675_692 ist die Leitung 671-675.
%   T3  Die Leitungen sind reine Laengs-RL-Zweige, keine Ladekapazitaet.
%   T4  Die Trafos sind Yg/Yg -> keine 30-Grad-Drehung; Umrechnen auf
%       4,16 kV ist eine reelle Skalierung mit ue^2 = 75.111.
%   T5  Der 1-kOhm-Zweig der Quelle liegt im Nullsystem -> hier ohne Wirkung.
%   T6  Dreieckslasten sind im 2x2-Modell nicht von Yg-Lasten gleicher
%       Leistung unterscheidbar (Unterschied nur im Nullsystem).
%
%  MODELLENTSCHEIDUNGEN
%   M1  Die beiden Kondensatorbaenke 611 Yg Z und 675 Yg Z nehmen keine
%       Wirkleistung auf und sind keine Maschinen. Sie gehen als C_kk in
%       Y_bkk ein - genau der Platz, den Gl.(17) dafuer vorsieht.
%       -> Y_b enthaelt damit 13 Lasten in Z_m und 2 Baenke in Y_b.
%   M2  633 und der oberspannungsseitige Punkt von XFM-680 tragen nichts.
%       Sie werden nicht als Knoten gefuehrt; stattdessen sind Leitung und
%       Trafostreuimpedanz zu EINEM Serien-RL-Zweig zusammengefasst:
%          632-634 = l633_644 + XFM-634,   671-680 = l671_680 + XFM-680
%       Das erhaelt die Form von Gl.(16) mit genau einem R_kl und L_kl.
%   M3  Damit M2 exakt ein Serien-RL ergibt, entfallen die Magnetisierungs-
%       zweige (opt.xfmMag = false). Preis: |Y_mag|/|Y_serie| = 6,5e-5,
%       Abweichung der Treiberimpedanzen <= 5e-4 relativ (0,5..1000 Hz,
%       gegen das exakte 14-Knoten-Modell nachgerechnet).
%       Exakte Referenz: opt.mergeSeriesNodes = false und opt.xfmMag = true
%       -> 14 Knoten. Gegengerechnet wurde beides gegen diese Referenz;
%       die Abweichung des 12-Knoten-Modells bleibt unter 5e-4.
%   M4  Der LV-Knoten hinter XFM-680 heisst wie im Modell "680". Bei 675
%       tragen beide Seiten im Modell dasselbe Label, deshalb hier
%       675 (OS, mit Bank und Last) und 675L (NS, mit Last und VSC).
%   M5  Die DEG ist eine Maschine wie jede andere - ideale Quelle hinter
%       Rg + (s +- j w0) Lg - und sitzt am Knoten 632.
%
%   E1  ERWEITERUNG zu Gl.(17): der Term 1/((s +- j w0) L_kk) ist ergaenzt.
%       Gl.(17) kennt nur R_kk und C_kk; ohne die Erweiterung ginge jede
%       induktive Querkomponente verloren (relevant nur, wenn Lasten per
%       opt.loadsInYb nach Y_b gelegt werden).
%
%  ------------------------------------------------------------------------
%  ERZEUGTE VARIABLEN (symbolisch, class sym)
%    s, w0, nodeName, R_kl, L_kl, R_kk, C_kk, L_kk
%    Ybkl {NxN cell, 2x2}, Ybkk {Nx1 cell, 2x2}, Yb (2N x 2N)
%    load_  Struct je Last mit R/L/C und den 2x2-Bloecken .Ym / .Zm
%    Ym_node, Zm_node {Nx1 cell, 2x2},  Y_m (2N x 2N, blockdiagonal)
%    Yb_plus_Ym = Y_b + Y_m = Zhat^-1   (NICHT Yhat, siehe assemble_Ysys.m)
%    par .name/.sym/.val,  Yb_num, YbYm_num, Yb_fun, YbYm_fun
%
%  VORAUSSETZUNG
%    init_para_with_symmtrc_loads.m wurde ausgefuehrt (para_grid, fg, Rg, Lg).
% ========================================================================


%% ------------------------------------------------------------- 0) Optionen
opt = struct();

opt.w0_symbolic      = false;  % true: w0 bleibt Symbol
opt.mergeSeriesNodes = true;   % M2: 633 und 680_OS nicht als Knoten fuehren
opt.xfmMag           = false;  % M3: Magnetisierungszweige beruecksichtigen
opt.capBanksInYb     = true;   % M1: Kondensatorbaenke als C_kk in Y_b
opt.shuntL           = true;   % Erweiterung E1
opt.DEG              = 'machine632';  % 'machine632'|'shunt632'|'none'
opt.deltaLoadV       = 'LL';   % V_LN_rms-Feld der D-Lasten enthaelt Vab
opt.exactSym         = false;  % true -> sym(x,'f')
opt.makeFun          = true;
opt.verbose          = true;

% --- Transformatoren XFM-634 / -671 / -675 / -680 (alle gleich) ---------
%  Aus den Blockmasken. Die Maske steht auf SI-Einheiten: R1/L1 in Ohm bzw.
%  Henry, nicht in pu. Gegenprobe: R1 = 0.011 pu, X1 = 0.020 pu, Rm = Xm = 500 pu.
opt.xfm.Sn = 500e3;  opt.xfm.V1 = 4160;  opt.xfm.V2 = 480;
opt.xfm.R1 = 0.38072;    opt.xfm.L1 = 0.0018362;
opt.xfm.R2 = 0;          opt.xfm.L2 = 0;
opt.xfm.Rm = 17306;      opt.xfm.Lm = 45.905;
fprintf('-> Schritt 0: Optionen\n');
fprintf('-> Paramter aller Transformatoren aus Simulinkmodell gelesen und hartkodiert.\n');

%% --------------------------------------------- 1) Vorbedingungen pruefen
fprintf('\n=== Y_b / Z_m nach Gu et al., Gl.(14)-(17) ===\n');
fprintf('-> Schritt 1: Vorbedingungen pruefen\n');

assert(exist('para_grid','var')==1, ...
    'para_grid fehlt. Bitte zuerst init_para_with_symmtrc_loads.m ausfuehren.');
assert(exist('fg','var')==1, 'fg (Nennfrequenz) fehlt im Workspace.');
if ~strcmpi(opt.DEG,'none')
    assert(exist('Rg','var')==1 && exist('Lg','var')==1, ...
        'Rg/Lg (DEG-Kurzschlussimpedanz) fehlen im Workspace.');
end
if opt.mergeSeriesNodes && opt.xfmMag
    warning('Yb:mag', ['opt.mergeSeriesNodes = true und opt.xfmMag = true: ' ...
        'die Magnetisierung von XFM-634 und XFM-680 haette keinen Knoten. ' ...
        'Sie wird fuer diese beiden weggelassen.']);
end

s = sym('s');
if opt.w0_symbolic, w0 = sym('w0','positive');
else,               w0 = 2*sym(pi)*sym(fg);
end
w0n   = 2*pi*fg;
V_ref = 4160;

Zb_x = opt.xfm.V1^2/opt.xfm.Sn;
fprintf(['   Trafo-Gegenprobe: Zbase = %.4f Ohm | R1 = %.4f pu | X1 = %.4f pu' ...
         ' | Rm = %.1f pu\n'], Zb_x, opt.xfm.R1/Zb_x, w0n*opt.xfm.L1/Zb_x, opt.xfm.Rm/Zb_x);


%% ------------------------------------------------------- 2) Knotenliste
fprintf('-> Schritt 2: Knotenliste aufbauen\n');

if opt.mergeSeriesNodes
    %      Name    Vbase[V]  Bemerkung
    nodeTab = { ...
        '632' , 4160 , 'DEG + Last; enthaelt 644 (T1)'            ; ...
        '645' , 4160 , ''                                         ; ...
        '646' , 4160 , 'Stichleitungsende'                        ; ...
        '671' , 4160 , 'enthaelt 692 (T2), drei Lasten'           ; ...
        '684' , 4160 , 'Verzweigung, KEINE Maschine'              ; ...
        '611' , 4160 , 'Last + Kondensatorbank (C_kk)'            ; ...
        '652' , 4160 , 'Stichleitungsende'                        ; ...
        '675' , 4160 , 'Last + Kondensatorbank (C_kk)'            ; ...
        '634' ,  480 , 'NS hinter XFM-634, Last + VSC'            ; ...
        '672' ,  480 , 'NS hinter XFM-671, Last + VSC'            ; ...
        '675L',  480 , 'NS hinter XFM-675, Last + VSC'            ; ... % hier ist 675 gemeint
        '680' ,  480 , 'NS hinter XFM-680, Last + VSC'            };
    LVbei680 = '680';
else
    nodeTab = { ...
        '632' , 4160 , 'DEG + Last; enthaelt 644 (T1)'            ; ...
        '633' , 4160 , 'OS-Punkt vor XFM-634, keine Maschine'     ; ...
        '645' , 4160 , ''                                         ; ...
        '646' , 4160 , 'Stichleitungsende'                        ; ...
        '671' , 4160 , 'enthaelt 692 (T2), drei Lasten'           ; ...
        '684' , 4160 , 'Verzweigung, KEINE Maschine'              ; ...
        '611' , 4160 , 'Last + Kondensatorbank (C_kk)'            ; ...
        '652' , 4160 , 'Stichleitungsende'                        ; ...
        '675' , 4160 , 'Last + Kondensatorbank (C_kk)'            ; ...
        '680' , 4160 , 'OS-Punkt vor XFM-680, keine Maschine'     ; ...
        '634' ,  480 , 'NS hinter XFM-634, Last + VSC'            ; ...
        '672' ,  480 , 'NS hinter XFM-671, Last + VSC'            ; ...
        '675L',  480 , 'NS hinter XFM-675, Last + VSC'            ; ...
        '680L',  480 , 'NS hinter XFM-680, Last + VSC'            };
    LVbei680 = '680L';
end

nodeName = nodeTab(:,1).';
nodeVb   = cell2mat(nodeTab(:,2)).'; % Knotenspannung in Matrix legen
N_node        = numel(nodeName);             % N = Anzahl der Knoten
kz       = (V_ref./nodeVb).^2;     % Impedanzfaktor Knotenebene -> V_ref % K =
                                   % siehe eq. circuit of trafo 
                               %  K = I1/I2 = V2/V1 = V_node/V_ref = kz^-1

fprintf('   N = %d Knoten  ->  Y_b ist %dx%d Bloecke (2x2), ausgeschrieben %dx%d\n', ...
        N_node, N_node, N_node, 2*N_node, 2*N_node);
fprintf('   %s\n', strjoin(nodeName, ', '));


%% ---------------- 3) Zusammenfassung der elektrischen Elemente je Knotenpaar
fprintf('-> Schritt 3: Laengszweige aufbauen\n');
% Widerstand und Kapazität des Ersatzschaltbilds des Trafos
Rx = opt.xfm.R1 + opt.xfm.R2*(opt.xfm.V1/opt.xfm.V2)^2;   % auf OS bezogen
Lx = opt.xfm.L1 + opt.xfm.L2*(opt.xfm.V1/opt.xfm.V2)^2;

internodeParams = struct('from',{},'to',{},'R',{},'L',{},'tag',{});
add = @(f,t,R,L,tg) struct('from',f,'to',t,'R',R,'L',L,'tag',tg);
pl  = @(f) para_grid.(f); % shortcut für die Übersichtlichkeit

% Prüft ob ... im struct para_grid existiert
for f = {'l644_645','l645_646','l633_644','l632_671','l671_684', ...
         'l652_684','l611_684','l675_692','l671_680'}
    assert(isfield(para_grid,f{1}), 'para_grid.%s fehlt.', f{1});
end

% Erstellt struct mit den Feldern für Gesamtwiderstände und Gesamtinduktivitäten zwischen den Knoten
internodeParams(end+1) = add('632','645', pl('l644_645').R, pl('l644_645').L, 'Leitung 632-645 (#603, 500 ft)');
internodeParams(end+1) = add('645','646', pl('l645_646').R, pl('l645_646').L, 'Leitung 645-646 (#603, 300 ft)');
internodeParams(end+1) = add('632','671', pl('l632_671').R, pl('l632_671').L, 'Leitung 632-671 (#601, 2000 ft)');
internodeParams(end+1) = add('671','684', pl('l671_684').R, pl('l671_684').L, 'Leitung 671-684 (#604, 300 ft)');
internodeParams(end+1) = add('684','652', pl('l652_684').R, pl('l652_684').L, 'Leitung 684-652 (#607, 800 ft)');
internodeParams(end+1) = add('684','611', pl('l611_684').R, pl('l611_684').L, 'Leitung 684-611 (#605, 300 ft)');
internodeParams(end+1) = add('671','675', pl('l675_692').R, pl('l675_692').L, 'Leitung 671-675 (#606, 500 ft)');
internodeParams(end+1) = add('671','672', Rx, Lx, 'XFM-671 Streuimpedanz'); % zw. 671 und 672 nur ein Trafo
internodeParams(end+1) = add('675','675L', Rx, Lx, 'XFM-675 Streuimpedanz'); % vor und hinter Trafo XFM-675

if opt.mergeSeriesNodes
    % Fasst 633-644 zusammen mit Trafo
    internodeParams(end+1) = add('632','634', pl('l633_644').R + Rx, pl('l633_644').L + Lx, ...
                    'Leitung 632-633 + XFM-634, zusammengefasst (M2)');
    % Fasst 671-680 zusammen mit Trafo
    internodeParams(end+1) = add('671','680', pl('l671_680').R + Rx, pl('l671_680').L + Lx, ...
                    'Leitung 671-680 + XFM-680, zusammengefasst (M2)');
else
    internodeParams(end+1) = add('632','633', pl('l633_644').R, pl('l633_644').L, 'Leitung 632-633 (#602, 500 ft)');
    internodeParams(end+1) = add('633','634', Rx, Lx, 'XFM-634 Streuimpedanz');
    internodeParams(end+1) = add('671','680', pl('l671_680').R, pl('l671_680').L, 'Leitung 671-680 (#601, 1000 ft)');
    internodeParams(end+1) = add('680','680L', Rx, Lx, 'XFM-680 Streuimpedanz');
end

if strcmpi(opt.DEG,'node650')
    error('Yb:deg','opt.DEG = ''node650'' wird in dieser Fassung nicht mehr gefuehrt.');
end

nBr = numel(internodeParams);
fprintf('   %d Zweige bei %d Knoten', nBr, N_node);
if nBr == N_node-1
    fprintf('  ->  N-1, also ein Baum (radiales Netz), wie erwartet\n');
else
    fprintf('  ->  %d unabhaengige Maschen\n', nBr-N_node+1);
end


%% ------------------- 4) Netzknotenelemente: passive Elemente (-> Y_b) und Maschinen und Lasten (-> Z_m)
fprintf('-> Schritt 4: Querglieder bestimmen\n');

%  Die 13 Lasten, die Maschinen sind.  Spalte: Knoten | Feld | Schaltung | Name
loadTab = { ...
    '632' ,'n632YgPQ'  ,'Yg','632 Yg PQ'                 ; ...
    '645' ,'n645YgPQ'  ,'Yg','645 Yg PQ'                 ; ...
    '646' ,'n646YgPQ'  ,'Yg','646 Yg PQ'                 ; ...
    '652' ,'n652YgZ'   ,'Yg','652 Yg Z'                  ; ...
    '611' ,'n611YgL'   ,'Yg','611 Yg L'                  ; ...
    '671' ,'n671YgPQ'  ,'Yg','671 Yg PQ'                 ; ...
    '671' ,'n671DPQ'   ,'D' ,'671 D PQ'                  ; ...
    '671' ,'n692DL'    ,'D' ,'692 D L (692==671)'        ; ...
    '675' ,'n675YgPQa' ,'Yg','675 Yg PQa'                ; ...
    '634' ,'n634YgPQ'  ,'Yg','634 Yg PQ 1 (schaltbar)'   ; ...
    '672' ,'n672YgPQ'  ,'Yg','672 Yg PQ   (schaltbar)'   ; ...
    '675L','n675YgPQb' ,'Yg','675 Yg PQb  (schaltbar)'   ; ...
    LVbei680,'n680YgPQ','Yg','680 Yg PQ   (schaltbar)'   };

%  Die 2 Kondensatorbaenke: keine Maschinen, sie gehoeren ins Netz (M1)
capTab = { ...
    '611' ,'n611YgZ'   ,'Yg','611 Yg Z (Kondensatorbank)' ; ...
    '675' ,'n675YgZ'   ,'Yg','675 Yg Z (Kondensatorbank)' };

if ~opt.capBanksInYb          % zurueck in die Maschinenmatrix
    loadTab = [loadTab; capTab];  capTab = cell(0,4);
end

% --- 4a) Lasten -> Maschinen

% Schaut ob die Lasten an den Knoten aus loadTab in nodeName gefunden
% werden können und ordnet einer Last die Knoten-ID nk zu
% 
% Berechnet die Parameter RLC der Reihen-RLC-Last aus den Leistungswerten 
% und speichert diese in load_ ab

load_ = struct('name',{},'field',{},'node',{},'config',{},'nk',{}, ...
               'R',{},'L',{},'C',{},'P3',{},'Q3',{},'QC3',{},'pf',{});

for i = 1:size(loadTab,1)
    nk = nodeIdx(nodeName, loadTab{i,1});
    if isempty(nk)
        fprintf('   [uebersprungen] %-22s (Knoten %s fehlt)\n', loadTab{i,4}, loadTab{i,1});
        continue
    end
    fn = loadTab{i,2};
    if ~isfield(para_grid,fn)
        warning('Yb:load','para_grid.%s fehlt - "%s" ignoriert.', fn, loadTab{i,4});
        continue
    end
    useLL = strcmp(loadTab{i,3},'D') && strcmpi(opt.deltaLoadV,'LL');
    % pqc gibt Wirkleistung, induktive und kapazitive Blindleistung zurück 
    [Ri,Li,Ci,pqc] = rlc_from_pq(para_grid.(fn), w0n, useLL, loadTab{i,4});
    S = hypot(3*pqc(1), 3*(pqc(2)-pqc(3)));
    load_(end+1) = struct('name',loadTab{i,4},'field',fn,'node',loadTab{i,1}, ...
        'config',loadTab{i,3},'nk',nk, ...
        'R',Ri*kz(nk),'L',Li*kz(nk),'C',Ci/kz(nk), ...
        'P3',3*pqc(1),'Q3',3*pqc(2),'QC3',3*pqc(3), ...
        'pf',ternary(S>0, 3*pqc(1)/max(S,eps), NaN));                  %#ok<SAGROW>
end

% --- 4b) Netz-Querglieder: Kondensatorbaenke (+ optional Magnetisierung)
Gnet = zeros(1,N_node); iLnet = zeros(1,N_node); Cnet = zeros(1,N_node);
netParts = cell(1,N_node);  for k = 1:N_node, netParts{k} = {}; end

% Sucht 
% Schleife macht (whrs.) 2 Durchläufe
for i = 1:size(capTab,1)
    nk = nodeIdx(nodeName, capTab{i,1});
    if isempty(nk) || ~isfield(para_grid, capTab{i,2}), continue; end
    % Lädt aus dem para_grid-Struct die Felder V_LN_rms, p, QL, QC und
    % berechnet daraus die Parameter einer (RLC-) Last R, L und C
    [Ri,Li,Ci] = rlc_from_pq(para_grid.(capTab{i,2}), w0n, false, capTab{i,4});
    assert(Ri==0 && Li==0 && isfinite(Ci), ['"%s" ist kein reiner ' ...
        'Kondensator (P oder QL ungleich null) und passt dann nicht in ' ...
        'C_kk nach Gl.(17). Bitte opt.capBanksInYb = false setzen.'], capTab{i,4});
    % kz ist ein Faktor um die Kapazität Ci auf das Spannungsiveau von
    % V_ref = 4160 zu bringen. Umrechnung gilt für Reaktanzen, aber auch für Kapazitäten?
    % --> kz() ist 1 in beiden Fällen, sodass bei der Platzierung der
    % Kapazität auf Spannungsebenbe V_ref = 4160 V kein Rechenfehler
    % entsteht
    fprintf('(Potentielle Fehlerquelle: Umrechnung Kapazität)\n')
    Cnet(nk)  = Cnet(nk)  + Ci/kz(nk);
    netParts{nk}{end+1} = capTab{i,4};
end


% Berücksichtigung der Magnetisierungszweige; default -> "nein"
if opt.xfmMag
    magTab = {'671','XFM-671'; '675','XFM-675'};
    if ~opt.mergeSeriesNodes
        magTab = [magTab; {'633','XFM-634'}; {'680','XFM-680'}];
    end
    for i = 1:size(magTab,1)
        nk = nodeIdx(nodeName, magTab{i,1});
        if isempty(nk), continue; end
        Gnet(nk)  = Gnet(nk)  + 1/opt.xfm.Rm;
        iLnet(nk) = iLnet(nk) + 1/opt.xfm.Lm;
        netParts{nk}{end+1} = [magTab{i,2} ' Magnetisierung'];
    end
end

% --- 4c) Lasten gehen als eigene 2x2-Bloecke nach Z_m (Abschnitt 6b).
%     Eine Reihen-RLC-Last laesst sich nicht einfach in ein gemeinsames R_kk/L_kk/
%     C_kk je Knoten zusammenfassen wie in Gl. (17) (Gu et al.)
%     Es muss dein Umweg über die 
%
%     Bei defaults: capBanksInYb = true, xfmMag = false
%     --> Gnet = 0, iLnet = 0
%     --> hasR = false, hasL = false
%     --> flag hasC ist zweimal true, für Knoten 611 und 675, wo
%     Kondensatorbänke liegen
Gsh = Gnet;  iLsh = iLnet;  Csh = Cnet;
fprintf('   %d Lasten -> Z_m, %d Kondensatorbaenke -> Y_b\n', ...
        numel(load_), size(capTab,1));

Rkk_n = 1./Gsh;   Rkk_n(Gsh==0)  = Inf;
Lkk_n = 1./iLsh;  Lkk_n(iLsh==0) = Inf;
Ckk_n = Csh;
hasR = isfinite(Rkk_n);
hasC = (Ckk_n ~= 0);
hasL = isfinite(Lkk_n) & opt.shuntL;


%% ---------------------------------------- 5) Symbole anlegen (class sym)
fprintf('-> Schritt 5: symbolische Parameter anlegen\n');

R_kl = sym(zeros(N_node,N_node));   hasRl = false(N_node,N_node);
L_kl = sym(zeros(N_node,N_node));   hasLl = false(N_node,N_node);
R_kk = sym(zeros(N_node,1));  C_kk = sym(zeros(N_node,1));  L_kk = sym(zeros(N_node,1));
% par enthält
par  = struct('name',{{}},'sym',sym([]),'val',[]);


% R_kl ist ein gute Variable für einen optische "quick check"
for i = 1:numel(internodeParams)
    k = nodeIdx(nodeName, internodeParams(i).from);
    l = nodeIdx(nodeName, internodeParams(i).to);
    if isempty(k) || isempty(l), continue; end
    if internodeParams(i).R ~= 0
        smbR = sym(sprintf('R_%s_%s', nodeName{k}, nodeName{l}));
        R_kl(k,l) = smbR; R_kl(l,k) = smbR; hasRl(k,l) = true; hasRl(l,k) = true;
        par = addpar(par, smbR, internodeParams(i).R);
    end
    if internodeParams(i).L ~= 0
        smbL = sym(sprintf('L_%s_%s', nodeName{k}, nodeName{l}));
        L_kl(k,l) = smbL; L_kl(l,k) = smbL; hasLl(k,l) = true; hasLl(l,k) = true;
        par = addpar(par, smbL, internodeParams(i).L);
    end
end
hasBr = hasRl | hasLl;

for k = 1:N_node
    if hasR(k)
        smbR = sym(sprintf('R_%s_%s', nodeName{k}, nodeName{k}));
        R_kk(k) = smbR;  par = addpar(par, smbR, Rkk_n(k));
    end
    if hasC(k)
        smbC = sym(sprintf('C_%s_%s', nodeName{k}, nodeName{k}));
        C_kk(k) = smbC;  par = addpar(par, smbC, Ckk_n(k));
    end
    if hasL(k)
        smbL = sym(sprintf('L_%s_%s', nodeName{k}, nodeName{k}));
        L_kk(k) = smbL;  par = addpar(par, smbL, Lkk_n(k));
    end
end
fprintf('   %d Netzparameter\n', numel(par.name));


%% ----------------------------- 6) Y_bkl und Y_bkk aufbauen (Gl. 16 / 17)
fprintf('-> Schritt 6: Y_bkl (Gl.16) und Y_bkk (Gl.17) aufbauen\n');

Ybkl = cell(N_node,N_node);
for k = 1:N_node
    for l = 1:N_node
        if k==l || ~hasBr(k,l), continue; end
        Zp  = R_kl(k,l) + (s + 1i*w0)*L_kl(k,l);
        Zm_ = R_kl(k,l) + (s - 1i*w0)*L_kl(k,l);
        Ybkl{k,l} = -[ 1/Zp , 0 ; 0 , 1/Zm_ ];                    % Gl. (16)
    end
end

Ybkk = cell(N_node,1);
for k = 1:N_node
    acc = sym(zeros(2,2));
    for l = 1:N_node
        if l==k || isempty(Ybkl{k,l}), continue; end
        acc = acc - Ybkl{k,l};
    end
    yp = sym(0);  ym = sym(0);
    if hasR(k), yp = yp + 1/R_kk(k);             ym = ym + 1/R_kk(k);          end
    if hasC(k), yp = yp + (s+1i*w0)*C_kk(k);     ym = ym + (s-1i*w0)*C_kk(k);  end
    if hasL(k), yp = yp + 1/((s+1i*w0)*L_kk(k)); ym = ym + 1/((s-1i*w0)*L_kk(k)); end
    acc = acc + [ yp , 0 ; 0 , ym ];                              % Gl. (17)
    Ybkk{k} = acc;  Ybkl{k,k} = acc;
end

% Macht aus der N x N cell array (mit 2x2 Matrizen je cell) eine
% 2N x 2N - Matrix
Yb = blk2full(Ybkl);      % ausgeschriebene Form, nur fuer die Numerik


%% ------------------------------- 6b) Z_m der Maschinen (Gl. 14)
% Besetzt die Diagonalelemente der 2x2 Matrix jeder Maschinenimpedanz
% Erstellt Q'D'-Matrix für die Lasten aus ihren zuvor berechneten Werten
% für Widerstand, Induktivität und Kapazität
fprintf('-> Schritt 6b: Z_m aufbauen (Lasten, DEG, VSC)\n');

for i = 1:numel(load_)
    zp = sym(0);  zm = sym(0);  tagn = load_(i).field;   % Reihenimpedanz
    if load_(i).R ~= 0
        smbR = sym(['R_' tagn]);  par = addpar(par, smbR, load_(i).R);
        zp = zp + smbR;                zm = zm + smbR;                load_(i).symR = smbR;
    end
    if load_(i).L ~= 0
        smbL = sym(['L_' tagn]);  par = addpar(par, smbL, load_(i).L);
        zp = zp + (s+1i*w0)*smbL;      zm = zm + (s-1i*w0)*smbL;      load_(i).symL = smbL;
    end
    % für alle i ist load_(i).C = inf
    if isfinite(load_(i).C)
        smbC = sym(['C_' tagn]);  par = addpar(par, smbC, load_(i).C);
        zp = zp + 1/((s+1i*w0)*smbC);  zm = zm + 1/((s-1i*w0)*smbC);  load_(i).symC = smbC;
    end
    % fügt jedem Eintrag (jeder Last) in load_ einen Impedanzmatrix und
    % eine Admittanzmatrix an
    load_(i).Zm = [zp 0; 0 zm];      % Gl.(14): Z_mk entsteht direkt
    load_(i).Ym = [1/zp 0; 0 1/zm];  % das, was in Y_m und damit in Gl.(21) geht
end

Ym_node = cell(N_node,1);  for k = 1:N_node, Ym_node{k} = sym(zeros(2,2)); end
Zm_node = cell(N_node,1);  for k = 1:N_node, Zm_node{k} = sym(zeros(2,2)); end
machParts = cell(1,N_node);  for k = 1:N_node, machParts{k} = {}; end

% an manchen Knoten sind mehr als eine Last/Maschine angebunden
% Mit der Schleife, werden je Knoten mit der ID nk die Admittanzen/Impedanzen der
% Lasten und die der VSC gesammelt
for i = 1:numel(load_)
    k = load_(i).nk;
    Ym_node{k} = Ym_node{k} + load_(i).Ym;
    Zm_node{k} = Zm_node{k} + load_(i).Zm;
    machParts{k}{end+1} = load_(i).name;
end

% Abbildung des dynamic equivalent of extended grid durch eine Admittanz
if strcmpi(opt.DEG,'machine632') || strcmpi(opt.DEG,'shunt632')
    k = nodeIdx(nodeName,'632');
    R_DEG = sym('R_DEG');  L_DEG = sym('L_DEG');
    par = addpar(par, R_DEG, Rg);  par = addpar(par, L_DEG, Lg);
    Ydeg = [ 1/(R_DEG+(s+1i*w0)*L_DEG) ,              0            ; ...
                          0            , 1/(R_DEG+(s-1i*w0)*L_DEG) ];
    Zdeg = [ Ydeg(1,1)^-1 ,              0            ; ...
                          0            , Ydeg(2,2)^-1 ];
    if strcmpi(opt.DEG,'shunt632')
        % diese Verzweigung wurde nicht redigiert
        Ybkk{k} = Ybkk{k} + Ydeg;  Ybkl{k,k} = Ybkk{k};
        Yb(2*k-1:2*k, 2*k-1:2*k) = Ybkk{k};
        fprintf('   DEG als Serien-RL-Querglied in Y_b (fuer die SPS-Validierung)\n');
    else
        Ym_node{k} = Ym_node{k} + Ydeg;
        % Fehler: die Impedanzen können je Knoten nicht einfach addiert
        % werden
        %Zm_node{k} = Zm_node{k} + Zdeg;
        machParts{k}{end+1} = 'DEG (Rg, Lg)';
        fprintf('   DEG als Maschine Z_m am Knoten 632\n');
    end
end

% --- VSC: 2x2-Platzhalter. Feldnamen mit 'n'-Praefix, weil MATLAB-Struktur-
%     felder nicht mit einer Ziffer beginnen duerfen. assemble_Ysys.m
%     benutzt dieselbe Konvention.
vscTab = {'634','n634','Y_VSC_634','Z_VSC_634'; '672','n672','Y_VSC_672','Z_VSC_672'; ...
          '675L','n675L','Y_VSC_675L','Z_VSC_675L'; LVbei680,['n' LVbei680],'Y_VSC_680','Z_VSC_680'};
Y_VSC = struct();  Y_VSCnode = struct();
Z_VSC = struct();  Z_VSCnode = struct();
for i = 1:size(vscTab,1)
    k = nodeIdx(nodeName, vscTab{i,1});
    if isempty(k), continue; end
    Ypl = sym(zeros(2,2));
    Zpl = sym(zeros(2,2));
    Ypl(1,1) = sym([vscTab{i,3} '_pp']);  Ypl(1,2) = sym([vscTab{i,3} '_pm']);
    Ypl(2,1) = sym([vscTab{i,3} '_mp']);  Ypl(2,2) = sym([vscTab{i,3} '_mm']);
    Zpl(1,1) = sym([vscTab{i,4} '_pp']);  Zpl(1,2) = sym([vscTab{i,4} '_pm']);
    Zpl(2,1) = sym([vscTab{i,4} '_mp']);  Zpl(2,2) = sym([vscTab{i,4} '_mm']);
    Y_VSC.(vscTab{i,2})     = Ypl;
    Z_VSC.(vscTab{i,2})     = Zpl;
    Y_VSCnode.(vscTab{i,2}) = k;
    Z_VSCnode.(vscTab{i,2}) = k;
    % NICHT automatisch addiert - sonst stehen freie Symbole in Yb_plus_Ym.
    Ym_node{k} = Ym_node{k} + Ypl;
end
fprintf('   VSC-Platzhalter: %s (eingehaengt)\n', ...
        strjoin(fieldnames(Y_VSC).', ', '));

% Hinweis: nachfolgender Block nicht bebraucht, falls Zm_node vorher
% symbolisch besetzt wird
%
% for k = 1:N_node
%     if isAlwaysZero(Ym_node{k}), Zm_node{k} = sym(inf(2,2));
%     else, Zm_node{k} = [1/Ym_node{k}(1,1) 0; 0 1/Ym_node{k}(2,2)];
%     end
% end

% "Übersetzt" die Admittanzen der Maschinen an den 12 Knoten
% in eine Diagonalmatrix, der Maschinenmatrix (mit Admittanzen)
Y_m = sym(zeros(2*N_node,2*N_node));
for k = 1:N_node
    Y_m(2*k-1:2*k, 2*k-1:2*k) = Ym_node{k}; 
end
Yb_plus_Ym = Yb + Y_m;


%% ------------------------------------------ 7) Zahlenwerte einsetzen
% Anmerkung:
% exprnew = subs(expr,match,replacement) 
% replaces all occurrences of match in expr with replacement, evaluates the expression, and returns 
% the result as a symbolic expression
fprintf('-> Schritt 7: Zahlenwerte einsetzen\n');

if opt.exactSym
    vals = sym(zeros(1,numel(par.val)));
    for i = 1:numel(par.val), vals(i) = sym(par.val(i),'f'); end
else
    vals = par.val;
end

% Redigiert bis hier hin 19.10.2026

Yb_num = Yb;
%Ym_num = Ym;
fprintf('Die Matrix der Admittanzen der Lasten und VSC Ym liegt weiterhin nur symbolisch vor.\n');

YbYm_num = Yb_plus_Ym;
if ~isempty(par.val)
    Yb_num   = subs(Yb_num,   par.sym, vals);
    YbYm_num = subs(YbYm_num, par.sym, vals);
end
if opt.w0_symbolic
    Yb_num   = subs(Yb_num,   w0, 2*sym(pi)*sym(fg));
    YbYm_num = subs(YbYm_num, w0, 2*sym(pi)*sym(fg));
end

Yb_fun = [];  YbYm_fun = [];
if opt.makeFun
    t0 = tic;
    Yb_fun   = matlabFunction(Yb_num,  'Vars', s);
    %YbYm_fun = matlabFunction(YbYm_num,'Vars', s);
    fprintf('   matlabFunction erzeugt (%.1f s)\n', toc(t0));
end


%% HINWEIS: VERIFIKATION AUSGESETZT
%
% tbd (261019): definiere Gld, iLld und Cld
%
% %% ------------------------------------------------- 8) Verifikation
% fprintf('-> Schritt 8: Verifikation\n');
% 
% maxRes = 0;
% for k = 1:N_node
%     rs = sym(zeros(2,2));
%     % Zeile 2k-1  ist die +Komponente von Knoten k
%     % Zeile 2k    ist die -Komponente von Knoten k
% 
%     % bei k=1: Y(1:2,1:2) + Y(1:2,3:4) + Y(1:2,5:6) + ...
%     % 
%     for l = 1:N_node
% 
%         rs = rs + Yb_num(2*k-1:2*k, 2*l-1:2*l); % Argumente liefern aus Yb_num eine 2x2-Matrix zurück
%     end
%     yp = 0; ym = 0;
%     if hasR(k), yp = yp + 1/Rkk_n(k);          ym = ym + 1/Rkk_n(k);          end
%     if hasC(k), yp = yp + 1i*w0n*Ckk_n(k);     ym = ym - 1i*w0n*Ckk_n(k);     end
%     if hasL(k), yp = yp + 1/(1i*w0n*Lkk_n(k)); ym = ym + 1/(-1i*w0n*Lkk_n(k));end
%     ref = [yp 0; 0 ym];
%     if strcmpi(opt.DEG,'shunt632') && strcmp(nodeName{k},'632')
%         ref = ref + [1/(Rg+1i*w0n*Lg) 0; 0 1/(Rg-1i*w0n*Lg)];
%     end
%     maxRes = max(maxRes, max(abs(reshape(double(subs(rs,s,0))-ref,[],1))));
% end
% fprintf('   (a) Blockzeilensumme Y_b (Gl.17)    : max. Abw. %.3e\n', maxRes);
% 
% Yph  = zeros(N_node,N_node);
% Gall = Gnet + Gld;  iLall = iLnet + iLld;  Call = Cnet + Cld;
% for i = 1:numel(internodeParams)
%     k = nodeIdx(nodeName, internodeParams(i).from);  l = nodeIdx(nodeName, internodeParams(i).to);
%     if isempty(k) || isempty(l), continue; end
%     y = 1/(internodeParams(i).R + 1i*w0n*internodeParams(i).L);
%     Yph(k,l) = Yph(k,l) - y;  Yph(l,k) = Yph(l,k) - y;
%     Yph(k,k) = Yph(k,k) + y;  Yph(l,l) = Yph(l,l) + y;
% end
% for k = 1:N_node
%     if Gall(k)  ~= 0, Yph(k,k) = Yph(k,k) + Gall(k);            end
%     if Call(k)  ~= 0, Yph(k,k) = Yph(k,k) + 1i*w0n*Call(k);     end
%     if iLall(k) ~= 0 && opt.shuntL
%                       Yph(k,k) = Yph(k,k) + iLall(k)/(1i*w0n);  end
% end
% if ~strcmpi(opt.DEG,'none')
%     k = nodeIdx(nodeName,'632');  Yph(k,k) = Yph(k,k) + 1/(Rg+1i*w0n*Lg);
% end
% Y0 = double(subs(YbYm_num, s, 0));
% Yp = Y0(1:2:end,1:2:end);  Ymn = Y0(2:2:end,2:2:end);
% fprintf('   (b) Y_b+Y_m (+) bei s=0 vs. Phasor-Y: max. Abw. %.3e\n', max(abs(Yp(:)-Yph(:))));
% fprintf('   (c) Y_b+Y_m (-) bei s=0 vs. conj    : max. Abw. %.3e\n', max(abs(Ymn(:)-conj(Yph(:)))));
% fprintf('   (d) +/- -Kopplung                   : %.3e\n', max(max(abs(Y0(1:2:end,2:2:end)))));
% 
% ftest = 137;  wt = 2*pi*ftest;
% Ysh = double(subs(YbYm_num, s, 1i*(wt-w0n)));  Ysh = Ysh(1:2:end,1:2:end);
% Yst = zeros(N_node,N_node);
% for i = 1:numel(internodeParams)
%     k = nodeIdx(nodeName, internodeParams(i).from);  l = nodeIdx(nodeName, internodeParams(i).to);
%     if isempty(k) || isempty(l), continue; end
%     y = 1/(internodeParams(i).R + 1i*wt*internodeParams(i).L);
%     Yst(k,l) = Yst(k,l)-y;  Yst(l,k) = Yst(l,k)-y;
%     Yst(k,k) = Yst(k,k)+y;  Yst(l,l) = Yst(l,l)+y;
% end
% for k = 1:N_node
%     if Gall(k)  ~= 0, Yst(k,k) = Yst(k,k) + Gall(k);           end
%     if Call(k)  ~= 0, Yst(k,k) = Yst(k,k) + 1i*wt*Call(k);     end
%     if iLall(k) ~= 0 && opt.shuntL
%                       Yst(k,k) = Yst(k,k) + iLall(k)/(1i*wt);  end
% end
% if ~strcmpi(opt.DEG,'none')
%     k = nodeIdx(nodeName,'632');  Yst(k,k) = Yst(k,k) + 1/(Rg+1i*wt*Lg);
% end
% fprintf('   (e) Y_(+)(j(w-w0)) vs. Y_stat(jw), f=%g Hz: %.3e\n', ftest, max(abs(Ysh(:)-Yst(:))));
% 
% Yb0 = double(subs(Yb_num, s, 0));
% fprintf('   (f) Y_b symmetrisch (Reziprozitaet) : %.3e\n', max(max(abs(Yb0-Yb0.'))));
% fprintf('   (g) Rang Y_b(0) = %d/%d, Rang Y_b+Y_m(0) = %d/%d\n', ...
%         rank(Yb0), 2*N_node, rank(Y0), 2*N_node);
% 
% 
% %% ------------------------------------------------------- 9) Report
% if opt.verbose
%     fprintf('\n--- Laengszweige (R_kl, L_kl) -------------------------------\n');
%     fprintf('%-6s %-6s %13s %14s %8s   %s\n','von','nach','R_kl [Ohm]','L_kl [H]','X/R','Bezeichnung');
%     for i = 1:numel(internodeParams)
%         if internodeParams(i).R ~= 0, xr = w0n*internodeParams(i).L/internodeParams(i).R; else, xr = Inf; end
%         fprintf('%-6s %-6s %13.6g %14.6g %8.3f   %s\n', ...
%                 internodeParams(i).from, internodeParams(i).to, internodeParams(i).R, internodeParams(i).L, xr, internodeParams(i).tag);
%     end
% 
%     fprintf('\n--- Netz-Querglieder in Y_b (R_kk, L_kk, C_kk) --------------\n');
%     if ~any(hasR|hasL|hasC), fprintf('   (keine)\n'); end
%     for k = 1:N_node
%         if hasR(k) || hasL(k) || hasC(k)
%             fprintf('%-6s %13.6g %14.6g %14.6g   %s\n', nodeName{k}, ...
%                     Rkk_n(k), Lkk_n(k), Ckk_n(k), strjoin(netParts{k},' + '));
%         end
%     end
% 
%     fprintf('\n--- Maschinen je Knoten -------------------------------------\n');
%     nMach = 0;
%     for k = 1:N_node
%         nv = '';
%         if isfield(Y_VSCnode, ['n' nodeName{k}]), nv = ' + VSC (Platzhalter)'; end
%         if isempty(machParts{k}) && isempty(nv)
%             fprintf('%-6s  -  keine Maschine (reiner Verzweigungsknoten)\n', nodeName{k});
%         else
%             fprintf('%-6s  %s%s\n', nodeName{k}, strjoin(initmachParts{k},' + '), nv);
%             nMach = nMach + 1;
%         end
%     end
% 
%     fprintf('\n--- Die %d Lasten in Z_m, bezogen auf %.0f V -----------------\n', ...
%             numel(load_), V_ref);
%     fprintf('%-22s %-5s %-3s %12s %13s %8s\n','Last','Kn.','Sch','R [Ohm]','L [H]','cos(phi)');
%     for i = 1:numel(load_)
%         fl = '';
%         if ~isnan(load_(i).pf) && load_(i).pf < 0.5 && load_(i).P3 > 0, fl = '  <== pruefen'; end
%         fprintf('%-22s %-5s %-3s %12.6g %13.6g %8.3f%s\n', load_(i).name, ...
%                 load_(i).node, load_(i).config, load_(i).R, load_(i).L, load_(i).pf, fl);
%     end
%     fprintf('SUMME  P = %.0f kW, QL = %.0f kvar, QC = %.0f kvar\n', ...
%             sum([load_.P3])/1e3, sum([load_.Q3])/1e3, sum([load_.QC3])/1e3);
% 
%     fprintf('\nN = %d Knoten, %d davon mit Maschine, %d Zweige.\n', N_node, nMach, numel(internodeParams));
%     fprintf('Y_b = %dx%d Bloecke (2x2), ausgeschrieben %dx%d, %d Parameter.\n', ...
%             N_node, N_node, 2*N_node, 2*N_node, numel(par.name));
%     fprintf('Weiter mit assemble_Ysys.m (Gl.20/21 dual, xi, VSC).\n\n');
% end


%% ========================================================================
%  Lokale Funktionen
%% ========================================================================

function k = nodeIdx(nodeName, name)
k = find(strcmp(nodeName, name), 1);
end

function par = addpar(par, symv, val)
par.name{end+1} = char(symv);
par.sym(end+1)  = symv;
par.val(end+1)  = val;
end

function y = ternary(c,a,b)
if c, y = a; else, y = b; end
end

function answer = isAlwaysZero(M)
answer = isempty(symvar(M)) && all(double(M(:)) == 0);
end

function [R,L,C,pqc] = rlc_from_pq(ld, w0n, useLL, tag)
% Reihen-RLC-Ersatz einer symmetrischen SPS-Last aus P, QL, QC.
% Die Lastbloecke sind "Three-Phase Series RLC Load": R, L und C liegen in
% REIHE, derselbe Strom fliesst durch alle drei. Mit Q = QL - QC und
% |I|^2 = (P^2+Q^2)/V^2 folgt aus P=|I|^2 R, QL=|I|^2 w0 L, QC=|I|^2/(w0 C):
%   R = P *V^2/(P^2+Q^2)   L = QL*V^2/(w0*(P^2+Q^2))   C = (P^2+Q^2)/(w0*QC*V^2)
% Fehlendes Element ist ein KURZSCHLUSS: R=0, L=0, C=Inf (parallel waere es
% ein Leerlauf: Inf, Inf, 0). Bei f0 und Nennspannung liefern beide Varianten
% dieselbe Impedanz; bei Lasten mit nur einem Element sind sie identisch.
V  = getfld(ld,'V_LN_rms',NaN);
P  = getfld(ld,'P' ,0);
QL = getfld(ld,'QL',0);
QC = getfld(ld,'QC',0);

isasym = @(x) numel(x)>1 && max(abs(x-x(1))) > 1e-9*max(1,max(abs(x)));
if isasym(V) || isasym(P) || isasym(QL) || isasym(QC)
    warning('Yb:asym', ['Last "%s" ist unsymmetrisch. Es wird nur Phase a ' ...
        'verwendet - das 2x2-Modell kann Unsymmetrie nicht abbilden.'], tag);
end
V = V(1);  P = P(1);  QL = QL(1);  QC = QC(1);
if useLL, V = V/sqrt(3); end

assert(isfinite(V) && V>0, 'Last "%s": V_LN_rms fehlt oder ist null.', tag);
assert(P~=0 || QL~=0 || QC~=0, 'Last "%s": P, QL und QC sind alle null.', tag);

Q = QL - QC;  d = P^2 + Q^2;
if P  ~= 0, R = P *V^2/d;       else, R = 0;   end
if QL ~= 0, L = QL*V^2/(w0n*d); else, L = 0;   end
if QC ~= 0, C = d/(w0n*QC*V^2); else, C = Inf; end
pqc = [P QL QC];
end

function v = getfld(st, f, dflt)
if isfield(st,f) && ~isempty(st.(f)), v = st.(f); else, v = dflt; end
end
