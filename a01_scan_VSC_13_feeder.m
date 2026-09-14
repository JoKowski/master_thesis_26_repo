%% ========================================================================
%  a_01_impedanzscan.m
%
%  Impedanzscan eines GFL-Umrichters, qd-Domaenen-Injektion mit
%  Injektionswinkel nach Gl. (28)/(29) der Referenz, jedoch als
%  SPANNUNGS-Stoerung (nicht Shunt-Strom).
%
%      dv_q = amp_ptb * sin(omega_ptb*t + phs0_ptb) * cos(phi_ptb)
%      dv_d = amp_ptb * sin(omega_ptb*t + phs0_ptb) * sin(phi_ptb)
%
%  Beide Achsen teilen sich EINEN Oszillator -> die Stoerung ist eine
%  gerade Linie in der qd-Ebene unter dem Winkel phi_ptb, keine Ellipse.
%  Deshalb im Modell ein Sine-Wave-Block und zwei Gain-Bloecke, nicht
%  zwei getrennte Sinusquellen.
%
%      phi_ptb  : RAEUMLICHER Injektionswinkel in der qd-Ebene.
%                 Unterscheidet die beiden Laeufe, bestimmt das Ergebnis.
%      phs0_ptb : ZEITLICHE Phase des Sinus. In beiden Laeufen gleich,
%                 kuerzt sich in Y = I*V^-1 heraus. Einziger Zweck: der
%                 Sinus soll bei t_start auf null einsetzen.
%
%  Pro Frequenz zwei Laeufe, die sich NUR in phi_ptb unterscheiden.
%  Gemessen bei f_p:  V = [v^1 v^2], I = [i^1 i^2]  (Spalten = Laeufe)
%      Y_qd = I * V^-1        Z_qd = Y_qd^-1
%
%  KONDITIONIERUNG:
%      v^(k) = A*[cos(phi_k); sin(phi_k)]   ->   det(V) = A^2*sin(dphi)
%      cond(V) = cot(dphi/2)
%          dphi = 30 Grad (+15/-15)  ->  cond = 3.73
%          dphi = 60 Grad (+30/-30)  ->  cond = 1.73
%          dphi = 70 Grad (+35/-35)  ->  cond = 1.43
%          dphi = 90 Grad (+45/-45)  ->  cond = 1.00   (Optimum)
%      Der ABSOLUTE Winkel ist irrelevant - eine gemeinsame Drehung beider
%      Laeufe wird von V^-1 absorbiert. Nur die Differenz zaehlt.
%
%  REIHENFOLGE DER LAEUFE:
%      Von der HOECHSTEN zur NIEDRIGSTEN Stoerfrequenz. Die kurzen Fenster
%      laufen zuerst -> Fehler in der Messkette fallen nach Sekunden auf,
%      nicht erst nach der ersten 50-Sekunden-Simulation. Gespeichert wird
%      weiterhin nach aufsteigender Frequenz sortiert.
%
%  EINSCHWINGEN NACH INJEKTIONSBEGINN:
%      Der Zuschalttransient klingt mit den Eigenwerten des Systems ab, NICHT
%      mit 1/f_p. Deshalb:
%          t_ptb = min( t_ptb_max, max( t_ptb_min, n_ptb_per/f_p ) )
%      Die Periodenkopplung deckt den Fall ab, dass f_p nahe einer schwach
%      gedaempften Resonanz liegt; t_ptb_max verhindert, dass sie bei kleinen
%      f_p die Laufzeit dominiert.
%
%  Kohaerenzbedingungen (zwingend, sonst Leckage in der DFT):
%      T_win  = k/f0     (k ganz)     -> f0 exakt auf einem DFT-Bin
%      f_p    = m/T_win  (m ganz)     -> f_p exakt auf einem Bin
%      => (f0 +- f_p)*T_win = k +- m  -> beide abc-Toene ebenfalls exakt
%      Ts     = 1/(f0*N)              -> T_win/Ts = k*N ganzzahlig
%      D_log | N                      -> dezimiertes Raster bleibt exakt
%      t_start*f0, t_ptb*f0 ganz      -> Fensterstart faellt aufs Raster
%
%  MODELLAUFBAU IM STOERPFAD:
%      Sine Wave:  Amplitude = amp_ptb
%                  Frequency = omega_ptb     [rad/s]
%                  Phase     = phs0_ptb      [rad]
%         |-- Gain  cos(phi_ptb) --> dv_q --|
%         |-- Gain  sin(phi_ptb) --> dv_d --|--> Mux [dv_q; dv_d; 0]
%                          Constant 0 ------|
%                                            --> qd0 to abc  <- theta_sys
%                                            --> 3x Controlled Voltage Source
%      Freigabe ab t_start durch Multiplikation mit einem Step.
%
%  WEITERE VORAUSSETZUNGEN:
%      - powergui: Discrete, Sample time Ts
%      - alle PID/Regler-Bloecke diskret mit Sample time Ts
%      - zwei To-Workspace-Bloecke, Format Timeseries, Sample time Ts,
%        Decimation D_log, mit dem qd0-Vektor von Spannung und Strom am PCC
%      - HINWEIS: 'ode5' ist ein EXPLIZITER Fixed-Step-Solver. Damit gilt die
%        Stabilitaetsgrenze h < 2.8/|lambda_max|. Sind alle Zustaende diskret,
%        ist 'FixedStepDiscrete' schneller und hat keine solche Grenze.
%
%  Literatur:
%      Francis, Diss. Virginia Tech, 2010              (zwei unabh. Anregungen)
%      Wen et al.,     IEEE TPEL   31(1):675-687, 2016 (2x2-dq-Impedanz GFL)
%      Rygg et al.,    IEEE JESTPE  4(4):1383-1396, 2016 (qd <-> Sequenz)
%      Roinila et al., IEEE TPEL   33(5):4015-4022, 2018 (Konditionierung)
%      Thakar & Ramasubramanian, arXiv:2207.13641, 2022  (Fehlerquellen)
% =========================================================================
%clear; close all; clc;

%% ---- 0) Konfiguration ---------------------------------------------------
% Ergebnisdatei: Basisname + Startzeitstempel + Skriptversion
% -> impedanzscan_2026-09-02_14-35_vers_1.1.mat
SCRIPT_VERSION = '1.2';        % <<< bei jeder inhaltlichen Aenderung erhoehen %Ts jetzt ungleich T_ctrl
RUN_STAMP      = char(datetime('now','Format','yyyy-MM-dd_HH-mm'));
SAVE_BASE      = 'impedanzscan';
SAVE_FILE      = sprintf('%s_%s_vers_%s.mat', SAVE_BASE, RUN_STAMP, SCRIPT_VERSION);

mdl      = 'PV_Geno_VSC_GFL_634';   % <<< ANPASSEN: Modellname
SOLVER   = 'FixedStepDiscrete';                  % 'ode5' | 'FixedStepDiscrete'

fg       = 60;                       % Netzfrequenz [Hz]
Vp       = sqrt(2)*480/sqrt(3);      % Scheitelwert Phasenspannung [V]
psi_0    = 0;                        % Rahmenoffset der Clock-Rampe [rad]

t_settle  = 0.5;                     % Einschwingen bis Injektionsbeginn [s]
t_start   = t_settle;                % t_start*fg muss ganzzahlig sein
t_ptb_min = 0.2;                     % Mindest-Einschwingen NACH Injektion [s]
t_ptb_max = 1.0;                     % Deckel dafuer [s]  <-- neu
n_ptb_per = 8;                       % zusaetzlich mind. so viele Stoerperioden

amp_rel   = 0.01;                    % Stoeramplitude relativ zu Vp
amp_ptb   = amp_rel*Vp;              % = I_m der Referenz, hier als Spannung

phi_ptb_deg = [35, -35];             % Injektionswinkel je Lauf [Grad]
                                     % [45,-45] -> cond(V) = 1 (Optimum)

% Indizes im geloggten qd0-Vektor  <<< ANPASSEN an deine Busreihenfolge
IQ = 1;                              % q-Komponente
ID = 2;                              % d-Komponente

% ---- Signal Logging (out.logsout) ---------------------------------------
% Die Reihenfolge in logsout ist nicht garantiert, deshalb wird der Index
% beim ersten sim-Aufruf einmalig aufgeloest. ncol dient der Plausipruefung.
LOG = struct();
LOG.v_qd0_PCC = struct('name','Vqd0_634_pv','numInLogdata',[],'ncol',3);
LOG.i_qd0_PCC = struct('name','Iqd0_634_pv','numInLogdata',[],'ncol',3);
%LOG.v_qd0_ptb = struct('name','v_qd0_ptb', 'numInLogdata',[], 'ncol',3);
LOG.phi_clk   = struct('name','phi_clk',   'numInLogdata',[], 'ncol',1);

log_resolved = false;

% ---- Frequenzraster -----------------------------------------------------
%   'lin'    : aequidistant ue ber [f_lo, f_hi]. Wenige subsynchrone Punkte.
%   'log'    : logarithmisch. Haelfte der Punkte unterhalb f0.
%   'hybrid' : n_pts_log logarithmische Punkte bis f_split, danach linear.
GRID_MODE = 'log';
f_lo      = 0.1;                     % [Hz]
f_hi      = 5e3;                     % [Hz] % Nyquist-Frequenz bei 1/2 * (1/T_ctrl) = 10 kHz
n_pts     = 40;
f_split   = 60;                      % nur 'hybrid': Grenze log <-> lin
n_pts_log = 8;                       % nur 'hybrid': Punkte unterhalb f_split

M_min   = 5;                         % Mindestzahl Stoerperioden im Fenster
tol     = 0.02;                      % zul. rel. Abweichung vom Wunschraster
k_max   = 6000;                      % Obergrenze (T_win <= 100 s)
SHOW_GRID = true;                    % Rastertabelle im Preflight ausgeben

% Diskretisierungsguete
SPP_min  = 30;                       % Mindest-Samples/Periode bei f0+f_hi
SPP_log  = 16;                       % Mindest-Samples/Stoerperiode nach Dezimierung
f_alias  = 10*2*fg;                  % dezimierte Abtastrate mind. so hoch
T_ctrl   = 5e-5;

%% ---- 0b) Abgeleitete Groessen -------------------------------------------
f_abc_max = fg + f_hi;
Ts_max    = 1/(SPP_min*f_abc_max);  % obere Schranke aus der Diskretisierungsguete

% n_ctrl = Solverschritte je Reglerschritt. Kleinstes n, das beide
% Teilbarkeiten erfuellt: Ts <= Ts_max  und  T0/Ts ganzzahlig.
n_ctrl = ceil(T_ctrl/Ts_max);
while abs( (1/fg)*n_ctrl/T_ctrl - round((1/fg)*n_ctrl/T_ctrl) ) > 1e-9
    n_ctrl = n_ctrl + 1;
    assert(n_ctrl < 1e5, 'Keine passende Schrittweite gefunden.');
end
Ts   = T_ctrl/n_ctrl;
N_Ts = round(1/(fg*Ts));

assert(abs(T_ctrl/Ts - round(T_ctrl/Ts)) < 1e-12, 'T_ctrl kein Vielfaches von Ts.');
assert(abs(1/(fg*Ts) - N_Ts)             < 1e-9,  'T0 kein Vielfaches von Ts.');
phi_ptb_vec = phi_ptb_deg*pi/180;
dphi        = abs(diff(phi_ptb_vec));
cond_theo   = abs(cot(dphi/2));

assert(numel(phi_ptb_deg)==2, 'phi_ptb_deg muss genau zwei Winkel enthalten.');
assert(abs(sin(dphi)) > 1e-6, ...
       'Die beiden Injektionswinkel sind (nahezu) kollinear -> V singulaer.');
assert(abs(t_start*fg - round(t_start*fg)) < 1e-12, ...
       't_start*fg muss ganzzahlig sein.');

%% ---- 1) Frequenzraster mit Kohaerenzbedingungen -------------------------
switch lower(GRID_MODE)
    case 'lin'
        f_wish = linspace(f_lo, f_hi, n_pts);
    case 'log'
        f_wish = logspace(log10(f_lo), log10(f_hi), n_pts);
    case 'hybrid'
        f_wish = unique([ logspace(log10(f_lo), log10(f_split), n_pts_log), ...
                          linspace(f_split, f_hi, n_pts - n_pts_log + 1) ]);
    otherwise
        error('GRID_MODE unbekannt: %s', GRID_MODE);
end
n_pts = numel(f_wish);

k_vec = nan(1,n_pts);
m_vec = nan(1,n_pts);
f_p   = nan(1,n_pts);

for ii = 1:n_pts
    best_err = inf;
    k0 = max(1, ceil(M_min*fg/f_wish(ii)));
    for kk = 3*ceil(k0/3):3:k_max  %T_ctrol muss Grundperiode mit Netzfrequenz teilen
        T_run = kk/fg;
        m     = max(M_min, round(f_wish(ii)*T_run));
        err   = abs(m/T_run - f_wish(ii))/f_wish(ii);
        if err < best_err
            best_err = err;  k_vec(ii) = kk;  m_vec(ii) = m;  f_p(ii) = m/T_run;
        end
        if err <= tol, break; end
    end
    if best_err > tol
        warning('a_01:raster','f_wish = %.4f Hz: nur %.2f %% erreicht (k_max = %d).', ...
                f_wish(ii), 100*best_err, k_max);
    end
end

[f_p, iu] = unique(f_p,'stable');
k_vec = k_vec(iu);  m_vec = m_vec(iu);
n_f   = numel(f_p);
T_win = k_vec/fg;

% Einschwingzeit nach Injektionsbeginn, gedeckelt, auf f0-Raster aufgerundet
t_ptb  = min(t_ptb_max, max(t_ptb_min, n_ptb_per./f_p));
t_ptb  = ceil(t_ptb*fg)/fg;
t_win0 = t_start + t_ptb;

% Dezimierung je Punkt: groesster Teiler von N_Ts, der beide Grenzen haelt
div_N = divisors_of(N_Ts);
D_vec = ones(1,n_f);
for ii = 1:n_f
    fs_need   = max([SPP_log*f_p(ii), f_alias]);
    cand      = div_N( 1./(div_N*Ts) >= fs_need );
    D_vec(ii) = max([1, cand]);
end
Ts_log = D_vec*Ts;
N_win  = round(T_win./Ts_log);

% Simulationsreihenfolge: hoechste Frequenz zuerst
[~, run_order] = sort(f_p,'descend');

%% ---- 1b) Preflight ------------------------------------------------------
fprintf('\n=== Preflight ===\n');
fprintf('Skript    = a_01_impedanzscan  Version %s\n', SCRIPT_VERSION);
fprintf('Lauf      = %s\n', RUN_STAMP);
fprintf('Ergebnis  = %s\n', SAVE_FILE);
fprintf('Modell    = %s   Solver = %s\n', mdl, SOLVER);
fprintf('Vp        = %.4f V\n', Vp);
fprintf('amp_ptb   = %.2f %% von Vp = %.4f V\n', 100*amp_rel, amp_ptb);
fprintf('phi_ptb   = [%+.1f, %+.1f] Grad   ->  dphi = %.1f Grad\n', ...
        phi_ptb_deg(1), phi_ptb_deg(2), dphi*180/pi);
fprintf('erwartetes cond(V) = cot(dphi/2) = %.3f', cond_theo);
if cond_theo > 3
    fprintf('   <- mit [+45,-45] waere cond = 1.00\n');
else
    fprintf('\n');
end
fprintf('N_Ts      = %d   ->  Ts = %.6e s\n', N_Ts, Ts);
fprintf('Samples/Periode bei f0+f_hi = %.0f Hz:  %.1f\n', f_abc_max, 1/(f_abc_max*Ts));
assert(abs((1/fg)/Ts - round((1/fg)/Ts)) < 1e-9, 'T0/Ts nicht ganzzahlig.');

% Verteilung der Messpunkte sichtbar machen
fprintf('\nRaster    = ''%s'',  %d Punkte von %.4f bis %.1f Hz\n', ...
        GRID_MODE, n_f, min(f_p), max(f_p));
fprintf('davon subsynchron (< %.0f Hz): %d   ueber f0: %d\n', ...
        fg, sum(f_p<fg), sum(f_p>=fg));
e = 10.^(floor(log10(min(f_p))) : ceil(log10(max(f_p))));
c = histcounts(f_p, e);
fprintf('pro Dekade: ');
for jj = 1:numel(c)
    fprintf('%g-%g Hz: %d  ', e(jj), e(jj+1), c(jj));
end
fprintf('\n');

if SHOW_GRID
    fprintf('\n   #     f_p [Hz]      k      m   T_win [s]  t_ptb [s]  D_log   Samples\n');
    fprintf(  '  ---------------------------------------------------------------------\n');
    for jj = 1:n_f
        ii = run_order(jj);
        fprintf('  %3d  %11.4f  %5d  %5d  %9.4f  %8.4f  %5d  %8d\n', ...
                jj, f_p(ii), k_vec(ii), m_vec(ii), T_win(ii), t_ptb(ii), ...
                D_vec(ii), N_win(ii));
    end
    fprintf('\n');
end

ok = true;
for ii = 1:n_f
    c1 = abs(fg*T_win(ii)          - round(fg*T_win(ii)))          < 1e-9;
    c2 = abs(f_p(ii)*T_win(ii)     - round(f_p(ii)*T_win(ii)))     < 1e-9;
    c3 = abs(T_win(ii)/Ts_log(ii)  - round(T_win(ii)/Ts_log(ii)))  < 1e-6;
    c4 = abs(t_win0(ii)/Ts_log(ii) - round(t_win0(ii)/Ts_log(ii))) < 1e-6;
    c5 = f_p(ii) < 1/(2*Ts_log(ii));
    if ~(c1&&c2&&c3&&c4&&c5)
        ok = false;
        fprintf('  VERLETZT bei f_p = %.4f Hz: [%d %d %d %d %d]\n', f_p(ii),c1,c2,c3,c4,c5);
    end
end
assert(ok,'Kohaerenzbedingungen verletzt.');

t_sim_total = sum(2*(t_win0 + T_win));
fprintf('%d Frequenzpunkte, %d Simulationen\n', n_f, 2*n_f);
fprintf('laengstes Fenster: %.2f s   Einschwingen nach Injektion: %.2f ... %.2f s\n', ...
        max(T_win), min(t_ptb), max(t_ptb));
fprintf('groesster Logvektor: %d Samples\n', max(N_win));
fprintf('Summe simulierte Zeit: %.1f s  ->  %.2e Solverschritte\n', ...
        t_sim_total, t_sim_total/Ts);
fprintf('=================\n\n');

%% ---- 2) Simulationsschleife (hoechste Frequenz zuerst) ------------------
V_ph = zeros(2,2,n_f);   % V_ph(:,run,ii) = [v_q ; v_d]
I_ph = zeros(2,2,n_f);   % I_ph(:,run,ii) = [i_q ; i_d]
done = false(1,n_f);

fprintf('Start model ''%s'' ...\n', mdl);
t_load = tic;
load_system(mdl);
fprintf('  model loaded (%.1f s)\n\n', toc(t_load));

% fprintf('Testweise werden nur zwei Frequenzpunkte simuliert.')

t_all = tic;
for jj = 1:2 %--> Testweise werden nur zwei Frequenzpunkte simuliert
    ii = run_order(jj);                 % absteigende Frequenz

    T_win_run = T_win(ii);
    fp        = f_p(ii);
    w_ptb     = 2*pi*fp;
    StopT     = t_win0(ii) + T_win_run;

    fprintf('[%2d/%2d]  f_p = %10.4f Hz   T_win = %8.4f s   (k = %d, m = %d)\n', ...
            jj, n_f, fp, T_win_run, k_vec(ii), m_vec(ii));

    for run = 1:2
        phi_ptb  = phi_ptb_vec(run);
        phs0 = mod(-w_ptb*t_start, 2*pi);   % Sinus beginnt bei t_start auf null

        fprintf(['  Run %d/2  phi_ptb = %+6.1f deg: perturbation at %.6f Hz ' ...
                 '(omega_ptb = %.4f rad/s)\n' ...
                 '           amplitude %.2f %% of Vp (= %.4f V) = %.4f V, ' ...
                 'dv_q = %+.4f V, dv_d = %+.4f V\n' ...
                 '           injected at t_start = %.3f s, settling %.4f s, ' ...
                 'window %.4f s from %.4f s\n' ...
                 '           D_log = %d  ->  %d samples\n'], ...
                run, phi_ptb_deg(run), fp, w_ptb, ...
                100*amp_rel, Vp, amp_ptb, amp_ptb*cos(phi_ptb), amp_ptb*sin(phi_ptb), ...
                t_start, t_ptb(ii), T_win_run, t_win0(ii), D_vec(ii), N_win(ii));

        in = Simulink.SimulationInput(mdl);
        in = in.setModelParameter('StopTime',  num2str(StopT,16));
        in = in.setModelParameter('SolverType','Fixed-step');
        in = in.setModelParameter('Solver',    SOLVER);
        in = in.setModelParameter('FixedStep', num2str(Ts,16));

        D_log =     D_vec(ii);

        %in = in.setVariable('Ts',        Ts);
        %in = in.setVariable('t_start',   t_start);
        %in = in.setVariable('amp_ptb',   amp_ptb);
        %in = in.setVariable('omega_ptb', w_ptb);      % rad/s
        %in = in.setVariable('phi_ptb',   phi_ptb);        % rad
        %in = in.setVariable('phs0_ptb',  phs0);       % rad

        fprintf('           simulating, StopTime = %.4f s ... ', StopT);
        t_run = tic;
        out   = sim(in);
        fprintf('done in %.1f s\n', toc(t_run));
        if ~log_resolved
            nEl = out.logsout.numElements;
            ds_names = strings(1,nEl);
            for kk = 1:nEl
                ds_names(kk) = string(out.logsout{kk}.Name);
            end
            fprintf('\n  logsout: %d Signale -> %s\n', nEl, ...
                    strjoin(cellstr(ds_names), ', '));

            fn = fieldnames(LOG);
            for kk = 1:numel(fn)
                s   = LOG.(fn{kk});
                pos = find(ds_names == string(s.name));
                assert(~isempty(pos), ...
                    'Signal "%s" fehlt in logsout. Vorhanden: %s', ...
                    s.name, strjoin(cellstr(ds_names),', '));
                assert(isscalar(pos), ...
                    'Signal "%s" ist %d-fach in logsout - Namen eindeutig machen.', ...
                    s.name, numel(pos));
                nc = size(out.logsout{pos}.Values.Data, 2);
                assert(nc == s.ncol, ...
                    'Signal "%s": %d Spalten erwartet, %d gefunden.', ...
                    s.name, s.ncol, nc);
                LOG.(fn{kk}).numInLogdata = pos;
                fprintf('    %-14s -> logsout{%d}  (%d Spalten)\n', s.name, pos, nc);
            end
            fprintf('\n');
            log_resolved = true;
        end
        
        % Nächste Zuordnung Log-Daten
        [t_v_dat, v_qd0_PCC_dat] = getlog(out, LOG.v_qd0_PCC);
        [~,       i_qd0_PCC_dat] = getlog(out, LOG.i_qd0_PCC);
        %[~,       v_qd0_ptb_dat] = getlog(out, LOG.v_qd0_ptb);
        [~,       phi_clk_dat  ] = getlog(out, LOG.phi_clk);

        % Tatsaechliches Log-Raster aus den Daten bestimmen, nicht annehmen
        Ts_raw = median(diff(t_v_dat));
        step   = round(Ts_log(ii)/Ts_raw);          % 1 falls Simulink dezimiert
        assert(step >= 1 && abs(step*Ts_raw - Ts_log(ii)) < 1e-12, ...
               ['Log-Raster %.6e passt nicht zu Ts_log = %.6e ' ...
                '(Verhaeltnis %.4f).'], Ts_raw, Ts_log(ii), Ts_log(ii)/Ts_raw);

        n0  = round(t_win0(ii)/Ts_raw) + 1;
        idx = n0 : step : n0 + (N_win(ii)-1)*step;
        assert(idx(end) <= numel(t_v_dat), ...
               'Fenster laenger als die Simulation (f_p = %g Hz).', fp);

        tt  = (0:N_win(ii)-1).' * Ts_log(ii);
        % run ist einer der beiden Läufe je Frequenzpunkt; es wird mit zwei
        % linear unabhänigen Perturbationsspannungsphasoren angeregt
        % ii ist der Frequenzpunkt
        % : -> in Zeile 1 und 2 für Spalte run geht der Spaltenvektor [vq;vd] 
        V_ph(:,run,ii) = [ dft1(v_qd0_PCC_dat(idx,IQ),fp,tt) ; dft1(v_qd0_PCC_dat(idx,ID),fp,tt) ];
        I_ph(:,run,ii) = [ dft1(i_qd0_PCC_dat(idx,IQ),fp,tt) ; dft1(i_qd0_PCC_dat(idx,ID),fp,tt) ];
        
    end


    cV = cond(V_ph(:,:,ii));
    done(ii) = true;
    if cV > 2*cond_theo
        fprintf(2,'  -> cond(V) = %.3f   ACHTUNG: erwartet waren %.3f\n\n', cV, cond_theo);
    else
        fprintf('  -> cond(V) = %.3f (erwartet %.3f)   |v^1| = %.4f V   |v^2| = %.4f V\n\n', ...
                cV, cond_theo, norm(V_ph(:,1,ii)), norm(V_ph(:,2,ii)));
    end
end
fprintf('All simulations finished in %.1f min.\n\n', toc(t_all)/60);

%% ---- 3) Auswertung ------------------------------------------------------
idx_done = find(done);
n_d      = numel(idx_done);
assert(n_d > 0, 'Keine Frequenzpunkte gemessen.');
if n_d < n_f
    warning('a_01:teilmenge','Nur %d von %d Frequenzpunkten gemessen.', n_d, n_f);
end

Ydat  = nan(2,2,n_f);   Zdat  = nan(2,2,n_f);
condV = nan(1,n_f);     condI = nan(1,n_f);   res_Z = nan(1,n_f);

for ii = idx_done                    %                       run1    run2
    V = V_ph(:,:,ii);            % ii = Laeufe                [vq1] [vq2]
                                 %                            [vd1] [vd2]
    I = I_ph(:,:,ii);

    condV(ii) = cond(V);
    condI(ii) = cond(I);

    % V ist eingepraegt und damit gut konditioniert -> V invertieren, nicht I
    % x = A\B ist die Lösung für die Gleichung Ax = B
    % Y = I/V ist die Lösung für die Gleichung YV = I.

    Y = I / V;
    Z = inv(Y);

    % Gegenprobe: V/I ist mathematisch dasselbe. Abweichung = Konditionsproblem
    res_Z(ii) = norm(Z - V/I, 'fro') / norm(Z, 'fro');

    Ydat(:,:,ii) = Y;
    Zdat(:,:,ii) = Z;
end

[rmax, kmax] = max(res_Z(idx_done));
if rmax > 1e-8
    warning('a_01:kondition', ...
        'Groesstes Residuum inv(I/V) vs V/I: %.2e bei f_p = %.4f Hz (cond(I) = %.1f).', ...
        rmax, f_p(idx_done(kmax)), condI(idx_done(kmax)));
end

Yqq = squeeze(Ydat(1,1,:)).';  Yqd = squeeze(Ydat(1,2,:)).';
Ydq = squeeze(Ydat(2,1,:)).';  Ydd = squeeze(Ydat(2,2,:)).';

% FRD-Objekte nur aus den tatsaechlich gemessenen Punkten
f_sel = f_p(idx_done);
w_rad = 2*pi*f_sel(:);
Y_frd = frd(Ydat(:,:,idx_done), w_rad, 'FrequencyUnit','rad/s');
Z_frd = frd(Zdat(:,:,idx_done), w_rad, 'FrequencyUnit','rad/s');
Y_frd.OutputName = {'i_q','i_d'};  Y_frd.InputName = {'v_q','v_d'};
Z_frd.OutputName = {'v_q','v_d'};  Z_frd.InputName = {'i_q','i_d'};

% +/- Basis (qd-Reihenfolge):  u_+ = u_d + j*u_q,  u_- = u_d - j*u_q
Tpm = [1j 1; -1j 1];
Zpm = nan(2,2,n_f);
for ii = idx_done
    Zpm(:,:,ii) = Tpm * Zdat(:,:,ii) / Tpm;
end

cfg = struct( ...
    'version',     SCRIPT_VERSION, ...
    'run_stamp',   RUN_STAMP, ...
    'finished',    char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss')), ...
    'model',       mdl,         'solver',    SOLVER, ...
    'fg',          fg,          'Vp',        Vp, ...
    'Ts',          Ts,          'N_Ts',      N_Ts, ...
    't_start',     t_start,     't_ptb_min', t_ptb_min, ...
    't_ptb_max',   t_ptb_max,   'n_ptb_per', n_ptb_per, ...
    'amp_rel',     amp_rel,     'amp_ptb',   amp_ptb, ...
    'phi_ptb_deg', phi_ptb_deg, 'cond_theo', cond_theo, ...
    'grid_mode',   GRID_MODE, ...
    'f_lo',        f_lo,        'f_hi',      f_hi, ...
    'n_pts',       n_pts,       'M_min',     M_min,   'tol', tol, ...
    'SPP_min',     SPP_min,     'SPP_log',   SPP_log, ...
    'IQ',          IQ,          'ID',        ID, ...
    'n_measured',  n_d,         'n_total',   n_f);

save(SAVE_FILE,'cfg','f_p','done','idx_done','T_win','k_vec','m_vec','D_vec', ...
               'Ts','Ts_log','fg','amp_ptb','amp_rel','phi_ptb_deg', ...
               't_start','t_ptb','t_win0','V_ph','I_ph', ...
               'Ydat','Zdat','Zpm','condV','condI','res_Z', ...
               'Yqq','Yqd','Ydq','Ydd');
fprintf('Ergebnis gespeichert: %s  (%d von %d Punkten)\n', SAVE_FILE, n_d, n_f);

%% ---- 4) Plots -----------------------------------------------------------

% 4a) Anregungsvektoren: Betrag und Winkel je Lauf
figure('Name','Anregung','Color','w');
subplot(2,1,1)
loglog(f_sel, abs(squeeze(V_ph(1,1,:))),'o-','LineWidth',1.2); hold on
loglog(f_sel, abs(squeeze(V_ph(2,1,:))),'o--','LineWidth',1.0);
loglog(f_sel, abs(squeeze(V_ph(1,2,:))),'s--','LineWidth',1.0);
loglog(f_sel, abs(squeeze(V_ph(2,2,:))),'s-','LineWidth',1.2);
grid on; ylabel('|v|  [V]');
legend('|v_q| Lauf 1','|v_d| Lauf 1','|v_q| Lauf 2','|v_d| Lauf 2','Location','best');
title('Kanalspannungen je Lauf');
subplot(2,1,2)
ang1 = atan2(real(squeeze(V_ph(2,1,:))), real(squeeze(V_ph(1,1,:))))*180/pi;
ang2 = atan2(real(squeeze(V_ph(2,2,:))), real(squeeze(V_ph(1,2,:))))*180/pi;
semilogx(f_sel, ang1,'o-', f_sel, ang2,'s-','LineWidth',1.2); grid on; hold on
yline(phi_ptb_deg(1),'k:'); yline(phi_ptb_deg(2),'k:');
xlabel('f_sel  [Hz]'); ylabel('Anregungswinkel  [Grad]');
legend('Lauf 1','Lauf 2','Location','best');
title('gemessener Injektionswinkel gegen Sollwert');

% 4b) Konditionierung
figure('Name','Konditionierung','Color','w');
semilogx(f_sel, condV,'o-','LineWidth',1.3); grid on; hold on
yline(cond_theo,'g--','theoretisch');
yline(10,'r--','Warnschwelle');
xlabel('f_sel  [Hz]'); ylabel('cond(V)');
title(sprintf('Konditionszahl  (d\\phi = %.0f Grad)', dphi*180/pi));

% 4c) Z-Matrix
lbl = {'Z_{qq}','Z_{qd}';'Z_{dq}','Z_{dd}'};
figure('Name','Z_qd Betrag','Color','w');
for r=1:2, for c=1:2
    subplot(2,2,(r-1)*2+c);
    loglog(f_sel, abs(squeeze(Zdat(r,c,idx_done))),'o-','LineWidth',1.2); grid on
    xlabel('f  [Hz]'); ylabel('|Z|  [\Omega]'); title(lbl{r,c});
end, end

figure('Name','Z_qd Phase','Color','w');
for r=1:2, for c=1:2
    subplot(2,2,(r-1)*2+c);
    semilogx(f_sel, unwrap(angle(squeeze(Zdat(r,c,idx_done))))*180/pi,'o-','LineWidth',1.2);
    grid on; yline(0,'k:');
    xlabel('f  [Hz]'); ylabel('\angle Z  [Grad]'); title(lbl{r,c});
end, end

% 4d) Re{Z_qq}: PLL-induzierter Negativwiderstandsbereich (Wen et al. 2016)
figure('Name','Re{Z_qq}','Color','w');
semilogx(f_sel, real(squeeze(Zdat(1,1,:))),'o-','LineWidth',1.4); grid on; hold on
yline(0,'r--','LineWidth',1.2);
xlabel('f  [Hz]'); ylabel('Re\{Z_{qq}\}  [\Omega]');
title('Negativwiderstandsbereich');

%% ---- Hilfsfunktionen ----------------------------------------------------
function X = dft1(x, f, t)
% Einzelbin-DFT bei f. Der Vorfaktor kuerzt sich in Y = I/V heraus,
% solange er fuer Spannung und Strom identisch ist.
    X = mean( x(:) .* exp(-1j*2*pi*f*t(:)) );
end

function d = divisors_of(n)
    d = 1:n;  d = d(mod(n,d)==0);
end

function [t, d] = getlog(out, entry)
% Zeit und Daten eines per Signal Logging erfassten Signals.
    v = out.logsout{entry.numInLogdata}.Values;
    t = v.Time;
    d = v.Data;
end