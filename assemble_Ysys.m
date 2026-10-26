%% ========================================================================
%  assemble_Ysys.m
%  ------------------------------------------------------------------------
%  Hinweis (Stand: 20.10.2026): Beschreibungstext nicht redigiert
%
%  ------------------------------------------------------------------------
%
%
%  ------------------------------------------------------------------------
%  Schliesst die Schleife nach Gu et al., Gl. (18)-(21), OHNE Kron-Reduktion,
%  und regelt die Frame-Umrechnung (lokaler q'd'-Frame -> globaler D'Q').
%
%  Voraussetzung: build_Yb_13feeder.m wurde ausgefuehrt.
%
%  ------------------------------------------------------------------------
%  WARUM HIER DIE DUALE FORM (20)/(21) BENUTZT WIRD
%
%  Knoten 684 ist eine reine Verzweigung und traegt keine Maschine, also
%      Y_m684 = 0   ->   rang(Y_m) = 2N-2   ->   Z_m = Y_m^-1 EXISTIERT NICHT.
%  Die primale Form braucht Z_m und ist damit nicht auswertbar:
%      Gl.(18)  Zhat = Z_m (I + Y_b Z_m)^-1
%      Gl.(19)  Yhat = (I + Y_b Z_m)^-1 Y_b
%
%  Gu et al. geben fuer genau diesen Fall die duale Form an ("more
%  convenient when the machine admittance Y_m can be more readily obtained
%  than the machine impedance Z_m"):
%      Gl.(20)  Y_m = Z_m^-1 ,  Z_b = Y_b^-1
%      Gl.(21)  Yhat = Y_m (I + Z_b Y_m)^-1 ,  Zhat = (I + Z_b Y_m)^-1 Z_b
%
%  Die braucht Z_b statt Z_m. Und Y_b IST invertierbar, seit die beiden
%  Kondensatorbaenke als C_kk im Netz stehen: sie sind die einzigen
%  Querglieder in Y_b und machen es regulaer (cond(Y_b) = 1e2 ... 2e6).
%  Ohne sie waere auch Y_b singulaer und keine der beiden Formen ginge.
%
%  Zusammenhang, exakt und ohne jede Inversion von Y_m:
%      Zhat = (I + Z_b Y_m)^-1 Z_b = (Y_b + Y_m)^-1
%      Yhat = Y_m (I + Z_b Y_m)^-1 = Y_m (Y_b + Y_m)^-1 Y_b
%  Die zweite Form nutzt die Identitaet A(A+B)^-1 B = B(A+B)^-1 A, weshalb
%  Gl.(19) und Gl.(21) dasselbe Yhat liefern (numerisch bestaetigt).
%
%  STABILITAET
%      det(I + Y_b Z_m) = det(Y_b+Y_m)/det(Y_m)     (primale Schleife)
%      det(I + Z_b Y_m) = det(Y_b+Y_m)/det(Y_b)     (duale Schleife)
%  Beide haben dieselben Nullstellen - die von det(Y_b+Y_m). Sie
%  unterscheiden sich nur darin, welche offenen Pole herausdividiert
%  werden. Fuer das verallgemeinerte Nyquist-Kriterium ist das der
%  Unterschied in der Zaehlung N = Z - P: die duale Schleife L' = Z_b Y_m
%  hat die Pole von Z_b, die primale die von Z_m. Hier ist nur die duale
%  verfuegbar.
%
%  ------------------------------------------------------------------------
%  FRAME-KONVENTION (explizit)
%    Raumzeiger      x_s = 2/3 (x_a + a x_b + a^2 x_c),  a = exp(j 2pi/3)
%    Frame-Winkel    theta_n(t) = w0 t + xi_n
%    komplexe Groesse   u = x_s * exp(-j theta_n)
%    'dq': d-Achse auf theta, q eilt 90 Grad vor -> [x_d; x_q], T = [1 j; 1 -j]
%    'qd': q-Achse auf theta, d eilt 90 Grad nach -> [x_q; x_d], T = conj(...)
%    u ist in beiden Faellen dieselbe komplexe Groesse, daher sind Z_++ und
%    Z_-- konventionsunabhaengig; nur die reelle 2x2-Darstellung dreht die
%    Nebendiagonalen um.
%
%    Y_b, alle Lastbloecke, die Kondensatorbaenke und die DEG sind in der
%    +- -Basis DIAGONAL und damit invariant unter jeder Frame-Drehung.
%    Nur der VSC hat wegen der PLL Nebendiagonalen - nur dort wirkt xi_n.
% ========================================================================


%% ------------------------------------------------------------- 0) Optionen
o2 = struct();
o2.form     = 'auto';      % 'auto' | 'dual' (Gl.20/21) | 'primal' (Gl.18/19)
o2.vscConv  = 'qd';        % Konvention des VSC-Modells: 'qd' oder 'dq'
o2.xiSource = 'manual';  % 'loadflow' | 'manual' | 'zero'
o2.xiManual = [];       % der Simulation des Simulink models entnommen
o2.verbose  = true;

% VSC-Impedanzen: Handle auf eine REELLE 2x2-Matrix im LOKALEN Frame.
% Leer lassen, solange das Konvertermodell noch nicht steht.
%   Beispiel:  vsc.n634 = @(sv) [Zqq(sv) Zqd(sv); Zdq(sv) Zdd(sv)];
% Feldnamen mit 'n'-Praefix (MATLAB-Felder duerfen nicht mit Ziffer beginnen).
vsc = struct();


%% --------------------------------------------- 1) Vorbedingungen
fprintf('\n=== Schleifenformulierung nach Gu et al., Gl.(18)-(21) ===\n');
fprintf('-> Schritt 1: Vorbedingungen pruefen\n');
% Prüft ob die u.g. Variablen existieren
for v = {'Yb_num','Y_m','par','nodeName','nodeVb','fg','s','Rg','Lg'}
    assert(exist(v{1},'var')==1, '%s fehlt - bitte build_Yb_13feeder.m ausfuehren.', v{1});
end
% Prüft ob die u.g. Funktionen existieren
for fn = {'real2pm','pm2real','applyXi','blk2full','full2blk'}
    assert(exist(fn{1},'file')==2, '%s.m nicht im Pfad.', fn{1});
end
N = numel(nodeName);  w0n = 2*pi*fg;
if exist('U_rms','var')~=1, U_rms = 4160; end

% vscNodes ist ein cell array, welches die Knotennamen, der Knoten mit U=480 V enthält 
vscNodes = nodeName(nodeVb == 480);
% fülle vsc mit den um "n" erweiterten Knotennamen
for i = 1:numel(vscNodes)
    fn = ['n' vscNodes{i}];
    if ~isfield(vsc, fn), vsc.(fn) = []; end
end

% Parametervariablen in Y_m werden durch deren Werte ersetzt
Ym_num = subs(Y_m, par.sym, par.val);
% falls Yb noch nicht als Funktion s vorlag, wird eine erstellt

%% ------------------------------------------------------------------------
% Redaktion 20.10.2026: sollte die Matrix der Übertragsungsglieder nicht
%                      in Impedanzen vorliegen (vgl. Glg. (21))
%% ------------------------------------------------------------------------
if exist('Yb_fun','var')==1 && ~isempty(Yb_fun), Yb_h = Yb_fun;
else,                                            Yb_h = matlabFunction(Yb_num,'Vars',s);
end
Ym_h = matlabFunction(Ym_num, 'Vars', s);
fprintf('   %d Knoten -> Y_b ist %dx%d Bloecke (2x2), ausgeschrieben %dx%d\n', ...
        N, N, N, 2*N, 2*N);


%% --------------------------------------------- 2) Selbsttest der Frame-Tools
fprintf('-> Schritt 2: Selbsttest der Frame-Transformationen\n');
fprintf('\t\t\t Auskommentiert.\n');

% Hinweis: (Stand 20.10.2026) auskommentiert, nicht redigiert 

% Lt = 0.01;  g = 1/(w0n*Lt);
% Ypm_L = diag([1/(1i*w0n*Lt), 1/(-1i*w0n*Lt)]);
% e1 = max(max(abs(pm2real(Ypm_L,'dq') - [0 g; -g 0])));
% e2 = max(max(abs(pm2real(Ypm_L,'qd') - [0 -g; g 0])));
% e3 = max(max(abs(real2pm(pm2real(Ypm_L,'qd'),'qd') - Ypm_L)));
% Zt = [3+2i, 0.4-0.9i; 0.15+0.3i, 2.5-1.1i];
% e4 = max(max(abs(applyXi(Zt, pi/2) - [Zt(1,1) -Zt(1,2); -Zt(2,1) Zt(2,2)])));
% e5 = max(max(abs(applyXi(diag([1+2i,3-1i]), 0.37) - diag([1+2i,3-1i]))));
% e6 = max(max(abs(blk2full(full2blk(Yb_h(0))) - Yb_h(0))));
% fprintf('   dq %.1e | qd %.1e | Rueckweg %.1e | xi=90 %.1e | diag inv %.1e | blk %.1e\n', ...
%         e1,e2,e3,e4,e5,e6);
% assert(max([e1 e2 e3 e4 e5 e6]) < 1e-10, 'Frame-/Blockwerkzeuge inkonsistent.');


%% --------------------------------------------- 3) Lastfluss -> xi_n
fprintf('-> Schritt 3: Arbeitspunkt und xi_n bestimmen\n');
switch lower(o2.xiSource)
    case 'zero'
        xi = zeros(1,N);  V0 = [];
    case 'manual'
        assert(numel(o2.xiManual)==N, 'o2.xiManual braucht %d Werte.', N);
        xi = o2.xiManual(:).';  V0 = [];
    otherwise
        % Hinweis KI-generiert
        % Alle Lasten sind konstante Impedanz, die DEG ist eine Quelle hinter
        % Rg+jXg -> der 60-Hz-Arbeitspunkt ist LINEAR, keine Iteration noetig.
        % Gilt exakt, solange die VSC nicht einspeisen.
        Y0 = Yb_h(0) + Ym_h(0);
        Yp = Y0(1:2:end, 1:2:end);
        Iinj = zeros(N,1);
        k632 = find(strcmp(nodeName,'632'),1);
        Iinj(k632) = (U_rms/sqrt(3))/(Rg + 1i*w0n*Lg);
        V0 = Yp\Iinj;
        xi = angle(V0).';
        fprintf('   |V|/V_n : %s\n', strjoin(arrayfun(@(k) ...
            sprintf('%s=%.3f', nodeName{k}, abs(V0(k))/(4160/sqrt(3))), 1:N,'uni',0),'  '));
        fprintf('   xi [Grad]: %s\n', strjoin(arrayfun(@(k) ...
            sprintf('%s=%+.2f', nodeName{k}, rad2deg(xi(k))), 1:N,'uni',0),'  '));
        fprintf(['   Genauer waere der Arbeitspunkt aus powergui -> ' ...
                 '"Steady-State Voltages and Currents",\n   dann ' ...
                 'o2.xiSource = ''manual''.\n']);
end
% Eine globale Drehung aller xi_n ist eine Aehnlichkeitstransformation der
% Schleife und aendert det(I + Z_b Y_m) nicht - der Bezugswinkel ist frei.


%% --------------------------------------------- 4) VSC einhaengen
fprintf('-> Schritt 4: VSC-Beitraege aufbereiten\n');
vscIdx = [];  vscFun = {};
for i = 1:numel(vscNodes)
    k = find(strcmp(nodeName, vscNodes{i}),1);
    f = vsc.(['n' vscNodes{i}]);
    if isempty(k) || isempty(f), continue; end
    % lokaler reeller Frame -> +- -Basis -> globaler D'Q'-Frame -> Admittanz
    vscIdx(end+1) = k;                                                 %#ok<SAGROW>
    vscFun{end+1} = @(sv) inv(applyXi(real2pm(f(sv), o2.vscConv), xi(k))); %#ok<SAGROW>
end
if isempty(vscIdx)
    fprintf(['   Keine VSC-Modelle gesetzt. Sobald sie stehen:\n' ...
             '     vsc.n634 = @(sv) [Zqq Zqd; Zdq Zdd];   %% reell, lokaler Frame\n' ...
             '     -> real2pm(.,''qd'') -> applyXi(., xi(k)) -> inv(.)\n']);
else
    fprintf('   %d VSC an: %s\n', numel(vscIdx), strjoin(nodeName(vscIdx), ', '));
end
Ymtot_h = @(sv) addVsc(Ym_h(sv), sv, vscIdx, vscFun);


%% --------------------------------------------- 5) Form waehlen
fprintf('-> Schritt 5: Form waehlen (primal / dual)\n');
sv0 = 0.7 + 1i*31.4;
Ymt = Ymtot_h(sv0);  Ybt = Yb_h(sv0);
noMach = {};
for k = 1:N
    if max(max(abs(Ymt(2*k-1:2*k, 2*k-1:2*k)))) == 0, noMach{end+1} = nodeName{k}; end %#ok<SAGROW>
end
rYm = rank(Ymt);  rYb = rank(Ybt);
fprintf('   Knoten ohne Maschine : %s\n', ternary(isempty(noMach),'(keine)',strjoin(noMach,', ')));
fprintf('   rang(Y_m) = %d/%d  -> Z_m %s\n', rYm, 2*N, ternary(rYm==2*N,'existiert','existiert NICHT'));
fprintf('   rang(Y_b) = %d/%d  -> Z_b %s,  cond(Y_b) = %.2e\n', ...
        rYb, 2*N, ternary(rYb==2*N,'existiert','existiert NICHT'), cond(Ybt));

form = o2.form;
if strcmpi(form,'auto')
    if rYm == 2*N, form = 'primal'; else, form = 'dual'; end
end
assert(~(strcmpi(form,'primal') && rYm < 2*N), ...
    ['o2.form = ''primal'' verlangt ein invertierbares Y_m. Knoten ohne ' ...
     'Maschine: %s. Bitte o2.form = ''dual'' benutzen.'], strjoin(noMach,', '));
assert(~(strcmpi(form,'dual') && rYb < 2*N), ...
    'o2.form = ''dual'' verlangt ein invertierbares Y_b.');
fprintf('   -> Gl.%s wird benutzt\n', ternary(strcmpi(form,'dual'),'(20)/(21), dual','(18)/(19), primal'));


%% --------------------------------------------- 6) Gu-Groessen
fprintf('-> Schritt 6: Zhat, Yhat, Schleifenverstaerkung\n');

Ifull  = eye(2*N);
Zb_h   = @(sv) inv(Yb_h(sv));                               % Gl. (20)
Zhat_h = @(sv) inv(Yb_h(sv) + Ymtot_h(sv));                 % = Gl.(18) = Gl.(21)
Yhat_h = @(sv) Ymtot_h(sv)*Zhat_h(sv)*Yb_h(sv);             % = Gl.(19) = Gl.(21)

if strcmpi(form,'dual')
    L_h  = @(sv) Zb_h(sv)*Ymtot_h(sv);                      % duale Schleife
    lbl  = 'L = Z_b Y_m';
else
    L_h  = @(sv) Yb_h(sv)/Ymtot_h(sv);                      % = Y_b Z_m
    lbl  = 'L = Y_b Z_m';
end
RD_h  = @(sv) Ifull + L_h(sv);
nyq_h = @(sv) eig(L_h(sv));

% Blockweise lesen: Zhat_blk{k,k} ist die 2x2-Treiberimpedanz am Knoten k
Zhat_blk_h = @(sv) full2blk(Zhat_h(sv));


%% --------------------------------------------- 7) Verifikation
fprintf('-> Schritt 7: Verifikation\n');
sv = sv0;  Yb = Yb_h(sv);  Ym = Ymtot_h(sv);  Zb = Zb_h(sv);
Zsum   = inv(Yb + Ym);
Zhat21 = (Ifull + Zb*Ym)\Zb;                                % Gl. (21) woertlich
Yhat21 = Ym/(Ifull + Zb*Ym);                                % Gl. (21) woertlich
fprintf('   (a) Gl.(21) Zhat vs. (Y_b+Y_m)^-1     : %.3e (rel.)\n', ...
        max(max(abs(Zhat21-Zsum)))/max(max(abs(Zsum))));
fprintf('   (b) Gl.(21) Yhat vs. Y_m (Y_b+Y_m)^-1 Y_b : %.3e (rel.)\n', ...
        max(max(abs(Yhat21-Ym*Zsum*Yb)))/max(max(abs(Yhat21))));
fprintf('   (c) A(A+B)^-1 B = B(A+B)^-1 A          : %.3e (rel.)\n', ...
        max(max(abs(Yb*Zsum*Ym - Ym*Zsum*Yb)))/max(max(abs(Ym*Zsum*Yb))));
[~,l1] = logdet(RD_h(sv));
[~,l2] = logdet(Yb+Ym);
[~,l3] = logdet(ternary(strcmpi(form,'dual'), Yb, Ym));
fprintf('   (d) log|det(I+L)| vs. log|det(Y_b+Y_m)/det(%s)| : %.3e (rel.)\n', ...
        ternary(strcmpi(form,'dual'),'Y_b','Y_m'), abs(l1-(l2-l3))/max(abs(l1),eps));
fprintf('   (e) Blockweg blk2full(full2blk(.))     : %.3e\n', ...
        max(max(abs(blk2full(full2blk(Zsum))-Zsum))));

if o2.verbose
    fprintf('\nKnoten (Reihenfolge der Bloecke): %s\n', strjoin(nodeName,', '));
    fprintf('Handles : Zhat_h(s), Yhat_h(s), Zb_h(s), L_h(s) [%s], RD_h(s), nyq_h(s)\n', lbl);
    fprintf('Blockweise lesen : Zk = Zhat_blk_h(s); Zk{k,k}  -> 2x2 am Knoten k\n');
    fprintf('Stationaere Frequenz f -> s = 1i*2*pi*(f - %g)\n\n', fg);
end


%% ========================================================================
%  Lokale Funktionen
%% ========================================================================

function Y = addVsc(Y, sv, idxv, funv)
for i = 1:numel(idxv)
    k = idxv(i);
    Y(2*k-1:2*k, 2*k-1:2*k) = Y(2*k-1:2*k, 2*k-1:2*k) + funv{i}(sv);
end
end

function y = ternary(c,a,b)
if c, y = a; else, y = b; end
end

function [sg, ld] = logdet(A)
[~,U,P] = lu(A);
d  = diag(U);
sg = det(P)*prod(sign(d));
ld = sum(log(abs(d)));
end
