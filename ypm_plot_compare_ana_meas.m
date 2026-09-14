function out = ypm_aus_scan(sys, FILES, swap_pm, f_max_plot)
%YPM_AUS_SCAN  Gemessene Scans in die +/-Basis bringen und gegen ein Modell legen.
%
%   out = ypm_aus_scan(sys)                       % Dateien per Dialog waehlen
%   out = ypm_aus_scan(sys, FILES)
%   out = ypm_aus_scan(sys, FILES, true)          % +/- Kanaele tauschen
%   out = ypm_aus_scan(sys, FILES, false, 2000)   % Darstellung bis +/-2000 Hz
%
%   sys         2x2 tf in der +/-Basis, z.B. Y_st_tf
%   FILES       Dateiname oder Zellarray der Ergebnisdateien aus a_01_impedanzscan.m
%   swap_pm     optional, tauscht + und - der Messung (Konventionsabgleich)
%   f_max_plot  optional, Darstellungsbereich in Hz, Vorgabe 1000 -> [-1000, 1000]
%
%   Der Frequenzbereich des Modells folgt f_max_plot. Messpunkte ausserhalb
%   des Bereichs bleiben in out.scan(a) erhalten, werden aber nicht gezeichnet.
%
%   Schriftgroessen von Achsen, Legende und Knopf skalieren mit der
%   Fenstergroesse (SizeChangedFcn). Knopf "Vektor-Export..." fuer PDF/SVG/EPS/EMF.
%
%   Rueckgabe:
%     out.f_mod, out.Y_mod      analytische Kurve, 2x2xnw, Frequenz in Hz
%     out.scan(a).f, .Ypm       Messung je Datei in der +/-Basis
%     out.scan(a).offdiag       max|Nebendiag| / max|Diag| je Frequenzpunkt
%     out.fig                   Handle der Abbildung
%
%   Grundlage: Y_+- = T * Y_qd / T ist eine exakte Aehnlichkeitstransformation.
%   Diagonal wird das Ergebnis nur, wenn Y_qq = Y_dd und Y_qd = -Y_dq gilt.
%   Siehe Harnefors, IEEE TIE 54(4), 2007 und Rygg et al., IEEE JESTPE, 2016.

if nargin < 3 || isempty(swap_pm),    swap_pm    = false; end
if nargin < 4 || isempty(f_max_plot), f_max_plot = 1000;  end

wrange = {0, 2*pi*f_max_plot};

%% ---- 0) Analytisches Modell auswerten ---------------------------------
[mag, phase, wout] = bode(sys, wrange);
out.f_mod = wout(:)/(2*pi);
out.Y_mod = mag .* exp(1j*phase*pi/180);

lab_mod = 'Y_{analytisch}';

%% ---- 1) Dateien laden --------------------------------------------------
if nargin < 2 || isempty(FILES)
    [fn, fp_dir] = uigetfile('*.mat','Ergebnisdatei(en) waehlen','MultiSelect','on');
    if isequal(fn,0), error('Abgebrochen.'); end
    if ischar(fn), fn = {fn}; end
    FILES = fullfile(fp_dir, fn);
elseif ischar(FILES)
    FILES = {FILES};
end
nF = numel(FILES);

Tpm = [1j 1; -1j 1];          % wirkt auf [x_q; x_d]
out.Tpm = Tpm;

fprintf('=== Geladene Ergebnisse ===\n');
for a = 1:nF
    S = load(FILES{a});
    need = {'f_p','Zdat','Ydat','V_ph','I_ph'};
    miss = need(~isfield(S, need));
    assert(isempty(miss), 'In %s fehlen: %s', FILES{a}, strjoin(miss,', '));

    if isfield(S,'idx_done')
        sel = S.idx_done(:).';
    elseif isfield(S,'done')
        sel = find(S.done(:).');
    else
        sel = find(~isnan(squeeze(S.Zdat(1,1,:))).');
    end

    f = S.f_p(sel);
    Y = S.Ydat(:,:,sel);

    % aufsteigend sortieren, sonst zickzackt die Verbindungslinie
    [f, ord] = sort(f(:).');
    Y = Y(:,:,ord);

    % ---- in die +/-Basis drehen
    Ypm = nan(2,2,numel(sel));
    off = nan(1,numel(sel));
    for q = 1:numel(sel)
        M = Tpm * Y(:,:,q) / Tpm;
        if swap_pm, M = M([2 1],[2 1]); end
        Ypm(:,:,q) = M;
        off(q) = max(abs([M(1,2) M(2,1)])) / max(abs([M(1,1) M(2,2)]));
    end

    if isfield(S,'cfg')
        c = S.cfg;
        lab = sprintf('%s v%s | amp %.1f%%', getfield_or(c,'run_stamp','?'), ...
                      getfield_or(c,'version','?'), 100*getfield_or(c,'amp_rel',NaN));
    else
        [~,bn] = fileparts(FILES{a});
        lab = bn;
    end
    % Legende bleibt nackt; die Zuordnung Index -> Datei steht in der
    % Konsolenausgabe weiter unten.
    if nF > 1
        labTex = sprintf('Y_{Mess,%d}', a);
    else
        labTex = 'Y_{Mess}';
    end

    % ---- negativen Ast ergaenzen:  Y(-f) = J * conj(Y(f)) * J
    % Gilt exakt, weil die gemessene dq-Matrix reelle Koeffizienten hat.
    J  = [0 1; 1 0];
    nQ = numel(sel);
    f_all   = [-fliplr(f(:).'), f(:).'];
    Ypm_all = nan(2,2,2*nQ);
    for q = 1:nQ
        Ypm_all(:,:,nQ+q)   = Ypm(:,:,q);
        Ypm_all(:,:,nQ+1-q) = J * conj(Ypm(:,:,q)) * J;
    end

    % ---- Diagnose: Messung gegen Modell am gesamten Messraster
    Hmod = freqresp(sys, 2*pi*f(:).');
    rat  = squeeze(abs(Ypm(1,1,:))) ./ squeeze(abs(Hmod(1,1,:)));

    out.scan(a).file     = FILES{a};
    out.scan(a).label    = lab;        % roh, fuer Konsole
    out.scan(a).labelTex = labTex;     % fuer Legende
    out.scan(a).f        = f;
    out.scan(a).Y        = Y;
    out.scan(a).Ypm      = Ypm;
    out.scan(a).offdiag  = off;
    out.scan(a).f_all    = f_all;
    out.scan(a).Ypm_all  = Ypm_all;
    out.scan(a).ratio11  = rat(:).';

    fprintf('  %s = %s\n     %d Punkte, %.4f ... %.1f Hz, Nebendiagonale %.1f ... %.1f %%\n', ...
            strrep(labTex,'\_','_'), lab, numel(sel), min(f), max(f), ...
            100*min(off), 100*max(off));
    fprintf('     |Y_mess/Y_modell| auf (+,+): min %.3f, median %.3f, max %.3f\n', ...
            min(rat), median(rat), max(rat));
    if max(f) > f_max_plot
        fprintf('     Hinweis: %d Punkte oberhalb %.0f Hz werden nicht gezeichnet.\n', ...
                sum(f > f_max_plot), f_max_plot);
    end
end
fprintf('\n');

%% ---- 2) Vergleich zeichnen --------------------------------------------
col  = lines(max(nF,3));
% Geschweifte Klammern gruppieren nur, sie erscheinen nicht im Text. Noetig,
% weil der TeX-Interpreter '\Deltav' als ein unbekanntes Schluesselwort liest.
inN  = {'{\Delta}v_+','{\Delta}v_-'};
outN = {'{\Delta}i_+','{\Delta}i_-'};

fh = figure('Name','Y_{+-}: Modell und Messung','Color','w', ...
            'Renderer','painters','Position',[100 80 900 760]);
out.fig = fh;
hLeg  = gobjects(1, nF+1);          % [Modell, Messung 1..nF]
axAll = gobjects(1, 8);
k = 0;

mmod = out.f_mod <= f_max_plot;     % Modell nur im Darstellungsbereich

for r = 1:2
    for c = 1:2
        % ---------- Betrag ----------
        k = k+1;
        axAll(k) = subplot(4,2,(r-1)*4+c);
        ymod = squeeze(out.Y_mod(r,c,:));
        h0 = plot(out.f_mod(mmod), 20*log10(abs(ymod(mmod))), ...
                  'k-','LineWidth',1.2); hold on
        for a = 1:nF
            [fa, ya] = band_mit_luecke(out.scan(a).f_all, ...
                           squeeze(out.scan(a).Ypm_all(r,c,:)).', f_max_plot);
            ha = plot(fa, 20*log10(abs(ya)), '-o', 'Color',col(a,:), ...
                      'LineWidth',0.9, 'MarkerSize',4);
            if r == 1 && c == 1, hLeg(a+1) = ha; end
        end
        if r == 1 && c == 1, hLeg(1) = h0; end
        grid on; xlim([-f_max_plot f_max_plot]); set(gca,'XTickLabel',[]);
        if r == 1, title(inN{c}); end
        if c == 1, ylabel({outN{r},'|Y|  [dB]'}); end

        % ---------- Phase ----------
        k = k+1;
        axAll(k) = subplot(4,2,(r-1)*4+2+c);
        plot(out.f_mod(mmod), angle(ymod(mmod))*180/pi, 'k-','LineWidth',1.2); hold on
        for a = 1:nF
            [fa, ya] = band_mit_luecke(out.scan(a).f_all, ...
                           squeeze(out.scan(a).Ypm_all(r,c,:)).', f_max_plot);
            plot(fa, angle(ya)*180/pi, '-o', 'Color',col(a,:), ...
                 'LineWidth',0.9, 'MarkerSize',4);
        end
        grid on; xlim([-f_max_plot f_max_plot]);
        if c == 1, ylabel('\angle Y  [Grad]'); end
        if r == 2, xlabel('f  [Hz]'); else, set(gca,'XTickLabel',[]); end
    end
end

% ---- Achsen stauchen: unten Platz fuer Legende und Knopf, oben fuer "von"
for k = 1:numel(axAll)
    p = axAll(k).Position;
    axAll(k).Position = [p(1), 0.065 + 0.88*p(2), p(3), 0.88*p(4)];
end

labs = [{lab_mod}, {out.scan.labelTex}];
lgd  = legend(axAll(1), hLeg, labs, 'Orientation','horizontal', ...
              'NumColumns', min(nF+1, 3), 'Tag','lgd_haupt');

% ---- "von" und "nach" je einmal zentral, in einer unsichtbaren Vollbildachse.
% Bezug sind die tatsaechlichen Achsenpositionen, nicht feste Zahlen.
% PickableParts 'none', damit Zoom und Datentipps der Subplots erreichbar bleiben.
P  = cell2mat(get(axAll(:),'Position'));
x0 = min(P(:,1));  x1 = max(P(:,1)+P(:,3));
y0 = min(P(:,2));  y1 = max(P(:,2)+P(:,4));

axLbl = axes('Parent',fh, 'Units','normalized', 'Position',[0 0 1 1], ...
             'XLim',[0 1], 'YLim',[0 1], 'Visible','off', ...
             'HitTest','off', 'PickableParts','none', 'Tag','ax_labels');
text(axLbl, (x0+x1)/2, min(y1+0.055, 0.965), 'von', ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontWeight','normal', 'HitTest','off', 'Tag','txt_von');
text(axLbl, 0.022, (y0+y1)/2, 'nach', 'Rotation',90, ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontWeight','normal', 'HitTest','off', 'Tag','txt_nach');

%% ---- 3) Export-Bedienung ----------------------------------------------
vorschlag = sprintf('Ypm_vergleich_%.0fHz', f_max_plot);

uicontrol(fh, 'Style','pushbutton', 'String','Vektor-Export...', ...
          'Tag','btn_export', 'Units','normalized', ...
          'Position',[0.855 0.004 0.140 0.042], ...
          'Callback', @(~,~) vektor_export(fh, vorschlag));

m = uimenu(fh, 'Label','Export');
uimenu(m, 'Label','Vektorgrafik speichern...', ...
          'Callback', @(~,~) vektor_export(fh, vorschlag));

%% ---- 4) Mitskalieren bei Fenstergroessenaenderung ----------------------
set(fh, 'SizeChangedFcn', @(s,~) fig_skalieren(s));
fig_skalieren(fh);

end

%% ------------------------------------------------------------------------
function v = getfield_or(S, name, dflt)
if isstruct(S) && isfield(S, name), v = S.(name); else, v = dflt; end
end

function [x, y] = band_mit_luecke(x, y, fmax)
%BAND_MIT_LUECKE  Auf |f| <= fmax beschneiden, dann NaN an der Nahtstelle.
% Die Maske ist symmetrisch, die Haelften bleiben also gleich lang. Das NaN
% verhindert, dass die Verbindungslinie ueber den datenfreien Bereich um
% f = 0 gezogen wird.
x = x(:).';  y = y(:).';
m = abs(x) <= fmax;
x = x(m);  y = y(m);
n = numel(x)/2;
if mod(numel(x),2) == 0 && n >= 1
    x = [x(1:n), NaN, x(n+1:end)];
    y = [y(1:n), NaN, y(n+1:end)];
end
end

function fig_skalieren(fh)
%FIG_SKALIEREN  Schriftgroessen an die Fenstergroesse koppeln.
% Bezug ist die MATLAB-Standardfigur (560 x 420 px) mit 10 pt. Das drawnow
% zwingt die Legende, ihre Box nach der Schriftaenderung neu zu vermessen,
% bevor sie mittig gesetzt wird; der Wiedereintrittsschutz verhindert, dass
% das drawnow die SizeChangedFcn erneut ausloest.
if ~isgraphics(fh), return; end
if isappdata(fh,'inScale') && getappdata(fh,'inScale'), return; end
setappdata(fh,'inScale',true);
schutz = onCleanup(@() setappdata(fh,'inScale',false));   %#ok<NASGU>

u = get(fh,'Units');  set(fh,'Units','pixels');
p = get(fh,'Position');  set(fh,'Units',u);

s  = min(p(3)/560, p(4)/420);
s  = max(0.75, min(s, 3));          % begrenzen, sonst wird es absurd
fs = 10*s;

set(findobj(fh,'Type','axes'), 'FontSize', fs);   % Titel und Labels folgen ueber
                                                  % Title-/LabelFontSizeMultiplier
btn = findobj(fh,'Tag','btn_export');
if ~isempty(btn), set(btn, 'FontSize', max(8, 0.9*fs)); end

% "von" und "nach" folgen nicht automatisch, sie sind eigene text-Objekte
for tg = {'txt_von','txt_nach'}
    t = findobj(fh,'Tag',tg{1});
    if ~isempty(t), set(t, 'FontSize', 1.1*fs, 'FontWeight','normal'); end
end

lgd = findobj(fh,'Tag','lgd_haupt');
if ~isempty(lgd)
    set(lgd, 'FontSize', fs, 'Units','normalized');
    drawnow limitrate;                       % Box neu vermessen lassen
    lp = get(lgd,'Position');
    set(lgd,'Position',[0.5 - lp(3)/2, 0.005, lp(3), lp(4)]);
end
end

function vektor_export(fh, vorschlag)
%VEKTOR_EXPORT  Speicherdialog und vektorieller Export der Abbildung.
filt = { '*.pdf', 'PDF, vektoriell (*.pdf)'; ...
         '*.svg', 'SVG (*.svg)'; ...
         '*.eps', 'EPS Level 3, Farbe (*.eps)'; ...
         '*.emf', 'EMF, nur Windows (*.emf)' };

[fn, fp] = uiputfile(filt, 'Abbildung vektoriell exportieren', vorschlag);
if isequal(fn,0), return; end
ziel = fullfile(fp, fn);
[~,~,ext] = fileparts(ziel);
if isempty(ext), ext = '.pdf'; ziel = [ziel ext]; end

% Bedienelemente aus dem Export heraushalten; wird auch bei Fehler zurueckgesetzt
btn = findobj(fh, 'Tag','btn_export');
set(btn, 'Visible','off');
zurueck = onCleanup(@() set(btn, 'Visible','on'));   %#ok<NASGU>
drawnow;

set(fh, 'Renderer','painters');   % erzwingt Vektorausgabe statt Rasterung

switch lower(ext)
    case '.svg'
        % exportgraphics kann SVG erst in neueren Releases, print immer
        print(fh, '-dsvg', '-painters', ziel);

    case {'.pdf','.eps','.emf'}
        if exist('exportgraphics','file') == 2
            exportgraphics(fh, ziel, 'ContentType','vector', 'BackgroundColor','white');
        else
            set(fh, 'PaperPositionMode','auto');
            pp = get(fh, 'PaperPosition');
            set(fh, 'PaperSize', pp(3:4));
            switch lower(ext)
                case '.pdf', print(fh, '-dpdf',   '-painters', ziel);
                case '.eps', print(fh, '-depsc2', '-painters', ziel);
                case '.emf', print(fh, '-dmeta',  '-painters', ziel);
            end
        end

    otherwise
        error('Nicht unterstuetztes Format: %s', ext);
end

fprintf('Abbildung gespeichert: %s\n', ziel);
end