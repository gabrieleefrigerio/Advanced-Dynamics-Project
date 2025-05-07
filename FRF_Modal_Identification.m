function FRF_Modal_Identification()
clc
close all
clear 

FRF = []; f = []; FRF_range = []; f_range = []; locs = []; minAmpField = []; ax = []; locs_idx = []; peaks = []; locs_pos = [];
 % fieldf0 = []; fieldA = []; fieldRh = []; fieldRl = []; fieldXi = []; 
% ===IMPORT FRF ===
[filename, pathname] = uigetfile('*.mat', 'Seleziona un file MAT');
data = load(fullfile(pathname, filename)); % caricato dentro 'data'

% Estrai numero massimo accelerometri
    num_acc = size(data.frf, 2);

    % === CREA GUI PER SCELTA ACC ===
    fig = uifigure('Name', 'Seleziona Accelerometro', 'Position', [500 400 500 200]);

    % Etichetta centrata
    uilabel(fig, ...
        'Text', sprintf('Numero massimo accelerometri disponibili: %d', num_acc), ...
        'Position', [100 140 300 30], ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    % Dropdown centrato
    dd = uidropdown(fig, ...
        'Items', arrayfun(@(i) sprintf('Accelerometro %d', i), 1:num_acc, 'UniformOutput', false), ...
        'Position', [100 90 300 30]);

    % Bottone per confermare selezione
    uibutton(fig, ...
        'Text', 'Conferma', ...
        'Position', [200 40 100 30], ...
        'ButtonPushedFcn', @(btn, ~) confermaCallback());





%% GUI PER SCEGLIERE SE SELEZIONARE RANGE
% choice = questdlg('Vuoi selezionare un range specifico della FRF?', ...
%     'Selezione range', 'Sì','No','No');
% 
% if strcmp(choice, 'Sì')
%     figure('Name', 'Seleziona Range FRF', 'Position', [100 100 1600 900]);
%     semilogy(f, abs(FRF)); grid on; grid minor;
%     title('Seleziona due punti per il range');
%     [xsel, ~] = ginput(2); close;
%     fmin = min(xsel); fmax = max(xsel);
% else
%     fmax = max(f);
%     if min(f) == 0
%         fmin = realmin;
%     else
%         fmin = min(f); 
%     end
% end
% 
% idx = f >= fmin & f <= fmax;
% f_range = f(idx);
% FRF_range = FRF(idx);
% 
% %% GUI PER SELEZIONARE PICCHI
% fig = uifigure('Name', 'Peak Finder', 'Position', [100 100 900 550]);
% 
% uilabel(fig, 'Text', 'Ampiezza minima del picco:', 'Position', [30, 430, 200, 22]);
% defaultMinAmp = round(0.15 * max(abs(FRF_range)));
% minAmpField = uieditfield(fig, 'numeric', 'Value', defaultMinAmp, 'Position', [30, 400, 150, 22]);
% 
% btnFind = uibutton(fig, 'Text', 'Find Peaks', 'Position', [30, 360, 150, 30]);
% btnProceedToFit = uibutton(fig, 'Text', 'Proceed', 'Position', [30, 310, 150, 30]);
% 
% ax = uiaxes(fig, 'Position', [200, 50, 670, 450]);
% title(ax, 'FRF con Picchi');
% xlabel(ax, 'Frequenza (Hz)');
% ylabel(ax, '|FRF|');
% grid(ax, 'on'); grid(ax, 'minor');
% 
% peaks = []; locs = []; locs_pos = [];
% 
% btnFind.ButtonPushedFcn = @(~,~) findPeaksCallback();
% btnProceedToFit.ButtonPushedFcn = @(~,~) startFitting();
% 
% % Trova picchi appena parte
% findPeaksCallback();

%% FUNZIONE TASTO PER TORVARE I PICCHI CON NUOVO VALORE MINIMO DEL PICCO
    function confermaCallback()
    % Ottieni indice dell'accelerometro selezionato
    selectedIndex = dd.Value;
    accIndex = sscanf(selectedIndex, 'Accelerometro %d');

    % Chiudi la finestra di selezione accelerometro
    close(fig);

    % Estrai FRF dell'accelerometro selezionato
    FRF = data.frf(:, accIndex);
    f = data.freq;

    % Salva FRF e f in base workspace per accesso nel resto dello script
    assignin('base', 'FRF', FRF);
    assignin('base', 'f', f);

  
% Avvia selezione range (emulando la logica successiva già presente)
    choice = questdlg('Vuoi selezionare un range specifico della FRF?', ...
        'Selezione range', 'Sì','No','No');

    if strcmp(choice, 'Sì')
        figure('Name', 'Seleziona Range FRF', 'Position', [100 100 1600 900]);
        semilogy(f, abs(FRF)); grid on; grid minor;
        title('Seleziona due punti per il range');
        [xsel, ~] = ginput(2); close;
        fmin = min(xsel); fmax = max(xsel);
    else
        fmax = max(f);
        if min(f) == 0
            fmin = realmin;
        else
            fmin = min(f);
        end
    end

    % Filtra i dati per il range selezionato
    idx = f >= fmin & f <= fmax;
    f_range = f(idx);
    FRF_range = FRF(idx);

    % Salva anche questi nel workspace base per uso successivo
    assignin('base', 'f_range', f_range);
    assignin('base', 'FRF_range', FRF_range);

    % Avvia la GUI per la selezione dei picchi
    runPeakSelectionGUI();
end

    function runPeakSelectionGUI()
    fig = uifigure('Name', 'Peak Finder', 'Position', [100 100 900 550]);

    uilabel(fig, 'Text', 'Ampiezza minima del picco:', 'Position', [30, 430, 200, 22]);
    defaultMinAmp = round(0.15 * max(abs(FRF_range)));
    minAmpField = uieditfield(fig, 'numeric', 'Value', defaultMinAmp, 'Position', [30, 400, 150, 22]);

    btnFind = uibutton(fig, 'Text', 'Find Peaks', 'Position', [30, 360, 150, 30]);
    btnProceedToFit = uibutton(fig, 'Text', 'Proceed', 'Position', [30, 310, 150, 30]);

    ax = uiaxes(fig, 'Position', [200, 50, 670, 450]);
    title(ax, 'FRF con Picchi');
    xlabel(ax, 'Frequenza (Hz)');
    ylabel(ax, '|FRF|');
    grid(ax, 'on'); grid(ax, 'minor');


    btnFind.ButtonPushedFcn = @(~,~) findPeaksCallback();
    btnProceedToFit.ButtonPushedFcn = @(~,~) startFitting();

    findPeaksCallback(); % Trova picchi inizialmente
    end
    
    
    function findPeaksCallback()

        % aggiorno il valore minimo del picco preso dalla gui
        minAmp = minAmpField.Value;

        % aggionro il plot con i nuovi picchi trovati
        cla(ax);
        semilogy(ax, f_range, abs(FRF_range)); hold(ax, 'on');
        [pks, locsFound] = findpeaks(abs(FRF_range), f_range, 'MinPeakHeight', minAmp);
        [~, locsPos] = findpeaks(abs(FRF_range), 'MinPeakHeight', minAmp);
        semilogy(ax, locsFound, pks, 'ro', 'MarkerFaceColor', 'r');
        peaks = pks;
        locs = locsFound;
        locs_pos = locsPos;
        legend(ax, 'FRF', 'Picchi trovati');
        hold(ax, 'off');
    end


    function startFitting()
        % verifico che ci sia almeno un picco trovato oppure termino
        % l'applicazione
        if isempty(peaks)
            uialert(fig, 'Trova prima i picchi!', 'Errore');
            return;
        end
        % chiudo la figure per trovare i picchi
        close(fig);
        % faccio partire la gui per ottimizzare i singoli picchi
        fitSingleMode(1, []);
    end

%% FUNZIONE PER OTTIMIZZARE I SINGOLI MODI
    function fitSingleMode(modeIndex, allParams)
        % se ho ottimizzato tutti i picchi apro la gui con il plot finale
        % di confronto
        if modeIndex > length(peaks)
            showFinalResult(allParams);
            return;
        end

        % creo la gui per ottimizzare i singoli picchi
        figFit = uifigure('Name', ['Modo ', num2str(modeIndex)], 'Position', [100 100 900 600]);

        % come guess delle frequenze proprie uso la posizione dei picchi
        % ricavati prima
        f0 = locs(modeIndex);
        [~, locs_idx] = min(abs(f_range - f0));

        % Half-power method per avere una guess del damping adimensionale
        mag_target = abs(FRF_range(locs_idx))/sqrt(2);
        left_idx = find(abs(FRF_range(1:locs_idx)) <= mag_target, 1, 'last');
        right_idx = find(abs(FRF_range(locs_idx:end)) <= mag_target, 1, 'first') + locs_idx - 1;


if isempty(left_idx) || isempty(right_idx)
            xi0 = 0.01;
        else
            f1 = f_range(left_idx);
            f2 = f_range(right_idx);
            xi0 = (f2 - f1) / (2*f0);
        end

        % definisco le altre guess da cui partire
        omega0 = 2*pi*f0;
        A0 = real( FRF_range(locs_idx) * (2i * xi0 * omega0^2) );

        Rh0 = 0;
        Rl0 = 0;

        % Tasti e caselle editabili GUI
        uilabel(figFit, 'Text', 'f0 (Hz):', 'Position', [30 450 100 22]);
        fieldf0 = uieditfield(figFit, 'numeric', 'Value', f0, 'Position', [150 450 100 22]);

        uilabel(figFit, 'Text', 'xi:', 'Position', [30 410 100 22]);
        fieldXi = uieditfield(figFit, 'numeric', 'Value', xi0, 'Position', [150 410 100 22]);

        uilabel(figFit, 'Text', 'A:', 'Position', [30 370 100 22]);
        fieldA = uieditfield(figFit, 'numeric', 'Value', A0, 'Position', [150 370 100 22]);


        uilabel(figFit, 'Text', 'Rh:', 'Position', [30 330 100 22]);
        fieldRh = uieditfield(figFit, 'numeric', 'Value', Rh0, 'Position', [150 330 100 22]);

        uilabel(figFit, 'Text', 'Rl:', 'Position', [30 290 100 22]);
        fieldRl = uieditfield(figFit, 'numeric', 'Value', Rl0, 'Position', [150 290 100 22]);

        btnOptimize = uibutton(figFit, 'Text', 'Ottimizza', ...
            'Position', [30 200 220 30], ...
            'ButtonPushedFcn', @(~,~) runOptimization());

        btnNext = uibutton(figFit, 'Text', 'Avanti', ...
            'Position', [30 160 220 30], ...
            'ButtonPushedFcn', @(~,~) proceedToNext());


        % Assi per ampiezza (in alto)
        axAmp = uiaxes(figFit, 'Position', [300, 350, 580, 250]);
        title(axAmp, 'FRF - Ampiezza');
        xlim(axAmp,[f0 - 50 f0 + 50]) % plotto il grafico solo nell'intorno del picco
        xlabel(axAmp, 'Frequenza (Hz)');
        ylabel(axAmp, '|FRF|');
        grid(axAmp, 'on'); grid(axAmp, 'minor');

        % Assi per fase (in basso)
        axPhase = uiaxes(figFit, 'Position', [300, 50, 580, 250]);
        title(axPhase, 'FRF - Fase');
        xlim(axPhase,[f0 - 50  f0 + 50]) % plotto il grafico solo nell'intorno del picco
        xlabel(axPhase, 'Frequenza (Hz)');
        ylabel(axPhase, 'Fase [rad]');
        grid(axPhase, 'on'); grid(axPhase, 'minor');

        optimized = [];

        % Appena apro, lancio subito l'ottimizzazione
        runOptimization();

function runOptimization()
            % vettore dei parametri da ottimizzare
            p0 = [fieldf0.Value*2*pi, fieldXi.Value, fieldA.Value, fieldRh.Value, fieldRl.Value];
            
            % === Selezione dinamica dell'intervallo attorno al picco ===
            
            % Frequenza centrale del picco corrente
            f_central = locs(modeIndex);
            
            % Calcolo della larghezza della finestra: metà distanza tra picchi adiacenti
            if modeIndex == 1
                % Primo picco: guarda solo verso il prossimo
                df = (locs(2) - locs(1)) / 15;
            elseif modeIndex == length(locs)
                % Ultimo picco: guarda solo verso il precedente
                df = (locs(end) - locs(end-1)) / 15;
            else
                % Picchi centrali: usa la media tra le due distanze adiacenti
                df1 = locs(modeIndex) - locs(modeIndex - 1);
                df2 = locs(modeIndex + 1) - locs(modeIndex);
                df = min(df1, df2) / 15; % più conservativo
            end
            
            % Estendi di un piccolo fattore (es. 20%) per sicurezza
            f_low = f_central -  df;
            f_high = f_central +  df;
            
            % Trova gli indici nell'intervallo
            idx_min = find(f_range >= f_low, 1, 'first');
            idx_max = find(f_range <= f_high, 1, 'last');
            
            % Fallback agli estremi se gli indici non vengono trovati
            if isempty(idx_min), idx_min = 1; end
            if isempty(idx_max), idx_max = length(f_range); end
            
            % Vettori finali da usare
            omega_vec = 2 * pi * f_range(idx_min:idx_max);
            G_exp = FRF_range(idx_min:idx_max).';


            % funzione di trasferimento numerica in formato anonymous
            modelFun = @(p, omega_vec) p(3)./ (-omega_vec.^2 + 2j*p(2)*p(1).*omega_vec + p(1)^2) + p(4) + ( p(5)./omega_vec.^2);
            % funzione con parametri scalati
            scale = [1, 1, 1, 1e-3, 1e-6];
            modelFun_scaled = @(p, omega_vec) modelFun(p .* scale, omega_vec);
            
            % cost function da minimizzare
            residui = @(p) sum( real(G_exp - modelFun_scaled(p, omega_vec)).^2 + imag(G_exp - modelFun_scaled(p, omega_vec)).^2  );
            
            % setting lsqnonlin
            opts = optimoptions('lsqnonlin','Display','off');

            % effettuo l'ottimizzazione
            [popt_scaled, ~] = lsqnonlin(residui, p0, [], [], opts);

            popt = popt_scaled .* scale;
            
            % salvo i valori ottimizzati
            optimized = popt;

            % calcolo la funzione di trasf numerica con i parmaetri
            % ottimizzati
            G_fit = modelFun(popt, omega_vec);

            % Aggiorno il plot
            % ---- Ampiezza ----
            cla(axAmp);
            semilogy(axAmp, f_range(idx_min:idx_max), abs(G_exp), 'b'); hold(axAmp, 'on'); grid(axAmp, 'minor');% grid(axAmp, 'on');
            xlim(axAmp,[fieldf0.Value-50 fieldf0.Value+50]) % plotto il grafico solo nell'intorno del picco
            semilogy(axAmp, f_range(idx_min:idx_max), abs(G_fit), 'r--', 'LineWidth', 1.5);
            legend(axAmp, 'FRF', 'Fit');
            grid(axAmp, 'minor'); %grid(axAmp, 'on');
            hold(axAmp, 'off');

            % ---- Fase ----
            cla(axPhase);
            plot(axPhase, f_range(idx_min:idx_max), angle(G_exp), 'b'); hold(axPhase, 'on');  grid(axPhase, 'minor'); % grid(axPhase, 'on');
            xlim(axPhase,[fieldf0.Value-50 fieldf0.Value+50]) % plotto il grafico solo nell'intorno del picco
            plot(axPhase, f_range(idx_min:idx_max), angle(G_fit), 'r--', 'LineWidth', 1.5);
            legend(axPhase, 'FRF', 'Fit');
            grid(axPhase, 'minor'); %grid(axPhase, 'on');
            hold(axPhase, 'off');

            % Aggiorna i campi nella gui
            fieldf0.Value = popt(1)/2/pi;
            fieldXi.Value = popt(2);
            fieldA.Value = popt(3);
            fieldRh.Value = popt(4);
            fieldRl.Value = popt(5);
        end

%% FUNZIONE TASTO PER PROCEDERE AL OTTIMIZZAZIONE DOPO
        function proceedToNext()
            if isempty(optimized)
                runOptimization();
            end
            % chiudo questa figure
            close(figFit);

            % apro la figure per il modo seguente
            fitSingleMode(modeIndex + 1, [allParams; optimized]);
        end
    end

%% GUI PER IL PLOT FINALE
    function showFinalResult(allParams)
        figFinal = figure('Name', 'FRF SDOF e Avvio MultiDOF', 'Position', [100 100 1600 900]);

        % Variabili locali
        omega_vec = 2*pi*f_range;
        G_total = zeros(size(omega_vec));
        modes = cell(size(allParams,1), 6);

        % Calcolo G_total e salvataggio parametri modali
        for i = 1:size(allParams,1)
            p = allParams(i,:);
            A = p(3);
            G_total = G_total + A ./ (-omega_vec.^2 + 2j*p(2)*p(1)*omega_vec + p(1)^2);
            modes(i,:) = {p(1)/(2*pi), p(2), p(1), p(3), p(4), p(5)};
        end
        ottimizzaResidui()
        G_total =  G_total + allParams(size(allParams,1), 4) + allParams(size(allParams,1), 5)./omega_vec.^2; %sommo i residui dell'ultimo modo
        modeTable = cell2table(modes, 'VariableNames', {'Frequenza (Hz)', 'Damping (xi)', 'omega0', 'A', 'Rh', 'Rl'});
            
        % Tabella a sinistra (20% larghezza)
        uit = uitable(figFinal, ...
            'Units', 'normalized', ...
            'Position', [0.05 0.15 0.20 0.75], ...
            'Data', modeTable{:,:}, ...
            'ColumnName', modeTable.Properties.VariableNames);

        % Grafico ampiezza (centrato in alto a destra)
        axAmplitude = axes(figFinal, ...
            'Units', 'normalized', ...
            'Position', [0.35 0.55 0.6 0.4]);
        semilogy(axAmplitude, f_range, abs(FRF_range), 'b'); hold on;
        semilogy(axAmplitude, f_range, abs(G_total), 'r--', 'LineWidth', 1.5);
        xlabel(axAmplitude, 'Frequenza (Hz)'); ylabel(axAmplitude, 'Ampiezza');
        title(axAmplitude, 'Ampiezza');
        legend(axAmplitude, 'FRF', 'Somma SDOF'); grid(axAmplitude, 'on'); grid(axAmplitude, 'minor');

        % Grafico fase (sotto il grafico ampiezza)
        axPhase = axes(figFinal, ...
            'Units', 'normalized', ...
            'Position', [0.35 0.15 0.6 0.3]);
        plot(axPhase, f_range, angle(FRF_range), 'b'); hold on;
        plot(axPhase, f_range, angle(G_total), 'r--', 'LineWidth', 1.5);
        xlabel(axPhase, 'Frequenza (Hz)'); ylabel(axPhase, 'Fase (rad)');
        title(axPhase, 'Fase');
        legend(axPhase, 'FRF', 'Somma SDOF'); grid(axPhase, 'on'); grid(axPhase, 'minor');

        % Pulsante per l'ottimizzazione MultiDOF
        uicontrol(figFinal, 'Style', 'pushbutton', 'String', 'Ottimizza MultiDOF', ...
            'Position', [460 10 200 30], ...
            'Callback', @(~,~) multiDOF_GUI(allParams, []));

        % Pulsante per esportare i dati
        uicontrol(figFinal, 'Style', 'pushbutton', 'String', 'Esporta Dati', ...
            'Position', [700 10 200 30], ...
            'Callback', @(~,~) exportData());

        % Pulsante Quit per uscire dal programma
        uicontrol(figFinal, 'Style', 'pushbutton', 'String', 'Quit', ...
            'Position', [940 10 200 30], ...
            'Callback', @(~,~) close(figFinal));

%% FUNZIONE OTTIMIZZAZIONE RESIDUI FINALE
   function ottimizzaResidui()
        omega_vec = 2*pi*f_range;

        % Parametri modali (senza residui)
        G_base = zeros(size(omega_vec));
        for i = 1:size(allParams,1)
            p = allParams(i,:);
            A = p(3);
            G_base = G_base + A ./ (-omega_vec.^2 + 2j*p(2)*p(1)*omega_vec + p(1)^2);
        end
         G_base = G_base.';
        % Funzione obiettivo per ottimizzazione Rh e Rl
        r0 = allParams(end,4:5);
        opts = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun',1e-12, 'TolX',1e-12);
        obj = @(res) sum(real(G_base + res(1) + res(2)./omega_vec.^2 - FRF_range).^2 + imag(G_base + res(1) + res(2)./omega_vec.^2 - FRF_range).^2);
        % Ottimizzazione (puoi usare fmincon se vuoi vincoli)
        res_ottimizzati = lsqnonlin(obj, r0, [], [], opts);
        % Aggiorna i parametri
        allParams(end,4) = res_ottimizzati(1);
        allParams(end,5) = res_ottimizzati(2);
        modes(end,5:6) = {res_ottimizzati(1), res_ottimizzati(2)}
    end
        %% FUNZIONE TASTO EXPORT
        function exportData()

    % Crea cartella Results se non esiste
    resultsDir = 'Results';
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end

    % === Salva FRF (ampiezza e fase) ===
    frfData = table(f_range(:), FRF_range(:), ...
        'VariableNames', {'Frequenza_Hz', 'FRF'});

    % Percorsi file FRF
    frfCSVname = fullfile(resultsDir, 'FRF_SDOF_Optimize.csv');
    frfMATname = fullfile(resultsDir, 'FRF_SDOF_Optimize.mat');

    % Salva CSV e MAT
    writetable(frfData, frfCSVname);
    save(frfMATname, 'frfData');

    % === Salva Tabella dei Modi ===
    modeTable = cell2table(modes, ...
        'VariableNames', {'Frequenza_Hz', 'Damping_xi', 'omega0', 'A', 'Rh', 'Rl'});

    % Percorsi file Modi
    modesCSVname = fullfile(resultsDir, 'ModeTable_FRF.csv');
    modesMATname = fullfile(resultsDir, 'ModeTable_FRF.mat');


% Salva CSV e MAT
    writetable(modeTable, modesCSVname);
    save(modesMATname, 'modeTable');

    % Conferma
    msgbox({'Dati esportati con successo!'}, 'Esportazione completata');

end


    end

    function multiDOF_GUI(initialParams, G_total_SDOF)
        figMDOF = uifigure('Name', 'Ottimizzazione MultiDOF', 'Position', [100 100 1000 550]);

        numModes = size(initialParams, 1);
        paramLabels = {'omega0', 'xi', 'A'};
        numParams = numel(paramLabels);

        fields = gobjects(numModes, numParams);
        spacingX = 200;
        startX = 50;
        startY = 500;
        spacingY = 70;

        % Crea i campi per ogni parametro
        for i = 1:numModes
            y0 = startY - (i - 1) * spacingY;
            for j = 1:numParams
                x0 = startX + (j - 1) * spacingX;
                uilabel(figMDOF, 'Text', sprintf(paramLabels{j}), ...
                    'Position', [x0 y0 80 22]);
                fields(i, j) = uieditfield(figMDOF, 'numeric', ...
                    'Value', initialParams(i, j), ...
                    'Position', [x0 + 85 y0 80 22]);
            end
        end

        % Bottone ottimizzazione centrato
        btnOptimize = uibutton(figMDOF, 'Text', 'Ottimizza MultiDOF', ...
            'Position', [figMDOF.Position(3)/2 - 75, 20, 150, 30], ...
            'ButtonPushedFcn', @(~, ~) runMultiDOF(G_total_SDOF));

        % Esegui subito l’ottimizzazione al lancio
        runMultiDOF(G_total_SDOF);

        function runMultiDOF(G_total_SDOF)
            p0 = zeros(numModes * numParams, 1);
            for i = 1:numModes
                for j = 1:numParams
                    p0(numParams * (i - 1) + j) = fields(i, j).Value;
                end
            end

            omega_vec = 2 * pi * f_range;
            G_exp = FRF_range;

            modelFun = @(p, w) sum( ...
                (p(3:4:end) + 1i * p(4:4:end)) ./ ...
                (-w.^2 + 2j * p(2:4:end) .* p(1:4:end) .* w + (p(1:4:end)).^2), 2);

            residui = @(p) [real(G_exp - modelFun(p, omega_vec)); ...
                imag(G_exp - modelFun(p, omega_vec))];

            opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 1e5);
            [popt, ~] = lsqnonlin(residui, p0, [], [], opts);

            showFinalMultiPlot(popt, G_total_SDOF);
        end

        function showFinalMultiPlot(popt, G_total_SDOF)
            figCompare = figure('Name', 'FRF Completo', 'Position', [100 100 1200 500]);

            omega_vec = 2 * pi * f_range;
            G_total_MDOF = sum( ...
                (popt(3:4:end) + 1i * popt(4:4:end)) ./ ...
                (-omega_vec.^2 + 2j * popt(2:4:end) .* popt(1:4:end) .* omega_vec + (popt(1:4:end)).^2), 2);

            subplot(1, 2, 1)
            plot(f_range, abs(FRF_range), 'k', 'DisplayName', 'FRF'); hold on;
            plot(f_range, abs(G_total_SDOF), 'b--', 'DisplayName', 'SDOF');
            plot(f_range, abs(G_total_MDOF), 'r-.', 'DisplayName', 'MDOF');
            xlabel('Frequenza (Hz)'); ylabel('|FRF|');
            title('Modulo FRF');
            legend; grid on; grid minor;

            subplot(1, 2, 2)
            plot(f_range, angle(FRF_range), 'k', 'DisplayName', 'FRF'); hold on;
            plot(f_range, angle(G_total_SDOF), 'b--', 'DisplayName', 'SDOF');
            plot(f_range, angle(G_total_MDOF), 'r-.', 'DisplayName', 'MDOF');
            xlabel('Frequenza (Hz)'); ylabel('Fase (rad)');
            title('Fase FRF');
            legend; grid on; grid minor;

            % Bottone per esportare
            uicontrol('Style', 'pushbutton', 'String', 'Esporta dati', ...
                'Position', [1050 20 100 30], ...
                'Callback', @(~, ~) exportResults(popt, G_total_SDOF, G_total_MDOF));
        end

        function exportResults(popt, G_total_SDOF, G_total_MDOF)
            [file, path] = uiputfile({'*.mat'; '*.csv'}, 'Salva risultati');
            if isequal(file, 0), return; end


omega_vec = 2 * pi * f_range;
            FRFdata = table(f_range.', abs(FRF_range).', abs(G_total_SDOF).', abs(G_total_MDOF).', ...
                'VariableNames', {'Freq', 'FRF', 'SDOF', 'MDOF'});

            if endsWith(file, '.mat')
                save(fullfile(path, file), 'FRFdata', 'popt');
            else
                writetable(FRFdata, fullfile(path, file));
            end
            msgbox('Dati esportati con successo!');
        end
    end

end

