% === Inizializzazione ===
clc;  % Pulisce la finestra dei comandi
clear;  % Rimuove tutte le variabili dalla memoria
close all;  % Chiude tutte le figure

% Definisce la frequenza angolare (omega) per la simulazione
load("Results\Analytic\freq.mat");
freq = freq';
omega = 2 * pi * freq;  % Frequenze da 0 a 200 Hz

% carico la matirce dei modi analitici
load("Results\Analytic\PhiMatrixco-located_1.20m.mat");

% creo il vettore x che discretizza la lungheza della trave
L = 1.2; %lunghezza trave
x = linspace(0,L, length(modes_shapes));



% Seleziona la cartella contenente i file .mat e .txt
cartella = uigetdir('Seleziona la cartella contenente i file .mat e .txt');

% Lista dei file .mat nella cartella
file_mat = dir(fullfile(cartella, '*.mat'));

% Carica il file di informazioni sulle posizioni
info_file_path = fullfile(cartella, 'FRF_analytical_co-located_1.20m_info.txt');
if ~exist(info_file_path, 'file')
    error('File di informazioni sulle posizioni non trovato.');
end

% Leggi il file delle informazioni sulle posizioni per ottenere l'indice dell'accelerometro co-locato
[acc_Index, acc_positions]  = loadPositionInfo(info_file_path);

% === Elaborazione dei file ===
n_file = length(file_mat);  % Numero di file .mat
pos_acc = zeros(1, n_file);  % Vettore delle posizioni degli accelerometri
FRF_cell = cell(1, n_file);  % Cella per memorizzare le FRF
mode_data = cell(1, n_file);  % Cella per i dati dei modi
co_located_index = [];  % Variabile per l'indice del file co-locato

% Ciclo attraverso i file .mat per estrarre i dati e calcolare le FRF
for i = 1:n_file
    nome_file = file_mat(i).name;  % Nome del file corrente
    file_path = fullfile(cartella, nome_file);  % Percorso completo del file

    % Estrazione delle posizioni dal nome del file (martellata e accelerometro)
    [xk, xj] = extractPositionsFromFilename(nome_file);
    
    % Salva la posizione dell'accelerometro
    pos_acc(i) = xj; 

    % Verifica se l'accelerometro è quello co-locato (confronta la posizione)
    if xj == acc_Index  
        co_located_index = i;  % Memorizza l'indice del file co-locato
    end

    % Caricamento dei dati dal file .mat
    load(file_path);  % Carica la variabile modeTable dal file .mat
    mode_data{i} = modeTable;  % Salva i dati della tabella dei modi
    n_modi = height(modeTable);  % Numero di modi nel file

    % Calcolo delle FRF per ogni modo
    FRF_modi = cell(1, n_modi);
    jw = 1j * omega;  % Definisce la parte immaginaria per la FRF

    for m = 1:n_modi
        wn = modeTable{m, 3};  % Pulsazione del modo [rad/s]
        xi = modeTable{m, 2};  % Smorzamento del modo
        A  = modeTable{m, 4};  % Ampiezza del modo
        R_hf = modeTable{m, 5};  % Residuo alta frequenza
        R_lf = modeTable{m, 6};  % Residuo bassa frequenza

        % Calcolo della FRF per il modo m (contributo modale + residui)
        H = A ./ (wn^2 - omega.^2 + jw * 2 * xi * wn) + R_hf + R_lf ./ jw;
        FRF_modi{m} = H;  % Memorizza la FRF del modo m
    end

    % Memorizza la cella delle FRF per il file corrente
    FRF_cell{i} = FRF_modi;
end

% === Verifica la presenza del file co-locato ===
if isempty(co_located_index)
    error('Nessuna FRF co-locata trovata.');
end

% === Calcolo delle ampiezze modali φ_ij ===
T_coll = mode_data{co_located_index};  % Tabella dei modi per il file co-locato
n_modi = height(T_coll);  % Numero di modi nel file co-locato
n_acc = n_file;  % Numero di accelerometri

% Vettori per le ampiezze modali
phi_coll = zeros(n_modi, 1);  % φ_coll per ogni modo
phi_matrix = zeros(n_acc, n_modi);  % φ(i,j) = modo j sull'accelerometro i

% Ciclo per calcolare le ampiezze modali
for m = 1:n_modi
    wn = T_coll{m, 3};  % Pulsazione del modo [rad/s]
    xi = T_coll{m, 2};  % Smorzamento del modo
    w0_idx = find(abs(omega - wn) == min(abs(omega - wn)), 1);  % Indice vicino a wn

    % FRF co-locata nel modo m
    H_coll = FRF_cell{co_located_index}{m};
    FRF_w0 = H_coll(w0_idx);  % Valore della FRF nel punto w0

    % Calcolo dell'ampiezza φ_coll (modal participation factor)
    phi_m = sqrt(abs(FRF_w0 * 2j * xi * wn^2));
    phi_coll(m) = phi_m;

    % Calcolo delle ampiezze su ogni accelerometro
    for i = 1:n_acc
        H_acc = FRF_cell{i}{m};
        FRF_acc = H_acc(w0_idx);
        phi_ij = (FRF_acc * 2j * xi * wn^2) / phi_m;  % Calcola φ_ij
        phi_matrix(i, m) = real(phi_ij);  % Memorizza la parte reale di φ_ij
    end
end

% === Plot dei Modi ===
phi_matrix(:,1) = -phi_matrix(:,1);
phi_matrix(:,3) = -phi_matrix(:,3);

% Ciclo per ogni modo modale
for m = 1:n_modi
    deformata = phi_matrix(:, m);  % Deformata del modo m
    
    % Plot della deformata del modo
    figure('Name', ['Modo ', num2str(m)], 'NumberTitle', 'off');
    plot(x, modes_shapes(m,:), 'b-', 'LineWidth', 1.5); hold on;
    plot(acc_positions, deformata', 'ro', 'MarkerFaceColor', 'r');  % Posizioni reali degli accelerometri
    xlabel('Posizione lungo la trave (m)');
    ylabel(['Ampiezza del modo ', num2str(m)]);
    title(['Deformata Modale - Modo ', num2str(m)]);
    grid on;
end

% === Funzioni di Supporto ===

function [acc_index, acc_positions] = loadPositionInfo(info_file_path)
    % Apri il file di testo per leggere le informazioni
    fileID = fopen(info_file_path, 'rt');
    if fileID == -1
        error('Impossibile aprire il file di informazioni.');
    end
    file_contents = fread(fileID, '*char')';  % Leggi tutto il file
    fclose(fileID);

    % Estrai la posizione della martellata (xk) dal file
    pattern_xk = 'Posizione forza \(martellata\): xk = (\d+\.\d+) m';
    tokens_xk = regexp(file_contents, pattern_xk, 'tokens');
    if ~isempty(tokens_xk)
        xk_file = str2double(tokens_xk{1}{1});  % Estrai la posizione xk
    else
        error('Posizione forza (martellata) non trovata nel file.');
    end

    % Estrai tutte le posizioni degli accelerometri
    pattern_acc = 'Accelerometro in xj = (\d+\.\d+) m';
    tokens_acc = regexp(file_contents, pattern_acc, 'tokens');

    % Converti le posizioni in numeri e salvale in un array
    acc_positions = zeros(1, length(tokens_acc));
    for i = 1:length(tokens_acc)
        acc_positions(i) = str2double(tokens_acc{i}{1});
    end

    % Trova l'indice dell'accelerometro co-locato (quello che ha la stessa posizione di xk)
    acc_index = find(abs(acc_positions - xk_file) < 1e-6, 1);

    % Se non si trova l'accelerometro co-locato, errore
    if isempty(acc_index)
        error('Accelerometro co-locato non trovato nel file.');
    end
end

% Funzione per estrarre le posizioni dal nome del file
function [xk, xj] = extractPositionsFromFilename(nome_file)
    tokens = regexp(nome_file, 'ModeTable_FRF_hammer_(\d+\.\d+)m_acc_(\d+)\.mat', 'tokens');
    if isempty(tokens)
        error('Nome file non conforme: %s', nome_file);
    end
    xk = str2double(tokens{1}{1});  % Posizione forza (martellata)
    xj = str2double(tokens{1}{2});  % Posizione accelerometro
end


% % === Plot interpolato per ogni modo ===
% % === (Opzionale) Polar plot simmetrico ===
% for m = 1:n_modi
% 
%     if usa_polar_plot
%         figure(n_modi + m); % numerazione figure successive
%         set(gcf, 'Name', ['Modo ', num2str(m), ' - Polar plot']);
% 
%         % Normalizza ampiezze
%         y_norm = A(m,:) ./ A_collocata(m);
%         validi = ~isnan(y_norm);
%         x_val = acc_pos(validi);
%         y_val = y_norm(validi);
% 
%         % Simmetria: specchia rispetto al centro (metà lunghezza rotaia)
%         L = max(acc_pos);         % lunghezza totale stimata
%         x_specchiati = L - x_val; % coordinate specchiate
%         y_specchiati = y_val;     % ampiezza uguale per simmetria
% 
%         % Unisci originali + specchiati
%         x_all = [x_val, x_specchiati];
%         y_all = [y_val, y_specchiati];
% 
%         % Ordina per angolo (converti posizione lineare in angolo su cerchio)
%         [x_all_sorted, idx_sort] = sort(x_all);
%         y_all_sorted = y_all(idx_sort);
%         angoli = linspace(0, pi, length(x_all_sorted)); % semi-cerchio
% 
%         % Polar plot
%         polarplot(angoli, y_all_sorted, 'b-', 'LineWidth', 1.5); hold on;
%         polarplot(angoli, y_all_sorted, 'ro', 'MarkerFaceColor', 'r'); % punti
%         title(['Modo ', num2str(m), ' (Polar plot simmetrico)']);
%     else
%         figure(n_modi); % numerazione figure successive
%         set(gcf, 'Name', ['Modo ', num2str(m), ' - Polar plot']);
% 
%             % Punti noti
%             x = acc_pos;
%             y = A(m,:) ./ A_collocata(m);
% 
%             % Rimuovi NaN
%             validi = ~isnan(y);
%             x_val = x(validi);
%             y_val = y(validi);
% 
%             % Interpolazione spline
%             x_interp = linspace(min(x_val), max(x_val), 200);
%             y_interp = spline(x_val, y_val, x_interp);
% 
%             % Plot
%             plot(x_interp, y_interp, 'b-', 'LineWidth', 1.5); hold on;
%             plot(x_val, y_val, 'ro', 'MarkerFaceColor', 'r'); hold off;
% 
%             title(['Modo ', num2str(m)]);
%             xlabel('Posizione accelerometro [m]');
%             ylabel('Ampiezza normalizzata');
%             grid on;
% 
% 
%     end
%  end