clc
clear
close all

%% IMPORT

% load freq vector
load("Results\Rail Optimized\FRF_matrix.mat");

% Seleziona la cartella contenente i file .mat e .txt
% cartella = uigetdir('Seleziona la cartella contenente i file .mat e .txt');
cartella = [ pwd '\Results\Rail Optimized\Table Modi'];
% Leggi e ordina i file .mat numericamente
file_struct = dir(fullfile(cartella, '*.mat'));
file_names = {file_struct.name};

% Estrai i numeri dai nomi dei file (assumendo che contengano numeri)
file_numbers = regexp(file_names, '\d+', 'match');
file_numbers = cellfun(@(x) str2double(x{end}), file_numbers);  % Usa l'ultimo numero nel nome

% Ordina i file in base ai numeri
[~, sort_idx] = sort(file_numbers);
file_mat = file_struct(sort_idx);  % Ora file_mat è ordinato correttamente


% === Elaborazione dei file ===
deltadeg = 15;

omega = freq * 2 * pi;

% === Elaborazione dei file ===
n_file = length(file_mat);  % Numero di file .mat

% Nuova definizione: posizione degli accelerometri lungo la circonferenza
% Accelerometri distanziati di 15° (convertiti in radianti), acc 1 = 0°
pos_acc = zeros(1, n_file);  % Angoli in radianti
FRF_cell = cell(1, n_file);  % Cella per memorizzare le FRF
mode_data = cell(1, n_file);  % Cella per i dati dei modi
co_located_index = 12;  % Variabile per l'indice del file co-locato

% Ciclo attraverso i file .mat per estrarre i dati e calcolare le FRF
for i = 1:n_file
    nome_file = file_mat(i).name;  % Nome del file corrente
    file_path = fullfile(cartella, nome_file);  % Percorso completo del file

    % Assegna la posizione angolare dell'accelerometro in radianti
    pos_acc(i) = deg2rad((i - 1) * 15);  % 0°, 15°, 30°, ... → in radianti


    % Caricamento dei dati dal file .mat
    load(file_path);  % Carica la variabile modeTable dal file .mat
    mode_data{i} = modeTable;  % Salva i dati della tabella dei modi
    n_modi = height(modeTable);  % Numero di modi nel file

    % Calcolo delle FRF per ogni modo
    FRF_modi = cell(1, n_modi);
    jw = 1j * omega;  % Parte immaginaria per la FRF

    for m = 1:n_modi
        wn = modeTable{m, 3};  % Pulsazione [rad/s]
        xi = modeTable{m, 2};  % Smorzamento
        A  = modeTable{m, 4};  % Ampiezza
        R_hf = modeTable{m, 5};  % Residuo alta frequenza
        R_lf = modeTable{m, 6};  % Residuo bassa frequenza

        % Calcolo della FRF (componente modale + residui)
        H = A ./ (wn^2 - omega.^2 + jw * 2 * xi * wn) + R_hf + R_lf ./ jw;
        FRF_modi{m} = H;
    end

    % Salva la FRF per il file corrente
    FRF_cell{i} = FRF_modi;
end


% === Calcolo delle ampiezze modali φ_ij ===
T_coll = mode_data{co_located_index};
n_modi = height(T_coll);
n_acc = n_file;

% Vettori per le ampiezze modali
phi_coll = zeros(n_modi, 1);  % φ_coll per ogni modo
phi_matrix = zeros(n_acc, n_modi);  % φ(i,j) = modo j sull'accelerometro i
fittato = zeros(n_acc,1);
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
    if length(FRF_cell{i}) < n_modi && fittato(i, 1) == 0
        warning("⚠️ Accelerometro %d ha meno modi di %d. Salto.", i, 2);
        dati_mode = mode_data{i};
        if dati_mode{1, 3} <= 1000*2*pi
           phi_matrix(i, m+1) = 0;
           fittato(i, 1) = 1;
        %continue;
        elseif dati_mode{1, 3} > 1000*2*pi
            if m == 1
               phi_matrix(i, m) = 0;
               continue;
            else
                H_acc = FRF_cell{i}{1};  % Vettore FRF del modo m per l'accelerometro i
                FRF_acc = H_acc(w0_idx);  % Estrai valore della FRF a ω₀
                phi_ij = (FRF_acc * 2j * xi * wn^2) / phi_m;  % Calcola φ_ij
                phi_matrix(i, m) = real(phi_ij);
                continue;
            end
        end
    elseif fittato(i, 1) == 1
        continue;
    end

    H_acc = FRF_cell{i}{m};  % Vettore FRF del modo m per l'accelerometro i
    FRF_acc = H_acc(w0_idx);  % Estrai valore della FRF a ω₀
    phi_ij = (FRF_acc * 2j * xi * wn^2) / phi_m;  % Calcola φ_ij
    phi_matrix(i, m) = real(phi_ij);
end

end

% === Polar plot dettagliato del modo finale ===
modo = 1;
phi = phi_matrix(:, modo);                   % Dati modali identificati
theta = deg2rad(0:15:165) - pi/2;            % Angoli in radianti + rotazione 90° oraria

% Raggio di riferimento
R = 10000;                                   % Raggio base per la ciambella

% Calcola raggi deformati
R_deformato = R + phi;

% Cerchio non deformato (per riferimento)
theta_circle = deg2rad(0:1:360);     % Anche il cerchio ruotato
circle = ones(size(theta_circle')) * R;

% === Plot ===
figure;
polarplot(theta_circle, circle, 'k--', 'LineWidth', 1.0);  % Cerchio di riferimento
hold on;
theta_fine = linspace(min(theta), max(theta), 360);  % Più punti angolari
R_interp = interp1(theta, R_deformato, theta_fine, 'spline');  % Interpolazione spline
polarplot(theta_fine, R_interp, 'r-', 'LineWidth', 1);      % Forma modale
hold on;
polarplot(-theta_fine, -R_interp, 'b-', 'LineWidth', 1);     % Forma riflessa
polarplot(theta, R_deformato, 'ro', 'LineWidth', 2);      % Forma modale
polarplot(-theta, -R_deformato, 'bo', 'LineWidth', 2);     % Forma riflessa
rlim([0, R + max(abs(phi))*1.1]); 
title('Modo deformato - rappresentazione polare');

% Etichette angolari (accelerometri)
for i = 1:length(theta)
    r_text = R + max(abs(phi)) * 1.05;  % Posizione leggermente esterna
    text(theta(i), r_text, sprintf('%d', i), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 10);
end

% Titolo e legenda
title(sprintf('%d^{\\rm th} Axial Mode Shape', modo), 'FontSize', 14);
legend({'Undeformed', 'Axial mode shape (identified)', 'Axial mode shape (symmetry)'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal');
hold off;

