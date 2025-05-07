% === Script per caricare FRF da più file .mat ===
clear
clc
% Directory contenente i file .mat
folderPath = 'Results\Rail Optimized';  
% Trova tutti i file .mat nella cartella
matFiles = dir(fullfile(folderPath, '*.mat'));
% load("Results\Rail Optimized\freq.mat");

% Inizializza cella per contenere le FRF temporaneamente
frfList = {};
% === Leggi e ordina i file numericamente ===
matFiles = dir(fullfile(folderPath, '*.mat'));
fileNames = {matFiles.name};

% Estrai i numeri dai nomi (assumendo che siano nel formato 'FRF_numero.mat')
fileNumbers = regexp(fileNames, '\d+', 'match');
fileNumbers = cellfun(@(x) str2double(x{end}), fileNumbers);  % Prende l'ultimo numero trovato

% Ordina i file in base ai numeri
[~, sortIdx] = sort(fileNumbers);
matFiles = matFiles(sortIdx);  % Riordina la struct dei file

for k = 1:length(matFiles)
    % Costruisci il path completo del file
    filePath = fullfile(folderPath, matFiles(k).name);
    
    % Carica il file
    data = load(filePath);
    
    % Verifica che la variabile esista
    if isfield(data, 'frfData')
        frf = data.frfData;
        
        % Assicurati che sia un vettore colonna
        if isrow(frf)
            frf = frf.';
        end
        
        frfList{end+1} = frf;
    else
        warning('File %s non contiene la variabile FRF.', matFiles(k).name);
    end
    if k == length(matFiles)
        freq = data.freqData;
    end
end

% Verifica che tutte le FRF abbiano la stessa lunghezza
frfLengths = cellfun(@length, frfList);


% Costruzione della matrice finale (ogni colonna è una FRF)
FRF_matrix = cell2mat(frfList);
frf= FRF_matrix;
save(fullfile(folderPath, 'FRF_matrix.mat'), 'frf', 'freq');
