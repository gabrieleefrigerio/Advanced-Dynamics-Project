clear 
close all
clc

%% IMPORT 

% carico la FRF analitica
load("Results\Analytic\FRF_analytical_co-located_1.20m.mat");

n_acc = 8;  % Scegli l'accelerometro (1–8)

% Directory base
dataDir = "Results/Beam Optimized/";

% Costruzione del nome del file in base a n_acc con switch
switch n_acc
    case 1
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_1.mat";
        xj = 0.1;
    case 2
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_2.mat";
        xj = 0.2;
    case 3
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_3.mat";
        xj = 0.3;
    case 4
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_4.mat";
        xj = 0.5;
    case 5
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_5.mat";
        xj = 0.7;
    case 6
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_6.mat";
        xj = 0.8;
    case 7
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_7.mat";
        xj = 1;
    case 8
        fileName = "FRF_SDOF_Optimize_hammer_1.20m_acc_8.mat"; % collocata
        xj = 1.2;
    otherwise
        error("Valore di n_acc non valido. Deve essere tra 1 e 8.");
end

% Caricamento del file corrispondente
filePath = fullfile(dataDir, fileName);
load(filePath);

xk = extractPositionsFromFilename(fileName);


%%


% Plot di confronto con ampiezza sopra e fase sotto
figure('Color', 'w', 'Name', 'Comparison Plot', 'Position', [100, 100, 1700, 900]);

% Titolo dinamico con posizioni
titolo = sprintf("Input at x_k = %.2f m, Output at x_j = %.1f m", xk, xj);


% --- Subplot 1: Ampiezza (modulo)
subplot(2,1,1);
semilogy(freq, abs(frf(:,n_acc)), 'b', 'LineWidth', 1.5);
hold on;
semilogy(freqData, abs(frfData), '--r', 'LineWidth', 1.5);
legend("FRF Analytic", "FRF Optimized");
xlabel("Frequenza [Hz]");
ylabel("|FRF|");
title(["Amplitude FRF", titolo], 'FontWeight', 'bold');
grid on;

% --- Subplot 2: Fase
subplot(2,1,2);
plot(freq, unwrap(angle(frf(:,n_acc))), 'b', 'LineWidth', 1.5);
hold on;
plot(freqData, unwrap(angle(frfData)), '--r', 'LineWidth', 1.5);
legend("FRF Analytic", "FRF Optimized");
xlabel("Frequenza [Hz]");
ylabel("Fase [rad]");
title(["Phase FRF", titolo], 'FontWeight', 'bold');
grid on;


function [xk, xj] = extractPositionsFromFilename(nome_file)
    tokens = regexp(nome_file, 'FRF_SDOF_Optimize_hammer_(\d+\.\d+)m_acc_(\d+)\.mat', 'tokens');
    if isempty(tokens)
        error('Nome file non conforme: %s', nome_file);
    end
    xk = str2double(tokens{1}{1});  % Posizione forza (martellata)
    xj = str2double(tokens{1}{2});  % Posizione accelerometro
end