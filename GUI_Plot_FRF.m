clc
clear
close all

%% === IMPORT DATI ===
[filename, pathname] = uigetfile('*.mat', 'Seleziona un file MAT');
load(fullfile(pathname, filename)); % caricato dentro 'data'

%% === SETUP DATI ===
num_accelerometri = size(frf, 2);
frequency = freq;
H1 = arrayfun(@(i) frf(:, i), 1:num_accelerometri, 'UniformOutput', false);
hasCoherence = exist('cohe', 'var') == 1;

ax_coh = []; 
if hasCoherence
    coherence = arrayfun(@(i) cohe(:,i), 1:num_accelerometri, 'UniformOutput', false);
else
    coherence = [];
end
colors = lines(num_accelerometri);

%% === CREAZIONE GUI ===
f_gui = uifigure('Name', 'FRF Media', 'Position', [100 50 1500 950]);
tg = uitabgroup(f_gui, 'Position', [10 10 1480 930]);

%% === TABS SINGOLI PER ACCELEROMETRO ===
for i = 1:num_accelerometri
    t = uitab(tg, 'Title', ['Acc. ' num2str(i)]);

    ax1 = uiaxes(t, 'Position', [60 650 1350 250]);
    semilogy(ax1, frequency, abs(H1{i}), 'Color', colors(i,:), 'LineWidth', 1.5);
    grid(ax1, 'on')
    xlabel(ax1, 'Frequenza [Hz]')
    ylabel(ax1, 'Ampiezza')
    title(ax1, ['Modulo FRF - Accelerometro ' num2str(i)])

    ax2 = uiaxes(t, 'Position', [60 360 1350 250]);
    plot(ax2, frequency, angle(H1{i}), 'Color', colors(i,:), 'LineWidth', 1.5);
    grid(ax2, 'on')
    xlabel(ax2, 'Frequenza [Hz]')
    ylabel(ax2, 'Fase [rad]')
    title(ax2, ['Fase FRF - Accelerometro ' num2str(i)])

    if hasCoherence
        ax3 = uiaxes(t, 'Position', [60 70 1350 250]);
        plot(ax3, frequency, coherence{i}, 'Color', colors(i,:), 'LineWidth', 1.5);
        grid(ax3, 'on')
        xlabel(ax3, 'Frequenza [Hz]')
        ylabel(ax3, 'Coerenza')
        title(ax3, ['Coerenza - Accelerometro ' num2str(i)])
    end
end

%% === TAB MULTI-ACCELEROMETRO ===
t_multi = uitab(tg, 'Title', 'Multi-FRF');

% Pannello checkbox a sinistra
panel_cb = uipanel(t_multi, 'Title', 'Seleziona Accelerometri', ...
    'Position', [10 10 160 880]);

cb = gobjects(num_accelerometri,1);
color_patch = gobjects(num_accelerometri,1);

for i = 1:num_accelerometri
    cb(i) = uicheckbox(panel_cb, 'Text', ['Acc. ' num2str(i)], ...
        'Position', [30 830 - (i-1)*30 120 22], 'Value', true);

    color_patch(i) = uilabel(panel_cb, ...
        'Position', [10 832 - (i-1)*30 15 15], ...
        'BackgroundColor', colors(i,:), ...
        'Text', '');
end

% Assi grafici a destra
ax_amp = uiaxes(t_multi, 'Position', [200 640 1260 240]); ax_amp.YScale = 'log';
title(ax_amp, 'Modulo FRF')
xlabel(ax_amp, 'Frequenza [Hz]')
ylabel(ax_amp, 'Ampiezza')
grid(ax_amp, 'on'); hold(ax_amp, 'on')

ax_phase = uiaxes(t_multi, 'Position', [200 360 1260 240]);
title(ax_phase, 'Fase FRF')
xlabel(ax_phase, 'Frequenza [Hz]')
ylabel(ax_phase, 'Fase [rad]')
grid(ax_phase, 'on'); hold(ax_phase, 'on')

if hasCoherence
    ax_coh = uiaxes(t_multi, 'Position', [200 80 1260 240]);
    title(ax_coh, 'Coerenza')
    xlabel(ax_coh, 'Frequenza [Hz]')
    ylabel(ax_coh, 'Coerenza')
    grid(ax_coh, 'on'); hold(ax_coh, 'on')
end

%% === FUNZIONE DI AGGIORNAMENTO ===
function update_multi(ax_amp, ax_phase, ax_coh, H1, coherence, frequency, cb, colors)
    cla(ax_amp); cla(ax_phase);
    hold(ax_amp, 'on'); hold(ax_phase, 'on');
    
    if ~isempty(ax_coh)
        cla(ax_coh); hold(ax_coh, 'on');
    end


    for j = 1:numel(cb)
        if cb(j).Value
            semilogy(ax_amp, frequency, abs(H1{j}), 'Color', colors(j,:), 'LineWidth', 1.2);
            plot(ax_phase, frequency, angle(H1{j}), 'Color', colors(j,:), 'LineWidth', 1.2);
            if ~isempty(ax_coh)
                plot(ax_coh, frequency, coherence{j}, 'Color', colors(j,:), 'LineWidth', 1.2);
            end

        end
    end
end

%% === CALLBACKS ===
for i = 1:num_accelerometri
   cb(i).ValueChangedFcn = @(src, event) ...
    update_multi(ax_amp, ax_phase, ax_coh, H1, coherence, frequency, cb, colors);
end

% Primo disegno
update_multi(ax_amp, ax_phase, ax_coh, H1, coherence, frequency, cb, colors);
