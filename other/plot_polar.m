clc
clear



phi = [
-1288.35887015845
-915.743086804604
95.8794561147716
0
1357.40359540156
940.479927944467
0
-1010.95550148070
-1350.74123045745
0
0
983.314718735323
];



% Raggio di riferimento
R = 10000;  % puoi modificarlo per scalare la ciambella

% Calcola raggi deformati
R_deformato = R + phi;

% Angoli in radianti (12 accelerometri ogni 15°)
theta = deg2rad(0:15:165);
phi_symmetry = R_deformato;

% % Chiudi il cerchio per continuità
% R_deformato = [R_deformato; R_deformato(1)];
% theta = [theta, theta(1)];
theta_circle = deg2rad(0:15:360);
circle = ones(size(theta_circle')) * R;

% === Plot ===
figure;
polarplot(theta_circle, circle, 'k--', 'LineWidth', 1.0); 
hold on;
theta_fine = linspace(min(theta), max(theta), 360);  % Più punti angolari
R_interp = interp1(theta, R_deformato, theta_fine, 'spline');  % Interpolazione spline

% Traccia curva continua
polarplot(theta_fine, R_interp, 'r-', 'LineWidth', 2);
hold on;
polarplot(-theta_fine, R_interp, 'b-o', 'LineWidth', 2);
rlim([0, R + max(abs(phi))*1.1]);
title('Modo deformato - rappresentazione polare');