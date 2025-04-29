clear all
close all
clc

%% BEAM 
%% Definition of the mechanical properties of the system

h = 0.008;               % thickness [m]
b = 0.04;                % width [m]
rho = 2700;             % density [kg/^3]
L = 1.2;                % beam length [m]  
E = 68e9;               % Young Modulus [Pa]
J = b*h^3/12;           % Inertia moment [m^4]
V = L*b*h;              % volume [m^3]
M = rho*V;
m = M/L; % mass [Kg]
xsi = 0.01;             % smorzamento adimensionale

% Setting the frequency range
fmax = 200;                        %[Hz]
% resolutions
n_points = 10000; 
% create frequency vectors
f=linspace(0,fmax,n_points);    % [Hz]
omega=2*pi*f;                   %[rad/s]


%% Calculation
% Building the matrix of the coefficients from the BCs for a cantilever
% beam
H=@(omega) [            1                                       0                                           1                                           0    ;
                        0                                       1                                           0                                           1    ;
              -cos(L*(m*omega^2/(E*J))^(1/4))        -sin(L*(m*omega^2/(E*J))^(1/4))            cosh(L*(m*omega^2/(E*J))^(1/4))             sinh(L*(m*omega^2/(E*J))^(1/4));
              sin(L*(m*omega^2/(E*J))^(1/4))         -cos(L*(m*omega^2/(E*J))^(1/4))            sinh(L*(m*omega^2/(E*J))^(1/4))             cosh(L*(m*omega^2/(E*J))^(1/4));];

% inizializzo il vettore dove inserirò il valore del determinante
dets = zeros(length(omega),1);
% calculate H matrix determinant in the frequency range
for i=1:length(omega)
    dets(i)=det(H(omega(i)));
end

% plot determinant
figure, box on
semilogy(f,abs(dets),'-b')
hold on;
grid on;
grid minor;
ylabel('det(H)');
xlabel('f [Hz]');
title("H Matrix determinant")


%% Imposing that the determinant is null
% inizializzo il vettore dove metto gli indici per cui si annulla il
% determinate (frequenze proprie)
i_nat=[];
% trovo i minimi locali del determinante e salvo gli indici in i_nat
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end
% stampo a schermo le frequenze proprie
fprintf('Natural frequencies [Hz]:\n ');
disp(f(i_nat));
% aggiungo nel plot precedente dei pallini dove ho le frequenze proprie
plot(f(i_nat),abs(dets(i_nat)),'or')

%% Solving the reduced system
% Ora sappiamo i valori di omega (frequenze proprie) per cui il sistema è singolare, 
% quindi possiamo risolvere il sistema ridotto per trovare i modi.

% inizializzo la matrice dove in ogni colonna metterò i modi per diverse
% frequenze proprie
X_hat = zeros(4, length(i_nat)); 

% Salvo i valori dei coefficienti quando omega è pari alla frequenza
% propria (in pratica salvo i modi) nella matrice C_hat
for i_mode = 1:length(i_nat)
    % trovo omega0 (frequenza propria)
    omega_i = omega(i_nat(i_mode));
    
    % calcolo la matrice H quando omega = omega_i (frequenza propria)
    Hi = H(omega_i);
    
    % estraggo la parte ridotta della matrice Hi (2:4, 2:4)
    Hi_hat = Hi(2:4, 2:4);  % 3x3 matrix
    Ni_hat = Hi(2:4, 1);    % 3x1 vector
    
    % Risolvo il sistema Hi_hat * Xi_hat = -Ni_hat
    % Trovo il vettore del modo di vibrare e lo aggiungo nella matrce dei modi
    X_hat(:, i_mode) = [1; -Hi_hat\Ni_hat];  % Il primo componente è 1 e gli altri sono risolti
    
end

%% Modal shapes computation
% Distretizzo l'asse della lunghezza della trave
x = linspace(0, L, n_points);  % vettore di lunghezza della trave
dx = x(2);  % distanza fra i punti

% Calcolo dei modi per ogni frequenza propria
modes_shapes = zeros(length(i_nat), length(x));  % inizializzo la matrice dei modi

for i_mode = 1:length(i_nat)
    omega_i = omega(i_nat(i_mode));  % Frequenza propria
    gamma_i = (L * m * omega_i^2 / (E * J))^(1/4);  % parametro di vibrazione gamma_i
    
    % Calcolo del modo di vibrazione per ogni punto dell'asse x
    modes_shapes(i_mode, :) = X_hat(1, i_mode) * cos(gamma_i * x) ...
                    + X_hat(2, i_mode) * sin(gamma_i * x) ...
                    + X_hat(3, i_mode) * cosh(gamma_i * x) ...
                    + X_hat(4, i_mode) * sinh(gamma_i * x);
end

% Normalizzo i modi per far sì che il massimo valore sia 1
normaliz = max(abs(modes_shapes), [], 2);  % Normalizzo per ciascun modo (per ogni riga)
modes_shapes = modes_shapes ./ normaliz;  % Normalizzazione dei modi

% Se vuoi normalizzare ogni modo separatamente per avere il primo valore pari a 1:
% phi = phi ./ phi(:,1);  % Normalizza ogni modo rispetto al primo valore



% Animation
mode = 1;

colors_p = lines(length(i_nat));  % palette colori

figure('Position', [100, 100, 1200, 500]); hold on; grid on;
title(['Mode ', num2str(mode)])
plot(x, modes_shapes(mode,:), ':k', 'LineWidth', 2)
h1 = plot(x, zeros(size(x)), 'LineWidth', 2, 'color', colors_p(mode, :));
xlabel('Posizione [m]')
ylabel('Deformazione (modale)')
ylim([-1, 1] * max(abs(modes_shapes(mode,:))) * 1.1)

% Parametri dell’animazione
T = 1 / f(i_nat(mode));     % Periodo reale della frequenza naturale
n_cycles = 4;               % Numero di cicli da animare
n_frames = 1000;             % Numero di frame totali

t = linspace(0, n_cycles*T, n_frames);  % tempo continuo

for k = 1:length(t)
    if ishandle(h1)
        w1 = modes_shapes(mode,:) * cos(2*pi*f(i_nat(mode)) * t(k));  % uso 2πf per l'argomento del coseno
        h1.YData = w1;
        pause(T/n_frames);
    else
        return
    end
end


%% Transfer Function

% posizione acceleroemtro
xj = 0.2; %[m]
% posizione martellata
xk = 1.2; %[m]

% Trova l'indice di xj nel vettore x
[~, pos_xj] = min(abs(x - xj));

% Trova l'indice di xk nel vettore x
[~, pos_xk] = min(abs(x - xk));

% modal mass vector initialization
modal_mass = zeros(1,length(i_nat));

% modal mass calculation for every mode
for i = 1:length(i_nat)
    modal_mass(i,1) = m * trapz(x,modes_shapes(i,:).^2);
end


% transfer function initialization
FRF = zeros(1,length(omega));

% transfer function calculation
for k = 1:length(omega)
    FRF(k) = 0; % inizializza a zero
    for i = 1:length(i_nat)
        FRF(k) = FRF(k) + ( ( modes_shapes(i,pos_xj) * modes_shapes(i,pos_xk) ) / modal_mass(i) ) / ...
                     ( -omega(k)^2 + 2i * xsi * omega(k) * omega(i_nat(i)) + omega(i_nat(i))^2 );
    end
end

% % transfer function calculation for single unit impulse force
% for i = 1:length(i_nat)
%     num = modes_shapes(i, pos_xj) * modes_shapes(i, pos_xk);
%     den = -omega.^2 + 2i * xsi * omega * omega(i_nat(i)) + omega(i_nat(i))^2;
%     FRF = FRF + (num / modal_mass(i)) ./ den;
% end

%% FRF plot 
% Calcolo ampiezza e fase della FRF
FRF_amp = abs(FRF);             % modulo
FRF_phase = angle(FRF);         % fase in radianti

% Crea la figura con due subplot
figure;

% ---- Ampiezza (semilogaritmica) ----
subplot(2,1,1);
semilogy(f, FRF_amp, 'b', 'LineWidth', 2);
grid on;
xlabel('Frequenza [Hz]');
ylabel('|G(j\omega)|');
title('Funzione di Trasferimento - Ampiezza');
xlim([0 max(f)]);

% ---- Fase ----
subplot(2,1,2);
plot(f, FRF_phase, 'r', 'LineWidth', 2);
grid on;
xlabel('Frequenza [Hz]');
ylabel('Fase [rad]');
title('Funzione di Trasferimento - Fase');
xlim([0 max(f)]);

%% Export FRF (numero complesso) con nome file in base alla posizione

% Crea la cartella "Results" se non esiste
if ~exist('Results', 'dir')
    mkdir('Results');
end

% Determina nome file in base alle posizioni
if abs(xj - xk) < 1e-6
    suffix = sprintf('collocata_in_%.2fm', xj);
else
    suffix = sprintf('xj_%.2fm_xk_%.2fm', xj, xk);
end

% Rimuovi caratteri non validi per nome file
suffix = strrep(suffix, '.', '_');  % cambia punto in underscore

% Crea una matrice con frequenza [Hz], parte reale e parte immaginaria
FRF_data = [f.' real(FRF).' imag(FRF).'];

% Percorsi completi dei file
csv_path = fullfile('Results', ['FRF_export_' suffix '.csv']);
mat_path = fullfile('Results', ['FRF_export_' suffix '.mat']);

% Salva su file CSV
writematrix(FRF_data, csv_path);

% Salva anche su file .mat
save(mat_path, 'f', 'FRF');

% Messaggio di conferma
disp(['FRF esportata in: ', csv_path, ' e ', mat_path]);
