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
fmax=200;                        %[Hz]
% resolutions
n_points = 10000; 
% create frequency vectors
f=linspace(0,fmax,n_points);    % [Hz]
omega=2*pi*f;                   %[rad/s]
% Setting the space domain along the beam
x=linspace(0,L,n_points);

%% Calculation
% Building the matrix of the coefficients from the BCs for a cantilever
% beam
H=@(omega) [            1                                   0                               1                               0    ;
                        0                                   1                               0                               1    ;
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
C_hat = zeros(4, length(i_nat)); 

% Salvo i valori dei coefficienti quando omega è pari alla frequenza
% propria (in pratica salvo i modi) nella matrice C_hat
for i_mode = 1:length(i_nat)
    % trovo omega0 (frequenza propria)
    omega_i = omega(i_nat(i_mode));
    
    % calcolo la matrice H quando omega = omega_i (frequenza propria)
    Hi = H(omega_i);
    
    % estraggo la parte ridotta della matrice Hi (2:4, 2:4)
    Hi_hat = Hi(2:4, 2:4);  % 3x3 matrix
    Ei_hat = Hi(2:4, 1);    % 3x1 vector
    
    % Risolvo il sistema Hi_hat * Ci_hat = -Ei_hat
    % Trovo i modi di vibrazione per la frequenza propria
    Ci_hat = [1; -Hi_hat\Ei_hat];  % Il primo componente è 1 e gli altri sono risolti
    
    % Salvo il vettore dei modi nella matrice C_hat
    C_hat(:, i_mode) = Ci_hat;
end

%% Modal shapes computation
% Distretizzo l'asse della lunghezza della trave
x = linspace(0, L, 1000);  % vettore di lunghezza della trave
dx = x(2);  % distanza fra i punti

% Calcolo dei modi per ogni frequenza propria
phi = zeros(length(i_nat), length(x));  % inizializzo la matrice dei modi
for i_mode = 1:length(i_nat)
    omega_i = omega(i_nat(i_mode));  % Frequenza propria
    gamma_i = (L * m * omega_i^2 / (E * J))^(1/4);  % parametro di vibrazione gamma_i
    
    % Calcolo del modo di vibrazione per ogni punto dell'asse x
    phi(i_mode, :) = C_hat(1, i_mode) * cos(gamma_i * x) ...
                    + C_hat(2, i_mode) * sin(gamma_i * x) ...
                    + C_hat(3, i_mode) * cosh(gamma_i * x) ...
                    + C_hat(4, i_mode) * sinh(gamma_i * x);
end

% Normalizzo i modi per far sì che il massimo valore sia 1
normaliz = max(abs(phi), [], 2);  % Normalizzo per ciascun modo (per ogni riga)
phi = phi ./ normaliz;  % Normalizzazione dei modi

% Se vuoi normalizzare ogni modo separatamente per avere il primo valore pari a 1:
% phi = phi ./ phi(:,1);  % Normalizza ogni modo rispetto al primo valore



%% Animation
% Seleziona il modo di cui vuoi fare l'animazione
mode = 1;

% Se i colori non sono definiti, possiamo crearne uno
colors_p = lines(length(i_nat));  % Usa una palette di colori predefiniti

figure, hold on, grid on
title(['Mode ', num2str(mode)])
plot(x, phi(mode,:), ':k', 'LineWidth', 2)
h1 = plot(x, zeros(size(x)), 'LineWidth', 2, 'color', colors_p(mode, :));

% Crea l'animazione
for t = linspace(0, 5 / f(i_nat(mode)), 400)
    if ishandle(h1)
        % Calcola la forma modale oscillante nel tempo (sinusoidale)
        w1 = phi(mode,:) * cos(omega(i_nat(mode)) * t);  % oscillazione coseno
        h1.YData = w1;  % Aggiorna i dati Y della curva
        pause(0.005)  % Pausa per il frame dell'animazione
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
    modal_mass(i,1) = m * trapz(phi(i,:).^2,x);
end

% omega vector
w = f * 2 * pi;

% transfer function initialization
FRF = zeros(1,length(w));

% transfer function calculation
for k = 1:length(w)
    FRF(k) = 0; % inizializza a zero
    for i = 1:length(i_nat)
        FRF(k) = FRF(k) + ( ( phi(i,pos_xj) * phi(i,pos_xk) ) / modal_mass(i) ) / ...
                     ( -w(k)^2 + 2i * xsi * w(k) * omega(i_nat(i)) + omega(i_nat(i))^2 );
    end
end

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

% Posizione accelerometro e martellata
xj = 0.2;  % [m]
xk = 1.2;  % [m]

% Crea la cartella "Results" se non esiste
if ~exist('Results', 'dir')
    mkdir('Results');
end

% Determina nome file in base alle posizioni
if abs(xj - xk) < 1e-6  % tolleranza per confronti floating point
    suffix = sprintf('collocata_in_%.2fm', xj);
else
    suffix = sprintf('xj_%.2fm_xk_%.2fm', xj, xk);
end

% Sostituisci il punto con la virgola nel nome file (opzionale, stile EU)
suffix = strrep(suffix, '.', ',');

% Crea una matrice con frequenza [Hz], parte reale e parte immaginaria
FRF_data = [f.' real(FRF).' imag(FRF).'];

% Percorsi completi dei file
csv_path = fullfile('Results', ['FRF_export_' suffix '.csv']);
mat_path = fullfile('Results', ['FRF_export_' suffix '.mat']);

% Salva su file CSV
writematrix(FRF_data, csv_path);

% Salva anche su file .mat (mantiene il numero complesso)
save(mat_path, 'f', 'FRF');

% Messaggio di conferma
disp(['FRF esportata in: ', csv_path, ' e ', mat_path]);
