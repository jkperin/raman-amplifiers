%% Numerical comparison between conventional dispersion compensating fibers 
%% and photonic crystal fibers as lumped Raman amplifiers
clear, close all

dBm2Watt = @(x) 1e-3*10^(x/10);

addpath fibras\ -end
addpath Analitico\ -end
addpath Numerico_CW\ -end

%% Dados do sinal - 20 canais WDM na banda C
sinal.lambda = linspace(1530.33, 1560.61, 20);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 0.1;
sinal.P = dBm2Watt(-10);

%% Caso A

% Dados do bombeio
% DCF
bombeioDCF.lambda = [1413.2 1446.7];
bombeioDCF.N_Bombeios = length(bombeioDCF.lambda);
bombeioDCF.P = [0.3193 0.3827];
bombeioDCF.Bwp = 0.1;
bombeioDCF.FPL = 2;
% PCF, alfa = -5dB/km
bombeioPCF5.lambda = [1421.5 1451.6];
bombeioPCF5.N_Bombeios = length(bombeioPCF5.lambda);
bombeioPCF5.P = [0.45 0.45];
bombeioPCF5.Bwp = 0.1;
bombeioPCF5.FPL = 2;
% PCF, alfa = -4 dB/km
bombeioPCF4.lambda = [1421.3 1451.6];
bombeioPCF4.N_Bombeios = length(bombeioPCF4.lambda);
bombeioPCF4.P = [0.4497 0.45];
bombeioPCF4.Bwp = 0.1;
bombeioPCF4.FPL = 2;
% PCF, alfa = -3 dB/km
bombeioPCF3.lambda = [1421.5 1451.4];
bombeioPCF3.N_Bombeios = length(bombeioPCF3.lambda);
bombeioPCF3.P = [0.4157 0.4487];
bombeioPCF3.Bwp = 0.1;
bombeioPCF3.FPL = 2;

% fibra
% DCF
fibraDCF = DCF(bombeioDCF, sinal);
fibraDCF.L = 12.15e3;

% PCF alfa = -5dB/km
fibraPCF5 = PCF(bombeioPCF5, sinal);
fibraPCF5.L = 1e3;

% PCF alfa = -4dB/km
fibraPCF4 = PCF(bombeioPCF4, sinal, 4);
fibraPCF4.L = 1e3;

% PCF alfa = -3dB/km
fibraPCF3 = PCF(bombeioPCF3, sinal, 3);
fibraPCF3.L = 1e3;

% Roda modelos
disp('DCF');
[DCF.num, DCF.ana] = DRA_Numerico(bombeioDCF, sinal, fibraDCF);
fprintf('Ganho Médio = %.2f dB [NUM], %.2f dB [ANA]\n', DCF.num.Ganho_Medio, DCF.ana.Ganho_Medio);
fprintf('Ripple = %.2f dB [NUM], %.2f dB [ANA]\n', DCF.num.ripple, DCF.ana.ripple);

disp('PCF5');
[PCF5.num, PCF5.ana] = DRA_Numerico(bombeioPCF5, sinal, fibraPCF5);
fprintf('Ganho Médio = %.2f dB [NUM], %.2f dB [ANA]\n', PCF5.num.Ganho_Medio, PCF5.ana.Ganho_Medio);
fprintf('Ripple = %.2f dB [NUM], %.2f dB [ANA]\n', PCF5.num.ripple, PCF5.ana.ripple);

disp('PCF4');
[PCF4.num, PCF4.ana] = DRA_Numerico(bombeioPCF4, sinal, fibraPCF4);
fprintf('Ganho Médio = %.2f dB [NUM], %.2f dB [ANA]\n', PCF4.num.Ganho_Medio, PCF4.ana.Ganho_Medio);
fprintf('Ripple = %.2f dB [NUM], %.2f dB [ANA]\n', PCF4.num.ripple, PCF4.ana.ripple);

disp('PCF3');
[PCF3.num, PCF3.ana] = DRA_Numerico(bombeioPCF3, sinal, fibraPCF3);
fprintf('Ganho Médio = %.2f dB [NUM], %.2f dB [ANA]\n', PCF3.num.Ganho_Medio, PCF3.ana.Ganho_Medio);
fprintf('Ripple = %.2f dB [NUM], %.2f dB [ANA]\n', PCF3.num.ripple, PCF3.ana.ripple);

%% Gráficos
plot(sinal.lambda, DCF.num.Ganho_on_off, '-ok', sinal.lambda, DCF.ana.Ganho_on_off, '-or')
hold on
plot(sinal.lambda, PCF5.num.Ganho_on_off, '-sk', sinal.lambda, PCF5.ana.Ganho_on_off, '-sr')
plot(sinal.lambda, PCF4.num.Ganho_on_off, '-*k', sinal.lambda, PCF4.ana.Ganho_on_off, '-*r')
plot(sinal.lambda, PCF3.num.Ganho_on_off, '-dk', sinal.lambda, PCF3.ana.Ganho_on_off, '-dr')
legend('DCF num', 'DCF ana', 'PCF5 num', 'PCF5 ana', 'PCF4 num', 'PCF4 ana', 'PCF3 num', 'PCF3 ana');



