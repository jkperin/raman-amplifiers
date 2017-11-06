clear
clc

addpath fibras\ -end 
addpath Analitico\ -end
addpath Numerico_CW\ -end

dBm2Watt = inline('1e-3*10.^(P/10)');

%% Dados do bombeio
bombeio.lambda = [1421.5 1451.4];
bombeio.N_Bombeios = length(bombeio.lambda);
bombeio.P = [0.4157 0.4487];
bombeio.Bwp = 0.2000;
bombeio.FPL = 2;

%% Dados do sinal
sinal.lambda = linspace(1530.33,1560.61,20);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 0.1;
sinal.P = dBm2Watt(-20);

%% Dados da fibra
fibra = SMF(bombeio,sinal);
fibra.L = 1e3;

tic
[resultados_num, resultados_an] = DRA_Numerico(bombeio,sinal,fibra);
toc

plot(sinal.lambda,resultados_num.Ganho_on_off,'o-',sinal.lambda,resultados_an.Ganho_on_off,'o-r')
xlabel('Comprimento de Onda (nm)');
ylabel('Ganho ON/OFF (dB)');
legend('Modelo Numérico', 'Modelo Analitico');
axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(resultados_an.Ganho_on_off)])