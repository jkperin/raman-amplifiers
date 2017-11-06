% % 1522:1:1562 40 canais
% % 1520-1595 75 canais
% % 
% % -10 dBm/Channel
% % 
% % 25 Km SMF
% % 
% % A) 1420, 1450, 1480 
% % 174, 300, 26
% % 40 canais
% % 
% % B) 1420, 1450, 1480 
% % 391,600,9
% % 40 canais
% % 
% % 
% % C) 1420, 1440.25, 1480
% % 166.6, 166.6, 166.6
% % 75canais

clear
clc

%% Dados Experimentais
a = [1522, 3.75; 1524, 4.1; 1526, 4.25; 1528, 4.2; 1530, 4.1; 1532, 4; 1534, 3.9;...
    1536, 3.6; 1538, 3.6; 1540, 3.9; 1542, 4; 1544, 4.1; 1546, 4.2; 1548, 4.2;...
    1550, 4.25; 1552, 4.25; 1554, 4.1; 1556, 4.1; 1558, 4.3; 1560, 4.4; 1562, 3.6];


%% Simulação

addpath ..\fibras\ -end 
addpath ..\Analitico\ -end
addpath ..\Numerico_CW\ -end

dBm2Watt = inline('1e-3*10.^(P/10)');

%% Dados do sinal - 40 canais WDM na banda C
sinal.lambda = 1522:1:1561;
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 1;
sinal.P = dBm2Watt(-10);

%% Caso A

% Dados do bombeio
bombeioA.lambda = [1420 1450 1480];
bombeioA.N_Bombeios = length(bombeioA.lambda);
bombeioA.P = [0.174, 0.300, 0.026];
bombeioA.Bwp = 1;
bombeioA.FPL = 2;

% Dados da fibra
fibraA = SMF(bombeioA,sinal);
fibraA.L = 20e3;


%% Resultados númericos e analíticos
disp('Caso A')
tic
[A.resultados_num,A.resultados_an] = DRA_Numerico(bombeioA, sinal, fibraA);
toc
% disp('Caso B')
% tic
% [B.resultados_num,B.resultados_an] = DRA_Numerico(bombeioB, sinal, fibraB);
% toc
% disp('Caso C')
% tic
% [C.resultados_num,C.resultados_an] = DRA_Numerico(bombeioC, sinal, fibraC);
% toc

%% Gráficos

figure(1)
plot(sinal.lambda,A.resultados_num.Ganho_on_off,'.-',sinal.lambda,A.resultados_an.Ganho_on_off,'.-r')
hold on
plot(a(:,1),a(:,2),'ko')
title('Caso A')
legend('Numérico CW', 'Analítico', 'Experimental');
xlabel('Comprimento de Onda (nm)');
ylabel('Ganho ON/OFF (dB)');
axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(A.resultados_an.Ganho_on_off))])

% figure(2)
% plot(sinal.lambda,B.resultados_num.Ganho_on_off,'.-',sinal.lambda,B.resultados_an.Ganho_on_off,'.-r')
% hold on
% plot(lambda,Ganho_B + 0.23*13*ones(length(lambda),1),'ko')
% title('Caso B')
% legend('Numérico CW', 'Analítico', 'Experimental');
% xlabel('Comprimento de Onda (nm)');
% ylabel('Ganho ON/OFF (dB)');
% axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(B.resultados_an.Ganho_on_off))])
% 
% figure(3)
% plot(sinal.lambda,C.resultados_num.Ganho_on_off,'.-',sinal.lambda,C.resultados_an.Ganho_on_off,'.-r')
% hold on
% plot(lambda,Ganho_C + 0.23*13*ones(length(lambda),1),'ko')
% title('Caso C')
% legend('Numérico CW', 'Analítico', 'Experimental');
% xlabel('Comprimento de Onda (nm)');
% ylabel('Ganho ON/OFF (dB)');
% axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(C.resultados_an.Ganho_on_off))])
