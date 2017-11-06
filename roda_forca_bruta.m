clear, clc, close all

addpath fibras\ -end  
addpath Analitico\ -end
addpath Numerico_CW\ -end

dBm2Watt = @(x) 1e-3*10^(x/10);

%% Bombeios
%lambda(nm)   Pmin(mW)   Pmax(mW)
lambdas = 1410:1460;
Pmin = 100*ones(length(lambdas),1);
Pmax = 450*ones(length(lambdas),1);
Tabela_Bombeios = [lambdas' Pmin Pmax];

Precisao_P = 20; % Variação da potênica dos bombeios (mW) 
bombeio.N_Bombeios = 4;
bombeio.Bws = 0.1000;
bombeio.FPL = 2;

%% Sinais
%% Dados do sinal - 20 canais WDM na banda C
sinal.lambda = linspace(1530.33, 1560.61, 20);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 0.1;
sinal.P = dBm2Watt(-10);

%% fibra
f_fibra  = @(bombeio,sinal) PCF(bombeio,sinal);
L_fibra    = 1e3; %Comprimento da fibra (m)
  

%%
forca_bruta(bombeio, sinal, f_fibra, L_fibra, Tabela_Bombeios,Precisao_P);
 