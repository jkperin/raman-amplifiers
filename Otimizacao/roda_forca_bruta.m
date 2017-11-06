clear
clc

addpath ..\fibras\ -end  
addpath ..\Analitico\ -end
addpath ..\Numerico_CW\ -end


%% Bombeios
%lambda(nm)   Pmin(mW)   Pmax(mW)
lambdas = 1420:1460;
Pmin = 400*ones(length(lambdas),1);
Pmax = 450*ones(length(lambdas),1);
Tabela_Bombeios = [lambdas' Pmin Pmax];

Precisao_P = 10; % Variação da potênica dos bombeios (mW) 
bombeio.N_Bombeios = 2;
bombeio.Bws = 0.2000;
bombeio.FPL = 2;

%% Sinais
% sinal.lambda = [1528.778174,1530.338948,1531.902913,1533.470077,1535.040451,...
%     1536.614044,1538.190867,1539.770930,1541.354242,1542.940813,1544.530654,...
%     1546.123775,1547.720186,1549.319897,1550.922918,1552.529259,1554.138932,...
%     1555.751946,1557.368312,1558.988040,1560.611140,1562.237624,1563.867501,...
%     1565.500783,1567.137480,1568.777603,1570.421163,1572.068170,1573.718635,...
%     1575.372570,1577.029984,1578.690890,1580.355298,1582.023219,1583.694665,...
%     1585.369646,1587.048174,1588.730260,1590.415915,1592.105151,1593.797980,...
%     1595.494412,1597.194459,1598.898133,1600.605446,1602.316408,1604.031033,...
%     1605.749330,1607.471314,1609.196994,1610.926384,1612.659494];

sinal.lambda = linspace(1530.33,1560.61,20);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 0.1;

%% fibra
    % 1 - SMF
    % 2 - SMF2
    % 3 - SMF_Padtec
    % 4 - TW_NZDSF
    % 5 - AllWave
    % 6 - LEAF
    % 7 - Truewave_RS
    % 8 - Corning_NZDSF
    % 9 - NDSF
    % 10 - Corning_DSF
f_fibra  = @(bombeio,sinal) TW_NZDSF(bombeio,sinal);
L_fibra    = 1e3; %Comprimento da fibra (m)
  

%%
forca_bruta(bombeio, sinal, f_fibra, L_fibra, Tabela_Bombeios,Precisao_P);
 