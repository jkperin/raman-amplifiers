%% Validação dos códigos análitico e numérico CW baseando-se nos dados
%% experimentais do artigo "Analysis of a multiple-pump Raman ampli?er" por
%% M. Achtenhagen

clear, clc, close all

%% Dados Experimentais
p0 = [134 622];  % 1520 0

c=[191 312;200 284;203 290;210 288;218 278;226 267;233 265;240 254;253 257;...
    248 249;268 244;261 230;276 220;281 208;288 199;297 196;303 189;308 178;...
    315 170;322 171;330 166;344 157;353 150;358 147;365 157;371 147;379 136;...
    385 128;393 125;403 110;407 114;414 122;421 111;428 92;435 116;442 108;...
    450 111;456 125;463 125;471 122;485 141;499 151;512 166;527 166;541 165;...
    555 171;544 195];

b=[191 486;198 470;205 475;209 460;218 459;226 447;234 449;240 445;246 439;...
    253 458;268 454;261 427;275 429;281 419;288 411;295 425;302 409;309 406;...
    315 410;322 403;328 408;343 412;350 403;357 403;363 409;371 417;378 398;...
    385 397;391 407;398 397;406 394;413 408;419 389;426 387;435 410;440 398;...
    450 417;463 426;456 433;469 436;484 458;498 466;511 470;526 479;541 491;...
    558 498;566 495];

a=[191 517;198 499;204 493;209 503;219 500;225 493;233 490;239 485;248 482;...
    254 491;267 486;260 471;274 466;281 459;288 454;295 459;302 450;309 442;...
    316 438;322 437;329 439;342 429;350 425;358 428;364 437;371 426;379 422;...
    385 420;393 421;400 404;406 413;413 416;420 407;427 402;436 405;442 411;...
    450 410;456 415;464 415;470 406;484 422;499 427;513 429;526 431;541 432;...
    555 435;565 438];

%    Eixo lambda

l1     = [134 223 309 572];
l2     = [1520 1530 1540 1570];
[pl,s] = polyfit(l1,l2,1);
lambda = polyval(pl,c(:,1));

%   Eixo Ganho

G1     = [622 487 354];
G2     = [0 5 10];
[pG,s] = polyfit(G1,G2,1);

Ganho_C = polyval(pG,c(:,2));
Ganho_B = polyval(pG,b(:,2));
Ganho_A = polyval(pG,a(:,2));


%% Simulação

addpath fibras\ -end 
addpath Analitico\ -end

dBm2Watt = inline('1e-3*10.^(P/10)');

%% Dados do sinal - 47 canais WDM na banda C
sinal.lambda = linspace(1526.3,1568.9,47);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 1;
sinal.P = dBm2Watt(-20);

%% Caso A

% Dados do bombeio
bombeioA.lambda = [1443 1465];
bombeioA.N_Bombeios = length(bombeioA.lambda);
bombeioA.P = [0.15 0.3];
bombeioA.Bwp = 1;
bombeioA.FPL = 1;

% Dados da fibra
fibraA = TW_NZDSF(bombeioA,sinal);
fibraA.L = 13e3;

%% Caso B

% Dados do bombeio
bombeioB.lambda = [1443 1465];
bombeioB.N_Bombeios = length(bombeioB.lambda);
bombeioB.P = [0.3 0.15];
bombeioB.Bwp = 1;
bombeioB.FPL = 1;

% Dados da fibra
fibraB = TW_NZDSF(bombeioB,sinal);
fibraB.L = 13e3;

%% Caso C

% Dados do bombeio
bombeioC.lambda = [1443 1455 1465];
bombeioC.N_Bombeios = length(bombeioC.lambda);
bombeioC.P = [0.3 0.3 0.3];
bombeioC.Bwp = 1;
bombeioC.FPL = 1;

% Dados da fibra
fibraC = TW_NZDSF(bombeioC,sinal);
fibraC.L = 13e3;

%% Resultados númericos e analíticos
disp('Caso A')
tic
[A.ripple,A.Ganho_Medio,A.Ganho_on_off_medio,A.Ganho_on_off, A.Ppump0] = DRA_Analitico_num(bombeioA,sinal,fibraA);
[Am.ripple,Am.Ganho_Medio,Am.Ganho_on_off_medio,Am.Ganho_on_off, Am.Ppump0] = DRA_Analitico_num_m(bombeioA,sinal,fibraA);
toc
disp('Caso B')
tic
[B.ripple,B.Ganho_Medio,B.Ganho_on_off_medio,B.Ganho_on_off, B.Ppump0] = DRA_Analitico_num(bombeioB,sinal,fibraB);
[Bm.ripple,Bm.Ganho_Medio,Bm.Ganho_on_off_medio,Bm.Ganho_on_off, Bm.Ppump0] = DRA_Analitico_num_m(bombeioB,sinal,fibraB);
toc
disp('Caso C')
tic
[C.ripple,C.Ganho_Medio,C.Ganho_on_off_medio,C.Ganho_on_off, C.Ppump0] = DRA_Analitico_num(bombeioC,sinal,fibraC);
[Cm.ripple,Cm.Ganho_Medio,Cm.Ganho_on_off_medio,Cm.Ganho_on_off, Cm.Ppump0] = DRA_Analitico_num_m(bombeioC,sinal,fibraC);
toc

%% Gráficos

figure(1)
plot(sinal.lambda,A.Ganho_on_off,'.-', sinal.lambda,Am.Ganho_on_off,'.r')
hold on
plot(lambda,Ganho_A + 0.23*13*ones(length(lambda),1),'ko')
title('Caso A')
legend('Analitico', 'Analítico M', 'Experimental');
xlabel('Comprimento de Onda (nm)');
ylabel('Ganho ON/OFF (dB)');
axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(A.Ganho_on_off))])

figure(2)
plot(sinal.lambda,B.Ganho_on_off,'.-', sinal.lambda,Bm.Ganho_on_off,'.r')
hold on
plot(lambda,Ganho_B + 0.23*13*ones(length(lambda),1),'ko')
title('Caso B')
legend('Analitico', 'Analítico M', 'Experimental');
xlabel('Comprimento de Onda (nm)');
ylabel('Ganho ON/OFF (dB)');
axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(B.Ganho_on_off))])

figure(3)
plot(sinal.lambda,C.Ganho_on_off,'.-', sinal.lambda,Cm.Ganho_on_off,'.r')
hold on
plot(lambda,Ganho_C + 0.23*13*ones(length(lambda),1),'ko')
title('Caso C')
legend('Analitico', 'Analítico M', 'Experimental');
xlabel('Comprimento de Onda (nm)');
ylabel('Ganho ON/OFF (dB)');
axis([sinal.lambda(1) sinal.lambda(end) 0 1.2*max(max(C.Ganho_on_off))])
