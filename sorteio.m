clear, clc, close all

addpath fibras\ -end  
addpath Analitico\ -end
addpath Numerico_CW\ -end

dBm2Watt = @(x) 1e-3*10^(x/10);

N = 150e3;

%% Dados do sinal - 20 canais WDM na banda C
sinal.lambda = linspace(1530.33, 1560.61, 20);
sinal.N_Sinais = length(sinal.lambda);
sinal.Bws = 0.1;
sinal.P = dBm2Watt(-10);

bombeio.N_Bombeios = 4;
bombeio.Bws = 0.1000;
bombeio.FPL = 2;

Pmin = 100; % mW
Pmax = 450; % mW

lmin = 1410; %nm
lmax = 1460; %nm

lambdas = zeros(4, N);
P = zeros(4, N);
ripple = zeros(1,N);
Ganho_medio = zeros(1,N);

Aguarde = waitbar(0, '', 'Name', 'Aguarde...', 'Visible', 'on');
kk = 0;

for k = 1:N
    bombeio.lambda = lmin + rand(1,bombeio.N_Bombeios)*(lmax - lmin);
    bombeio.P = 1e-3*(Pmin + rand(1,bombeio.N_Bombeios)*(Pmax - Pmin));
    
    lambdas(:,k) = bombeio.lambda;
    P(:,k) = bombeio.P;
    
    fibra = PCF(bombeio, sinal);
    fibra.L = 1e3;
    
    [ripple(k),Ganho_medio(k),~, ~] = DRA_Analitico(bombeio,sinal,fibra);
    
    if kk == 999 || k == N
        plot(ripple(k-kk:k),Ganho_medio(k-kk:k),'o')
        hold on
        kk = -1;
    end
   
    waitbar(k/N, Aguarde, sprintf('Calculando Combinação %d de %.0f', k, N))
    kk = kk + 1;
end

close(Aguarde)

xlabel('Ripple [dB]', 'FontSize', 14);
ylabel('Ganho Médio [dB]', 'FontSize', 14);
set(gca, 'FontSize', 14);
saveas(gca, 'Sorteio_150e3_4pump', 'fig')
saveas(gca, 'Sorteio_150e3_4pump', 'bmp')
save Sorteio_150e3_4pump
 