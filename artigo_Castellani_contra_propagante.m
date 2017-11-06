%% Validação de analise de sinais contra propagante baseando-se nos dados
%% do artigo "Design methodology for multi-pumped discrete Raman amplifiers: case-study employing photonic crystal fibers.

clear
clc

addpath Analise_de_sinais\fibras\ -end 
addpath numerico_CW\ -end
addpath Analise_de_sinais\ -end

dBm2Watt = inline('1e-3*10.^(P/10)');

PsinaldBm = [0 -5 -10];
Caso = {'primeiro', 'segundo', 'terceiro'};
Caso2 = {'P = 0dBm','P = -5dBm','P = -10dBm'};

for kk = 1:3
    disp('##############################')
    disp(Caso2{kk});
    disp('##############################')
    % Dados do bombeio
    bombeio.lambda = [1422.4 1451.9];
    bombeio.N_Bombeios = length(bombeio.lambda);
    bombeio.P = [0.2389 0.4352];
    bombeio.Bwp = 0.1;
    bombeio.FPL = 2.053;

    %% Dados do sinal - 47 canais WDM na banda C
    sinal.lambda = [1530.33 1531.9 1533.47 1535.04 1536.61 1538.19 1539.77 1541.35 1542.94 1544.53 1546.12 1547.72 1549.32 1550.92 1552.52 1554.13 1555.75 1557.36 1558.98 1560.61];
    sinal.N_Sinais = length(sinal.lambda);
    sinal.Bws = 0.1;
    sinal.P = dBm2Watt(PsinaldBm(kk));

    % Dados da fibra
    fibraSMF = SMF(bombeio,sinal);
    fibraSMF.L = 70e3;

    fibraPCF = PCF(bombeio,sinal);
    fibraPCF.L = 1e3;

    %% Resultados númericos e analíticos
    tic
    %% Características da transmissão

    Trans.BitRate       = 10e9;                                                % Taxa de bits  [Bits/s]            
    Trans.Nbits         = 64;                                                  % Número de bits (potencia de 2)
    Trans.N             = 2^14;                                                % Número de amostras
    Trans.T0            = 1/Trans.BitRate;                                     % Tempo do bit    
    Trans.Tf            = Trans.T0*Trans.Nbits;                                % Tempo de transmissão

    %% Sinais
    sinal.P             = sinal.P*ones(1,sinal.N_Sinais);                      % Potência dos sinais em z = 0[W]
    sinal.pulso         = 'SuperGaussiano';                                    % Formato dos pulsos dos sinais                                                   
    sinal.modo          = 'semirand';                                          % Modo de geração de sinais 
    for i = 1:sinal.N_Sinais 
        sinal.y(i,:)    = gerando_bits(sinal,Trans,sinal.P(i));                % Gera Sinais 
    end

    sinal_SMF_in = sinal.y;

    %% Bombeios
    bombeio.pulso       = 'SuperGaussiano';                                    % Formato dos pulsos dos bombeios                                                  
    bombeio.modo        = 'CW';                                                % Modo de geração de bombeios
    for i = 1:bombeio.N_Bombeios
        bombeio.y(i,:)  = gerando_bits(bombeio,Trans,bombeio.P(i));            % Gera bombeios
    end

    disp('Propaga somente sinais na fibra SMF');
    sinal_SMF_out = propaga_sb_novo(sinal,fibraSMF,Trans);
    
    sinal.y = sinal_SMF_out;
    sinal.P = (sum(abs(sinal.y).^2,2)/size(sinal.y,2))';

    [resultados_AS,resultados_num, resultados_an] = propaga_contra_propagante(sinal,bombeio,fibraPCF,Trans);
    tempoa = toc;
    stra = myhhmmss(tempoa);
    fprintf('Tempo Gasto = %s\n', stra);

    np = bombeio.N_Bombeios;
    Azs = resultados_AS.Az(np+1:end,:);

    for k=1:sinal.N_Sinais
        Pin(k) = sum(abs(sinal_SMF_out(k,:)).^2)/length(sinal_SMF_out(k,:));


        Pout(k) = sum(abs(Azs(k,:)).^2)/length(Azs(k,:));


        Patenua(k) = fibraPCF.alfasdBkm(k)*fibraPCF.L/1e3;
    end
        PindBm = 10*log10(Pin) + 30;
        PoutdBm = 10*log10(Pout) + 30;
        Ganho_on_off = (PoutdBm - PindBm);      %%% ganho on/off (net gain)
        NG = PoutdBm - PindBm;
        
    load(Caso{kk})
    figure(kk);
    plot(sinal.lambda,Ganho_on_off,'o-k',sinal.lambda,Gof,'o-r');
    set(gca,'FontSize', 14);
    a = [1:4:20 20];
    set(gca, 'xTick', sinal.lambda(a));
    ylabel('Ganho On/Off (dB)', 'FontSize', 14);
    xlabel('Comprimento de Onda (nm)', 'FontSize', 14);
    axis([sinal.lambda(1) sinal.lambda(end) 0 12]);
    title(Caso2{kk})
    legend('Simulado', 'Artigo');
    grid on
    1;
end

