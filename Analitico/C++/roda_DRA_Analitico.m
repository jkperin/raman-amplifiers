clear all
clc
load sample
for k = 1:5
    fprintf('-----------compilado-----------\n')
    tic
    for j = 1:200
        [compilado.ripple,compilado.Ganho_Medio,compilado.Ganho_on_off_medio,compilado.GA_sinaldB] = DRA_Analitico(bombeio,sinal,fibra); 
    end
    toc
        
    fprintf('-----------matlab-----------\n');
    tic
    for i = 1:200
        [matlab.ripple,matlab.Ganho_Medio,matlab.Ganho_on_off_medio,matlab.GA_sinaldB] = DRA_Analitico_m(bombeio,sinal,fibra); 
    end
    toc
end

