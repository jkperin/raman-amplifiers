function err = minimiza(P,bombeio,sinal)

if any(P < 0.4 | P > 0.45)
    err = Inf;
else


%% Dados do bombeio
bombeio.P = P;

%% Dados da fibra
fibra = PCF(bombeio,sinal);
fibra.L = 1e3;

[ripple,Ganho_Medio,Ganho_on_off_medio,Ganho_on_off] = DRA_Analitico(bombeio,sinal,fibra);

err = ripple;
end