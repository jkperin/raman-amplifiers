clear
clc

lambda = 1410:1460;
N = length(lambda);

bombeio.N_Bombeios = 2;
bombeio.Bws = 0.2000;
bombeio.FPL = 2;


%% Dados do sinal
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
sinal.Bws = 1;

options = optimset('MaxFunEvals',10e3,'MaxIter',10e3);
P0 = 0.41*ones(1,bombeio.N_Bombeios);

I_pumps = combnk(1:length(lambda),bombeio.N_Bombeios);
I_pumps_size = size(I_pumps, 1);

for k = 1:I_pumps_size
    bombeio.lambda = lambda(I_pumps(k,:));
	[Popt(k,:), fval(k), exitflag] = fminsearch(@(P)minimiza(P,bombeio,sinal),P0,options);
    switch (exitflag)
        case 1
            disp('fminsearch converged to a solution x.');
        case 0
            disp('Maximum number of function evaluations or iterations was reached');
        case -1
            disp('Algorithm was terminated by the output function');
    end
end

[err_min, ind] = min(fval);
disp('Resultados!!!!!!!!!!!!!!!!!');
bombeio.lambda = lambda(I_pumps(ind,:));
bombeio.P = Popt(ind,:)

fibra = PCF(bombeio,sinal);
fibra.L = 1e3;
[ripple,Ganho_Medio,Ganho_on_off_medio,Ganho_on_off] = DRA_Analitico(bombeio,sinal,fibra)
