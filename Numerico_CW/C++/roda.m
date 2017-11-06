clear
clc

mex numerico_ODE.cpp

load ode_teste

for i = 1:5
    disp('Compilado Otimizado')
    tic
    for j = 1:200
        dydx_c = numerico_ODE(x,y,bombeio,sinal,fibra,Param);
    end
    toc
    disp('----------------------------------------------');
    
    disp('Matlab');
    tic
    for k = 1:200
        dydx_m = numerico_ODE_m(x,y,bombeio,sinal,fibra,Param);
    end
    toc
    disp('----------------------------------------------');
end
