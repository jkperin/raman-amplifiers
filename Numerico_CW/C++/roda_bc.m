clear
clc
load bc_teste

mex numerico_BC.cpp

for i = 1:5
    disp('Compilado Otimizado')
    tic
    for j = 1:200
        res_c = numerico_BC(ya,yb,bombeio,sinal);
    end
    toc
    disp('----------------------------------------------');
    
    disp('Matlab');
    tic
    for k = 1:200
        res_m = numerico_BC_m(ya,yb,bombeio,sinal);
    end
    toc
    disp('----------------------------------------------');
end
