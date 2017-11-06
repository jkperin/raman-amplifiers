%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                         sinal =           
%    N_Bombeios: int(m)                     lambda: [1xn] 
%        lambda: [1xm]                           P: [1xn]
%             P: [1xm]                    N_Sinais: int(n)
%           Bwp: double                        Bws: double
%           FPL: double                        
%Sa�da:
% fibra =
%           L: double
%       alfas: [1xn]
%       alfap: [1xm]
%   alfasdBkm: [1xn]
%   alfapdBkm: [1xm]
% Crpicosinal: [1xn]
%  Crpicopump: [1xm]
%    Crnormal: [1xK]
%     sepfreq: [1xK]
%       Aeffs: [1xn]
%       Aeffp: [1xm]       
%          KR: double % Constante que depende do material dopante da fibra
%          NA: double % Abertura num�rica da fibra
%          no: double % �ndice de Refra��o do n�cleo da fibra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fibra = NDSF(bombeio, sinal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               Atenua��o                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sepfreq = 1440:20:1660;
% 
% NDSF = [0.268 0.240 0.220 0.208 0.201 0.193 0.188  0.191 0.192 0.202 ...
%     0.218 0.246];
% 
% polyNDSF = polyfit(sepfreq,NDSF,4);

polyNDSF = [1.916703088963640e-10,-1.182803597107636e-06,0.002740208637416,-2.824781186650773,1.093541489158918e+03];

fibra.alfasdBkm = polyval(polyNDSF, sinal.lambda);
fibra.alfapdBkm = polyval(polyNDSF, bombeio.lambda);

fibra.alfas = fibra.alfasdBkm*log(10)*1e-4;   %(neper/m)
fibra.alfap = fibra.alfapdBkm*log(10)*1e-4;   %(neper/m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Calculo do ganho de Raman                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdap_ref    = 1455;
Crpico_refp    = 0.39;  %(W^-1*km^-1) Coeficiente de Raman de pico 
Crpico_up_refp = Crpico_refp*1e-3;  % (1/Wm) unidade padrao
fibra.Crpicopump     = lambdap_ref*Crpico_up_refp./bombeio.lambda; % Corrige desn�vel no pico entre lambdas
%
lambdas_ref    = 1455;
Crpico_refs    = 0.39;  %(W^-1*km^-1) Coeficiente de Raman de pico 
Crpico_up_refs = Crpico_refs*1e-3;  % (1/Wm) unidade padrao
fibra.Crpicosinal = lambdas_ref*Crpico_up_refs./sinal.lambda;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Pontos da curva de ganho normalizada                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fibra.Crnormal=[0.00000 0.04235 0.08747 0.12658 0.18108 0.22152 0.24789 0.27602 ...
                0.29008 0.31646 0.34810 0.36920 0.40963 0.46159 0.52836 0.63001 ...
                0.71755 0.78734 0.86769 0.90014 0.95566 1.00000 0.98412 0.93179 ...
                0.97850 0.93687 0.52954 0.30158 0.17941 0.18388 0.20439 0.22679 ...
                0.14241 0.09318 0.07560 0.07384 0.07560 0.07736 0.09845 0.14768 ...
                0.17053 0.15999 0.12131 0.07736 0.05098 0.04395 0.03868 0.03165 ...
                0.03165 0.02989 0.03516 0.04571 0.05978 0.06505 0.05450 0.03868 ...
                0.02813 0.02637 0.02813 0.03692 0.03692 0.02637 0.02110 0.01231 ...
                0.00879 0.01055 0.00176 0.00000];
        
fibra.sepfreq= [0.00000000 6.26866E11 1.25373E12 1.88060E12 2.50746E12 3.13433E12 ...
                3.76120E12 4.38806E12 5.01493E12 5.64179E12 6.26866E12 6.89553E12 ...
                7.52239E12 8.14926E12 8.77612E12 9.40299E12 1.00299E13 1.06567E13 ...
                1.12836E13 1.19105E13 1.25373E13 1.31642E13 1.37911E13 1.44179E13 ...
                1.50448E13 1.56717E13 1.62985E13 1.69254E13 1.75522E13 1.81791E13 ...
                1.88060E13 1.94328E13 2.00597E13 2.06866E13 2.13134E13 2.19403E13 ...
                2.25672E13 2.31940E13 2.38209E13 2.44478E13 2.50746E13 2.57015E13 ...
                2.63284E13 2.69552E13 2.75821E13 2.82090E13 2.88358E13 2.94627E13 ...
                3.00896E13 3.07164E13 3.13433E13 3.19702E13 3.25970E13 3.32239E13 ...
                3.38508E13 3.44776E13 3.51045E13 3.57314E13 3.63582E13 3.69851E13 ...
                3.76120E13 3.82388E13 3.88657E13 3.94926E13 4.01194E13 4.07463E13 ...
                4.13732E13 4.2E13]; 
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             �rea Efetiva                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lambda de Refer�ncia = 1550nm
lamb  = 1400:1649;   

Aeff_lambda_ref = 80; %Livro Mohammed Islam, pg 115
Aeff_ref = [0.86858 0.86937 0.87015 0.87094 0.87172 0.87251 0.87330 ...
            0.87408 0.87487 0.87566 0.87644 0.87723 0.87801 0.87880 ...
            0.87972 0.88050 0.88129 0.88208 0.88286 0.88365 0.88443 ...
            0.88535 0.88614 0.88692 0.88771 0.88850 0.88941 0.89020 ...
            0.89099 0.89177 0.89269 0.89347 0.89426 0.89518 0.89596 ...
            0.89675 0.89754 0.89845 0.89924 0.90003 0.90094 0.90173 ...
            0.90252 0.90343 0.90422 0.90514 0.90592 0.90671 0.90763 ...
            0.90841 0.90933 0.91012 0.91103 0.91182 0.91274 0.91352 ...
            0.91444 0.91523 0.91614 0.91693 0.91785 0.91863 0.91955 ...
            0.92034 0.92125 0.92204 0.92296 0.92387 0.92466 0.92558 ...
            0.92636 0.92728 0.92820 0.92898 0.92990 0.93082 0.93160 ...
            0.93252 0.93344 0.93436 0.93514 0.93606 0.93698 0.93776 ...
            0.93868 0.93960 0.94051 0.94143 0.94222 0.94313 0.94405 ...
            0.94497 0.94589 0.94667 0.94759 0.94851 0.94942 0.95034 ...
            0.95126 0.95218 0.95309 0.95401 0.95493 0.95571 0.95663 ...
            0.95755 0.95846 0.95938 0.96030 0.96122 0.96213 0.96305 ...
            0.96397 0.96488 0.96580 0.96685 0.96777 0.96868 0.96960 ...
            0.97052 0.97144 0.97235 0.97327 0.97419 0.97524 0.97615 ...
            0.97707 0.97799 0.97890 0.97995 0.98087 0.98179 0.98270 ...
            0.98375 0.98467 0.98559 0.98650 0.98755 0.98847 0.98939 ...
            0.99044 0.99135 0.99227 0.99332 0.99423 0.99515 0.99620 ...
            0.99712 0.99817 0.99908 1.00000 1.00105 1.00197 1.00301 ...
            1.00393 1.00498 1.00590 1.00694 1.00786 1.00891 1.00983 ...
            1.01088 1.01192 1.01284 1.01389 1.01481 1.01585 1.01690 ...
            1.01782 1.01887 1.01979 1.02083 1.02188 1.02280 1.02385 ...
            1.02490 1.02594 1.02686 1.02791 1.02896 1.03001 1.03092 ...
            1.03197 1.03302 1.03407 1.03512 1.03616 1.03708 1.03813 ...
            1.03918 1.04023 1.04127 1.04232 1.04337 1.04442 1.04547 ...
            1.04651 1.04756 1.04861 1.04966 1.05071 1.05176 1.05280 ...
            1.05385 1.05490 1.05595 1.05700 1.05805 1.05909 1.06014 ...
            1.06119 1.06224 1.06342 1.06447 1.06551 1.06656 1.06761 ...
            1.06879 1.06984 1.07089 1.07193 1.07311 1.07416 1.07521 ...
            1.07626 1.07744 1.07849 1.07953 1.08071 1.08176 1.08281 ...
            1.08399 1.08504 1.08622 1.08726 1.08844 1.08949 1.09054 ...
            1.09172 1.09277 1.09395 1.09499 1.09617 1.09735 1.09840 ...
            1.09958 1.10063 1.10181 1.10286 1.10404];
 
poly_Aeff = polyfit(lamb, Aeff_ref*Aeff_lambda_ref, 2);   

fibra.Aeffs = 1e-12*polyval(poly_Aeff, sinal.lambda);
fibra.Aeffp = 1e-12*polyval(poly_Aeff, bombeio.lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Caracter�sticas da fibra para o c�lculo do coeficiente de Rayleigh  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fibra.KR = 60;   % (dB/km) constante que depende do material dopante da fibra (Assumido)
fibra.NA = 0.14; % abertura numerica da fibra (Assumido)
fibra.no = 1.468;  %NA^2/(2*deltaindice);  % indice de refracao do nucleo da fibra (Assumido)
