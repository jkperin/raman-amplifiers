%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                         sinal =           
%    N_Bombeios: int(m)                     lambda: [1xn] 
%        lambda: [1xm]                    N_Sinais: int(n)
%             P: [1xm]                         Bws: double
%           Bwp: double                        
%           FPL: double                        
%Saída:
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
%          NA: double % Abertura numérica da fibra
%          no: double % Índice de Refração do núcleo da fibra
%      beta1s: [1xn]  *
%      beta2s: [1xn]  *
%      beta3s: [1xn]  *
%      beta1p: [1xm]  *
%      beta2p: [1xm]  *
%      beta3p: [1xm]  *
%       gamas: [1xn]  *
%       gamap: [1xm]  *

%
% * Utilizado somente no código de análise de sinais.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fibra = DCF(bombeio,sinal)

alfa_sinal_DCF 		= 0.5; %atenuação dB/km dos sinais
alfa_bombeio_DCF 	= 0.6; %atenuação dB/km dos bombeios

fibra.alfasdBkm = alfa_sinal_DCF*ones(1,sinal.N_Sinais);            %atenuação(dB/km) dos sinais 
fibra.alfapdBkm = alfa_bombeio_DCF*ones(1,bombeio.N_Bombeios);      %atenuação(dB/km)do pump 

fibra.alfas = fibra.alfasdBkm*log(10)*1e-4;     %(neper/m)
fibra.alfap = fibra.alfapdBkm*log(10)*1e-4;     %(neper/m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Pontos da curva de ganho normalizada                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CRPICO_DCF = 3.2;
fibra.Crpicopump = CRPICO_DCF*1e-3*ones(1,bombeio.N_Bombeios);
fibra.Crpicosinal = CRPICO_DCF*1e-3*ones(1,sinal.N_Sinais);

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
%%                             Área Efetiva                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=1530;         %nm
x2=1565;

y1=19.5e-12;      % (micro*m)^2
y2=24e-12;

fibra.Aeffs = y1 + ((sinal.lambda-x1)/(x2-x1))*(y2-y1);     % Curva de Aeff obtida pela ref. do Gaarde
fibra.Aeffp = 15e-12*ones(1,bombeio.N_Bombeios); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Características da fibra para o cálculo do coeficiente de Rayleigh  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fibra.KR= 60;                              % (dB/km) constante que depende do material dopante da fibra
fibra.NA= 0.14;                            % abertura numerica da fibra
fibra.no= 1.468;  %NA^2/(2*deltaindice);    % indice de refracao do nucleo da fibra   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Características da dispersão cromática da fibra                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e5; %km/s

% Dados retirados do datasheet da fibra
lambda0 = 1550;
S0 = -0.3229;                %Slope de fibra (ps/nm^2-km) em lambda0 
D0 = -98;                   %Dispersao em 1550nm (ps/nm-km) em lambda0 
ng = 1.46;                   %Índice de refração efetivo do grupo em lambda0

%Sinais
% Aproximação linear da curva de dispersão por um ponto e inclinação
Ds = (D0+S0*(sinal.lambda-lambda0));   %(ps/nm-km) eq. tirada do ITU-T G652
% Expansão da série de Taylor para beta1
fibra.beta1s = ng/3e8 + 1e-15*D0*(sinal.lambda - lambda0) + 1e-15*S0*(sinal.lambda - lambda0).^2;
% Parâmetro GVD (Eq. 1.2.11 - Agrawal 4ªed)
fibra.beta2s = -((sinal.lambda.^2).*Ds/(2*pi*c))*1e-27;
% GVD de alta ordem (Desenvolver beta3 = dbeta2/dw)
fibra.beta3s = 1e-39*((sinal.lambda.^3)/(2*pi*c)^2).*(2*Ds + sinal.lambda*S0);

%Bombeios
% Aproximação linear da curva de dispersão por um ponto e inclinação
Dp = (D0+S0*(bombeio.lambda-lambda0));   %(ps/nm-km) eq. tirada do ITU-T G652
% Expansão da série de Taylor para beta1
fibra.beta1p = ng/3e8 + 1e-15*D0*(bombeio.lambda - lambda0) + 1e-15*S0*(bombeio.lambda - lambda0).^2;
% Parâmetro GVD (Eq. 1.2.11 - Agrawal 4ªed)
fibra.beta2p = -((bombeio.lambda.^2).*Dp/(2*pi*c))*1e-27;
% GVD de alta ordem (Desenvolver beta3 = dbeta2/dw)
fibra.beta3p = 1e-39*((bombeio.lambda.^3)/(2*pi*c)^2).*(2*Dp + bombeio.lambda*S0);

% Parâmetro não linear Eq. 2.3.29 - Nonlinear fiber optics, 4ªed - Agrawal
n2 = 2.45*10^(-20);
fibra.gamas = 2*pi*n2./(sinal.lambda.*fibra.Aeffs);
fibra.gamap = 2*pi*n2./(bombeio.lambda.*fibra.Aeffp);
