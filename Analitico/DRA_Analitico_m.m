% Arquivo .m id�ntico a vers�o compilada de DRA_Analitico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                  sinal =                  fibra =
%    N_Bombeios: int(m)          lambda: [1xn]                  L: double 
%        lambda: [1xm]         N_Sinais: int(n)              lamb: [1xX] (NU)
%             P: [1xm]              Bws: double (NU)        alfas: [1xn] 
%           Bwp: double (NU)                                alfap: [1xm]  
%           FPL: double                                 alfasdBkm: [1xn] (NU)
%                                                       alfapdBkm: [1xm] (NU)
%                                                     Crpicosinal: [1xn] (NU)
%                                                      Crpicopump: [1xm]
%                                                        Crnormal: [1xK]
%                                                         sepfreq: [1xK]
%                                                            Aeff: [1xX] (NU)
%
% (NU) -> N�o Utilizado pelo modelo anal�tico, as structs podem ser
% passadas sem esses par�metros
%
%Sa�das:
%               ripple: double %ripple com atenua��o
%          Ganho_Medio: double %Com atenua��o
%   Ganho_on_off_medio: double %Sem atenua��o
%         Ganho_on_off: [1xn]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [ripple,Ganho_Medio,Ganho_on_off_medio,Ganho_on_off] = DRA_Analitico_m(bombeio,sinal,fibra)
    

    alfap = mean(fibra.alfap);
    pontos = length(fibra.Crnormal);

    fp = (2.99792458e8./bombeio.lambda)*10^9;  %Hz
    fs = (2.99792458e8./sinal.lambda)*10^9;  %Hz
    wp=2*pi*fp;
    sinal.N_Sinais = length(sinal.lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Crnovo entre pumps e sinais                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltafps = zeros(bombeio.N_Bombeios,sinal.N_Sinais);
    Crnovops = zeros(bombeio.N_Bombeios,sinal.N_Sinais);
    Cr_upps = zeros(bombeio.N_Bombeios,sinal.N_Sinais);
    for k=1:sinal.N_Sinais
        for m=1:bombeio.N_Bombeios
            deltafps(m,k)=abs(fp(m)-fs(k));
             for i=1:pontos-1
                 if (fibra.sepfreq(i)<=(abs(deltafps(m,k))))&&(fibra.sepfreq(i+1)>=(abs(deltafps(m,k))))
                    Crnovops(m,k)=fibra.Crnormal(i)+(((abs(deltafps(m,k)))-fibra.sepfreq(i))...
                        *(fibra.Crnormal(i+1)-fibra.Crnormal(i))/(fibra.sepfreq(i+1)-fibra.sepfreq(i)));
                    break;
                end
            end
            Cr_upps(m,k)=fibra.Crpicopump(m)*Crnovops(m,k);  %(m/W)
        end
    end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Crnovo entre pumps                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltafpp = zeros(bombeio.N_Bombeios);
    Crnovopp = zeros(bombeio.N_Bombeios);
    Cr_uppp = zeros(bombeio.N_Bombeios);
    for k=1:bombeio.N_Bombeios
        for m=1:bombeio.N_Bombeios
            deltafpp(m,k)=abs(fp(m)-fp(k));
             for i=1:pontos-1
                 if (fibra.sepfreq(i)<=(abs(deltafpp(m,k))))&&(fibra.sepfreq(i+1)>=(abs(deltafpp(m,k))))
                    Crnovopp(m,k)=fibra.Crnormal(i)+(((abs(deltafpp(m,k)))-fibra.sepfreq(i))...
                        *(fibra.Crnormal(i+1)-fibra.Crnormal(i))/(fibra.sepfreq(i+1)-fibra.sepfreq(i)));
                    break;
                end
            end
           Cr_uppp(m,k)= fibra.Crpicopump(m)*Crnovopp(m,k);
       end
    end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               EQUACAO DOS PUMPS COM ITERACAO ENTRE ELES               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      z=linspace(0,fibra.L,100);

%%%%%%%%%%%%%%%%%%%%% PROPAGA�AO DOS BOMBEIOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pp_lp=0;   Pp_gp=0;    
    Np_lp=0;   Np_gp=0;        

    for j=1:bombeio.N_Bombeios
        Pp_p=0;  Pg_p=0;  
        for i=1:bombeio.N_Bombeios
            if i>j
                for k=1:bombeio.N_Bombeios
                    if k>i
                        Pp_lp=Pp_lp+(wp(i)/wp(k))*Cr_uppp(i,k)*bombeio.P(k);
                    end
                    if k<i
                        Pp_gp=Pp_gp+Cr_uppp(i,k)*bombeio.P(k);
                    end
                end
                Pp_g=Pp_gp;
                Pp_l=Pp_lp;
                Pp_p=Pp_p-((wp(j)/wp(i))*Cr_uppp(j,i)*bombeio.P(i).*(1-exp(-1/(alfap*bombeio.FPL)...
                    *(1-exp(-alfap*(fibra.L-z)))*(-Pp_g+Pp_l))))./(-Pp_g+Pp_l);
                  
                Pp_lp=0;   Pp_gp=0;   Pp_g=0;    Pp_l=0;
            end  %end if i>j

            if i<j
                for k=1:bombeio.N_Bombeios
                    if k>i
                        Np_lp=Np_lp+(wp(i)/wp(k))*Cr_uppp(i,k)*bombeio.P(k);
                    end
                    if k<i
                        Np_gp=Np_gp+Cr_uppp(i,k)*bombeio.P(k);
                    end
                end

                Np_g=Np_gp;
                Np_l=Np_lp;
                Pg_p=Pg_p+(Cr_uppp(j,i)*bombeio.P(i).*(1-exp(-1/(alfap*bombeio.FPL)...
                    *(1-exp(-alfap*(fibra.L-z)))*(-Np_g+Np_l))))./(-Np_g+Np_l);

                Np_lp=0;  Np_gp=0;  Np_g=0;  Np_l=0;
            end %end if i<j
        end  %for i=1:bombeio.N_Bombeios

        clear k

         Pp(j,:)=bombeio.P(j).*exp(-alfap.*(fibra.L-z)).*exp(Pp_p+Pg_p);
     end  %for j=1:bombeio.N_bombeios

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      GANHO ANALITICO OBTIDO COM A SEGUNDA ITERA�AO DOS BOMBEIOS       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear j i
    Ganho_sem_RamandB = zeros(sinal.N_Sinais,1);
    k = zeros(sinal.N_Sinais,1);
    GA_sinaldB = zeros(sinal.N_Sinais,1);
    Ganho_on_off = zeros(sinal.N_Sinais,1);

    for j=1:sinal.N_Sinais
         Ganho_pump=0;
         for i=1:bombeio.N_Bombeios
            k(i,j)=(Cr_upps(i,j)/bombeio.FPL);
            Ganho_pump = Ganho_pump+k(i,j)*trapz(z,Pp(i,:));  
         end

         Ganho_sem_RamandB(j)=10*log10(exp(-fibra.alfas(j)*fibra.L));

         GA_sinaldB(j)=10*log10(exp(Ganho_pump - fibra.alfas(j)*fibra.L));

         Ganho_on_off(j)=GA_sinaldB(j) - Ganho_sem_RamandB(j);
    end
    
    GA_sinaldB = GA_sinaldB';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              SAIDA                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ripple                  = max(GA_sinaldB)-min(GA_sinaldB);
    Ganho_Medio             = mean(GA_sinaldB);
    Ganho_on_off_medio      = mean(Ganho_on_off);

end
