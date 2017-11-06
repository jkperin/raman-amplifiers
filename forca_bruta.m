function forca_bruta(bombeio, sinal, f_fibra,L_fibra, Tabela_Bombeios, Precisao_P)
figure(1)
%Gera todas combinações possíveis entre os comprimentos de onda
%combnk(v,k) returns all combinations of the n elements in v 
%taken k at a time.
I_pumps = combnk(1:length(Tabela_Bombeios(:,1)),bombeio.N_Bombeios);
I_pumps_size = size(I_pumps, 1);  

Bombeios_P_max = Tabela_Bombeios(:,3);
Bombeios_P_min = Tabela_Bombeios(:,2);
delta_P = Precisao_P;

%Número de combinações de potências para cada bombeio
P_N_comb = floor((Bombeios_P_max - Bombeios_P_min)/delta_P) + 1; 

%Número total de combinações comprimento de onda/potência possível
Total_comb = 0;
for i = 1:size(I_pumps,1)
    Total_comb = Total_comb + det(diag(P_N_comb(I_pumps(i,:))));
end
fprintf('Combinações = %d\n', Total_comb);
if Total_comb > 1e10
    fprintf('Número de combinações é muito grande!!!!\n');
    return
end
%% Força Bruta   
    count = 1;
    sair = false;

    MAX_IT = 500; %Máximo número de iterações para estimar o tempo total    
    
    t_amostra = tic;
    Aguarde = waitbar(count,'', 'Name', 'Aguarde...', 'Visible', 'off');
    if Total_comb > MAX_IT %Para muitas combinações amostra tempo
        for j = 1:I_pumps_size %Percorre combinações de bombeios 
            comb_Pot = Tabela_Bombeios(I_pumps(j,1),2):Precisao_P:Tabela_Bombeios(I_pumps(j,1),3);
            for i = 2:bombeio.N_Bombeios
                comb_Pot = combvec(comb_Pot, Tabela_Bombeios(I_pumps(j,i),2):...
                    Precisao_P:Tabela_Bombeios(I_pumps(j,i),3));
            end
            bombeio.lambda = Tabela_Bombeios(I_pumps(j,:),1)';
            
            tmp_ripple = zeros(1,size(comb_Pot,2));
            tmp_Ganho_on_off_medio = zeros(1,size(comb_Pot,2));
            for k = 1:size(comb_Pot,2)
                bombeio.P = 1e-3*comb_Pot(:,k)';
                fibra = f_fibra(bombeio,sinal);
                fibra.L = L_fibra;

                [tmp_ripple(k),tmp_Ganho_Medio,tmp_Ganho_on_off_medio(k),...
                    tmp_Ganho_on_off] = DRA_Analitico(bombeio,sinal,fibra);
                             
                waitbar(count/Total_comb, Aguarde, sprintf('Calculando Combinação %d de %.0f', count, Total_comb))
                count = count + 1;
                if count == MAX_IT
                    amostra_tempo = toc(t_amostra);
                    c = ceil(amostra_tempo*Total_comb/MAX_IT);
                    tempo_estimado_str = hhmmss(c); %Converte tempo em s para formato hh:mm:ss
                    if(c > 5) % Se gastar pouco tempo(5s) nem pergunta
                        alerta = ['O tempo estimado de solução é de ' tempo_estimado_str '. Deseja continuar?'];         
                        choice = questdlg(alerta, 'Atenção!!', 'Sim', 'Não', 'Não');
                        tic;
                        if strcmp(choice,'Não')
                            sair = true;
                            break;
                        else
                            set(Aguarde,'Visible', 'on');
                        end
                    end
                end
            end
            
            plot(tmp_ripple,tmp_Ganho_on_off_medio,'o')
            hold on
            
            if sair
                break;
            end
        end
        t_total = toc;
        str = hhmmss(amostra_tempo + t_total);
        fprintf('Tempo total = %s\n',str);
        close(Aguarde);
    else
        for j = 1:size(I_pumps,1) %Percorre combinações de bombeios 
            comb_Pot = Tabela_Bombeios(I_pumps(j,1),2):Precisao_P:Tabela_Bombeios(I_pumps(j,1),3);
            for i = 2:bombeio.N_Bombeios
                comb_Pot = combvec(comb_Pot,Tabela_Bombeios(I_pumps(j,i),2):...
                    Precisao_P:Tabela_Bombeios(I_pumps(j,i),3));
            end
            bombeio.lambda = Tabela_Bombeios(I_pumps(j,:),1)';
            tmp_ripple = zeros(1,size(comb_Pot,2));
            tmp_Ganho_on_off_medio = zeros(1,size(comb_Pot,2));
            
            for k = 1:size(comb_Pot,2)
                bombeio.P = 1e-3*comb_Pot(:,k)';
                fibra = f_fibra(bombeio,sinal);
                fibra.L = L_fibra;

                [tmp_ripple(k),tmp_Ganho_Medio,tmp_Ganho_on_off_medio(k),...
                    tmp_Ganho_on_off] = DRA_Analitico(bombeio,sinal,fibra);
            end
            plot(tmp_ripple,tmp_Ganho_on_off_medio,'o')
            hold on
        end
        toc(t_amostra);
    end
    
    if ~sair
        sound(1+square(1:100));
        set(gca,'Visible','on');
        grid('on'),hold('off'),zoom('on')
        xlabel('ripple [dB]'), ylabel('Ganho On/Off Médio [dB]')
    else
        close all
    end

%% Função auxiliar    
function str = hhmmss(c)
    if c < 60 %segundos
        str = sprintf('%.0f segundos', c);
    elseif c < 3600
        mm = floor(c/60);
        ss = mod(c,60);
        str = sprintf('%d min %d seg', mm, ss);
    else
        hh = floor(c/3600);
        mm = floor(mod(c,3600)/60);
        ss = mod(mm, 60);
        str = sprintf('%d h %d min %d seg', hh, mm, ss);
    end