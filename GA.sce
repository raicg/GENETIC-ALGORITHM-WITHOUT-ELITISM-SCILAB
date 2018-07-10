RateMutacao=0.01;//Taxa de Mutação
RateCrossOver=0.60;//Taxa de Cruzamento
Geracoes=40;//Numero de Gerações
MaxX=500;//Valor máximo para X
MinX=-500;//Valor mínimo para X
MaxY=500;//Valor máximo para Y
MinY=-500;//Valor mínimo para Y
PopSize=50;//Tamanho da População

Bits=14;

//função que calcula o valor na função dada
function valor=ValoresFuncao(x, y)
    z=-x.*sin(sqrt(abs(x)))-y.*sin(sqrt(abs(y)));
    x=x/250;
    valor=x.*z;
endfunction

//função que normaliza os valores da função para que o menor valor seja a maior fitness
function normalizado=Normalizar(valorFuncao)
    max_old=max(valorFuncao);
    normalizado=valorFuncao-max_old; //subtrai todos os valores pelo maior, para que todos se tornem negativos
    normalizado=(-1)*normalizado; //o mais negativo é o que possui melhor fitness ao inverter o sinal
    normalizado=(normalizado+0.1*max_old)*10; //+0.1*max_old para não ter fitness 0, *10 para aumentar a diferença do melhor para o pior 
endfunction

//função que cria a população inicial
function[popX, popY]=CriarPopulacao(PopSize)
    popX=int(2*rand(PopSize,Bits));//em binário
    popY=int(2*rand(PopSize,Bits));//em binário
endfunction

//função do cross over
function[filho1, filho2]=CrossOver(pai1, pai2, RateCrossOver)
    if rand()<=RateCrossOver then
        corte=int((Bits-1)*rand())+1;
        pai1bin=pai1;
        pai2bin=pai2;
        //cortando o pai 1
        aux1=pai1bin(1:corte)
        aux2=pai1bin((corte+1):Bits)
        //cortando o pai 2
        aux3=pai2bin(1:corte)
        aux4=pai2bin((corte+1):Bits)
        //Unindo os cortes formando os novos filhos
        filho1=[aux1 aux4]
        filho2=[aux3 aux2]
    else
    filho1=pai1;
    filho2=pai2;
    end
endfunction

//função da mutação               
function novofilho=Mutacao(filho, RateMutacao)
    for k=1:Bits
        if rand()<=RateMutacao then //se acontecer a mutação, inverte o bit na posição K
            aux=k;
            invertido=filho(aux); 
            if invertido==0 then
                invertido=1;
            else
                invertido=0;
            end
            filho(aux)=invertido;
        end
    end
    novofilho=filho
endfunction


function [Menores_Valores, Media_Valores]=GA(PopSize, RateMutacao, RateCrossOver, Geracoes, MinX, MinY, MaxX, MaxY)
    [Xb,Yb]=CriarPopulacao(PopSize)
    for k=1:1:Geracoes
        for n=1:1:PopSize
            X(k,n)=MinX+(MaxX-MinX)*(sum(Xb(n,:).*(2.0.^[(Bits-1):-1:0]))/(2^Bits-1)); //transforma o número binário X em decimal entre o seu MAX e MIN
            Y(k,n)=MinY+(MaxX-MinX)*(sum(Yb(n,:).*(2.0.^[(Bits-1):-1:0]))/(2^Bits-1)); //transforma o número binário Y em decimal entre o seu MAX e MIN
        end

        Valores(k,:)=ValoresFuncao(X(k,:),Y(k,:)); //adquire os valores da população de X e Y na função 
        fitness=Normalizar(Valores(k,:)); //normaliza os valores para que os menores se tornem os maiores e todos sejam positivos

        soma=sum(fitness); //soma todos os normalizados
        if soma==0 then //soma não pode ser 0
            disp(soma);
            clf;
        end
        probabilidadesRoleta=fitness/soma; //a soma de todas as probabilidades das fitness deve ser 1
        for l=2:PopSize //faz com que cada probabilidade tenha um intervalo diferente, todas entre 0 e 1
            probabilidadesRoleta(l) = probabilidadesRoleta(l-1) + probabilidadesRoleta(l);
        end
        cross_Over_X=[];
        cross_Over_Y=[];

        //----------------------Roleta ----------------------------
        for i=1:1:PopSize
            escolhido=rand();
            posicao=0;
            acabou=0;
            while acabou==0 //enquanto não achar quem foi o escolhido
                posicao=posicao+1; //posição do escolhido
                if escolhido<=probabilidadesRoleta(posicao); //se achar o escolhido
                    acabou=1; 
                end
            end
            cross_Over_X(i,:)=Xb(posicao,:); //os pais X escolhidos
            cross_Over_Y(i,:)=Yb(posicao,:); //os pais Y escolhidos
        end


        // -------------------Cross Over------------------------------
        for j=1:2:PopSize
            //Cross Over em X
            Pai1_X=cross_Over_X(j,:);
            Pai2_X=cross_Over_X(j+1,:);
            [Filho1_X,Filho2_X]=CrossOver(Pai1_X,Pai2_X,RateCrossOver); //faz o cross over dos 2 pais X gerando 2 filhos X
            Xb(j,:)=Mutacao(Filho1_X,RateMutacao); //verifica a probabilidade de mutação do filho 1 X e, se for menor que o RateMutação, a faz e o inclui em Xb
            Xb(j+1,:)=Mutacao(Filho2_X,RateMutacao); //verifica a probabilidade de mutação do filho 2 X  e, se for menor que o RateMutação, a faz e o inclui em Xb

            //CrossOver em Y
            Pai1_Y=cross_Over_Y(j,:);
            Pai2_Y=cross_Over_Y(j+1,:);
            [Filho1_Y,Filho2_Y]=CrossOver(Pai1_Y,Pai2_Y,RateCrossOver); //faz o cross over dos 2 pais Y gerando 2 filhos Y
            Yb(j,:)=Mutacao(Filho1_Y,RateMutacao);  //verifica a probabilidade de mutação do filho 1 Y  e, se for menor que o RateMutação, a faz e o inclui em Yb
            Yb(j+1,:)=Mutacao(Filho2_Y,RateMutacao);  //verifica a probabilidade de mutação do filho 2 Y  e, se for menor que o RateMutação, a faz e o inclui em Yb
        end

        Menor_Valor=min(Valores(k,:)); //adquire o menor valor na função da atual população
        Menores_Valores(k)=Menor_Valor; //guarda
        media=sum(Valores(k,:))/PopSize;  //adquire a media dos valores na função da atual população
        Media_Valores(k)=media; //guarda
    end
    xtitle("G.A.");
    plot(1:Geracoes,Menores_Valores,'b');
    plot(1:Geracoes,Media_Valores,'r');
    xlabel('Geraçoes');
    ylabel('Valor Mínimo por geração');
    legend('MENOR VALOR','MEDIA POPULAÇÃO')
    xgrid();
    //disp(Menor_Valor);
endfunction

testes = 1;
for i=1:1:testes
    [Menores_Valores2, Media_Valores2]=GA(PopSize, RateMutacao, RateCrossOver, Geracoes, MinX, MinY, MaxX, MaxY);
    Menor(i) = Menores_Valores2(Geracoes);
    disp(Menor(i));    
end

//plot2d2(Menor);
//xtitle("MENOR VALOR POR TENTATIVA");
//xlabel('Tentativa');
//ylabel('Menor valor por tentativa');
//square(0, -1500, testes, -600)
disp(min(Menor));
