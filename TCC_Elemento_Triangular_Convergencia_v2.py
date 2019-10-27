import numpy as np
from matplotlib import pyplot as plt
import math
'''
Esse programa descreve o grafico de convergencia de uma viga atraves do metodo dos elementos finitos. Muitos dos processos
descritos aqui aparecem no 'TCC_Elemento_Triangular_v5', mas são inseridos nesse em uma tentativa de deixa o codigo menos extenso e objetivo.
'''


#O corpo a ser estudado eh uma viga com uma extremidade engastada.
#Descricao da viga:
P = -5 # Força, em KN, na extremidade direita da viga.
L = 200 # Largura da viga em cm.
h0 = 10 # Altura da viga em cm.
M = int(L/h0) # Relacao largura por altura, importante para determinar fatores como quantidade de elementos, nodulos etc...
E = 2.07*10**4 # Modulo de elasticidade da viga em Kn/cm^2
v = 0.3 # Coeficiente de Poisson
t = 2 # Espessura da viga em cm.
I = t*h0**3/12 # Inercia em cm^4.
Deslocamentos = []
Tensoes = []
Precisoes_D = []
Precisoes_T = []
#Elementos:
'''
Um bloco de elemento definido como dois elementos triangulares retangulos isoceles, Sendo seus catetos de dimensao h.
Seu esboco pode ser visto a seguir.
k _______i
 |      /|j
 |  2 /  |
 |  /  1 |
j|/______|
  i       k
Note que, apesar de semelhantes, (1) e (2) sao elementos distintos e devem ter caracteristicas como sua matriz de rigidez([k]) , matriz de
transformacao ([Lmd]) e sua matriz [B] diferentes. Note que as vertices de cada elemento eh rotulada por 'i','j','k'.
'''
n_max = 20     # Relevante para a geracao da malha, "n" eh o numero de subdivicoes da altura do elemento.
              #Mais especificamente, "n" define a altura do elemento e impacta na quantidade de elementos.
              # Quanto maior o "n" maior a precisao, porem, maior o tempo computacional.
'''
Exemplo:
Se n = 2 e h0 = 10cm (Altura da viga) entao temos que um elemento da malha vai ter altura:
h = h0/n = 10/2 = 5 cm
Se a viga tem L = 100cm (Largura da viga) entao sera necessario que a malha tenha:
L/h = 100/5 = 20 blocos de elementos na horizontal. Como n = 2, entao temos que a quantidade de elementos que compoe a malha eh descrita por:
(L/h)*2*n = 20*(numero de blocos)*2(elementos em cada bloco)*2(divisoes da viga) = 80 elementos.
Vale notar que, em geral, o numero de elementos de uma dada malha pode ser calculada como:
(L/h)*2*n = [L/(h0/n)]*2*n = 2*(L/h0)*n^2 = 2*M*n^2.
'''
if 2*M*n_max**2 > 12000: #Avaliacao do numero de elementos.
    print("Atencao, o numero de elementos na malha eh muito grande. Podendo causar problemas de desempenho no computador.")
    print("Numero de elementos: ",2*M*n_max**2)
    Certeza = input("Insira 's' para continuar qualquer outra coisa para cancelar.")
    if Certeza == "s" or Certeza == "S":
        print("Prosseguindo.")
    else:
        exit()
'''
'M' tem que ser um numero inteiro para o programa funcionar, portanto, se M =! L/h0 entao nao eh possivel realizar uma malha corretamente.
'''
if M/(L/h0) == 1: 
    print("Ok")
    print("Relacao L/h0 = ", M)
    if M< 5 :
        print("A relacao L/h0 escolhida esta baixa.")
else:
    '''
    Se a relacao nao der um numero inteiro os elementos não emglobam todo o corpo.
    '''
    print("Cuidado, a relacao L/h0 nao eh um numero inteiro e, consequentemente, nao apresenta resultados confiaveis.") 

for n in range(1,n_max+1):                                                                                                                        
    print("Numero de elementos da malha com n=",n,":",2*M*n**2)
    #Descricao dos elementos:
    h = (h0/n) # Altura e largura de cada elemento.
    # Dimensoes dos dois elementos em um bloco qualquer.
    #Note que as variaveis a baixo sao tupler onde sua primeira posicao pertence ao elemento tipo (1) e a segunda posicao ao elemento tipo (2).
    Xi = (0,h) #posicao no eixo X na vertice i do elemento.
    Yi = (0,h) #posicao no eixo Y na vertice i do elemento.
    Xj = (h,0) #posicao no eixo X na vertice j do elemento.
    Yj = (h,0) #posicao no eixo Y na vertice j do elemento.
    Xk = (h,0) #posicao no eixo X na vertice k do elemento.
    Yk = (0,h) #posicao no eixo Y na vertice k do elemento.
    A = h**2/2 # Area do elemento

        
    #Matriz deslocamento-deformacao [B]. Capitulo 10.2, pag 381 do livro (Rao,2018); Capitulo 5.3.3 ,pag 143 (Alves,2014).
    B = []
    for e in range(2):
        #Matriz [B] no elemento 'e'.
        Bs = (1/(A*2))*np.array([[Yk[e]-Yj[e],0,-(Yk[e]-Yi[e]),0,Yj[e]-Yi[e],0],
                        [0,-(Xk[e]-Xj[e]),0,Xk[e]-Xi[e],0,-(Xj[e]-Xi[e])],
                        [-(Xk[e]-Xj[e]),Yk[e]-Yj[e],Xk[e]-Xi[e],-(Yk[e]-Yi[e]),-(Xj[e]-Xi[e]),Yj[e]-Yi[e]]])
        B.append(Bs) #Lista contendo as duas matrizes [B]
    #Inversa da matriz de coeficientes elastico, matriz [D]. Capitulo 10.2, pag 381 do livro (Rao,2018).
    D = E/(1-v**2)*np.array([[1,v,0],
                              [v,1,0],
                        [0,0,(1-v)/2]])
    '''
    Matriz de rigidez local [k] e global [K]. Capitulo 10.2, pag 381 do livro (Rao,2018).
     [k] = integral([B].T * [D] * [B])Dv = t*A*([B].T * [D] * [B])
    '''
    k = []
    for e in range(2):
        
        k.append(t*A*np.dot(np.dot(B[e].T,D),B[e]))

    '''
    Com a matriz global ainda eh necessario fazer a ordenação dos 
    nodulo para associa-los com cada [K]. Sendo que, a enumeracao dos nodulo eh feita da esquerda para direita, de baixo para cima.

    Exemplo : Se n = 1 e M = 1, para o elemento 1:
    pontas  i j k     
    nó      1 4 3
    (2)_______(4)
      |      /|
      |  2 /  |
      |  /  1 |
      |/______|
    (1)      (3)
    '''
    Ordem = []
    Ordem_K = []
    elem = 0
    for J in range(M*n):
        a = 0
        b = 0
        for Col in range(1,2*n+1):
            elem += 1
            if elem % 2 != 0: # Os elementos impares são diferentes dos pares.
                a += 1
                Ordem.append([J*n+J+a,(J+1)*n+J+2+a,(J+1)*n+J+1+a]) #Criando lista que orienta a posicao de cada nodulo em cada elemento impar.
                Ordem_K.append([2*Ordem[elem-1][0]-1,2*Ordem[elem-1][0],2*Ordem[elem-1][1]-1,2*Ordem[elem-1][1],2*Ordem[elem-1][2]-1,2*Ordem[elem-1][2]]) #Usando lista acima para orientar as matrizes [K] impar. 
            else:
                b += 1
                Ordem.append([(J+1)*n+J+2+b,J*n+J+b,J*n+J+1+b])#Criando lista que orienta a posicao de cada nodulo em cada elemento par.
                Ordem_K.append([2*Ordem[elem-1][0]-1,2*Ordem[elem-1][0],2*Ordem[elem-1][1]-1,2*Ordem[elem-1][1],2*Ordem[elem-1][2]-1,2*Ordem[elem-1][2]]) #Usando lista acima para orientar as matrizes [K] par.

    #Juncao das matrizes de rigidez [_K_].
    # (20*n**2+21*n+1) é o numero de nodulo. Como cada nodulo tem uma deformacao na direcao x e em y, [_K_] tem de ser de dimensao ((M*n*(n+1)+n+1)*2,(n*M*(n+1)+n+1)*2)).
    _K_ = np.zeros(((M*n*(n+1)+n+1)*2,(n*M*(n+1)+n+1)*2))
    for e in range(1,2*M*n**2+1):
    #Usando lista Criada 'Ordem_K' para garantir que cada vertice e coluna de _K_ dos elementos impares esteja no local certo.
        if e % 2 !=0:
            for i in range(6):
                for j in range(6):
                    _K_[Ordem_K[e-1][i]-1][Ordem_K[e-1][j]-1] += k[0][i][j] 
    #Usando lista Criada 'Ordem_K' para garantir que cada vertice e coluna de _K_ dos elementos pares esteja no local certo.
        else:
            for i in range(6):
                for j in range(6):
                    _K_[Ordem_K[e-1][i]-1][Ordem_K[e-1][j]-1] += k[1][i][j]    
          
    #Vetor de forças F e Vetor de deslocamento W
    '''
    Note que a formula que descreve o numero de nodulo da malha eh dada por:
    Num_No = M*n*(n+1)+n+1
    Como a cada nodulo se acrescenta 2 graus de liberdade, temos o vetor de forças tem 2*(M*n*(n+1)+n+1) termos.
    '''

    F = np.zeros(((M*n*(n+1)+n+1)*2,1)) #Criando vetor forca com o numero certo de termos.
    F[-1][0] += P #Aplicando carga pontual em sua extremidade.
    # Definindo os deslocamentos na extremidade esquerda como zero.
    for num in range(2*(n+1)):
        F = np.delete(F,1,0) #Retirando os termos irrelevantes para o calculo.
    # Excluir condicoes de contorno do problema.
    for i in range(2*(n+1)):
        _K_ = np.delete(_K_,0,0) # Excluindo as linhas com deslocamento nulo de _K_
        _K_ = np.delete(_K_,0,1) # Excluindo as colunas com deslocamento nulo de _K_
    # Calculo do deslocamento global [_K_] * W = F.
    Q_desordenado = np.linalg.solve(_K_,F) #Resolvendo a equação linear, descobrindo os deslocamento dos elementos em coordenadas globais, [Q].
    #Vale notar que Q_desordenado nao apresenta os zeros da condicao de contorno, e esta ordenado na mesma sequencia que os nodulos. 
    O = np.zeros((2*(n+1),1)) # Matriz com zeros para acrescentar de volta os deslocamentos zero do engaste.
    Q_desordenado = np.append(O,Q_desordenado,0) # Acrescimo dos zeros no deslocamento.
    # Organizar o deslocamento por elemento.
    '''
    Para melhor explicar Q_desordenado e Q, segue o exemplo com n = 1 e M = 1:
    Q_desordenado = {deslocamento no nodulo 1 em x}                  Q = {deslocamento no nodulo 1 em x}  
                    {deslocamento no nodulo 1 em y}                      {deslocamento no nodulo 1 em y}  
                    {deslocamento no nodulo 2 em x}                      {deslocamento no nodulo 4 em x}  
                              .                                          {deslocamento no nodulo 4 em y}            
                              .                                                        .
                              .                                                        .
                    {deslocamento no nodulo 8 em y}
    '''
    Q = []
    Q_1 = []

    for e in range(1,2*M*n**2+1):
        for i in range(6):
            Q_1.append(Q_desordenado[Ordem_K[e-1][i]-1][0]) #Organizando [Q] na ordem dos elementos.
        Q.append([Q_1])
        Q_1 = []
    Q = np.asarray(Q, dtype=np.float32) #Transformando [Q] de List para um numpy array.
    
    #Calculo dos deslocamento e Tensões, e a guardando essa informação para a geração do grafico.
    Desl_comparativo = 100*Q[-1].T[1]/((-P/6*L**3+P/2*L**3)/(E*I))
    Precisoes_D.append(Desl_comparativo)
    Deslocamentos.append(Q[-1].T[1])
    k_m = int((math.ceil(L/h))/2-1)
    D_B = np.dot(D,B[0])
    T1 = np.dot(D_B,Q[2*k_m*n].T) # Calculo da Tensão no meio da viga no ínicio da seção.Capitulo 10.3, pag 392 do livro (Rao, 2018). 
    D_B = np.dot(D,B[1])
    T2 = (np.dot(D_B,Q[2*k_m*n+1].T)) 
    T = (T1+T2)[0]/2
    Tensoes.append(T)
    Tens_comparativo = T/(h0/2*P*(L-math.ceil(L/h)/2*h)/I)*100
    Precisoes_T.append(Tens_comparativo)
    
'''
geracao dos graficos: 
'''

print("Gráfico de relação (Valor obtido)/(Valor real) para a Tensão e o deslocamento.")
     
x = [n for n in range(n_max)]
plt.xlabel('Numero de subdivisões "n".')
plt.ylabel("(Valor obtido)/(Valor Real) em porcentagem.")
plt.plot(x,Precisoes_D,label='Deslocamento y')
plt.plot(x,Precisoes_T,label='Tensão x')
plt.legend()
plt.title('Gráfico da precisão percentual por subdivisões')
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()

print('Gráfico de deslocamento no final da viga')
plt.xlabel('Número de subdivisões "n".')
plt.ylabel("Deslocamento em cm")
plt.plot(x,Deslocamentos)
plt.title('Gráfico do deslocamento por subdivisões')
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()

print('Inclinação (Derivada discreta) do gráfico de deslocamento.')
x_i = [n for n in range(1,n_max)]
plt.xlabel('Número de subdivisões "n".')
plt.ylabel("Modulo das diferenças entre os deslocamento em cm")
Delta_D = [abs(Deslocamentos[x+1]-Deslocamentos[x]) for x in range(n_max-1)]
plt.plot(x_i,Delta_D)
plt.title('Inclinação da função deslocamento')
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()


print('Gráfico da tensão no começo da viga (parte inferior)')
plt.xlabel('Número de subdivisões "n".')
plt.ylabel("Tensão em KN/cm^2")
plt.plot(x,Tensoes)
plt.title('Gráfico da tensão por subdivisões')
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()

print('Inclinação (Derivada discreta) do gráfico de tensões.')
plt.xlabel('Número de subdivisões "n".')
plt.ylabel("Modulo das diferenças entre as tensões em KN/cm^2")
Delta_T = [abs(Tensoes[x+1]-Tensoes[x]) for x in range(n_max-1)]
plt.plot(x_i,Delta_T)
plt.title('Inclinação da função Tensão')
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()

