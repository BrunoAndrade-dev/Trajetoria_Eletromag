# Simulação da força de Lorentz

Projeto didático em Python: uma partícula carregada em campos elétrico e magnético
uniformes e constantes. Compara integração numérica com solução analítica,
verifica energia e produz gráficos, dados e animação 3D.

## Como executar

Na pasta do projeto:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python simulacao.py
```

São abertas duas figuras: painel de diagnóstico e animação. Em Linux, a abertura
de janelas requer um backend gráfico do Matplotlib, como Tk ou Qt. Para executar
em um ambiente sem interface gráfica, use:

```bash
python simulacao.py --no-show
python simulacao.py --no-show --gif
```

Resultados na pasta `output/`:

- `diagnosticos.png`: trajetória 3D, projeção XY, posições, rapidez, energia e erro.
- `trajetoria.csv`: tempo, posição, velocidade, energias por massa e erro.
- `resumo.json`: parâmetros, período, deriva e métricas de precisão.
- `animacao.gif`: animação, quando solicitado com `--gif`.

O GIF tem até 360 quadros e 30 quadros por segundo. O tempo de reprodução é
independente do tempo físico, exibido em cada quadro. A integração e o CSV
preservam todas as amostras, mesmo quando a animação usa menos quadros.

## Experimentos

```bash
python simulacao.py --cenario livre --saida output/livre
python simulacao.py --cenario eletrico --saida output/eletrico
python simulacao.py --cenario magnetico --saida output/magnetico
python simulacao.py --cenario cruzados --saida output/cruzados
python simulacao.py --q-m -1 --saida output/carga_negativa
python simulacao.py --E 0 0 0.2 --B 0 0 1 --saida output/paralelos
python simulacao.py --B 0 0 2 --duracao 15 --metodo DOP853
python -m unittest -v
```

`--E` e `--B` substituem os campos do cenário. Use `--v0`, `--r0`, `--q-m`,
`--duracao` e `--pontos` para modificar as condições. Unidades sempre no SI:
metros, segundos, m/s, V/m, teslas e C/kg. `--pontos` é um mínimo; o programa
aumenta a amostragem para representar cada volta com pelo menos 80 pontos.
Existe um limite de 200 mil amostras para evitar uma execução acidental enorme.

## A teoria como uma história

Imagine uma bolinha muito pequena que tem uma propriedade chamada carga
elétrica. A bolinha pode ter carga positiva ou negativa. Ela também tem massa:
quanto maior a massa, mais difícil é mudar seu movimento.

Um campo é como uma regra invisível presente em cada ponto do espaço. O campo
elétrico diz como empurrar uma carga; o campo magnético diz como desviar uma
carga que já está andando. No nosso modelo, essas regras são iguais em todos
os lugares e não mudam com o tempo.

O campo elétrico funciona como um empurrão numa direção fixa. Uma carga positiva
recebe o empurrão no sentido do campo; uma negativa, no sentido contrário.
Ele pode acelerar ou frear a bolinha e, portanto, alterar sua energia cinética.

O campo magnético funciona como um empurrão de lado, sempre perpendicular à
velocidade. Sozinho, ele muda a direção sem mudar a rapidez. Lembra uma bolinha
presa a uma corda girando: a corda puxa para o centro. A analogia é apenas para
a direção da força; não existe corda nem atrito na simulação.

Se a bolinha estiver parada, o campo magnético sozinho não a coloca em
movimento. Se ela estiver andando exatamente na direção do campo magnético,
também não haverá força magnética. O desvio aparece na parte do movimento
perpendicular ao campo.

Quando existe movimento perpendicular e também paralelo ao campo, a trajetória
vira uma hélice, como uma mola: a bolinha gira e ao mesmo tempo avança.

## A regra matemática

Força de Lorentz:

    F = q (E + v × B)

Segunda lei de Newton:

    a = F/m = (q/m) (E + v × B)

`×` significa produto vetorial, não multiplicação comum. Sua expansão é:

    v × B = (vy Bz - vz By, vz Bx - vx Bz, vx By - vy Bx)

O produto aponta perpendicularmente aos dois vetores. A regra da mão direita
indica seu sentido para carga positiva; carga negativa inverte a força.
O fator `q/m` indica a facilidade de desviar ou acelerar a partícula.

O estado guardado pelo programa é `[x, y, z, vx, vy, vz]`. Sua derivada é
`[vx, vy, vz, ax, ay, az]`. A função `derivadas` calcula essa derivada e o
`solve_ivp` resolve o sistema com passos adaptativos. O RK45 compara fórmulas
de ordens diferentes para estimar o erro local e ajustar o passo.

## O que acontece com os números do relatório

Usamos `q/m = 1`, `E = (0, 0.2, 0)`, `B = (0, 0, 1)`, posição inicial zero
e velocidade inicial `(1, 0, 0.2)`.

No começo, `v × B = (0, -1, 0)`. Somando o campo elétrico, a aceleração inicial
é `(0, -0.8, 0)`: apesar de o campo elétrico apontar para +Y, a força magnética
é maior e o movimento começa curvando para -Y.

A bolinha gira no plano XY, enquanto o centro desse giro se desloca para +X.
Ela também avança em Z com velocidade constante de 0.2 m/s.

    velocidade de deriva = (E × B)/|B|² = (0.2, 0, 0) m/s
    frequência angular ciclotrônica = |q/m| |B| = 1 rad/s
    período = 2π / frequência angular ≈ 6.283 s

No referencial que acompanha a deriva, a velocidade perpendicular inicial vale
`1 - 0.2 = 0.8 m/s`. O raio é `0.8 / 1 = 0.8 m` e permanece constante.
A projeção XY é uma trocoide com laços para estas condições iniciais.

Os ciclos não crescem continuamente, ao contrário da interpretação da tabela
do relatório. A rapidez no laboratório oscila porque o campo elétrico realiza
trabalho positivo ou negativo ao longo da órbita. Depois de uma volta completa,
a rapidez volta ao mesmo valor.

Para esse caso específico, a solução exata em SI é:

    x(t) = 0.2 t + 0.8 sin(t)
    y(t) = 0.8 (cos(t) - 1)
    z(t) = 0.2 t

O argumento trigonométrico representa `omega*t`, com `omega = 1 rad/s`.
Essa fórmula é a régua usada para conferir o cálculo numérico.

## Por que os gráficos provam algo

O painel compara as posições calculadas com as posições exatas. O erro é a
distância entre elas. Um erro pequeno indica boa concordância; não significa
que o computador calculou com precisão infinita.

A energia cinética por massa é `K/m = |v|²/2`. Para o campo elétrico uniforme,
escolhemos a energia potencial por massa `U/m = -(q/m) E·r`, com zero na origem.
O campo magnético não realiza trabalho porque `(v × B)·v = 0`. Assim:

    H/m = |v|²/2 - (q/m) E·r = constante

O gráfico de energia deve mostrar trocas entre cinética e potencial, mantendo
a soma aproximadamente constante. Potencial negativo é permitido: seu zero é
uma escolha de referência. Como informamos somente q/m, usamos J/kg. Para
calcular energias em joules é preciso informar a massa separadamente.

Os testes verificam movimento livre, aceleração elétrica constante, conservação
da rapidez sob campo magnético, período, sinal da carga, partícula neutra,
campos oblíquos e concordância com a solução exata.

## Organização e otimização

- `Parametros`: condições físicas e validação das entradas.
- `derivadas`: lei de movimento em seis componentes.
- `simular`: integração adaptativa; `rtol=1e-8`, `atol=1e-10`.
- `solucao_exata`: referência independente da integração, para campos uniformes.
- `diagnosticos`: energia e erro de posição.
- `criar_painel` e `criar_animacao`: apresentação dos resultados já calculados.
- `main`: opções de linha de comando e exportação.

A integração acontece uma vez. Os quadros apenas reutilizam os dados e atualizam
os objetos existentes. A animação pode ter menos quadros que a integração.
Os três eixos usam a mesma escala espacial para evitar deformações da órbita.
O passo interno é limitado a um vigésimo do período ciclotrônico, além do
controle adaptativo. `t_eval` define amostras de saída, não os passos internos.
As tolerâncias controlam erro local; a comparação analítica mede o erro global
observado. Para outras escalas físicas, revise tolerâncias e unidades ou
adimensionalize o problema.

## Limites do modelo

Uma partícula de teste, velocidades muito menores que a velocidade da luz,
campos constantes e uniformes, sem colisões, radiação, gravidade, interação
entre partículas ou reação da partícula sobre os campos. Não simula um tokamak
completo nem a aceleração de um cíclotron com campo elétrico alternado.
Os parâmetros padrão são didáticos e não representam um elétron real.

Para campos não uniformes ou dependentes do tempo, altere a função das forças,
substitua a referência analítica e revise a lei de energia usada na validação.

Referência da API: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
