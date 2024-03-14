# Avaliador de interação de partículas
Estudo sobre a performance de um problema de sistemas dinâmicos, usando yatuner.

## Instalando os requisitos:
* Python 3.6 ou mais recente
* GCC 9 ou mais recente
* Instalar o framework yatuner (https://github.com/JuniMay/yatuner)
* libgsl-dev
* pip3 install requirements.txt 
* Dica: crie um ambiente virtual com python (https://docs.python.org/3/tutorial/venv.html)

## Instruções para execução do programa de simulação

Caso queira executar somente o arquivo de simulação:
> gcc downloading_wall_restitution_GPU.c -o downloading_wall_restitution_GPU -std=c11 -lgsl -lgslcblas -lm -pg -O3

Executar:

> ./downloading_wall_restitution_GPU 0.0001 8 10 0.09
 


* *./downloading_wall_restitution_GPU* (tempo total) (grãos na vertical) (grãos na horizontal) (parâmetro de afilamento)

## Execução do arquivo de Tunning

> python tune_size.py
 
Ao executar o tunning da aplicação um diretório `yatuner-db/` será formado, este diretório contém os parâmetros selecionados pelo framework para melhor performance do algorítmo, além de algumas outras peças para análise.

