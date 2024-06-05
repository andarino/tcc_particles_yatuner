#!/usr/bin/bash

install(){
        cp ../CONFIG_FILES/utils.py yatuner/utils.py
        python3 setup.py install
}

install_conf_one(){
        cp CONFIG_FILES/conf_one/tuner.py yatuner/tuner.py
        python3 setup.py install
}

install_conf_two(){
        cp CONFIG_FILES/conf_two/tuner.py yatuner/tuner.py
        python3 setup.py install
}

install_conf_three(){
        cp CONFIG_FILES/conf_three/tuner.py yatuner/tuner.py
        python3 setup.py install
}

echo "***Iniciando script MASTER para rodar meu fucking TCC!!!!***\n"
#config do pc de testes
apt install -y neofetch
neofetch > neofetch

#instalando os trem para rodar os arquivetos

apt-get -y install linux-tools-generic linux-cloud-tools-generic
apt-get -y install libgsl-dev

#requisitos do programa em python3 (env, yatuner, requirements)
apt install python3.8-venv -y
apt-get install python3-pip -y
apt-get install python3-setuptools -y
apt install python3.12-dev -y

python3 -m venv yatuner_vagrant && \\
source yatuner_vagrant/bin/activate && \\
pip install setuptools && \\
pip install --upgrade setuptools && \\
pip install -r requirements.txt

#se ha o diretorio tuner -> remove
<<LOW
if test -d "yatuner"; then 
        rm -rf yatuner/
fi
mkdir yatuner
mkdir RESULTADOS
git clone https://github.com/JuniMay/yatuner.git yatuner 

cd yatuner
cp ../CONFIG_FILES/utils.py yatuner/utils.py

# rodando os arquivos para database_one
install_conf_one
for i in {1..10}; do 
        echo "RODANDO SCRIPT 1 -> exec N_${i}"
        python3 ../TCC/stat_time_one.py >> ../RESULTADOS/term_output_one_${i}
done

# rodando arquivos para database_two
install_conf_two
for i in {1..10}; do 
        echo "RODANDO SCRIPT 2..."
        python3 ./stat_time_two.py >> term_output_two
done

# rodando arquivos para database_three
install_conf_three
for i in {1..10}; do 
        echo "RODANDO SCRIPT 3..."
        python3 ./stat_time_three.py >> term_output_three
done

# rodando arquivos com o manual
install
for i in {1..10}; do 
        echo "RODANDO SCRIPT SEM PARAMENTROS"
        python3 ./stat_no_time.py >> term_no_param
done

LOW

    