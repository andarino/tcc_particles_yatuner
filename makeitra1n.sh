#!/usr/bin/bash

install(){
	mv ../CONFIG_FILES/tuner.py yatuner/utils
	python setup.py install
}

install_conf_one(){
mv ../CONFIG_FILES/conf_one/tuner.py yatuner/utils

python setup.py install
}

install_conf_two(){
mv ../CONFIG_FILES/conf_two/tuner.py yatuner/utils

python setup.py install
}

install_conf_three(){
mv ../CONFIG_FILES/conf_three/tuner.py yatuner/utils

python setup.py install
}

echo "***Iniciando script MASTER para rodar meu fucking TCC!!!!***\n"

#config do pc de testes
apt install -y neofetch
neofetch > neofetch

#instalando os trem para rodar os arquivetos

apt-get -y install linux-tools-generic linux-cloud-tools-generic
apt-get -y install libgsl-dev

#requisitos do programa em python (env, yatuner, requirements)

python -m venv yatuner_vagrant
source yatuner_vagrant/bin/activate
pip install -r TCC/requirements.txt

git clone https://github.com/JuniMay/yatuner.git
cd yatuner/
	#yatuner/yatuner/utils
mv ../CONFIG_FILES/utils.py yatuner/utils

# rodando os arquivos para database_one
install_conf_one
for i in {1..10}; do 
	echo "RODANDO SCRIPT 1..."
	python ./stat_time_one.py >> term_output_one
done

# rodando arquivos para database_two
install_conf_two
for i in {1..10}; do 
	echo "RODANDO SCRIPT 2..."
       	python ./stat_time_two.py >> term_output_two
done

# rodando arquivos para database_three
install_conf_three
for i in {1..10}; do 
	echo "RODANDO SCRIPT 3..."
       	python ./stat_time_three.py >> term_output_three
done

# rodando arquivos com o manual
install
for i in {1..10}; do 
	echo "RODANDO SCRIPT SEM PARAMENTROS"
       	python ./stat_no_time.py >> term_no_param
done

