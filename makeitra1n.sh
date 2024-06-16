#!/usr/bin/bash

install(){
	cd yatuner
        cp ../CONFIG_FILES/utils.py yatuner/utils.py
        python3.8 setup.py install
}

install_conf_one(){
	cd yatuner
        cp ../CONFIG_FILES/conf_one/tuner.py yatuner/tuner.py
        python3.8 setup.py install
}

install_conf_two(){
	cd ../yatuner
        cp ../CONFIG_FILES/conf_two/tuner.py yatuner/tuner.py
        python3.8 setup.py install
}

install_conf_three(){
	cd ../yatuner
        cp ../CONFIG_FILES/conf_three/tuner.py yatuner/tuner.py
        python3.8 setup.py install
}

echo "***Iniciando script MASTER para rodar meu fucking TCC!!!!***\n"
#config do pc de testes
apt install -y neofetch
neofetch > neofetch

#instalando os trem para rodar os arquivetos
apt install gcc-10 -y
apt-get -y install linux-tools-generic linux-cloud-tools-generic
apt-get -y install libgsl-dev

#requisitos do programa em python3.8 (env, yatuner, requirements)
apt install python3.8-venv python3.8-dev -y
apt-get install python3.8-pip -y
apt-get install python3.8-setuptools -y
apt install python3.8-dev -y

python3.8 -m venv yatuner_vagrant && 
source yatuner_vagrant/bin/activate &&
pip install setuptools &&
pip install --upgrade setuptools &&
pip install -r requirements.txt 2> /dev/null

mkdir yatuner
mkdir RESULTADOS
git clone https://github.com/JuniMay/yatuner.git yatuner 

chown -R $USER ../tcc_particles_yatuner/yatuner/
chown -R $USER ../tcc_particles_yatuner/yatuner_vagrant/
chown -R $USER CONFIG_FILES/

#cd yatuner
#cp ../CONFIG_FILES/utils.py yatuner/utils.py

# rodando os arquivos para database_one
sysctl kernel.perf_event_paranoid=-1

install_conf_one
for i in {1..10}; do 
        echo "RODANDO SCRIPT 1 -> exec N_${i}"
	cd ../TCC/ && python3.8 stat_time_one.py > ../RESULTADOS/term_output_one_${i}
done

# rodando arquivos para database_two
install_conf_two
for i in {1..10}; do 
        echo "RODANDO SCRIPT 2 -> exec N_${i}"
        cd ../TCC/ && python3.8 stat_time_two.py > ../RESULTADOS/term_output_two_${i}
done

# rodando arquivos para database_three
install_conf_three
for i in {1..10}; do
	echo "RODANDO SCRIPT 3 -> exec N_${i}"
	cd ../TCC/ && python3.8 stat_time_three.py > ../RESULTADOS/term_output_three_${i}
done

# rodando arquivos com o manual
install
for i in {1..10}; do 
        echo "RODANDO SCRIPT SEM PARAMETROS -> exec N_${i}"
        cd ../TCC/ && python3.8 ./stat_no_time.py > ../RESULTADOS/term_no_param_${i}
done

