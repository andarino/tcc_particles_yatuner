
case $1 in
"-h") echo "./$0 [param]"; echo "ex: ./$0 { -std=c11  -lgsl -lgslcblas -lm -pg }"
exit 0
;;
"-v") echo "hilarious version"
exit 0
;;
esac

if [ -z $1 ]; then 
	echo "[!] Informe parametros para compilacao [!]"
	echo "ex: ./autotestes.sh -std=c11 -lgsl -lgslcblas -lm -pg"
	exit 0
fi 

if [ -e downloading_wall_restitution0 ]
then 
	echo "[!] removendo arquivo downloading_wall_restitution1 existente "
	rm downloading_wall_restitution0

elif [ -e downloading_wall_restitution1 ]
then 
	echo "[!] removendo arquivo downloading_wall_restitution1 existente "
	rm downloading_wall_restitution1

elif [ -e downloading_wall_restitution2 ] 
then 
	echo "[!] removendo arquivo downloading_wall_restitution2 existente "
	rm downloading_wall_restitution2

elif [ -e downloading_wall_restitution3 ] 
then 
	echo "[!] removendo arquivo downloading_wall_restitution3 existente "
	rm downloading_wall_restitution3
else 
	echo  "[x] Nenhum executável para ser removido !!!"
fi

echo "[!] compilando -O0, -O1, -O2, -O3 com os parametro: $*"

OPT0=`gcc downloading_wall_restitution.c -o downloading_wall_restitutionOFast $* -Ofast`

OPT1=`gcc downloading_wall_restitution.c -o downloading_wall_restitutionOs $* -Os`

#OPT2=`gcc downloading_wall_restitution.c -o downloading_wall_restitution2 $* -O2`

#OPT3=`gcc downloading_wall_restitution.c -o downloading_wall_restitution3 $* -O3`

test -d "time" || mkdir -m 700 time

#***************************************************************
echo "[!] executando o programa 10 vezes com a flag 'Ofast' "
echo "[!] Aguarde ... "

echo " "
echo "PARAMETROS USADOS: $*" > time/gcctimeO0.txt
for i in $(seq 10);
 do
	perf stat -x , -e 'duration_time' ./downloading_wall_restitutionOFast 0.0001 8 10 0.09 >> time/gcctimeO0.txt
	echo "=============================================================" >> time/gcctimeO0.txt
  sleep 5
done
echo "[v] tempos salvos em 'time/gcctimeO0'"
#***********************************************************
echo "[!] executando o programa 10 vezes com a flag 'Os' "
echo "[!] Aguarde ... "

echo " "
echo "PARAMETROS USADOS: $*" > time/gcctimeO1.txt

for i in $(seq 10);
 do
  	perf stat -x , -e 'duration_time' ./downloading_wall_restitutionOs 0.0001 8 10 0.09 >> time/gcctimeO1.txt
	echo "=============================================================" >> time/gcctimeO1.txt
  sleep 5
done
echo "[v] tempos salvos em 'time/gcctimeO1'"

#***************************************************************
'''
echo "[!] executando o programa 10 vezes com flag 'O2' "
echo "[!] Aguarde ... "
echo " "

echo "PARAMETROS USADOS: $*" > time/gcctimeO2.txt
echo " "
for i in $(seq 10);
 do
   perf stat -x , -e 'duration_time' ./downloading_wall_restitution2 0.0001 8 10 0.09 >> time/gcctimeO2.txt
   echo "=============================================================" >> time/gcctimeO2.txt
   	sleep 5
done
echo "[v] tempos salvos em 'time/gcctimeO2'"

#***************************************************************
echo "[!] executando o programa 10 vezes com flag 'O3' "
echo "[!] Aguarde ... "

echo "PARAMETROS USADOS: $*" > time/gcctimeO3.txt
echo " "
for i in $(seq 10);
 do
   perf stat -x , -e 'duration_time' ./downloading_wall_restitution3 0.0001 8 10 0.09 >> time/gcctimeO3.txt
   echo "=============================================================" >> time/gcctimeO3.txt
 	sleep 5
 done

echo "[v] tempos salvos em 'time/gcctimeO3'"

'''