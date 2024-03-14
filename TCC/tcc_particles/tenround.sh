for i in $(seq 10);
 do
   perf stat -x , -e 'duration_time' ./downloading_wall_restitution 0.0001 8 10 0.09
 done


