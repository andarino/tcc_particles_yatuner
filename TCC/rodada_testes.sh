for i in {1..10}; do
	perf stat -x , -e 'duration_time' ./tcc_particles/downloading_wall_restitution_GPU 0.0001 8 10 0.09
done
