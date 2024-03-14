# -*- coding: utf-8 -*-

import os
import yatuner
import logging

cc = 'gcc-10'
src = '/home/thiago/TCC/tcc_particles/downloading_wall_restitution.c '
parametros_exec = ' 0.0001 8 10 0.09'
out = 'build/timePerf/downloading_wall_restitution'
binario = 'downloading_wall_restitution'
base = '-O3'
metric = 'duration_time'

if not os.path.isdir('./buildPerf/time/'):
    try:
        os.makedirs('./buildPerf/time/')
    except OSError as error:
        print(error)

optimizers = set(yatuner.utils.fetch_gcc_optimizers(cc=cc)).difference(
    yatuner.utils.fetch_gcc_enabled_optimizers(options=base))
optimizers.remove('-fipa-pta')

parameters = yatuner.utils.fetch_gcc_parameters(cc=cc)


def comp(optimizers, parameters, additional):
    options = '-std=c11 -lgsl -lgslcblas -lm -pg '
    if additional is not None:
        options += f'{additional} '
    else:
        options += f'{base} '

    if optimizers is not None:
        for optimizer in optimizers:
            options += f'{optimizer} '

    if parameters is not None:
        for parameter, val in parameters.items():
            options += f'--param={parameter}={val} '

    res = yatuner.utils.execute(
        f'{cc} {src} {options} -o {out}')
    
    if res['returncode'] != 0:
        raise RuntimeError(res['stderr'])


def run():
    return yatuner.utils.fetch_perf_stat(binario + parametros_exec)[metric] / 1000


def perf():
    return yatuner.utils.fetch_perf_stat(binario + parametros_exec)


tuner = yatuner.Tuner(comp,
                    run,
                    optimizers,
                    parameters,
                    call_perf=perf,
                    norm_range=0.99,
                    deterministic=False)

# Tuning process!
tuner.initialize()
tuner.test_run(num_samples=50, warmup=5)#200, 10
tuner.hypotest_optimizers(num_samples=5, num_epochs=30)#1, 30 ()
tuner.hypotest_parameters(num_samples=5)
tuner.optimize(num_samples=20, num_epochs=100) #Using bayesian optimization 
tuner.optimize_linUCB(alpha=0.25,
                      num_bins=30,
                      num_epochs=20, 
                      nth_choice=4,
                      metric=metric) #using linUCB
tuner.run(num_samples=5) #50 -Ofast, -Os, -O0, (10 - 30)
tuner.plot_data()