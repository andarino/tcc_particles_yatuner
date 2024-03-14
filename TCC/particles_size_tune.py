# -*- coding: utf-8 -*-

import os
import yatuner
import logging

cc = 'gcc-10'
src = './tcc_particles/downloading_wall_restitution.c'
parametros_exec = '0.0001 8 10 0.09'
out = './build/sizeTune/downloading_wall_restitution'
base = '-Os'

if not os.path.isdir('./build/sizeTune/'):
    try:
        os.makedirs('./build/sizeTune/')
    except OSError as error:
        print(error)

optimizers = set(yatuner.utils.fetch_gcc_optimizers(cc=cc)).difference(
        yatuner.utils.fetch_gcc_enabled_optimizers(options=base))
optimizers.remove('-fipa-pta')
parameters = yatuner.utils.fetch_gcc_parameters(cc=cc)

#./autotestes.sh -std=c11 -lgsl -lgslcblas -lm -pg"

logger = logging.getLogger('Particles')
logger.setLevel(logging.INFO)

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
    return yatuner.utils.fetch_size(f'{out}')


def perf():
    return NotImplementedError()

tuner = yatuner.Tuner(call_compile = comp,
                      call_running =run,
                      optimizers = optimizers,
                      parameters = parameters,
                      call_perf=perf,
                      norm_range=0.99,
                      workspace='yatuner_size.db',
                      log_level=logging.DEBUG,
                      deterministic=True)

logger.info(f'Performing Optimization /build/sizeTune')

tuner.initialize()
#tuner.test_run(num_samples=1, warmup=0)
tuner.hypotest_optimizers(num_samples=1, num_epochs=30)
tuner.hypotest_parameters(num_samples=1)
tuner.optimize(num_samples=1, num_epochs=60) # using bayesian optimization
tuner.run(num_samples=1)