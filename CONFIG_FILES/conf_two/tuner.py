# -*- coding: utf-8 -*-

# Copyright (c) 2022 Synodic Month, Juni May
# yaTuner is licensed under Mulan PSL v2.
# You can use this software according to the terms and conditions of the Mulan PSL v2.
# You may obtain a copy of Mulan PSL v2 at:
#          http://license.coscl.org.cn/MulanPSL2
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
# EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
# MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
# See the Mulan PSL v2 for more details.

from math import inf, log10
import os
from subprocess import TimeoutExpired
from typing import Callable, Dict, Mapping, Sequence, Tuple, Any

import seaborn as sns
import pandas as pd
import numpy as np

import yatuner
import logging
import GPyOpt

from yatuner import LinUCB
from rich.logging import RichHandler
from rich.table import Table
from rich.console import Console
from rich.progress import track
from scipy import stats
from scipy import optimize
from matplotlib import pyplot as plt


class Tuner:

    def __init__(self,
                 call_compile: Callable[[Sequence[str], Mapping, str], None],
                 call_running: Callable[[], float],
                 optimizers: Sequence[str],
                 parameters: Mapping[str, Tuple],
                 call_perf: Callable[[], Dict[str, Any]] = None,
                 workspace='yatuner.db_three',
                 log_level=logging.DEBUG,
                 norm_range=None,
                 deterministic=False,
                 u=None,
                 z=None,
                 std=None,
                 exec_data=None) -> None:
        """A tuner.

        Args:
            call_compile ((Sequence[str], Mapping, str) -> None): A function runing compilation process.
            call_running (() -> float): A function fetching result of target program.
            optimizers (Sequence[str]): List of on/of options.
            parameters (Mapping[str, Tuple]): List of parameters, in format of `param: (min, max, default)`
            call_perf (() -> Dict[str, Any], optional): Fectch feature of given program, used in linUCB. Defaults to None.
            workspace (str, optional): Directory to store output files. Defaults to 'yatuner.db'.
            log_level (optional): Log level. Defaults to logging.DEBUG.
            norm_range (float, optional): Cut the data of test run to get more accurate result, None if doing a symmetrization. Defaults to None.
            deterministic (bool): False if the result of `call_compile` is random to a certain extent, otherwise True. Defaults to False.
        """
        logging.basicConfig(format='[ %(name)s ] %(message)s',
                            handlers=[
                                RichHandler(level=log_level,
                                            markup=True,
                                            show_path=False)
                            ])
        self.logger = logging.getLogger('yatuner')
        self.logger.setLevel(log_level)

        self.exec_data=exec_data
        self.u=u
        self.z=z   
        self.std=std #desvio padrao

        self.call_compile = call_compile
        self.call_running = call_running
        self.call_perf = call_perf
        self.optimizers = optimizers
        self.parameters = parameters
        self.workspace = workspace
        if norm_range is None:
            self.symmetrization = True
            self.norm_range = None
        else:
            self.symmetrization = False
            self.norm_range = norm_range
        self.deterministic = deterministic

    def initialize(self):
        if not os.path.isdir(self.workspace):
            self.logger.info("workspace is not detected, creating one.")
            os.mkdir(self.workspace)

    def test_run(self, num_samples=200, warmup=50):
        """Doing a test run with no options indicated.

        Args:
            num_samples (int, optional): Times to run. Defaults to 200.
            warmup (int, optional): Times to warmup. Defaults to 50.
        """

        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if os.path.exists(self.workspace + '/test_run.csv'):
            self.logger.info("using existing test run result.")
            pd_data = pd.read_csv(self.workspace + '/test_run.csv')
            self.exec_data = pd_data.to_numpy().transpose()[0].tolist()

        else:
            self.call_compile(None, None, None)
            self.exec_data = []

            if self.deterministic:
                res = self.call_running()
                self.exec_data.append(res)
            else:
                for i in track(range(warmup), description='  warmup'):
                    res = self.call_running()
                    self.logger.debug(f"warmup {i}/{warmup} result: {res}")

                for _ in track(range(num_samples), description='test run'):
                    res = self.call_running()
                    self.exec_data.append(res)

            pd_data = pd.DataFrame({'test_run': self.exec_data})
            pd_data.to_csv(self.workspace + '/test_run.csv', index=0)

        if not self.deterministic:
            if self.symmetrization:
                self.exec_data.sort()
                kernel = stats.gaussian_kde(self.exec_data)
                maxima = optimize.minimize_scalar(
                    lambda x: -kernel(x),
                    bounds=(np.min(self.exec_data), np.max(self.exec_data)),
                    method='bounded').x[0]
                data_cnt = len(self.exec_data)
                for i in range(data_cnt):
                    self.exec_data.append(maxima + maxima - self.exec_data[i])
            else:
                self.exec_data = np.sort(
                    self.exec_data
                )[:int(len(self.exec_data) * self.norm_range)]
                '''
                calcular a média (ou a média aritmética) de um conjunto de dados. 
                A média é uma medida estatística que representa o valor central de um conjunto de números.   
                '''
            self.u = np.mean(self.exec_data)
            '''
            o desvio padrão de um conjunto de dados. O desvio padrão é uma medida de dispersão 
            que indica o quão dispersos estão os valores em relação à média.
            '''
            self.std = np.std(self.exec_data)
            #---truque para fazer a linha 156 compilar
            au = np.array([self.exec_data.argmax()])
            amostra_float = au.astype(float)
            #---
            kstest = stats.kstest(amostra_float, 'norm', (self.u, self.std))

            self.logger.info(
                f"test run finished with u: {self.u:.2f}, std: {self.std:.2f}")

            if kstest.pvalue < 0.05:
                self.logger.info("test run ... [bold green]OK[/]")
            else:
                self.logger.warning(
                    "[red]execution time disobey normal distribution[/]")

            plt.clf()
            plt.hist(self.exec_data, bins=50, density=True, label='test run')
            plt.legend()
            plt.savefig(self.workspace + "/test_run_distribution.png")
        else:
            self.u = self.exec_data[0]
            self.logger.info(
                f"teste run finished with res: {self.exec_data[0]:.2f}")

    def hypotest_optimizers(self,
                            num_samples=10,
                            z_threshold=0.05,
                            t_threshold=0.05,
                            num_epochs=30
                            ):
        """Hypothesis test for on/of options.

        Args:
            num_samples (int, optional): Sampling times for each option. Defaults to 10.
            z_threshold (float, optional): Z threshold. Defaults to 0.05.
            t_threshold (float, optional): T threshold. Defaults to 0.05.
            num_epochs (int, optional): optimization epoches after selection. Defaults to 30.
        """
        print("\ndef hypotest_optimizers | o otimizador escolhido é quando: p < z_threshold=0.05 or t < t_threshold=0.05")
        print("p: Função de sobrevivência de uma distribuição normal (gaussiana).")
        print("t: indica se a média de uma única amostra é significativamente diferente de uma média populacional conhecida ou hipotética.\n")
        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if os.path.exists(self.workspace + '/selected_optimizers.txt'):
            self.logger.info("using existing selected optimizers.")
            return

        self.selected_optimizers = []
        hypotest_exec_data = []
        for i, optimizer in enumerate(self.optimizers):
            try:
                self.call_compile([optimizer], None, None)
            except RuntimeError as err:
                self.logger.error(f"[red]compile error with {optimizer}[/]")
                self.logger.exception(err)
                continue

            samples = np.zeros(num_samples)
            err = False
            if self.deterministic:
                try:
                    res = self.call_running()
                except RuntimeError as err:
                    self.logger.error(f"[red]runtime error for {optimizer}[/]")
                    self.logger.exception(err)
                    err = True
                    break

                samples[0] = res
                hypotest_exec_data.append(res)

            else:
                for j in range(num_samples):
                    try:
                    	res = self.call_running()
                    except RuntimeError as err:
                        self.logger.error(
                            f"[red]runtime error for {optimizer}[/]")
                        self.logger.exception(err)
                        err = True
                        break

                    samples[j] = res
                    hypotest_exec_data.append(res)

            if err: 
                continue

            if not self.deterministic:
                samples_mean = samples.mean()
                z = (samples_mean - self.u) / (self.std / np.sqrt(len(samples)))
                '''
                A função stats.norm.sf é parte do módulo scipy.stats em Python e é utilizada 
                para calcular a função de sobrevivência (também conhecida como complemento da 
                função de distribuição acumulativa) de uma distribuição normal (gaussiana).

                A função de sobrevivência em um ponto x de uma distribuição é a probabilidade de que a 
                variável aleatória seja maior que x. Matematicamente, para uma distribuição normal, 
                isso pode ser expresso como:
                '''
                p = 2 * stats.norm.sf(abs(z))

                '''
                Ela é usada para realizar um teste t de uma amostra. 
                Um teste t de uma amostra é um teste estatístico usado para determinar 
                se a média de uma única amostra é significativamente diferente de uma 
                média populacional conhecida ou hipotética.
                '''
                t = stats.ttest_1samp(samples, self.u).pvalue 

                self.logger.debug(f"{i}/{len(self.optimizers)} {optimizer} "
                                  f"u: {self.u:.2f} -> {samples_mean:.2f}, "
                                  f"p: {p:.2f}, t: {t:.2f}")

                if (p < z_threshold or t < t_threshold) and z < 0:
                    self.selected_optimizers.append(optimizer)
                    self.logger.info(f"[green]{optimizer} is selected[/]")
            else:
                self.logger.debug(f"{i}/{len(self.optimizers)} {optimizer} "
                                  f"res: {samples[0]:.2f}")
                if samples[0] < self.u:
                    self.selected_optimizers.append(optimizer)
                    self.logger.info(f"[green]{optimizer} is selected[/]")

        self.logger.info(
            f"{len(self.selected_optimizers)} optimizers selected")
        self.logger.info(f"optimizing optimizers")
        cnt = 0

        def step(vals) -> float:

            nonlocal cnt

            step_optimizers = []
            for i, opt in enumerate(self.selected_optimizers):
                if (vals[0][i]):
                    step_optimizers.append(opt)
            try:
                self.call_compile(step_optimizers, None, None)
            except RuntimeError as err:
                self.logger.error("[red]compile error[/]")
                self.logger.exception(err)
                return inf

            res = 0
            if self.deterministic:
                res = self.call_running()
            else:
                for _ in track(range(num_samples), f'epoch {cnt:<5}'):
                    res += self.call_running()
                res /= num_samples

            self.logger.debug(f'{cnt}/{num_epochs} result: {res:.2f}')

            cnt += 1

            return res

        bounds = [{
            'name': opt,
            'type': 'discrete',
            'domain': (0, 1)
        } for opt in self.selected_optimizers]

        method = GPyOpt.methods.BayesianOptimization(step,
                                                     domain=bounds,
                                                     acquisition_type='EI',
                                                     acquisition_weight=2,
                                                     initial_design_numdata=10)
        method.run_optimization(max_iter=num_epochs, eps=0.0)
        # method.plot_convergence(self.workspace + '/convergence.png')

        self.logger.info(f"best result: {method.fx_opt.flatten()[0]}")
        self.logger.info(f"best option: {method.x_opt}")
        new_optimizers = []
        for idx, x in enumerate(method.x_opt):
            if x:
                new_optimizers.append(self.selected_optimizers[idx])
        self.selected_optimizers = new_optimizers

        with open(self.workspace + '/selected_optimizers.txt',
                  'w',
                  encoding='utf-8') as f:
            f.writelines(
                [optimizer + '\n' for optimizer in self.selected_optimizers])

        if not self.deterministic:
            plt.clf()
            # bin_min = min(np.min(self.exec_data), np.min(hypotest_exec_data))
            # bin_max = max(np.max(self.exec_data), np.max(hypotest_exec_data))
            # bins = np.arange(bin_min, bin_max + 500, 500)
            plt.hist(
                self.exec_data,
                #bins=bins,
                density=True,
                alpha=0.5,
                label='test run')
            plt.hist(
                hypotest_exec_data,
                #  bins=bins,
                density=True,
                alpha=0.5,
                label='hypotest-optimizers')
            plt.legend()
            plt.savefig(self.workspace +
                        "/hypotest_optimizers_distribution.png")

    def hypotest_parameters(self, num_samples=10, t_threshold=0.05):
        """Hypothesis test for parameters.

        Args:
            num_samples (int, optional): Sampling times for each parameter. Defaults to 10.
            t_threshold (float, optional): T threshold. Defaults to 0.05.
        """
        print("\ndef hypotest_parameters: ('p') Este teste é usado para verificar se há uma diferença significativa entre as médias das duas amostras.")
        print("O otimizador é escolhido quando: p < t_threshold=0.05\n")
        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if os.path.exists(self.workspace + '/selected_parameters.txt'):
            self.logger.info("using existing selected parameters")
            return

        if os.path.exists(self.workspace + '/selected_optimizers.txt'):
            self.selected_optimizers = []
            with open(self.workspace + '/selected_optimizers.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_optimizers = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_optimizers)} optimizers")
        else:
            self.logger.error(
                "no selected optimizers detected, "
                "please run hypothesis test for optimizers first.")
            return

        self.selected_parameters = []
        i = 0
        for parameter, r in self.parameters.items():
            r_min, r_max, default = r

            samples_min = np.zeros(num_samples)
            try:
                self.call_compile(self.selected_optimizers, {parameter: r_min},
                                  None)
            except RuntimeError as err:
                self.logger.error(f"[red]compile error with {parameter}[/]")
                self.logger.exception(err)
                continue
            except TimeoutExpired:
                self.logger.warning(
                    f"[red]compile timeout with {parameter}[/]")

            if self.deterministic:
                res = self.call_running()
                samples_min[0] = res
            else:
                for j in range(num_samples):
                    res = self.call_running()
                    samples_min[j] = res

            samples_max = np.zeros(num_samples)
            try:
                self.call_compile(self.selected_optimizers, {parameter: r_max},
                                  None)
            except RuntimeError as err:
                self.logger.error(f"[red]compile error with {parameter}[/]")
                self.logger.exception(err)
                continue
            except TimeoutExpired:
                self.logger.warning(
                    f"[red]compile timeout with {parameter}[/]")

            if self.deterministic:
                res = self.call_running()
                samples_max[0] = res
            else:
                for j in range(num_samples):
                    res = self.call_running()
                    samples_max[j] = res

            if not self.deterministic:
                '''
                O teste de levene é um teste para determinar se as variâncias 
                ou os desvios padrão de dois ou mais grupos são diferentes. 

                Variância é uma medida de dispersão que mostra o quão distante 
                cada valor desse conjunto está do valor central (médio).
                '''
                l = stats.levene(samples_min, samples_max).pvalue

                '''
                Este teste é usado para verificar se há uma diferença significativa entre as médias das duas amostras.

                equal_var indica se as variâncias das amostras são assumidas como iguais

                '''
                p = stats.ttest_ind(samples_min,
                                    samples_max,
                                    equal_var=(l > 0.05)).pvalue

                self.logger.debug(f"{i}/{len(self.parameters)} {parameter} "
                                  f"min: {np.mean(samples_min):.2f}, "
                                  f"max: {np.mean(samples_max):.2f} "
                                  f"l: {l:.2f}, p: {p:.2f}")

                # if p < t_threshold=0.05
                if p < t_threshold:
                    self.selected_parameters.append(parameter)
                    self.logger.info(f"[green]{parameter} is selected[/]")
            else:
                self.logger.debug(f"{i}/{len(self.parameters)} {parameter} "
                                  f"min: {samples_min[0]:.2f} "
                                  f"max: {samples_max[0]:.2f}")
                if samples_min[0] != samples_max[0]:
                    self.selected_parameters.append(parameter)
                    self.logger.info(f"[green]{parameter} is selected[/]")

            i += 1

        with open(self.workspace + '/selected_parameters.txt',
                  'w',
                  encoding='utf-8') as f:
            f.writelines(
                [parameter + '\n' for parameter in self.selected_parameters])

    def optimize_linUCB(self,
                        alpha=0.5,
                        num_epochs=200,
                        num_samples=10,
                        num_bins=25,
                        nth_choice=3,
                        metric='duration_time') -> None:
        """Optimize selected parameters with linUCB."""
        print("\ndef optimize_linUCB: \n")
        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if self.deterministic:
            self.logger.info("using bayesian instead of LinUCB.")
            self.optimize(num_samples=num_samples, num_epochs=num_samples)
            return

        if self.call_perf == None:
            self.logger.error("call_perf not found")
            return

        if num_bins > num_epochs:
            self.logger.error("num_epochs needs to be larger than num_bins")
            return

        if os.path.exists(self.workspace + '/optimized_parameters.txt'):
            self.logger.info("using existing optimized parameters.")
            return

        if os.path.exists(self.workspace + '/selected_optimizers.txt'):
            self.selected_optimizers = []
            with open(self.workspace + '/selected_optimizers.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_optimizers = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_optimizers)} optimizers")

        try:
            self.call_compile(self.selected_optimizers, None, None)
        except RuntimeError as err:
            self.logger.error(f"[red]compile error[/]")
            self.logger.exception(err)
            return

        baseline = 0
        for _ in track(range(num_samples), "before optimization"):
            baseline += self.call_running()
        baseline /= num_samples
        self.logger.info(f"execution result before optimize: {baseline}")

        if os.path.exists(self.workspace + '/selected_parameters.txt'):
            self.selected_parameters = []
            with open(self.workspace + '/selected_parameters.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_parameters = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_parameters)} parameters")

        print("\nAplicando LinUCB\n")

        dim = len(self.call_perf())
        features = [log10(1 + x) for x in self.call_perf().values()]
        ucbs = [
            LinUCB.LinUCB(dim,
                          np.linspace(self.parameters[param][0],
                                      self.parameters[param][1],
                                      num_bins,
                                      endpoint=True,
                                      dtype='int'),
                          alpha=alpha,
                          nth_choice=nth_choice)
            for param in self.selected_parameters
        ]
        for ucb in ucbs:
            ucb.init()
        choices = np.zeros(len(self.selected_parameters), dtype=int)
        timearr = np.zeros(num_epochs)
        best_choices = np.zeros(len(self.selected_parameters), dtype=int)
        best_time = float('inf')
        for i in track(range(num_epochs), description='optimizing'):
            for idx, ucb in enumerate(ucbs):
                choices[idx] = ucb.recommend(features)
            step_parameters = {}
            for idx, param in enumerate(self.selected_parameters):
                step_parameters[param] = choices[idx]
            # print(step_parameters)
            try:
                self.call_compile(self.selected_optimizers, step_parameters,
                                  None)
                new_perf = self.call_perf()
                features = [log10(1 + x) for x in new_perf.values()]
                new_time = new_perf[metric] / 1000
                reward = (baseline - new_time) / 100
                self.logger.debug(
                    f"r={reward:.2f} t={new_time:.2f} baseline={baseline:.2f}")
                for ucb in ucbs:
                    ucb.update(reward)
                timearr[i] = new_time
                if new_time < best_time:
                    best_time = new_time
                    best_choices = choices.copy()
            except TimeoutExpired:
                self.logger.warning(f"[red]compile timeout[/]")
        plt.clf()
        plt.plot(timearr)
        plt.savefig(self.workspace + '/convergence_linUCB.png')
        print(choices)
        with open(self.workspace + '/optimized_parameters.txt',
                  'w',
                  encoding='utf-8') as file:
            for idx, param in enumerate(self.selected_parameters):
                file.write(param + " " + str(best_choices[idx]) + '\n')

    def optimize(self, num_samples=10, num_epochs=60) -> None:
        """Optimize selected parameters with bayesian."""
        print("\ndef optimize\n")
        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if os.path.exists(self.workspace + '/optimized_parameters.txt'):
            self.logger.info("using existing optimized parameters.")
            return

        if os.path.exists(self.workspace + '/selected_optimizers.txt'):
            self.selected_optimizers = []
            with open(self.workspace + '/selected_optimizers.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_optimizers = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_optimizers)} optimizers")

        try:
            self.call_compile(self.selected_optimizers, None, None)
        except RuntimeError as err:
            self.logger.error(f"[red]compile error[/]")
            self.logger.exception(err)
            return
        except TimeoutExpired:
            self.logger.error(f"[red]compile timeout[/]")

        res = 0
        if self.deterministic:
            res = self.call_running()
        else:
            for _ in track(range(num_samples), "before optimization"):
                res += self.call_running()
            res /= num_samples
        self.logger.info(f"execution result before optimize: {res}")

        if os.path.exists(self.workspace + '/selected_parameters.txt'):
            self.selected_parameters = []
            with open(self.workspace + '/selected_parameters.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_parameters = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_parameters)} parameters")

        cnt = 0

        def step(vals) -> float:

            nonlocal cnt

            step_parameters = {}
            for i, parameter in enumerate(self.selected_parameters):
                v = round(vals[0][i] * (self.parameters[parameter][1] -
                                        self.parameters[parameter][0]) +
                          self.parameters[parameter][0])
                step_parameters[parameter] = int(v)

            try:
                self.call_compile(self.selected_optimizers, step_parameters,
                                  None)
            except RuntimeError as err:
                self.logger.error("[red]compile error[/]")
                self.logger.exception(err)
                return inf

            res = 0
            if self.deterministic:
                res = self.call_running()
            else:
                for _ in track(range(num_samples), f'epoch {cnt:<5}'):
                    res += self.call_running()
                res /= num_samples

            self.logger.debug(f'{cnt}/{num_epochs} result: {res:.2f}')

            cnt += 1

            return res

        bounds = [{
            'name': parameter,
            'type': 'continuous',
            'domain': (0, 1)
        } for parameter in self.selected_parameters]

        '''
        A convergência da otimização bayesiana refere-se ao processo pelo qual o algoritmo melhora sua estimativa do valor ótimo da função objetivo ao longo das iterações.

        O gráfico que mostra a convergência da otimização ao longo das iterações. Isso pode ser útil para visualizar como a otimização melhora a estimativa do valor ótimo ao longo do tempo.
        O gráfico de convergência mostrará como o valor da função objetivo estimado melhora com o número de iterações.

        '''

        method = GPyOpt.methods.BayesianOptimization(step,
                                                     domain=bounds,
                                                     acquisition_type='LCB',
                                                     acquisition_weight=0.2)
        method.run_optimization(max_iter=num_epochs)
        method.plot_convergence(self.workspace + '/convergence.png')

        self.logger.info(f"best result: {method.fx_opt.flatten()[0]}")
        self.logger.info(f"best option: {method.x_opt}")

        vals = method.x_opt
        with open(self.workspace + '/optimized_parameters.txt',
                  'w',
                  encoding='utf-8') as file:

            for i, parameter in enumerate(self.selected_parameters):
                v = round(vals[i] * (self.parameters[parameter][1] -
                                     self.parameters[parameter][0]) +
                          self.parameters[parameter][0])

                file.write(f'{parameter} {int(v)}\n')

    def run(self, num_samples=10):
        """Run result and existing options to compare.

        Args:
            num_samples (int, optional): Sampling times. Defaults to 10.
        """
        if self.deterministic and num_samples != 1:
            self.logger.warning(f"num_samples of {num_samples} "
                                f"is not necessary for deterministic goal.")
            num_samples = 1

        if os.path.exists(self.workspace + '/result.csv'):
            self.logger.info("aready done.")
            return

        samples_ofast = np.zeros(num_samples)
        samples_os = np.zeros(num_samples)
        samples_o0 = np.zeros(num_samples)
        samples_o1 = np.zeros(num_samples)
        samples_o2 = np.zeros(num_samples)
        samples_o3 = np.zeros(num_samples)
        samples_optimizers = np.zeros(num_samples)
        samples_parameters = np.zeros(num_samples)

        self.call_compile(None, None, '-Ofast')
        if self.deterministic:
            res = self.call_running()
            samples_ofast[0] = res
        else:
            for i in track(range(num_samples), description='-Ofast'):
                res = self.call_running()
                samples_ofast[i] = res

        self.call_compile(None, None, '-Os')
        if self.deterministic:
            res = self.call_running()
            samples_os[0] = res
        else:
            for i in track(range(num_samples), description='   -Os'):
                res = self.call_running()
                samples_os[i] = res

        self.call_compile(None, None, '-O0')
        if self.deterministic:
            res = self.call_running()
            samples_o0[0] = res
        else:
            for i in track(range(num_samples), description='   -O0'):
                res = self.call_running()
                samples_o0[i] = res

        self.call_compile(None, None, '-O1')
        if self.deterministic:
            res = self.call_running()
            samples_o1[0] = res
        else:
            for i in track(range(num_samples), description='   -O1'):
                res = self.call_running()
                samples_o1[i] = res

        self.call_compile(None, None, '-O2')
        if self.deterministic:
            res = self.call_running()
            samples_o2[0] = res
        else:
            for i in track(range(num_samples), description='   -O2'):
                res = self.call_running()
                samples_o2[i] = res

        self.call_compile(None, None, '-O3')
        if self.deterministic:
            res = self.call_running()
            samples_o3[0] = res
        else:
            for i in track(range(num_samples), description='   -O3'):
                res = self.call_running()
                samples_o3[i] = res

        if os.path.exists(self.workspace + '/selected_optimizers.txt'):
            self.selected_optimizers = []
            with open(self.workspace + '/selected_optimizers.txt',
                      'r',
                      encoding='utf-8') as file:
                self.selected_optimizers = [
                    x.strip() for x in file.readlines()
                ]

            self.logger.info(
                f"loaded {len(self.selected_optimizers)} optimizers")
        else:
            pass

        self.call_compile(self.selected_optimizers, None, None)
        if self.deterministic:
            res = self.call_running()
            samples_optimizers[0] = res
        else:
            for i in track(range(num_samples), description='optimizers'):
                res = self.call_running()
                samples_optimizers[i] = res

        if os.path.exists(self.workspace + '/optimized_parameters.txt'):
            self.optimized_parameters = {}
            with open(self.workspace + '/optimized_parameters.txt',
                      'r',
                      encoding='utf-8') as file:
                for line in file.readlines():
                    parameter, val = line.strip().split(' ')
                    self.optimized_parameters[parameter] = val

            self.logger.info(
                f"loaded {len(self.optimized_parameters)} parameters")
        else:
            pass

        self.call_compile(self.selected_optimizers, self.optimized_parameters,
                          None)
        if self.deterministic:
            res = self.call_running()
            samples_parameters[0] = res
        else:
            for i in track(range(num_samples), description='parameters'):
                res = self.call_running()
                samples_parameters[i] = res

        if not self.deterministic:
            bin_min = min(np.min(samples_os), np.min(samples_o0),
                          np.min(samples_o1), np.min(samples_o2),
                          np.min(samples_o3), np.min(samples_ofast),
                          np.min(samples_optimizers),
                          np.min(samples_parameters))
            bin_max = max(np.max(samples_os), np.max(samples_o0),
                          np.max(samples_o1), np.max(samples_o2),
                          np.max(samples_o3), np.max(samples_ofast),
                          np.max(samples_optimizers),
                          np.max(samples_parameters))
            bins = np.arange(bin_min, bin_max + 500, 500)

            plt.clf()
            plt.hist(samples_o0,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='O0')
            plt.hist(samples_o1,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='O1')
            plt.hist(samples_o2,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='O2')
            plt.hist(samples_o3,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='O3')
            plt.hist(samples_os,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='Os')
            plt.hist(samples_ofast,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='Ofast')
            plt.hist(samples_optimizers,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='Optimizers')
            plt.hist(samples_parameters,
                     bins=bins,
                     density=True,
                     alpha=0.5,
                     label='Parameters')
            plt.xlabel("Time/Tick")
            plt.ylabel("Density")
            plt.legend()
            plt.savefig(self.workspace + "/result.png")

        pd_data = pd.DataFrame({
            'O0': samples_o0,
            'O1': samples_o1,
            'O2': samples_o2,
            'O3': samples_o3,
            'Os': samples_os,
            'Ofast': samples_ofast,
            'Optimizers': samples_optimizers,
            'Parameters': samples_parameters
        })
        pd_data.to_csv(self.workspace + "/result.csv", index=0)

        mean_ofast = samples_ofast.mean()
        mean_os = samples_os.mean()
        mean_o0 = samples_o0.mean()
        mean_o1 = samples_o1.mean()
        mean_o2 = samples_o2.mean()
        mean_o3 = samples_o3.mean()
        mean_optimizers = samples_optimizers.mean()
        mean_parameters = samples_parameters.mean()

        minimal = min(mean_ofast, mean_os, mean_o0, mean_o1, mean_o2, mean_o3,
                      mean_optimizers, mean_parameters)

        score_ofast = 100 * minimal / mean_ofast
        print("\n\n")
        print("samples_os: ", samples_os)
        print("minimal: ", minimal)
        print("mean_os: ", mean_os)
        print("\n\n")
        score_os = 100 * minimal / mean_os
        score_o0 = 100 * minimal / mean_o0
        score_o1 = 100 * minimal / mean_o1
        score_o2 = 100 * minimal / mean_o2
        score_o3 = 100 * minimal / mean_o3
        score_optimizers = 100 * minimal / mean_optimizers
        score_parameters = 100 * minimal / mean_parameters

        delta_optimizers = (score_optimizers - score_o2) / score_o2 * 100
        delta_parameters = (score_parameters - score_o2) / score_o2 * 100

        table = Table(title="Result")
        table.add_column("Method")
        table.add_column(f"Result", style="cyan")
        table.add_column("Score", style="green")
        table.add_column("Delta", style="green")
        table.add_row("O0", f"{mean_o0:.2f}", f"{score_o0:.2f}", "")
        table.add_row("O1", f"{mean_o1:.2f}", f"{score_o1:.2f}", "")
        table.add_row("O2", f"{mean_o2:.2f}", f"{score_o2:.2f}", "")
        table.add_row("O3", f"{mean_o3:.2f}", f"{score_o3:.2f}", "")
        table.add_row("Ofast", f"{mean_ofast:.2f}", f"{score_ofast:.2f}", "")
        table.add_row("Os", f"{mean_os:.2f}", f"{score_os:.2f}", "")
        table.add_row("Optimizers", f"{mean_optimizers:.2f}",
                      f"{score_optimizers:.2f}", f"{delta_optimizers:.2f}%")
        table.add_row("Parameters", f"{mean_parameters:.2f}",
                      f"{score_parameters:.2f}", f"{delta_parameters:.2f}%")

        console = Console()
        console.print(table)

    def plot_data(self) -> None:
        """Plot result in violin graph."""
        if self.deterministic:
            self.logger.warning("skipping plotting data.")
            return

        if not os.path.exists(self.workspace + '/result.csv'):
            self.logger.error("No data found.")
            return

        plt.style.use('seaborn')
        pd_data = pd.read_csv(self.workspace + "/result.csv")
        plt.clf()
        plt.figure(figsize=(10, 5), dpi=100)
        sns.violinplot(
            data=pd_data,
            orient='horizontal',
            palette='Set3',
            width=0.9,
        )
        plt.title('Time Comparison')
        plt.ylabel('Optimization_methods')
        plt.xlabel('Time/Tick - Lower the better')
        plt.grid()
        plt.savefig(self.workspace + "/result_violin.png")
