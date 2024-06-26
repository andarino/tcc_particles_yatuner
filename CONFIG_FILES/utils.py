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

import subprocess
import re
import os
import numpy as np
import platform
import ast
from typing import List, Tuple
from typing import Any, Dict


def execute(command) -> Dict[str, Any]:
    """Execute given command.

    Args:
        command: Command to be executed.

    """

    p = subprocess.Popen(command,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(timeout=600)
    returncode = p.returncode

    p.terminate()

    return {
        'returncode': returncode,
        'stdout': stdout.decode(),
        'stderr': stderr.decode(),
    }


def fetch_perf_stat(command) -> Dict[str, Any]:
    perf_command = (
        'perf stat -x,'
        ' -e branch-instructions '
        ' -e branch-misses '
        ' -e bus-cycles '
        ' -e cache-misses '
        ' -e cache-references '
        ' -e duration_time'
        ' -e cpu-cycles '
        ' -e instructions '
        ' -e ref-cycles '
        ' -e alignment-faults '
        ' -e bpf-output '
        ' -e context-switches '
        ' -e cpu-clock '
        ' -e cpu-migrations '
        ' -e dummy '
        ' -e emulation-faults '
        ' -e major-faults '
        ' -e minor-faults '
        ' -e page-faults '
        ' -e task-clock '
        ' -e duration_time '
        + command)

    res = execute(perf_command)
    
    if res['returncode'] != 0:
    	raise RuntimeError(res['stderr'])
    else:
        counter_tuples = []
        for line in res['stderr'].splitlines():
            stat = line.split(',')
            counter_tuples.append((stat[2], stat[0]))
        # print(counter_tuples)
        counter_dict = dict(
            (x, 0 if y in ['<not supported>', '<not counted>'] else float(y))
            for (x, y) in counter_tuples)
        return counter_dict


def ir2vec(llvm_ir, outfile, mode='fa', vocab='vocab.txt', level='p'):
    """Call ir2vec tool and fetch the result."""
    if (os.path.exists(outfile)):
        os.remove(outfile)
    command = (f'ir2vec -{mode} -vocab {vocab}'
               f' -o {outfile} -level {level} {llvm_ir}')
    res = execute(command)
    if (res['returncode'] != 0):
        raise RuntimeError(res)
    vec = np.genfromtxt(outfile)
    return vec


def generate_llvm(src, llvm_ir, llvm_compiler='clang') -> None:
    """Generate llvm ir for source code."""
    command = f'{llvm_compiler} -S -emit-llvm -o {llvm_ir} {src}'
    res = execute(command)
    if (res['returncode'] != 0):
        raise RuntimeError(res)


def fetch_src_feature(src: str,
                      llvm_ir=None,
                      outfile=None,
                      clean=False) -> np.array:
    """Fetch source code feature by using ir2vec."""
    if llvm_ir is None:
        llvm_ir = src + '.ll'
    if outfile is None:
        outfile = src + '.csv'

    llvm_compiler = ''

    if src.endswith(('.c')):
        llvm_compiler = 'clang'
    elif src.endswith(('.F')):
        llvm_compiler = 'flang'
    elif src.endswith(('.cpp', '.cc', '.cxx')):
        llvm_compiler = 'clang++'

    try:
        generate_llvm(src, llvm_ir, llvm_compiler)
    except RuntimeError as err:
        print('[ ERROR ] llvm-ir: ', err)

    try:
        vec = ir2vec(llvm_ir, outfile)
    except RuntimeError as err:
        print('[ ERROR ] ir2vec: ', err)

    if (clean):
        if (os.path.exists(llvm_ir)):
            os.remove(llvm_ir)

        if (os.path.exists(outfile)):
            os.remove(outfile)

    return vec


def fetch_arch() -> str:
    """Fetch machine architecture information.

    Returns:
        str: The architecture.
    """
    return platform.machine()


def fetch_gcc_optimizers(cc='gcc') -> List[str]:
    """Fetch gcc optimizers

    Args:
        cc (str, optional): Specified gcc binary. Defaults to 'gcc'.

    Returns:
        List[str]: Optimizers.
    """
    optimizers = []

    opt = execute(f'{cc} --help=optimizers')['stdout']

    optimizers = re.findall(r'^  (-f[a-z0-9-]+) ', opt, re.M)

    return optimizers


def fetch_gcc_version(cc='gcc') -> Tuple[int, int, int]:
    """Fetch gcc version.

    Args:
        cc (str, optional): Specified gcc binary. Defaults to 'gcc'.

    Returns:
        Tuple[int, int, int]: gcc version in triplet.
    """
    version = None
    ver = re.search(r'([0-9]+)[.]([0-9]+)[.]([0-9]+)',
                    execute(f'{cc} --version')['stdout'])
    if ver:
        version = tuple(map(int, ver.group(1, 2, 3)))

    return version


def fetch_gcc_parameters(cc='gcc',
                         params_def: str = None
                         ) -> Dict[str, Tuple[int, int, int]]:
    """Fetch gcc parameters and its range

    Args:
        cc (str, optional): Specified gcc binary. Defaults to 'gcc'.
        params_def (str, optional): `params.def` of gcc source for version earlier then 10. Defaults to None.

    Raises:
        RuntimeError: If `params.def` is not given for gcc earlier than 10.

    Returns:
        Dict[str, Tuple[int, int, int]]: Parameters in format of `param: (min, max, default)`.
    """

    version = fetch_gcc_version(cc)

    params = {}

    if version[0] > 9:
        params_str = execute(f'{cc} -Q --help=params ')['stdout']
        regex = r'^  --param=([a-z-]+)=(<[0-9]+,[0-9]+>)?\s+([0-9]+)'

        for param, r, default in re.findall(regex, params_str, re.M):
            if r == '':
                params[param] = (0, int(default) * 10, int(default))
            else:
                r_min = int(r[1:-1].split(',')[0])
                r_max = int(r[1:-1].split(',')[1])
                params[param] = (r_min, r_max, int(default))

    else:
        if params_def is None:
            raise RuntimeError(
                "Params definition file must be given for gcc-9 or earlier")
        params_str = execute(f'{cc} --help=params ')['stdout']

        regex = r'^  --param=([a-z-]+)'
        raw_params = re.findall(regex, params_str, re.M)

        regex = r'^  ([a-z-]+)'
        params_def_file = open(params_def)
        raw_params_def = params_def_file.read()
        for m in re.finditer(r'DEFPARAM *\((([^")]|"[^"]*")*)\)',
                             raw_params_def):
            param_def_str = (
                m.group(1)  #
                .replace('GGC_MIN_EXPAND_DEFAULT', '30')  #
                .replace('GGC_MIN_HEAPSIZE_DEFAULT', '4096')  #
                .replace('50 * 1024 * 1024', '52428800')  #
                .replace('128 * 1024 * 1024', '134217728')  #
                .replace('INT_MAX', '2147483647'))  #

            param, _, default, r_min, r_max = ast.literal_eval(
                '[' + param_def_str.split(',', 1)[1] + ']')

            if param not in raw_params:
                continue

            if r_max == 0:
                r_max = min(default * 10, 2147483647)

            params[param] = (r_min, r_max, default)

        params_def_file.close()

    return params


def fetch_gcc_enabled_optimizers(cc='gcc', options='-O3') -> List[str]:
    """Fetch enabled optimizers of certain options

    Args:
        cc (str, optional): Specified gcc binary. Defaults to 'gcc'.
        options (str, optional): Options to look into. Defaults to '-O3'.

    Returns:
        List[str]: Enabled optimizers.
    """
    raw = execute(f'{cc} {options} -Q --help=optimizers')['stdout']
    regex = r'(-f[a-z0-9-]+)\s+(\[enabled\]|\[disabled\])'
    optimizers = []
    for optimizer, status in re.findall(regex, raw):
        if status == '[enabled]':
            optimizers.append(optimizer)

    return optimizers


def fetch_size(filename: str) -> int:
    """Fetch file size

    Args:
        filename (str): Path to file.

    Returns:
        int: File size in bytes.
    """
    return os.path.getsize(filename)
