from setuptools import setup

setup(name='yatuner',
      version='0.0.1',
      author='Synodic Month, Juni May',
      author_email='synodic_month@163.com, juni_may@outlook.com',
      description='Yet another auto tuner for compilers.',
      long_description='README.md',
      packages=['yatuner'],
      requires=[
          'GPyOpt', 'cython', 'GPy', 'numpy', 'matplotlib', 'scipy', 'rich', 'seaborn'
      ],
      license='Mulan PSL v2',
      entry_points={'console_scripts': ['yatuner=yatuner.__main__:main']})
