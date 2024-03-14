import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

data = pd.read_csv('../TESTES_REAIS/yatuner.db_1nd_linUCB/result.csv')

O0 = data['O0'][::]
O1 = data['O1'][::]
O2 = data['O2'][::]
O3 = data['O3'][::]
Os = data['Os'][::]
Ofast = data['Ofast'][::]
Optimizers = data['Optimizers'][::]
Parameters = data['Parameters'][::]

 
df = pd.DataFrame([[min(O0), min(O1), min(O2), min(O3), min(Os), min(Ofast), min(Optimizers), min(Parameters)], 
	[max(O0), max(O1), max(O2), max(O3), max(Os), max(Ofast), max(Optimizers), max(Parameters)]], 
	columns = ['O0', 'O1', 'O2', 'O3', 'Os', 'Ofast', 'Optimizers', 'Parameters'])
# writing data frame to a CSV file
df.to_csv('primeiro/tabelaum.csv')	



plt.style.use('seaborn')
pd_data = pd.read_csv("primeiro/tabelaum.csv")
plt.clf()
plt.figure(figsize=(10, 5), dpi=100)
sns.violinplot(
            data=pd_data,
            orient='vertical',
            palette='Set3',
            width=0.9)
plt.title('Time Comparison')
plt.ylabel('Optimization_methods')


'''
eixo_y = [3712654.733,2652017.221,2573221.866,
2546447.99,2576195.949,2444390.514,2377303.433,
1351549.485,3818220.612,2689851.039,2647972.404,
2717234.943,2617238.142,2856473.672,2472749.779,1377074.064]

plt.yticks(eixo_y)
'''

plt.xlabel('Time/Tick - Lower the better')
plt.grid()
plt.savefig("primeiro/1stgraph.png")