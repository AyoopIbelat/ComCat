import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy import stats
from ase.units import Bohr, Hartree, Rydberg, fs
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from decimal import Decimal
from math import ceil, floor
from scipy import stats

# def rm_outlierz(df, zscore=2):
#     mask=np.zeros((df.shape[0],1))==0
#     for i in ['CO','CHO','OH','H']:
#         print(i)
#         # print(np.abs(stats.zscore(df[i])) < zscore)
#         df=df[(np.abs(stats.zscore(df[i])) < zscore)]
#         print(df.shape)
#     return df
def rm_outlier(df):
    for i in ['CO','CHO']:
        print(i)
        df=df[(df[i]<1) & (df[i]>-4)]
    return df
df100=pd.read_csv('100xy.csv',index_col=0)
df100=rm_outlier(df100)
df110=pd.read_csv('110xy.csv',index_col=0)
df110=rm_outlier(df110)

gco = np.linspace(-4.0, 1, 1000)
gcho = np.linspace(-4.0, 4, 1000)
GCO, GCHO = np.meshgrid(gco,gcho)
OP = (GCO + GCHO)

for i in range(0,len(gco)):
    for j in range(0,len(gcho)):
        bco = GCO[i,j]
        if bco < 0:
            OP[i,j] = GCHO[i,j]-bco
        if bco > 0 or bco == 0:
            OP[i,j] = GCHO[i,j]
OP=-OP

cp = plt.contourf(GCO, GCHO, OP,cmap = 'jet',levels = np.arange(-3, 2, 0.01),extend = 'both')
cbar = plt.colorbar(cp)
cbar.set_label(r'$U_{L}$ / V vs. RHE', fontsize=17)
v = np.linspace(-1, 0, 10, endpoint=True)

plt.xlabel("$n_1$")
plt.xlabel(r'$\Delta G_{CO}\ / \ eV$',fontsize=20)#{\epsilon$}')
plt.ylabel(r'$\Delta G_{CHO}\ / \ eV$',fontsize=20)#{\epsilon$}')

GCO_data = [0.004803165,-0.492671966,-0.773346328,-0.912723963,-0.990530916,-1.02196395,-1.107329304,-1.20046096,-1.192854582,-1.123581567] #Put CO binding energies here
GCHO_data  = [0.991206968,0.669408745,-0.140116072,0.169172151,0.061990269,-0.211758855,-0.164237195,-0.116722965,-0.235723975,-0.699790026] #Put CHO binding energies here


plt.scatter(GCO_data,GCHO_data,marker='^', c='r', alpha=0.8)
plt.scatter(df100.CO,df100.CHO, alpha=0.8,s=15,label='TMNs (100)')
plt.scatter(df110.CO,df110.CHO, alpha=0.8,s= 15,label='TMNs (110)')
plt.legend()

plt.savefig('Volcano.png',dpi = 400, bbox_inches='tight')
plt.show()
