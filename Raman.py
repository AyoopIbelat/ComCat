import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
###DATA INPUT###
df1=pd.read_excel('Calf PSF 50 - 4.xlsx',index_col=None)
df2=pd.read_excel('CALF-20-old.xlsx')
df3=pd.read_excel('PSF Pellet.xlsx')
df1.columns=['x','y']
df2.columns=['x','y']
df3.columns=['x','y']

x1=np.array(df1.x)[102:1700]
y1=StandardScaler().fit_transform(df1)[:,1][102:1700]
x2=np.array(df2.x)[102:1700]
y2=StandardScaler().fit_transform(df2)[:,1][102:1700]
x3=np.array(df3.x)[102:1700]
y3=StandardScaler().fit_transform(df3)[:,1][102:1700]
###################################

ax=plt.gca()
ax.plot(x1,y1+20,'C0', label='CALF-20+10wt.% PSF')
ax.plot(x2,y2+10,'C1', label='CALF-20')
ax.plot(x3,y3,'C2', label='PSF')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax.set_yticks([])
ax.legend(frameon=False)
plt.xlim([249.1,3328.9])
plt.ylim([-1, 30])
plt.xlabel('Raman Shift $\mathregular{(cm^{-1})}$',fontweight='bold')
plt.ylabel('Intensity (a.u.)',fontweight='bold')
plt.savefig('CALF20_Raman.png',dpi=720)
plt.show()

