import pandas as pd
import numpy as np

def ads_energy(mp, ads, energy,termination):
    energy=float(energy)
    global ref
    global clean
    [co]=(ref['CO'])
    [h2]=(ref['H2'])
    [h2o]=(ref['H2O'])
    cl=float(clean[(clean['mp-id']==mp) & (clean['Termination']==termination)]['Final Energy'])
    if ads=='CHO':
        x= (energy-co-h2/2-cl)
    elif ads=='CO':
        x= energy-co-cl
    elif ads=='OH':
        x= energy+h2/2-h2o-cl
    elif ads=='H':
        x= energy-h2/2-cl
    else:
        x= np.nan
    return x
    
ref=pd.read_csv('./result/free_ads.csv',index_col=0).transpose()
data=pd.read_csv('Data.csv',index_col=0)
clean=data[data['Adsorbate']=='clean'][['mp-id','Termination','Final Energy']]

data['Adsorption Energy']=data.apply(lambda row: ads_energy(str(row['mp-id']),
                                                            str(row['Adsorbate']),
                                                            str(row['Final Energy']),
                                                            str(row['Termination'])),
                                      axis=1)
data=data[['mp-id', 'Adsorbate', 'Termination','Site#', 'Convergence','Final Energy', 
           'Dissociation', 'Adsorbate-Surface Connection', 'Adsorption Energy']]
data.set_index(['mp-id','Termination','Adsorbate'], inplace=True)
data.to_excel('Data_with_adsE.xlsx')