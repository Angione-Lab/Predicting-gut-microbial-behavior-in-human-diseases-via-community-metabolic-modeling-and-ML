# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:33:39 2023

@author: U0034546
"""

import pandas as pd
import glob
import os

covid = glob.glob('Covid/*.csv')
healthy = glob.glob('Healthy/*.csv')

samples = covid + healthy
exRxns = {}

for i in samples:
    df = pd.read_csv(i, index_col=0).T
    df['label'] = i.split("\\")[0]
    df.iloc[0,:].to_dict()
    exRxns[os.path.basename(i).split('_')[0]] = df.iloc[0,:].to_dict()


exRxn_df = pd.DataFrame.from_dict(exRxns).T.dropna(axis = 1)

exRxn_df.to_csv('ExchangeReactions.csv')


#%%
def extract_interaction_type(cohert, label):
    covid_samples = {}
    j = 0
    for i in cohert:
        print(i)
        sample = pd.read_csv(i, sep='\t', header=0, index_col=0)
        
        sample = pd.read_csv(i, sep='\t', header=0, index_col=0)
        
        interactions = (sample[' TypeOfInteraction'].value_counts()/sample.shape[0]).to_dict()
        
        mean_flux = sample[[' GRSpeciesAFull ',' GRSpeciesBFull ', ' GRASolo ', ' GRBSolo ']].mean(axis = 0).to_dict()
        mean_flux['label'] = label
        interactions.update(mean_flux)
        item = {j:interactions}
        covid_samples = {**covid_samples, **item}
        
        j +=1
    return covid_samples



covid_cohert = glob.glob('Covidpolished/*/*.tsv')    
healthy_cohert = glob.glob('hcdpolished/*/*.tsv')       

covid_item = pd.DataFrame(extract_interaction_type(covid_cohert, 'C')).T
healthy_item = pd.DataFrame(extract_interaction_type(healthy_cohert, 'H')).T

all_samples = pd.concat([covid_item, healthy_item], axis = 0).fillna(0).reset_index(drop=True)

all_samples.to_csv('covid_healthy.csv')
