#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 13:00:27 2023

@author: patricfernandez
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:52:19 2023

@author: patricfernandez
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 06:51:33 2023

@author: rz73
"""

# %% Import libraries

import pandas as pd
import numpy as np
from os import listdir
#import random
#import math
#from pysam import FastaFile
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy import signal
#import os
#from plotnine import *
import re

# %% Establishing path to the folder containing raw .csv data files

path_data_NGS = "/Users/patricfernandez/Documents/python/Pu-Seq_data_senataxin_run2/"

path_data_NGS = "/Users/patricfernandez/Documents/python/episilon figure try/"
path_data_NGS = "/Users/patricfernandez/Documents/python/epsilon2 figure try/"


path_data_NGS = "/Users/patricfernandez/Documents/python/idfromdelta/"
# %%

def call_peak(data, bin_size = 300):
   
    data = data.reset_index()
    data['peak_no'] = 1
    for i in range(1,len(data)):
        pos = data.iloc[i]["pos"]
        bpos = data.iloc[i-1]["pos"]
       # print(len(data))
        if i != len(data) and pos == bpos + bin_size:
            data.loc[i, 'peak_no'] = data.loc[i-1, 'peak_no']
               
        elif i == len(data) and pos == bpos + bin_size:
            data.loc[i, 'peak_no'] = data.loc[i-1, 'peak_no']
                   
        elif i == len(data) and pos != bpos + bin_size:
            data.loc[i, 'peak_no'] = data.loc[i-1, 'peak_no'] + 1
                       
        else:
            data.loc[i, 'peak_no'] = data.loc[i-1, 'peak_no'] + 1

    return data

# %% Importing raw Pu-Seq

list_files = list(filter(lambda file : ".csv" in file, listdir(path_data_NGS)))

data_raw = pd.DataFrame()

for file in list_files :
   
    tmp_data = pd.read_csv(path_data_NGS + file, sep = ",")
    tmp_data["strain"] = re.split('-|_|\\.', file)[0]
    tmp_data["condition"] = re.split('-|_|\\.', file)[1]
    tmp_data["pol"] = re.split('-|_|\\.', file)[2]
    tmp_data["strand"] = re.split('-|_|\\.', file)[4]
   
    data_raw = pd.concat([data_raw, tmp_data])
    del([tmp_data, file])
   
data_raw["experiment"] = np.where(data_raw["strain"].isin(["RZ259", "RZ261", "RZ263", "RZ265", "RZ267", "RZ269", "RZ271", "RZ273"]), 1,
                                np.where(data_raw["strain"].isin(["RZ260", "RZ262", "RZ264", "RZ266", "RZ268", "RZ270", "RZ272", "RZ274"]), 2, "error"))

data_raw["genotype"] = np.where(data_raw["strain"].isin(["RZ259", "RZ260", "RZ267", "RZ268"]), "WT",
                                np.where(data_raw["strain"].isin(["RZ261", "RZ262", "RZ269", "RZ270"]), "sen1_d",
                                         np.where(data_raw["strain"].isin(["RZ263", "RZ264", "RZ271", "RZ272"]), "dbl8_d",
                                                  np.where(data_raw["strain"].isin(["RZ265", "RZ266", "RZ273", "RZ274"]), "sen1dbl8_dd", "error"))))

data_raw["count"] = np.where(data_raw["count"] == 0, 1, data_raw["count"])

data_raw["count_norm"] = data_raw.groupby(
    ["experiment", "strain", "genotype", "condition", "pol", "strand"])["count"].transform(
        lambda x : x / sum(x))
       
# %% Calculating polymerase tracks

data_pol_track = data_raw.copy()

data_pol_track = data_pol_track.pivot_table(
   
    index = ["experiment", "strain", "genotype", "condition", "pol", "chro", "pos"],
    columns = "strand",
    values = "count").reset_index()

data_pol_track["track"] = (data_pol_track["f"] - data_pol_track["r"]) / (data_pol_track["f"] + data_pol_track["r"])

data_pol_track["track_smooth"] = data_pol_track.groupby(
    ["experiment", "strain", "genotype", "condition", "pol", "chro"])["track"].transform(
        lambda x : x.rolling(7, center = True).mean())

data_pol_track["track_diff"] = data_pol_track.groupby(
    ["experiment", "strain", "genotype", "condition", "pol", "chro"])["track_smooth"].transform(
        lambda x : -np.diff(x, append = 1, n = 1))
       
data_pol_track["diff_smooth"] = data_pol_track.groupby(
    ["experiment", "strain", "genotype", "condition", "pol", "chro"])["track_diff"].transform(
        lambda x : x.rolling(7, center = True).mean())

x = data_pol_track.loc[data_pol_track["pol"] == "e"]



# %% Averaging two biological repeats

data_pol_track_mean = data_pol_track.copy().groupby(
    ["genotype", "condition", "pol", "chro", "pos"],
    as_index = False).agg(
        track_mean = ("track_smooth", "mean"))
       
data_pol_track_mean["track_diff"] = data_pol_track_mean.groupby(
    ["genotype", "condition", "pol", "chro"])["track_mean"].transform(
        lambda x : -np.diff(x, append = 1, n = 1))
           
y = data_pol_track_mean.loc[data_pol_track_mean["pol"] == "e"]

# %% Calculating differences in termination

data_termination_change = data_pol_track_mean.copy().loc[data_pol_track_mean["pol"] == "e"]
data_termination_change["termination"] = np.where(data_termination_change["track_diff"] <= 0, 0, data_termination_change["track_diff"])
       
data_termination_change = data_termination_change.pivot(
   
    index = ["condition", "pol", "chro", "pos"],
    columns = "genotype",
    values = "termination").reset_index()

data_termination_change["sen1D"] = data_termination_change["sen1_d"] - data_termination_change["WT"]
data_termination_change["dbl8D"] = data_termination_change["dbl8_d"] - data_termination_change["WT"]
data_termination_change["sen1dbl8DD"] = data_termination_change["sen1dbl8_dd"] - data_termination_change["WT"]

data_termination_change = data_termination_change.drop(["WT", "sen1_d", "dbl8_d", "sen1dbl8_dd"], axis = 1)

data_termination_change = data_termination_change.melt(
   
    id_vars = ["condition", "pol", "chro", "pos"],
    var_name = "genotype",
    value_name = "termination_change")

data_termination_change["termination_change_smooth"] = data_termination_change.groupby(
    ["condition", "genotype", "pol", "chro"])["termination_change"].transform(
        lambda x : x.rolling(3, center = True).mean())

data_termination_change = data_termination_change.fillna(0)
       
data_termination_change = data_termination_change.sort_values(
            by = ["condition", "genotype", "chro", "pos"]).reset_index()


# Removing telomeres (chr1)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr1") & (data_termination_change["pos"] < 100000), 0,
                                                                np.where((data_termination_change["chro"] == "chr1") & (data_termination_change["pos"] > 5579133 - 100000), 0, data_termination_change["termination_change_smooth"]))

# Removing telomeres (chr2)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr2") & (data_termination_change["pos"] < 100000), 0,
                                                                np.where((data_termination_change["chro"] == "chr2") & (data_termination_change["pos"] > 4539804 - 100000), 0, data_termination_change["termination_change_smooth"]))

# Removing telomeres (chr3)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr3") & (data_termination_change["pos"] < 100000), 0,
                                                                np.where((data_termination_change["chro"] == "chr3") & (data_termination_change["pos"] > 2452883 - 100000), 0, data_termination_change["termination_change_smooth"]))

# Removing centromere (chr1)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr1") & (data_termination_change["pos"] < 3789421) & (data_termination_change["pos"] > 3753687), 0, data_termination_change["termination_change_smooth"])

# Removing centromere (chr2)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr2") & (data_termination_change["pos"] < 1644747) & (data_termination_change["pos"] > 1602418), 0, data_termination_change["termination_change_smooth"])

# Removing centromere (chr3)    
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr3") & (data_termination_change["pos"] < 1137003) & (data_termination_change["pos"] > 1070904), 0, data_termination_change["termination_change_smooth"])

# Removing sen1 locus
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr1") & (data_termination_change["pos"] < 3262405 + 1800) & (data_termination_change["pos"] > 3262405 - 1800), 0, data_termination_change["termination_change_smooth"])

# Removing dbl8 locus
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr2") & (data_termination_change["pos"] < 2557150 + 1800) & (data_termination_change["pos"] > 2557150 - 1800), 0, data_termination_change["termination_change_smooth"])

# Removing mating type locus
data_termination_change["termination_change_smooth"] = np.where((data_termination_change["chro"] == "chr2") & (data_termination_change["pos"] < 2140000) & (data_termination_change["pos"] > 2110000), 0, data_termination_change["termination_change_smooth"])

# %% Peak calling

data_termination_change_peaks = data_termination_change.copy().loc[data_termination_change["termination_change_smooth"] > 0]

data_termination_change_peaks = data_termination_change_peaks.groupby(
    ["condition", "genotype", "chro"], as_index = False).apply(lambda x: call_peak(x))

# %% Filtering Peaks

data_termination_change_peaks_filtered = data_termination_change_peaks.copy().groupby(
    ["condition", "genotype", "chro", "peak_no"], as_index = False).agg(      
        chro = ("chro", "first"),
        start_pos = ("pos", "first"),
        end_pos = ("pos", "last"),
        peak_area = ("termination_change_smooth", "sum"),
        peak_max = ("termination_change_smooth", "max"))
       
data_termination_change_peaks_filtered = data_termination_change_peaks_filtered.loc[
    data_termination_change_peaks_filtered["condition"] == "RhArepr"]
       
data_termination_change_peaks_filtered = data_termination_change_peaks_filtered.loc[
   
    (data_termination_change_peaks_filtered["peak_area"] > np.quantile(data_termination_change_peaks_filtered["peak_area"], 0.98)) &
    (data_termination_change_peaks_filtered["peak_max"] > np.quantile(data_termination_change_peaks_filtered["peak_max"], 0.98))]


#%%
'''
pos_start_series = data_termination_change_peaks.groupby('peak_no')['pos'].min()
pos_end_series = data_termination_change_peaks.groupby('peak_no')['pos'].max()


# Map the values to df based on 'peak_no'
data_termination_change_peaks_filtered['pos_start'] = data_termination_change_peaks_filtered['peak_no'].map(pos_start_series)
data_termination_change_peaks_filtered['pos_end'] = data_termination_change_peaks_filtered['peak_no'].map(pos_end_series)

# Print the result
print(data_termination_change_peaks_filtered)

'''





#%%

sensen = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'sen1D')]
dbldbl = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'dbl8D')]
dsds = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'sen1dbl8DD')]

sensen1 = sensen[(sensen['chro'] == 'I')]
sensen2 = sensen[(sensen['chro'] == 'II')]
sensen3 = sensen[(sensen['chro'] == 'III')]

dbldbl1 = dbldbl[(dbldbl['chro'] == 'I')]
dbldbl2 = dbldbl[(dbldbl['chro'] == 'II')]
dbldbl3 = dbldbl[(dbldbl['chro'] == 'III')]

dsds1 = dsds[(dsds['chro'] == 'I')]
dsds2 = dsds[(dsds['chro'] == 'II')]
dsds3 = dsds[(dsds['chro'] == 'III')]

dtcp1 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr1')]
dtcp2 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr2')]
dtcp3 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr3')]



#%%
#plotting here 

data_pol_track_I = data_pol_track[(data_pol_track['chro'] == 'chr1')]
#data_pol_track_Ii = data_pol_track[(data_pol_track['chro'] == 'chr1')]
data_pol_track_II = data_pol_track[(data_pol_track['chro'] == 'chr2')]
data_pol_track_III = data_pol_track[(data_pol_track['chro'] == 'chr3')]

data_pol_track_I_r = data_pol_track_I[(data_pol_track_I['condition'] == 'RhArepr')]

data_pol_track_I_wt = data_pol_track_I_r[(data_pol_track_I_r['genotype'] == 'WT')]
data_pol_track_I_sen = data_pol_track_I_r[(data_pol_track_I_r['genotype'] == 'sen1_d')]
data_pol_track_I_dbl = data_pol_track_I_r[(data_pol_track_I_r['genotype'] == 'dbl8_d')]
data_pol_track_I_ds = data_pol_track_I_r[(data_pol_track_I_r['genotype'] == 'sen1dbl8_dd')]



data_termination_change_I = data_termination_change[(data_termination_change['chro'] == 'chr1')]
data_termination_change_Ir = data_termination_change_I[(data_termination_change_I['condition'] == 'RhArepr')]

data_termination_change_I_sen = data_termination_change_Ir[(data_termination_change_Ir['genotype'] == 'sen1D')]
data_termination_change_I_dbl = data_termination_change_Ir[(data_termination_change_Ir['genotype'] == 'dbl8D')]
data_termination_change_I_ds = data_termination_change_Ir[(data_termination_change_Ir['genotype'] == 'sen1dbl8DD')]



def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    genesfor['Sbinpos'] = genesfor['start']/300
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*300 +150
    genesfor['Ebinpos'] = genesfor['end']/300
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*300 +150


    genesrev['Sbinposr'] = genesrev['end']/300
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*300 +150
    genesrev['Ebinposr'] = genesrev['start']/300
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*300 +150

    return genesfor, genesrev, genes

genesfor, genesrev, ggenes = Find("dbl8_stall_sites_direction.txt")
#controlfor, controlrev, ccontrol = Find('control_genes.txt')

gene1 = ggenes[(ggenes['chro'] == 'chr1')]
gene2 = ggenes[(ggenes['chro'] == 'chr2')]
gene3 = ggenes[(ggenes['chro'] == 'chr3')]


def Chromosome_plot_for (data1,data2, data3, data4, featurex, data2c, data3c, data4c, peak1):
    ff, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6,1, sharex=True)

    ax1.set_ylabel('f counts (e)')
    
    #ax2.set_title('sen1')
    ax1.plot(data1['pos'], data1['f'], color ='black', alpha=0.8, linewidth = 1)
    ax1.plot(data2['pos'], data2['f'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax1.plot(data3['pos'], data3['f'], color ='orange', alpha=0.8, linewidth = 1)
    ax1.plot(data4['pos'], data4['f'], color ='tomato', alpha=0.8, linewidth = 1)
    
    ax2.plot(data1['pos'], data1['track_smooth'], color ='black', alpha=0.8, linewidth = 1)
    ax2.plot(data2['pos'], data2['track_smooth'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.plot(data3['pos'], data3['track_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax2.plot(data4['pos'], data4['track_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    
    ax3.plot(data1['pos'], data1['diff_smooth'], color ='black', alpha=0.8, linewidth = 1)
    ax3.plot(data2['pos'], data2['diff_smooth'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax3.plot(data3['pos'], data3['diff_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax3.plot(data4['pos'], data4['diff_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    
    
  #  ax4.plot(data1c['pos'], data1c['termination_change'], color ='black', alpha=0.8)
    ax4.plot(data2c['pos'], data2c['termination_change_smooth'], color ='steelblue', alpha=0.8,linewidth = 1)
    ax4.plot(data3c['pos'], data3c['termination_change_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax4.plot(data4c['pos'], data4c['termination_change_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    ax4.set_ylim(-0.01,0.2)
    ax5.set_ylim(-0.01,0.2)
    #ax2.set_ylim(-5,5)
  #  ax2.set_ylabel('differential values (e)')
   # ax2.set_ylim(0,200)

    for fe in featurex.itertuples(index=False, name=None):
      #  ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                
            #    ax5.annotate(fe[0], xy = [fe[2],-0.15])    
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                ax6.set_ylabel('Gene annotations')
    
    for index, row in peak1.iterrows():
        if row['genotype'] == 'sen1D':

            start_pos = row['start_pos']
            end_pos = row['end_pos']
    
        # Filter df2 based on the conditions
            filtered_data = data2c[(data2c['pos'] >= start_pos-300) & (data2c['pos'] <= end_pos+300)]
         #   filtered_datadb = data3c[(data3c['pos'] >= start_pos-300) & (data3c['pos'] <= end_pos+300)]
          #  filtered_datads = data4c[(data4c['pos'] >= start_pos-300) & (data4c['pos'] <= end_pos+300)]
    
        # Plot the filtered data
            ax5.plot(filtered_data['pos'], filtered_data['termination_change_smooth'],color ='steelblue',)
         #   ax5.plot(filtered_datadb['pos'], filtered_datadb['termination_change_smooth'],color ='orange',)
          #  ax5.plot(filtered_datads['pos'], filtered_datads['termination_change_smooth'],color ='tomato', label=f'Segment {index + 1}')


        if row['genotype'] == 'dbl8D':
           start_pos = row['start_pos']
           end_pos = row['end_pos']
    
        # Filter df2 based on the conditions
           filtered_datab = data3c[(data3c['pos'] >= start_pos-600) & (data3c['pos'] <= end_pos+600)]
    
        # Plot the filtered data
           ax5.plot(filtered_datab['pos'], filtered_datab['termination_change_smooth'],color ='orange', label=f'Segment {index + 1}')

        if row['genotype'] == 'sen1dbl8DD':
           start_pos = row['start_pos']
           end_pos = row['end_pos']
    
        # Filter df2 based on the conditions
           filtered_datads = data4c[(data4c['pos'] >= start_pos-600) & (data4c['pos'] <= end_pos+600)]
    
        # Plot the filtered data
           ax5.plot(filtered_datads['pos'], filtered_datads['termination_change_smooth'],color ='tomato', label=f'Segment {index + 1}')
        

   
   
    return ff

chr1 = Chromosome_plot_for (data_pol_track_I_wt, data_pol_track_I_sen, data_pol_track_I_dbl, data_pol_track_I_ds, gene1,
                            data_termination_change_I_sen,data_termination_change_I_dbl, data_termination_change_I_ds, dtcp1)
    
#%% 


etcp1 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr1')]
etcp2 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr2')]
etcp3 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'chr3')]
 
def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 'chr1'
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 'chr2'
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 'chr3'
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]
    
    featrev.loc[featrev['Start position'] < 25, 'Start position'] = 25
     
    featfor['Sbinpos'] = featfor['Start position']/50
    featfor['Sbinpos'] = featfor['Sbinpos'].astype(int)
    featfor['Sbinpos'] = featfor['Sbinpos']*50 + 25
    featfor['Ebinpos'] = featfor['End position']/50
    featfor['Ebinpos'] = featfor['Ebinpos'].astype(int)
    featfor['Ebinpos'] = featfor['Ebinpos']*50 +25
    
    featrev['Sbinpos'] = featrev['End position']/50 
    featrev['Sbinpos'] = featrev['Sbinpos'].astype(int)
    featrev['Sbinpos'] = featrev['Sbinpos']*50 +25
    featrev['Ebinpos'] = featrev['Start position']/50
    featrev['Ebinpos'] = featrev['Ebinpos'].astype(int)
    featrev['Ebinpos'] = featrev['Ebinpos']*50 +25
    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')

feat1 = ffeat[(ffeat['chro'] == 'chr1')]
feat2 = ffeat[(ffeat['chro'] == 'chr2')]
feat3 = ffeat[(ffeat['chro'] == 'chr3')]


def assign_id_to_deltastall_peak(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['Start position'] <= replicon_row['start_pos'] and
                gene_row['End position'] >= replicon_row['start_pos']):
                matching_replicon_ids.append(replicon_row['peak_no'])
                gene_list.append(gene_row['Systematic ID'])
            elif (gene_row['Start position'] <= replicon_row['end_pos'] and
                gene_row['End position'] >= replicon_row['end_pos']):
                matching_replicon_ids.append(replicon_row['peak_no'])
                gene_list.append(gene_row['Systematic ID'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'Systematic ID': gene_list})
    replicon_df = replicon_df.merge(id_and_ID, on='peak_no', how='left')
  #  duplicate_counts = replicon_df['Systematic ID'].value_counts()
    df_no_duplicates = replicon_df.drop_duplicates(subset='peak_no', keep='first')
    

    return df_no_duplicates

etcp1 = assign_id_to_deltastall_peak(etcp1, feat1)
etcp2 = assign_id_to_deltastall_peak(etcp2, feat2)
etcp3 = assign_id_to_deltastall_peak(etcp3, feat3)





sensen = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'sen1D')]
dbldbl = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'dbl8D')]
dsds = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['genotype'] == 'sen1dbl8DD')]

sensen1 = sensen[(sensen['chro'] == 'chr1')]
sensen2 = sensen[(sensen['chro'] == 'chr2')]
sensen3 = sensen[(sensen['chro'] == 'chr3')]

dbldbl1 = dbldbl[(dbldbl['chro'] == 'chr1')]
dbldbl2 = dbldbl[(dbldbl['chro'] == 'chr2')]
dbldbl3 = dbldbl[(dbldbl['chro'] == 'chr3')]

dsds1 = dsds[(dsds['chro'] == 'chr1')]
dsds2 = dsds[(dsds['chro'] == 'chr2')]
dsds3 = dsds[(dsds['chro'] == 'chr3')]

dtcp1 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'I')]
dtcp2 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'II')]
dtcp3 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'III')]

#etcp1.to_csv('/Users/patricfernandez/Documents/python/epi_stalls_chr1.csv', index=False)
#etcp2.to_csv('/Users/patricfernandez/Documents/python/epi_stalls_chr2.csv', index=False)
#etcp3.to_csv('/Users/patricfernandez/Documents/python/epi_stalls_chr3.csv', index=False)

#so that hasn't worked, i will need to go by replicon, and common replicons 









