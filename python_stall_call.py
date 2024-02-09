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

@author: rz73 and jf383
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
        lambda x : np.diff(x, append = 1, n = 1))
       
data_pol_track["diff_smooth"] = data_pol_track.groupby(
    ["experiment", "strain", "genotype", "condition", "pol", "chro"])["track_diff"].transform(
        lambda x : x.rolling(7, center = True).mean())

x = data_pol_track.loc[data_pol_track["pol"] == "d"]

#%%
#plotting here 

data_pol_track_I = data_pol_track[(data_pol_track['chro'] == 'I')]
data_pol_track_II = data_pol_track[(data_pol_track['chro'] == 'II')]
data_pol_track_III = data_pol_track[(data_pol_track['chro'] == 'III')]

data_pol_track_I_wt = data_pol_track_I[(data_pol_track_I['genotype'] == 'WT')]
data_pol_track_I_sen = data_pol_track_I[(data_pol_track_I['genotype'] == 'sen1_d')]
data_pol_track_II_dbl = data_pol_track_I[(data_pol_track_I['genotype'] == 'dbl8_d')]
data_pol_track_III_ds = data_pol_track_I[(data_pol_track_I['genotype'] == 'sen1dbl8_dd')]



data_termination_change_I = data_termination_change[(data_termination_change['chro'] == 'I')]
data_termination_change_I_sen = data_termination_change_I[(data_termination_change_I['genotype'] == 'sen1D')]
data_termination_change_I_dbl = data_termination_change_I[(data_termination_change_I['genotype'] == 'dbl8D')]
data_termination_change_I_ds = data_termination_change_I[(data_termination_change_I['genotype'] == 'sen1dbl8DD')]


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

    ax1.set_ylabel('f counts (d)')
    
    #ax2.set_title('sen1')
    ax1.plot(data1['pos'], data1['f'], color ='black', alpha=0.8, linewidth = 1)
    ax1.plot(data2['pos'], data2['f'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax1.plot(data3['pos'], data3['f'], color ='orange', alpha=0.8, linewidth = 1)
    ax1.plot(data4['pos'], data4['f'], color ='tomato', alpha=0.8, linewidth = 1)
    
    ax2.plot(data1['pos'], data1['track_smooth'], color ='black', alpha=0.8, linewidth = 1)
    ax2.plot(data2['pos'], data2['track_smooth'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.plot(data3['pos'], data3['track_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax2.plot(data4['pos'], data4['track_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    ax2.set_ylabel('fork direction (d)')
    
    ax3.plot(data1['pos'], data1['diff_smooth'], color ='black', alpha=0.8, linewidth = 1)
    ax3.plot(data2['pos'], data2['diff_smooth'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax3.plot(data3['pos'], data3['diff_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax3.plot(data4['pos'], data4['diff_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    ax3.set_ylabel('smoo_diff(d)')
    
  #  ax4.plot(data1c['pos'], data1c['termination_change'], color ='black', alpha=0.8)
    ax4.plot(data2c['pos'], data2c['termination_change_smooth'], color ='steelblue', alpha=0.8,linewidth = 1)
    ax4.plot(data3c['pos'], data3c['termination_change_smooth'], color ='orange', alpha=0.8, linewidth = 1)
    ax4.plot(data4c['pos'], data4c['termination_change_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    ax4.set_ylim(-0.01,0.2)
    ax5.set_ylim(-0.01,0.2)
    ax4.set_ylabel('termination change (d)')
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
            filtered_data = data2c[(data2c['pos'] >= start_pos-600) & (data2c['pos'] <= end_pos+600)]
            filtered_datadb = data3c[(data3c['pos'] >= start_pos-300) & (data3c['pos'] <= end_pos+300)]
            filtered_datads = data4c[(data4c['pos'] >= start_pos-300) & (data4c['pos'] <= end_pos+300)]
    
        # Plot the filtered data
            ax5.plot(filtered_data['pos'], filtered_data['termination_change_smooth'],color ='steelblue',)
          #  ax5.plot(filtered_datadb['pos'], filtered_datadb['termination_change_smooth'],color ='orange',)
           # ax5.plot(filtered_datads['pos'], filtered_datads['termination_change_smooth'],color ='tomato', label=f'Segment {index + 1}')


        elif row['genotype'] == 'dbl8D':
            start_pos = row['start_pos']
            end_pos = row['end_pos']
    
        # Filter df2 based on the conditions
            filtered_datab = data3c[(data3c['pos'] >= start_pos-300) & (data3c['pos'] <= end_pos+300)]
    
        # Plot the filtered data
            ax5.plot(filtered_datab['pos'], filtered_datab['termination_change_smooth'],color ='orange', label=f'Segment {index + 1}')

        elif row['genotype'] == 'sen1dbl8DD':
            start_pos = row['start_pos']
            end_pos = row['end_pos']
        

            filtered_datad = data4c[(data4c['pos'] >= start_pos-300) & (data4c['pos'] <= end_pos+300)]
    
        # Plot the filtered data
            ax5.plot(filtered_datad['pos'], filtered_datad['termination_change_smooth'],color ='tomato', label=f'Segment {index + 1}')

   
   
    return ff

chr1 = Chromosome_plot_for (data_pol_track_I_wt, data_pol_track_I_sen, data_pol_track_II_dbl, data_pol_track_III_ds, gene1,
                            data_termination_change_I_sen,data_termination_change_I_dbl, data_termination_change_I_ds, dtcp1)
    
    
    

# %% Averaging two biological repeats

data_pol_track_mean = data_pol_track.copy().groupby(
    ["genotype", "condition", "pol", "chro", "pos"],
    as_index = False).agg(
        track_mean = ("track_smooth", "mean"))
       
data_pol_track_mean["track_diff"] = data_pol_track_mean.groupby(
    ["genotype", "condition", "pol", "chro"])["track_mean"].transform(
        lambda x : np.diff(x, append = 1, n = 1))
           
y = data_pol_track_mean.loc[data_pol_track_mean["pol"] == "d"]

# %% Calculating differences in termination

data_termination_change = data_pol_track_mean.copy().loc[data_pol_track_mean["pol"] == "d"]
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

dtcp1 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'I')]
dtcp2 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'II')]
dtcp3 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'III')]



#%%%

 
def Create_df(df,dr,ef,er,w,df2,dr2,ef2,er2):
    df = pd.read_csv(df) # Reads a csv file
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file
    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef = pd.read_csv(ef, usecols=[2])
    ef['count'].replace(0,1, inplace = True)   
    ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
    er = pd.read_csv(er, usecols=[2])
    er['count'].replace(0,1, inplace = True)
    er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    
    df2 = pd.read_csv(df2) # Reads a csv file
    df2['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df2.rename(columns = {"count" : "df2_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr2 = pd.read_csv(dr2, usecols=[2]) # Read only the counts column from the next file
    dr2['count'].replace(0,1, inplace = True)
    dr2.rename(columns = {'count' : 'dr2_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef2 = pd.read_csv(ef2, usecols=[2])
    ef2['count'].replace(0,1, inplace = True)   
    ef2.rename(columns = {'count' : 'ef2_count'}, inplace = True)
    
    er2 = pd.read_csv(er2, usecols=[2])
    er2['count'].replace(0,1, inplace = True)
    er2.rename(columns = {'count' : 'er2_count'}, inplace = True)
    
    all_data = pd.concat([df, dr, ef, er], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes
    all_data2 = pd.concat([df2, dr2, ef2, er2], axis=1, join='outer')
  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
    all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
    all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
    all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
    all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])
    
    all_data2['norm_df'] = all_data2['df2_count']/all_data2['df2_count'].sum()
    all_data2['norm_dr'] = all_data2['dr2_count']/all_data2['dr2_count'].sum()
    all_data2['norm_ef'] = all_data2['ef2_count']/all_data2['ef2_count'].sum()
    all_data2['norm_er'] = all_data2['er2_count']/all_data2['er2_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data2['ratio_delta_f'] = all_data2['norm_df']/(all_data2['norm_df'] + all_data2['norm_ef'])
    all_data2['ratio_delta_r'] = all_data2['norm_dr']/(all_data2['norm_dr'] + all_data2['norm_er'])
    all_data2['ratio_epsilon_f'] = all_data2['norm_ef']/(all_data2['norm_ef'] + all_data2['norm_df'])
    all_data2['ratio_epsilon_r'] = all_data2['norm_er']/(all_data2['norm_er'] + all_data2['norm_dr'])


    
#INSTEAD, I will try okazaki method for rightward forks
    all_data['right_forks'] = (all_data['norm_ef'] - all_data['norm_er'])/ (all_data['norm_ef'] + all_data['norm_er'])
    all_data2['right_forks'] = (all_data2['norm_ef'] - all_data2['norm_er'])/ (all_data2['norm_ef'] + all_data2['norm_er'])
    all_data['av_right_forks'] = (all_data['right_forks'] + all_data2['right_forks'])/2
    
    all_data['smoo_right_forks'] = all_data['av_right_forks'].rolling(window = w, center=True).mean()
    #what i need to do now is try this without the smoothing 
   # all_data['smoo_right_forks'] = all_data['right_forks']

    #INSTEAD
    ser = all_data['av_right_forks']
    ser = ser.diff()


    ser = ser*-1
    


    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser
    all_data['smoo_differentials'] = all_data['differentials'].rolling(window = w, center =True).mean()
    all_data['smoo_differentials'].clip(0, 1, inplace=True)

# create three separate dataframes for the three chromosomes and return these too
 
    chrI = all_data.loc[all_data['chro']=='chr1']
    chrII = all_data.loc[all_data['chro']=='chr2']
    chrIII = all_data.loc[all_data['chro']=='chr3']
    return chrI, chrII, chrIII, all_data


 
ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 3, 'RZ260-RhArepr_d.e1.f-w300.count.csv', 'RZ260-RhArepr_d.e1.r-w300.count.csv', 'RZ268-RhArepr_e.e1.f-w300.count.csv', 'RZ268-RhArepr_e.e1.r-w300.count.csv')
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, seen1chr3, esen1 = Create_df ('RZ261-RhArepr-d.e1.f-w300.count.csv', 'RZ261-RhArepr-d.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 3, 'RZ262-RhArepr_d.e1.f-w300.count.csv', 'RZ262-RhArepr_d.e1.r-w300.count.csv', 'RZ270-RhArepr_e.e1.f-w300.count.csv', 'RZ270-RhArepr_e.e1.r-w300.count.csv')
edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('RZ263-RhArepr-d.e1.f-w300.count.csv', 'RZ263-RhArepr-d.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 3, 'RZ264-RhArepr_d.e1.f-w300.count.csv', 'RZ264-RhArepr_d.e1.r-w300.count.csv', 'RZ272-RhArepr_e.e1.f-w300.count.csv', 'RZ272-RhArepr_e.e1.r-w300.count.csv')
edschr1, edschr2, edschr3, eds = Create_df('RZ265-RhArepr-d.e1.f-w300.count.csv', 'RZ265-RhArepr-d.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 3, 'RZ266-RhArepr_d.e1.f-w300.count.csv', 'RZ266-RhArepr_d.e1.r-w300.count.csv', 'RZ274-RhArepr_e.e1.f-w300.count.csv', 'RZ274-RhArepr_e.e1.r-w300.count.csv')



def difbaby (data, wt):
    data['smoo_differentials'] = data['smoo_differentials'] - wt['smoo_differentials']
  #  data['smoo_differentials'].clip(0, 1, inplace=True)
    data = data[(data['smoo_differentials'] > 0)]
    
    data = data.dropna()
    
   # data = data.loc[data['differentials'] > 0]
   # data = data.dropna()
    #templist = np.array([])
    
 #   for i in range(1,len(data['differentials'])):
  #      if data['differentials'].iloc[i] != 0: 
   #         templist = np.append(templist, data['differentials'].iloc[i])
    
 #   cutoff = np.quantile(templist, 0.98)
  #  for i in range(1, len(data['differentials'])-1):
   #     if data['differentials'].iloc[i] != 0:
    #        if data['differentials'].iloc[i]< cutoff:
     #           data['differentials'].iloc[i] = 0
    return data

esen1chr1 = difbaby(esen1chr1, ewtchrI)
esen1chr2 = difbaby(esen1chr2, ewtchrII)
esen1chr3 = difbaby(seen1chr3, ewtchrIII)

edbl8chr1 = difbaby(edbl8chr1, ewtchrI)
edbl8chr2 = difbaby(edbl8chr2, ewtchrII)
edbl8chr3 = difbaby(edbl8chr3, ewtchrIII)

edschr1 = difbaby(edschr1, ewtchrI)
edschr2 = difbaby(edschr2, ewtchrII)
edschr3 = difbaby(edschr3, ewtchrIII)

####let's do delta forks
 
def Create_df(df,dr,w):
    df = pd.read_csv(df) # Reads a csv file
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file
    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
 #   ef = pd.read_csv(ef, usecols=[2])
  #  ef['count'].replace(0,1, inplace = True)   
  #  ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
 #   er = pd.read_csv(er, usecols=[2])
  #  er['count'].replace(0,1, inplace = True)
  #  er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    all_data = pd.concat([df, dr], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes
    
  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    print(all_data['df_count'].sum())
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    print(all_data['dr_count'].sum())
 #   all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
  #  all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()


# Next we calculate the ratios for each strand and assign to a new colunm
 #   all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
  #  all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
  #  all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
   # all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])


# Now we have column for pol delta useage for the duplex
   # all_data['d_usage'] = (all_data['ratio_delta_f'] + all_data['ratio_delta_r']) / 2

# now we a column for the percentage of right-moving forks
   # all_data['right_forks']  = all_data['ratio_epsilon_f']*2 - 1
    
#INSTEAD, I will try okazaki method for rightward forks
    all_data['right_forks'] = (all_data['norm_df'] - all_data['norm_dr'])/ (all_data['norm_df'] + all_data['norm_dr'])
# now we will produce a new colum for each sliding window average for each of the calculated columns
# Note: centre = True, means that the data is summed from both left and right. False means its the sum of the last of the number of values.
    
    #all_data['smoo_ratio_d_f'] = all_data['ratio_delta_f'].rolling(window = w, center=True).mean()
   # all_data['smoo_ratio_d_r'] = all_data['ratio_delta_r'].rolling(window = w, center=True).mean()
#    all_data['smoo_ratio_e_f'] = all_data['ratio_epsilon_f'].rolling(window = w, center=True).mean()
#    all_data['smoo_ratio_e_r'] = all_data['ratio_epsilon_r'].rolling(window = w, center=True).mean()
   # all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = w, center=True).mean()
   
    all_data['smoo_right_forks'] = all_data['right_forks'].rolling(window = w, center=True).mean()
    #what i need to do now is try this without the smoothing 
   # all_data['smoo_right_forks'] = all_data['right_forks']
    

# now create a differential collunm and a discreate orgin (centre = 1 bin) column.
# Note the script removes differentials not present in both strands and which are singular 
    
    # take unsmoothed pol usage data into two pandas arrays
 #   ser_etop = all_data['ratio_epsilon_f']
  #  ser_ebot = all_data['ratio_epsilon_r']
    
    #INSTEAD
    ser = all_data['smoo_right_forks']
    ser = ser.diff()

    # Calculate the differentials 
   # ser_topd = ser_etop.diff()
   # ser_botd = ser_ebot.diff()
    # Reverse the sign on the bottom strand differentials to match top strand
   # ser_botd = ser_botd*-1

    # curently dont use a roling average, but here they are if needed
    #ser_topd = ser_topd.rolling(window = 3, center=True).mean()
    #ser_botd = ser_botd.rolling(window = 3, center=True).mean()
    #ser = ser.rolling(window = 3, center=True).mean() ###################
    ser = ser*-1

    # Removing all the negative valuse to zero
  #  ser_topd.clip(0, 1, inplace=True)
  #  ser_botd.clip(0, 1, inplace=True)
    #INSTEAD
    #ser.clip(0, 1, inplace=True)
    
    # If the value in one is zero, then seting it in zero in both datasets (this removes a lot of noise)
 #   for i in range(len(ser_topd)):
 #       if ser_topd.iloc[i] == 0:
#            ser_botd.iloc[i] = 0
#    for i in range(len(ser_botd)):
 #       if ser_botd.iloc[i] == 0:
  #          ser_topd.iloc[i] = 0
    
    # Now we add the two things together and divide by two - i.e we have an average origin activity
  #  ser_cumd = (ser_topd + ser_botd)/2     

    # Now we want to calculate the quantile of all the positive data.
#    templist = np.array([])
  #  for i in range(1,len(ser_cumd)):
 #       if ser_cumd.iloc[i] != 0: 
 #           templist = np.append(templist, ser_cumd.iloc[i])
    # set a cutoff threshold at of the top 10% of all positive values
 #   cutoff = np.quantile(templist, 0.9)

    # now if a single +ve value (i.e at least one zero either side) is below this cutoff, then set to zero. This again removes noise.
 #   for i in range(1, len(ser_cumd)-1):
#        if ser_cumd.iloc[i] != 0:
 #           if ser_cumd.iloc[i-1] == 0 and ser_cumd.iloc[i+1] == 0 and ser_cumd.iloc[i]< cutoff:
 #               ser_cumd.iloc[i] = 0

    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser

    # Next we want to position "zonal" origins to a single point (bin) and give them a cumulative value (efficiency of the zone)
#    start=0
 #   stop=0
 #   origins = ser_cumd.clip(0, 0) # creates a new pd.series of zero values same length of ser_cumd
 #   for i in range(len(ser_cumd)):
 #       if i <= stop: # simply prevents itterating over the non-zero block
  #          continue # continue goes back to the for loop
   #     else:
    #        if ser_cumd.iloc[i] != 0:
     #           start = i        
      #          for z in range(start+1, len(ser_cumd)): # find the ned of the non-zero block
       #             if ser_cumd.iloc[z] == 0:
        #                stop = z
         #               tot = 0
          #              tem = 0
           #             tem_loc = 0
            #            for v in range(i,z): # adds the non-zero block numbers together and identifies the location of the highest value
             #               tot = tot + ser_cumd.iloc[v]
              #              if ser_cumd.iloc[v] > tem:
               #                 tem = ser_cumd.iloc[v]
                #                tem_loc = v
                 #       origins.iloc[tem_loc] = tot # adds the total of the non-zero block to the bin with the highest individula differential
                  #      break

  #  all_data['origin'] = origins

# create three separate dataframes for the three chromosomes and return these too

    all_data.loc[all_data['chro'] == "chr1", 'chro'] = 'I'
    all_data.loc[all_data['chro'] == "chr2", 'chro'] = 'II'
    all_data.loc[all_data['chro'] == "chr3", 'chro'] = 'III'
 
    chrI = all_data.loc[all_data['chro']=='I']
    chrII = all_data.loc[all_data['chro']=='II']
    chrIII = all_data.loc[all_data['chro']=='III']
    return chrI, chrII, chrIII, all_data



#wtchrI, wtchrII, wtchrIII = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)



#    mat_loc = all_data.loc[all_data['chro']=='mating_type_region']
#    mito = all_data.loc[all_data['chro']=='mitochondrial']
#    tel_gap = all_data.loc[all_data['chro']=='chr_II_telomeric_gap']
    
# Note: the %s sign allows us to asiign the variable outputname within a string for multiple file names
#    all_data.to_pickle('./%s.pkl' %(outputname))
#    chrI.to_pickle('./%s_chrI.pkl' %(outputname))
#    chrII.to_pickle('./%s_chrII.pkl' %(outputname))
#    chrIII.to_pickle('./%s_chrIII.pkl' %(outputname))
#    mat_loc.to_pickle('./%s_mat_locus.pkl' %(outputname))
#    mito.to_pickle('./%s_mito.pkl' %(outputname))
#    tel_gap.to_pickle('./%s_tel_gap.pkl' %(outputname))
    
#delta
ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df('Jo_wt_YE.e1.f-w300.count.csv', 'Jo_wt_YE.e1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, esen1chr3, esen1 = Create_df ('Jo_sen1d_YE.e1.f-w300.count.csv', 'Jo_sen1d_YE.e1.r-w300.count.csv', 3)
edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('Jo_dbl8_YE.e1.f-w300.count.csv', 'Jo_dbl8_YE.e1.r-w300.count.csv', 3)
edschr1, edschr2, edschr3, eds = Create_df('Jo_sen1_dbl8_YE.e1.f-w300.count.csv', 'Jo_sen1_dbl8_YE.e1.r-w300.count.csv', 3)

#epsilon
ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df( 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, seen1chr3, esen1 = Create_df ('RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 3)
edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 3)
edschr1, edschr2, edschr3, eds = Create_df('RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 3)

#epi two 
ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df( 'RZ268-RhArepr_e.e1.f-w300.count.csv', 'RZ268-RhArepr_e.e1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, seen1chr3, esen1 = Create_df ('RZ270-RhArepr_e.e1.f-w300.count.csv', 'RZ270-RhArepr_e.e1.r-w300.count.csv', 3)
edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('RZ272-RhArepr_e.e1.f-w300.count.csv', 'RZ272-RhArepr_e.e1.r-w300.count.csv', 3)
edschr1, edschr2, edschr3, eds = Create_df('RZ274-RhArepr_e.e1.f-w300.count.csv', 'RZ274-RhArepr_e.e1.r-w300.count.csv', 3)

 
ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, seen1chr3, esen1 = Create_df ('RZ261-RhArepr-d.e1.f-w300.count.csv', 'RZ261-RhArepr-d.e1.r-w300.count.csv', 3)
edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('RZ263-RhArepr-d.e1.f-w300.count.csv', 'RZ263-RhArepr-d.e1.r-w300.count.csv', 3)
edschr1, edschr2, edschr3, eds = Create_df('RZ265-RhArepr-d.e1.f-w300.count.csv', 'RZ265-RhArepr-d.e1.r-w300.count.csv', 3)


 
#ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df('Jo_wt_YE.d1.f-w300.count.csv', 'Jo_wt_YE.d1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
#esen1chr1, esen1chr2, esen1chr3, esen1 = Create_df ('Jo_sen1d_YE.d1.f-w300.count.csv', 'Jo_sen1d_YE.d1.r-w300.count.csv', 3)
#edbl8chr1, edbl8chr2, edbl8chr3, edbl8 = Create_df('Jo_dbl8_YE.d1.f-w300.count.csv', 'Jo_dbl8_YE.d1.r-w300.count.csv', 3)
#edschr1, edschr2, edschr3, eds = Create_df('Jo_sen1_dbl8_YE.d1.f-w300.count.csv', 'Jo_sen1_dbl8_YE.d1.r-w300.count.csv', 3)


def difbaby (data, wt):
    data['differentials'] = -data['differentials']
    data = data[(data['differentials'] > 0)]
    data['differentials'] = data['differentials'] - wt['differentials']
    
  #  data = data[(data['differentials'] > 0)]
    
    data = data.dropna()
    return data

esen1chr1 = difbaby(esen1chr1, ewtchrI)
esen1chr2 = difbaby(esen1chr2, ewtchrII)
esen1chr3 = difbaby(esen1chr3, ewtchrIII)

edbl8chr1 = difbaby(edbl8chr1, ewtchrI)
edbl8chr2 = difbaby(edbl8chr2, ewtchrII)
edbl8chr3 = difbaby(edbl8chr3, ewtchrIII)

edschr1 = difbaby(edschr1, ewtchrI)
edschr2 = difbaby(edschr2, ewtchrII)
edschr3 = difbaby(edschr3, ewtchrIII)



def fill_missing_positions(df, p, pos_column='pos', columns_to_fill_zero=None, bin_size=300):
    if columns_to_fill_zero is None:
        columns_to_fill_zero = []

    # Create a DataFrame with the complete range of positions
    complete_range = pd.DataFrame({pos_column: range(df[pos_column].min(), df[pos_column].max() + 1, bin_size)})

    # Merge the original DataFrame with the complete range DataFrame
    merged_df = pd.merge(complete_range, df, on=pos_column, how='left')

    # Fill NaN values with 0 for the specified columns
    merged_df[columns_to_fill_zero] = merged_df[columns_to_fill_zero].fillna(0)
    merged_df['chro'].fillna(p, inplace=True)

    return merged_df



# Example usage:
# Assuming df is your DataFrame
esen1chr1 = fill_missing_positions(esen1chr1, "I", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
esen1chr2 = fill_missing_positions(esen1chr2, "II", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
esen1chr3 = fill_missing_positions(esen1chr3, "III", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])

edbl8chr1 = fill_missing_positions(edbl8chr1, "I", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
edbl8chr2 = fill_missing_positions(edbl8chr2, "II", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
edbl8chr3 = fill_missing_positions(edbl8chr3, "III", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])

edschr1 = fill_missing_positions(edschr1, "I", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
edschr2 = fill_missing_positions(edschr2, "II", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
edschr3 = fill_missing_positions(edschr3, "III", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])


ewtchrI = fill_missing_positions(ewtchrI, "I", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
ewtchrII = fill_missing_positions(ewtchrII, "II", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])
ewtchrIII = fill_missing_positions(ewtchrIII, "III", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'norm_df', 'norm_dr',
                                                                               'right_forks', 'smoo_right_forks', 'differentials'])







sen1l = [esen1chr1, esen1chr2, esen1chr3]
esen1 = pd.concat(sen1l,ignore_index=True)

dbl8l = [edbl8chr1, edbl8chr2, edbl8chr3]
edbl8 = pd.concat(dbl8l ,ignore_index=True)


dsl = [edschr1, edschr2, edschr3]
eds = pd.concat(dsl,ignore_index=True)


ewtchrI = ewt[(ewt['chro'] == 'I')]
ewtchrII = ewt[(ewt['chro'] == 'II')]
ewtchrIII = ewt[(ewt['chro'] == 'III')]

esen1chr1 = esen1[(esen1['chro'] == 'I')]
esen1chr2 = esen1[(esen1['chro'] == 'II')]
esen1chr3 = esen1[(esen1['chro'] == 'III')]

edbl8chr1 = edbl8[(edbl8['chro'] == 'I')]
edbl8chr2 = edbl8[(edbl8['chro'] == 'II')]
edbl8chr3 = edbl8[(edbl8['chro'] == 'III')]

edschr1 = eds[(eds['chro'] == 'I')]
edschr2 = eds[(eds['chro'] == 'II')]
edschr3 = eds[(eds['chro'] == 'III')]






def Create_df(df,dr,ef,er,w):
    df = pd.read_csv(df) # Reads a csv file

    df['chro'] = df['chro'].replace({"I": "chr1", "II": "chr2", "III": "chr3", 'mitochondrial': 'mito'})
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file

    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef = pd.read_csv(ef, usecols=[2])
    ef['count'].replace(0,1, inplace = True)   
    ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
    er = pd.read_csv(er, usecols=[2])
    er['count'].replace(0,1, inplace = True)
    er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    all_data = pd.concat([df, dr, ef, er], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes
    
    

    

  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
    all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
    all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
    all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
    all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])


# Now we have column for pol delta useage for the duplex
    all_data['d_usage'] = (all_data['ratio_delta_f'] + all_data['ratio_delta_r']) / 2

# now we a column for the percentage of right-moving forks
    all_data['right_forks']  = all_data['ratio_epsilon_f']*2 - 1
    
    
# now we will produce a new colum for each sliding window average for each of the calculated columns
# Note: centre = True, means that the data is summed from both left and right. False means its the sum of the last of the number of values.
    
    all_data['smoo_ratio_d_f'] = all_data['ratio_delta_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_d_r'] = all_data['ratio_delta_r'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_f'] = all_data['ratio_epsilon_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_r'] = all_data['ratio_epsilon_r'].rolling(window = w, center=True).mean()
    all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = w, center=True).mean()
    all_data['smoo_right_forks'] = all_data['right_forks'].rolling(window = w, center=True).mean()
    

# now create a differential collunm and a discreate orgin (centre = 1 bin) column.
# Note the script removes differentials not present in both strands and which are singular 
    
    # take unsmoothed pol usage data into two pandas arrays
    ser_etop = all_data['ratio_epsilon_f']
    ser_ebot = all_data['ratio_epsilon_r']

    # Calculate the differentials 
    ser_topd = ser_etop.diff()
    ser_botd = ser_ebot.diff()
    # Reverse the sign on the bottom strand differentials to match top strand
    ser_botd = ser_botd*-1

    # curently dont use a roling average, but here they are if needed
    #ser_topd = ser_topd.rolling(window = 3, center=True).mean()
    #ser_botd = ser_botd.rolling(window = 3, center=True).mean()

    # Removing all the negative valuse to zero
    ser_topd.clip(0, 1, inplace=True)
    ser_botd.clip(0, 1, inplace=True)
    
    # If the value in one is zero, then seting it in zero in both datasets (this removes a lot of noise)
    for i in range(len(ser_topd)):
        if ser_topd.iloc[i] == 0:
            ser_botd.iloc[i] = 0
    for i in range(len(ser_botd)):
        if ser_botd.iloc[i] == 0:
            ser_topd.iloc[i] = 0
    
    # Now we add the two things together and divide by two - i.e we have an average origin activity
    ser_cumd = (ser_topd + ser_botd)/2     

    # Now we want to calculate the quantile of all the positive data.
    templist = np.array([])
    for i in range(1,len(ser_cumd)):
        if ser_cumd.iloc[i] != 0: 
            templist = np.append(templist, ser_cumd.iloc[i])
    # set a cutoff threshold at of the top 10% of all positive values
    cutoff = np.quantile(templist, 0.9)

    # now if a single +ve value (i.e at least one zero either side) is below this cutoff, then set to zero. This again removes noise.
    for i in range(1, len(ser_cumd)-1):
        if ser_cumd.iloc[i] != 0:
            if ser_cumd.iloc[i-1] == 0 and ser_cumd.iloc[i+1] == 0 and ser_cumd.iloc[i]< cutoff:
                ser_cumd.iloc[i] = 0

    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser_cumd

    # Next we want to position "zonal" origins to a single point (bin) and give them a cumulative value (efficiency of the zone)
    start=0
    stop=0
    origins = ser_cumd.clip(0, 0) # creates a new pd.series of zero values same length of ser_cumd
    for i in range(len(ser_cumd)):
        if i <= stop: # simply prevents itterating over the non-zero block
            continue # continue goes back to the for loop
        else:
            if ser_cumd.iloc[i] != 0:
                start = i        
                for z in range(start+1, len(ser_cumd)): # find the ned of the non-zero block
                    if ser_cumd.iloc[z] == 0:
                        stop = z
                        tot = 0
                        tem = 0
                        tem_loc = 0
                        for v in range(i,z): # adds the non-zero block numbers together and identifies the location of the highest value
                            tot = tot + ser_cumd.iloc[v]
                            if ser_cumd.iloc[v] > tem:
                                tem = ser_cumd.iloc[v]
                                tem_loc = v
                        origins.iloc[tem_loc] = tot # adds the total of the non-zero block to the bin with the highest individula differential
                        break

    all_data['origin'] = origins

# create three separate dataframes for the three chromosomes and return these too
 
    chrI = all_data.loc[all_data['chro']=='chr1']
    chrII = all_data.loc[all_data['chro']=='chr2']
    chrIII = all_data.loc[all_data['chro']=='chr3']
    return chrI, chrII, chrIII



#wtchrI, wtchrII, wtchrIII = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)



#    mat_loc = all_data.loc[all_data['chro']=='mating_type_region']
#    mito = all_data.loc[all_data['chro']=='mitochondrial']
#    tel_gap = all_data.loc[all_data['chro']=='chr_II_telomeric_gap']
    
# Note: the %s sign allows us to asiign the variable outputname within a string for multiple file names
#    all_data.to_pickle('./%s.pkl' %(outputname))
#    chrI.to_pickle('./%s_chrI.pkl' %(outputname))
#    chrII.to_pickle('./%s_chrII.pkl' %(outputname))
#    chrIII.to_pickle('./%s_chrIII.pkl' %(outputname))
#    mat_loc.to_pickle('./%s_mat_locus.pkl' %(outputname))
#    mito.to_pickle('./%s_mito.pkl' %(outputname))
#    tel_gap.to_pickle('./%s_tel_gap.pkl' %(outputname))
    

owtchrI, owtchrII, owtchrIII = Create_df('Jo_wt_YE.e1.f-w300.count.csv', 'Jo_wt_YE.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)
#sen1chrI, sen1chrII, sen1chrIII 
osen1chr1, osen1chr2, osen1chr3 = Create_df ('Jo_sen1d_YE.e1.f-w300.count.csv', 'Jo_sen1d_YE.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 5)
odbl8chr1, odbl8chr2, odbl8chr3 = Create_df('Jo_dbl8_YE.e1.f-w300.count.csv', 'Jo_dbl8_YE.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 5)
odschr1, odschr2, odschr3 = Create_df('Jo_sen1_dbl8_YE.e1.f-w300.count.csv', 'Jo_sen1_dbl8_YE.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 5)

#%%



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


con1 = ccontrol[(ccontrol['chro'] == 'chr1')]
con2 = ccontrol[(ccontrol['chro'] == 'chr2')]
con3 = ccontrol[(ccontrol['chro'] == 'chr3')]


def Find(file):
    genes = pd.read_csv(file, delimiter=",")
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


newcontrolfor, newcontrolrev, newccontrol = Find('new_control.csv')

con1 = newccontrol[(newccontrol['chro'] == 'chr1')]
con2 = newccontrol[(newccontrol['chro'] == 'chr2')]
con3 = newccontrol[(newccontrol['chro'] == 'chr3')]


sen1stall = ggenes[(ggenes['genotype'] == 'sen1D')]
dbl8stall = ggenes[(ggenes['genotype'] == 'dbl8D')]
doublestall = ggenes[(ggenes['genotype'] == 'sen1dbl8DD_unique')]


sen1stalls = sen1stall['ID'].to_list()
dbl8stalls = dbl8stall['ID'].to_list()
dsstalls = doublestall['ID'].to_list()

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



gene1 = ggenes[(ggenes['chro'] == 'chr1')]
gene2 = ggenes[(ggenes['chro'] == 'chr2')]
gene3 = ggenes[(ggenes['chro'] == 'chr3')]



aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'Chromosome': [1,2,3]}
comcentro = pd.DataFrame(data=aaaah)
comcentro['Sbinpos'] = comcentro['start']/50
comcentro['Sbinpos'] = comcentro['Sbinpos'].astype(int)
comcentro['Sbinpos'] = comcentro['Sbinpos']*50 +25
comcentro['Ebinpos'] = comcentro['end']/50
comcentro['Ebinpos'] = comcentro['Ebinpos'].astype(int)
comcentro['Ebinpos'] = comcentro['Ebinpos']*50 +25


d = {'start': [3753687], 'end': [3789421]}
centro1 = pd.DataFrame(data=d)

dd = {'start': [1602264], 'end': [1644747]}
centro2 = pd.DataFrame(data=dd)

ddd = {'start': [1070904], 'end': [1137003]}
centro3 = pd.DataFrame(data=ddd)


te1 = {'start': [1, 5554844], 'end': [29663,5579133]}
telo1 = pd.DataFrame(data=te1)
te2 ={'start': [1, 4500619], 'end': [39186,4539800]}
telo2 = pd.DataFrame(data=te2)
te3 ={'start': [], 'end': []}
telo3 = pd.DataFrame(data=te3)

def Chromosome_plot (data1, data2, data3, data4, featurex, centro, telo, genee, c, con,data5, data6, datat):
    ff, (ax1, ax2, ax3, ax5) = plt.subplots(4,1, sharex=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['right_forks'], color ='black', alpha=0.8)
    ax1.plot(data2['pos'], data2['right_forks'], color ='steelblue', alpha=0.8)
    ax1.plot(data3['pos'], data3['right_forks'], color ='orange', alpha=0.8)
    ax1.plot(data4['pos'], data4['right_forks'], color ='tomato', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('rightward forks (e)')
    
    #ax2.set_title('sen1')
    #ax2.plot(data1['pos'], data1['differentials'], color ='saddlebrown', alpha=0.8)
    ax2.plot(data2['pos'], data2['smoo_differentials'], color ='steelblue', alpha=0.8)
    ax2.plot(data3['pos'], data3['smoo_differentials'], color ='orange', alpha=0.8)
    ax2.plot(data4['pos'], data4['smoo_differentials'], color ='tomato', alpha=0.8)
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('differential values (e)')
   # ax2.set_ylim(0,200)

    ax3.set_ylabel('filter dif') 
    
    #ax3.set_title('dbl8')
    #ax3.plot(data3['pos'], data3['ratio_epsilon_f'], color ='teal', alpha=0.8)
   # ax3.plot(data3['pos'], data3['ratio_delta_f'], color ='tomato', alpha=0.8)
    #ax3.set_ylim(-5,5)
    #ax3.set_ylabel('Dbl8∆')
    #ax3.set_ylim(0,200)
    ax5.set_xlabel('Chromosome position')
#    ax4.plot(cc['Pos'], cc['Wnorm'], color ='firebrick', alpha=0.8)
 #   ax4.plot(cc['Pos'], cc['Cnorm'], color ='steelblue', alpha=0.8)
  #  ax4.set_ylim(-10,24)
   # ax4.set_ylabel('WT (HpM)')
    #ax4.plot(data4['pos'], data4['ratio_epsilon_f'], color ='teal', alpha=0.8)
    #ax4.plot(data4['pos'], data4['ratio_delta_f'], color ='tomato', alpha=0.8)
    #ax4.set_ylabel('Sen1∆Dbl8∆')
    #ax4.set_ylim(0,200)
  #  for be in odata1.itertuples(index=False, name = None):
   #     if be[23] > 0.2:
    #        ax1.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
     #       ax2.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            
                  
    for fe in featurex.itertuples(index=False, name=None):
      #  ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                
            #    ax5.annotate(fe[0], xy = [fe[2],-0.15])    
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                ax5.set_ylabel('Gene annotations')
                
    for c in centro.itertuples(index=False, name=None):
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax5.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax5.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                
    for sug in data5.itertuples(index=False, name=None):
        ax3.axvspan(sug[4],sug[4],0.1,0.9,color="blue",alpha=0.3)
        
    for dug in data6.itertuples(index=False, name=None):
        ax3.axvspan(dug[4],dug[4],0.1,0.9,color="orange",alpha=0.3)                

    for dsug in datat.itertuples(index=False, name=None):
        ax3.axvspan(dsug[4],dsug[4],0.1,0.9,color="red",alpha=0.3)
                
  #  for co in con.itertuples(index=False, name=None):
        
   #     if co[4] == 'forward':
    #        ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
     #   elif co[4] == 'reverse':
      #      ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(ewtchrI, esen1chr1, edbl8chr1, edschr1, gene1, centro1, telo1, feat1, 'I', con1, sensen1, dbldbl1, dsds1)
chromosome2 = Chromosome_plot(ewtchrII, esen1chr2, edbl8chr2, edschr2, gene2, centro2, telo2, feat2, 'II', con2, sensen2, dbldbl2, dsds2)
chromosome2 = Chromosome_plot(ewtchrIII, esen1chr3, edbl8chr3, edschr3, gene3, centro3, telo3, feat3, 'III', con3, sensen3, dbldbl3, dsds3)

   
    ax3.bar(data5['pos'], data5['peak_area'], color ='steelblue')
    ax3.bar(data6['pos'], data6['peak_area'], color ='orange')
    ax3.bar(datat['pos'], datat['peak_area'], color ='tomato')







####what i need to add to this plot is the origin positions 
def Chromosome_plot (data1, data2, data3, data4, featurex, centro, telo, genee, c, con,data5, data6, datat):
    ff, (ax1, ax2, ax3, ax5) = plt.subplots(4,1, sharex=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['smoo_right_forks'], color ='black', alpha=0.8)
    ax1.plot(data2['pos'], data2['smoo_right_forks'], color ='steelblue', alpha=0.8)
    ax1.plot(data3['pos'], data3['smoo_right_forks'], color ='orange', alpha=0.8)
    ax1.plot(data4['pos'], data4['smoo_right_forks'], color ='tomato', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('rightward forks (e)')
    
    #ax2.set_title('sen1')
    #ax2.plot(data1['pos'], data1['differentials'], color ='saddlebrown', alpha=0.8)
    ax2.plot(data2['pos'], data2['differentials'], color ='steelblue', alpha=0.8)
    ax2.plot(data3['pos'], data3['differentials'], color ='orange', alpha=0.8)
    ax2.plot(data4['pos'], data4['differentials'], color ='tomato', alpha=0.8)
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('differential values (e)')
   # ax2.set_ylim(0,200)

    ax3.set_ylabel('filter dif') 
    

    ax5.set_xlabel('Chromosome position')

  #  for be in odata1.itertuples(index=False, name = None):
   #     if be[23] > 0.2:
    #        ax1.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
     #       ax3.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
          # ax2.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
 
                  
    for fe in featurex.itertuples(index=False, name=None):
      #  ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                
            #    ax5.annotate(fe[0], xy = [fe[2],-0.15])    
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                ax5.set_ylabel('Gene annotations')
                
    for c in centro.itertuples(index=False, name=None):
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax5.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax5.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
                ax5.annotate(ge[0], xy = [ge[4],0.15])
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                ax5.annotate(ge[0], xy = [ge[4],0.1])
             
                
    for sug in data5.itertuples(index=False, name=None):
        ax3.axvspan(sug[4],sug[5],0.1,0.9,color="steelblue",alpha=0.9)
        
    for dug in data6.itertuples(index=False, name=None):
        ax3.axvspan(dug[4],dug[5],0.1,0.9,color="orange",alpha=0.9)                

    for dsug in datat.itertuples(index=False, name=None):
        ax3.axvspan(dsug[4],dsug[5],0.1,0.9,color="red",alpha=0.9)
                
    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(ewtchrI, esen1chr1, edbl8chr1, edschr1, gene1, centro1, telo1, feat1, 'I', con1, sensen1, dbldbl1, dsds1)
chromosome2 = Chromosome_plot(ewtchrII, esen1chr2, edbl8chr2, edschr2, gene2, centro2, telo2, feat2, 'II', con2, sensen2, dbldbl2, dsds2)
chromosome3 = Chromosome_plot(ewtchrIII, esen1chr3, edbl8chr3, edschr3, gene3, centro3, telo3, feat3, 'III', con3, sensen3, dbldbl3, dsds3)


#%%
###
#now let's asign genes to each peak 

dtcp1 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'I')]
dtcp2 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'II')]
dtcp3 = data_termination_change_peaks_filtered[(data_termination_change_peaks_filtered['chro'] == 'III')]



#this sort of worked yeah 
#i need to do it all manually !!! FUCK ! 
#export the dataframe dcp1 to excel and then go through everything one by one 

#I would need to specifically enter 
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

dtcp1 = assign_id_to_deltastall_peak(dtcp1, feat1)
dtcp2 = assign_id_to_deltastall_peak(dtcp2, feat2)
dtcp3 = assign_id_to_deltastall_peak(dtcp3, feat3)


#dtcp1.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/5jandelta_stalls_chr1.csv', index=False)
#dtcp2.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/5jandelta_stalls_chr2.csv', index=False)
#dtcp3.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/5jandelta_stalls_chr3.csv', index=False)

####so i then manually checked and added missing IDs, so import again now

def Find_pos(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)


    return genes

dtcp1yay = Find_pos("5jandelta_stalls_chr1.csv")


dtcp2yay = Find_pos("5jandelta_stalls_chr2.csv")


dtcp3yay = Find_pos("5jandelta_stalls_chr3.csv")



# Assuming dtcp1yay and feat1 are your dataframes
columns_to_add = ['Start position', 'End position', 'Strand']

# Merge only the specified columns from feat1
dtcp1yay = pd.merge(dtcp1yay, feat1[columns_to_add + ['Systematic ID']], left_on='Systematic ID', right_on='Systematic ID', how='left')
dtcp2yay = pd.merge(dtcp2yay, feat2[columns_to_add + ['Systematic ID']], left_on='Systematic ID', right_on='Systematic ID', how='left')
dtcp3yay = pd.merge(dtcp3yay, feat3[columns_to_add + ['Systematic ID']], left_on='Systematic ID', right_on='Systematic ID', how='left')




dtcp1yay['peak_midpoint'] = (dtcp1yay['start_pos'] + dtcp1yay['end_pos'])/2
dtcp2yay['peak_midpoint'] = (dtcp2yay['start_pos'] + dtcp2yay['end_pos'])/2
dtcp3yay['peak_midpoint'] = (dtcp3yay['start_pos'] + dtcp3yay['end_pos'])/2



def score(file, frame1, frame2, frame3):
    for i in range(len(file)):
        midpoint = file.iloc[i]["peak_midpoint"]
        start = file.iloc[i]["start_pos"]
        stop = file.iloc[i]["end_pos"]
     #   print(premid)
        genotype = file.iloc[i]["genotype"]
        
        if genotype == 'sen1D':
           direction = frame1.loc[(frame1['pos'] >= start) & (frame1['pos'] <= stop)]
           
           file.loc[file.index[i], 'average_forks'] = direction['smoo_right_forks'].sum() /len(direction)

        if genotype == 'dbl8D':
           direction2 = frame2.loc[(frame2['pos'] >= start) & (frame2['pos'] <= stop)]
           file.loc[file.index[i], 'average_forks'] = direction2['smoo_right_forks'].sum() /len(direction2)
                                  

        if genotype == 'sen1dbl8DD':
           direction3 = frame3.loc[(frame3['pos'] >= start) & (frame3['pos'] <= stop)]
           file.loc[file.index[i], 'average_forks'] = direction3['smoo_right_forks'].sum() /len(direction3)



        for i in range(len(file)):
            forkyboi = file.iloc[i]["average_forks"]
          #  if -0.1 < forkyboi < 0.1: 
           #     file.loc[file.index[i], 'stalled_fork'] = 'equal zone'
                
            if forkyboi > 0 :
                   file.loc[file.index[i], 'stalled_fork'] = 'rightward'
            elif forkyboi < 0 :
                   file.loc[file.index[i], 'stalled_fork'] = 'leftward'

                   


    return file

dtcp1yay = score(dtcp1yay, esen1chr1, edbl8chr1, edschr1)
dtcp2yay = score(dtcp2yay, esen1chr2, edbl8chr2, edschr2)
dtcp3yay = score(dtcp3yay, esen1chr3, edbl8chr3, edschr3)

#%%

#i think for the venn diagram i need to do that but to recreate rob's stall diagram
#we can just use 5kb either side of the start_pos 
#so basically adjust for bins being 300bp mid point
#and then use the tranacription pipeline 

dtcp1yay_a = dtcp1yay.dropna(subset=['condition'])

def asign_orientation(file):
    for i in range(len(file)):
        fork = file.iloc[i]["stalled_fork"]
        gene = file.iloc[i]["Strand"]    
        print(fork)
        print(gene)
    
        if fork == 'rightward' and gene == 'forward':
            file.loc[file.index[i], 'orientation'] = 'Head-Tail'
        elif fork == 'rightward' and gene == 'reverse':
            file.loc[file.index[i], 'orientation'] = 'Head-Head'
    
        if fork == "leftward" and gene == 'forward':
            file.loc[file.index[i], 'orientation'] = 'Head-Head'
        elif fork == "leftward" and gene == 'reverse':
            file.loc[file.index[i], 'orientation'] = 'Head-Tail'                 
        
    return file

dtcp1yay_a = asign_orientation(dtcp1yay_a)
dtcp2yay_a = asign_orientation(dtcp2yay)
dtcp3yay_a = asign_orientation(dtcp3yay)


#%%

def Create_df(df,dr,ef,er,w):
    df = pd.read_csv(df) # Reads a csv file

    df['chro'] = df['chro'].replace({"I": "chr1", "II": "chr2", "III": "chr3", 'mitochondrial': 'mito'})
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file

    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef = pd.read_csv(ef, usecols=[2])
    ef['count'].replace(0,1, inplace = True)   
    ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
    er = pd.read_csv(er, usecols=[2])
    er['count'].replace(0,1, inplace = True)
    er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    all_data = pd.concat([df, dr, ef, er], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes
    
    

    

  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
    all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
    all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
    all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
    all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])


# Now we have column for pol delta useage for the duplex
    all_data['d_usage'] = (all_data['ratio_delta_f'] + all_data['ratio_delta_r']) / 2

# now we a column for the percentage of right-moving forks
   # all_data['right_forks']  = all_data['ratio_epsilon_f']*2 - 1
    all_data['right_forks'] = (all_data['norm_ef'] - all_data['norm_er'])/ (all_data['norm_ef'] + all_data['norm_er'])
    
    
# now we will produce a new colum for each sliding window average for each of the calculated columns
# Note: centre = True, means that the data is summed from both left and right. False means its the sum of the last of the number of values.
    
    all_data['smoo_ratio_d_f'] = all_data['ratio_delta_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_d_r'] = all_data['ratio_delta_r'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_f'] = all_data['ratio_epsilon_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_r'] = all_data['ratio_epsilon_r'].rolling(window = w, center=True).mean()
    all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = w, center=True).mean()
    all_data['smoo_right_forks'] = all_data['right_forks'].rolling(window = w, center=True).mean()
    

# now create a differential collunm and a discreate orgin (centre = 1 bin) column.
# Note the script removes differentials not present in both strands and which are singular 
    
    # take unsmoothed pol usage data into two pandas arrays
    ser_etop = all_data['ratio_epsilon_f']
    ser_ebot = all_data['ratio_epsilon_r']

    # Calculate the differentials 
    ser_topd = ser_etop.diff()
    ser_botd = ser_ebot.diff()
    # Reverse the sign on the bottom strand differentials to match top strand
    ser_botd = ser_botd*-1

    # curently dont use a roling average, but here they are if needed
    #ser_topd = ser_topd.rolling(window = 3, center=True).mean()
    #ser_botd = ser_botd.rolling(window = 3, center=True).mean()

    # Removing all the negative valuse to zero
    ser_topd.clip(0, 1, inplace=True)
    ser_botd.clip(0, 1, inplace=True)
    
    # If the value in one is zero, then seting it in zero in both datasets (this removes a lot of noise)
    for i in range(len(ser_topd)):
        if ser_topd.iloc[i] == 0:
            ser_botd.iloc[i] = 0
    for i in range(len(ser_botd)):
        if ser_botd.iloc[i] == 0:
            ser_topd.iloc[i] = 0
    
    # Now we add the two things together and divide by two - i.e we have an average origin activity
    ser_cumd = (ser_topd + ser_botd)/2     

    # Now we want to calculate the quantile of all the positive data.
    templist = np.array([])
    for i in range(1,len(ser_cumd)):
        if ser_cumd.iloc[i] != 0: 
            templist = np.append(templist, ser_cumd.iloc[i])
    # set a cutoff threshold at of the top 10% of all positive values
    cutoff = np.quantile(templist, 0.9)

    # now if a single +ve value (i.e at least one zero either side) is below this cutoff, then set to zero. This again removes noise.
    for i in range(1, len(ser_cumd)-1):
        if ser_cumd.iloc[i] != 0:
            if ser_cumd.iloc[i-1] == 0 and ser_cumd.iloc[i+1] == 0 and ser_cumd.iloc[i]< cutoff:
                ser_cumd.iloc[i] = 0

    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser_cumd

    # Next we want to position "zonal" origins to a single point (bin) and give them a cumulative value (efficiency of the zone)
    start=0
    stop=0
    origins = ser_cumd.clip(0, 0) # creates a new pd.series of zero values same length of ser_cumd
    for i in range(len(ser_cumd)):
        if i <= stop: # simply prevents itterating over the non-zero block
            continue # continue goes back to the for loop
        else:
            if ser_cumd.iloc[i] != 0:
                start = i        
                for z in range(start+1, len(ser_cumd)): # find the ned of the non-zero block
                    if ser_cumd.iloc[z] == 0:
                        stop = z
                        tot = 0
                        tem = 0
                        tem_loc = 0
                        for v in range(i,z): # adds the non-zero block numbers together and identifies the location of the highest value
                            tot = tot + ser_cumd.iloc[v]
                            if ser_cumd.iloc[v] > tem:
                                tem = ser_cumd.iloc[v]
                                tem_loc = v
                        origins.iloc[tem_loc] = tot # adds the total of the non-zero block to the bin with the highest individula differential
                        break

    all_data['origin'] = origins

# create three separate dataframes for the three chromosomes and return these too
 
    chrI = all_data.loc[all_data['chro']=='chr1']
    chrII = all_data.loc[all_data['chro']=='chr2']
    chrIII = all_data.loc[all_data['chro']=='chr3']
    return chrI, chrII, chrIII, all_data



#wtchrI, wtchrII, wtchrIII = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)



#    mat_loc = all_data.loc[all_data['chro']=='mating_type_region']
#    mito = all_data.loc[all_data['chro']=='mitochondrial']
#    tel_gap = all_data.loc[all_data['chro']=='chr_II_telomeric_gap']
    
# Note: the %s sign allows us to asiign the variable outputname within a string for multiple file names
#    all_data.to_pickle('./%s.pkl' %(outputname))
#    chrI.to_pickle('./%s_chrI.pkl' %(outputname))
#    chrII.to_pickle('./%s_chrII.pkl' %(outputname))
#    chrIII.to_pickle('./%s_chrIII.pkl' %(outputname))
#    mat_loc.to_pickle('./%s_mat_locus.pkl' %(outputname))
#    mito.to_pickle('./%s_mito.pkl' %(outputname))
#    tel_gap.to_pickle('./%s_tel_gap.pkl' %(outputname))
    

owtchrI, owtchrII, owtchrIII, owt = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)
#sen1chrI, sen1chrII, sen1chrIII 
osen1chr1, osen1chr2, osen1chr3 = Create_df ('RZ261-RhArepr-d.e1.f-w300.count.csv', 'RZ261-RhArepr-d.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 5)
odbl8chr1, odbl8chr2, odbl8chr3 = Create_df('RZ263-RhArepr-d.e1.f-w300.count.csv', 'RZ263-RhArepr-d.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 5)
odschr1, odschr2, odschr3 = Create_df('RZ265-RhArepr-d.e1.f-w300.count.csv', 'RZ265-RhArepr-d.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 5)


owt1 = owt[(owt['chro'] == 'chr1')]
owt2 = owt[(owt['chro'] == 'chr2')]
owt3 = owt[(owt['chro'] == 'chr3')]

owt1n = owt1[owt1['origin'] > 0.229]
owt2n = owt2[owt2['origin'] > 0.229]
owt3n = owt3[owt3['origin'] > 0.229]



def create_replicon_dataframe(df):
    # Sort the DataFrame by the 'pos' column
    df = df.sort_values(by='pos')

    # Initialize lists to store replicon information
    start_positions = []
    stop_positions = []
    replicon_ids = []

    # Iterate through the sorted DataFrame to find replicons
    for i in range(len(df) - 1):
        start_positions.append(df.iloc[i]['pos'])
        stop_positions.append(df.iloc[i + 1]['pos'])
        replicon_ids.append(i + 1)  # ID starts from 1

    # Create a new DataFrame from the replicon information
    replicon_df = pd.DataFrame({
        'start': start_positions,
        'stop': stop_positions,
        'id': replicon_ids
    })

    return replicon_df

repliconchr1 = create_replicon_dataframe(owt1n)
repliconchr2 = create_replicon_dataframe(owt2n)
repliconchr3 = create_replicon_dataframe(owt3n)



#assign replicon id to identified gene from peaks, so which replicon does it occur in 
#this is an updated version from the delta bias peaks

def assign_id_to_delta(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['Start position'] >= replicon_row['start'] and
                gene_row['End position'] <= replicon_row['stop']):
                matching_replicon_ids.append(gene_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return gene_df

dtcp1yay_aa = assign_id_to_delta(repliconchr1, dtcp1yay_a)
dtcp2yay_aa = assign_id_to_delta(repliconchr2, dtcp2yay_a)
dtcp3yay_aa = assign_id_to_delta(repliconchr3, dtcp3yay_a)


dtcp3yay_aa['id'].fillna(55, inplace=True)


def assign_id_to_stalls(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start'] >= replicon_row['start'] and
                gene_row['end'] <= replicon_row['stop']):
                matching_replicon_ids.append(replicon_row['id'])
                gene_list.append(gene_row['ID'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'id': matching_replicon_ids, 'ID': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='ID', how='left')
    print(id_and_ID)
    
    

   # gene_df = gene_df.append({'id':matching_replicon_ids}, ignore_index=True)
  #  matching_replicon_ids = pd.DataFrame(matching_replicon_ids)
   # print(matching_replicon_ids)
    #gene_df = gene_df.concat({'id':matching_replicon_ids}, ignore_index=True)
    #filtered_replicon_df = replicon_df[replicon_df['id'].isin(matching_replicon_ids)]

    return gene_df

gene1 = assign_id_to_stalls(repliconchr1, gene1)
gene2 = assign_id_to_stalls(repliconchr2, gene2)
gene3 = assign_id_to_stalls(repliconchr3, gene3)
    
#%%
####now import your venn diagram

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3 as venn3
from matplotlib_venn import venn3_circles as venn3_circles
from matplotlib import pyplot as plt


venn3(subsets = (0, 2, 0, 108, 7, 26, 2), set_labels = ('sen1∆', 'dbl8∆', 'sen1∆dbl8∆'), set_colors=('b', 'g', 'r'), alpha = 0.3);
plt.savefig('foo.png', bbox_inches='tight')



#%%

###this is my code copied from untitled 2 (rnap2 data script)
def getmenamesofthebastards(file, col, x):
    threshold = 2
    subset = file[(file[col]) > threshold]
    subset = subset.sort_values(x, ascending = False)
    listi = subset['ID'].to_list()
    return subset, listi

senlist ,senover2list = getmenamesofthebastards(joinedd, 'sen1Z', 'sen1ratio')
dbl8list ,dbl8over2list = getmenamesofthebastards(joinedddbl, 'dbl8Z', 'dbl8ratio')
dslist ,dsover2list = getmenamesofthebastards(joineddds, 'dsZ', 'dsratio')

dsdlist, dsdover2list = getmenamesofthebastards(joinedc, 'dsdZ', 'dsdratio')
sddlist, sddover2list = getmenamesofthebastards(joineds, 'sddZ', 'sddratio')
bslist, bsover2list = getmenamesofthebastards(joinedbs, 'bsZ', 'bsratio')


def getmenamesofthebastardsbelow(file, col, x):
    threshold = -2
    subset = file[(file[col]) < threshold]
    subset = subset.sort_values(x, ascending = False)
    listi = subset['ID'].to_list()
    return subset, listi


senlistlow ,senover2listlow = getmenamesofthebastardsbelow(joinedd, 'sen1Z', 'sen1ratio')
dbl8listlow ,dbl8over2listlow = getmenamesofthebastardsbelow(joinedddbl, 'dbl8Z', 'dbl8ratio')
dslistlow ,dsover2listlow = getmenamesofthebastardsbelow(joineddds, 'dsZ', 'dsratio')


dsdlistlow, dsdover2listlow = getmenamesofthebastardsbelow(joinedc, 'dsdZ', 'dsdratio')
sddlistlow, sddover2listlow = getmenamesofthebastardsbelow(joineds, 'sddZ', 'sddratio')
bslistlow, bsover2listlow = getmenamesofthebastardsbelow(joinedbs, 'bsZ', 'bsratio')

##think about plotting differences or similarites between dbl8 and ds 

SEN = senover2list.extend(senover2listlow)
DBL = dbl8over2list.extend(dbl8over2listlow)
DS = dsover2list.extend(dsover2listlow)



control = pd.read_csv('control_genes.txt', delimiter="\t") 
controlname = control['ID'].to_list()

#%%
#this next bit checks for overlaps between delta stalls and delta bias.
def Find_fix(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)


    return genes

dweltause1 = Find_fix("delta_usage_chr1.csv")


dweltause2 = Find_fix("delta_usage_chr2.csv")


dweltause3 = Find_fix("delta_usage_chr3.csv")



'''

#pull replicons with stalls
def filter_replicons_by_genes(replicon_df, gene_df):
    matching_replicon_ids = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start'] >= replicon_row['start'] and
                gene_row['end'] <= replicon_row['stop']):
                matching_replicon_ids.append(replicon_row['id'])

    filtered_replicon_df = replicon_df[replicon_df['id'].isin(matching_replicon_ids)]

    return filtered_replicon_df

repliconchr1 = filter_replicons_by_genes(repliconchr1, gene1)
repliconchr2 = filter_replicons_by_genes(repliconchr2, gene2)
repliconchr3 = filter_replicons_by_genes(repliconchr3, gene3)

#what i need is pull replicon with stall 
# give also deltas, and make new col in delta
#

        temp = pd.DataFrame({'ID': [key],
                             'Status': ['YES' if interest.loc['G'] == interest.max() or interest.loc['C'] == interest.max() else 'NO'],
                             'Rank': [interest.max()]})

        gene_status_df = pd.concat([gene_status_df, temp])



def building_layers(replicon_df, gene_df, delta):
    matching_replicon_ids = []
    
    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            temp = pd.DataFrame({'ID': replicon_row['id'],
                                 'stall_containing': ['YES' if gene_row['start'] >= replicon_row['start'] and gene_row['end'] <= replicon_row['stop'] else 'NO']})
            
            replicon_df = pd.concat([replicon_df, temp])
    return replicon_df

repliconchr1 = building_layers(repliconchr1, gene1, deltachr1)
 '''             
    

#OpenAI. (2022). GPT-3.5 (ChatGPT 3.5). Retrieved from https://www.openai.com/
def building_layers(replicon_df, gene_df, delta):
    
    matching_replicon_ids = []
    mathcing_rep_and_d_id = []
    topdeltapeaks = []
  #  matching_gene_ids = []  # New list to store gene IDs
    
    # Iterate over rows in gene_df
    #this line of code asks if a stall gene occurs within a given replicon and iterates, 
    #creates boolean mask
    for _, gene_row in gene_df.iterrows():
        # Check the condition using vectorized operations
        mask = (replicon_df['start'] <= gene_row['start']) & (gene_row['end'] <= replicon_df['stop'])
   #     print(mask)
        
        # Extract matching replicon IDs and add to the list
        #this takes which replicon_id the stall occurs in 
        
        matching_replicon_ids.extend(replicon_df.loc[mask, 'id'].tolist())
        

    # Create a DataFrame with matching IDs and 'YES' for 'stall_containing'
    temp = pd.DataFrame({'id': matching_replicon_ids, 'stall_containing': 'YES'})
    
    
    # Merge the temporary DataFrame with replicon_df
    replicon_df = replicon_df.merge(temp, on='id', how='left')

    # Fill NaN values in 'stall_containing' with 'NO'
    replicon_df['stall_containing'] = replicon_df['stall_containing'].fillna('NO')
    
       
    delta['peak_no'] = pd.to_numeric(delta['peak_no'], errors='coerce')
    for _, delta_row in delta.iterrows():
        
        #so now here this asks if a 2% top delta peak occurs in a replicon 
        mask2 = (replicon_df['start'] <= delta_row['start_pos']) & (delta_row['end_pos'] <= replicon_df['stop'])
       # print(f"Mask2: {mask2}")
      #  print(f"Delta columns: {delta.columns}")
      #  if any(mask2):
       #     print("At least one True value in mask2")
        #takes the replicon id where they is a delta peak 
        mathcing_rep_and_d_id.extend(replicon_df.loc[mask2, 'id'].tolist())
      #  topdeltapeaks.extend(delta.loc[mask2, 'peak_no'].tolist())
       # topdeltapeaks.append(delta_row['peak_no'])
       # print(topdeltapeaks)
        
  #  temp2 = pd.DataFrame({'id': mathcing_rep_and_d_id, 'delta_containing': 'YES', 'delta_peak_no' : topdeltapeaks})
    temp2 = pd.DataFrame({'id': mathcing_rep_and_d_id, 'delta_containing': 'YES'})
    replicon_df = replicon_df.merge(temp2, on='id', how='left')
    replicon_df['delta_containing'] = replicon_df['delta_containing'].fillna('NO')
    
    
   
    #what i need to do now is add another later, whereby if mask = true, i also add peak_no to list
    
    df_no_duplicates = replicon_df.drop_duplicates(subset='id', keep='first')

        


    return df_no_duplicates

# Example usage
repliconchr1 = building_layers(repliconchr1, gene1, deltachr1)
repliconchr2 = building_layers(repliconchr2, gene2, deltachr2)
repliconchr3 = building_layers(repliconchr3, gene3, deltachr3)



'''

def assign_peak_no_to_stalls(delta_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, delta_row in delta_df.iterrows():
            if (gene_row['start'] >= delta_row['start_pos'] and
                gene_row['end'] <= delta_row['end_pos']):
                matching_replicon_ids.append(delta_row['peak_no'])
                gene_list.append(gene_row['ID'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'ID': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='ID', how='left')
    print(id_and_ID)
    
    
    return gene_df

gene1 = assign_peak_no_to_stalls(deltachr1, gene1)
gene2 = assign_id_to_stalls(repliconchr2, gene2)
gene3 = assign_id_to_stalls(repliconchr3, gene3)



#assign replicon id to delta peak, so which replicon does it occur in 
def assign_id_to_stalls(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start_pos'] >= replicon_row['start'] and
                gene_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(replicon_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    
    

   # gene_df = gene_df.append({'id':matching_replicon_ids}, ignore_index=True)
  #  matching_replicon_ids = pd.DataFrame(matching_replicon_ids)
   # print(matching_replicon_ids)
    #gene_df = gene_df.concat({'id':matching_replicon_ids}, ignore_index=True)
    #filtered_replicon_df = replicon_df[replicon_df['id'].isin(matching_replicon_ids)]

    return gene_df

gene1 = assign_id_to_stalls(repliconchr1, deltachr1)
gene2 = assign_id_to_stalls(repliconchr2, deltachr2)
gene3 = assign_id_to_stalls(repliconchr3, deltachr3)
'''    


#%%%

gene1 = gene1.merge(repliconchr1, on='id', how='left')
gene2 = gene2.merge(repliconchr2, on='id', how='left')
gene3 = gene3.merge(repliconchr3, on='id', how='left')

gene1 = gene1.merge(deltachr1, on='id', how='left')
assert all(gene1.loc[gene1['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"
gene2 = gene2.merge(deltachr2, on='id', how='left')
assert all(gene2.loc[gene2['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"
gene3 = gene3.merge(deltachr3, on='id', how='left')
assert all(gene3.loc[gene3['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"

allthegoodg = [gene1, gene2, gene3]
allthestalls = pd.concat(allthegoodg)

#%%

deltastallswithrepliconidsblg = [dtcp1yay_aa, dtcp2yay_aa, dtcp3yay_aa]
deltastallswithrepliconids = pd.concat(deltastallswithrepliconidsblg)
deltastallswithrepliconids.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/5thjan_deltastallswithrepliconids.csv', index=False)

dtcp1yay_aa_list = dtcp1yay_aa['id'].to_list()
dtcp2yay_aa_list = dtcp2yay_aa['id'].to_list()
dtcp3yay_aa_list = dtcp3yay_aa['id'].to_list()

gene1_list = gene1['id'].to_list()
gene2_list = gene2['id'].to_list()
gene3_list = gene3['id'].to_list()


dweltause1l = dweltause1['id'].to_list()
dweltause2l = dweltause2['id'].to_list()
dweltause3l = dweltause3['id'].to_list()

deltastallswithrepliconids_uc = gene3['id'].nunique() #90




gene3_list = set(gene3_list)
gene3_list = list(gene3_list)

#gene1 - 52
#gene 2= 46
#gene 3 = 26
#epsilon stalls 124 unique replicons


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

d1vg1 = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2 = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3 = intersection (dtcp3yay_aa_list, gene3_list)


d1vd1 = intersection (dtcp1yay_aa_list, dweltause1l)
d2vd2 = intersection (dtcp2yay_aa_list, dweltause2l)
d3vd3 = intersection (dtcp3yay_aa_list, dweltause3l)


def fucker(lst1, lst2):
    lst3 = [value for value in lst1 if not value in lst2]
    return lst3

d1vg1_NOT = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2_NOT = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3_NOT = intersection (dtcp3yay_aa_list, gene3_list)



#%%

#overlappingdelta_and_epsi = [d1vg1, d2vg2, d3vg3]
# Assuming df is your DataFrame and 'id' is the column name
delta_filter_foroverlap1 = dtcp1yay_aa[dtcp1yay_aa['id'].isin(d1vg1)]
delta_filter_foroverlap2 = dtcp2yay_aa[dtcp2yay_aa['id'].isin(d2vg2)]
delta_filter_foroverlap3 = dtcp3yay_aa[dtcp3yay_aa['id'].isin(d3vg3)]


##now do the same for the epsilon joey
epi_filter_foroverlap1 = gene1[gene1['id'].isin(d1vg1)]
epi_filter_foroverlap2 = gene2[gene2['id'].isin(d2vg2)]
epi_filter_foroverlap3 = gene3[gene3['id'].isin(d3vg3)]

def fuckfuckfuck(fileE, fileD):
    for i in range(len(fileE)):
        repliconE = fileE.iloc[i]["id"]
        startE =fileE.iloc[i]["start"]
        endE =fileE.iloc[i]["end"]
        
        for i in range(len(fileD)):
            repliconD = fileD.iloc[i]["id"]
            startD =fileD.iloc[i]["Start position"]
            endD = fileD.iloc[i]["End position"]
            
                
            if repliconE == repliconD:
                print('yes')
            
    return fileE, fileD
#epi_filter_foroverlap1, delta_filter_foroverlap1 = fuckfuckfuck( epi_filter_foroverlap1, delta_filter_foroverlap1)
        



def classify_positions(df_delta, df_epi):
    # Create a new column in df_delta to store the classification
    df_delta['position_class'] = None
    
    # Iterate over rows in df_delta
    for index, row in df_delta.iterrows():
        # Find matching rows in df_epi based on 'id'
        matching_rows = df_epi[df_epi['id'] == row['id']]
        
        # Check each matching row in df_epi
        for _, epi_row in matching_rows.iterrows():
            # Check conditions and assign classification
            if (row['End position'] < epi_row['start']) and (row['End position'] < epi_row['end']):
                df_delta.at[index, 'position_class'] = 'before'
                break  # Stop checking once a classification is made
            elif (row['Start position'] > epi_row['end']) and (row['Start position'] > epi_row['start']):
                df_delta.at[index, 'position_class'] = 'after'
                break  # Stop checking once a classification is made
            elif row['Systematic ID'] == epi_row['ID']:
                df_delta.at[index, 'position_class'] = 'same'
                break  # Stop checking once a classification is made
    
    return df_delta

# Example usage

delta_filter_foroverlap1 = classify_positions(delta_filter_foroverlap1, epi_filter_foroverlap1)
delta_filter_foroverlap2 = classify_positions(delta_filter_foroverlap2, epi_filter_foroverlap2)
delta_filter_foroverlap3 = classify_positions(delta_filter_foroverlap3, epi_filter_foroverlap3)



def classify_positions(df_delta, df_epi):
    # Create new columns in df_delta to store the classification and epsilon_ID
    df_delta['position_class'] = None
    df_delta['epsilon_ID'] = None
    df_delta['epsilon_coding_strand'] = None
    df_delta['epsilon_stalled_fork'] = None
    df_delta['epsilon_genotype'] = None
    
    # Iterate over rows in df_delta
    for index, row in df_delta.iterrows():
        # Find matching rows in df_epi based on 'id'
        matching_rows = df_epi[df_epi['id'] == row['id']]
        
        # Check each matching row in df_epi
        for _, epi_row in matching_rows.iterrows():
            # Check conditions and assign classification
            if (row['End position'] < epi_row['start']) and (row['End position'] < epi_row['end']):
                df_delta.at[index, 'position_class'] = 'before'
            elif (row['Start position'] > epi_row['end']) and (row['Start position'] > epi_row['start']):
                df_delta.at[index, 'position_class'] = 'after'
            elif row['Systematic ID'] == epi_row['ID']:
                df_delta.at[index, 'position_class'] = 'same'
                df_delta.at[index, 'epsilon_ID'] = epi_row['ID']
                df_delta.at[index, 'epsilon_coding_strand'] = epi_row['coding_strand']
                df_delta.at[index, 'epsilon_stalled_fork'] = epi_row['stalled_fork']
                df_delta.at[index, 'epsilon_genotype'] = epi_row['genotype']
                break  # Stop checking once a classification is made
    
    # Drop 'id' column from df_delta
    df_delta.drop('id', axis=1, inplace=True)
    
    # Rename columns with the prefix "epsilon_"
    df_delta.rename(columns={
        'coding_strand': 'epsilon_coding_strand',
        'stalled_fork': 'epsilon_stalled_fork',
        'genotype': 'epsilon_genotype'
    }, inplace=True)
    
    return df_delta


checkitbitch = classify_positions(delta_filter_foroverlap1, epi_filter_foroverlap1)

def classify_positions(df_delta, df_epi):
    # Create a new column in df_delta to store the classification
    df_delta['position_class'] = None
    df_delta['epsilon_ID'] = None
    df_delta['epsilon_coding_strand'] = None
    df_delta['epsilon_stalled_fork'] = None
    df_delta['epsilon_genotype'] = None
    
    # Iterate over rows in df_delta
    for index, row in df_delta.iterrows():
        # Find matching rows in df_epi based on 'id'
        matching_rows = df_epi[df_epi['id'] == row['id']]
        
        # Check each matching row in df_epi
        for _, epi_row in matching_rows.iterrows():
            # Check conditions and assign classification
            if (row['End position'] < epi_row['start']) and (row['End position'] < epi_row['end']):
                df_delta.at[index, 'position_class'] = 'before'
                df_delta.at[index, 'epsilon_ID'] = epi_row['ID']
                df_delta.at[index, 'epsilon_coding_strand'] = epi_row['coding_strand']
                df_delta.at[index, 'epsilon_stalled_fork'] = epi_row['stalled_fork']
                df_delta.at[index, 'epsilon_genotype'] = epi_row['genotype']
                
                break  # Stop checking once a classification is made
            elif (row['Start position'] > epi_row['end']) and (row['Start position'] > epi_row['start']):
                df_delta.at[index, 'position_class'] = 'after'
                df_delta.at[index, 'epsilon_ID'] = epi_row['ID']
                df_delta.at[index, 'epsilon_coding_strand'] = epi_row['coding_strand']
                df_delta.at[index, 'epsilon_stalled_fork'] = epi_row['stalled_fork']
                df_delta.at[index, 'epsilon_genotype'] = epi_row['genotype']
                
                
                break  # Stop checking once a classification is made
            elif row['Systematic ID'] == epi_row['ID']:
                df_delta.at[index, 'position_class'] = 'same'
                df_delta.at[index, 'epsilon_ID'] = epi_row['ID']
                df_delta.at[index, 'epsilon_coding_strand'] = epi_row['coding_strand']
                df_delta.at[index, 'epsilon_stalled_fork'] = epi_row['stalled_fork']
                df_delta.at[index, 'epsilon_genotype'] = epi_row['genotype']
                
                break  # Stop checking once a classification is made
    
    return df_delta
checkitbitch = classify_positions(delta_filter_foroverlap1, epi_filter_foroverlap1)

delta_filter_foroverlap1 = classify_positions(delta_filter_foroverlap1, epi_filter_foroverlap1)
delta_filter_foroverlap2 = classify_positions(delta_filter_foroverlap2, epi_filter_foroverlap2)
delta_filter_foroverlap3 = classify_positions(delta_filter_foroverlap3, epi_filter_foroverlap3)



ohyeahboo = pd.concat([delta_filter_foroverlap1, delta_filter_foroverlap2, delta_filter_foroverlap3], ignore_index=True)
ohyeahboo.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/5janoverlapsbetweendeltaandepsilon', index=False)


duplicate_counts = ohyeahboo['position_class'].value_counts()

#%%
set1 = {1, 2, 3, 4, 5}
set2 = {3, 4, 5, 6, 7}

# Create Venn diagram
venn_labels = {'100': 'Set 1', '010': 'Set 2', '110': 'Intersection'}
venn_diagram = venn2([set1, set2], set_labels=('Set 1', 'Set 2'))
venn_diagram.get_label_by_id('100').set_text(venn_labels['100'])
venn_diagram.get_label_by_id('010').set_text(venn_labels['010'])
venn_diagram.get_label_by_id('110').set_text(venn_labels['110'])

# Display the plot
plt.show()

#48 in common between e and d 
#145 from d, 147 fom d


from matplotlib_venn import venn2

delta_size = 111
epsilon_size = 124
overlap_size = 39

delta_color = 'steelblue'
epsilon_color = 'tomato'
overlap_color = 'purple'

# Create a Venn diagram with custom colors
venn2(subsets=(delta_size - overlap_size, epsilon_size - overlap_size, overlap_size),
      set_labels=('Delta stall replicons', 'Epsilon stall replicons'),
      set_colors=(delta_color, epsilon_color),
      alpha=0.5,  # Adjust transparency if needed
      normalize_to=1.0)  # Set to 1.0 for proportions, adjust if needed

# Display the plot
plt.show()
#wt



#####please note the next step is to seperate by genotype woooooo

dinerosen1 = dtcp1yay_aa[(dtcp1yay_aa['genotype'] == 'sen1D')]
dinerosen2 = dtcp2yay_aa[(dtcp2yay_aa['genotype'] == 'sen1D')]
dinerosen3 = dtcp3yay_aa[(dtcp3yay_aa['genotype'] == 'sen1D')]

ogsen1 = gene1[(gene1['genotype'] == 'dbl8D')]
ogsen2 = gene2[(gene2['genotype'] == 'dbl8D')]
ogsen3 = gene3[(gene3['genotype'] == 'dbl8D')]

#all to list
dinerosen1l = dinerosen1['id'].to_list()
dinerosen2l = dinerosen2['id'].to_list()
dinerosen3l = dinerosen3['id'].to_list()

ogsen1l = ogsen1['id'].to_list()
ogsen2l = ogsen2['id'].to_list()
ogsen3l = ogsen3['id'].to_list()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

d1vg1 = intersection (dinerosen1l, ogsen1l)
d2vg2 = intersection (dinerosen2l, ogsen2l)
d3vg3 = intersection (dinerosen3l, ogsen3l)
# 2 in common
#117 from the delta
#9 from the epi
len(ogsen1l) + len(ogsen2l) + len(ogsen3l)
len(dinerosen1l) + len(dinerosen2l) + len(dinerosen3l)
len(d1vg1) + len(d2vg2) + len(d3vg3)

delta_size = 117
epsilon_size = 30
overlap_size = 9

delta_color = 'blue'
epsilon_color = 'red'
overlap_color = 'purple'

# Create a Venn diagram with custom colors
venn2(subsets=(delta_size - overlap_size, epsilon_size - overlap_size, overlap_size),
      set_labels=('Delta ssen1D tall replicons', 'Epsilon dbl8D stall replicons'),
      set_colors=(delta_color, epsilon_color),
      alpha=0.3,  # Adjust transparency if needed
      normalize_to=1.0)  # Set to 1.0 for proportions, adjust if needed

# Display the plot
plt.show()

#####now that that's done omg
##just try to do a plot like a violin or swarm to show that the 
#intensity of peaks from epsilon are more than delta
#or maybe just try to adapt the python stall call to call from both delta and epsilon peaks 
#see how many survive that omg 

#%%
epistallscon = [etcp1, etcp2, etcp3]
epistallsall = pd.concat(epistallscon)



#%%


#using pol2 plot 
deltastalls = [dtcp1, dtcp2, dtcp3]
deltastallsall = pd.concat(deltastalls)

   

#what if we do again but with more detail,
# so divide by genotype but also orientation 


deltastalls = [dtcp1yay_a, dtcp2yay_a, dtcp3yay_a]
deltastallsall = pd.concat(deltastalls) 

deltastallswithrepliconids

epistallsall = Find_pos("epistallsall_frompythonpuseqcallsforepsilon.csv")
epistallsall['source'] = 'epsilon'
deltastallsall['source'] = 'delta'


allthestallsnow = [epistallsall, deltastallsall]
compiledepianddelta = pd.concat(allthestallsnow) 


compiledepianddelta['peak_size'] = (compiledepianddelta['end_pos'] - compiledepianddelta['start_pos'])/300

compiledepianddelta['peak_strength'] = (compiledepianddelta['peak_area'] / compiledepianddelta['peak_size'])*compiledepianddelta['peak_max']
compiledepianddelta['log2_peak_strength'] = np.log2(compiledepianddelta['peak_strength'])

compiledepianddelta['log2_peakmax'] = np.log2(compiledepianddelta['peak_max'])
compiledepianddelta['log2_peakarea'] = np.log2(compiledepianddelta['peak_area'])


custom_palette = ["#EF6262", "#1D5B79" ]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

sns.violinplot(
    data=compiledepianddelta, x="source", y="log2_peak_strength",
    hue="source"
)


sns.kdeplot(data=compiledepianddelta, x="log2_peakmax", y="log2_peakarea", hue="source")

sns.jointplot(data=compiledepianddelta, x="log2_peakmax", y="log2_peakarea", hue="source", kind="kde")




sns.kdeplot(
   data=compiledepianddelta, x="log2_peak_strength", hue="source",
   fill=True, common_norm=False,
   alpha=.5, linewidth=0,
)





selected_columns = compiledepianddelta[['peak_area', 'peak_max', 'source']]



# Map colors to unique values in 'source'
lut = dict(zip(compiledepianddelta['source'].unique(), "rbg"))
row_colors = compiledepianddelta['source'].map(lut)

# Create a clustermap
sns.clustermap(compiledepianddelta[['peak_area', 'peak_max']], row_colors=row_colors)


sns.clustermap(
    selected_columns,
    figsize=(7, 5),
    row_cluster=False,
    dendrogram_ratio=(.1, .2),
    cbar_pos=(0, .2, .03, .4)
)


#%%

#please note that this is redundant
epistallsall['log2_peakmax'] = np.log2(epistallsall['peak_max'])
epistallsall['log2_peakarea'] = np.log2(epistallsall['peak_area'])

sns.jointplot(data=epistallsall, x="log2_peakmax", y="log2_peakarea", hue="genotype", kind="kde")



epistallsall['peak_size'] = (epistallsall['end_pos'] - epistallsall['start_pos'])/300

epistallsall['peak_strength'] = (epistallsall['peak_area'] / epistallsall['peak_size'])*epistallsall['peak_max']
epistallsall['log2_peak_strength'] = np.log2(epistallsall['peak_strength'])

epistallsally = epistallsall.copy()
epistallsally["genotype"] = pd.Categorical(epistallsally["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD"])
custom_palette = ["steelblue", "orange", 'tomato' ]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

sns.violinplot(
    data=epistallsally, x="genotype", y="log2_peak_strength",
    hue="genotype", alpha= 0.8
)
plt.ylim(-13, -5)

# Show the plot
plt.show()

sns.lmplot(data=epistallsall, x="log2_peakarea", y="log2_peakmax", hue="genotype")



#####################
#and now for delta 
deltastallsall['log2_peakmax'] = np.log2(deltastallsall['peak_max'])
deltastallsall['log2_peakarea'] = np.log2(deltastallsall['peak_area'])

sns.jointplot(data=deltastallsall, x="log2_peakmax", y="log2_peakarea", hue="genotype", kind="kde")



deltastallsall['peak_size'] = (deltastallsall['end_pos'] - deltastallsall['start_pos'])/300

deltastallsall['peak_strength'] = (deltastallsall['peak_area'] / deltastallsall['peak_size'])*deltastallsall['peak_max']
deltastallsall['log2_peak_strength'] = np.log2(deltastallsall['peak_strength'])


custom_palette = ["steelblue", "orange", 'tomato' ]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

ax.set_ylim=(-5, -13)
sns.violinplot(
    data=deltastallsall, x="genotype", y="log2_peak_strength",
    hue="genotype", alpha= 0.8
)
plt.ylim(-13, -5)

# Show the plot
plt.show()
plt.show()

sns.lmplot(data=deltastallsall, x="log2_peakarea", y="log2_peakmax", hue="genotype")



epistallsall.to_csv('/Users/patricfernandez/Documents/python/1stjanepsilonstalls', index=False)
deltastallsall.to_csv('/Users/patricfernandez/Documents/python/1stjanDELTAstalls.csv', index=False)
#%%

#merged_df = pd.merge(complete_range, df, on=pos_column, how='left')
netho = pd.read_csv("NetSeqbin.csv", header=0, index_col=False)  

def Find_pos(file):

    HHst = file[(file['orientation'] == 'Head-Tail')]
    HTst = file[(file['orientation'] == 'Head-Head')]


    return HHst, HTst

HTstallsall, HHstallsall = Find_pos(deltastallsall)

sdelta = deltastallsall[(deltastallsall['genotype'] == 'sen1D')]
ddelta = deltastallsall[(deltastallsall['genotype'] == 'dbl8D')]
dsdelta = deltastallsall[(deltastallsall['genotype'] == 'sen1dbl8DD')]



#note that since we are interogating position of peak (not genes), 
#don't care if forward or reverse strand)

def RStart(genesfor, chr1, chr2, chr3, p, k):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[7] == k:
            if g[8] == p:
                if g[5] == 'chr1':
                    x = chr1.index[chr1['Pos'] == g[9]].tolist()
                    xx.append(x) 

                if g[5] == 'chr2':
                    x2 = chr2.index[chr2['Pos'] == g[9]].tolist()
                    xx.append(x2)

                if g[5] == 'chr3':
                    x3 = chr3.index[chr3['Pos'] == g[9]].tolist()
                    xx.append(x3)
                
    return xx
dsforxxR = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxR = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxL = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxL = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')


def REnd(genesfor, chr1, chr2, chr3, p, k):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[7] == k:
            if ge[8] == p:
                if ge[5] == 'chr1':
                    xe = chr1.index[chr1['Pos'] == ge[10]].tolist()
                    xxe.append(xe) 

                if ge[5] == 'chr2':
                    xe2 = chr2.index[chr2['Pos'] == ge[10]].tolist()
                    xxe.append(xe2)

                if ge[5] == 'chr3':
                    xe3 = chr3.index[chr3['Pos'] == ge[10]].tolist()
                    xxe.append(xe3)
    return xxe

dsforxxeR = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxeR = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxeL = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxeL = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')

#lets make two list, one is ds hh start 
# rev right and 
# for left 


#for wt , forward coding strand genes only 
#TSS
def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[2] == p:
            if g[4] == 'I':
                x = chr1.index[chr1['pos'] == g[5]].tolist()
                xx.append(x) 

            if g[4] == 'II':
                x2 = chr2.index[chr2['pos'] == g[5]].tolist()
                xx.append(x2)

            if g[4] == 'III':
                x3 = chr3.index[chr3['pos'] == g[5]].tolist()
                xx.append(x3)
                
    return xx

senforxx = Start(deltastallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
#senrevxx = Start(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxx = Start(deltastallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxR = Start(deltastallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')
#dbl8revxx = Start(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxx = Start(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxx = Start(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')

HTstallsall, HHstallsall
###
senforxx = Start(HTstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
senforxxHH = Start(HHstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')


dbl8forxx = Start(HTstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dbl8forxxHH = Start(HHstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')

dsforxxR = Start(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')
dsforxxRHH = Start(HHstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')




#TES
def End(file, chr1, chr2, chr3, p):
    xxe =[]
    for ge in file.itertuples():
        if ge[2] == p:
            if ge[4] == 'I':
                xe = chr1.index[chr1['pos'] == ge[6]].tolist()
                xxe.append(xe) 

            if ge[4] == 'II':
                xe2 = chr2.index[chr2['pos'] == ge[6]].tolist()
                xxe.append(xe2)

            if ge[4] == 'III':
                xe3 = chr3.index[chr3['pos'] == ge[6]].tolist()
                xxe.append(xe3)
    return xxe

senforxxe = End(deltastallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
#senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(deltastallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxeR = End(deltastallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')
#dbl8revxxe = End(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxxe = End(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxxe = End(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


HTstallsall, HHstallsall
###
senforxxe = End(HTstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
senforxxeHH = End(HHstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
#senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(HTstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dbl8forxxeHH = End(HHstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxeR = End(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')
dsforxxeRHH = End(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')

#%%
#specificall y for netho data
HTstallsall['sbinpos'] = HTstallsall['Start position'] -25
HTstallsall['ebinpos'] = HTstallsall['End position'] -25

HHstallsall['sbinpos'] = HHstallsall['Start position'] -25
HHstallsall['ebinpos'] = HHstallsall['End position'] - 25


    HTstallsall['Sbinpos'] = HTstallsall['Start position']/50
    HTstallsall['Sbinpos'] = HTstallsall['Sbinpos'].astype(int)
    HTstallsall['Sbinpos'] = HTstallsall['Sbinpos']*50 + 25
    HTstallsall['Ebinpos'] = HTstallsall['End position']/50
    HTstallsall['Ebinpos'] = HTstallsall['Ebinpos'].astype(int)
    HTstallsall['Ebinpos'] = HTstallsall['Ebinpos']*50 +25
    
    HHstallsall['Sbinpos'] = HHstallsall['Start position']/50
    HHstallsall['Sbinpos'] = HHstallsall['Sbinpos'].astype(int)
    HHstallsall['Sbinpos'] = HHstallsall['Sbinpos']*50 + 25
    HHstallsall['Ebinpos'] = HHstallsall['End position']/50
    HHstallsall['Ebinpos'] = HHstallsall['Ebinpos'].astype(int)
    HHstallsall['Ebinpos'] = HHstallsall['Ebinpos']*50 +25


netho = pd.read_csv("NetSeqbin.csv", header=0, index_col=False)  
netho = netho[(netho['coding_strand'] == 'reverse')] 

net1 = netho[(netho['chro'] == 'chr1')]
net2 = netho[(netho['chro'] == 'chr2')]
net3 = netho[(netho['chro'] == 'chr3')]

HTstallsallF = HTstallsall[(HTstallsall['Strand'] == 'forward')]
HTstallsallR = HTstallsall[(HTstallsall['Strand'] == 'reverse')]

HHstallsallF = HHstallsall[(HHstallsall['Strand'] == 'forward')]
HHstallsallR = HHstallsall[(HHstallsall['Strand'] == 'reverse')]

def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[2] == p:
            if g[4] == 'I':
                x = chr1.index[chr1['pos'] == g[25]].tolist()
                xx.append(x) 

            if g[4] == 'II':
                x2 = chr2.index[chr2['pos'] == g[25]].tolist()
                xx.append(x2)

            if g[4] == 'III':
                x3 = chr3.index[chr3['pos'] == g[25]].tolist()
                xx.append(x3)
                
    return xx



###
Fsenforxxx = Start(HTstallsallF, net1, net2, net3, 'sen1D')
FsenforxxHH = Start(HHstallsallF,  net1, net2, net3, 'sen1D')
Rsenforxxx = Start(HTstallsallR, net1, net2, net3, 'sen1D')
RsenforxxHH = Start(HHstallsallR,  net1, net2, net3, 'sen1D')


Fdbl8forxx = Start(HTstallsallF,  net1, net2, net3, 'dbl8D')
Fdbl8forxxHH = Start(HHstallsallF,  net1, net2, net3, 'dbl8D')
Rdbl8forxx = Start(HTstallsallR,  net1, net2, net3, 'dbl8D')
Rdbl8forxxHH = Start(HHstallsallR,  net1, net2, net3, 'dbl8D')


FdsforxxR = Start(HTstallsallF,  net1, net2, net3, 'sen1dbl8DD')
FdsforxxRHH = Start(HHstallsallF,  net1, net2, net3, 'sen1dbl8DD')
RdsforxxR = Start(HTstallsallR,  net1, net2, net3, 'sen1dbl8DD')
RdsforxxRHH = Start(HHstallsallR,  net1, net2, net3, 'sen1dbl8DD')






#TES
def End(file, chr1, chr2, chr3, p):
    xxe =[]
    for ge in file.itertuples():
        if ge[2] == p:
            if ge[4] == 'I':
                xe = chr1.index[chr1['pos'] == ge[26]].tolist()
                xxe.append(xe) 

            if ge[4] == 'II':
                xe2 = chr2.index[chr2['pos'] == ge[26]].tolist()
                xxe.append(xe2)

            if ge[4] == 'III':
                xe3 = chr3.index[chr3['pos'] == ge[26]].tolist()
                xxe.append(xe3)
    return xxe


###
Fsenforxxe = End(HTstallsallF,  net1, net2, net3, 'sen1D')
FsenforxxeHH = End(HHstallsallF,  net1, net2, net3, 'sen1D')
Rsenforxxe = End(HTstallsallR,  net1, net2, net3, 'sen1D')
RsenforxxeHH = End(HHstallsallR,  net1, net2, net3, 'sen1D')

#senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')

Fdbl8forxxe = End(HTstallsallF,  net1, net2, net3, 'dbl8D')
Fdbl8forxxeHH = End(HHstallsallF,  net1, net2, net3, 'dbl8D')
Rdbl8forxxe = End(HTstallsallR,  net1, net2, net3, 'dbl8D')
Rdbl8forxxeHH = End(HHstallsallR,  net1, net2, net3, 'dbl8D')

FdsforxxeR = End(HTstallsallF,  net1, net2, net3, 'sen1dbl8DD')
FdsforxxeRHH = End(HHstallsallF,  net1, net2, net3, 'sen1dbl8DD')
RdsforxxeR = End(HTstallsallR,  net1, net2, net3, 'sen1dbl8DD')
RdsforxxeRHH = End(HHstallsallR,  net1, net2, net3, 'sen1dbl8DD')



def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'density'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

#sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
Fsen1stallssen1 = Gene_bits(Fsenforxxx, Fsenforxxe, netho)
Rsen1stallssen1 = Gene_bits(Rsenforxxx, Rsenforxxe, netho)


#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)

Fdbl8stallsdbl8 = Gene_bits(Fdbl8forxx, Fdbl8forxxe, netho)
Rdbl8stallsdbl8 = Gene_bits(Rdbl8forxx, Rdbl8forxxe, netho)


#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)

FdsstallsdsR = Gene_bits(FdsforxxR, FdsforxxeR, netho)
RdsstallsdsR = Gene_bits(RdsforxxR, RdsforxxeR, netho)

######HEAD to HEAD
Fsen1stallssen1HH = Gene_bits(FsenforxxHH, FsenforxxeHH, netho)
Rsen1stallssen1HH = Gene_bits(RsenforxxHH, RsenforxxeHH, netho)


#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)

Fdbl8stallsdbl8HH = Gene_bits(Fdbl8forxxHH, Fdbl8forxxeHH, netho)
Rdbl8stallsdbl8HH = Gene_bits(Rdbl8forxxHH, Rdbl8forxxeHH, netho)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)

#FdsstallsdsRHH = Gene_bits(FdsforxxRHH, FdsforxxeRHH, netho)
RdsstallsdsRHH = Gene_bits(RdsforxxRHH, RdsforxxeRHH, netho)


Fsen1stallssen1
Rsen1stallssen1

Fdbl8stallsdbl8
Rdbl8stallsdbl8

FdsstallsdsR
RdsstallsdsR

Fsen1stallssen1HH
Rsen1stallssen1HH

Fdbl8stallsdbl8HH
Rdbl8stallsdbl8HH

FdsstallsdsRHH
RdsstallsdsRHH



def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 250/length
        print (num*length)

        te = np.arange(0,250,num, dtype = int)

        expand = np.zeros(250)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata


#sen1
#aw = Expand_genes(sen1stallswt)
assF = Expand_genes(Fsen1stallssen1)
assR = Expand_genes(Rsen1stallssen1)
#dbl8
dbdF = Expand_genes(Fdbl8stallsdbl8)
dbdR = Expand_genes(Rdbl8stallsdbl8)
#ds right 
dsdsRF = Expand_genes(FdsstallsdsR)
dsdsRR = Expand_genes(RdsstallsdsR)
###HEADTO HEAD
assHHF = Expand_genes(Fsen1stallssen1HH)
assHHR = Expand_genes(Rsen1stallssen1HH)
#dbl8
dbdHHF = Expand_genes(Fdbl8stallsdbl8HH)
dbdHHR = Expand_genes(Rdbl8stallsdbl8HH)
#ds right 
#dsdsRHHF = Expand_genes(FdsstallsdsRHH)
dsdsRHHR = Expand_genes(RdsstallsdsRHH)


fuckkkkkk, (ax4, ax6, ax8) = plt.subplots(3,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(ass, cmap = 'YlOrRd', ax=ax4, cbar=True, cbar_ax=cbar_ax,  vmin= -0, vmax=5)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(dbd, cmap = 'YlOrRd', ax=ax6,cbar=True, cbar_ax=cbar_ax, vmin= -0, vmax=5 )
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(dsdsR, cmap = 'YlOrRd', ax=ax8,cbar=True, cbar_ax=cbar_ax, vmin= -0, vmax=5)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,11, 31, 42])
ax8.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])





    fig = plt.figure() #create a figure
   # fig, ax = plt.subplots(figsize=(7, 11))
  #  fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14, ax15)) = plt.subplots(5, 3, figsize=(7, 11))
    gs = fig.add_gridspec(190, 29)  #sets the figure width and height (can only refer to these as integers)
   # cbar_ax = fig.add_axes(gs[10:5, 5:8])
    cbar_ax = fig.add_axes([.91, .3, .03, .4])

    ax1 = fig.add_subplot(gs[10:24, 5:10])
    sns.heatmap(assF, cmap = 'YlOrRd', ax=ax1, cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8 )
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    ax2 = fig.add_subplot(gs[26:63, 5:10])
    sns.heatmap(assR, cmap = 'YlOrRd', ax=ax2, cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8 )
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis


# ax3.set_title('difference right forks',loc = 'left', pad = 0) 
    ax4 = fig.add_subplot(gs[65:103, 5:10])
    sns.heatmap(assHHF, cmap = 'YlOrRd', ax=ax4, cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    ax5 = fig.add_subplot(gs[105:130, 5:10])
    sns.heatmap(assHHR, cmap = 'YlOrRd', ax=ax5, cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8)
    ax5.set_xticks([])
    ax5.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

 
    ax8 = fig.add_subplot(gs[135:140, 5:10])
    sns.heatmap(dbdF, cmap = 'YlOrRd', ax=ax8,  cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8)
    ax8.set_xticks([])
    ax8.set_yticks([])
    
    ax9 = fig.add_subplot(gs[142:149, 5:10])
    sns.heatmap(dbdR, cmap = 'YlOrRd', ax=ax9,  cbar=True, cbar_ax=cbar_ax,vmin= 0, vmax=8)
    ax9.set_xticks([])
    ax9.set_yticks([])
    

    ax11 = fig.add_subplot(gs[151:159, 5:10])
    sns.heatmap(dbdHHF, cmap = 'YlOrRd', ax=ax11,  cbar=True, cbar_ax=cbar_ax, vmin= 0, vmax=8)
    ax11.set_xticks([])
    ax11.set_yticks([])
    
    ax12 = fig.add_subplot(gs[161:169, 5:10])
    sns.heatmap(dbdHHR, cmap = 'YlOrRd', ax=ax12,  cbar=True, cbar_ax=cbar_ax, vmin= 0, vmax=8)
    ax12.set_xticks([])
    ax12.set_yticks([])    
    
    
    
    
   # ax15 = fig.add_subplot(gs[166:169, 5:10]) 
   # sns.heatmap(dsdsR, cmap = 'YlOrRd', ax=ax15,cbar=True, cbar_ax=cbar_ax, vmin= 0, vmax=5)
    
    
    
    
    ax12 = fig.add_subplot(gs[174:175, 5:10])
    sns.heatmap(dsdsRF, cmap = 'YlOrRd', ax=ax12,  cbar=True, cbar_ax=cbar_ax, vmin= 0, vmax=8)
    ax12.set_xticks([])
    ax12.set_yticks([])
    
    ax13 = fig.add_subplot(gs[177:180, 5:10])
    sns.heatmap(dsdsRR, cmap = 'YlOrRd', ax=ax13,  cbar=True, cbar_ax=cbar_ax, vmin= 0, vmax=8)
    ax13.set_xticks([])
    ax13.set_yticks([])
    
 #   ax15 = fig.add_subplot(gs[163:164, 5:10]) 
 #   sns.heatmap(dsdsRHHF, cmap = 'YlOrRd', ax=ax15,cbar=True, cbar_ax=cbar_ax,  vmin= 0, vmax=8)
 #   ax15.set_xticks([])
  #  ax15.set_yticks([])
    
    ax16 = fig.add_subplot(gs[182:186, 5:10]) 
    sns.heatmap(dsdsRHHR, cmap = 'YlOrRd', ax=ax16,cbar=True, cbar_ax=cbar_ax,  vmin= 0, vmax=8)


     
    ax16.set_xticks([0,250])
    ax16.set_xticklabels(['s', 'e'])
    ax16.set_yticks([])



#%%


#updated to take 3kb downstream, thankyou openai 
def Start(file, chr1, chr2, chr3, p):
    xx = []
    
    for g in file.itertuples():
        if g[2] == p:
            if g[4] == 'I':
                x = chr1.index[chr1['pos'] == g[5]].tolist()
                if x:
                    x = [max(0, pos - 10) for pos in x]
                    xx.append(x)

            elif g[4] == 'II':
                x2 = chr2.index[chr2['pos'] == g[5]].tolist()
                if x2:
                    x2 = [max(0, pos - 10) for pos in x2]
                    xx.append(x2)

            elif g[4] == 'III':
                x3 = chr3.index[chr3['pos'] == g[5]].tolist()
                if x3:
                    x3 = [max(0, pos - 10) for pos in x3]
                    xx.append(x3)
                
    return xx

senforxx = Start(deltastallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
dbl8forxx = Start(deltastallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxR = Start(deltastallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')

def End(file, chr1, chr2, chr3, p):
    xxe = []
    
    for ge in file.itertuples():
        if ge[2] == p:
            if ge[4] == 'I':
                xe = chr1.index[chr1['pos'] == ge[6]].tolist()
                if xe:
                    xe = [pos + 10 for pos in xe]
                    xxe.append(xe)

            elif ge[4] == 'II':
                xe2 = chr2.index[chr2['pos'] == ge[6]].tolist()
                if xe2:
                    xe2 = [pos + 10 for pos in xe2]
                    xxe.append(xe2)

            elif ge[4] == 'III':
                xe3 = chr3.index[chr3['pos'] == ge[6]].tolist()
                if xe3:
                    xe3 = [pos + 10 for pos in xe3]
                    xxe.append(xe3)
                
    return xxe

senforxxe = End(deltastallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
dbl8forxxe = End(deltastallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxeR = End(deltastallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD')



def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'differentials'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

#sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1 = Gene_bits(senforxx, senforxxe, esen1)
sen1stallsdbl8 = Gene_bits(senforxx, senforxxe, edbl8)
sen1stallsds = Gene_bits(senforxx, senforxxe, eds)

#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1 = Gene_bits(dbl8forxx, dbl8forxxe, esen1)
dbl8stallsdbl8 = Gene_bits(dbl8forxx, dbl8forxxe, edbl8)
dbl8stallsds = Gene_bits(dbl8forxx, dbl8forxxe, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1R = Gene_bits(dsforxxR, dsforxxeR, esen1)
dsstallsdbl8R = Gene_bits(dsforxxR, dsforxxeR, edbl8)
dsstallsdsR = Gene_bits(dsforxxR, dsforxxeR, eds)

######HEAD to HEAD
sen1stallssen1HH = Gene_bits(senforxxHH, senforxxeHH, esen1)
sen1stallsdbl8HH = Gene_bits(senforxxHH, senforxxeHH, edbl8)
sen1stallsdsHH = Gene_bits(senforxxHH, senforxxeHH, eds)

#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1HH = Gene_bits(dbl8forxxHH, dbl8forxxeHH, esen1)
dbl8stallsdbl8HH = Gene_bits(dbl8forxxHH, dbl8forxxeHH, edbl8)
dbl8stallsdsHH = Gene_bits(dbl8forxxHH, dbl8forxxeHH, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1RHH = Gene_bits(dsforxxRHH, dsforxxeRHH, esen1)
dsstallsdbl8RHH = Gene_bits(dsforxxRHH, dsforxxeRHH, edbl8)
dsstallsdsRHH = Gene_bits(dsforxxRHH, dsforxxeRHH, eds)




def ahead_Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]-10:z[0],'differentials'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

sen1stallssen1before = ahead_Gene_bits(senforxx, senforxxe, esen1)
sen1stallsdbl8before = ahead_Gene_bits(senforxx, senforxxe, edbl8)
sen1stallsdsbefore = ahead_Gene_bits(senforxx, senforxxe, eds)

dbl8stallssen1before = ahead_Gene_bits(dbl8forxx, dbl8forxxe, esen1)
dbl8stallsdbl8before = ahead_Gene_bits(dbl8forxx, dbl8forxxe, edbl8)
dbl8stallsdsbefore = ahead_Gene_bits(dbl8forxx, dbl8forxxe, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1Rbefore = ahead_Gene_bits(dsforxxR, dsforxxeR, esen1)
dsstallsdbl8Rbefore = ahead_Gene_bits(dsforxxR, dsforxxeR, edbl8)
dsstallsdsRbefore = ahead_Gene_bits(dsforxxR, dsforxxeR, eds)

###HEAD TO HEAD
sen1stallssen1beforeHH = ahead_Gene_bits(senforxxHH, senforxxeHH, esen1)
sen1stallsdbl8beforeHH = ahead_Gene_bits(senforxxHH, senforxxeHH, edbl8)
sen1stallsdsbeforeHH = ahead_Gene_bits(senforxxHH, senforxxeHH, eds)

dbl8stallssen1beforeHH = ahead_Gene_bits(dbl8forxxHH, dbl8forxxeHH, esen1)
dbl8stallsdbl8beforeHH = ahead_Gene_bits(dbl8forxxHH, dbl8forxxeHH, edbl8)
dbl8stallsdsbeforeHH = ahead_Gene_bits(dbl8forxxHH, dbl8forxxeHH, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1RbeforeHH = ahead_Gene_bits(dsforxxRHH, dsforxxeRHH, esen1)
dsstallsdbl8RbeforeHH = ahead_Gene_bits(dsforxxRHH, dsforxxeRHH, edbl8)
dsstallsdsRbeforeHH = ahead_Gene_bits(dsforxxRHH, dsforxxeRHH, eds)



def behind_Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[ze[0]:ze[0]+10,'differentials'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

#sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1behind = behind_Gene_bits(senforxx, senforxxe, esen1)
sen1stallsdbl8behind = behind_Gene_bits(senforxx, senforxxe, edbl8)
sen1stallsdsbehind = behind_Gene_bits(senforxx, senforxxe, eds)

#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1behind = behind_Gene_bits(dbl8forxx, dbl8forxxe, esen1)
dbl8stallsdbl8behind = behind_Gene_bits(dbl8forxx, dbl8forxxe, edbl8)
dbl8stallsdsbehind = behind_Gene_bits(dbl8forxx, dbl8forxxe, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1Rbehind = behind_Gene_bits(dsforxxR, dsforxxeR, esen1)
dsstallsdbl8Rbehind = behind_Gene_bits(dsforxxR, dsforxxeR, edbl8)
dsstallsdsRbehind = behind_Gene_bits(dsforxxR, dsforxxeR, eds)

#HEAD TO HEAD 
#sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1behindHH = behind_Gene_bits(senforxxHH, senforxxeHH, esen1)
sen1stallsdbl8behindHH = behind_Gene_bits(senforxxHH, senforxxeHH, edbl8)
sen1stallsdsbehindHH = behind_Gene_bits(senforxxHH, senforxxeHH, eds)

#dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1behindHH = behind_Gene_bits(dbl8forxxHH, dbl8forxxeHH, esen1)
dbl8stallsdbl8behindHH = behind_Gene_bits(dbl8forxxHH, dbl8forxxeHH, edbl8)
dbl8stallsdsbehindHH = behind_Gene_bits(dbl8forxxHH, dbl8forxxeHH, eds)

#for right double 
#dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1RbehindHH = behind_Gene_bits(dsforxxRHH, dsforxxeRHH, esen1)
dsstallsdbl8RbehindHH = behind_Gene_bits(dsforxxRHH, dsforxxeRHH, edbl8)
dsstallsdsRbehindHH = behind_Gene_bits(dsforxxRHH, dsforxxeRHH, eds)



#%%
'''
def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 1000/length
        print (num*length)

        te = np.arange(0,1000,num, dtype = int)

        expand = np.zeros(1000)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata

aad = Expand_genes(sen1stalls)
'''
def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 20/length
        print (num*length)

        te = np.arange(0,20,num, dtype = int)

        expand = np.zeros(20)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata


#sen1
#aw = Expand_genes(sen1stallswt)
ass = Expand_genes(sen1stallssen1)
ad = Expand_genes(sen1stallsdbl8)
ads = Expand_genes(sen1stallsds)


#dbl8
#dbw = Expand_genes(dbl8stallswt)
dbss = Expand_genes(dbl8stallssen1)
dbd = Expand_genes(dbl8stallsdbl8)
dbds = Expand_genes(dbl8stallsds)


#ds right 
#dswR = Expand_genes(dsstallswtR)
dssR = Expand_genes(dsstallssen1R)
dsbR = Expand_genes(dsstallsdbl8R)
dsdsR = Expand_genes(dsstallsdsR)



###HEADTO HEAD

assHH = Expand_genes(sen1stallssen1HH)
adHH = Expand_genes(sen1stallsdbl8HH)
adsHH = Expand_genes(sen1stallsdsHH)


#dbl8
#dbw = Expand_genes(dbl8stallswt)
dbssHH = Expand_genes(dbl8stallssen1HH)
dbdHH = Expand_genes(dbl8stallsdbl8HH)
dbdsHH = Expand_genes(dbl8stallsdsHH)


#ds right 
#dswR = Expand_genes(dsstallswtR)
dssRHH = Expand_genes(dsstallssen1RHH)
dsbRHH = Expand_genes(dsstallsdbl8RHH)
dsdsRHH = Expand_genes(dsstallsdsRHH)





#%%



def concati(before, peaky, behind):
    
    peakys = np.stack(peaky, axis=0)
    peakys_df = pd.DataFrame(peakys)
    
    new = pd.concat([before, peakys_df, behind], axis=1)

    return new

#reverse LEFT, forward right 

HTs = concati(dsstallssen1Rbefore, dssR, dsstallssen1Rbehind)
HTb = concati(dsstallsdbl8Rbefore, dsbR, dsstallsdbl8Rbehind)
HTds = concati (dsstallsdsRbefore ,dsdsR, dsstallsdsRbehind)

sas = concati(sen1stallssen1before,ass,sen1stallssen1behind)
sab = concati(sen1stallsdbl8before,ad, sen1stallsdbl8behind)
sads = concati(sen1stallsdsbefore,ads, sen1stallsdsbehind)

dbl8s = concati(dbl8stallssen1before, dbss, dbl8stallssen1behind)
dbl8b = concati(dbl8stallsdbl8before, dbd, dbl8stallsdbl8behind)
dbl8ds = concati(dbl8stallsdsbefore, dbds, dbl8stallsdsbehind)


sasHH = concati(sen1stallssen1beforeHH,assHH,sen1stallssen1behindHH)
sabHH = concati(sen1stallsdbl8beforeHH,adHH, sen1stallsdbl8behindHH)
sadsHH = concati(sen1stallsdsbeforeHH,adsHH, sen1stallsdsbehindHH)

dbl8sHH = concati(dbl8stallssen1beforeHH, dbssHH, dbl8stallssen1behindHH)
dbl8bHH = concati(dbl8stallsdbl8beforeHH, dbdHH, dbl8stallsdbl8behindHH)
dbl8dsHH = concati(dbl8stallsdsbeforeHH, dbdsHH, dbl8stallsdsbehindHH)






#%%

import scipy.stats as stats
from scipy.stats import t

x = pd.DataFrame(sa)
m = x[0].mean()

s = x[0].std() 
dof = len(x[0])-1 
confidence = 0.95
t_crit = np.abs(t.ppf((1-confidence)/2,dof))
print(t_crit)
(m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x))) 


def conint(array):
    confidence = 0.95
    trythis = []
    x = pd.DataFrame(array)
    #print(x)
    for column in x:

        m = (x[column]).mean()
        s = x[column].std()

        dof = len(x[column])-1 
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        interval = (m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x)))
        
        trythis.append(interval)
        saint = pd.DataFrame(trythis, columns=['lower', 'upper'])
        #oft = saint.T
    return saint

saCON = conint(sa)
sasCON = conint(sas)
sabCON = conint(sab)
sadsCON = conint(sads)

bwCON = conint(dbl8w)
bsCON = conint(dbl8s)
bbCON = conint(dbl8b)
bdsCON = conint(dbl8ds)

HHwCON = conint(HHw)
HHsCON = conint(HHs)
HHbCON = conint(HHb)
HHdsCON = conint(HHds)


HTwCON = conint(HTw)
HTsCON = conint(HTs)
HTbCON = conint(HTb)
HTdsCON = conint(HTds)

cawCON = conint(caw)
casCON = conint(cas)
cabCON = conint(cab)
cadsCON = conint(cads)

sewCON = conint(sew)
sesCON = conint(ses)
sebCON = conint(seb)
sedsCON = conint(seds)


        
    
    
def pileup (thing):
    #thing = -thing
    #print(thing)
    #print(-thing)
    whack = np.stack( thing, axis=0 )
    stuff = whack.mean(axis=0)
    wow = pd.DataFrame(stuff)
    
    return wow
#sen compiled


#hh ht 
HHwline = pileup(HHw)
HHsline = pileup(HHs)
HHbline = pileup(HHb)
HHdsline = pileup(HHds)
HTwline = pileup(HTw)
HTsline = pileup(HTs)
HTbline = pileup(HTb)
HTdsline = pileup(HTds)

#selection
swline = pileup(sew)
ssline = pileup(ses)
sbline = pileup(seb)
sdsline = pileup(seds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)

rcwline = pileup(cwr)
rcsline = pileup(csr)
rcbline = pileup(cbr)
rcdsline = pileup(cdsr)


#sen
awline = pileup(sa)
assline = pileup(sas)
adline = pileup(sab)
adsline= pileup(sads)


#dbl8
dbwline = pileup(dbl8w)
dbssline = pileup(dbl8s)
dbdline = pileup(dbl8b)
dbdsline= pileup(dbl8ds)
#%%

import seaborn as sns


fuckkkkkk, (ax4, ax6, ax8) = plt.subplots(3,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.71, .3, .03, .4])
sns.heatmap(HTs, cmap = 'coolwarm', ax=ax4, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(HTb, cmap = 'coolwarm', ax=ax6,cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(HTds, cmap = 'coolwarm', ax=ax8,cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,11, 31, 42])
ax8.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])




#DBL8 heatmap !!!!!
fxb, (( ax3, ax5,ax7)) = plt.subplots(3,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])

#ax1.set_title('Forward Strand')
#ax1.set_ylabel('WT')
sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax3,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.5)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.5)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax7,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.5)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,11, 31, 42])
ax7.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])

#SEN1 heatmap bitches 
fx, (( ax3, ax5,ax7)) = plt.subplots(3,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.91, .3, .03, .4])


sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.7)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.7)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.7)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')


ax7.set_xticks([0,11, 31, 42])
ax7.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])




HTs = concati(dsstallssen1Rbefore, dssR, dsstallssen1Rbehind)
HTb = concati(dsstallsdbl8Rbefore, dsbR, dsstallsdbl8Rbehind)
HTds = concati (dsstallsdsRbefore ,dsdsR, dsstallsdsRbehind)

sas = concati(sen1stallssen1before,ass,sen1stallssen1behind)
sab = concati(sen1stallsdbl8before,ad, sen1stallsdbl8behind)
sads = concati(sen1stallsdsbefore,ads, sen1stallsdsbehind)

dbl8s = concati(dbl8stallssen1before, dbss, dbl8stallssen1behind)
dbl8b = concati(dbl8stallsdbl8before, dbd, dbl8stallsdbl8behind)
dbl8ds = concati(dbl8stallsdsbefore, dbds, dbl8stallsdsbehind)


sasHH = concati(sen1stallssen1beforeHH,assHH,sen1stallssen1behindHH)
sabHH = concati(sen1stallsdbl8beforeHH,adHH, sen1stallsdbl8behindHH)
sadsHH = concati(sen1stallsdsbeforeHH,adsHH, sen1stallsdsbehindHH)

dbl8sHH = concati(dbl8stallssen1beforeHH, dbssHH, dbl8stallssen1behindHH)
dbl8bHH = concati(dbl8stallsdbl8beforeHH, dbdHH, dbl8stallsdbl8behindHH)
dbl8dsHH = concati(dbl8stallsdsbeforeHH, dbdsHH, dbl8stallsdsbehindHH)


#%%

arr_range = np.ptp(arr) 

    fig = plt.figure() #create a figure
   # fig, ax = plt.subplots(figsize=(7, 11))
  #  fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14, ax15)) = plt.subplots(5, 3, figsize=(7, 11))
    gs = fig.add_gridspec(190, 29)  #sets the figure width and height (can only refer to these as integers)
#    cbar_ax = fig.add_axes(gs[0:5, 5:8])

    ax1 = fig.add_subplot(gs[10:51, 5:10])
    sns.heatmap(sas, cmap = 'coolwarm', ax=ax1, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax1.set_xticks([])
    ax1.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:51, 12:17]) 
    sns.heatmap(sab, cmap = 'coolwarm', ax=ax2,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:51, 19:24]) 
    sns.heatmap(sads, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax3.set_xticks([])
    ax3.set_yticks([])

# ax3.set_title('difference right forks',loc = 'left', pad = 0) 
    ax4 = fig.add_subplot(gs[56:119, 5:10])
    sns.heatmap(sasHH, cmap = 'coolwarm', ax=ax4, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax4.set_xticks([])
    ax4.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax5 = fig.add_subplot(gs[56:119, 12:17]) 
    sns.heatmap(sabHH, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax5.set_xticks([])
    ax5.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax6 = fig.add_subplot(gs[56:119, 19:24]) 
    sns.heatmap(sadsHH, cmap = 'coolwarm', ax=ax6, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax6.set_xticks([])
    ax6.set_yticks([])   
   
   
   
    
    ax7 = fig.add_subplot(gs[124:136, 5:10])
    sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax7,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax7.set_xticks([])
    ax7.set_yticks([])
    
    ax8 = fig.add_subplot(gs[124:136, 12:17])
    sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax8,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax8.set_xticks([])
    ax8.set_yticks([])
    
    ax9 = fig.add_subplot(gs[124:136, 19:24])
    sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax9,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax9.set_xticks([])
    ax9.set_yticks([])
    
    
    
    ax10 = fig.add_subplot(gs[141:157, 5:10])
    sns.heatmap(dbl8sHH, cmap = 'coolwarm', ax=ax10,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax10.set_xticks([])
    ax10.set_yticks([])
    
    ax11 = fig.add_subplot(gs[141:157, 12:17])
    sns.heatmap(dbl8bHH, cmap = 'coolwarm', ax=ax11,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax11.set_xticks([])
    ax11.set_yticks([])
    
    ax12 = fig.add_subplot(gs[141:157, 19:24])
    sns.heatmap(dbl8dsHH, cmap = 'coolwarm', ax=ax12,  cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax12.set_xticks([])
    ax12.set_yticks([])
        
    
    
    
    ax13 = fig.add_subplot(gs[162:170, 5:10])
    sns.heatmap(HTs, cmap = 'coolwarm', ax=ax13, cbar=True, cbar_ax=cbar_ax, vmin= -0.3, vmax=0.7)
    ax13.set_xticks([0,11, 31, 42])
    ax13.set_xticklabels(['-','s', 'e', '+'])
    ax13.set_yticks([])
    
    
    ax14 = fig.add_subplot(gs[162:170, 12:17]) 
    sns.heatmap(HTb, cmap = 'coolwarm', ax=ax14,cbar=True, cbar_ax=cbar_ax,  vmin= -0.3, vmax=0.7)
    ax14.set_xticks([0,11, 31, 42])
    ax14.set_xticklabels(['-','s', 'e', '+'])
    ax14.set_yticks([])
    
    ax15 = fig.add_subplot(gs[162:170, 19:24]) 
    sns.heatmap(HTds, cmap = 'coolwarm', ax=ax15,cbar=True, cbar_ax=cbar_ax,  vmin= -0.3, vmax=0.7)
    ax15.set_xticks([0,11, 31, 42])
    ax15.set_xticklabels(['-','s', 'e', '+'])
    ax15.set_yticks([])








#%%

#so this is for from before, when i worked without the HH and HT 

    fig = plt.figure(constrained_layout=True) #create a figure
    gs = fig.add_gridspec(185, 29)  #sets the figure width and height (can only refer to these as integers)




# Now add as many subplots using the gridspace (gs) to position them precisely
# the first position is the height (ie here from 1 to 3 on the grid)
# the second position is the width (i.e here from 1 to 10 on the grid)
    ax1 = fig.add_subplot(gs[10:127, 5:10])
    sns.heatmap(sas, cmap = 'coolwarm', ax=ax1, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax1.set_xticks([])
    ax1.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:127, 12:17]) 
    sns.heatmap(sab, cmap = 'coolwarm', ax=ax2,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:127, 19:24]) 
    sns.heatmap(sads, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax3.set_xticks([])
    ax3.set_yticks([])

   # ax3.set_title('difference right forks',loc = 'left', pad = 0) 
   
   
   
   
    
    ax4 = fig.add_subplot(gs[132:152, 5:10])
    sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax4,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    ax5 = fig.add_subplot(gs[132:152, 12:17])
    sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    ax6 = fig.add_subplot(gs[132:152, 19:24])
    sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax6,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax6.set_xticks([])
    ax6.set_yticks([])
    
    
    
    
    
    ax7 = fig.add_subplot(gs[157:165, 5:10])
    sns.heatmap(HTs, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
    ax7.set_xticks([0,11, 31, 42])
    ax7.set_xticklabels(['-','s', 'e', '+'])
    ax7.set_yticks([])
    
    
    ax8 = fig.add_subplot(gs[157:165, 12:17]) 
    sns.heatmap(HTb, cmap = 'coolwarm', ax=ax8,cbar=True, cbar_ax=cbar_ax,  vmin= -0.1, vmax=0.1)
    ax8.set_xticks([0,11, 31, 42])
    ax8.set_xticklabels(['-','s', 'e', '+'])
    ax8.set_yticks([])
    
    ax9 = fig.add_subplot(gs[157:165, 19:24]) 
    sns.heatmap(HTds, cmap = 'coolwarm', ax=ax9,cbar=True, cbar_ax=cbar_ax,  vmin= -0.1, vmax=0.1)
    ax9.set_xticks([0,11, 31, 42])
    ax9.set_xticklabels(['-','s', 'e', '+'])
    ax9.set_yticks([])
    
    
#%%
deltastallsall






