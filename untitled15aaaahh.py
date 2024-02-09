#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 15:42:07 2023

@author: patricfernandez
"""

#let's call this trying to adapt stall call for calling delta bias 


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
    all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = 7, center=True).mean()
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
    

owtchrI, owtchrII, owtchrIII, owt = Create_df('Jo_wt_YE.e1.f-w300.count.csv', 'Jo_wt_YE.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 3)
#sen1chrI, sen1chrII, sen1chrIII 
osen1chr1, osen1chr2, osen1chr3, osen = Create_df ('Jo_sen1d_YE.e1.f-w300.count.csv', 'Jo_sen1d_YE.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 3)
odbl8chr1, odbl8chr2, odbl8chr3, odbl = Create_df('Jo_dbl8_YE.e1.f-w300.count.csv', 'Jo_dbl8_YE.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 3)
odschr1, odschr2, odschr3, ods = Create_df('Jo_sen1_dbl8_YE.e1.f-w300.count.csv', 'Jo_sen1_dbl8_YE.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 3)


owtchrI, owtchrII, owtchrIII, owt = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)
#sen1chrI, sen1chrII, sen1chrIII 
osen1chr1, osen1chr2, osen1chr3 = Create_df ('RZ261-RhArepr-d.e1.f-w300.count.csv', 'RZ261-RhArepr-d.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 5)
odbl8chr1, odbl8chr2, odbl8chr3 = Create_df('RZ263-RhArepr-d.e1.f-w300.count.csv', 'RZ263-RhArepr-d.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 5)
odschr1, odschr2, odschr3 = Create_df('RZ265-RhArepr-d.e1.f-w300.count.csv', 'RZ265-RhArepr-d.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 5)



def dusebaby (data, wt):
    data['d_usage_norm'] = data['d_usage'] - wt['d_usage']
    data = data.loc[data['d_usage_norm'] > 0]
   # data['d_usage_norm'].clip(0, 1, inplace=True)
    #data = data.dropna()
 #   templist = np.array([])
    
  #  for i in range(1,len(data['d_usage_norm'])):
   #     if data['d_usage_norm'].iloc[i] != 0: 
    #        templist = np.append(templist, data['d_usage_norm'].iloc[i])
    
   # cutoff = np.quantile(templist, 0.98)
 #   for i in range(1, len(data['d_usage_norm'])-1):
  #      if data['d_usage_norm'].iloc[i] != 0:
   #         if data['d_usage_norm'].iloc[i]< cutoff:
    #            data['d_usage_norm'].iloc[i] = 0
    return data

osen1chr1 = dusebaby(osen1chr1, owtchrI)
osen1chr2 = dusebaby(osen1chr2, owtchrII)
osen1chr3 = dusebaby(osen1chr3, owtchrIII)

odbl8chr1 = dusebaby(odbl8chr1, owtchrI)
odbl8chr2 = dusebaby(odbl8chr2, owtchrII)
odbl8chr3 = dusebaby(odbl8chr3, owtchrIII)

odschr1 = dusebaby(odschr1, owtchrI)
odschr2 = dusebaby(odschr2, owtchrII)
odschr3 = dusebaby(odschr3, owtchrIII)




def remove_telo (data, start, end):
    
    data = data[(data[('pos')] > start) & (data[('pos')] < end)]
    
    return data


owtchrI = remove_telo(owtchrI, 29663,5554844)
osen1chr1 = remove_telo(osen1chr1, 29663,5554844)
odbl8chr1 = remove_telo(odbl8chr1, 29663,5554844)
odschr1 = remove_telo(odschr1, 29663,5554844)

owtchrII = remove_telo(owtchrII, 39186,4500619)
osen1chr2 = remove_telo(osen1chr2, 39186,4500619)
odbl8chr2 = remove_telo(odbl8chr2, 39186,4500619)
odschr2 = remove_telo(odschr2, 39186,4500619)

owtchrIII = remove_telo(owtchrIII, 13722,2439540)
osen1chr3 = remove_telo(osen1chr3, 13722,2439540)
odbl8chr3 = remove_telo(odbl8chr3, 13722,2439540)
odschr3 = remove_telo(odschr3, 13722,2439540)



def remove_mat_loc (data, start, end):
    
    data = data[(data[('pos')] < start) | (data[('pos')] > end)]
    
    return data

#mat loc
owtchrII = remove_mat_loc(owtchrII, 2114008, 2115135)
osen1chr2 = remove_mat_loc(osen1chr2, 2114008,2115135)
odbl8chr2 = remove_mat_loc(odbl8chr2, 2114008,2115135)
odschr2 = remove_mat_loc(odschr2, 2114008,2115135)

#remove dbl8
owtchrII = remove_mat_loc(owtchrII, 2555350, 2555350)
osen1chr2 = remove_mat_loc(osen1chr2, 2555350,2558950)
odbl8chr2 = remove_mat_loc(odbl8chr2, 2555350,2558950)
odschr2 = remove_mat_loc(odschr2, 2555350,2558950)

#remove cen2
owtchrII = remove_mat_loc(owtchrII, 1602264, 1644747)
osen1chr2 = remove_mat_loc(osen1chr2, 1602264,1644747)
odbl8chr2 = remove_mat_loc(odbl8chr2, 1602264,1644747)
odschr2 = remove_mat_loc(odschr2, 1602264,1644747)


def remove_chr1 (data, start, end):
    
    data = data[(data[('pos')] < start) | (data[('pos')] > end)]
    
    return data

owtchrI = remove_chr1(owtchrI, 3260605,3264205)
osen1chr1 = remove_chr1(osen1chr1, 3260605,3264205)
odbl8chr1 = remove_chr1(odbl8chr1, 3260605,3264205)
odschr1 = remove_chr1(odschr1, 3260605,3264205)

#centromwere 1]
owtchrI = remove_chr1(owtchrI, 3753687,3789421)
osen1chr1 = remove_chr1(osen1chr1, 3753687,3789421)
odbl8chr1 = remove_chr1(odbl8chr1, 3753687,3789421)
odschr1 = remove_chr1(odschr1, 3753687,3789421)



def remove_chr3 (data, start, end):
    
    data = data[(data[('pos')] < start) | (data[('pos')] > end)]
    
    return data

osen1chr3 = remove_chr3(osen1chr3, 1070904,1137003)
odbl8chr3 = remove_chr3(odbl8chr3, 1070904,1137003)
odschr3 = remove_chr3(odschr3, 1070904,1137003)

owtchrIII = remove_chr3(owtchrIII, 1070904,1137003)

#%%
#please note:
    #there is an issue with origin calling using YE and EMM data together 
    #note that some are shifted/don't read out properly. but they are correct (?) positions
#first thing's first, let's make sure origins of all the genotypes align
#
percentile_value = ods['origin'].quantile(0.98)
owt.loc[owt['origin'] < percentile_value, 'origin'] = 0


osen.loc[osen['origin'] < percentile_value, 'origin'] = 0
odbl.loc[odbl['origin'] < percentile_value, 'origin'] = 0
ods.loc[ods['origin'] < percentile_value, 'origin'] = 0


owt1 = owt[(owt['chro'] == 'chr1')]
owt2 = owt[(owt['chro'] == 'chr2')]
owt3 = owt[(owt['chro'] == 'chr3')]

osen1 = osen[(osen['chro'] == 'chr1')]
osen2 = osen[(osen['chro'] == 'chr2')]
osen3 = osen[(osen['chro'] == 'chr3')]

odbl1 = odbl[(odbl['chro'] == 'chr1')]
odbl2 = odbl[(odbl['chro'] == 'chr2')]
odbl3 = odbl[(odbl['chro'] == 'chr3')]

ods1 = ods[(ods['chro'] == 'chr1')]
ods2 = ods[(ods['chro'] == 'chr2')]
ods3 = ods[(ods['chro'] == 'chr3')]




def Chromosome_plot (data1, data2, data3, data4, featurex, centro, telo, genee, c, con,data5, data6, datat, odata1, odata2, odata3, odata4):
    ff, (ax1, ax3, ax5) = plt.subplots(3,1, sharex=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['smoo_right_forks'], color ='black', alpha=0.8)
    ax1.plot(data2['pos'], data2['smoo_right_forks'], color ='steelblue', alpha=0.8)
    ax1.plot(data3['pos'], data3['smoo_right_forks'], color ='orange', alpha=0.8)
    ax1.plot(data4['pos'], data4['smoo_right_forks'], color ='tomato', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('rightward forks (e)')


    ax5.set_xlabel('Chromosome position')

    for be in odata1.itertuples(index=False, name = None):
        if be[23] > 0.2:
            ax1.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            ax3.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
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
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                
    for sug in data5.itertuples(index=False, name=None):
        if sug[23] > 0.1:
            ax3.axvspan(sug[1],sug[1],0.1,0.9,color="steelblue",alpha=0.9)
        
    for dug in data6.itertuples(index=False, name=None):
        if dug[23] > 0.1:
            ax3.axvspan(dug[1],dug[1],0.1,0.9,color="orange",alpha=0.9)                
        
    for dsug in datat.itertuples(index=False, name=None):
        if dsug[23] > 0.1:
            ax3.axvspan(dsug[1],dsug[1],0.1,0.9,color="red",alpha=0.9)
                
    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(owtchrI, osen1chr1, odbl8chr1, odschr1, gene1, centro1, telo1, feat1, 'I', con1, osen1, odbl1, ods1, owtchrI, osen1chr1, odbl8chr1, odschr1)



# Create a function to find unique positions where 'origin' > 0 for each DataFrame
def unique_positions(df):
    return set(df.loc[df['origin'] > 0, 'pos'])

# Get unique positions for each DataFrame
unique_positions_df = unique_positions(owt)
unique_positions_df1 = unique_positions(osen)
unique_positions_df2 = unique_positions(odbl)
unique_positions_df3 = unique_positions(ods)

# Check if all unique positions are the same for all DataFrames
if (
    unique_positions_df == unique_positions_df1
    and unique_positions_df == unique_positions_df2
    and unique_positions_df == unique_positions_df3
):
    print("All unique positions are the same across the four genotypes.")
else:
    print("Unique positions differ among the four genotypes.")




# Find unique positions that differ among the four genotypes
differ_positions = (
    unique_positions_df
    .symmetric_difference(unique_positions_df1)
    .symmetric_difference(unique_positions_df2)
    .symmetric_difference(unique_positions_df3)
)

# Create a DataFrame with the differing positions
differ_df = pd.DataFrame({'pos': list(differ_positions)})

# Print or use the 'differ_df' DataFrame
print(differ_df)



def unique_positions(df):
    return set(df.loc[df['origin'] > 0, ['pos', 'chro']].apply(tuple, axis=1))

# Get unique positions for each DataFrame
unique_positions_df = unique_positions(owt)
unique_positions_df1 = unique_positions(osen)
unique_positions_df2 = unique_positions(odbl)
unique_positions_df3 = unique_positions(ods)

# Find unique positions that differ among the four genotypes
differ_positions = (
    unique_positions_df
    .symmetric_difference(unique_positions_df1)
    .symmetric_difference(unique_positions_df2)
    .symmetric_difference(unique_positions_df3)
)

# Create a DataFrame with the differing positions and 'chro' values
differ_df = pd.DataFrame(list(differ_positions), columns=['pos', 'chro'])

# Print or use the 'differ_df' DataFrame
print(differ_df)


differ1 = differ_df[(differ_df['chro'] == 'chr1')]
differ2 = differ_df[(differ_df['chro'] == 'chr2')]
differ3 = differ_df[(differ_df['chro'] == 'chr3')]



#%%


frameseny = [osen1chr1, osen1chr2, osen1chr3]
sencut = pd.concat(frameseny)

framesdbly = [odbl8chr1,odbl8chr2,odbl8chr3]
dblcut = pd.concat(framesdbly)

framesds = [odschr1, odschr2, odschr3]
dscut = pd.concat(framesds)

sencut['source'] = 'sen1'
dblcut['source'] = 'dbl8'
dscut['source'] = 'ds'

# Concatenate the DataFrames
concatall = pd.concat([sencut, dblcut, dscut], ignore_index=True)
#%%

def Chromosome_plot (datat, data5, data6, data7):
    ff, (ax5, ax6) = plt.subplots(2,1, sharex=True)
    
    ax5.plot(data5['pos'], data5['d_usage'], color= 'blue', alpha=0.5)
    ax5.plot(data6['pos'], data6['d_usage'], color= 'orange', alpha=0.5)
    ax5.plot(datat['pos'], datat['d_usage'], color= 'black', alpha=0.5)
    ax5.plot(data7['pos'], data7['d_usage'], color= 'red', alpha=0.5)
    
    ax6.plot(data5['pos'], data5['d_usage_norm'], color= 'blue', alpha=0.5)
    ax6.plot(data6['pos'], data6['d_usage_norm'], color= 'orange', alpha=0.5)
    ax6.plot(data7['pos'], data7['d_usage_norm'], color= 'red', alpha=0.5)
    
    return ff

chro1plot = Chromosome_plot(owtchrI, osen1chr1, odbl8chr1, odschr1)


def call_peak(data, bin_size):
    
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
            
        
            data_output = pd.DataFrame()
            data_output = data.groupby('peak_no').agg({'source':'first', 'chro':'first' ,'pos':'first', 'd_usage_norm':['max', 'sum']})
            print(data_output)
            
                        
    return data, data_output

concatall, econcatall = call_peak(concatall, 300)



delta_peaks_maybe = concatall.copy().groupby(
    ["source", "chro", "peak_no"], as_index = False).agg(      
        chro = ("chro", "first"),
        start_pos = ("pos", "first"),
        end_pos = ("pos", "last"),
        peak_area = ("d_usage_norm", "sum"),
        peak_max = ("d_usage_norm", "max"))
       

delta_peaks_maybe_filt = delta_peaks_maybe.loc[
   
    (delta_peaks_maybe["peak_area"] > np.quantile(delta_peaks_maybe["peak_area"], 0.98)) &
    (delta_peaks_maybe["peak_max"] > np.quantile(delta_peaks_maybe["peak_max"], 0.98))]



sensen = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'sen1')]
dbldbl = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'dbl8')]
dsds = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'ds')]


deltachr1 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr1')]
deltachr2 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr2')]
deltachr3 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr3')]

sensen1 = sensen[(sensen['chro'] == 'chr1')]
sensen2 = sensen[(sensen['chro'] == 'chr2')]
sensen3 = sensen[(sensen['chro'] == 'chr3')]

dbldbl1 = dbldbl[(dbldbl['chro'] == 'chr1')]
dbldbl2 = dbldbl[(dbldbl['chro'] == 'chr2')]
dbldbl3 = dbldbl[(dbldbl['chro'] == 'chr3')]

dsds1 = dsds[(dsds['chro'] == 'chr1')]
dsds2 = dsds[(dsds['chro'] == 'chr2')]
dsds3 = dsds[(dsds['chro'] == 'chr3')]

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
    return chrI, chrII, chrIII, all_data



    

ewtchrI, ewtchrII, ewtchrIII, ewt = Create_df('Jo_wt_YE.e1.f-w300.count.csv', 'Jo_wt_YE.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)
#sen1chrI, sen1chrII, sen1chrIII 
esen1chr1, esen1chr2, esen1chr3, esen = Create_df ('Jo_sen1d_YE.e1.f-w300.count.csv', 'Jo_sen1d_YE.e1.r-w300.count.csv', 'RZ269-RhArepr-e.e1.f-w300.count.csv', 'RZ269-RhArepr-e.e1.r-w300.count.csv', 5)
edbl8chr1, edbl8chr2, edbl8chr3, edbl = Create_df('Jo_dbl8_YE.e1.f-w300.count.csv', 'Jo_dbl8_YE.e1.r-w300.count.csv', 'RZ271-RhArepr-e.e1.f-w300.count.csv', 'RZ271-RhArepr-e.e1.r-w300.count.csv', 5)
edschr1, edschr2, edschr3, eds = Create_df('Jo_sen1_dbl8_YE.e1.f-w300.count.csv', 'Jo_sen1_dbl8_YE.e1.r-w300.count.csv', 'RZ273-RhArepr-e.e1.f-w300.count.csv', 'RZ273-RhArepr-e.e1.r-w300.count.csv', 5)


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

#%%
def Chromosome_plot (data1, data2, data3, data4, featurex, centro, telo, genee, c, con,data5, data6, datat, odata1, odata2, odata3, odata4):
    ff, (ax1, ax2, axb, ax3, ax5) = plt.subplots(5,1, sharex=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['smoo_right_forks'], color ='black', alpha=0.8)
    ax1.plot(data2['pos'], data2['smoo_right_forks'], color ='steelblue', alpha=0.8)
    ax1.plot(data3['pos'], data3['smoo_right_forks'], color ='orange', alpha=0.8)
    ax1.plot(data4['pos'], data4['smoo_right_forks'], color ='tomato', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('rightward forks (e)')
    
    #ax2.set_title('sen1')
    ax2.plot(data1['pos'], data1['smoo_d_usage'], color ='black', alpha=0.8)
    ax2.plot(data2['pos'], data2['smoo_d_usage'], color ='steelblue', alpha=0.8)
    ax2.plot(data3['pos'], data3['smoo_d_usage'], color ='orange', alpha=0.8)
    ax2.plot(data4['pos'], data4['smoo_d_usage'], color ='tomato', alpha=0.8)
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('drightward forks (e)')
   # ax2.set_ylim(0,200)

    ax3.set_ylabel('filter dif') 
   # ax2.plot(data1['pos'], data1['d_usage_norm'], color ='black', alpha=0.8)
    axb.plot(odata2['pos'], odata2['d_usage_norm'], color ='steelblue', alpha=0.8)
    axb.plot(odata3['pos'], odata3['d_usage_norm'], color ='orange', alpha=0.8)
    axb.plot(odata4['pos'], odata4['d_usage_norm'], color ='tomato', alpha=0.8)
    

    ax5.set_xlabel('Chromosome position')

    for be in odata1.itertuples(index=False, name = None):
        if be[23] > 0.23:
            ax1.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            ax3.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
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
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                
    for sug in data5.itertuples(index=False, name=None):
        ax3.axvspan(sug[3],sug[4],0.1,0.9,color="steelblue",alpha=0.9)
        
    for dug in data6.itertuples(index=False, name=None):
        ax3.axvspan(dug[3],dug[4],0.1,0.9,color="orange",alpha=0.9)                

    for dsug in datat.itertuples(index=False, name=None):
        ax3.axvspan(dsug[3],dsug[4],0.1,0.9,color="red",alpha=0.9)
                
    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(ewtchrI, esen1chr1, edbl8chr1, edschr1, gene1, centro1, telo1, feat1, 'I', con1, sensen1, dbldbl1, dsds1, owtchrI, osen1chr1, odbl8chr1, odschr1)
chromosome2 = Chromosome_plot(ewtchrII, esen1chr2, edbl8chr2, edschr2, gene2, centro2, telo2, feat2, 'II', con2, sensen2, dbldbl2, dsds2, owtchrII, osen1chr2, odbl8chr2, odschr2)
chromosome3 = Chromosome_plot(ewtchrIII, esen1chr3, edbl8chr3, edschr3, gene3, centro3, telo3, feat3, 'III', con3, sensen3, dbldbl3, dsds3, owtchrIII, osen1chr3, odbl8chr3, odschr3)




#%%
#%%
#%%

#so now i want to see which replicon has delta bias in it ]
# lets adapt the function so that we know which replicon it occurs in
#if overlap with stall gene : class as overlap 
#if in same replicom: same 
#if not: independent 



fullgenesandstalls = [newccontrol, ggenes]
gggenes = pd.concat(fullgenesandstalls)
gene1 = gggenes[(gggenes['chro'] == 'chr1')]
gene2 = gggenes[(gggenes['chro'] == 'chr2')]
gene3 = gggenes[(gggenes['chro'] == 'chr3')]
    

#we do this for both control and stall 

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

#this is with rob only origins 
#rrepliconchr1 = create_replicon_dataframe(owt1n)
#rrepliconchr2 = create_replicon_dataframe(owt2n)
#rrepliconchr3 = create_replicon_dataframe(owt3n)

###tremeber to change origins to rob only 
repliconchr1 = create_replicon_dataframe(owt1n)
repliconchr2 = create_replicon_dataframe(owt2n)
repliconchr3 = create_replicon_dataframe(owt3n)

#assign replicon id to delta peak, so which replicon does it occur in 
#so nans are occuring here 
#this is because some delta peaks extend over more than one replicon wtf 
#check that wiht plotting please 
#also please note that i only need to do this function once, there is no need to do it 
#again with the control data 
def assign_id_to_delta(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start_pos'] >= replicon_row['start'] and
                gene_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(gene_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return gene_df

deltachr1 = assign_id_to_delta(repliconchr1, deltachr1)
deltachr2 = assign_id_to_delta(repliconchr2, deltachr2)
deltachr3 = assign_id_to_delta(repliconchr3, deltachr3)




####not this one 

#def assign_id_to_delta(replicon_df, delta_df):
    matching_replicon_ids = []
    gene_list = []

    for _, delta_row in delta_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (delta_row['start_pos'] >= replicon_row['start'] and
                delta_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(delta_row['peak_no'])
                gene_list.append(replicon_row['id'])
            elif (delta_row['start_pos'] >= replicon_row['start'] and
                delta_row['end_pos'] >= replicon_row['stop']):
                matching_replicon_ids.append(delta_row['peak_no'])
                gene_list.append(replicon_row['id'])
                
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    delta_df = delta_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return delta_df

#deltachr1 = assign_id_to_delta(repliconchr1, deltachr1)
#deltachr2 = assign_id_to_delta(repliconchr2, deltachr2)
#deltachr3 = assign_id_to_delta(repliconchr3, deltachr3)


#deltachr1.to_csv('/Users/patricfernandez/Documents/python/ALLdelta_usage_chr1.csv', index=False)
#deltachr2.to_csv('/Users/patricfernandez/Documents/python/ALLdelta_usage_chr2.csv', index=False)
#deltachr3.to_csv('/Users/patricfernandez/Documents/python/ALLdelta_usage_chr3.csv', index=False)  



def Chromosome_plot (featurex, centro, telo, genee, c, con, odata1, rep, delta):
    ff, (ax5) = plt.subplots(1,1, sharex=True)
    
    
    for bed in rep.itertuples(index=False, name = None):
        ax5.axvspan(bed[0], bed[0], 0.0, 1, color = 'black', alpha =0.5)
        ax5.axvspan(bed[1], bed[1], 0.0, 1, color = 'black', alpha =0.5)
        ax5.annotate(bed[2], xy = [bed[0],0.45])  
     
    for bedd in delta.itertuples(index=False, name = None):
         ax5.axvspan(bedd[3], bedd[4], 0.0, 1, color = 'purple', alpha =0.2)
         ax5.annotate(bed[2], xy = [bed[1],0.45]) 
         
    
    
    for fe in featurex.itertuples(index=False, name=None):
        ax5.annotate(fe[0], xy = [fe[2],0.35])  
       # ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[4] == 'reverse':
            if fe[9] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[9] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[9] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                 
        
        elif fe[4] == 'forward':
            if fe[9] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[9] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[9] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                ax5.set_ylabel('Gene annotations')
    
    for c in centro.itertuples(index=False, name=None):
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
    
    for t in telo.itertuples(index=False, name=None):
            ax5.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
    
    
    
    return ff

chromosome1 = Chromosome_plot(gene1, centro1, telo1, feat1, 'I', con1, owtchrI, repliconchr1, deltachr1)
chromosome2 = Chromosome_plot(gene2, centro2, telo2, feat2, 'I', con2, owtchrII, repliconchr2, deltachr2)
chromosome3 = Chromosome_plot(gene3, centro3, telo3, feat3, 'I', con3, owtchrIII, repliconchr3, deltachr3)


    


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



#gene1.to_csv('/Users/patricfernandez/Documents/python/delta_usage_chr1.csv', index=False)
#gene2.to_csv('/Users/patricfernandez/Documents/python/delta_usage_chr2.csv', index=False)
#gene3.to_csv('/Users/patricfernandez/Documents/python/delta_usage_chr3.csv', index=False)  

###I think i may have to manually deal with these 
#so they were semi automated, let's import them back in 

def Find_fix(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)


    return genes

#gene1 = Find_fix("delta_usage_chr1.csv")


#gene2 = Find_fix("delta_usage_chr2.csv")


#gene3 = Find_fix("delta_usage_chr3.csv")




deltachr1 = Find_fix("ALLdelta_usage_chr1.csv")


deltachr2 = Find_fix("ALLdelta_usage_chr2.csv")


deltachr3 = Find_fix("ALLdelta_usage_chr3.csv")




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



agahgah = [deltachr1, deltachr2, deltachr3]
deltabias_and_ids = pd.concat(agahgah)   

#deltabias_and_ids.to_csv('/Users/patricfernandez/Documents/python/fixingtheworldstartswithmeandyourmum/deltabias_df_6thjan)with_ids', index=False)

#deltabiasdeltabias = Find_fix("deltabias_df_6thjan)with_ids.csv")


#deltachr1 = deltabiasdeltabias[(deltabiasdeltabias['DBIAS_chro'] == 'chr1')]
#deltachr2 = deltabiasdeltabias[(deltabiasdeltabias['chro'] == 'chr2')]
#deltachr3 = deltabiasdeltabias[(deltabiasdeltabias['chro'] == 'chr3')]

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


#this was my dumb way to fix a bug which doesnt seem to be happening anymore 6/1/24

#gene1['delta_containing'] = np.where(~gene1['source'].isna(), 'YES', gene1['delta_containing'])
#gene2['delta_containing'] = np.where(~gene2['source'].isna(), 'YES', gene2['delta_containing'])
#gene3['delta_containing'] = np.where(~gene3['source'].isna(), 'YES', gene3['delta_containing'])



allthegoodg = [gene1, gene2, gene3]
allthestalls = pd.concat(allthegoodg)
    


####please note that there are more than we initially 
#thought becayse some loci are delta bias +ve in more than one loci 
#so maybe for the simple task of asking who is delta +ve (loci/replicon wise) 
#we can just remove duplicates and keep first (31/12/23)

##1/1/24 this has been fixed by dropping duplicate ID


controlaaaaal = allthestalls[(allthestalls['type'] == 'mRNA')]
controlaaaaal.drop_duplicates(subset=['ID'], keep='first', inplace=True)



sen1stall = allthestalls[(allthestalls['genotype'] == 'sen1D')]
sen1stall.drop_duplicates(subset=['ID'], keep='first', inplace=True)

dbl8stall = allthestalls[(allthestalls['genotype'] == 'dbl8D')]
dbl8stall.drop_duplicates(subset=['ID'], keep='first', inplace=True)

doublestall = allthestalls[(allthestalls['genotype'] == 'sen1dbl8DD_unique')]
doublestall.drop_duplicates(subset=['ID'], keep='first', inplace=True)

condition = (
    (doublestall['genotype'] == 'sen1dbl8DD_unique') & 
    ((doublestall['coding_strand'] == 'reverse') & (doublestall['stalled_fork'] == 'leftward') | 
     (doublestall['coding_strand'] == 'forward') & (doublestall['stalled_fork'] == 'rightward'))
)

# Apply the condition to filter the DataFrame
HTstall = doublestall.loc[condition]


condition2 = (
    (doublestall['genotype'] == 'sen1dbl8DD_unique') & 
    ((doublestall['coding_strand'] == 'reverse') & (doublestall['stalled_fork'] == 'rightward') | 
     (doublestall['coding_strand'] == 'forward') & (doublestall['stalled_fork'] == 'leftward'))
)

# Apply the condition to filter the DataFrame
HHstall = doublestall.loc[condition2]





#%%
import matplotlib.pyplot as plt
import seaborn as sns

count = 0

for x in controlaaaaal.itertuples():
    if x[15] == 'YES':
        count += 1
        print(count)
        


####new plot, for stalls and controls, are they delta +ve? 
groups = ['sen1D stall', 'dbl8D stall ', 'HT sddd', 'HH sddd', 'control loci']
values1 = [8, 16, 17, 31, 60]
values2 = [1, 14, 21, 39, 129]

groups = ['sen1D stall', 'dbl8D stall ', 'HT sddd', 'HH sddd']
values1 = [8/9*100, 16/30*100, 17/38*100, 31/70*100, 60/189*100]
values2 = [1/9*100, 14/30*100, 21/38*100, 39/70*100, 129/189*100]


fig, ax = plt.subplots()

# Stacked bar chart
ax.bar(groups, values1, color = "#1D5B79", alpha = 0.8, label = "delta bias positive")
ax.bar(groups, values2, bottom = values1, color = "black", alpha = 0.8,  label = "delta bias negative")

for bar in ax.patches:
  ax.text(bar.get_x() + bar.get_width() / 2,
          bar.get_height() / 2 + bar.get_y(),
          round(bar.get_height()), ha = 'center',
          color = 'w', weight = 'bold', size = 10)

total_values = np.add(values1, values2)

# Total values labels
for i, total in enumerate(total_values):
  ax.text(i, total + 0.5, round(total),
          ha = 'center', weight = 'bold', color = 'black')


#ax.legend()
plt.show()

#%%
###note that we updates to deal with controls and genes in the same earleir function 
#and now for the controls that have things in it 

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

repliconchr1c = create_replicon_dataframe(owt1n)
repliconchr2c = create_replicon_dataframe(owt2n)
repliconchr3c = create_replicon_dataframe(owt3n)



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

con1 = assign_id_to_stalls(repliconchr1c, con1)
con2 = assign_id_to_stalls(repliconchr2c, con2)
con3 = assign_id_to_stalls(repliconchr3c, con3)
    


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
  #  matching_gene_ids = []  # New list to store gene IDs
    
    # Iterate over rows in gene_df
    for _, gene_row in gene_df.iterrows():
        # Check the condition using vectorized operations
        mask = (replicon_df['start'] <= gene_row['start']) & (gene_row['end'] <= replicon_df['stop'])
        print(mask)
        
        # Extract matching replicon IDs and add to the list
        matching_replicon_ids.extend(replicon_df.loc[mask, 'id'].tolist())
        

    # Create a DataFrame with matching IDs and 'YES' for 'stall_containing'
    temp = pd.DataFrame({'id': matching_replicon_ids, 'stall_containing': 'YES'})
    
    
    # Merge the temporary DataFrame with replicon_df
    replicon_df = replicon_df.merge(temp, on='id', how='left')

    # Fill NaN values in 'stall_containing' with 'NO'
    replicon_df['stall_containing'] = replicon_df['stall_containing'].fillna('NO')
    
       
    
    for _, delta_row in delta.iterrows():
        
        mask2 = (replicon_df['start'] <= delta_row['start_pos']) & (delta_row['end_pos'] <= replicon_df['stop'])
        mathcing_rep_and_d_id.extend(replicon_df.loc[mask2, 'id'].tolist())
        
    temp2 = pd.DataFrame({'id': mathcing_rep_and_d_id, 'delta_containing': 'YES'})
    replicon_df = replicon_df.merge(temp2, on='id', how='left')
    replicon_df['delta_containing'] = replicon_df['delta_containing'].fillna('NO')
    
    df_no_duplicates = replicon_df.drop_duplicates(subset='id', keep='first')

        

    return df_no_duplicates

# Example usage
repliconchr1c = building_layers(repliconchr1c, con1, deltachr1)
repliconchr2c = building_layers(repliconchr2c, con2, deltachr2)
repliconchr3c = building_layers(repliconchr3c, con3, deltachr3)

'''
def assign_id_to_delta(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start_pos'] >= replicon_row['start'] and
                gene_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(gene_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return gene_df

deltachr1 = assign_id_to_delta(repliconchr1c, deltachr1)
deltachr2 = assign_id_to_delta(repliconchr2c, deltachr2)
deltachr3 = assign_id_to_delta(repliconchr3c, deltachr3)
'''


con1 = con1.merge(repliconchr1c, on='id', how='left')
con2 = con2.merge(repliconchr2c, on='id', how='left')
con3 = con3.merge(repliconchr3c, on='id', how='left')



con1 = con1.merge(deltachr1, on='id', how='left')
assert all(con1.loc[con1['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"
con2 = con2.merge(deltachr2, on='id', how='left')
assert all(con2.loc[con2['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"
con3 = con3.merge(deltachr3, on='id', how='left')
assert all(con3.loc[con3['delta_containing'] == 'NO', 'peak_no'].isna()), "Assertion failed: col2 should be NaN when col1 is 'NO'"


allthegdg = [con1, con2, con3]
allthecon = pd.concat(allthegdg)
    

count = 0


for x in allthecon.itertuples():
    if x[12] == 'YES':
        count += 1
        print(count)
        


####new plot, for stalls only
groups = ['sen1D stall', 'dbl8D stall ', 'HT sddd', 'HH sddd', 'control']
values1 = [8, 15, 16, 28, 49]
values2 = [1, 15, 22, 42, 140]


groups = ['sen1D stall', 'dbl8D stall ', 'HT sddd', 'HH sddd', 'control']
values1 = [89, 50, 42, 40, 10]
values2 = [11, 50, 58, 60, 90]


fig, ax = plt.subplots()

# Stacked bar chart
ax.bar(groups, values1, color = "#1D5B79", alpha = 0.8, label = "delta bias positive")
ax.bar(groups, values2, bottom = values1, color = "black", alpha = 0.8,  label = "delta bias negative")

for bar in ax.patches:
  ax.text(bar.get_x() + bar.get_width() / 2,
          bar.get_height() / 2 + bar.get_y(),
          round(bar.get_height()), ha = 'center',
          color = 'w', weight = 'bold', size = 10)

total_values = np.add(values1, values2)

# Total values labels
for i, total in enumerate(total_values):
  ax.text(i, total + 0.5, round(total),
          ha = 'center', weight = 'bold', color = 'black')


#ax.legend()
plt.show()



#%%



#this is now the df without duplicates. 

allthegoodg_nodups = [sen1stall, dbl8stall, doublestall, controlaaaaal]
allthestalls_nodups = pd.concat(allthegoodg_nodups)

  



#I guess now what is important is to find out for each stall that has delta bias, 
#what genotype is it id's from
#note that sen1d had the most delta bias peaks 
def score(file, frame1, frame2, frame3):
    for i in range(len(file)):
        #this is replicon start and stop
        tempstart = file.iloc[i]["start_y"]
        tempend = file.iloc[i]['stop']
        tempchro = file.iloc[i]["chro_x"]
       # rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['pos'] >= tempstart) & (frame1['pos'] <= tempend) & (frame1['chro'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['pos'] >= tempstart) & (frame2['pos'] <= tempend) & (frame2['chro'] == tempchro)]
        tempsubsetm = frame3.loc[(frame3['pos'] >= tempstart) & (frame3['pos'] <= tempend) & (frame3['chro'] == tempchro)]
        

        file.loc[file.index[i], 'thing'] = tempsubset['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_sen'] = tempsubset['d_usage_norm'].sum() / len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_dbl'] = tempsubsett['d_usage_norm'].sum() / len(tempsubsett)
        
        file.loc[file.index[i], 'thingm'] = tempsubsetm['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_ds'] = tempsubsetm['d_usage_norm'].sum() / len(tempsubsetm)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb_sen'] = np.log2(file.loc[file.index[i], 'epb_sen'])
        file.loc[file.index[i], 'log2_epb_dbl'] = np.log2(file.loc[file.index[i], 'epb_dbl'])
        file.loc[file.index[i], 'log2_epb_ds'] = np.log2(file.loc[file.index[i], 'epb_ds'])

        # Calculate the length of the region and add log2_length column
   #     length = tempend - tempstart
    #    file.loc[file.index[i], 'log2_length'] = np.log2(length)
     #   file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

allthestalls_extra = score(allthestalls_nodups, sencut, dblcut, dscut)
#allthecon_extra = score(allthecon, sencut, dblcut, dscut)

#tmb = [allthestalls_extra, allthecon_extra]
#try_me_bitch = pd.concat(tmb)
#try_me_bitch['type'] = try_me_bitch['type'].fillna('STALL')



allthestalls_extra['type'] = allthestalls_extra['type'].fillna('STALL')
allthestalls_extra['source'].fillna('non-delta', inplace=True)

#

melted_df = allthestalls_extra.melt(id_vars=['ID', 'gene', 'start_x', 'end', 'chro_x',
                                       'coding_strand', 'stalled_fork', 'genotype', 'id',
                                       'start_y', 'stop', 'stall_containing', 'delta_containing', 
                                       "peak_no", "chro_y", "start_pos", "end_pos", "peak_area",
                                       "peak_max", 'thing', 'epb_sen', 'thingt',
                                       'epb_dbl', 'thingm', 'epb_ds', 'cpc', 'type', 'source'],
                               var_name='source_e',
                               value_vars= ['log2_epb_sen', "log2_epb_dbl", "log2_epb_ds"],
                               value_name='log2_epb')


#melted_df = try_me_bitch.melt(id_vars=['ID', 'gene', 'start_x', 'end', 'chro_x',
 #                                      'coding_strand', 'stalled_fork', 'genotype', 'id',
  #                                     'start_y', 'stop', 'stall_containing', 'delta_containing', 
   #                                    "peak_no", "chro_y", "start_pos", "end_pos", "peak_area",
    #                                   "peak_max", 'thing', 'epb_sen', 'thingt',
     #                                  'epb_dbl', 'thingm', 'epb_ds', 'cpc', 'type'],
      #                         var_name='source',
       #                        value_name='log2_epb')

# Extract the suffix from the 'source' column
melted_df['source_e'] = melted_df['source_e'].str.split('_').str[-1]

# Print the result
print(melted_df)
melted_df['type'] = melted_df['type'].astype(str)
melted_df['log2_epb'] = pd.to_numeric(melted_df['log2_epb'], errors='coerce')



custom_palette = ["steelblue", "orange", 'tomato' ]  

sns.set_palette(custom_palette)



sns.swarmplot(
    data=melted_df, x="type", y="log2_epb",
    hue="source"
)

plt.figure(figsize=(7, 11)) 
sns.violinplot(
    data=melted_df, x="genotype", y="log2_epb",
    hue="source_e", alpha =0.8
)
plt.ylim(-12, 1)

sns.violinplot(
    data=allthestalls_extra, x="type", y="log2_epb_sen",
    hue="delta_containing"
)

fx= sns.pointplot(data=melted_df, x="type", y="log2_epb", hue="source_e", dodge=True)
fx.show()

fig =sns.lmplot(
    data=melted_df, x="peak_max", y="peak_area",
    hue="type", col="source_e", height=4,
)


#%%
#now i want to rank by intensity of delta bias, so let's make a scoring metric 
#normalise area by 


melted_df['log2_peakmax'] = np.log2(melted_df['peak_max'])
melted_df['log2_peakarea'] = np.log2(melted_df['peak_area'])
sns.lmplot(
    data=melted_df, x="log2_peakmax", y="log2_peakarea",
    hue="type", col="source_e", height=4,
)

melted_df['peak_size'] = (melted_df['end_pos'] - melted_df['start_pos'])/300

melted_df['d_peak_strength'] = (melted_df['peak_area'] / melted_df['peak_size'])*melted_df['peak_max']
melted_df['log2_d_peak_strength'] = np.log2(melted_df['d_peak_strength'])




####so really for this one what i need to do is correlate the strenegth of the stall vs the delta bias
##we can do this simply with epistallall and deltastall all dfs from python stall call. just import 

sns.lmplot(
    data=melted_df, x="log2_d_peak_strength", y="log2_epb",
    hue="type", col="source_e", height=4,
)

#why are there duplicates. i need to find out

#%%
def Find_pos(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)


    return genes

epsilonforyourplot = Find_pos("1stjanepsilonstalls.csv")
#deltaforyourplot = Find_pos("1stjanDELTAstalls.csv")
#deltaforyourplot = Find_fix("4thjan_deltastallswithrepliconids.csv")
deltaforyourplot = Find_pos('5thjan_deltastallswithrepliconids.csv')


deltaforyourplot['peak_size'] = (deltaforyourplot['end_pos'] - deltaforyourplot['start_pos'])/300

deltaforyourplot['peak_strength'] = (deltaforyourplot['peak_area'] / deltaforyourplot['peak_size'])*deltaforyourplot['peak_max']
deltaforyourplot['log2_peak_strength'] = np.log2(deltaforyourplot['peak_strength'])

deltaforyourplot['log2_peakmax'] = np.log2(deltaforyourplot['peak_max'])
deltaforyourplot['log2_peakarea'] = np.log2(deltaforyourplot['peak_area'])



epsilonforyourplot.rename(columns={'log2_peak_strength': 'STALLlog2_peak_strength'}, inplace=True)
deltaforyourplot.rename(columns={'log2_peak_strength': 'STALLlog2_peak_strength'}, inplace=True)
allthestalls_extra.rename(columns={'ID': 'Systematic ID'}, inplace=True)


#so now i want to merge on SID but only take STALLlog2_peak_strength


anexperiment = pd.merge(allthestalls_extra, epsilonforyourplot[['Systematic ID', 'STALLlog2_peak_strength']], on='Systematic ID', how='left')
anexperiment.drop_duplicates(subset=['Systematic ID'], keep='first', inplace=True)

#an_experiment_stall_only =anexperiment[(anexperiment[''])]


an_experiment_stall_only = anexperiment[(anexperiment['type'] == 'STALL')]
an_experiment_stall_only_nona = an_experiment_stall_only.dropna(subset=['STALLlog2_peak_strength'])


an_experiment_stall_only['peak_size'] = (an_experiment_stall_only['end_pos'] - an_experiment_stall_only['start_pos'])/300

an_experiment_stall_only['d_peak_strength'] = (an_experiment_stall_only['peak_area'] / an_experiment_stall_only['peak_size'])*an_experiment_stall_only['peak_max']
an_experiment_stall_only['log2_d_peak_strength'] = np.log2(an_experiment_stall_only['d_peak_strength'])



sns.lmplot(
    data=an_experiment_stall_only, x="log2_d_peak_strength", y="STALLlog2_peak_strength",
    hue="type", col="source", height=4,
)



#%%

#why this hasn't worked is because systematic ID does not match 
#so what i should do is see if i can go by replicon ID
#epsilonforyourplot needs a replicon ID

#####JOANNA PAY ATTENTION
####LOOK HERE. NEED TO DO WITH deltaforyourplot TO GET DELTA VERSION
epsilonforyourplot1 = epsilonforyourplot[(epsilonforyourplot['chro'] == 'chr1')]
epsilonforyourplot2 = epsilonforyourplot[(epsilonforyourplot['chro'] == 'chr2')]
epsilonforyourplot3 = epsilonforyourplot[(epsilonforyourplot['chro'] == 'chr3')]
                                           

def assign_id_to_delta(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start_pos'] >= replicon_row['start'] and
                gene_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(gene_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return gene_df

epsilonforyourplot11 = assign_id_to_delta(repliconchr1, epsilonforyourplot1)
epsilonforyourplot2 = assign_id_to_delta(repliconchr2, epsilonforyourplot2)
epsilonforyourplot3 = assign_id_to_delta(repliconchr3, epsilonforyourplot3)


epsilonforyourplotewlh = [epsilonforyourplot11, epsilonforyourplot2, epsilonforyourplot3]
epsilonforyourplot = pd.concat(epsilonforyourplotewlh)

#epsilonforyourplot.rename(columns={'log2_peak_strength': 'STALLlog2_peak_strength'}, inplace=True)
#allthestalls_extra.rename(columns={'ID': 'Systematic ID'}, inplace=True)


#so now i want to merge on SID but only take STALLlog2_peak_strength


anexperiment = pd.merge(allthestalls_extra, epsilonforyourplot[['id', 'STALLlog2_peak_strength']], on='id', how='left')
anexperiment.drop_duplicates(subset=['Systematic ID'], keep='first', inplace=True)

#an_experiment_stall_only =anexperiment[(anexperiment[''])]


an_experiment_stall_only = anexperiment[(anexperiment['type'] == 'STALL')]
an_experiment_stall_only_nona = an_experiment_stall_only.dropna(subset=['STALLlog2_peak_strength'])


an_experiment_stall_only['peak_size'] = (an_experiment_stall_only['end_pos'] - an_experiment_stall_only['start_pos'])/300

an_experiment_stall_only['d_peak_strength'] = (an_experiment_stall_only['peak_area'] / an_experiment_stall_only['peak_size'])*an_experiment_stall_only['peak_max']
an_experiment_stall_only['log2_d_peak_strength'] = np.log2(an_experiment_stall_only['d_peak_strength'])




an_experiment_stall_only_DELTAYES = an_experiment_stall_only[(an_experiment_stall_only['delta_containing'] == 'YES')]



#so there are only 58 stall genes that have deltabias
#thus why i can't plot more :)
an_experiment_stall_only['log2_peakmax'] = np.log2(an_experiment_stall_only['peak_max'])
an_experiment_stall_only['log2_peakarea'] = np.log2(an_experiment_stall_only['peak_area'])
sns.lmplot(
    data=an_experiment_stall_only, x="log2_d_peak_strength", y="log2_peakmax",
    hue="genotype", height=4
)





        
an_experiment_stall_onlySEN = an_experiment_stall_only[(an_experiment_stall_only['genotype'] == 'dbl8D')]
an_experiment_stall_onlySEN= an_experiment_stall_onlySEN.dropna(subset=['peak_no'])



ffeatNO = ffeat[(ffeat['DStop3'] == 'NO')]
ffeatNO = ffeatNO.dropna(subset=['Rank_DS'])
ffeatNO

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(an_experiment_stall_onlySEN['log2_d_peak_strength']).reshape((-1, 1))
yy = np.array(an_experiment_stall_onlySEN['STALLlog2_peak_strength'])

# Add a constant term to the independent variable
xx = sm.add_constant(xx)

# Create a linear regression model
model = sm.OLS(yy, xx).fit()

# Get the p-values for the coefficients
p_values = model.pvalues

# Print the p-values
print('P-values:', p_values)

# Print the equation of the regression line
slope = model.params[1]
intercept = model.params[0]
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')



sns.lmplot(
    data=an_experiment_stall_only_DELTAYES, x="log2_d_peak_strength", y="STALLlog2_peak_strength",
    col="genotype", height=4,
)



###i need to melt an_experiment
#want to plot log2 (delta_bias) versus STALLlog2_peak_strength. hue yes or no for delta bias

experiment_melted_df = an_experiment_stall_only.melt(id_vars=['Systematic ID', 'chro_x', 'start_x', 'end', 'coding_strand', 'cpc',
       'type', 'gene', 'stalled_fork', 'genotype', 'id', 'start_y', 'stop',
       'stall_containing', 'delta_containing', 'peak_no', 'chro_y',
       'start_pos', 'end_pos', 'peak_area', 'peak_max', 'thing', 'epb_sen',
       'thingt', 'epb_dbl', 'thingm', 'epb_ds', 'STALLlog2_peak_strength', 'peak_size',
       'd_peak_strength', 'log2_d_peak_strength'],
                               var_name='source',
                               value_vars= ['log2_epb_sen', "log2_epb_dbl", "log2_epb_ds"],
                               value_name='log2_epb')


#melted_df = try_me_bitch.melt(id_vars=['ID', 'gene', 'start_x', 'end', 'chro_x',
 #                                      'coding_strand', 'stalled_fork', 'genotype', 'id',
  #                                     'start_y', 'stop', 'stall_containing', 'delta_containing', 
   #                                    "peak_no", "chro_y", "start_pos", "end_pos", "peak_area",
    #                                   "peak_max", 'thing', 'epb_sen', 'thingt',
     #                                  'epb_dbl', 'thingm', 'epb_ds', 'cpc', 'type'],
      #                         var_name='source',
       #                        value_name='log2_epb')

# Extract the suffix from the 'source' column
experiment_melted_df['source'] = experiment_melted_df['source'].str.split('_').str[-1]

# Print the result
print(experiment_melted_df)
experiment_melted_df['type'] = experiment_melted_df['type'].astype(str)
experiment_melted_df['log2_epb'] = pd.to_numeric(experiment_melted_df['log2_epb'], errors='coerce')



custom_palette = ["steelblue", "orange", 'tomato' ]  
custom_palette = ["#1D5B79", "#EF6262" ]
custom_palette = ["steelblue", "black"] 
sns.set_palette(custom_palette)


sns.kdeplot(
    data=experiment_melted_df, x="STALLlog2_peak_strength", y="log2_epb",
    hue="source"
)

sns.violinplot(
    data=experiment_melted_df, x="delta_containing", y="STALLlog2_peak_strength",
    hue="delta_containing", alpha =0.8
)
#plt.ylim(, -2)

sns.kdeplot(
    data=experiment_melted_df,  x="STALLlog2_peak_strength",
    hue="delta_containing", alpha =0.8
)


fx= sns.pointplot(data=experiment_melted_df, x="delta_containing", y="STALLlog2_peak_strength", hue="genotype", dodge=True)
fx.show()




fig =sns.lmplot(
    data=experiment_melted_df, x="STALLlog2_peak_strength", y="log2_epb",
    hue="delta_containing", col="genotype", height=4,
)



count = 0

for x in an_experiment_stall_only.itertuples():
    if x[15] == 'YES':
        count += 1
        print(count)

#%%

#####JOANNA PAY ATTENTION
####LOOK HERE. NEED TO DO WITH deltaforyourplot TO GET DELTA VERSION



#%%
#so now i want to merge on SID but only take STALLlog2_peak_strength

    deltaforyourplot.loc[deltaforyourplot['chro'] == "I", 'chro'] = 'chr1'
    deltaforyourplot.loc[deltaforyourplot['chro'] == "II", 'chro'] = 'chr2'
    deltaforyourplot.loc[deltaforyourplot['chro'] == "III", 'chro'] = 'chr3'

deltaforyourplot
def score(file, frame1, frame2, frame3):
    for i in range(len(file)):
        #this is replicon start and stop
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]['End position']
        tempchro = file.iloc[i]["chro"]
       # rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['pos'] >= tempstart) & (frame1['pos'] <= tempend) & (frame1['chro'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['pos'] >= tempstart) & (frame2['pos'] <= tempend) & (frame2['chro'] == tempchro)]
        tempsubsetm = frame3.loc[(frame3['pos'] >= tempstart) & (frame3['pos'] <= tempend) & (frame3['chro'] == tempchro)]
        

        file.loc[file.index[i], 'thing'] = tempsubset['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_sen'] = tempsubset['d_usage_norm'].sum() / len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_dbl'] = tempsubsett['d_usage_norm'].sum() / len(tempsubsett)
        
        file.loc[file.index[i], 'thingm'] = tempsubsetm['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_ds'] = tempsubsetm['d_usage_norm'].sum() / len(tempsubsetm)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb_sen'] = np.log2(file.loc[file.index[i], 'epb_sen'])
        file.loc[file.index[i], 'log2_epb_dbl'] = np.log2(file.loc[file.index[i], 'epb_dbl'])
        file.loc[file.index[i], 'log2_epb_ds'] = np.log2(file.loc[file.index[i], 'epb_ds'])

        # Calculate the length of the region and add log2_length column
   #     length = tempend - tempstart
    #    file.loc[file.index[i], 'log2_length'] = np.log2(length)
     #   file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

allthestalls_extra = score(deltaforyourplot, sencut, dblcut, dscut)


deltaforyourplot1 = deltaforyourplot[(deltaforyourplot['chro'] == 'chr1')]
deltaforyourplot2 = deltaforyourplot[(deltaforyourplot['chro'] == 'chr2')]
deltaforyourplot3 = deltaforyourplot[(deltaforyourplot['chro'] == 'chr3')]

#Ithink you need to merege deltachr1 with repliconchr1 so that i can have replicon 
#start and stop sites
#or i can manually go through and check they're in the right replicon 

###please note that deltachr1,2,3  is delta bias                                         
'''
def assign_id_to_delta(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start_pos'] >= replicon_row['start'] and
                gene_row['end_pos'] <= replicon_row['stop']):
                matching_replicon_ids.append(gene_row['peak_no'])
                gene_list.append(replicon_row['id'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'peak_no': matching_replicon_ids, 'id': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='peak_no', how='left')
    print(id_and_ID)
    

    return gene_df

deltaforyourplot11 = assign_id_to_delta(repliconchr1, deltaforyourplot1)
deltaforyourplot2 = assign_id_to_delta(repliconchr2, deltaforyourplot2)
deltaforyourplot3 = assign_id_to_delta(repliconchr3, deltaforyourplot3)


deltaforyourplotwlh = [deltaforyourplot11, deltaforyourplot2, deltaforyourplot3]
deltaforyourplot = pd.concat(deltaforyourplotwlh)

#epsilonforyourplot.rename(columns={'log2_peak_strength': 'STALLlog2_peak_strength'}, inplace=True)
#allthestalls_extra.rename(columns={'ID': 'Systematic ID'}, inplace=True)

'''


columns_to_rename = {'source': 'DBIAS_source', 'peak_no': 'DBIAS_peak_no', 'chro': 'DBIAS_chro',
                     'start_pos': 'DBIAS_start_pos', 'end_pos': 'DBIAS_end_pos',
                     'peak_area': 'DBIAS_peak_area', 'peak_max': 'DBIAS_peak_max',
                     'Unnamed: 8': 'DBIAS_Unnamed_8'}

# Rename the specified columns
deltachr1.rename(columns=columns_to_rename, inplace=True)
deltachr2.rename(columns=columns_to_rename, inplace=True)
deltachr3.rename(columns=columns_to_rename, inplace=True)


# Display the modified DataFrame
print(deltachr1)



oooooooooooowe = deltaforyourplot1.merge(deltachr1, on='id', how='left')
oooooooooooowe2 = deltaforyourplot2.merge(deltachr2, on='id', how='left')
oooooooooooowe3 = deltaforyourplot3.merge(deltachr3, on='id', how='left')


owe = [oooooooooooowe, oooooooooooowe2, oooooooooooowe3]
deltastalldbias = pd.concat(owe,ignore_index=True)
#deltastalldbias_nodupy = deltastalldbias.drop_duplicates(subset=deltastalldbias.columns)
#this means it appears to be duplicates but they are actaually all unique across the board

deltastalldbias_nodupy = deltastalldbias.drop_duplicates(subset=['Systematic ID', 'genotype'], keep='first')

#####alternatively to drop na
deltastalldbias_nodupy['DBIAS_source'] = deltastalldbias_nodupy['DBIAS_source'].fillna('dbias_NO')


dan_experiment_dpositiveonly = deltastalldbias.dropna(subset=['DBIAS_source'])
##61 delta stalls have delta bias
#######
dan_experiment_stall_only = deltastalldbias_nodupy.copy()

dan_experiment_stall_only['DBIASpeak_size'] = (dan_experiment_stall_only['DBIAS_end_pos'] - dan_experiment_stall_only['DBIAS_start_pos'])/300

dan_experiment_stall_only['d_peak_strength'] = (dan_experiment_stall_only['DBIAS_peak_area'] / dan_experiment_stall_only['DBIASpeak_size'])*dan_experiment_stall_only['DBIAS_peak_max']
dan_experiment_stall_only['log2_d_peak_strength'] = np.log2(dan_experiment_stall_only['d_peak_strength'])


#an_experiment_stall_only_DELTAYES = an_experiment_stall_only[(anexperiment['delta_containing'] == 'YES')]



#so there are only 58 stall genes that have deltabias
#thus why i can't plot more :)

sns.lmplot(
    data=dan_experiment_stall_only, x="log2_d_peak_strength", y="STALLlog2_peak_strength",
    hue="genotype", height=4,
)




an_experiment_stall_onlySEN = an_experiment_stall_only[(an_experiment_stall_only['genotype'] == 'dbl8D')]
an_experiment_stall_onlySEN= an_experiment_stall_onlySEN.dropna(subset=['peak_no'])



ffeatNO = ffeat[(ffeat['DStop3'] == 'NO')]
ffeatNO = ffeatNO.dropna(subset=['Rank_DS'])
ffeatNO

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(an_experiment_stall_onlySEN['log2_d_peak_strength']).reshape((-1, 1))
yy = np.array(an_experiment_stall_onlySEN['STALLlog2_peak_strength'])

# Add a constant term to the independent variable
xx = sm.add_constant(xx)

# Create a linear regression model
model = sm.OLS(yy, xx).fit()

# Get the p-values for the coefficients
p_values = model.pvalues

# Print the p-values
print('P-values:', p_values)

# Print the equation of the regression line
slope = model.params[1]
intercept = model.params[0]
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')







sns.lmplot(
    data=dan_experiment_stall_only, x="log2_d_peak_strength", y="log2_peakarea",
    hue="genotype", height=4,
)




custom_palette = ["orange", "steelblue", 'tomato' ]  
custom_palette = ["#1D5B79", "#EF6262" ]
custom_palette = ["steelblue", "black"] 
sns.set_palette(custom_palette)



###i need to melt an_experiment
#want to plot log2 (delta_bias) versus STALLlog2_peak_strength. hue yes or no for delta bias

dexperiment_melted_df = dan_experiment_stall_only.melt(id_vars=['condition', 'genotype', 'peak_no', 'chro', 'start_pos', 'end_pos',
       'peak_area', 'peak_max', 'Systematic ID', 'could be', 'Start position',
       'End position', 'Strand', 'peak_midpoint', 'average_forks',
       'stalled_fork', 'orientation', 'id', 'Unnamed: 9', 'peak_size',
       'peak_strength', 'STALLlog2_peak_strength', 'log2_peakmax',
       'log2_peakarea', 'thing', 'epb_sen', 'thingt', 'epb_dbl', 'thingm',
       'epb_ds', 'DBIAS_source',
       'DBIAS_peak_no', 'DBIAS_chro', 'DBIAS_start_pos', 'DBIAS_end_pos',
       'DBIAS_peak_area', 'DBIAS_peak_max',
       'DBIASpeak_size', 'd_peak_strength', 'log2_d_peak_strength'],
                               var_name= 'source',
                               value_vars= ['log2_epb_sen', "log2_epb_dbl", "log2_epb_ds"],
                               value_name='log2_epb')


#melted_df = try_me_bitch.melt(id_vars=['ID', 'gene', 'start_x', 'end', 'chro_x',
 #                                      'coding_strand', 'stalled_fork', 'genotype', 'id',
  #                                     'start_y', 'stop', 'stall_containing', 'delta_containing', 
   #                                    "peak_no", "chro_y", "start_pos", "end_pos", "peak_area",
    #                                   "peak_max", 'thing', 'epb_sen', 'thingt',
     #                                  'epb_dbl', 'thingm', 'epb_ds', 'cpc', 'type'],
      #                         var_name='source',
       #                        value_name='log2_epb')

# Extract the suffix from the 'source' column
dexperiment_melted_df['source'] = dexperiment_melted_df['source'].str.split('_').str[-1]

# Print the result
print(dexperiment_melted_df)
#dexperiment_melted_df['type'] = dexperiment_melted_df['type'].astype(str)
dexperiment_melted_df['log2_epb'] = pd.to_numeric(dexperiment_melted_df['log2_epb'], errors='coerce')



custom_palette = ["steelblue", "orange", 'tomato' ]  
custom_palette = ["#1D5B79", "#EF6262" ]
custom_palette = ["steelblue", "black"] 
sns.set_palette(custom_palette)


    dexperiment_melted_df.loc[dexperiment_melted_df['DBIAS_source'] == "sen1", 'DBIAS_source'] = 'd_STALL'
    dexperiment_melted_df.loc[dexperiment_melted_df['DBIAS_source'] == "dbl8", 'DBIAS_source'] = 'd_STALL'
    dexperiment_melted_df.loc[dexperiment_melted_df['DBIAS_source'] == "ds", 'DBIAS_source'] = 'd_STALL'


dexperiment_melted_df["DBIAS_source"] = pd.Categorical(dexperiment_melted_df["DBIAS_source"],
                                                 ordered = True,
                                                 categories = ["d_STALL", "dbias_NO"])

sns.violinplot(
    data=dexperiment_melted_df, x="DBIAS_source", y="STALLlog2_peak_strength",
    hue="DBIAS_source", alpha =0.8
)

plt.ylim(-12.2, -6.4)




sns.kdeplot(
    data=experiment_melted_df, x="STALLlog2_peak_strength", y="log2_epb",
    hue="source"
)

sns.violinplot(
    data=dexperiment_melted_df, x="DBIAS_source", y="STALLlog2_peak_strength",
    hue="DBIAS_source", alpha =0.8
)
plt.ylim(-12, -2)

sns.violinplot(
    data=dexperiment_melted_df, x="DBIAS_source", y="log2_epb",
    hue="DBIAS_source", alpha =0.8
)

sns.kdeplot(
    data=experiment_melted_df,  x="STALLlog2_peak_strength",
    hue="DBIAS_source", alpha =0.8
)


fx= sns.pointplot(data=experiment_melted_df, x="delta_containing", y="STALLlog2_peak_strength", hue="genotype", dodge=True)
fx.show()




fig =sns.lmplot(
    data=experiment_melted_df, x="STALLlog2_peak_strength", y="log2_epb",
    hue="delta_containing", col="genotype", height=4,
)



count = 0

for x in dan_experiment_stall_only.itertuples():
    if x[34] == 'dbias_NO':
        count += 1
        print(count)
        
        
        
        
dan_experiment_dpositiveonly





dddexperiment_melted_df = dan_experiment_dpositiveonly.melt(id_vars=['condition', 'genotype', 'peak_no', 'chro', 'start_pos', 'end_pos',
       'peak_area', 'peak_max', 'Systematic ID', 'could be', 'Start position',
       'End position', 'Strand', 'peak_midpoint', 'average_forks',
       'stalled_fork', 'orientation', 'id', 'Unnamed: 9', 'peak_size',
       'peak_strength', 'STALLlog2_peak_strength', 'log2_peakmax',
       'log2_peakarea', 'thing', 'epb_sen', 'thingt', 'epb_dbl', 'thingm',
       'epb_ds', 'DBIAS_source',
       'DBIAS_peak_no', 'DBIAS_chro', 'DBIAS_start_pos', 'DBIAS_end_pos',
       'DBIAS_peak_area', 'DBIAS_peak_max'],
                               var_name= 'source',
                               value_vars= ['log2_epb_sen', "log2_epb_dbl", "log2_epb_ds"],
                               value_name='log2_epb')



sns.violinplot(
    data=dddexperiment_melted_df, x="genotype", y="log2_epb",
    hue="source", alpha =0.8, legend=False
)
plt.ylim(-12, 1)

dddexperiment_melted_df['type'] = 'delta_stall'

sns.violinplot(
    data=dddexperiment_melted_df, x="condition", y="log2_epb",
    hue="genotype", alpha =0.8
)

    


dddexperiment_melted_df['source_e'] = dddexperiment_melted_df['source'].str.split('_').str[-1]

plt.ylim(-12, 0)

blahblahblah = pd.concat([dddexperiment_melted_df, melted_df], axis=0, ignore_index=True)



sns.stripplot(
    data=blahblahblah, x="type", y="log2_epb", hue="source_e",
    dodge=True, alpha=.2, legend=False,
)
sns.pointplot(
    data=blahblahblah, x="type", y="log2_epb", alpha =0.8, hue='source_e', dodge=True
)
plt.ylim(7, 3)
plt.figure(figsize=(7, 3)) 
        

melted_df['type'].unique()
#%%





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
esen1chr1 = fill_missing_positions(osen1chr1, "chr1", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
esen1chr2 = fill_missing_positions(osen1chr2, "chr2", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
esen1chr3 = fill_missing_positions(osen1chr3, "chr3", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])





edbl8chr1 = fill_missing_positions(odbl8chr1, "chr1", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
edbl8chr2 = fill_missing_positions(odbl8chr2, "chr2", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
edbl8chr3 = fill_missing_positions(odbl8chr3, "chr3", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])

edschr1 = fill_missing_positions(odschr1, "chr1", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
edschr2 = fill_missing_positions(odschr2, "chr2", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])
edschr3 = fill_missing_positions(odschr3, "chr3", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin', 'd_usage_norm'])


ewtchrI = fill_missing_positions(owtchrI, "chr1", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin'])
ewtchrII = fill_missing_positions(owtchrII, "chr2", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin'])
ewtchrIII = fill_missing_positions(owtchrIII, "chr3", pos_column='pos', columns_to_fill_zero=['df_count', 'dr_count', 'ef_count', 'er_count',
       'norm_df', 'norm_dr', 'norm_ef', 'norm_er', 'ratio_delta_f',
       'ratio_delta_r', 'ratio_epsilon_f', 'ratio_epsilon_r', 'd_usage',
       'right_forks', 'smoo_ratio_d_f', 'smoo_ratio_d_r', 'smoo_ratio_e_f',
       'smoo_ratio_e_r', 'smoo_d_usage', 'smoo_right_forks', 'differentials',
       'origin'])







sen1l = [esen1chr1, esen1chr2, esen1chr3]
esen1 = pd.concat(sen1l,ignore_index=True)

dbl8l = [edbl8chr1, edbl8chr2, edbl8chr3]
edbl8 = pd.concat(dbl8l ,ignore_index=True)


dsl = [edschr1, edschr2, edschr3]
eds = pd.concat(dsl,ignore_index=True)


ewtchrI = ewt[(ewt['chro'] == 'chr1')]
ewtchrII = ewt[(ewt['chro'] == 'chr2')]
ewtchrIII = ewt[(ewt['chro'] == 'chr3')]

esen1chr1 = esen1[(esen1['chro'] == 'chr1')]
esen1chr2 = esen1[(esen1['chro'] == 'chr2')]
esen1chr3 = esen1[(esen1['chro'] == 'chr3')]

edbl8chr1 = edbl8[(edbl8['chro'] == 'chr1')]
edbl8chr2 = edbl8[(edbl8['chro'] == 'chr2')]
edbl8chr3 = edbl8[(edbl8['chro'] == 'chr3')]

edschr1 = eds[(eds['chro'] == 'chr1')]
edschr2 = eds[(eds['chro'] == 'chr2')]
edschr3 = eds[(eds['chro'] == 'chr3')]


an_experiment_stall_only
10 genotype
11 start_y
12 stop

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

an_experiment_stall_only

an_experiment_stall_only_for = an_experiment_stall_only[(an_experiment_stall_only['coding_strand'] == 'forward')]
an_experiment_stall_only_rev = an_experiment_stall_only[(an_experiment_stall_only['coding_strand'] == 'reverse')]

an_experiment_stall_only['sbinpos'] = an_experiment_stall_only['start_x']/300
an_experiment_stall_only['sbinpos'] = an_experiment_stall_only['sbinpos'].astype(int)
an_experiment_stall_only['sbinpos'] = an_experiment_stall_only['sbinpos']*300 +150

an_experiment_stall_only['ebinpos'] = an_experiment_stall_only['end']/300
an_experiment_stall_only['ebinpos'] = an_experiment_stall_only['ebinpos'].astype(int)
an_experiment_stall_only['ebinpos'] = an_experiment_stall_only['ebinpos']*300 +150



def asign_orientation(file):
    for i in range(len(file)):
        fork = file.iloc[i]["stalled_fork"]
        gene = file.iloc[i]["coding_strand"]    
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
an_experiment_stall_only = asign_orientation(an_experiment_stall_only)


HTstallsall = an_experiment_stall_only[(an_experiment_stall_only['orientation'] == 'Head-Tail')]
HHstallsall = an_experiment_stall_only[(an_experiment_stall_only['orientation'] == 'Head-Head')]




10 genotype
12 start_y
13 stop

#for wt , forward coding strand genes only 
#TSS
def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[10] == p:
            if g[2] == 'chr1':
                x = chr1.index[chr1['pos'] == g[36]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[36]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[36]].tolist()
                xx.append(x3)
                
    return xx

senforxx = Start(an_experiment_stall_only, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
#senrevxx = Start(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxx = Start(an_experiment_stall_only, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxR = Start(an_experiment_stall_only, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')
#dbl8revxx = Start(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxx = Start(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxx = Start(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')

HTstallsall, HHstallsall
###
senforxx = Start(HTstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
senforxxHH = Start(HHstallsall, esen1chr1, esen1chr2, esen1chr3, 'sen1D')


dbl8forxx = Start(HTstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dbl8forxxHH = Start(HHstallsall, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')

dsforxxR = Start(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')
dsforxxRHH = Start(HHstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')




#TES
def End(file, chr1, chr2, chr3, p):
    xxe =[]
    for ge in file.itertuples():
        if ge[10] == p:
            if ge[2] == 'chr1':
                xe = chr1.index[chr1['pos'] == ge[37]].tolist()
                xxe.append(xe) 

            if ge[2] == 'chr2':
                xe2 = chr2.index[chr2['pos'] == ge[37]].tolist()
                xxe.append(xe2)

            if ge[2] == 'chr3':
                xe3 = chr3.index[chr3['pos'] == ge[37]].tolist()
                xxe.append(xe3)
    return xxe

senforxxe = End(an_experiment_stall_only, esen1chr1, esen1chr2, esen1chr3, 'sen1D')
#senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(an_experiment_stall_only, edbl8chr1, edbl8chr2, edbl8chr3, 'dbl8D')
dsforxxeR = End(an_experiment_stall_only, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')
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

dsforxxeR = End(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')
dsforxxeRHH = End(HTstallsall, edschr1, edschr2, edschr3, 'sen1dbl8DD_unique')






def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'d_usage_norm'].tolist()
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
        yy = frame.loc[z[0]-10:z[0],'d_usage_norm'].tolist()
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
        yy = frame.loc[ze[0]:ze[0]+10,'d_usage_norm'].tolist()
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
        num = 50/length
        print (num*length)

        te = np.arange(0,50,num, dtype = int)

        expand = np.zeros(50)

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


fuckkkkkk, (ax4, ax6, ax8) = plt.subplots(3,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(dssR, cmap = 'coolwarm', ax=ax4, cbar=True, cbar_ax=cbar_ax)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(dsbR, cmap = 'coolwarm', ax=ax6,cbar=True, cbar_ax=cbar_ax)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(dsdsR, cmap = 'coolwarm', ax=ax8,cbar=True, cbar_ax=cbar_ax)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,300])
ax8.set_xticklabels(['peak start', 'peak end'])




#DBL8 heatmap !!!!!
fxb, (( ax3, ax5,ax7)) = plt.subplots(3,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])

ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax3,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax7,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,11, 31, 42])
ax7.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])

#SEN1 heatmap bitches 
fx, (( ax3, ax5,ax7)) = plt.subplots(3,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.81, .3, .03, .24])

ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax, vmin= -0.1, vmax=0.1)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')


ax7.set_xticks([0,11, 31, 42])
ax7.set_xticklabels(['-3Kb','peak start', 'peak end', '+3kb'])




#%%

#so this is for from before, when i worked without the HH and HT 

    fig = plt.figure(constrained_layout=True) #create a figure
    gs = fig.add_gridspec(177, 29)  #sets the figure width and height (can only refer to these as integers)




# Now add as many subplots using the gridspace (gs) to position them precisely
# the first position is the height (ie here from 1 to 3 on the grid)
# the second position is the width (i.e here from 1 to 10 on the grid)
    ax1 = fig.add_subplot(gs[10:19, 5:10])
    sns.heatmap(sas, cmap = 'coolwarm', ax=ax1, cbar=False, vmin=-0.1, vmax=0.1)
    ax1.set_xticks([])
    ax1.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:19, 12:17]) 
    sns.heatmap(sab, cmap = 'coolwarm', ax=ax2,  cbar=False, vmin=-0.1, vmax=0.1)
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:19, 19:24]) 
    sns.heatmap(sads, cmap = 'coolwarm', ax=ax3, cbar=False,vmin=-0.1, vmax=0.1)
    ax3.set_xticks([])
    ax3.set_yticks([])

   # ax3.set_title('difference right forks',loc = 'left', pad = 0) 
   
   
   
   
    
    ax4 = fig.add_subplot(gs[24:54, 5:10])
    sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax4,  cbar=False, vmin=-0.1,vmax=0.1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    ax5 = fig.add_subplot(gs[24:54, 12:17])
    sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5,  cbar=False,vmin=-0.1, vmax=0.1)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    ax6 = fig.add_subplot(gs[24:54, 19:24])
    sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax6,  cbar=False,vmin=-0.1, vmax=0.1)
    ax6.set_xticks([])
    ax6.set_yticks([])
    
    
    
    
    
    ax7 = fig.add_subplot(gs[59:167, 5:10])
    sns.heatmap(HTs, cmap = 'coolwarm', ax=ax7, cbar=False, vmin=-0.1,vmax=0.1)
    ax7.set_xticks([0,11, 61, 72])
    ax7.set_xticklabels(['-', 's', 'e', '+'])
    ax7.set_yticks([])
    
    
    ax8 = fig.add_subplot(gs[59:167, 12:17]) 
    sns.heatmap(HTb, cmap = 'coolwarm', ax=ax8,cbar=False,vmin=-0.1, vmax=0.1)
    ax8.set_xticks([0,11, 61, 72])
    ax8.set_xticklabels(['-', 's', 'e', '+'])
    ax8.set_yticks([])
    
    ax9 = fig.add_subplot(gs[59:167, 19:24]) 
    sns.heatmap(HTds, cmap = 'coolwarm', ax=ax9,cbar=False, vmin=-0.1,vmax=0.1)
    ax9.set_xticks([0,11, 61, 72])
    ax9.set_xticklabels(['-', 's', 'e', '+'])
    ax9.set_yticks([])
    
 
    
 
    
 
    
    
    #####GENE BODY ONLY
    fig = plt.figure(constrained_layout=True) #create a figure
    gs = fig.add_gridspec(177, 29)  #sets the figure width and height (can only refer to these as integers)


# Now add as many subplots using the gridspace (gs) to position them precisely
# the first position is the height (ie here from 1 to 3 on the grid)
# the second position is the width (i.e here from 1 to 10 on the grid)
    ax1 = fig.add_subplot(gs[10:19, 5:10])
    sns.heatmap(ass, cmap = 'coolwarm', ax=ax1, cbar=False, vmin=-0.1, vmax=0.1)
    ax1.set_xticks([])
    ax1.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:19, 12:17]) 
    sns.heatmap(ad, cmap = 'coolwarm', ax=ax2,  cbar=False, vmin=-0.1, vmax=0.1)
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:19, 19:24]) 
    sns.heatmap(ads, cmap = 'coolwarm', ax=ax3, cbar=False,vmin=-0.1, vmax=0.1)
    ax3.set_xticks([])
    ax3.set_yticks([])

   # ax3.set_title('difference right forks',loc = 'left', pad = 0) 
   
   
   
   
    
    ax4 = fig.add_subplot(gs[24:54, 5:10])
    sns.heatmap(dbss, cmap = 'coolwarm', ax=ax4,  cbar=False, vmin=-0.1,vmax=0.1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    ax5 = fig.add_subplot(gs[24:54, 12:17])
    sns.heatmap(dbd, cmap = 'coolwarm', ax=ax5,  cbar=False,vmin=-0.1, vmax=0.1)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    ax6 = fig.add_subplot(gs[24:54, 19:24])
    sns.heatmap(dbds, cmap = 'coolwarm', ax=ax6,  cbar=False,vmin=-0.1, vmax=0.1)
    ax6.set_xticks([])
    ax6.set_yticks([])
    
    
    
    
    
    ax7 = fig.add_subplot(gs[59:167, 5:10])
    sns.heatmap(dssR, cmap = 'coolwarm', ax=ax7, cbar=False, vmin=-0.1,vmax=0.1)
    ax7.set_xticks([0,50])
    ax7.set_xticklabels(['-','+'])
    ax7.set_yticks([])
    
    
    ax8 = fig.add_subplot(gs[59:167, 12:17]) 
    sns.heatmap(dsbR, cmap = 'coolwarm', ax=ax8,cbar=False,vmin=-0.1, vmax=0.1)
    ax8.set_xticks([0,50])
    ax8.set_xticklabels(['-', '+'])
    ax8.set_yticks([])
    
    ax9 = fig.add_subplot(gs[59:167, 19:24]) 
    sns.heatmap(dsdsR, cmap = 'coolwarm', ax=ax9,cbar=False, vmin=-0.1,vmax=0.1)
    ax9.set_xticks([0,50])
    ax9.set_xticklabels(['-','+'])
    ax9.set_yticks([])
    
    
    
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

oldgenesfor, oldgenesrev, oldggenes = Find("dbl8_stall_sites_direction.txt")
#controlfor, controlrev, ccontrol = Find('control_genes.txt')
oldgene1 = oldggenes[(oldggenes['chro'] == 'chr1')]
oldgene2 = oldggenes[(oldggenes['chro'] == 'chr2')]
oldgene3 = oldggenes[(oldggenes['chro'] == 'chr3')]



sensen = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'sen1')]
dbldbl = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'dbl8')]
dsds = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['source'] == 'ds')]


deltachr1 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr1')]
deltachr2 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr2')]
deltachr3 = delta_peaks_maybe_filt[(delta_peaks_maybe_filt['chro'] == 'chr3')]

sensen1 = sensen[(sensen['chro'] == 'chr1')]
sensen2 = sensen[(sensen['chro'] == 'chr2')]
sensen3 = sensen[(sensen['chro'] == 'chr3')]

dbldbl1 = dbldbl[(dbldbl['chro'] == 'chr1')]
dbldbl2 = dbldbl[(dbldbl['chro'] == 'chr2')]
dbldbl3 = dbldbl[(dbldbl['chro'] == 'chr3')]

dsds1 = dsds[(dsds['chro'] == 'chr1')]
dsds2 = dsds[(dsds['chro'] == 'chr2')]
dsds3 = dsds[(dsds['chro'] == 'chr3')]


def Chromosome_plot_for (data1,data2, data3, data4, featurex, data2c, data3c, data4c, peak1):
    ff, (ax1, ax2, ax3, ax5, ax6) = plt.subplots(5,1, sharex=True)

    ax1.set_ylabel('f counts (e)')
    
    #ax2.set_title('sen1')
    ax1.plot(data1['pos'], data1['d_usage'], color ='black', alpha=0.8, linewidth = 1)
    ax1.plot(data2['pos'], data2['d_usage'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax1.plot(data3['pos'], data3['d_usage'], color ='orange', alpha=0.8, linewidth = 1)
    ax1.plot(data4['pos'], data4['d_usage'], color ='tomato', alpha=0.8, linewidth = 1)
    
    ax2.plot(data1['pos'], data1['smoo_d_usage'], color ='black', alpha=0.8, linewidth = 1)
    ax2.plot(data2['pos'], data2['smoo_d_usage'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.plot(data3['pos'], data3['smoo_d_usage'], color ='orange', alpha=0.8, linewidth = 1)
    ax2.plot(data4['pos'], data4['smoo_d_usage'], color ='tomato', alpha=0.8, linewidth = 1)
    
   # ax3.plot(data1['pos'], data1['d_usage_norm'], color ='black', alpha=0.8, linewidth = 1)
    ax3.plot(data2['pos'], data2['d_usage_norm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax3.plot(data3['pos'], data3['d_usage_norm'], color ='orange', alpha=0.8, linewidth = 1)
    ax3.plot(data4['pos'], data4['d_usage_norm'], color ='tomato', alpha=0.8, linewidth = 1)
    
    
  #  ax4.plot(data1c['pos'], data1c['termination_change'], color ='black', alpha=0.8)
 #   ax4.plot(data2c['pos'], data2c['termination_change_smooth'], color ='steelblue', alpha=0.8,linewidth = 1)
  #  ax4.plot(data3c['pos'], data3c['termination_change_smooth'], color ='orange', alpha=0.8, linewidth = 1)
   # ax4.plot(data4c['pos'], data4c['termination_change_smooth'], color ='tomato', alpha=0.8, linewidth = 1)
    ax5.set_ylim(0.01,0.34)
  #  ax5.set_ylim(-0.01,0.2)
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
        if row['source'] == 'sen1':

            start_pos = row['start_pos']
            end_pos = row['end_pos']
    
        # Filter df2 based on the conditions
            filtered_data = data2[(data2['pos'] >= start_pos-300) & (data2['pos'] <= end_pos+300)]
         #   filtered_datadb = data3c[(data3c['pos'] >= start_pos-300) & (data3c['pos'] <= end_pos+300)]
          #  filtered_datads = data4c[(data4c['pos'] >= start_pos-300) & (data4c['pos'] <= end_pos+300)]
    
        # Plot the filtered data
            ax5.plot(filtered_data['pos'], filtered_data['d_usage_norm'],color ='steelblue',)
         #   ax5.plot(filtered_datadb['pos'], filtered_datadb['termination_change_smooth'],color ='orange',)
          #  ax5.plot(filtered_datads['pos'], filtered_datads['termination_change_smooth'],color ='tomato', label=f'Segment {index + 1}')


        if row['source'] == 'dbl8':
           start_posb = row['start_pos']
           end_posb = row['end_pos']
    
        # Filter df2 based on the conditions
           filtered_datab = data3[(data3['pos'] >= start_posb-300) & (data3['pos'] <= end_posb+300)]
    
        # Plot the filtered data
           ax5.plot(filtered_datab['pos'], filtered_datab['d_usage_norm'],color ='orange', label=f'Segment {index + 1}')

        if row['source'] == 'ds':
           start_posd = row['start_pos']
           end_posd = row['end_pos']
    
        # Filter df2 based on the conditions
           filtered_datads = data4[(data4['pos'] >= start_posd-300) & (data4['pos'] <= end_posd+300)]
    
        # Plot the filtered data
           ax5.plot(filtered_datads['pos'], filtered_datads['d_usage_norm'],color ='tomato', label=f'Segment {index + 1}')
        

   
   
    return ff

chr1 = Chromosome_plot_for (owtchrI, osen1chr1, odbl8chr1, odschr1, oldgene1,
                            sensen1,dbldbl1, dsds1, deltachr1)


#%%

sensenl = sensen[''].to_list()
dtcp2yay_aa_list = dtcp2yay_aa['id'].to_list()
dtcp3yay_aa_list = dtcp3yay_aa['id'].to_list()

gene1_list = gene1['id'].to_list()
gene2_list = gene2['id'].to_list()
gene3_list = gene3['id'].to_list()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

d1vg1 = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2 = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3 = intersection (dtcp3yay_aa_list, gene3_list)


def fucker(lst1, lst2):
    lst3 = [value for value in lst1 if not value in lst2]
    return lst3

d1vg1_NOT = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2_NOT = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3_NOT = intersection (dtcp3yay_aa_list, gene3_list)

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

#%%

repliconchr1['chro'] = 'chr1'
repliconchr2['chro'] = 'chr2'
repliconchr3['chro'] = 'chr3'

eatme = [repliconchr1, repliconchr2, repliconchr3]
replicondf = pd.concat(eatme ,ignore_index=True)
#now i'm going to try something quickly which is a little new, and an idea i had at 2am

#sum delta bias in each replicon using score function
#distribution and zscore if normal dist 
def score(file, frame1, frame2, frame3):
    for i in range(len(file)):
        #this is replicon start and stop
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]["End position"]
        tempchro = file.iloc[i]["chro"]
       # rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['pos'] >= tempstart) & (frame1['pos'] <= tempend) & (frame1['chro'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['pos'] >= tempstart) & (frame2['pos'] <= tempend) & (frame2['chro'] == tempchro)]
        tempsubsetm = frame3.loc[(frame3['pos'] >= tempstart) & (frame3['pos'] <= tempend) & (frame3['chro'] == tempchro)]
        

       # file.loc[file.index[i], 'thing'] = tempsubset['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_sen'] = tempsubset['d_usage_norm'].sum() / len(tempsubset)
        
      #  file.loc[file.index[i], 'thingt'] = tempsubsett['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_dbl'] = tempsubsett['d_usage_norm'].sum() / len(tempsubsett)
        
     #   file.loc[file.index[i], 'thingm'] = tempsubsetm['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_ds'] = tempsubsetm['d_usage_norm'].sum() / len(tempsubsetm)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb_sen'] = np.log2(file.loc[file.index[i], 'epb_sen'])
        file.loc[file.index[i], 'log2_epb_dbl'] = np.log2(file.loc[file.index[i], 'epb_dbl'])
        file.loc[file.index[i], 'log2_epb_ds'] = np.log2(file.loc[file.index[i], 'epb_ds'])

        # Calculate the length of the region and add log2_length column
   #     length = tempend - tempstart
    #    file.loc[file.index[i], 'log2_length'] = np.log2(length)
     #   file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

#replicondfyy = score(replicondf, sencut, dblcut, dscut)
ffeatyy = score(ffeat, sencut, dblcut, dscut)


def bastard(file, x):
    mean = file[x].mean()
    print(mean)
    sd = file[x].std()
    print(sd)
    threshold = 2
    temp = file[x]
    file['zscores'] = [(temp - mean) / sd for temp in temp]
    subset = file[(file['zscores']) > threshold]
    subset = subset.sort_values(x, ascending = False)
    return file, subset

replicondfyy, replicondfyySEN = bastard(replicondfyy, 'epb_sen')
replicondfyy, replicondfyyDBL = bastard(replicondfyy, 'epb_dbl')
replicondfyy, replicondfyyDS = bastard(replicondfyy, 'epb_ds')

ffeatyy, ffeatyySEN = bastard(ffeatyy, 'epb_sen')
ffeatyy, ffeatyyDBL = bastard(ffeatyy, 'epb_dbl')
ffeatyy, ffeatyyDS = bastard(ffeatyy, 'epb_ds')




meltyrepli = replicondfyy.melt(id_vars=['start', 'stop', 'id', 'stall_containing', 'delta_containing', 'chro',
        'log2_epb_sen', 'log2_epb_dbl','log2_epb_ds', 'zscores'],
                               var_name= 'source',
                               value_vars= ['epb_sen', 'epb_dbl', 'epb_ds',],
                               value_name='EPB')



sns.kdeplot(
   data=meltyrepli, x="EPB", hue="source",
   fill=True, common_norm=False,
   alpha=.5, linewidth=0,
)

#x_coords = [ 0.053405066731942744, 0.07234119865049557, value3]  # Replace with your actual x-axis values

annotations = [(0.053405066731942744, 'steelblue'), (0.07234119865049557, 'orange'), (0.06171829134069566, 'tomato')]  

# Annotate each point
for x_coord, color in annotations:
    plt.axvline(x_coord, color=color, linestyle='--')  # Vertical line


# Show the plot
plt.show()

#sen1
#mean 0.031347596591847895
#sd 0.011028735070047424
#dbl8
#0.039838707776592035
#0.016251245436951767
#ds
#0.03226025159006109
#0.014729019875317285



meltyfeat = ffeatyy.melt(id_vars=['Systematic ID', 'Gene name', 'Product description', 'Feature type',
       'Start position', 'End position', 'Chromosome', 'Strand', 'chro',
       'log2_epb_sen', 'log2_epb_dbl',
       'log2_epb_ds', 'zscores'],
                               var_name= 'source',
                               value_vars= ['epb_sen', 'epb_dbl', 'epb_ds',],
                               value_name='EPB')


annotations = [(0.06246107305746121, 'steelblue'), (0.07162792074789012, 'orange'), (0.06758004003664128, 'tomato')]
sns.kdeplot(data=meltyfeat, x="EPB", hue="source")
for x_coord, color in annotations:
    plt.axvline(x_coord, color=color, linestyle='--')  # Vertical line
plt.show()


sns.kdeplot(
   data=meltyfeat, x="EPB", hue="source",
   fill=True, common_norm=False,
   alpha=.5, linewidth=0,
)

0.02880805043405443
0.016826511311703388
0.032089037927077735
0.01976944141040619
0.029448627009987207
0.019065706513327035








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

oldgenesfor, oldgenesrev, oldggenes = Find("dbl8_stall_sites_direction.txt")
#controlfor, controlrev, ccontrol = Find('control_genes.txt')
oldgene1 = oldggenes[(oldggenes['chro'] == 'chr1')]
oldgene2 = oldggenes[(oldggenes['chro'] == 'chr2')]
oldgene3 = oldggenes[(oldggenes['chro'] == 'chr3')]






ffeatyy1 = ffeatyy[(ffeatyy['chro'] == 'chr1')]
ffeatyy2 = ffeatyy[(ffeatyy['chro'] == 'chr2')]
ffeatyy3 = ffeatyy[(ffeatyy['chro'] == 'chr3')]


ffeatyy1aa_list = ffeatyy1['Systematic ID'].to_list()
ffeatyy2aa_list = ffeatyy2['Systematic ID'].to_list()
ffeatyy3aa_list = ffeatyy3['Systematic ID'].to_list()


gene1_list = oldgene1['ID'].to_list()
gene2_list = oldgene2['ID'].to_list()
gene3_list = oldgene3['ID'].to_list()

ffeatyySEN_list = ffeatyySEN['Systematic ID'].to_list()
ffeatyyDBL_list = ffeatyyDBL['Systematic ID'].to_list()
ffeatyyDS_list = ffeatyyDS['Systematic ID'].to_list()

oldggenes_list = oldggenes['ID'].to_list()
newccontrol_list = newccontrol['ID'].to_list()





def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

d1vg1 = intersection (ffeatyy1aa_list, gene1_list)
d2vg2 = intersection (ffeatyy2aa_list, gene2_list)
d3vg3 = intersection (ffeatyy3aa_list, gene3_list)

okaysen = intersection(ffeatyySEN_list, newccontrol_list)
okaydbl = intersection(ffeatyyDBL_list, newccontrol_list)
okayds = intersection(ffeatyyDS_list, newccontrol_list)


#%%


def find_overlapping_intervals(df):
    # Group by "chro"
    grouped = df.groupby("chro")

    # Initialize an empty DataFrame to store overlapping intervals
    overlapping_intervals = pd.DataFrame(columns=df.columns)

    # Iterate through groups
    for group_name, group_df in grouped:
        # Iterate through pairs of rows in the group
        for i in range(len(group_df)):
            for j in range(i + 1, len(group_df)):
                row_i = group_df.iloc[i]
                row_j = group_df.iloc[j]

                # Check for overlapping intervals
                if (row_i["start_pos"] <= row_j["end_pos"] and row_i["end_pos"] >= row_j["start_pos"]) or \
                   (row_j["start_pos"] <= row_i["end_pos"] and row_j["end_pos"] >= row_i["start_pos"]):
                    # Add overlapping rows to the result DataFrame
                    overlapping_intervals = pd.concat([overlapping_intervals, pd.DataFrame([row_i, row_j])])

    return overlapping_intervals

trythisbadbitchshit = find_overlapping_intervals(delta_peaks_maybe_filt)

trythisbadbitchshity =trythisbadbitchshit.drop_duplicates(keep='first')





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

