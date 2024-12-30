#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:53:03 2024

@author: ixn004
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 13:03:57 2023

@author: dsw002
"""

import numpy as np
import pandas as pd
from fcsextract import fcsextract
import os
import sys
from scipy.stats import ttest_ind, sem
import matplotlib.pyplot as plt
import csv

# Function reads .fcs files, returns the data as a matrix and column names
def read_fcs(filename):
    x = fcsextract(filename)
    names = []
    for i in range(int(x[0][b'$PAR'])):
        names.append(x[0].get(b"$P%dS"%(i+1),x[0][b"$P%dN"%(i+1)]))
    data = np.asarray(x[1])
    return data,names

# Performs one-sided t-tests to determine if an increase occurred, and if a decrease occurred, 
# for two datasets.  Additionally, tests whether increases or decreases occurred, ignoring 
# statistical significance.  Returns a boolean vector indicating which data columns are significantly
# less in data2 than in data1, a boolean vector indicating same but greater in data2 than data1,
# and a boolean vector indicating whether the mean of data2 is greater (no statistical significance)
def check_increase(data1,data2):
    mean1 = np.mean(data1,0)
    mean2 = np.mean(data2,0)
    decrease = np.zeros(np.size(data1,1),dtype='bool')
    increase = np.zeros(np.size(data1,1),dtype='bool')
    alpha = 0.05
    for idx,a in enumerate(increase):
        decrease[idx] = ttest_ind(data1[:,idx],data2[:,idx],equal_var=False,alternative='greater')[1]<alpha
        increase[idx] = ttest_ind(data1[:,idx],data2[:,idx],equal_var=False,alternative='less')[1]<alpha
    comparison_result = mean1<mean2
    return decrease,increase,comparison_result

# All "analyze_XXX" functions are designed to be applied to particular data files.  They all
# indicate "good_idx", which are indices corresponding to measurable/important proteins (excludes
# proteins like cell cycle proteins and housekeeping markers).  They all indicate the data to be read,
# then iterate through datasets to compare whether protein expressions increased or decreased between
# cell cycle stages at various times/receptor stimulations.  "a", "b", and "c" refer to the outputs
# of "check_increase", respectively.

# Primary NK cells, but this can be ignored for the later, more updated functions to analyze
# primary NK cells.  This only analyzes them at t=0.
def analyze_primary_data():
    good_idx = []
    stage_strings = ['G1','S','G2']
    data,name_list = read_fcs('./data/primary/stage_gates/debarcode_PBS_G1.fcs') #note: These are CD56 dim cells
    a = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    for protein_idx,name in enumerate(name_list):
        for i in range(len(stage_strings)-1):
            data1,names = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i]+'.fcs')
            data2,names = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i+1]+'.fcs')
            data1[data1<0] = 0
            data2[data2<0] = 0
            a[protein_idx,i],b[protein_idx,i],c[protein_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# NKL receptor comparison (i.e. Ly49H vs. NKG2D vs. IgG).  Input argument is a string of the
# receptor name (i.e. "Ly49H")
def analyze_NKL_data(receptor):
    good_idx = [0,1,2,5,7,8,9,10,11,12,13,14,15,17,20,21,22,23,24,25,26,28,35,37,38,39,40,41,42,43,44,45,46,47]
    time_strings = ['8','32','64','128','256']
    stage_strings = ['G1','S','G2']
    data,name_list = read_fcs('./data/NKLs/stage_gates/BC1_CT604_IgG_8min_G1.fcs')
    a = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    for time_idx,time in enumerate(time_strings):
        for protein_idx,name in enumerate(name_list):
            for i in range(len(stage_strings)-1):
                data1,names = read_fcs('./data/NKLs/stage_gates/BC1_CT604_'+receptor+'_'+time+'min_'+stage_strings[i]+'.fcs')
                data2,names = read_fcs('./data/NKLs/stage_gates/BC1_CT604_'+receptor+'_'+time+'min_'+stage_strings[i+1]+'.fcs')
                data1[data1<0] = 0
                data2[data2<0] = 0
                a[protein_idx,time_idx,i],b[protein_idx,time_idx,i],c[protein_idx,time_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# Data from CellCycleTRACER
def analyze_tracer_data():
    good_idx = [2,3,4,5,6,7,8,9,10,12,13,15,16,17,25,26,27,28,29,35,44,45,46,47,48,49,51,52,53,54,55]
    time_strings = ['D5','B11','C11','D11','E11','G8','G9']
    stage_strings = ['G1','S','G2']
    data,name_list = read_fcs('./data/TRACER/BaseFileName_B11_G1.fcs')
    a = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    for time_idx,time in enumerate(time_strings):
        for protein_idx,name in enumerate(name_list):
            for i in range(len(stage_strings)-1):
                data1,names = read_fcs('./data/TRACER/BaseFileName_'+time+'_'+stage_strings[i]+'.fcs')
                data2,names = read_fcs('./data/TRACER/BaseFileName_'+time+'_'+stage_strings[i+1]+'.fcs')
                data1[data1<0] = 0
                data2[data2<0] = 0
                a[protein_idx,time_idx,i],b[protein_idx,time_idx,i],c[protein_idx,time_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# Data from a separate NKL experiment with t0, corresponding to NKG2D stimulation.
def analyze_NKL_t0():
    good_idx = [0,1,2,3,6,7,8,9,10,12,13,14,15,17,20,21,22,23,24,25,26,27,28,35,38,39,40,41,42,43,44,45,46]
    time_strings = ['0','8','16','32','64','128']
    stage_strings = ['G1','S','G2']
    data,name_list = read_fcs('./data/NKLs/stage_gates/gated_NKG2D_0min_G1.fcs')
    a = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    for time_idx,time in enumerate(time_strings):
        for protein_idx,name in enumerate(name_list):
            for i in range(len(stage_strings)-1):
                data1,names = read_fcs('./data/NKLs/stage_gates/gated_NKG2D_'+time+'min_'+stage_strings[i]+'.fcs')
                data2,names = read_fcs('./data/NKLs/stage_gates/gated_NKG2D_'+time+'min_'+stage_strings[i+1]+'.fcs')
                data1[data1<0] = 0
                data2[data2<0] = 0
                a[protein_idx,time_idx,i],b[protein_idx,time_idx,i],c[protein_idx,time_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# Primary NK data with IL2 priming
def analyze_primary_data_IL2():
    stage_strings = ['G1','S','G2']
    data,name_list = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_G1.fcs') #note: These are CD56 dim cells
    a = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(stage_strings)-1),dtype='bool')
    for protein_idx,name in enumerate(name_list):
        for i in range(len(stage_strings)-1):
            data1,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i]+'.fcs')
            data2,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i+1]+'.fcs')
            data1[data1<0] = 0
            data2[data2<0] = 0
            a[protein_idx,i],b[protein_idx,i],c[protein_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list

# A plotting function that will plot average protein expression for a given protein for a 
# given dataset.  Provide the protein name as a string, and the dataset as a integer.
# I can provide dataset option context if needed; forgoing for now because datasets rely on
# my local path.
def plot_protein(protein_name,dataset):
    stage_strings = ['G1','S','G2']
    mean = np.zeros(len(stage_strings))
    std = np.zeros(len(stage_strings))
    if dataset==1:
        time_strings = ['B11','C11','D11','E11','G8','G9']
        mean = np.zeros((len(stage_strings),len(time_strings)))
        std = np.zeros((len(stage_strings),len(time_strings)))
        data,name_list = read_fcs('./data/TRACER/BaseFileName_B11_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            for j in range(len(time_strings)):
                data,names = read_fcs('./data/TRACER/BaseFileName_'+time_strings[j]+'_'+stage_strings[i]+'.fcs')
                data[data<0] = 0
                mean[i,j] = np.mean(data[:,idx])
                std[i,j] = sem(data[:,idx])
    if dataset==2:
        stage_strings = ['G1','S','G2']
        data,name_list = read_fcs('./data/NKLs/stage_gates/BC1_CT604_IgG_8min_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            data = read_fcs('./data/NKLs/stage_gates/BC1_CT604_'+'IgG_0min_'+stage_strings[i]+'.fcs')
    elif dataset==3:
        data,name_list = read_fcs('./data/primary/stage_gates/debarcode_PBS_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            data,names = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i]+'.fcs')
            data[data<0] = 0
            mean[i] = np.mean(data[:,idx])
            std[i] = sem(data[:,idx])
    elif dataset==5:
        data,name_list = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            data,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i]+'.fcs')
            data[data<0] = 0
            mean[i] = np.mean(data[:,idx])
            std[i] = sem(data[:,idx])
    elif dataset==6:
        time_strings = ['2','4','8','16','32','64','256']
        mean = np.zeros((len(stage_strings),len(time_strings)+1))
        std = np.zeros((len(stage_strings),len(time_strings)+1))
        data,name_list = read_fcs('./data/primary/stage_gates/debarcode_NKG2D_2min_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            data,name_list = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i]+'.fcs')
            data[data<0] = 0
            mean[i,0] = np.mean(data[:,idx])
            std[i,0] = sem(data[:,idx])
            for j in range(len(time_strings)):
                data,names = read_fcs('./data/primary/stage_gates/debarcode_NKG2D_'+time_strings[j]+'min_'+stage_strings[i]+'.fcs')
                data[data<0] = 0
                mean[i,j+1] = np.mean(data[:,idx])
                std[i,j+1] = sem(data[:,idx])
    elif dataset==7:
        time_strings = ['2','4','8','16','32','64','256']
        mean = np.zeros((len(stage_strings),len(time_strings)+1))
        std = np.zeros((len(stage_strings),len(time_strings)+1))
        data,name_list = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_G1.fcs')
        idx = name_list.index(protein_name)
        for i in range(len(stage_strings)):
            data,name_list = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i]+'.fcs')
            data[data<0] = 0
            mean[i,0] = np.mean(data[:,idx])
            std[i,0] = sem(data[:,idx])
            for j in range(len(time_strings)):
                data,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_NKG2D_'+time_strings[j]+'min_'+stage_strings[i]+'.fcs')
                data[data<0] = 0
                mean[i,j+1] = np.mean(data[:,idx])
                std[i,j+1] = sem(data[:,idx])
    plt.figure()
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['font.size'] = 18
    plt.xlabel('Stage')
    plt.ylabel('Average Protein Expression')
    if dataset==1:
        plt.title(protein_name.decode('ascii')[6:])
    else:
        plt.title(protein_name.decode('ascii'))
    ind = np.arange(3)
    if len(mean.shape) == 1:
        plt.bar(['G1','S','G2'],mean,yerr=std)
    else:
        width = np.min(np.diff(ind))/(len(mean[0,:])+1)
        for i in range(len(mean[0,:])):
            plt.bar(ind+width*i,mean[:,i],yerr=std[:,i],width=width)
            plt.xticks(ind + ((len(mean[0,:])-1)/2)*width , ('G1', 'S', 'G2'))
    #legend = plt.legend(('0min','2min','4min','8min','16min','32min','64min','256min'),framealpha=1.0)
    #export_legend(legend,"my_legend")
    
# unneeded; checks the quotient of the minimum and maximum of a dataset
def minmax():
    data1,names = read_fcs('./data/TRACER/BaseFileName_B11_G1.fcs')
    data2,names = read_fcs('./data/TRACER/BaseFileName_B11_G2.fcs')
    data1[data1<=0] = 0.00000001
    data2[data2<=0] = 0.00000001
    mymin = np.min(data1,0)
    mymax = np.max(data2,0)
    mult = mymax/mymin
    return mult

# If a legend is needed as a separate file (due to space constraints or redundancy), use this
def export_legend(legend, filename="legend.png"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

# Primary NK cells, but many times (not just t=0)
def analyze_primary_time():
    good_idx = [0,1,2,3,5,6,7,8,9,10,12,13,14,15,17,18,21,22,23,24,27,28,29,36,38,39,40,41,42,43,44,45,46,47,48]
    stage_strings = ['G1','S','G2']
    time_strings = ['2','4','8','16','32','64','256']
    data,name_list = read_fcs('./data/primary/stage_gates/debarcode_NKG2D_2min_G1.fcs') #note: These are CD56 dim cells
    a = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    for protein_idx,name in enumerate(name_list):
        for i in range(len(stage_strings)-1):
            data1,names = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i]+'.fcs')
            data2,names = read_fcs('./data/primary/stage_gates/debarcode_PBS_'+stage_strings[i+1]+'.fcs')
            data1[data1<0] = 0
            data2[data2<0] = 0
            a[protein_idx,0,i],b[protein_idx,0,i],c[protein_idx,0,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    for time_idx,time in enumerate(time_strings):
        for protein_idx,name in enumerate(name_list):
            for i in range(len(stage_strings)-1):
                data1,names = read_fcs('./data/primary/stage_gates/debarcode_NKG2D_'+time+'min_'+stage_strings[i]+'.fcs')
                data2,names = read_fcs('./data/primary/stage_gates/debarcode_NKG2D_'+time+'min_'+stage_strings[i+1]+'.fcs')
                data1[data1<0] = 0
                data2[data2<0] = 0
                a[protein_idx,time_idx+1,i],b[protein_idx,time_idx+1,i],c[protein_idx,time_idx+1,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# Primary NK cells with IL2 priming, but many times (not just t=0)
def analyze_primary_IL2_time():
    good_idx = [0,1,2,3,5,6,7,8,9,10,12,13,14,15,17,18,21,22,23,24,27,28,29,36,38,39,40,41,42,43,44,45,46,47,48]
    stage_strings = ['G1','S','G2']
    time_strings = ['2','4','8','16','32','64','256']
    data,name_list = read_fcs('./data/primary/stage_gates/debarcode_IL2_NKG2D_2min_G1.fcs') #note: These are CD56 dim cells
    a = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    b = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings)+1,len(stage_strings)-1),dtype='bool')
    for protein_idx,name in enumerate(name_list):
        for i in range(len(stage_strings)-1):
            data1,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i]+'.fcs')
            data2,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_PBS_'+stage_strings[i+1]+'.fcs')
            data1[data1<0] = 0
            data2[data2<0] = 0
            a[protein_idx,0,i],b[protein_idx,0,i],c[protein_idx,0,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    for time_idx,time in enumerate(time_strings):
        for protein_idx,name in enumerate(name_list):
            for i in range(len(stage_strings)-1):
                data1,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_NKG2D_'+time+'min_'+stage_strings[i]+'.fcs')
                data2,names = read_fcs('./data/primary/stage_gates/debarcode_IL2_NKG2D_'+time+'min_'+stage_strings[i+1]+'.fcs')
                data1[data1<0] = 0
                data2[data2<0] = 0
                a[protein_idx,time_idx+1,i],b[protein_idx,time_idx+1,i],c[protein_idx,time_idx+1,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list,good_idx

# Which of the 9 possible conditions do our booleans generate??  Output used for figures.
# Remember, "a" is "has decrease", "b" is "has increase", and "c" is a raw comparison (no 
# statistical significance).  "c" is unused here.
def do_tests(a,b,c):
	#example: test 1 is that G1->S has a increase and S->G2 has a increase
    test1 = b[:,:,0] & b[:,:,1]
    # test 2 is that G1->S has a increase and S->G2 does not have a decrease or a increase
    test2 = b[:,:,0] & ~(a[:,:,1] | b[:,:,1])
    # etc
    test3 = b[:,:,0] & a[:,:,1]
    test4 = a[:,:,0] & b[:,:,1]
    test5 = a[:,:,0] & ~(a[:,:,1] | b[:,:,1])
    test6 = a[:,:,0] & a[:,:,1]
    test7 = ~(a[:,:,0] | b[:,:,0]) & b[:,:,1]
    test8 = ~(a[:,:,0] | b[:,:,0]) & ~(a[:,:,1] | b[:,:,1])
    test9 = ~(a[:,:,0] | b[:,:,0]) & a[:,:,1]
    return test1,test2,test3,test4,test5,test6,test7,test8,test9

# This takes the output in the next function and converts it to a list of proteins in each
# category.  For example, if two decreases only happened for 4 proteins, it will return those
# protein names.  Not super useful anymore.
def give_names(test,names,good_idx):
    name_str = []
    for name in names:
        name_str.append(name.decode('ascii'))
    pass_test = []
    for res in test:
        pass_test.append(np.asarray(name_str)[good_idx][res[good_idx].astype(bool)])
    return pass_test

# converts the results of all those booleans to a index for the figure result.  For instance,
# two succesive increases (test1) has index "0", an increase followed by "no change" (test2) 
# has index "1", etc.  "out" is the output, and it's given as an output argument.  The written
# file indicates which proteins belong to each category of trends.
def write_class_file(a,b,c,names,good_idx,filename):
    tests = do_tests(a,b,c)
    out = np.zeros((len(tests[0][:,0]),len(tests[0][0,:])))
    for idx,test in enumerate(tests):
        out[test] = idx
    with open(filename,'w',newline='\n') as file:
        wr = csv.writer(file)
        for idx,test in enumerate(tests):
            res = give_names(test.T,names,good_idx)
            wr.writerow(res)
    return out

# if you "run" this python script, this executes.  It gives examples of how to run
# particular functions too, and where my results came from.
if __name__=='__main__':

    # a,b,c,names,good_idx = analyze_tracer_data()
    # test_idx = write_class_file(a,b,c,names,good_idx,'tracer.csv')
    
    # a,b,c,names,good_idx = analyze_NKL_data('IgG')
    # test_idx = write_class_file(a,b,c,names,good_idx,'NKL_IgG.csv')
    
    # a,b,c,names,good_idx = analyze_NKL_t0()
    # test_idx = write_class_file(a,b,c,names,good_idx,'NKL_NKG2D.csv')
    
    a,b,c,names,good_idx = analyze_primary_time()
    test_idx = write_class_file(a,b,c,names,good_idx,'Primary_untreated.csv')
    
    a,b,c,names,good_idx = analyze_primary_IL2_time()
    test_idx = write_class_file(a,b,c,names,good_idx,'Primary_IL2.csv')
    
    # a,b,c,names,good_idx = analyze_NKL_data('Ly49H')
    # test_idx = write_class_file(a,b,c,names,good_idx,'NKL_Ly49H.csv')
    
    
    #plot_protein(b'CD45',6)
    plot_protein(b'CD45',7)
    #plot_protein(b'CD57',6)
    plot_protein(b'CD57',7)
    #plot_protein(b'pCreb',6)
    plot_protein(b'pCreb',7)
    #plot_protein(b'CD107a',6)
    plot_protein(b'CD107a',7)
    #plot_protein(b'pNFkB',6)
    plot_protein(b'pNFkB',7)
    #plot_protein(b'Ki67',6)
    plot_protein(b'Ki67',7)
    plot_protein(b'CD16',7)
    plot_protein(b'142Nd_pSRC',1)
    plot_protein(b'169Tm_CAPZA1',1)
    #plot_protein(b'166Er_pSTAT3',1)
    #plot_protein(b'145Nd_pMAPKAPK2',1)
    plot_protein(b'150Sm_pNFKB',1)
    #plot_protein(b'154Gd_pERK',1)
    #plot_protein(b'152Gd_pAMPK',1)
    
    
    
    