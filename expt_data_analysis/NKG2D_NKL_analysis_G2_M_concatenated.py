#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:28:07 2025

@author: ixn004
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 11:39:20 2024

@author: ixn004
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:04:12 2024
@author:written by Darren Wethington, Indrani Nayak 
Data:  NKG2D stimulated data on NKL cell
This is done for various time points ['8','32','64','128','256']

1.This code is to perform one sided  t-test between protein data between  (a) G1 and S, (b) S and G2 at a time point.
2. The result of 1 is summarize in the color map
3. This code is to plot the protein abundances in terms of geometric mean





#The raw( untransformed) data looks asymmetric and distributed log-normally Aug 26
#data processing steps
#step-1 : choose the negative raw data and replace them with a small value
#step-2 : transform them to log transformation

"""
import fcsextract
import numpy as np
from scipy.stats import ttest_ind
import os
import sys
from scipy.stats import ttest_ind, sem
import matplotlib.pyplot as plt
import csv
import pandas as pd



##read .facs fies as an matrix
def load_fcs_as_matrix(filename):
 """
  Loads flow cytometry data from an FCS file using fcsextract.

  Returns:
      tuple: A tuple containing two elements:
          - dict: A dictionary with metadata extracted from the FCS header.
          - list: A 2D list containing the event data (each sub-list represents an event).
  """
  # Use fcsextract to parse the FCS file
 metadata, data = fcsextract.fcsextract(filename)
 

 param_names = []
 for i in range(int(metadata[b'$PAR'])):
      param_names.append(metadata.get(b"$P%dS"%(i+1),metadata[b"$P%dN"%(i+1)]))
 data = np.array(data)
 
 #This part is added to prune the unwanted data ( The data had data1 and data 2 are same, for example barcode at column 53) Nov 5,24
 param_names=param_names[3:49]
 data=data[:,3:49]
 
 
 
 #do the transformation of the raw data needed
 data[data<0]=1e-11
 data=np.log(data)
 
 return data,param_names

 
##analyse CD107 high population with gating for G1, S and G2
def analyze_protein_variation(path):
  
    time_strings = ['8','32','64','128','256']
 
    stage_strings = ['G1','S','G2']
    
    data,name_list = load_fcs_as_matrix(path+'8min_UMAP_G1.fcs')
    data,param_names=load_fcs_as_matrix(filename)

    print(data.shape)#(20 x 54)
    print(len(param_names))#54
    
    #54x6x3 size for a, b and c
    a = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool') #(54 x 6 x3)
    b = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    c = np.zeros((len(name_list),len(time_strings),len(stage_strings)-1),dtype='bool')
    
    print(a.shape)
    
    for time_idx,time in enumerate(time_strings): #[0,5]
        for protein_idx,name in enumerate(name_list): #[0,53]
            for i in range(len(stage_strings)-1):#2 [0,1] between [G1,S], [S,G2]
                data1,names = load_fcs_as_matrix(path+time+'min_UMAP_'+stage_strings[i]+'.fcs')
                data2,names = load_fcs_as_matrix(path+time+'min_UMAP_'+stage_strings[i+1]+'.fcs')
                # Concatenate G2 and M data
                #print(f'new case------before: {data1.shape},{data2.shape}')  
                    
                if (i+1) == 2:  # If comparing G2 and M
                    print('here priniting)')
                    print(data2.shape)
                    print(load_fcs_as_matrix(path+time+'min_UMAP_M.fcs')[0].shape)
                    data2 = np.concatenate((data2, load_fcs_as_matrix(path+time+'min_UMAP_M.fcs')[0]), axis=0)
                    
                    print(data2.shape)
                #print(f'after: {data1.shape},{data2.shape}')   
             
                #if data is log transformed we do not need this step, if data is raw then keep this step
                # data1[data1<0] = 0 #( x 54)
                # data2[data2<0] = 0 # ( x54)
                #---------------------------------#
              
                a[protein_idx,time_idx,i],b[protein_idx,time_idx,i],c[protein_idx,time_idx,i] = check_increase(data1[:,protein_idx:protein_idx+1],data2[:,protein_idx:protein_idx+1])
    return a,b,c,name_list
# Remember, "a" is "has decrease", "b" is "has increase", and "c" is a raw comparison (no 


def check_increase(data1,data2): # this returns increase or decrease of specific protein between two cell stages
    # print(np.shape(data1))#(nx1)
    # print(np.shape(data2)) 
    mean1 = np.mean(data1,0)   #mean calculation for a column, among the rows
    mean2 = np.mean(data2,0)
    
    #print(mean1,mean2)
    
    decrease = np.zeros(np.size(data1,1),dtype='bool') #just a vector with 1 element
    increase = np.zeros(np.size(data1,1),dtype='bool')#just a vector with 1 element
    
   # print(decrease, increase)
    
    
    # print(data1)
    # print(np.shape(data1))
    
    alpha = 0.05
    for idx,a in enumerate(increase):#54
    
        #This loop compares proteins between two stages if they have be decreased or increased
        ##calculate p value for decrease between (G1,S) or (S,G2),checks p value < alpha or TRUE, then data1 is significantly greater than data 2
        
        #alternative Hypothesis 1, data 1 is significantly higher than data 2
        p_decrease=ttest_ind(data1[:,idx],data2[:,idx],equal_var=False,alternative='greater')[1]
        decrease[idx] =p_decrease<alpha
        
        #alternative Hypothesis 2
        ##calculate p value for increase between (G1,S) or (S,G2), then data1 is significantly less than data 2
        p_increase=ttest_ind(data1[:,idx],data2[:,idx],equal_var=False,alternative='less')[1]
        increase[idx] = p_increase<alpha
        
        
    # print(p_decrease)
    # print(p_increase)
    comparison_result = mean1<mean2 # it should consist the same result as increase
    
   # print(decrease, increase, comparison_result)
    return decrease,increase,comparison_result # decrease = a, increase =b, c= increase=b returns [ True] [False] [False]



# A plotting function that will plot average protein expression for a given protein for a 
# given dataset.  Provide the protein name as a string, and the dataset as a integer.
# I can provide dataset option context if needed; forgoing for now because datasets rely on
# my local path.
def plot_protein(protein_name,path):
    stage_strings = ['G1','S','G2']
    mean = np.zeros(len(stage_strings))
    std = np.zeros(len(stage_strings))

    time_strings = ['8','32','64','128','256']
    gray_shades = [
    '#F5F5F5',  # White Smoke
    '#DCDCDC',  # Gainsboro
    '#C0C0C0',  # Silver
    '#A9A9A9',  # Dark Gray
    '#808080',  # Gray
    '#696969',  # Dim Gray
    '#505050',  # Darker Gray
    '#303030']   # Very Dark Gray]
    #time_strings = ['2','8','16','32','64','256']
    mean = np.zeros((len(stage_strings),len(time_strings)+1)) #3 x 6
    std = np.zeros((len(stage_strings),len(time_strings)+1)) # 3 x 6
    data,name_list = load_fcs_as_matrix(path+'8min_UMAP_G1.fcs')
    idx = name_list.index(protein_name)
    for i in range(len(stage_strings)): #G1, S, G2
        data,name_list = load_fcs_as_matrix(path+'8min_UMAP_'+stage_strings[i]+'.fcs')
        
        # if stage_strings[i] == 'G2':  # If stage is G2, combine with M
        #     data_m, _ = load_fcs_as_matrix(path +'8min_UMAP_M.fcs')
        #     data = np.concatenate((data, data_m), axis=0)
        
        #we do not need this step if the data is already log-transformed
        #data[data<0] = 0
        
        mean[i,0] = np.mean(data[:,idx])
        std[i,0] = sem(data[:,idx])
        
        print(f'mean={mean[i,0]}')
        print(len(time_strings))
        for j in range(len(time_strings)): #time points
            data,names = load_fcs_as_matrix(path+time_strings[j]+'min_UMAP_'+stage_strings[i]+'.fcs')
            if stage_strings[i] == 'G2':  # If stage is G2, combine with M
                data_m, _ = load_fcs_as_matrix(path + time_strings[j] + 'min_UMAP_M.fcs')
                data = np.concatenate((data, data_m), axis=0)
            #we do not need this step if the data is already log-transformed
            #data[data<0] = 0
            
            mean[i,j] = np.mean(data[:,idx])
            std[i,j] = sem(data[:,idx])
            
            
    print(mean)
    #print(std)
    
    # Set global font settings
    plt.rcParams.update({
    'font.family': 'Arial',  # Use Arial font
    'font.size': 14,         # Set font size
    })
 
   
    plt.figure()
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['font.size'] = 18
    #plt.xlabel('Stage')
    #plt.ylabel('Average Protein Expression')
    # Access the Axes object
    ax = plt.gca()

    # Set the width of the x and y axes (spines)
    axes_width_size=2
    ax.spines['bottom'].set_linewidth(axes_width_size)
    ax.spines['left'].set_linewidth(axes_width_size)
    ax.spines['top'].set_linewidth(axes_width_size)  # Optional: set top spine width to 0 if not needed
    ax.spines['right'].set_linewidth(axes_width_size) 
    # Set the tick parameters to adjust tick width and font size
    plt.tick_params(axis='both', which='major', width=3, labelsize=20)
   

    plt.yscale('log') #Block it if doing in normal scale
    #plt.ylim(0,52)
    #plt.title(protein_name.decode('ascii'))
    ind = np.arange(3)
    print(ind)
    print(mean.shape)
    print(mean)
  
    if len(mean.shape) == 1:
        plt.bar(['G1','S','G2-M'],mean,yerr=std)
    else:
        width = np.min(np.diff(ind))/(len(mean[0,:])+1)
        
        for i in range(len(mean[0,:])-1): # 6 time points
        
            # Y_mean = mean[:, i]  # sept 4, 2024
            # Y_std = std[:, i]
            # plt.bar(ind+width*i, Y_mean, yerr=Y_std, width=width,
            #         color=gray_shades[i % len(gray_shades)])

            #geometric mean Sept 4th 2024, you can block the lower part if you want to show data only in log-transformed scale
            
            
            #print(f"y_mean={Y_mean}")
            plt.bar(ind+width*i,np.exp(mean[:,i]),width=width,color=gray_shades[i % len(gray_shades)])
            plt.xticks(ind + ((len(mean[0,:])-1)/2)*width , ('G1', 'S', 'G2-M'))
            
            print(f'{i},{np.exp(mean[:,i])}')
    # Adjust legend size and position

    #plt.ylim(0,100)#to print the legend only keep this line on
    #legend = plt.legend(('0 min','2min','4min','8min', '16min', '32min', '64min', '256min'), loc='upper right', fontsize='16', framealpha=1)
    #legend = plt.legend(('8min', '32min', '64min','128min', '256min'), loc='upper right', fontsize='16', framealpha=1)

    
    
# Which of the 9 possible conditions do our booleans generate??  Output used for figures.
# Remember, "a" is "has decrease", "b" is "has increase", and "c" is a raw comparison (no 
# statistical significance).  "c" is unused here.
def do_tests(a,b,c):
	#example: test 1 is that G1->S has a increase and S->G2 has a increase
    test1 = b[:,:,0] & b[:,:,1]
    # test 2 is that G1->S has a increase and S->G2 does not have a decrease or a increase
    test2 = b[:,:,0] & ~(a[:,:,1] | b[:,:,1])
    # etc
    test3 = b[:,:,0] & a[:,:,1] # increase -> decrease
    test4 = a[:,:,0] & b[:,:,1] # decrease -> increase
    test5 = a[:,:,0] & ~(a[:,:,1] | b[:,:,1])# dcrease -> flat
    test6 = a[:,:,0] & a[:,:,1] # decraese -> decrease
    test7 = ~(a[:,:,0] | b[:,:,0]) & b[:,:,1] #flat -> increase
    test8 = ~(a[:,:,0] | b[:,:,0]) & ~(a[:,:,1] | b[:,:,1]) #flat -> flat
    test9 = ~(a[:,:,0] | b[:,:,0]) & a[:,:,1] # #flat -> decrease
    return test1,test2,test3,test4,test5,test6,test7,test8,test9    # return a tuple


# This takes the output in the next function and converts it to a list of proteins in each
# category.  For example, if two decreases only happened for 4 proteins, it will return those
# protein names.  Not super useful anymore.
def give_names(test,names,protein_idx):
    name_str = []
    for name in names:
        name_str.append(name.decode('ascii'))
    pass_test = []
    for res in test:
        pass_test.append(np.asarray(name_str)[protein_idx][res[protein_idx].astype(bool)])
    return pass_test

# converts the results of all those booleans to a index for the figure result.  For instance,
# two succesive increases (test1) has index "0", an increase followed by "no change" (test2) 
# has index "1", etc.  "out" is the output, and it's given as an output argument.  The written
# file indicates which proteins belong to each category of trends.
def write_class_file(a,b,c,names):
    tests = do_tests(a,b,c)
    print(len(tests)) #9
    out = np.zeros((len(tests[0][:,0]),len(tests[0][0,:])))
    for idx,test in enumerate(tests):
        out[test] = idx
    return out

    

# Performs one-sided t-tests to determine if an increase occurred, and if a decrease occurred, 
# for two datasets.  Additionally, tests whether increases or decreases occurred, ignoring 
# statistical significance.  Returns a boolean vector indicating which data columns are significantly
# less in data2 than in data1, a boolean vector indicating same but greater in data2 than data1,
# and a boolean vector indicating whether the mean of data2 is greater (no statistical significance)
def specific_protein_behavior(protein_name,timepoint, cell_cycle_stage):   
    # Find the index of the protein in the name_list
    protein_index = name_list.index(protein_name)
    print(protein_index)
    # Extract the corresponding values from a, b, and c
    a_values = a[protein_index,timepoint,cell_cycle_stage]
    print(a.shape)
    b_values = b[protein_index,timepoint,cell_cycle_stage]
    c_values = c[protein_index,timepoint,cell_cycle_stage]

    # Print the index and corresponding values for verification
    print(a_values,b_values,c_values)
    
    if (cell_cycle_stage==0):
        print("G1 to S:")
    else:
        print("S to G2:")
        
   # Check for significant increase or decrease
    if c_values.any():  # Check if any value in c_values is True
        if b_values.any():  # Check if any value in b_values is True
            print(f"There is a protein increase, with p < 0.05 {b_values}")
        else:
            print("There is a protein increase with no significance that means no change of protein between stages ")
    else:
       if a_values.any():  # Check if any value in a_values is True
           print(f"There is a protein decrease, with p < 0.05 {a_values}")
       else:
           print("There is a protein decrease with no significance that means no change of protein between stages ")

if __name__=='__main__':
    
    
    #greg's gating on NKLs
    #path='./data/NKLs/NKG2D_umap/NKG2D_'
    #path='./data/without_cisplatin/NKL/NKG2D/try2/NKG2D_'
    
    
    path='./data/NKL_greg/NKG2D/NKG2D_'

    filename= path+'8min_UMAP_G1.fcs' #note: These are CD56 dim cells
    #filename= './data/primary/CD107_cell_cycle_gates/debarcode_IL2_NKG2D_2min_G1.fcs' #note: These are CD56 dim cells
    data,param_name_list=load_fcs_as_matrix(filename)
    

    # print(param_names)
    # print(len(param_names))
    # print(data[0:5][:])
    
    #print(data.shape)
    # print(type(data))


    # a is a two dimesional array for all proteins, 1st column represents if proteins decrease between G1 and S and 2nd column represents the protein increase between S and G2 stages
    a,b,c,name_list=analyze_protein_variation(path)

    #print(a)
    
    #writing the file for colormap
    test_idx = write_class_file(a,b,c,name_list)
    
    print(name_list)
    # # Save to CSV file
    # np.savetxt('CD107_high.csv', test_idx, delimiter=',', fmt='%d')

    print(test_idx)
    print(test_idx.shape)

    
    # Ensure name_list is a list of strings
    name_list_str = [name.decode('ascii') if isinstance(name, bytes) else name for name in name_list]

    # Convert test_idx to a pandas DataFrame
    df_test_idx = pd.DataFrame(test_idx)

    # Add the name_list as the first column
    df_test_idx.insert(0, 'Protein Name', name_list_str)

    # Save the DataFrame to a CSV file
    #df_test_idx.to_csv('greg_umap_colormap_logtransformed_aug26.csv', index=False) #skip the previous file and look at the file now
    df_test_idx.to_csv('NKG2D_colormap_logtransformed_feb8_25_G2_M_combined_yes.csv', index=False) #skip the previous file and look at the file now

 
    # Specify the protein name you are interested in
    #protein_name = b'CD16'
    protein_name = b'pS6'
    timepoint=4 #(8 ,32, 64,128, 256)
    cell_cycle_stage=0 #{g1-S, S-G2]}
    
    plot_protein(protein_name,path)
    
    specific_protein_behavior(protein_name,timepoint, cell_cycle_stage)
    

