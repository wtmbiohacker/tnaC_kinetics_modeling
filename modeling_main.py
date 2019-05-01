# This script is the main script for tnaC kinetics simulation
# The following modules
# 1. Process the configure file, get the value for all constants, the variable to be studied, etc;
# 2. Given the constants and variables (tens of values to be tested), run simulation, get the result
# 3. Store and visualization of the result

# one input files:

# configure

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import math
import ConfigParser
plt.rcParams["font.family"] = "Times New Roman"

# The function processing the configure file
def configprocess(text):
    config=ConfigParser.ConfigParser()
    config.read(text)
    if len(config.sections())>1:
        print('Make a mistake in the configure file') and os._exit(0) 
    section=''.join(config.sections())
    varlist={option:config.get(section,option) for option in config.options(section)}
    return varlist

# input inspection
def metricJudge(varlist):
    f_ribo_attack    = True if float(varlist['f_ribo_attack'])>0 else False
    f_beyond_stemloop = True if float(varlist['f_beyond_stemloop'])>0 else False
    f_stalling_trp_no = True if float(varlist['f_stalling_trp_no'])>0 else False
    f_stalling_trp_yes = True if float(varlist['f_stalling_trp_yes'])>0 else False
    trp_ribosome = True if float(varlist['trp_ribosome'])>=0 else False
    f_rho_attack = True if float(varlist['f_rho_attack'])>0 else False
    f_terminator  = True if float(varlist['f_terminator'])>0 else False
    rnap_elongation  = True if float(varlist['rnap_elongation'])>0 else False
    ribosome_elongation  = True if float(varlist['ribosome_elongation'])>0 else False
    rho_translocation  = True if float(varlist['rho_translocation'])>0 else False
    variable  = True if varlist['variable'] in ['f_ribo_attack', 'f_beyond_stemloop', 'f_stalling_trp_no', 'f_stalling_trp_yes', 'trp_ribosome', 'f_rho_attack', 'f_terminator'] else False
    variable_range = True if len(varlist['variable_range'].split(','))==2 and float(varlist['variable_range'].split(',')[0])>0 else False
    increment_log = True if varlist['increment_log'] in ['Yes', 'No'] else False
    variable_number  = True if varlist['variable_number'].isdigit() else False
    simulation  = True if varlist['simulation'].isdigit() else False
    step  = True if float(varlist['step'])>0 and float(varlist['step'])<1 else False
    for item in varlist:
        if not eval(item) :
            print('The varaint %s is wrong, please input the correct value!!!!'%(item))
            print(varlist[item])
            os._exit(1)
    # ribosome must translate integer number of codons per step
    if float(varlist['ribosome_elongation'])*float(varlist['step'])%3 != 0:
        print ('ribosome must translate integer number of codons per step')
        os._exit(1)

# return the value list of variable to be studied uniformly distributed in either linear or log10 space
def get_variable_value(variable_range, increment_log, variable_number):
    valueLst = []
    upper_bound = float(variable_range.split(',')[1])
    lower_bound = float(variable_range.split(',')[0])
    if increment_log == 'Yes':
        upper_bound = math.log(upper_bound, 10)
        lower_bound = math.log(lower_bound, 10)
        increment = (upper_bound - lower_bound) / (variable_number - 1)
        for number in range(variable_number):
            valueLst.append(10 ** (lower_bound + increment * number))
    elif increment_log == 'No':
        increment = (upper_bound - lower_bound) / (variable_number - 1)
        for number in range(variable_number):
            valueLst.append(lower_bound + increment * number)
    return valueLst

# construct temp configure file for simulation usage, the variable to be studied is given a particular value
# For variables related to time scale, normalization by step; otherwise (trp_ribosome and simulation), do not perform this.
def construct_temp_configure(varlist, variable_value):
    os.system('cat /dev/null > temp_configure.txt')
    g = open ('temp_configure.txt', 'r+')
    g.write('[configdesign]\n')
    for item in varlist:
        if item == varlist['variable']:  # if it is the variables to be investigated
            if item == 'trp_ribosome':
                g.write('%s:%s\n'%(item, str(variable_value)))
            elif item in ['f_ribo_attack', 'f_beyond_stemloop', 'f_stalling_trp_no', 'f_stalling_trp_yes', 'f_rho_attack', 'f_terminator']:
                g.write('%s:%s\n'%(item, str( variable_value * float(varlist['step']) )))
        elif item in ['f_ribo_attack', 'f_beyond_stemloop', 'f_stalling_trp_no', 'f_stalling_trp_yes', 'f_rho_attack', 'f_terminator', 'rnap_elongation', 'ribosome_elongation', 'rho_translocation']: 
            value_stepwise = float(varlist[item])*float(varlist['step'])
            g.write('%s:%s\n'%(item, str(value_stepwise)))
        elif item in ['simulation', 'trp_ribosome']:
            g.write('%s:%s\n'%(item, str(varlist[item])))
    g.close()
    return 'temp_configure.txt'

# extract the value for one run of simulation
# the result file has one line, a float number, ratio of successful tnaA expression in all transcription
def get_result(result_file, simulation):
    f = open(result_file, 'r')
    result = 0.0
    for i, line in enumerate(f):
        if i == 0:
            if float(line.rstrip())>0:
                result = float(line.rstrip())
            else: # == 0, no expression
                result = 1/float(simulation)/2 # < 1 expression in N simulations
    f.close()
    os.system('rm %s'%(result_file))
    return result

# visualize the trend of expression vs. variable from simulation
def scatter_plot(increment_log, Xaxis, Yaxis, xlabel, ylabel, name, simulation):
    """
    funciton to draw scatter plot
    Parameters
    ________________
    numpy array: Xaixs, Yaxis
    output_dir, name
    """
    plt.scatter(Xaxis,Yaxis,s=30,color='#5DADE2')
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.xlim([min(Xaxis), max(Xaxis)])
    if increment_log:
        plt.xscale('log')
    #plt.ylim([min(Yaxis), max(Yaxis)])
    plt.ylim([1/float(simulation),1])
    plt.yscale('log')
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    plt.savefig('%s.png'%(name),dpi=200)
    plt.clf()
    return 0

# the main function
def main(configureFile, RNAfile):
    varlist = configprocess(configureFile)
    metricJudge(varlist)
    for item in varlist:
        print ('%s: %s'%(item, varlist[item]))
    valueLst = get_variable_value(varlist['variable_range'], varlist['increment_log'], int(varlist['variable_number']))
    resultLst = []
    output_dir = varlist['variable'] + '_' + varlist['variable_range']
    name = varlist['variable'] + '_' + varlist['variable_range']
    os.system('mkdir %s'%(output_dir))
    # used to store the final data
    data_Dic={}
    data_Dic[varlist['variable']]=valueLst
    data_Dic['expression']=[]
    for variable_value in valueLst:
        output_subdir = output_dir + '/' + str(variable_value)+'/'
        os.system('mkdir %s'%(output_subdir))
        temp_configure_file = construct_temp_configure(varlist, variable_value)
        prefix = varlist['variable'] + '_' + str(variable_value)
        os.system('python simulation_onerun.py %s %s %s'%(temp_configure_file, RNAfile, prefix))
        this_result = get_result('temp_result.txt', int(varlist['simulation']))
        resultLst.append(this_result)
        data_Dic['expression'].append(this_result)
        os.system('mv %s %s'%(temp_configure_file, output_subdir))
        os.system('mv *.png %s'%(output_subdir))
        print 'one simulation finalized successfully: %s %.5f'%(varlist['variable'], variable_value)
        print ('')
    os.system('cp %s %s'%(configureFile, output_dir))
    Xlabel = varlist['variable']
    scatter_plot(varlist['increment_log'], np.array(valueLst), np.array(resultLst), Xlabel, 'Expression (A.U.)', name, int(varlist['simulation']))
    os.system('mv %s.png %s'%(name, output_dir))
    data_DF = pd.DataFrame(data_Dic)
    data_DF.to_csv('%s/%s.csv'%(output_dir, name), sep='\t', index=None)
    os.system('rm temp* *.png')

if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2])
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)
