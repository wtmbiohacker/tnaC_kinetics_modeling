# run simulation of tnaC dynamics for given times

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np 
from numpy import random
import math
from random import shuffle
import ConfigParser
#plt.rcParams["font.family"] = "Times New Roman"

# The function processing the configure or RNA file
def configprocess(text):
    config=ConfigParser.ConfigParser()
    config.read(text)
    if len(config.sections())>1:
        print('Make a mistake in the configure file') and os._exit(0) 
    section=''.join(config.sections())
    varlist={option:config.get(section,option) for option in config.options(section)}
    return varlist

# judge the current state of RNAP, determine the available options in the next step of simulation
def judge_RNAP(RNAP_position, Ribosome_position, Rho_position, RNAlist):
    '''
    Parameters
    ----------
    RNAP_position: the current position of RNAP within tnaC transcript
    Ribosome_position: the current position of Ribosome within tnaC transcript, -1 before translation initiation; 0, ready for translation, (0,75), normal translation elongation; >=75 translation finalization
    Rho_position: the current position of Rho factor
    RNAlist: the structure file of tnaC transcript
    Return value
    ----------
    'attenuation', 'stemploop', 'terminator' or 'elongation', refers to premature transcription termination, stalling at stemloop, terminator or can conduct normal elongation, respectively
    '''
    # get the values for RNA structure
    stemloop = int(RNAlist['stemloop'])
    terminator = int(RNAlist['terminator'])
    # determine state
    if RNAP_position == Rho_position:
        return 'attenuation'
    elif RNAP_position != stemloop and RNAP_position != terminator:
        return 'elongation'
    elif RNAP_position == stemloop:
        if RNAP_position == Ribosome_position:  # the stemloop can be overcome given ribosome catches up with RNAP
            return 'elongation'
        else:
            return 'stemloop'
    elif RNAP_position == terminator:
        return 'terminator'

# determien the action of the next step for RNAP based on the current RNAP state
def act_RNAP(RNAP_position, varlist, RNAlist, RNAP_state):
    '''
    Parameters
    ----------
    RNAP_position: the current position of RNAP within tnaC transcript
    varlist: the file for inherent parameters of macromolecules
    RNAlist: the structure file of tnaC transcript
    RNAP_state: the current state of RNAP
    Return value
    ----------
    number of nucleotides transcribed by RNAP in this step
    '''
    transcribed_nt = 0
    # get the values for all variables
    F_beyond_stemloop =  float(varlist['f_beyond_stemloop'])
    F_terminator = float(varlist['f_terminator'])
    RNAP_elongation = float(varlist['rnap_elongation'])
    # get the values for RNA structure
    stemloop = int(RNAlist['stemloop'])
    terminator = int(RNAlist['terminator'])
    # determine action
    if RNAP_state == 'stemloop':
        transcribed_nt = RNAP_elongation if random.poisson(F_beyond_stemloop)>=1 else 0
    elif RNAP_state == 'terminator':
        transcribed_nt = RNAP_elongation if random.poisson(F_terminator)>=1 else 0
    elif RNAP_position < stemloop and (RNAP_position + RNAP_elongation) >= stemloop:
        transcribed_nt = stemloop - RNAP_position
    elif RNAP_position < terminator and (RNAP_position + RNAP_elongation) >= terminator:
        transcribed_nt = terminator - RNAP_position
    elif RNAP_state == 'elongation': # normal elongation
        transcribed_nt = RNAP_elongation
    return transcribed_nt

# judge the current state of ribosome, determine the available options in the next step of simulation
def judge_Ribosome(Ribosome_position, trp_ribosome, RNAlist):
    '''
    Parameters
    ----------
    Ribosome_position: the current position of Ribosome within tnaC transcript, -1 before translation initiation; 0, ready for translation, (0,75), normal translation elongation; >=75 translation finalization
    trp_ribosome: binary (0 or 1), whether the ribosome has two Trp molecules within the tunnel (1) or not (0)
    RNAlist: the structure file of tnaC transcript
    Return value
    ----------
    'before_initiation', 'trp_absence_stalling', 'Trp_presence_stalling', 'elongation' or 'finalization', refers to prior to translation initiation, stalling in the absence or presence of Trp, normal elongation or translation finalization, respectively
    '''
    # get the values for RNA structure
    Trp_absence_stalling = int(RNAlist['trp_absence_stalling'])
    Trp_presence_stalling = int(RNAlist['trp_presence_stalling'])
    stop_codon = int(RNAlist['stop_codon'])
    # determine the state
    if Ribosome_position == -1:
        return 'before_initiation'
    elif Ribosome_position >= stop_codon:
        # translation is finalized
        return 'finalizaion'
    elif Ribosome_position >= 0:
        if Ribosome_position not in [Trp_absence_stalling, Trp_presence_stalling]:
            return 'elongation'
        elif Ribosome_position == Trp_absence_stalling:
            return 'trp_absence_stalling'
        elif Ribosome_position == Trp_presence_stalling:
            if trp_ribosome == 0:
                return 'elongation'
            elif trp_ribosome == 1:
                return 'trp_presence_stalling'

# determien the action of the next step for Ribosome based on the current Ribosome state and RNAP position
def act_Ribosome(Ribosome_position, RNAP_position, varlist, RNAlist, Ribosome_state, trp_ribosome):
    '''
    Parameters
    ----------
    Ribosome_position: the current position of Ribosome within tnaC transcript
    RNAP_position: the current position of RNAP within tnaC transcript
    varlist: the file for inherent parameters of macromolecules
    RNAlist: the structure file of tnaC transcript
    Ribosome_state: the current state of Ribosome
    trp_ribosome: determine whether the ribosome used in this simulation run has trp within exit tunnel or not
    Return value
    ----------
    number of nucleotides translated by Ribosome in this step
    '''
    translated_nt = 0
    # get the values for all variables
    F_ribo_attack = float(varlist['f_ribo_attack'])
    F_stalling_Trp_no = float(varlist['f_stalling_trp_no'])
    F_stalling_Trp_yes = float(varlist['f_stalling_trp_yes'])
    Ribosome_elongation = float(varlist['ribosome_elongation'])
    # get the values for RNA structure
    Trp_absence_stalling = int(RNAlist['trp_absence_stalling'])
    Trp_presence_stalling = int(RNAlist['trp_presence_stalling'])
    stop_codon = int(RNAlist['stop_codon'])
    # determine the action
    if Ribosome_state == 'before_initiation':
        translated_nt = 1 if random.poisson(F_ribo_attack)>=1 else 0
    elif Ribosome_state == 'trp_absence_stalling':
        translated_nt = Ribosome_elongation if random.poisson(F_stalling_Trp_no)>=1 else 0
    elif Ribosome_state == 'trp_presence_stalling':
        translated_nt = Ribosome_elongation if random.poisson(F_stalling_Trp_yes)>=1 else 0
    elif Ribosome_state == 'elongation':
        if (Ribosome_position + Ribosome_elongation) > RNAP_position:
            translated_nt = (RNAP_position/3)*3 - Ribosome_position  # keep the ribosome at the first nucleotide of the codon that is the most proximal to RNAP
        else:
            translated_nt = Ribosome_elongation
    elif Ribosome_state == 'finalization':
        translated_nt = 0
    return translated_nt

# judge the current state of Rho factor, determine the available options in the next step of simulation
# Mechanically, we can regard the action of Ribosome and Rho factor as simutaneously
def judge_Rho(Rho_position, Ribosome_state, RNAP_position, RNAlist):
    '''
    Parameters
    ----------
    Rho_position: the current position of Rho factor, -1 before binding; rutA (defined in RNAlist, e.g. 110), ready for translation; (>rutA,75), normal translocation; >=75 translation finalization
    Ribosome_state: 'before_initiation', 'trp_absence_stalling', 'Trp_presence_stalling', 'elongation' or 'finalization', refers to prior to translation initiation, stalling in the absence or presence of Trp, normal elongation or translation finalization, respectively
    RNAP_position: the current position of RNAP within tnaC transcript
    RNAlist: the structure file of tnaC transcript
    Return value
    ----------
    Rho_state: 'cannot_initiation', 'before_initiation' or 'translocation', refers to cannot initiation due to unavailability of rutA (not transcribed yet or blocked by Ribosome), attacking rutA but not succeed yet and normal elongation, respectively
    '''
    # get the values for RNA structure
    rutA = int(RNAlist['ruta'])
    # determine the state
    if Rho_position == -1:
        if RNAP_position < rutA:
            return 'cannot_initiation'
        elif Ribosome_state in ['trp_absence_stalling', 'trp_presence_stalling']:
            return 'cannot_initiation'
        else:
            return 'before_initiation'
    elif Rho_position >= rutA:
        return 'elongation'

# determien the action of the next step for Ribosome based on the current Ribosome state and RNAP position
def act_Rho(Rho_position, RNAP_position, varlist, RNAlist, Rho_state):
    '''
    Parameters
    ----------
    Rho_position: the current position of Rho factor, -1 before binding; rutA (defined in RNAlist, e.g. 110), ready for translation; (>rutA,75), normal translocation; >=75 translation finalization
    RNAP_position: the current position of RNAP within tnaC transcript
    varlist: the file for inherent parameters of macromolecules
    RNAlist: the structure file of tnaC transcript
    Rho_state: 'cannot_initiation', 'before_initiation' or 'translocation', refers to cannot initiation due to unavailability of rutA (not transcribed yet or blocked by Ribosome), attacking rutA but not succeed yet and normal elongation, respectively
    Return value
    ----------
    number of nucleotides translocated by Rho in this step
    '''
    translocated_nt = 0
    # get the values for all variables
    F_rho_attack = float(varlist['f_rho_attack'])
    Rho_translocation = float(varlist['rho_translocation'])
    # get the values for RNA structure
    rutA = int(RNAlist['ruta'])
    # determine the action
    if Rho_state == 'cannot_initiation':
        translocated_nt = 0
    elif Rho_state == 'before_initiation':
        translocated_nt = rutA + 1 if random.poisson(F_rho_attack)>=1 else 0 # make the Rho binds rutA site if successful attack
    elif Rho_state == 'elongation':
        if (Rho_position + Rho_translocation) > RNAP_position:
            translocated_nt = RNAP_position - Rho_position  # keep the ribosome at the first nucleotide of the codon that is the most proximal to RNAP
        else:
            translocated_nt = Rho_translocation
    return translocated_nt

# one simulation given the values of variables
# algorithm:
# 1. extract all values for macromolecule kinetics and RNA milestone structure positions
# 2. for each step, use three judge functions to determine the action of the three macromolecules in this step
# 3.1 If RNAP goes beyond the terminator and is not caught up with by Rho, it gives rise to a successful tnaA transcription;
# 3.2 If RNAP is upstream of the terminator and is not caught up with by Rho, go on;
# 3.3 If RNAP is upstream of the terminator and caught up with by Rho, it gives rise to an unsuccessful tnaA transcription;
def onerun(varlist, RNAlist):
    '''
    Parameters
    ----------
    varlist: the file for inherent parameters of macromolecules
    RNAlist: the structure file of tnaC transcript
    Return value
    ----------
    tnaA_flag: 'Yes' or 'No' binary, successful transcription or not
    trp_ribosome: 1 or 0 binary, determine whether the ribosome used in this simulation run has trp within exit tunnel or not
    kinetics dataframe: column: simulation steps (step), positions for RNAP (RNAP), Ribosome (Ribosome) and Rho (Rho), respectively. Updated by each round of simulation 
    '''
    kinetics_Dic={}
    kinetics_Dic['step']=[]
    kinetics_Dic['RNAP']=[]
    kinetics_Dic['Ribosome']=[]
    kinetics_Dic['Rho']=[]
    RNAP_position = 0
    Ribosome_position = -1
    Rho_position = -1
    simulation_step = 0
    # get the values for RNA structure
    terminator = int(RNAlist['terminator'])
    # determine whether the ribosome used in this simulation run has trp within exit tunnel or not
    trp_ribosome = 1 if random.poisson(float(varlist['trp_ribosome']))>=1 else 0
    tnaA_flag = 'No' # whether tnaA is successfully transcribed or not
    kinetics_Dic={} # kinetics storage file
    kinetics_Dic['step']=[]
    kinetics_Dic['RNAP']=[]
    kinetics_Dic['Ribosome']=[]
    kinetics_Dic['Rho']=[]
    while RNAP_position <= terminator:
        simulation_step = simulation_step + 1
        # //////////////////////////////////////////////////////////////
        # test RNAP action in this step
        RNAP_state = judge_RNAP(RNAP_position, Ribosome_position, Rho_position, RNAlist)
        transcribed_nt = act_RNAP(RNAP_position, varlist, RNAlist, RNAP_state)
        RNAP_position = RNAP_position + transcribed_nt
        if RNAP_position > terminator:
            tnaA_flag = 'Yes'
            break
        # //////////////////////////////////////////////////////////////
        # test Ribosome action in this step
        Ribosome_state = judge_Ribosome(Ribosome_position, trp_ribosome, RNAlist)
        translated_nt = act_Ribosome(Ribosome_position, RNAP_position, varlist, RNAlist, Ribosome_state, trp_ribosome)
        Ribosome_position = Ribosome_position + translated_nt
        # //////////////////////////////////////////////////////////////
        # test Rho action in this step
        Rho_state = judge_Rho(Rho_position, Ribosome_state, RNAP_position, RNAlist)
        translocated_nt = act_Rho(Rho_position, RNAP_position, varlist, RNAlist, Rho_state)
        Rho_position = Rho_position + translocated_nt
        # //////////////////////////////////////////////////////////////
        # recheck whether RNAP is caught up with by Rho after Rho translocation in this step
        RNAP_state = judge_RNAP(RNAP_position, Ribosome_position, Rho_position, RNAlist)
        if RNAP_state == 'attenuation':
            break
        # //////////////////////////////////////////////////////////////
        # update the kinetics record after one round of simulation
        kinetics_Dic['step'].append(simulation_step)
        kinetics_Dic['RNAP'].append(RNAP_position)
        kinetics_Dic['Ribosome'].append(Ribosome_position)
        kinetics_Dic['Rho'].append(Rho_position)
    return tnaA_flag, trp_ribosome, pd.DataFrame(kinetics_Dic)

# plot the kinetics of the three macromolecules for intuition
def plot_kinetics(kinetics_DF, varlist, RNAlist, trp_ribosome, name):
    '''
    Parameters
    ----------
    kinetics_DF: column: simulation steps (step), positions for RNAP (RNAP), Ribosome (Ribosome) and Rho (Rho), respectively. Updated by each round of simulation 
    varlist: the file for inherent parameters of macromolecules
    RNAlist: the structure file of tnaC transcript
    trp_ribosome: determine whether the ribosome used in this simulation run has trp within exit tunnel or not
    name: name of the plot produced here
    Return value
    ----------
    1 if successful
    '''
    # get the values for all variables
    F_beyond_stemloop =  float(varlist['f_beyond_stemloop'])
    RNAP_elongation = float(varlist['rnap_elongation'])
    # get the values for RNA structure
    stemloop = int(RNAlist['stemloop'])
    Trp_absence_stalling = int(RNAlist['trp_absence_stalling'])
    Trp_presence_stalling = int(RNAlist['trp_presence_stalling'])
    stop_codon = int(RNAlist['stop_codon'])
    rutA = int(RNAlist['ruta'])
    terminator = int(RNAlist['terminator'])
    # prepare the data
    time_record = np.array(kinetics_DF['step'])
    RNAP_kinetics = np.array(kinetics_DF['RNAP'])
    Ribosome_kinetics = np.array(kinetics_DF['Ribosome'])
    Rho_kinetics = np.array(kinetics_DF['Rho'])
    # plot the milestone structures in RNA
    end_time = np.max(time_record)
    plt.figure(figsize=(5,7))
    plt.plot([0, end_time], [stemloop, stemloop], ls = '--', color = '#626567', alpha = 0.5, label = 'stemloop')
    plt.plot([0, end_time], [stop_codon, stop_codon], ls = ':', color = '#626567', alpha = 0.5, label = 'stop codon')
    plt.plot([0, end_time], [Trp_absence_stalling, Trp_absence_stalling], ls = '-.', color = '#48C9B0', alpha = 0.5, label = 'Trp($-$) stalling')
    plt.plot([0, end_time], [Trp_presence_stalling, Trp_presence_stalling], ls = '-.', color = '#0E6251', alpha = 0.5, label = 'Trp($+$) stalling')
    plt.plot([0, end_time], [rutA, rutA], ls = '--', color = '#873600', alpha = 0.5, label = '$rut$')
    plt.plot([0, end_time], [terminator, terminator], ls = ':', color = '#873600', alpha = 0.5, label = 'terminator')
    # plot the kinetics
    plt.plot(time_record, RNAP_kinetics, color = '#B03A2E', label = 'RNAP')
    plt.plot(time_record, Ribosome_kinetics, color = '#2874A6', label = 'Ribosome')
    plt.plot(time_record, Rho_kinetics, color = '#76448A', label = 'Rho')
    expected_end_time = 1/F_beyond_stemloop + terminator/RNAP_elongation
    plt.xlim([1, max(expected_end_time, end_time)])
    plt.ylim([0, terminator+5])
    plt.xlabel('Time (step)',fontsize=20)
    plt.ylabel('Nucleotide in transcript',fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper left', fontsize=14)
    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.95)
    plt.savefig('%s.png'%(name), dpi = 600)
    plt.clf()
    return 1

# the main function
# run simulations for determined times, collect the result, store the representative kinetics profile
def main(configureFile, RNAFile, prefix):
    varlist = configprocess(configureFile)
    RNAlist = configprocess(RNAFile)
    Fre_storage=200 # frequency to store kinetics
    successful_tnaA = 0
    # run the simulation by number of times according defined in configureFile
    for i in range(int(float(varlist['simulation']))):
        tnaA_flag, trp_ribosome, kinetics_DF = onerun(varlist, RNAlist)  # each simulation
        onerun_result = 1 if tnaA_flag == 'Yes' else 0
        successful_tnaA = successful_tnaA + onerun_result
        if i%Fre_storage == 0:  # store limited number of representative kinetics profiles
            plot_kinetics(kinetics_DF, varlist, RNAlist, trp_ribosome, 'kinetics%d'%(i/Fre_storage))
    successful_ratio = float(successful_tnaA) / float(varlist['simulation'])
    print ('tnaA expression: %.4f'%(successful_ratio))
    os.system('cat /dev/null > temp_result.txt')
    g=open('temp_result.txt', 'r+')
    g.write('%s\n'%(str(successful_ratio)))
    g.close()

if __name__ == '__main__':
    try:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)
