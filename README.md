# tnaC_kinetics_modeling

## What is this?
This python script collection is a stochatic model to simulate the dynamics of transcription and translation of bacterial tryptophan sensor (tnaC), which control the downstream gene expression in a tryptophan-dependent manner. The paper describing this program is currently under review. Please cite this paper if this program is useful to your work.

## General description of the algorithm
A comprehensive list of essential intermediate states of transcription-translation dynamics enables tnaC to regulate downstream gene expression in a Trp-dependent manner. In particular, the movements of three macromolecules, namely RNAP, the ribosome and Rho factor, are precisely spatiotemporally coordinated during the transcription of tnaC and the subsequent translation. To mimic this process, we developed this stochatis model [algorithm](./image/algorithm.png) to simulate the movement of RNAP, ribosome and Rho, during the transcription and translation of tnaC. Here, one round of modeling can be regarded as one tnaC transcription event; and its output is binary, depending on whether RNAP is able to transcribe through the entire intergenic region between tnaC and tnaA (success), or is halted by Rho factor before reaching tnaA (failure). Hence, the results of multiple rounds of dynamic modeling allowed us to calculate the ratio of successful tnaA expression compared with failed expression, which could be used to quantify the tnaC sensor response. The timescales of critical slow steps that are essential to the sensor response are regarded as the core variables of the model. For a specific bottleneck step, its timescale is regarded as the expectation of a Poisson distribution; and the time needed for the relevant macromolecule to go beyond this bottleneck step follows this Poisson distribution. Given these timescales as core parameters of the model, we can evolve the dynamics of RNAP, the ribosome and Rho factor (one tnaC transcription event), which determines whether tnaA is successfully expressed or not. Here we present two simulated macromolecular dynamics during tnaC transcription, corresponding to [successful expression](./image/Success.png) and [failed expression](./image/Failure.png).

Briefly, we defined a configure file that specifying the timescale for each bottleneck step of the macromolecule movement. Then, one of them is regarded as a variable to be tested. N rounds of simulations are executed to calculate the downstream gene expression, given different values of this variable, while all other bottleneck step timescales held constant. The response (downstream gene expression) vs. variable value is recorded and plotted.

## How to use it?
MacOS or Linux is supported environment to run this tool; for Windows users, some system shell commands in the scripts ('os.system(...)') need to be revised.

### Step 1：Installation and dependency
1. Install Python version 2.7
2. Install Scipy version 0.19.1 or above
3. Install Matplotlib version 2.0.2 or above
4. Install Numpy version 1.13.1 or above
5. Install Pandas version 0.18.1 or above

### Step 2：Prepare the configure file.
One example is shown below. In this case, "f_ribo_attack" is regarded as a variable to test its impact on sensor response. "variable_range", "increment_log" and "variable_number" together define the values set to this variable during simulation. "variable_range" specifies the minimum and maximum values; "variable_number" specifies the number of values set to the variable within this range; "increment_log" specifies to set the values in log10 ("Yes" argument) or linear ("No" argument) space. For the detailed description of these parameters and their corresponding units, see Table S2 of our paper. In this flat configure file, the name and value of the parameters are seperated by ":".

|parameter|value|
|---------|-----|
f_ribo_attack|0.25
f_beyond_stemloop|0.8
f_stalling_trp_no|0.3
f_stalling_trp_yes|0.06
trp_ribosome|100
f_rho_attack|5
f_terminator|1
rnap_elongation|40
ribosome_elongation|30
rho_translocation|100
variable|f_ribo_attack
variable_range|0.99,1
increment_log|Yes
variable_number|20
simulation|2000
step|0.1

### Step 3：Prepare the RNA structure file.
This file specifies the positions for important structures in tnaC transcript, which corresponds to the timescale values defined in the configure file. For the detailed description of these parameters and their corresponding units, see Table S1 of our paper. One example is shown below. In this flat configure file, the name and value of the parameters are seperated by ":".

|parameter|value|
|---------|-----|
stemloop|63
trp_absence_stalling|69
trp_presence_stalling|72
stop_codon|75
ruta|90
terminator|210

### Step 4：Run the pipeline
Put all necessary files mentioned above, as well as the data directory under the working directory.
Open the command line window, cd to the working directory and run the analysis pipeline.
cd path_to_your_working_directory
python HTsensor_main.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 
python HTsensor_main.py example_configure.txt

## Output files
The output files will be organized in the subdirectory whose name is specified by the variable name and its range. For instance, the example post here results in a . We term this subdirectory 'result directory' thereafter.

Under this result directory, N folders corresponding to N values of the variable, stores the representative simulated macromolecular dynamics when using these variable values.

In addition, one .csv file stores the sensor response for each simulated variable value, the main result of modeling. Here is an example.

Moreover, an expression vs. variable [plot](./image/result.png) is also presented.

