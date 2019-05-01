# tnaC_kinetics_modeling

## What is this?
This python script collection is a stochatic model to simulate the dynamics of transcription and translation of bacterial tryptophan sensor (tnaC), which control the downstream gene expression in a tryptophan-dependent manner. The paper describing this program is currently under review. Please cite this paper if this program is useful to your work.

## General description of the algorithm
A comprehensive list of essential intermediate states of transcription-translation dynamics enables tnaC to regulate downstream gene expression in a Trp-dependent manner. In particular, the movements of three macromolecules, namely RNAP, the ribosome and Rho factor, are precisely spatiotemporally coordinated during the transcription of tnaC and the subsequent translation. To mimic this process, we developed [this stochatis model](./image/FACS-seq_schematics.png) to simulate the movement of RNAP, ribosome and Rho, during the transcription and translation of tnaC. Here, one round of modeling can be regarded as one tnaC transcription event; and its output is binary, depending on whether RNAP is able to transcribe through the entire intergenic region between tnaC and tnaA (success), or is halted by Rho factor before reaching tnaA (failure). Hence, the results of multiple rounds of dynamic modeling allowed us to calculate the ratio of successful tnaA expression compared with failed expression, which could be used to quantify the tnaC sensor response. The timescales of critical slow steps that are essential to the sensor response are regarded as the core variables of the model. For a specific bottleneck step, its timescale is regarded as the expectation of a Poisson distribution; and the time needed for the relevant macromolecule to go beyond this bottleneck step follows this Poisson distribution. Given these timescales as core parameters of the model, we can evolve the dynamics of RNAP, the ribosome and Rho factor (one tnaC transcription event), which determines whether tnaA is successfully expressed or not. Two simulated macromolecular dynamics during tnaC transcription are [Here](./image/FACS-seq_schematics.png).

## How to use it?
MacOS or Linux is supported environment to run this tool; for Windows users, some system shell commands in the scripts ('os.system(...)') need to be revised.

### Step 1：Installation and dependency
1. Install Python version 2.7
2. Install Scipy version 0.19.1 or above
3. Install Matplotlib version 2.0.2 or above
4. Install Numpy version 1.13.1 or above
5. Install Pandas version 0.18.1 or above
