# EROS HUMAN
For animal version please swith to **"master"** branch.

##Experimental code for Event Realated Oscilation detection.

Developed software calculates advanced time-frequency quantities addressing evoked and induced EEG events. The software is designed to analyze both animal and human EEG data in order to provide translation between the observed phenomena. The analyzer is based on MATLAB platform and accepts LabChart and EDF data formats. More specifically, Event Related Oscillations are calculated to address induced oscillations. Phase Locking Index evaluates evoked oscillations over the trials and Phase Difference Locking Index measures functional connections between selected electrodes. The main outputs of the analyzer are the time-frequency characteristics of the quantities mentioned above.

##How to use EROS HUMAN
There are one executable scripts (in MATLAB) in repository:

```MATLAB
EROS_HUMAN.m
```

The script accepts EDF+ data format and reads EDF+ files with triggers (anotation).

####EROS_HUMAN.m

Calculates, visualizes, and stores Event Related Ocilations (ERO), Phase Locking Index (PLI), and Event Related Potencial (ERP). Additional output is Averaged Event Related Ocilations (AVG_ERO), which is the ERO obtained from ERP.

###Workflow

####Configure parallel computation:
For optimal performance set a number of parallel workers: Prallel->Manage Cluster Profiles->Cluster Profile->Edit->NumWorkers

For older MATLAB versions execute "matlabpool open" before calculation.

####Open your dataset:

####select multiple electrodes (use control+click or multiple selectin) the script evaluates all


* select two triggers (target/non-target) - the script operates on the fixed protocol
* select files to visualize - the script gives names to figures and subplots automatically
* MEAN SUBJECT structure contains averages values of ERP,PLI,ERP over all files




