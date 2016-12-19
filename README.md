#EROs, PLI, & DPLI Analyzer
##Authors:
* Vlastimil Koudelka
* Grygoriy Tsenov

This is for rodents, for human version switch to **"human"** branch, please. 

##Experimental code for Event Related Oscillation detection.
Developed software calculates advanced time-frequency quantities addressing evoked and induced EEG events. The software is designed to analyze both animal and human EEG data in order to provide translation between the observed phenomena. The analyzer is based on MATLAB platform and accepts LabChart and EDF data formats. More specifically, Event Related Oscillations are calculated to address induced oscillations. Phase Locking Index evaluates evoked oscillations over the trials and Phase Difference Locking Index measures functional connections between selected electrodes. The main outputs of the analyzer are the time-frequency characteristics of the quantities mentioned above.

##How to use EROS RODENTS

There are two executable scripts (in MATLAB) in repository:

```MATLAB
EROS_RODENTS.m
```
and
```MATLAB
PDLI_RODENTS.m
```
Both scripts accept .mat files structure exported by LabCart software. For EDF support, switch to the **"human"** branch.

####EROS_RODENTS.m

Calculates, visualizes, and stores Event Related Oscillations (ERO), Phase Locking Index (PLI), and Event Related Potential (ERP). Additional output is Averaged Event Related Oscillations (AVG_ERO), which is the ERO obtained from ERP.

####PDLI_RODENTS.m

Calculates, visualizes, and stores Phase Locking Index as a functional connectivity measure between electrodes. All combinations of electrodes are provided.

###Workflow

####Configure parallel computation:
For optimal performance set a number of parallel workers: Prallel->Manage Cluster Profiles->Cluster Profile->Edit->NumWorkers

For older MATLAB versions execute "matlabpool open" before calculation.

####Open your dataset:

![Open](https://github.com/VlastaKoudelka/EROs/blob/master/Doc/Open_rodents.png)

####Visualize:

![Visual](https://github.com/VlastaKoudelka/EROs/blob/master/Doc/Visual_rodents.png)

####Store your results:

Output data from **EROS_RODENTS.m** is the following array of subjects:
```MATLAB
subject(i)
```
ans = 
```MATLAB
          n_ch: 2                           %number of channels
    chan_label: {{1x1 cell}  {1x1 cell}}    %channel names
        f_name: 'P2_NT4T880 28MAY2014_2chan.mat'   %source file
        f_path: 'D:\Grisa\DPLI\'                   %source path
      triggers: [320x3 double]                     %event list
           ERO: {2x2 cell}                         %rows - channels
       AVG_ERO: {2x2 cell}                         %columns - events  
           ERP: {2x2 cell}                         %(1,2,3...) 
           PLI: {2x2 cell}                         %according to event list
```
Output data from **DPLI_RODENTS.m** is the following array of subjects:

```MATLAB
subject(1)
```
ans = 
```MATLAB
          n_ch: 2
      raw_data: {[1x5974400 double]  [1x5974400 double]}
    chan_label: {{1x1 cell}  {1x1 cell}}
        Fs_raw: 4000
        f_name: 'P2_NT4T880 28MAY2014_2chan.mat'
        f_path: 'D:\Grisa\DPLI\'
      triggers: [320x2 double]
        A_PDLI: {2x1 cell}              %DPLI for event A
        B_PDLI: {2x1 cell}              %DPLI for event B
             t: [1x326 double]
             f: [1x92 double]
```             
