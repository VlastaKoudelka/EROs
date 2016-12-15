#EROs, PLI, & DPLI Analyzer
##Authors:
* Vlastimil Koudelka
* Grygoriy Tsenov

This is for human, for animal version please swith to **"master"** branch.

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

![Open](https://github.com/VlastaKoudelka/EROs/blob/human/Doc/Open.png)

####Select multiple electrodes (use control+click or multiple selectin) for calculation:

![sel_el](https://github.com/VlastaKoudelka/EROs/blob/human/Doc/sel_el.png)

####Select triggers (target/non-target):

![sel_trg](https://github.com/VlastaKoudelka/EROs/blob/human/Doc/sel_trig.png)

####Select files to visualize, MEAN SUBJECT structure contains averages values of ERP, PLI, ERP over all files - the script gives names to figures and subplots automatically:

![sel_file](https://github.com/VlastaKoudelka/EROs/blob/human/Doc/sel_file.png)

####Vizualize:

![result](https://github.com/VlastaKoudelka/EROs/blob/human/Doc/result.png)

####Store your results:

Output date are stored into an array of structures:

```MATLAB
subject(1)
```

ans = 
```MATLAB
        f_name: 'Easrec_audio_ep-base_AM-v3_150616-0839.c0.edf'
        f_path: 'D:\EP_misa\adults\'
          n_ch: 1
    chan_label: {'O2'}
    trig_label: {'Trigger-31'  'Trigger-32'}
      triggers: [200x2 double]
           ERO: {2x1 cell}    %two cells, two triggers
           ERP: {2x1 cell}    
           PLI: {2x1 cell}
             f: [1x69 double]
```








