#EROS RODENTS

For human version swith to **"human"** branch, please. 

Developed software calculates advanced time-frequency quantities addressing evoked and induced EEG events. The software is designed to analyze both animal and human EEG data in order to provide translation between the observed phenomena. The analyzer is based on MATLAB platform and accepts LabChart and EDF data formats. More specifically, Event Related Oscillations are calculated to address induced oscillations. Phase Locking Index evaluates evoked oscillations over the trials and Phase Difference Locking Index measures functional connections between selected electrodes. The main outputs of the analyzer are the time-frequency characteristics of the quantities mentioned above.


output data is an array of subjects:

>> subject(i)

ans = 

          n_ch: 2                           %number of channels
    chan_label: {{1x1 cell}  {1x1 cell}}    %channel names
        f_name: 'P2_NT4T880 28MAY2014_2chan.mat'   %source file
        f_path: 'D:\Grisa\DPLI\'                   %source path
      triggers: [320x3 double]                     %event list
           ERO: {2x2 cell}                         %rows - channels
       AVG_ERO: {2x2 cell}                         %columns - events  
           ERP: {2x2 cell}                         %(1,2,3...) 
           PLI: {2x2 cell}                         %according to event list
    
AVG_ERO is ERO obtained from ERP signal.

