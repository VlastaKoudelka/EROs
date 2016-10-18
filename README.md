# EROS HUMAN
For animal version please swith to **"master"** branch.

##Experimental code for Event Realated Oscilation detection.

Developed software calculates advanced time-frequency quantities addressing evoked and induced EEG events. The software is designed to analyze both animal and human EEG data in order to provide translation between the observed phenomena. The analyzer is based on MATLAB platform and accepts LabChart and EDF data formats. More specifically, Event Related Oscillations are calculated to address induced oscillations. Phase Locking Index evaluates evoked oscillations over the trials and Phase Difference Locking Index measures functional connections between selected electrodes. The main outputs of the analyzer are the time-frequency characteristics of the quantities mentioned above.


ERO_HUMAN.m:

Is a piece of code for human eeg. The script reads EDF+ files with triggers.

* select multiple electrodes (use control+click or multiple selectin) the script evaluates all
* select two triggers (target/non-target) - the script operates on the fixed protocol
* select files to visualize - the script gives names to figures and subplots automatically
* MEAN SUBJECT structure contains averages values of ERP,PLI,ERP over all files




