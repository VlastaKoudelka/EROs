# EROs
Experimental code for Event Realated Oscilation detection

ERO_ST_BATCH_FILE_par.m:

Newly computes and stores EROs, ERPs, and PLIs:

MEAN_POW - mean energy value over all mice
MEAN_PLI - mean phase locking index over all mice

T_ERO 	- target EROs of all mice and channels
NOT_T_ERO 	- not-target EROs of all mice and channels

T_ERP 	- target evoked potencials of all mice and channels
NOT_T_ERP	- not-target evoked potencials of all mice and channels

T_PLI 	- target phase locking index of all mice and channels
NOT_T_PLI 	- not-target phase locking index of all mice and channels


This version contains validated batch parallel version of EROs calculation. Its validated against the previous versions of the algorithm.

This version accepts both options: one or two animals recorded in one input file.



EROS_ROI.m:

Still an experimental code for quantitative analyses of the evoked oscilations.



