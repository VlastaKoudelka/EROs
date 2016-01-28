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

