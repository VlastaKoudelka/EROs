function EROS_HUMAN
%                  _______  _______  _______  _______ 
%                 (  ____ \(  ____ )(  ___  )(  ____ \
%                 | (    \/| (    )|| (   ) || (    \/
%                 | (__    | (____)|| |   | || (_____ 
%                 |  __)   |     __)| |   | |(_____  )    
%                 | (      | (\ (   | |   | |      ) |   
%                 | (____/\| ) \ \__| (___) |/\____) |
%                 (_______/|/   \__/(_______)\_______)
%                                   
% modified> 27.1.2016                          coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
% 

%% Batch execution

close all
[f_name,f_path] = uigetfile('*.edf','Open the source file','MultiSelect', 'on');
tic
subject = par_manage(f_path,f_name);
subject(end + 1) = crt_mean_sbj(subject);
save ROI_in subject
toc
visualize_eros(subject);

end

%% Parallel management - data loading, computation, and completion
function subject = par_manage(f_path,f_name)

if ~(iscell(f_name))                    %in the case of one file
    f_name = {f_name};
end

no_files = length(f_name);              %a number of loaded files

for i = 1:length(f_name);               %subject in a sequence
    
    [data,header] = ReadEDF(fullfile(f_path,f_name{i}));       %read EDF 
    subject(i).f_name = f_name{i}; 
    subject(i).f_path = f_path;
    if i == 1                   %for the first subject only
        %---user input---BEGIN
        [el_idx,ok] = listdlg('PromptString','Select electrodes:','ListString',header.labels);   %select electrodes
        if length(el_idx) > 20
            el_idx = el_idx(1:20);
        end
        [ref_idx,new_ref_flag] = listdlg('PromptString','Select a new ref. or cancel:','ListString',header.labels);   %select reference          

        event = unique(header.annotation.event);
        trig_idx = listdlg('PromptString','Select triggers:','ListString',event);   %select triggers
        %---user input---END
    end
   
    if new_ref_flag
        new_ref = data{ref_idx(1)}/length(ref_idx);
        for j = 2:length(ref_idx)
            new_ref = new_ref + data{ref_idx(j)}/length(ref_idx);
        end
    end
    
    subject(i).n_ch = length(el_idx);
    for j = 1:length(el_idx)
        if new_ref_flag
            data{j} = data{el_idx(j)} - new_ref;
        else
            data{j} = data{el_idx(j)};
        end
        subject(i).chan_label{j} = header.labels{el_idx(j)};
    end
    
    for j = 1:length(trig_idx)
        subject(i).trig_label{j} = event{trig_idx(j)};
    end
    
    m = 1;
    for k = 1:length(header.annotation.event)
       for j = 1:length(subject(i).trig_label)
            if strcmp(subject(i).trig_label{j},header.annotation.event{k})
                subject(i).triggers(m,1) = j;              %event type
                subject(i).triggers(m,2) = header.annotation.starttime(k);
                m = m + 1;
            end
        end
    end
    subject(i).triggers(:,2) = round(subject(i).triggers(:,2)*1e3);
    %time [s] -> time [ms] -> sample indices
    
    %---parallel computation---BEGIN
    parfor j = 1:length(data)  %channels in parallel
        [A_rpow_ERO{j}, B_rpow_ERO{j}, A_ERP{j}, B_ERP{j}, A_PLI{j}, B_PLI{j}, f{j}, t{j}] = EROS_CALC(data{j},subject(i).triggers);
    end
    %---parallel computation---END
    
    k = 1;
    for j = 1:length(subject(i).chan_label) %over all channels
        subject(i).ERO{1,j} = A_rpow_ERO{k};
        subject(i).ERO{2,j} = B_rpow_ERO{k};
        subject(i).ERP{1,j} = A_ERP{k};
        subject(i).ERP{2,j} = B_ERP{k};        
        subject(i).PLI{1,j} = A_PLI{k};
        subject(i).PLI{2,j} = B_PLI{k};
        subject(i).f = f{1};
        subject(i).t = t{1};
        k = k + 1;
    end
end
end
%% EROS calculation
function [A_rpow_ERO, B_rpow_ERO, A_ERP, B_ERP, A_PLI, B_PLI, f, t] = EROS_CALC(data, flags)
t_pre = 600*1e-3;            %start trial before trigger [s]
t_post = 1100*1e-3;          %stop trial after trigger [s]
delay = 0;                   %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 40;                  %maximum frequency in spectogram [Hz]

Fs = 250;                               %down-sampled 1kHz -> 250Hz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial

%% Prefiltering & Down-sampling
load filter_DP_40_50.mat Num            %the same anti-aliasing filter for:

                                        %Fs=1kHz, fp=40Hz, fs=50Hz
                                        
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 1kHz -> 250Hz

load filter_HP_0_2.mat Num

data = filtfilt(Num,1,data);            %high pass filter

flags(:,2) = round(flags(:,2)/4);       %4xundersampled

%% Segmantation @ transformation @ averaging (saves memory)
if ((flags(end,2)- n_delay + n_post)>length(data))
    flags = flags(1:end-1,:);                %^the last segment overflow 
end

A_ERP_cum = 0; B_ERP_cum = 0; A_ERO_cum = 0; B_ERO_cum = 0; 
A_PLI_cum = 0; B_PLI_cum = 0;
no_A = 0; no_B = 0;

for i = 1:size(flags,1)                             %the first event
    if (flags(i,1) == 1)
        start_idx = flags(i,2) - n_delay - n_pre;   %begin idx. of trial
        stop_idx = flags(i,2) - n_delay + n_post;   %end idx. of trial
        [A_ST,t,f] = st_tuned(data(start_idx:stop_idx),0,f_max,T, f_res);  %S-transformation  
        
        A_ERP_cum = A_ERP_cum + data(start_idx:stop_idx);   %cumulative
        A_ERO_cum = A_ERO_cum + abs(A_ST);                  %operation
        A_PLI_cum = A_PLI_cum + A_ST./abs(A_ST);            %for average
        no_A = no_A + 1;        
    end    
    if (flags(i,1) == 2)                            %the second event
        start_idx = flags(i,2) - n_delay - n_pre;
        stop_idx = flags(i,2) - n_delay + n_post;
        [B_ST,t,f] = st_tuned(data(start_idx:stop_idx),0,f_max,T, f_res);  %S-transformation  
        
        B_ERP_cum = B_ERP_cum + data(start_idx:stop_idx);   %cumulative
        B_ERO_cum = B_ERO_cum + abs(B_ST);                  %operation
        B_PLI_cum = B_PLI_cum + B_ST./abs(B_ST);            %for average                
        no_B = no_B + 1;  
    end
end
clear data;

A_ERP = A_ERP_cum/no_A;     %cumulated devided by number of events
A_ERO = A_ERO_cum/no_A;
A_PLI = abs(A_PLI_cum/no_A);

B_ERP = B_ERP_cum/no_B;     %cumulated devided by number of events
B_ERO = B_ERO_cum/no_B;
B_PLI = abs(B_PLI_cum/no_B);

%extra S-transformation for averaged ERP
[A_AVG_ERO,t,f] = st_tuned(A_ERP,0,f_max,T, f_res);  %S-transformation
[B_AVG_ERO,t,f] = st_tuned(B_ERP,0,f_max,T, f_res);  %for energy
A_AVG_ERO = abs(A_AVG_ERO);
B_AVG_ERO = abs(B_AVG_ERO);

%% Normalization
t = (t - t_pre)*1e3;       %Time axis [ms]

% Normalize to base line (average from bROI for each frequency)

t_bROI(1) = -500;        %start time for base line
t_bROI(2) = -200;        %stop time for base line
t_vis(1) = -200;         %start time of visualization (exclude artefacts)
t_vis(2) = 1000;          %stop time of visualization (exclude artefacts)

[c t_idx(1)] = min(abs(t-t_bROI(1)));
[c t_idx(2)] = min(abs(t-t_bROI(2)));
[c t_vis_idx(1)] = min(abs(t-t_vis(1)));
[c t_vis_idx(2)] = min(abs(t-t_vis(2)));

t = t(t_vis_idx(1):t_vis_idx(2));   %a new time axis for visualization

% ERO
A_ERO = A_ERO.^2;
B_ERO = B_ERO.^2;

base_pow{1} = mean(A_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{1}] = meshgrid(1:size(A_ERO,2),base_pow{1}); 
base_pow{2} = mean(B_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{2}] = meshgrid(1:size(B_ERO,2),base_pow{2}); 

A_rpow_ERO = (A_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{1}(:,t_vis_idx(1):t_vis_idx(2)))... 
            ./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2))*100;
B_rpow_ERO = (B_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{2}(:,t_vis_idx(1):t_vis_idx(2)))...
            ./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2))*100;

% ERP based ERO
A_AVG_ERO = A_AVG_ERO.^2;
B_AVG_ERO = B_AVG_ERO.^2;

base_pow{1} = mean(A_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{1}] = meshgrid(1:size(A_AVG_ERO,2),base_pow{1}); 
base_pow{2} = mean(B_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{2}] = meshgrid(1:size(B_AVG_ERO,2),base_pow{2});

A_rpow_AVG_ERO = (A_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{1}(:,t_vis_idx(1):t_vis_idx(2)))...
                ./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2))*100;
B_rpow_AVG_ERO = (B_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{2}(:,t_vis_idx(1):t_vis_idx(2)))...
                ./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2))*100;

B_ERP = B_ERP(:,t_vis_idx(1):t_vis_idx(2));           %short the time series correspondingly
A_ERP = A_ERP(:,t_vis_idx(1):t_vis_idx(2));

A_PLI = A_PLI(:,t_vis_idx(1):t_vis_idx(2));
B_PLI = B_PLI(:,t_vis_idx(1):t_vis_idx(2));

% Normalize to the maximal power

% A_rpow_ERO = A_ERO.^2/max(max(A_ERO.^2));    %Relative spectral pow.
% B_rpow_ERO = B_ERO.^2/max(max(B_ERO.^2));
% 
% A_rpow_AVG_ERO = A_AVG_ERO.^2/max(max(A_AVG_ERO.^2));    %Relative spectral pow. Averaged
% B_rpow_AVG_ERO = B_AVG_ERO.^2/max(max(B_AVG_ERO.^2));

end

%% Data loading
function subject = read_data(path,name)



end

%% Average
function mean_sbj = crt_mean_sbj(subject)

n_sbj = length(subject);
n_evt = length(subject(1).trig_label);
n_ch = subject(1).n_ch;

for i = 1:n_evt
    for j = 1:n_ch    
        for k = 2:n_sbj             %first subj is cummulative
            subject(1).ERO{i,j} = subject(1).ERO{i,j} + subject(k).ERO{i,j};
            subject(1).PLI{i,j} = subject(1).PLI{i,j} + subject(k).PLI{i,j};
            subject(1).ERP{i,j} = subject(1).ERP{i,j} + subject(k).ERP{i,j};
        end
        subject(1).ERO{i,j} = subject(1).ERO{i,j} / n_sbj;
        subject(1).PLI{i,j} = subject(1).PLI{i,j} / n_sbj;
        subject(1).ERP{i,j} = subject(1).ERP{i,j} / n_sbj;
    end            
end

mean_sbj = subject(1);
mean_sbj.f_name = 'MEAN SUBJECT';
mean_sbj.f_path = 'N/A';
mean_sbj.chan_label{end + 1} = 'all-ch-avrg';

% Average across all channels

for i = 1:n_evt         %initiate cumulative channel
    mean_sbj.ERO{i,n_ch + 1} = mean_sbj.ERO{i,1};
    mean_sbj.PLI{i,n_ch + 1} = mean_sbj.PLI{i,1};
    mean_sbj.ERP{i,n_ch + 1} = mean_sbj.ERP{i,1};
end

for i = 1:n_evt         %averaging
    for j = 2:n_ch 
        mean_sbj.ERO{i,n_ch + 1} = mean_sbj.ERO{i,n_ch + 1} + mean_sbj.ERO{i,j};
        mean_sbj.PLI{i,n_ch + 1} = mean_sbj.PLI{i,n_ch + 1} + mean_sbj.PLI{i,j};
        mean_sbj.ERP{i,n_ch + 1} = mean_sbj.ERP{i,n_ch + 1} + mean_sbj.ERP{i,j};
    end
    mean_sbj.ERO{i,n_ch + 1} = mean_sbj.ERO{i,n_ch + 1}/n_ch;
    mean_sbj.PLI{i,n_ch + 1} = mean_sbj.PLI{i,n_ch + 1}/n_ch;
    mean_sbj.ERP{i,n_ch + 1} = mean_sbj.ERP{i,n_ch + 1}/n_ch;
end
    
end

%% Visuaslization
function visualize_eros(subject)

for i = 1:length(subject)
    files{i} = subject(i).f_name;
end
f_idx = listdlg('PromptString','Select files for visualization:','ListString',files);

for i = 1:length(f_idx)
    for j = 1:subject(f_idx(i)).n_ch + (1*isequal(f_idx(i),length(subject))) %an exception for the mean subject
        figure('name',subject(f_idx(i)).f_name)
        subplot(2,3,1)
        contourf(subject(f_idx(i)).t,subject(f_idx(i)).f,subject(f_idx(i)).ERO{1,j},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title(['ENERGY > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{1}])
    
        subplot(2,3,4)
        contourf(subject(f_idx(i)).t,subject(f_idx(i)).f,subject(f_idx(i)).ERO{2,j},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title(['ENERGY > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{2}])
        
        subplot(2,3,2)
        contourf(subject(f_idx(i)).t,subject(f_idx(i)).f,subject(f_idx(i)).PLI{1,j},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title(['PLI > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{1}])
        
        subplot(2,3,5)
        contourf(subject(f_idx(i)).t,subject(f_idx(i)).f,subject(f_idx(i)).PLI{2,j},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title(['PLI > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{2}])
        
        subplot(2,3,3)
        plot(subject(f_idx(i)).t,subject(f_idx(i)).ERP{1,j})
        xlabel('time [ms]')
        ylabel('voltage [uV]')
        title(['ERP > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{1}])
        
        subplot(2,3,6)
        plot(subject(f_idx(i)).t,subject(f_idx(i)).ERP{2,j})
        xlabel('time [ms]')
        ylabel('voltage [uV]')
        title(['ERP > ch: ' subject(f_idx(i)).chan_label{j} ' > evt: ' subject(f_idx(i)).trig_label{2}])
    end                          
end
end


