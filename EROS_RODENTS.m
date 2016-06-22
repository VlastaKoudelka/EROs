function EROS_RODENTS
%                  _______  _______  _______  _______ 
%                 (  ____ \(  ____ )(  ___  )(  ____ \
%                 | (    \/| (    )|| (   ) || (    \/
%                 | (__    | (____)|| |   | || (_____ 
%                 |  __)   |     __)| |   | |(_____  )    (BATCH_FILE)
%                 | (      | (\ (   | |   | |      ) |     (parallel)   
%                 | (____/\| ) \ \__| (___) |/\____) |
%                 (_______/|/   \__/(_______)\_______)
%                                   
%  modified> 22.6.2016                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
% 
% - for optimal performance set a number of parallel workers:
%    Prallel->Manage Cluster Profiles->Cluster Profile->Edit->NumWorkers
%
% - for older MATLAB versions execute "matlabpool open" before calculation


%% Batch execution
close all
[f_name,f_path] = uigetfile('*.mat','Open the source file','MultiSelect', 'on');
tic

subject = par_manage(f_path,f_name);
subject(end+1) = crt_mean_sbj(subject);
toc
visualize_eros(subject(end));
save ROI_in_new subject
end

%% Parallel computation - load, calculate in one function
function subject = par_manage(path,name)

if ~(iscell(name))                     %in the case of one file
    name = {name};
end

raw_data = {};       %will accumulate data from all subjects
triggers = {};       %will accumulate triggers from all subjects
for i = 1:length(name)                 %over all files
    load(fullfile(path,name{i}),'com', 'data','titles','samplerate');
    subject(i).n_ch = size(titles,1);
 
    for j = 1:subject(i).n_ch
        subject(i).chan_label{j} = {titles(j,:)}; %from char array to string
    end
    subject(i).Fs_raw = samplerate(1);                                                 
    subject(i).f_name = name{i};
    subject(i).f_path = path;
    subject(i).triggers = load_flags(com);
    
    data = vec2mat(data,length(data)/subject(i).n_ch); %channel per row
    raw_data = [raw_data; mat2cell(data,ones(1,subject(i).n_ch))]; %channel per cell
    
    for j = 1:subject(i).n_ch       %sliced triggers (one trigger table per ch.)
        triggers = [triggers; {subject(i).triggers}];
    end
    
end
clear data

parfor i = 1:length(raw_data)     %each channel one thread
    [A_ERO{i}, B_ERO{i}, A_AVG_ERO{i}, B_AVG_ERO{i}, A_ERP{i}, B_ERP{i}, A_PLI{i}, B_PLI{i}, f{i}, t{i}] = EROS_CALC(raw_data{i},triggers{i},subject(1).Fs_raw);
end

k = 1;
for i = 1:length(subject)     %sorting the outputs     
    for j = 1:subject(i).n_ch %over all cahnnels
        subject(i).ERO{1,j} = A_ERO{k}; 
        subject(i).ERO{2,j} = B_ERO{k};
        subject(i).AVG_ERO{1,j} = A_AVG_ERO{k}; 
        subject(i).AVG_ERO{2,j} = B_AVG_ERO{k};
        subject(i).ERP{1,j} = A_ERP{k}; 
        subject(i).ERP{2,j} = B_ERP{k};
        subject(i).PLI{1,j} = A_PLI{k}; 
        subject(i).PLI{2,j} = B_PLI{k};
        k = k + 1;
    end
    subject(i).t = t{i};
    subject(i).f = f{i}; 
end
end
%% EROS calculation
function [A_ERO, B_ERO, A_AVG_ERO, B_AVG_ERO, A_ERP, B_ERP, A_PLI, B_PLI, f, t] = EROS_CALC(data, flags,Fs_raw)
t_pre = 600*1e-3;            %start trial before trigger [s]
t_post = 1100*1e-3;          %stop trial after trigger [s]
delay = 0;                   %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 70;                  %maximum frequency in spectogram [Hz]

Fs = 4000;                               %down-sampled 4kHz -> 1kHz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial

%% Prefiltering & Down-sampling
load filters.mat Num                    %the same anti-aliasing filter for:
                                        %Fs=4kHz, fp=400Hz, fs=500Hz
                                        %Fs=1kHz, fp=100Hz, fs=125Hz
                                        
flags(:,2) = round(flags(:,2)/(Fs_raw/Fs));

if (Fs == 250)
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 4kHz -> 1kHz
end

if (Fs == 250) || ((Fs == 1000))
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 1kHz -> 250Hz
end

if ((Fs ~= 250) && (Fs ~= 1000) && (Fs ~= 4000))
    error('Wrong sampling rate, try 250Hz, 1000Hz, or 4000Hz!')
end

%% Segmantation
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

A_ERO = (A_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{1}(:,t_vis_idx(1):t_vis_idx(2)))...
        ./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2))*100;
        
B_ERO = (B_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{2}(:,t_vis_idx(1):t_vis_idx(2)))...
        ./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2))*100;
% ERP based ERO
A_AVG_ERO = A_AVG_ERO.^2;
B_AVG_ERO = B_AVG_ERO.^2;

base_pow{1} = mean(A_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{1}] = meshgrid(1:size(A_AVG_ERO,2),base_pow{1}); 
base_pow{2} = mean(B_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{2}] = meshgrid(1:size(B_AVG_ERO,2),base_pow{2});

A_AVG_ERO = (A_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{1}(:,t_vis_idx(1):t_vis_idx(2)))...
            ./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2))*100;
B_AVG_ERO = (B_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2)) - base_pow{2}(:,t_vis_idx(1):t_vis_idx(2)))...
            ./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2))*100;

A_ERP = A_ERP(t_vis_idx(1):t_vis_idx(2));           %shorten the time series correspondingly
B_ERP = B_ERP(t_vis_idx(1):t_vis_idx(2));

A_PLI = A_PLI(:,t_vis_idx(1):t_vis_idx(2));
B_PLI = B_PLI(:,t_vis_idx(1):t_vis_idx(2));

% Normalize to the maximal power

% A_rpow_ERO = A_ERO.^2/max(max(A_ERO.^2));    %Relative spectral pow.
% B_rpow_ERO = B_ERO.^2/max(max(B_ERO.^2));
% 
% A_rpow_AVG_ERO = A_AVG_ERO.^2/max(max(A_AVG_ERO.^2));    %Relative spectral pow. Averaged
% B_rpow_AVG_ERO = B_AVG_ERO.^2/max(max(B_AVG_ERO.^2));

end



%% Flag loading
function flags = load_flags(com)

for i = 1:size(com,1)/2
    temp(i,:) = com(2*i - 1:2*i,5);
end

if range(temp(:,1)) == 0    % > the first event means time tag
    time_flag_first = true;
else                        % > the second event means time tag
    time_flag_first = false;
end

if any(com(:,5) == 3)
    if time_flag_first == true;
        j = 1;
        for i = 1:size(com,1)
            if (com(i,5) == 2)              %non-target index
               flags(j,1) = 1;
               j = j + 1;
            end

            if (com(i,5) == 3)              %target index
               flags(j,1) = 2;
               j = j + 1;
            end

            if (com(i,5) == 1)              %sample index
               flags(j,2) = com(i,3);
            end
        end  
    else        
        j = 1;
        for i = 1:size(com,1)
            if (com(i,5) == 1)              %non-target index
               flags(j,1) = 1;
            end

            if (com(i,5) == 3)              %target index
               flags(j,1) = 2;
            end

            if (com(i,5) == 2)              %sample index
               flags(j,2) = com(i,3);
               j = j + 1;
            end
        end
    end
else
    flags(:,1) = com(:,5);
    flags(:,2) = com(:,3);   
end          
end

%% Visualization

function visualize_eros(subject)
%% EROS
figure
subplot(2,2,1)
contourf(subject.t,subject.f,subject.ERO{2,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH1 - target')

subplot(2,2,2)
contourf(subject.t,subject.f,subject.ERO{2,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH2 - target')

subplot(2,2,3)
contourf(subject.t,subject.f,subject.ERO{1,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH1 - non-target')

subplot(2,2,4)
contourf(subject.t,subject.f,subject.ERO{1,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH2 - non-target')

%% EROS AVERAGED (based on ERP)
figure
subplot(2,2,1)
contourf(subject.t,subject.f,subject.AVG_ERO{2,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH1 - target')

subplot(2,2,2)
contourf(subject.t,subject.f,subject.AVG_ERO{2,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH2 - target')

subplot(2,2,3)
contourf(subject.t,subject.f,subject.AVG_ERO{1,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH1 - non-target')

subplot(2,2,4)
contourf(subject.t,subject.f,subject.AVG_ERO{1,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH2 - non-target')

%% PLI
figure
subplot(2,2,1)
contourf(subject.t,subject.f,subject.PLI{2,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH1 - target')

subplot(2,2,2)
contourf(subject.t,subject.f,subject.PLI{2,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH2 - target')

subplot(2,2,3)
contourf(subject.t,subject.f,subject.PLI{1,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH1 - non-target')

subplot(2,2,4)
contourf(subject.t,subject.f,subject.PLI{1,2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH2 - non-target')

%% ERP
figure
subplot(2,2,1)
plot(subject.t,subject.ERP{2,1})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH1 - target')

subplot(2,2,2)
plot(subject.t,subject.ERP{2,2})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH2 - target')

subplot(2,2,3)
plot(subject.t,subject.ERP{1,1})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH1 - non-target')

subplot(2,2,4)
plot(subject.t,subject.ERP{1,2})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH2 - non-target')
end

%% Average
function mean_sbj = crt_mean_sbj(subject)

n_sbj = length(subject);
n_evt = 2;
n_ch = subject(1).n_ch;

for i = 1:n_evt
    for j = 1:n_ch    
        for k = 2:n_sbj             %first subj is cummulative
            subject(1).ERO{i,j} = subject(1).ERO{i,j} + subject(k).ERO{i,j};
            subject(1).PLI{i,j} = subject(1).PLI{i,j} + subject(k).PLI{i,j};
            subject(1).ERP{i,j} = subject(1).ERP{i,j} + subject(k).ERP{i,j};
            subject(1).AVG_ERO{i,j} = subject(1).AVG_ERO{i,j} + subject(k).AVG_ERO{i,j};
        end
        subject(1).ERO{i,j} = subject(1).ERO{i,j} / n_sbj;
        subject(1).PLI{i,j} = subject(1).PLI{i,j} / n_sbj;
        subject(1).ERP{i,j} = subject(1).ERP{i,j} / n_sbj;
        subject(1).AVG_ERO{i,j} = subject(1).AVG_ERO{i,j} / n_sbj;
    end            
end

mean_sbj = subject(1);
mean_sbj.f_name = 'MEAN SUBJECT';
mean_sbj.f_path = 'N/A';
    
end
