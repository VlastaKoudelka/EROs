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
%  modified> 30.9.2015                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
% 
% - for optimal performance set a number of parallel workers:
%    Prallel->Manage Cluster Profiles->Cluster Profile->Edit->NumWorkers
%
% - for older MATLAB versions execute "matlabpool open" before calculation


%% Batch execution
close all
[names,path] = uigetfile('*.mat','Open the source file','MultiSelect', 'on');
tic

if ~(iscell(names))                     %in the case of one file
    names = {names};
end

[raw_data, flags, labels] = read_data(path,names);

parfor i = 1:length(raw_data)   %over all channels
    [NOT_T_ERO{i},T_ERO{i}, NOT_T_ERP_ERO{i}, T_ERP_ERO{i}, NOT_T_ERP{i}, T_ERP{i}, NOT_T_PLI{i}, T_PLI{i}, f{i},t{i}] = EROS_CALC(raw_data{i},flags{i});
end

MEAN_POW{1} = zeros(size(NOT_T_ERO{1}));  %CH1 target
MEAN_POW{2} = zeros(size(NOT_T_ERO{1}));  %CH2 target
MEAN_POW{3} = zeros(size(NOT_T_ERO{1}));  %CH1 not-target
MEAN_POW{4} = zeros(size(NOT_T_ERO{1}));  %CH2 not_target
MEAN_ERP_POW{1} = zeros(size(NOT_T_ERO{1}));  %CH1 target
MEAN_ERP_POW{2} = zeros(size(NOT_T_ERO{1}));  %CH2 target
MEAN_ERP_POW{3} = zeros(size(NOT_T_ERO{1}));  %CH1 not-target
MEAN_ERP_POW{4} = zeros(size(NOT_T_ERO{1}));  %CH2 not_target
MEAN_PLI{1} = zeros(size(NOT_T_PLI{1}));  %CH1 target
MEAN_PLI{2} = zeros(size(NOT_T_PLI{1}));  %CH2 target
MEAN_PLI{3} = zeros(size(NOT_T_PLI{1}));  %CH1 not-target
MEAN_PLI{4} = zeros(size(NOT_T_PLI{1}));  %CH2 not_target
MEAN_ERP{1} = zeros(size(T_ERP{1}));  %CH1 target
MEAN_ERP{2} = zeros(size(T_ERP{1}));  %CH2 target
MEAN_ERP{3} = zeros(size(T_ERP{1}));  %CH1 not-target
MEAN_ERP{4} = zeros(size(T_ERP{1}));  %CH2 not_target


n_mice = length(raw_data)/2;
for i = 1:n_mice                    %average over all mice
    MEAN_POW{1} = MEAN_POW{1} + T_ERO{2*i - 1}/n_mice;     %CH1 target
    MEAN_POW{2} = MEAN_POW{2} + T_ERO{2*i}/n_mice;         %CH2 target
    MEAN_POW{3} = MEAN_POW{3} + NOT_T_ERO{2*i - 1}/n_mice; %CH1 not-target
    MEAN_POW{4} = MEAN_POW{4} + NOT_T_ERO{2*i}/n_mice;     %CH2 target
    
    MEAN_ERP_POW{1} = MEAN_ERP_POW{1} + T_ERP_ERO{2*i - 1}/n_mice;     %CH1 target
    MEAN_ERP_POW{2} = MEAN_ERP_POW{2} + T_ERP_ERO{2*i}/n_mice;         %CH2 target
    MEAN_ERP_POW{3} = MEAN_ERP_POW{3} + NOT_T_ERP_ERO{2*i - 1}/n_mice; %CH1 not-target
    MEAN_ERP_POW{4} = MEAN_ERP_POW{4} + NOT_T_ERP_ERO{2*i}/n_mice;     %CH2 target
    
    MEAN_PLI{1} = MEAN_PLI{1} + T_PLI{2*i - 1}/n_mice;     %CH1 target
    MEAN_PLI{2} = MEAN_PLI{2} + T_PLI{2*i}/n_mice;         %CH2 target
    MEAN_PLI{3} = MEAN_PLI{3} + NOT_T_PLI{2*i - 1}/n_mice; %CH1 not-target
    MEAN_PLI{4} = MEAN_PLI{4} + NOT_T_PLI{2*i}/n_mice;     %CH2 target
    
    MEAN_ERP{1} = MEAN_ERP{1} + T_ERP{2*i - 1}/n_mice;     %CH1 target
    MEAN_ERP{2} = MEAN_ERP{2} + T_ERP{2*i}/n_mice;         %CH2 target
    MEAN_ERP{3} = MEAN_ERP{3} + NOT_T_ERP{2*i - 1}/n_mice; %CH1 not-target
    MEAN_ERP{4} = MEAN_ERP{4} + NOT_T_ERP{2*i}/n_mice;     %CH2 target
end  
f = f{1};
t = t{1};
toc
visualize_eros(MEAN_POW, MEAN_ERP_POW, MEAN_PLI, MEAN_ERP, f, t);
save ROI_in MEAN_POW MEAN_ERP_POW MEAN_PLI MEAN_ERP T_ERO NOT_T_ERO T_ERP_ERO NOT_T_ERP_ERO T_ERP NOT_T_ERP T_PLI NOT_T_PLI labels f t
end

%% EROS calculation
function [A_rpow_ERO, B_rpow_ERO, A_rpow_AVG_ERO, B_rpow_AVG_ERO, NOT_T_ERP, T_ERP, A_PLI, B_PLI, f, t] = EROS_CALC(data, flags)
t_pre = 600*1e-3;            %start trial before trigger [s]
t_post = 1100*1e-3;          %stop trial after trigger [s]
delay = flags(1,3);          %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 70;                  %maximum frequency in spectogram [Hz]

Fs = 250;                               %down-sampled 4kHz -> 250Hz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial

%% Prefiltering & Down-sampling
load filters.mat Num                    %the same anti-aliasing filter for:
                                        %Fs=4kHz, fp=400Hz, fs=500Hz
                                        %Fs=1kHz, fp=100Hz, fs=125Hz
                                        
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 4kHz -> 1kHz
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 1kHz -> 250Hz

%% Segmantation
if ((flags(end,2)- n_delay + n_post)>length(data))
    flags = flags(1:end-1,:);                %^the last segment overflow 
end

j = 1;
k = 1;
for i = 1:size(flags,1)                             %the first event
    if (flags(i,1) == 1)
        start_idx = flags(i,2) - n_delay - n_pre;   %begin idx. of trial
        stop_idx = flags(i,2) - n_delay + n_post;   %end idx. of trial
        A(j,:) = data(start_idx:stop_idx);          %trial segment   
        j = j + 1;        
    end
    
    if (flags(i,1) == 2)                            %the second event
        start_idx = flags(i,2) - n_delay - n_pre;
        stop_idx = flags(i,2) - n_delay + n_post;
        B(k,:) = data(start_idx:stop_idx);
        k = k + 1;  
    end
end

%% Event related potencials
T_ERP = mean(B,1);              %target 
NOT_T_ERP = mean(A,1);          %not-target

%% Event related oscillations: Stockwell Transform 
for i = 1:size(A,1) 
    [A_ST{i},t,f] = st_tuned(A(i,:),0,f_max,T, f_res);  %S-transformation
end

for i = 1:size(B,1)
    [B_ST{i},t,f] = st_tuned(B(i,:),0,f_max,T, f_res);
end

[A_AVG_ERO,t,f] = st_tuned(NOT_T_ERP,0,f_max,T, f_res);  %S-transformation
[B_AVG_ERO,t,f] = st_tuned(T_ERP,0,f_max,T, f_res);      %for energy

A_AVG_ERO = abs(A_AVG_ERO);
B_AVG_ERO = abs(B_AVG_ERO);

%% Postprocessing
for i = 1:size(A_ST{1},1)           %mean value calculation
    for j = 1:size(A_ST{1},2)
        cum_abs = 0;
        cum_phase = 0;
        for k = 1:length(A_ST)
            cum_abs   = cum_abs + abs(A_ST{k}(i,j));
            cum_phase = cum_phase + A_ST{k}(i,j)/abs(A_ST{k}(i,j));
        end
        A_ERO(i,j) = cum_abs/length(A_ST);         %mean absolute value
        A_PLI(i,j)  = abs(cum_phase/length(A_ST));  %PLI
        cum_abs = 0;
        cum_phase = 0;
        for k = 1:length(B_ST)
            cum_abs   = cum_abs + abs(B_ST{k}(i,j));
            cum_phase = cum_phase + B_ST{k}(i,j)/abs(B_ST{k}(i,j));
        end
        B_ERO(i,j) = cum_abs/length(B_ST);
        B_PLI(i,j)  = abs(cum_phase/length(B_ST)); 
    end
end

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

A_rpow_ERO = A_ERO(:,t_vis_idx(1):t_vis_idx(2))./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2));
B_rpow_ERO = B_ERO(:,t_vis_idx(1):t_vis_idx(2))./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2));

% ERP based ERO
A_AVG_ERO = A_AVG_ERO.^2;
B_AVG_ERO = B_AVG_ERO.^2;

base_pow{1} = mean(A_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{1}] = meshgrid(1:size(A_AVG_ERO,2),base_pow{1}); 
base_pow{2} = mean(B_AVG_ERO(:,t_idx(1):t_idx(2)),2);
[c, base_pow{2}] = meshgrid(1:size(B_AVG_ERO,2),base_pow{2});

A_rpow_AVG_ERO = A_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2))./base_pow{1}(:,t_vis_idx(1):t_vis_idx(2));
B_rpow_AVG_ERO = B_AVG_ERO(:,t_vis_idx(1):t_vis_idx(2))./base_pow{2}(:,t_vis_idx(1):t_vis_idx(2));

T_ERP = T_ERP(:,t_vis_idx(1):t_vis_idx(2));           %short the time series correspondingly
NOT_T_ERP = NOT_T_ERP(:,t_vis_idx(1):t_vis_idx(2));

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
function [raw_data, flags, labels] = read_data(path,names)

n_sig = 0;                              %a number of signals
for i = 1:length(names)                 %over all files
    load(fullfile(path,names{i}),'com', 'data','titles');
    if (size(titles,1) == 2)
        raw_data{n_sig + 1} = data(1:length(data)/2);
        flags{n_sig + 1} = load_flags(com);
        labels{n_sig + 1} = fullfile(path,names{i});
        raw_data{n_sig + 2} = data(length(data)/2 + 1:end);
        flags{n_sig + 2} = load_flags(com); 
        labels{n_sig + 2} = fullfile(path,names{i});
        n_sig = n_sig + 2;
    else
        raw_data{n_sig + 1} = data(1:length(data)/4);
        flags{n_sig + 1} = load_flags(com);
        labels{n_sig + 1} = fullfile(path,names{i});
        raw_data{n_sig + 2} = data(length(data)/4 + 1:2*length(data)/4);
        flags{n_sig + 2} = load_flags(com);
        labels{n_sig + 2} = fullfile(path,names{i});
        raw_data{n_sig + 3} = data(2*length(data)/4 + 1:3*length(data)/4);
        flags{n_sig + 3} = load_flags(com);
        labels{n_sig + 3} = fullfile(path,names{i});
        raw_data{n_sig + 4} = data(3*length(data)/4 + 1:end);
        flags{n_sig + 4} = load_flags(com);
        labels{n_sig + 4} = fullfile(path,names{i});
        n_sig = n_sig + 4;
    end
end
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
        flags(1,3) = 0;                     %delay [s]    
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
        flags(1,3) = 0;                     %delay [s]
    end
else
    flags(:,1) = com(:,5);
    flags(:,2) = com(:,3);
    flags(1,3) = 70*1e-3;                    %delay [s] (uncorrected data)   
end
    
        

flags(:,2) = round(flags(:,2)/16);  %down-sampled
end

%% Visualization

function visualize_eros(ERO_vis, AVG_ERO_vis, PLI_vis, ERP_vis, f, t)
%% EROS
figure
subplot(2,2,1)
contourf(t,f,ERO_vis{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH1 - target')

subplot(2,2,2)
contourf(t,f,ERO_vis{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH2 - target')

subplot(2,2,3)
contourf(t,f,ERO_vis{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH1 - non-target')

subplot(2,2,4)
contourf(t,f,ERO_vis{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ENERGY > CH2 - non-target')

%% EROS AVERAGED (based on ERP)
figure
subplot(2,2,1)
contourf(t,f,AVG_ERO_vis{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH1 - target')

subplot(2,2,2)
contourf(t,f,AVG_ERO_vis{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH2 - target')

subplot(2,2,3)
contourf(t,f,AVG_ERO_vis{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH1 - non-target')

subplot(2,2,4)
contourf(t,f,AVG_ERO_vis{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('ERP ENERGY > CH2 - non-target')

%% PLI
figure
subplot(2,2,1)
contourf(t,f,PLI_vis{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH1 - target')

subplot(2,2,2)
contourf(t,f,PLI_vis{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH2 - target')

subplot(2,2,3)
contourf(t,f,PLI_vis{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH1 - non-target')

subplot(2,2,4)
contourf(t,f,PLI_vis{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('PLI > CH2 - non-target')

%% ERP
figure
subplot(2,2,1)
plot(t,ERP_vis{1})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH1 - target')

subplot(2,2,2)
plot(t,ERP_vis{2})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH2 - target')

subplot(2,2,3)
plot(t,ERP_vis{3})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH1 - non-target')

subplot(2,2,4)
plot(t,ERP_vis{4})
xlabel('time [ms]')
ylabel('voltage [uV]')
title('ERP > CH2 - non-target')
end


