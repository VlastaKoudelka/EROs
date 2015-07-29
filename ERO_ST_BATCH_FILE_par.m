function ERO_ST_BATCH_FILE_par
%                  _______  _______  _______  _______ 
%                 (  ____ \(  ____ )(  ___  )(  ____ \
%                 | (    \/| (    )|| (   ) || (    \/
%                 | (__    | (____)|| |   | || (_____ 
%                 |  __)   |     __)| |   | |(_____  )    (BATCH_FILE)
%                 | (      | (\ (   | |   | |      ) |     (parallel)   
%                 | (____/\| ) \ \__| (___) |/\____) |
%                 (_______/|/   \__/(_______)\_______)
%                                   
%  modified> 29.7.2015                         coded by> Vlastimil Koudelka
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
if ~(iscell(names))
    names = {names};
end

n_ch = 0;
for i = 1:length(names)                 %over all files
    load(fullfile(path,names{i}),'com', 'data','titles');
    if (size(titles,1) == 2)
        raw_data{n_ch + 1} = data(1:length(data)/2);
        flags{n_ch + 1} = load_flags(com); 
        raw_data{n_ch + 2} = data(length(data)/2 + 1:end);
        flags{n_ch + 2} = load_flags(com); 
        n_ch = n_ch + 2;
    else
        raw_data{n_ch + 1} = data(1:length(data)/4);
        flags{n_ch + 1} = load_flags(com);
        raw_data{n_ch + 2} = data(length(data)/4 + 1:2*length(data)/4);
        flags{n_ch + 2} = load_flags(com);
        raw_data{n_ch + 3} = data(2*length(data)/4 + 1:3*length(data)/4);
        flags{n_ch + 3} = load_flags(com);
        raw_data{n_ch + 4} = data(3*length(data)/4 + 1:end);
        flags{n_ch + 4} = load_flags(com);
        n_ch = n_ch + 4;
    end
end

parfor i = 1:length(raw_data)   %over all channels
    [NOT_TARGET{i},TARGET{i},f{i},t{i}] = EROS_CALC(raw_data{i},flags{i});
end

MEAN{1} = zeros(size(NOT_TARGET{1}));  %CH1 target
MEAN{2} = zeros(size(NOT_TARGET{1}));  %CH2 target
MEAN{3} = zeros(size(NOT_TARGET{1}));  %CH1 not-target
MEAN{4} = zeros(size(NOT_TARGET{1}));  %CH2 not_target

n_mice = n_ch/2;
for i = 1:n_mice                    %average over all mice
    MEAN{1} = MEAN{1} + TARGET{2*i - 1}/n_mice;
    MEAN{2} = MEAN{2} + TARGET{2*i}/n_mice;
    MEAN{3} = MEAN{3} + NOT_TARGET{2*i - 1}/n_mice;
    MEAN{4} = MEAN{4} + NOT_TARGET{2*i}/n_mice;
end  
f = f{1};
t = t{1};
toc
visualize_eros(MEAN, f, t);
save ROI_in MEAN f t
end

%% EROS calculation
function [A_rel_pow, B_rel_pow, f, t] = EROS_CALC(data, flags)
t_pre = 100*1e-3;            %start trial before trigger [s]
t_post = 900*1e-3;           %stop trial after trigger [s]
delay = 0*1e-3;              %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 50;                 %maximum frequency in spectogram [Hz]

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

%% Stockwell Transform
for i = 1:size(A,1) 
    [A_ST{i},t,f] = st(A(i,:),0,f_max,T, f_res);  %S-transformation
end

for i = 1:size(B,1)
    [B_ST{i},t,f] = st(B(i,:),0,f_max,T, f_res);
end

%% Postprocessing
for i = 1:size(A_ST{1},1)           %mean value calculation
    for j = 1:size(A_ST{1},2)
        cum = 0;
        for k = 1:length(A_ST)
            cum = cum + abs(A_ST{k}(i,j));
        end
        A_mean(i,j) = cum/length(A_ST);
        cum = 0;
        for k = 1:length(B_ST)
            cum = cum + abs(B_ST{k}(i,j));
        end
        B_mean(i,j) = cum/length(B_ST);
    end
end

A_rel_pow = A_mean.^2/max(max(A_mean.^2));    %Relative spectral pow.
B_rel_pow = B_mean.^2/max(max(B_mean.^2));

t = (t * (t_pre + t_post) - t_pre)*1e3;
end

%% Flag loading
function flags = load_flags(com)

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
flags(:,2) = round(flags(:,2)/16);  %down-sampled
end

%% Visuaslization

function visualize_eros(MEAN_vis, f, t)

figure
subplot(2,2,1)
contourf(t,f,MEAN_vis{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('CH1 - target')

subplot(2,2,2)
contourf(t,f,MEAN_vis{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('CH2 - target')

subplot(2,2,3)
contourf(t,f,MEAN_vis{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('CH1 - non-target')

subplot(2,2,4)
contourf(t,f,MEAN_vis{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('CH2 - non-target')
end


