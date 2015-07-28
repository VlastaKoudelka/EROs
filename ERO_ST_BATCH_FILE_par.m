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
%  modified> 28.7.2015                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
%
% for older MATLAB versions execute "matlabpool open" command at first

%% Batch execution
close all
[names,paths] = uigetfile('*.mat','Open the source file','MultiSelect', 'on');

if ~(iscell(names))
    names = {names};
    paths = {paths};
end

for i = 1:length(names)                 %over all files
    load(fullfile(paths{i},names{i}),'com', 'data');
    [NOT_TARGET{i},TARGET{i}, f, t] = EROS_CALC(data, com);   
end

for i = 1:length(names)                 %subject - wise organization
    MOUSE_NOT_TARGET{2*i-1} = {NOT_TARGET{i}{1:2}};
    MOUSE_NOT_TARGET{2*i} = {NOT_TARGET{i}{3:4}};
    MOUSE_TARGET{2*i-1} = {TARGET{i}{1:2}};
    MOUSE_TARGET{2*i} = {TARGET{i}{3:4}};
end

clear TARGET NOT_TARGET

MEAN{1} = zeros(size(MOUSE_TARGET{1,1}{1,1}));  %CH1 target
MEAN{2} = zeros(size(MOUSE_TARGET{1,1}{1,1}));  %CH2 target
MEAN{3} = zeros(size(MOUSE_TARGET{1,1}{1,1}));  %CH1 not-target
MEAN{4} = zeros(size(MOUSE_TARGET{1,1}{1,1}));  %CH2 not_target
n_mice = length(MOUSE_TARGET);                  %# mice

for i = 1:n_mice                    %average over all mice
    MEAN{1} = MEAN{1} + MOUSE_TARGET{1,i}{1,1}/n_mice;
    MEAN{2} = MEAN{2} + MOUSE_TARGET{1,i}{1,2}/n_mice;
    MEAN{3} = MEAN{3} + MOUSE_NOT_TARGET{1,i}{1,1}/n_mice;
    MEAN{4} = MEAN{4} + MOUSE_NOT_TARGET{1,i}{1,2}/n_mice;
end  

visualize_eros(MEAN, f, t);
save ROI_in MEAN f t
end

%% EROS calculation
function [A_rel_pow, B_rel_pow, f, t] = EROS_CALC(data, com)
t_pre = 100*1e-3;            %start trial before trigger [s]
t_post = 900*1e-3;           %stop trial after trigger [s]
delay = 0*1e-3;              %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 50;                  %maximum frequency in spectogram [Hz]
n_ch = 4;                    %#channels recorded in data vector

Fs = 250;                               %down-sampled 4kHz -> 250Hz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial
n_f = round(f_max/f_res) + 1;           % #frequency samples in spectogram
flags = load_flags(com);                % (trigger labels , sample idx.)

%% Data to channels
b = round(length(data)/n_ch);                     %data contains 4 channels
s{1} = data(1,1:b);
s{2} = data(1,b+1:2*b);
s{3} = data(1,2*b+1:3*b);
s{4} = data(1,3*b+1:4*b);

clear data

%% Prefiltering & Down-sampling
load filters.mat Num                    %the same anti-aliasing filter for:
                                        %Fs=4kHz, fp=400Hz, fs=500Hz
                                        %Fs=1kHz, fp=100Hz, fs=125Hz   
                                        
Num = Num;                            %extra declaration for parallel comp.                                       
                                      
parfor i = 1:n_ch                                        
    data{i} = filtfilt(Num,1,s{i});            %Zero phase filtering
    data{i} = downsample(data{i},4)';          %Fs 4kHz -> 1kHz
    data{i} = filtfilt(Num,1,data{i});         %Zero phase filtering
    data{i} = downsample(data{i},4)';          %Fs 1kHz -> 250Hz
end

%% Segmantation

parfor l = 1:n_ch                            %over all channels in parallel    
    j = 1;
    k = 1;
    for i = 1:size(flags,1)                            %the first event
        if (flags(i,1) == 1)
            start_idx = flags(i,2) - n_delay - n_pre;  %begin idx. of trial
            stop_idx = flags(i,2) - n_delay + n_post;  %end idx. of trial
            A{l}(j,:) = data{l}(start_idx:stop_idx);   %trial segment   
            j = j + 1;        
        end

        if (flags(i,1) == 2)                           %the second event
            start_idx = flags(i,2) - n_delay - n_pre;
            stop_idx = flags(i,2) - n_delay + n_post;
            B{l}(k,:) = data{l}(start_idx:stop_idx);
            k = k + 1;  
        end
    end
end

%% Stockwell Transform
n_non_t = size(A{1},1);         %#non-target trials
n_t =  size(B{1},1);            %#target trials

parfor j = 1:n_ch               %parallel FOR execution
    for i = 1:n_non_t 
        [A_ST{i,j},t(j,:),f(j,:)] = st(A{j}(i,:),0,f_max,T, f_res);  %S-transform
    end
end

parfor j = 1:n_ch               %parallel FOR execution
    for i = 1:n_t
        [B_ST{i,j},t(j,:),f(j,:)] = st(B{j}(i,:),0,f_max,T, f_res);  %S-transform
    end
end
%% Postprocessing

parfor l = 1:n_ch    
    for i = 1:length(f)           %mean value calculation
        for j = 1:length(t)
            cum = 0;
            for k = 1:n_non_t
                cum = cum + abs(A_ST{k,l}(i,j));
            end
            A_mean{l}(i,j) = cum/n_non_t;
            cum = 0;
            for k = 1:n_t
                cum = cum + abs(B_ST{k,l}(i,j));
            end
            B_mean{l}(i,j) = cum/n_t;
        end
    end
end

for i = 1:n_ch
    A_rel_pow{i} = A_mean{i}.^2/max(max(A_mean{i}.^2));    
    B_rel_pow{i} = B_mean{i}.^2/max(max(B_mean{i}.^2));
end                                                %^Relative spectral pow.


t = t(1,:);
f = f(1,:);
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


