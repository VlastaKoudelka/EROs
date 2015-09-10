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
%  modified> 10.9.2015                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
% 

%% Batch execution

close all
[f_names,f_path] = uigetfile('*.edf','Open the source file','MultiSelect', 'on');
tic

if ~(iscell(f_names))                     %in the case of one file
    f_names = {f_names};
end

subject = read_data(f_path,f_names);

for i = 1:length(subject)   %over all subject
    for j = 1:length(subject.labels)
        [subject(i).ERO{1,j},subject(i).ERO{2,j},subject(i).ERP{1,j}, ...
         subject(i).ERP{2,j},subject(i).PLI{1,j},subject(i).PLI{2,j}, ...
         subject(i).f,subject(i).t] = EROS_CALC(subject(i).raw_data{i},subject.triggers);
    end
end

toc
visualize_eros(subject);
% save ROI_in subject
end

%% EROS calculation
function [A_rel_pow, B_rel_pow,T_ERP, NOT_T_ERP, B_PLI, A_PLI, f, t] = EROS_CALC(data, flags)
t_pre = 200*1e-3;            %start trial before trigger [s]
t_post = 1000*1e-3;          %stop trial after trigger [s]
delay = 0;                   %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 50;                  %maximum frequency in spectogram [Hz]

Fs = 250;                               %down-sampled 1kHz -> 250Hz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial

%% Prefiltering & Down-sampling
load filters.mat Num                    %the same anti-aliasing filter for:

                                        %Fs=1kHz, fp=100Hz, fs=125Hz
                                        
data = filtfilt(Num,1,data);            %Zero phase filtering
data = downsample(data,4)';              %Fs 1kHz -> 250Hz

flags(:,2) = round(flags(:,2)/4);

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
    [A_ST{i},t,f] = st(A(i,:),0,f_max,T, f_res);  %S-transformation
end

for i = 1:size(B,1)
    [B_ST{i},t,f] = st(B(i,:),0,f_max,T, f_res);
end

%% Postprocessing
for i = 1:size(A_ST{1},1)           %mean value calculation
    for j = 1:size(A_ST{1},2)
        cum_abs = 0;
        cum_phase = 0;
        for k = 1:length(A_ST)
            cum_abs   = cum_abs + abs(A_ST{k}(i,j));
            cum_phase = cum_phase + A_ST{k}(i,j)/abs(A_ST{k}(i,j));
        end
        A_mean(i,j) = cum_abs/length(A_ST);         %mean absolute value
        A_PLI(i,j)  = abs(cum_phase/length(A_ST));  %PLI
        cum_abs = 0;
        cum_phase = 0;
        for k = 1:length(B_ST)
            cum_abs   = cum_abs + abs(B_ST{k}(i,j));
            cum_phase = cum_phase + B_ST{k}(i,j)/abs(B_ST{k}(i,j));
        end
        B_mean(i,j) = cum_abs/length(B_ST);
        B_PLI(i,j)  = abs(cum_phase/length(B_ST)); 
    end
end

A_rel_pow = A_mean.^2/max(max(A_mean.^2));    %Relative spectral pow.
B_rel_pow = B_mean.^2/max(max(B_mean.^2));

t = (t * (t_pre + t_post) - t_pre)*1e3;       %Time axis

end

%% Data loading
function subject = read_data(path,names)

for i = 1:length(names)                 %over all files
    [data,header] = ReadEDF(fullfile(path,names{i}));       %read EDF 
    if i == 1    
        [el_idx,ok] = listdlg('PromptString','Select electrodes:','ListString',header.labels);   %select electrodes
        if length(el_idx) > 2
            el_idx = el_idx(1:2);
        end
        
        event = unique(header.annotation.event);
        trig_idx = listdlg('PromptString','Select triggers:','ListString',event);   %select triggers
    end
    
    for j = 1:length(el_idx)
        subject(i).raw_data{j} = data{el_idx(j)};
        subject(i).labels{j} = header.labels{el_idx(j)};
    end
    
    for j = 1:length(trig_idx)
        subject(i).trig_label{j} = event{trig_idx(j)};
    end
    subject(i).filename = fullfile(path,names{i});
    
    m = 1;
    for k = 1:length(header.annotation.event)
       for j = 1:length(subject(i).trig_label)
            if strcmp(subject(i).trig_label{j},header.annotation.event{k})
                subject(i).triggers(m,1) = j;                     %event type
                subject(i).triggers(m,2) = header.annotation.starttime(k);
                m = m + 1;
            end
        end
    end
    subject(i).triggers(:,2) = round(subject(i).triggers(:,2)*1e3);
    %time [s] -> time [ms] -> sample indices
end
end

%% Visuaslization

function visualize_eros(subject)
%% EROS
figure
subplot(2,2,1)
contourf(subject.t,subject.f,subject.ERO{1,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title(['ENERGY > CH1 - ' subject.trig_label{1}])

subplot(2,2,2)
contourf(subject.t,subject.f,subject.ERO{2,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title(['ENERGY > CH2 - ' subject.trig_label{2}])

subplot(2,2,3)
contourf(subject.t,subject.f,subject.PLI{1,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title(['ENERGY > CH1 - '  subject.trig_label{1}])

subplot(2,2,4)
contourf(subject.t,subject.f,subject.PLI{2,1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title(['ENERGY > CH2 - '  subject.trig_label{2}])
end


