function PDLI_RODENTS
% 
% __/\\\\\\\\\\\\\____/\\\\\\\\\\\\_____/\\\___________/\\\\\\\\\\\_        
%  _\/\\\/////////\\\_\/\\\////////\\\__\/\\\__________\/////\\\///__       
%   _\/\\\_______\/\\\_\/\\\______\//\\\_\/\\\______________\/\\\_____      
%    _\/\\\\\\\\\\\\\/__\/\\\_______\/\\\_\/\\\______________\/\\\_____     
%     _\/\\\/////////____\/\\\_______\/\\\_\/\\\______________\/\\\_____    
%      _\/\\\_____________\/\\\_______\/\\\_\/\\\______________\/\\\_____   
%       _\/\\\_____________\/\\\_______/\\\__\/\\\______________\/\\\_____  
%        _\/\\\_____________\/\\\\\\\\\\\\/___\/\\\\\\\\\\\\__/\\\\\\\\\\\_ 
%         _\///______________\////////////_____\////////////__\///////////_
% 
%modified> 24.6.2016                           coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
%
% The function calculates  Phase Difference Locking Index between EEG chan.

%% Batch execution
close all
[names,path] = uigetfile('*.mat','Open the source file','MultiSelect', 'on');

if ~(iscell(names))                     %in the case of one file
    names = {names};
end

subject = read_data(path,names);
tic
 
parfor i = 1:length(subject)
    [A_PDLI{i},B_PDLI{i},t{i},f{i}] = PDLI_CALC(subject(i));
end
for i = 1:length(subject)
    subject(i).A_PDLI = A_PDLI{i};
    subject(i).B_PDLI = B_PDLI{i};
    subject(i).t = t{i};
    subject(i).f = f{i};    
end
clear A_PDLI B_PDLI t f
subject(end + 1) = crt_mean_sbj(subject);
toc
pdli_vis(subject)
save PDLI_out subject
end

%% Core calculation
function [A_PDLI,B_PDLI,t,f] = PDLI_CALC(subject)
t_pre = 400*1e-3;            %start trial before trigger [s]
t_post = 900*1e-3;          %stop trial after trigger [s]
delay = 0;                   %some delay of trigger flag [s]
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
                                        
subject.triggers(:,2) = round(subject.triggers(:,2)/(subject.Fs_raw/Fs));  %down-sampled flags                                        

for i = 1:subject.n_ch
    if (Fs == 250)
    subject.raw_data{i} = filtfilt(Num,1,subject.raw_data{i});            %Zero phase filtering
    subject.raw_data{i} = downsample(subject.raw_data{i},4)';              %Fs 4kHz -> 1kHz
    end
    
    if (Fs == 250) || ((Fs == 1000))
    subject.raw_data{i} = filtfilt(Num,1,subject.raw_data{i});            %Zero phase filtering
    subject.raw_data{i} = downsample(subject.raw_data{i},4)';              %Fs 1kHz -> 250Hz
    end
    if ((Fs ~= 250) && (Fs ~= 1000) && (Fs ~= 4000))
        error('Wrong sampling rate, try 250Hz, 1000Hz, or 4000Hz!')
        break
    end
end

%% Calculation
if ((subject.triggers(end,2)- n_delay + n_post)>length(subject.raw_data{1}))
    subject.triggers = subject.triggers(1:end-1,:);  %^the last segment overflow 
end

A_cum = 0;
B_cum = 0;
for i = 1:size(subject.triggers,1)          %the first event
    if (subject.triggers(i,1) == 1)
        for m = 1:subject.n_ch              %over all channels
            for n = 1:m-1                   %below the diagonal         
                start_idx = subject.triggers(i,2) - n_delay - n_pre;   %begin idx. of trial
                stop_idx = subject.triggers(i,2) - n_delay + n_post;   %end idx. of trial
                if A_cum == 0
                    [S_m,t,f] = st_tuned(subject.raw_data{m}(start_idx:stop_idx),0,f_max,T, f_res);
                    [S_n,t,f] = st_tuned(subject.raw_data{n}(start_idx:stop_idx),0,f_max,T, f_res);
                    A_PDLI{m,n} = S_m./S_n.*abs(S_n)./abs(S_m);
                else
                    S_m = st_tuned(subject.raw_data{m}(start_idx:stop_idx),0,f_max,T, f_res);
                    S_n = st_tuned(subject.raw_data{n}(start_idx:stop_idx),0,f_max,T, f_res);
                    A_PDLI{m,n} = A_PDLI{m,n} + S_m./S_n.*abs(S_n)./abs(S_m);
                end
            end
        end                                                                                                                                        
        A_cum = A_cum + 1;        
    end    
    if (subject.triggers(i,1) == 2)         %the second event
        for m = 1:subject.n_ch              %over all channels
            for n = 1:m-1                   %below the diagonal         
                start_idx = subject.triggers(i,2) - n_delay - n_pre;   %begin idx. of trial
                stop_idx = subject.triggers(i,2) - n_delay + n_post;   %end idx. of trial
                if B_cum == 0
                    S_m = st_tuned(subject.raw_data{m}(start_idx:stop_idx),0,f_max,T, f_res);
                    S_n = st_tuned(subject.raw_data{n}(start_idx:stop_idx),0,f_max,T, f_res);
                    B_PDLI{m,n} = S_m./S_n.*abs(S_n)./abs(S_m);
                else
                    S_m = st_tuned(subject.raw_data{m}(start_idx:stop_idx),0,f_max,T, f_res);
                    S_n = st_tuned(subject.raw_data{n}(start_idx:stop_idx),0,f_max,T, f_res);
                    B_PDLI{m,n} = B_PDLI{m,n} + S_m./S_n.*abs(S_n)./abs(S_m);
                end
            end
        end                                                                                                                                        
        B_cum = B_cum + 1; 
    end
end

for i = 1:subject.n_ch          %because I didn't know a number of stims.
    for j = 1:i-1
        A_PDLI{m,n} = abs(A_PDLI{m,n}/A_cum);
        B_PDLI{m,n} = abs(B_PDLI{m,n}/B_cum);
    end
end
t = (t - t_pre)*1e3;       %Time axis [ms]
subject.raw_data = 'erased by DPLI_CALC()';
end

%% Data loading
function subject = read_data(path,names)

n_sig = 0;                              %a number of signals
for i = 1:length(names)                 %over all files
    load(fullfile(path,names{i}),'com', 'data','titles','samplerate');
    subject(i).n_ch = size(titles,1);
    raw_data = vec2mat(data,length(data)/subject(i).n_ch);
    for j = 1:size(raw_data,1)  %mat2cell
        subject(i).raw_data{j} = raw_data(j,:);
    end
    for j = 1:subject(i).n_ch
        subject(i).chan_label{j} = {titles(j,:)}; %from char array to string
    end
    subject(i).Fs_raw = samplerate(1);
    subject(i).f_name = names{i};
    subject(i).f_path = path;
    subject(i).triggers = load_flags(com);
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

function mean_sbj = crt_mean_sbj(subject)

n_sbj = length(subject);
n_ch = subject(1).n_ch;

for m = 1:n_ch
    for n = 1:m-1    
        for k = 2:n_sbj             %first subj is cummulative
            subject(1).A_PDLI{m,n} = subject(1).A_PDLI{m,n} + subject(k).A_PDLI{m,n};
            subject(1).B_PDLI{m,n} = subject(1).B_PDLI{m,n} + subject(k).B_PDLI{m,n};
        end
        subject(1).A_PDLI{m,n} = subject(1).A_PDLI{m,n} / n_sbj;
        subject(1).B_PDLI{m,n} = subject(1).B_PDLI{m,n} / n_sbj;
    end
end


mean_sbj = subject(1);
mean_sbj.f_name = 'MEAN SUBJECT';
mean_sbj.f_path = 'N/A';
end

function pdli_vis(subject)

for m = 1:subject(1).n_ch
    for n = 1:m-1
        figure
        subplot(1,2,1)
        contourf(subject(end).t,subject(end).f,subject(end).A_PDLI{m,n},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title('MEAN PDLI A event')
        subplot(1,2,2)
        contourf(subject(end).t,subject(end).f,subject(end).B_PDLI{m,n},20,'LineStyle','none')
        xlabel('time [ms]')
        ylabel('frequency [Hz]')
        title('MEAN PDLI B event')    
    end
end
end