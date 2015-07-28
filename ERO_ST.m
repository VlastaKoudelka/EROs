function ERO_ST
%                  _______  _______  _______  _______ 
%                 (  ____ \(  ____ )(  ___  )(  ____ \
%                 | (    \/| (    )|| (   ) || (    \/
%                 | (__    | (____)|| |   | || (_____ 
%                 |  __)   |     __)| |   | |(_____  )
%                 | (      | (\ (   | |   | |      ) |
%                 | (____/\| ) \ \__| (___) |/\____) |
%                 (_______/|/   \__/(_______)\_______)
%                                   
%  modified> 8.7.2015                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
close all
load(uigetfile('*.mat','Open the source file'),'com', 'data', 'samplerate')

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
flags = load_flags(com);                % (trigger labels , sample idx.)

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

t = (t * (t_pre + t_post) - t_pre)*1e3;
A_rel_pow = A_mean.^2/max(max(A_mean.^2));    %Relative spectral pow.
B_rel_pow = B_mean.^2/max(max(B_mean.^2));

%% Visualization
figure
subplot(2,2,1)
contourf(t,f,A_rel_pow,20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('non-target')
% colormap(winter)
subplot(2,2,2)
contourf(t,f,B_rel_pow,20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('target')
% colormap(winter)

A_ERP = mean(A,1);
B_ERP = mean(B,1);
subplot(2,2,3)
plot(t,A_ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('non-target')
subplot(2,2,4)
plot(t,B_ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('target')

save EROS.mat A_mean B_mean A_ERP B_ERP f t 
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
