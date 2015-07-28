function ERO_ST_BATCH
%                  _______  _______  _______  _______ 
%                 (  ____ \(  ____ )(  ___  )(  ____ \
%                 | (    \/| (    )|| (   ) || (    \/
%                 | (__    | (____)|| |   | || (_____ 
%                 |  __)   |     __)| |   | |(_____  )      (BATCH)
%                 | (      | (\ (   | |   | |      ) |        
%                 | (____/\| ) \ \__| (___) |/\____) |
%                 (_______/|/   \__/(_______)\_______)
%                                   
%  modified> 8.7.2015                         coded by> Vlastimil Koudelka
%                                       used code by>Robert Glenn Stockwell
%
close all
load(uigetfile('*.mat','Open the source file'),'com', 'data', 'samplerate')

t_pre = 100*1e-3;            %start trial before trigger [s]
t_post = 900*1e-3;           %stop trial after trigger [s]
delay = 0*1e-3;              %some delay of trigger flag [s]
f_res = 1;                   %desired resolution in spectogram [Hz]
f_max = 50;                 %maximum frequency in spectogram [Hz]
n_ch = 4;                    %#channels recorded in data vector

Fs = 250;                               %down-sampled 4kHz -> 250Hz        
T = 1/Fs;                               %sample period
n_pre = round(t_pre*Fs);                % #samples before trigger
n_post = round(t_post*Fs);              % #samples after trigger
n_delay = round(delay*Fs);              % #samples of delay
N = n_pre + n_post + 1;                 % #samples within a trial
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
for i = 1:n_ch                                        
    data{i} = filtfilt(Num,1,s{i});            %Zero phase filtering
    data{i} = downsample(data{i},4)';              %Fs 4kHz -> 1kHz
    data{i} = filtfilt(Num,1,data{i});            %Zero phase filtering
    data{i} = downsample(data{i},4)';              %Fs 1kHz -> 250Hz
end

%% Segmantation

for l = 1:n_ch                                          %over all channels     
    j = 1;
    k = 1;
    for i = 1:size(flags,1)                             %the first event
        if (flags(i,1) == 1)
            start_idx = flags(i,2) - n_delay - n_pre;   %begin idx. of trial
            stop_idx = flags(i,2) - n_delay + n_post;   %end idx. of trial
            A{l}(j,:) = data{l}(start_idx:stop_idx);    %trial segment   
            j = j + 1;        
        end

        if (flags(i,1) == 2)                            %the second event
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

for j = 1:n_ch    
    for i = 1:n_non_t 
        [A_ST{i,j},t,f] = st(A{j}(i,:),0,f_max,T, f_res);  %S-transformation
    end

    for i = 1:n_t
        [B_ST{i,j},t,f] = st(B{j}(i,:),0,f_max,T, f_res);
    end
end
%% Postprocessing
for l = 1:n_ch    
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
    A_rel_pow{i} = A_mean{i}.^2/max(max(A_mean{i}.^2));    %Relative spectral pow.
    B_rel_pow{i} = B_mean{i}.^2/max(max(B_mean{i}.^2));
end


%% Visualization
t = (t * (t_pre + t_post) - t_pre)*1e3;

%% 
figure
subplot(2,2,1)
contourf(t,f,A_rel_pow{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 1 - CH1 - non-target')

subplot(2,2,2)
contourf(t,f,B_rel_pow{1},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 1 - CH1 - target')

subplot(2,2,3)
ERP = mean(A{1},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 1 - CH1 - non-target')

subplot(2,2,4)
ERP = mean(B{1},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 1 - CH1 - target')

%% 
figure
subplot(2,2,1)
contourf(t,f,A_rel_pow{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 1 - CH2 - non-target')

subplot(2,2,2)
contourf(t,f,B_rel_pow{2},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 1 - CH2 - target')

subplot(2,2,3)
ERP = mean(A{2},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 1 - CH2 - non-target')

subplot(2,2,4)
ERP = mean(B{2},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 1 - CH2 - target')

%% 
figure
subplot(2,2,1)
contourf(t,f,A_rel_pow{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 2 - CH1 - non-target')

subplot(2,2,2)
contourf(t,f,B_rel_pow{3},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 2 - CH1 - target')

subplot(2,2,3)
ERP = mean(A{3},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 2 - CH1 - non-target')

subplot(2,2,4)
ERP = mean(B{3},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 2 - CH1 - target')
%% 
figure
subplot(2,2,1)
contourf(t,f,A_rel_pow{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 2 - CH2 - non-target')

subplot(2,2,2)
contourf(t,f,B_rel_pow{4},20,'LineStyle','none')
xlabel('time [ms]')
ylabel('frequency [Hz]')
title('Mouse 2 - CH2 - target')

subplot(2,2,3)
ERP = mean(A{4},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 2 - CH2 - non-target')

subplot(2,2,4)
ERP = mean(B{4},1);
plot(t,ERP)
xlabel('time [ms]')
ylabel('<voltage> [uV]')
title('Mouse 2 - CH2 - target')

save EROS.mat A_mean B_mean f t 
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
