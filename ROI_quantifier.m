function ROI_quantifier

%              _____   ____ _____             ____  
%              |  __ \ / __ \_   _|           / __ \ 
%              | |__) | |  | || |    ______  | |  | |
%              |  _  /| |  | || |   |______| | |  | |
%              | | \ \| |__| || |_           | |__| |
%              |_|  \_\\____/_____|           \___\_\
%                                        
% 
% 
%  modified> 22.9.2015                         coded by> Vlastimil Koudelka

%% Load data
[f_name,f_path] = uigetfile('*.mat','Open the source file','MultiSelect', 'off');
load(fullfile(f_path,f_name),'subject');

%% ROI definitions %%

% DELTA % - for other frequencies change f1 and f2 in lines 68-86
n_roi = 4;          % No of ROIs
t = subject(1).t;
f = subject(1).f;

        %BASAL%
f1(1) = 1;          %start of the freq. interval
f2(1) = 5;          %end of the freq. interval
t1(1) = -200;        %start of the time basal interval
t2(1) = 0;          %end of the time basal interval
        %P1/N1%
f1(2) = 1;          %start of the freq. interval
f2(2) = 5;          %end of the freq. interval
t1(2) = 0;          %start of the time interval p1-n1
t2(2) = 50;         %end of the time interval p1-n1
         %N2/P3%
f1(3) = 1;          %start of the freq. interval
f2(3) = 5;          %end of the freq. interval
t1(3) = 50;         %start of the time interval n2-p3
t2(3) = 120;        %end of the time interval n2-p3
        %P3%
f1(4) = 1;          %start of the freq. interval
f2(4) = 5;          %end of the freq. interval
t1(4) = 120;         %start of the time interval p3
t2(4) = 260;        %end of the time interval p3

for i = 1:n_roi         %over all ROIs - find indexes based on "t" and "f" variables
    [c t1_idx(i)] = min(abs(t-t1(i)));
    [c t2_idx(i)] = min(abs(t-t2(i)));

    [c f1_idx(i)] = min(abs(f-f1(i)));
    [c f2_idx(i)] = min(abs(f-f2(i)));
end

%% ROI selection

for i = 1:n_roi                                     %over all ROIs
    for j = 1:length(subject)                       %over all subjects
        for k = 1:subject(j).n_ch                   %over all channels
            for l = 1:length(subject(j).trig_label) %over all events
                subject(j).ERO_ROI{l,k,i} = subject(j).ERO{l,k}(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
                subject(j).PLI_ROI{l,k,i} = subject(j).PLI{l,k}(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
                subject(j).ERO_Q{l,k,i} = sum(sum(subject(j).ERO_ROI{l,k,i}));
                subject(j).PLI_Q{l,k,i} = sum(sum(subject(j).PLI_ROI{l,k,i}));
            end
        end
    end
end

%%

save ROI_out subject