function EROS_ROI
%% Quantitavive analyses of EROS 
%  modified> 17.7.2015                         coded by> Vlastimil Koudelka

load ROI_in MEAN f t

t1(1) = -50;        %start of the time interval (base ROI)
t2(1) = 0;          %end of the time interval (base ROI)
f1(1) = 8;          %start of the freq. interval
f2(1) = 35;          %end of the freq. interval

t1(2) = 0;          %start of the time interval (1st evoked)
t2(2) = 50;         %end of the time interval (1st evoked)
f1(2) = 8;          %start of the freq. interval
f2(2) = 35;          %end of the freq. interval

n_roi = 2;          %# ROIs

for i = 1:n_roi         %over all ROIs - find indices
    [c t1_idx(i)] = min(abs(t-t1(i)));
    [c t2_idx(i)] = min(abs(t-t2(i)));

    [c f1_idx(i)] = min(abs(f-f1(i)));
    [c f2_idx(i)] = min(abs(f-f2(i)));
end

%% Average power

for i = 1:n_roi     %over all ROIs
    for j = 1:4     %over all CHs & flags (target,not-target)
        ROI_mean(i,j) = mean2(MEAN{j}(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i)));
        ROI_std(i,j) =std2(MEAN{j}(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i)));
    end
end

ROI_mean    %mean values of power within a particular ROI (rows) and Channels(columns)
ROI_std     %standard deviation values of power within a particular ROI (rows) and Channels(columns)


% ROI1:  CH1 target, CH2 target, CH1 not-target, CH2 not-target
% ROI2:  CH1 target, CH2 target, CH1 not-target, CH2 not-target

