%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   create ROI output for TOTAL POWER (ERO), EVOKED POWER (ERP-ERO), PhaseLockedIndex (PLI)   %
%                                   for Ch1 and Ch2 of TARGET and NOT-TARGET stimuli                                               %
%                       NEW VERSION                                                                                                                                                    %
%                                                                           modified 24-Feb-2016                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all; clear all; 
load ROI_in_bkgC57.mat %load your ROI_in.mat 

n_sbj = length(subject) -1;
%% DELTA ROI definitions %%

n_roi = 4;          % No of ROIs
        %BASAL%
f1(1) = 1;          %start of the freq. interval
f2(1) = 4;          %end of the freq. interval
t1(1) = -150;        %start of the time basal interval
t2(1) = -50;          %end of the time basal interval
        %P1/N1%
f1(2) = 1;          %start of the freq. interval
f2(2) = 4;          %end of the freq. interval
t1(2) = 0;          %start of the time interval p1-n1
t2(2) = 50;         %end of the time interval p1-n1
         %N2/P3%
f1(3) = 1;          %start of the freq. interval
f2(3) = 4;          %end of the freq. interval
t1(3) = 50;         %start of the time interval n2-p3
t2(3) = 250;        %end of the time interval n2-p3
        %P3%
f1(4) = 1;          %start of the freq. interval
f2(4) = 4;          %end of the freq. interval
t1(4) = 120;         %start of the time interval p3
t2(4) = 260;        %end of the time interval p3

for i = 1:n_roi         %over all ROIs - find indexes based on "t" and "f" variables
    [c t1_idx(i)] = min(abs(subject(1).t-t1(i)));
    [c t2_idx(i)] = min(abs(subject(1).t-t2(i)));

    [c f1_idx(i)] = min(abs(subject(1).f-f1(i)));
    [c f2_idx(i)] = min(abs(subject(1).f-f2(i)));
end

clear c f1 f2 i t1 t2;

% ERO %
for i=1:n_roi;
    for k= 1:n_sbj
        temp=subject(k).ERO{1,1};
        results.Ch1.ERO_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj
        temp=subject(k).ERO{2,1};
        results.Ch1.ERO_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj
        temp=subject(k).ERO{1,2};
        results.Ch2.ERO_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj
        temp=subject(k).ERO{2,2};
        results.Ch2.ERO_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% AVG_ERO %
for i=1:n_roi;
    for k= 1:n_sbj
        temp=subject(k).AVG_ERO{1,1};
        results.Ch1.AVG_ERO_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,1};
        results.Ch1.AVG_ERO_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,2};
        results.Ch2.AVG_ERO_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,2};
        results.Ch2.AVG_ERO_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% PLI %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,1};
        results.Ch1.PLI_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,1};
        results.Ch1.PLI_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,2};
        results.Ch2.PLI_Delta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,2};
        results.Ch2.PLI_Delta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% DELTA ROI MEAN calculations %
for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Delta_ERO_T(i,k) = mean2( results.Ch1.ERO_Delta_T{i,k}); %for target
        results.Ch1.mean_Delta_AVG_ERO_T(i,k) = mean2(results.Ch1.AVG_ERO_Delta_T{i,k}); 
        results.Ch1.mean_Delta_PLI_T(i,k) = mean2(results.Ch1.PLI_Delta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Delta_ERO_T(i,k) = mean2( results.Ch2.ERO_Delta_T{i,k}); %for target
        results.Ch2.mean_Delta_AVG_ERO_T(i,k) = mean2(results.Ch2.AVG_ERO_Delta_T{i,k}); 
        results.Ch2.mean_Delta_PLI_T(i,k) = mean2(results.Ch2.PLI_Delta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Delta_ERO_NT(i,k) = mean2( results.Ch1.ERO_Delta_NT{i,k}); %for target
        results.Ch1.mean_Delta_AVG_ERO_NT(i,k) = mean2(results.Ch1.AVG_ERO_Delta_NT{i,k}); 
        results.Ch1.mean_Delta_PLI_NT(i,k) = mean2(results.Ch1.PLI_Delta_NT{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Delta_ERO_NT(i,k) = mean2( results.Ch2.ERO_Delta_NT{i,k}); %for target
        results.Ch2.mean_Delta_AVG_ERO_NT(i,k) = mean2(results.Ch2.AVG_ERO_Delta_NT{i,k}); 
        results.Ch2.mean_Delta_PLI_NT(i,k) = mean2(results.Ch2.PLI_Delta_NT{i,k}); 
    end
end
clear i k;

% AVERAGE OVER ALL MICE Target and Non-target %
for i=1:n_roi;
results.Ch1.mean_Delta(i,1)=mean2(results.Ch1.mean_Delta_ERO_T(i,:)); %target ERO Ch1
results.Ch2.mean_Delta(i,1)=mean2(results.Ch2.mean_Delta_ERO_T(i,:)); %target ERO Ch2

results.Ch1.mean_Delta(i,2)=mean2(results.Ch1.mean_Delta_ERO_NT(i,:)); % not-target ERO Ch1
results.Ch2.mean_Delta(i,2)=mean2(results.Ch2.mean_Delta_ERO_NT(i,:)); % not-target ERO Ch2

results.Ch1.mean_Delta(i,3)=mean2(results.Ch1.mean_Delta_AVG_ERO_T(i,:)); %target ERP_ERO Ch1
results.Ch2.mean_Delta(i,3)=mean2(results.Ch2.mean_Delta_AVG_ERO_T(i,:)); %target ERP_ERO Ch2

results.Ch1.mean_Delta(i,4)=mean2(results.Ch1.mean_Delta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch1
results.Ch2.mean_Delta(i,4)=mean2(results.Ch2.mean_Delta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch2

results.Ch1.mean_Delta(i,5)=mean2(results.Ch1.mean_Delta_PLI_T(i,:)); %target PLI Ch1
results.Ch2.mean_Delta(i,5)=mean2(results.Ch2.mean_Delta_PLI_T(i,:)); %target PLI Ch2

results.Ch1.mean_Delta(i,6)=mean2(results.Ch1.mean_Delta_PLI_NT(i,:)); % not-target PLI Ch1
results.Ch2.mean_Delta(i,6)=mean2(results.Ch2.mean_Delta_PLI_NT(i,:)); % not-target PLI Ch2
end
clear i f1_idx f2_idx

% BARS PLOT %
figure; 
subplot (2,3,1); bar(results.Ch1.mean_Delta(:,1));
title('Delta ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch1.mean_Delta(:,2));
title('Delta ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch1.mean_Delta(:,3));
title('Delta ERP-ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch1.mean_Delta(:,4));
title('Delta ERP-ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch1.mean_Delta(:,5));
title('Delta PLI > Ch1 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch1.mean_Delta(:,6));
title('Delta PLI > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

figure; 
subplot (2,3,1); bar(results.Ch2.mean_Delta(:,1));
title('Delta ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch2.mean_Delta(:,2));
title('Delta ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch2.mean_Delta(:,3));
title('Delta ERP-ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch2.mean_Delta(:,4));
title('Delta ERP-ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch2.mean_Delta(:,5));
title('Delta PLI > Ch2 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch2.mean_Delta(:,6));
title('Delta PLI > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

%% THETA ROI definitions and calculation %%
        %BASAL%
f1(1) = 4;          %start of the freq. interval
f2(1) = 8;          %end of the freq. interval
        %P1/N1%
f1(2) = 4;          %start of the freq. interval
f2(2) = 8;          %end of the freq. interval
         %N2/P3%
f1(3) = 4;          %start of the freq. interval
f2(3) = 8;          %end of the freq. interval
        %P3%
f1(4) = 4;          %start of the freq. interval
f2(4) = 8;          %end of the freq. interval

for i = 1:n_roi         %over all ROIs - find indexes based on "f" variables
    [c f1_idx(i)] = min(abs(subject(1).f-f1(i)));
    [c f2_idx(i)] = min(abs(subject(1).f-f2(i)));
end

clear c f1 f2 i;
% ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,1};
        results.Ch1.ERO_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,1};
        results.Ch1.ERO_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,2};
        results.Ch2.ERO_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,2};
        results.Ch2.ERO_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% AVG_ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,1};
        results.Ch1.AVG_ERO_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,1};
        results.Ch1.AVG_ERO_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,2};
        results.Ch2.AVG_ERO_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,2};
        results.Ch2.AVG_ERO_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% PLI %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,1};
        results.Ch1.PLI_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,1};
        results.Ch1.PLI_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,2};
        results.Ch2.PLI_Theta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,2};
        results.Ch2.PLI_Theta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% Theta ROI MEAN calculations %
for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Theta_ERO_T(i,k) = mean2( results.Ch1.ERO_Theta_T{i,k}); %for target
        results.Ch1.mean_Theta_AVG_ERO_T(i,k) = mean2(results.Ch1.AVG_ERO_Theta_T{i,k}); 
        results.Ch1.mean_Theta_PLI_T(i,k) = mean2(results.Ch1.PLI_Theta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Theta_ERO_T(i,k) = mean2( results.Ch2.ERO_Theta_T{i,k}); %for target
        results.Ch2.mean_Theta_AVG_ERO_T(i,k) = mean2(results.Ch2.AVG_ERO_Theta_T{i,k}); 
        results.Ch2.mean_Theta_PLI_T(i,k) = mean2(results.Ch2.PLI_Theta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Theta_ERO_NT(i,k) = mean2( results.Ch1.ERO_Theta_NT{i,k}); %for target
        results.Ch1.mean_Theta_AVG_ERO_NT(i,k) = mean2(results.Ch1.AVG_ERO_Theta_NT{i,k}); 
        results.Ch1.mean_Theta_PLI_NT(i,k) = mean2(results.Ch1.PLI_Theta_NT{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Theta_ERO_NT(i,k) = mean2( results.Ch2.ERO_Theta_NT{i,k}); %for target
        results.Ch2.mean_Theta_AVG_ERO_NT(i,k) = mean2(results.Ch2.AVG_ERO_Theta_NT{i,k}); 
        results.Ch2.mean_Theta_PLI_NT(i,k) = mean2(results.Ch2.PLI_Theta_NT{i,k}); 
    end
end
clear i k;

% AVERAGE OVER ALL MICE Target and Non-target %
for i=1:n_roi;
results.Ch1.mean_Theta(i,1)=mean2(results.Ch1.mean_Theta_ERO_T(i,:)); %target ERO Ch1
results.Ch2.mean_Theta(i,1)=mean2(results.Ch2.mean_Theta_ERO_T(i,:)); %target ERO Ch2

results.Ch1.mean_Theta(i,2)=mean2(results.Ch1.mean_Theta_ERO_NT(i,:)); % not-target ERO Ch1
results.Ch2.mean_Theta(i,2)=mean2(results.Ch2.mean_Theta_ERO_NT(i,:)); % not-target ERO Ch2

results.Ch1.mean_Theta(i,3)=mean2(results.Ch1.mean_Theta_AVG_ERO_T(i,:)); %target ERP_ERO Ch1
results.Ch2.mean_Theta(i,3)=mean2(results.Ch2.mean_Theta_AVG_ERO_T(i,:)); %target ERP_ERO Ch2

results.Ch1.mean_Theta(i,4)=mean2(results.Ch1.mean_Theta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch1
results.Ch2.mean_Theta(i,4)=mean2(results.Ch2.mean_Theta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch2

results.Ch1.mean_Theta(i,5)=mean2(results.Ch1.mean_Theta_PLI_T(i,:)); %target PLI Ch1
results.Ch2.mean_Theta(i,5)=mean2(results.Ch2.mean_Theta_PLI_T(i,:)); %target PLI Ch2

results.Ch1.mean_Theta(i,6)=mean2(results.Ch1.mean_Theta_PLI_NT(i,:)); % not-target PLI Ch1
results.Ch2.mean_Theta(i,6)=mean2(results.Ch2.mean_Theta_PLI_NT(i,:)); % not-target PLI Ch2
end
clear i f1_idx f2_idx;

% BARS PLOT %
figure; 
subplot (2,3,1); bar(results.Ch1.mean_Theta(:,1));
title('Theta ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch1.mean_Theta(:,2));
title('Theta ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch1.mean_Theta(:,3));
title('Theta ERP-ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch1.mean_Theta(:,4));
title('Theta ERP-ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch1.mean_Theta(:,5));
title('Theta PLI > Ch1 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch1.mean_Theta(:,6));
title('Theta PLI > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

figure; 
subplot (2,3,1); bar(results.Ch2.mean_Theta(:,1));
title('Theta ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch2.mean_Theta(:,2));
title('Theta ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch2.mean_Theta(:,3));
title('Theta ERP-ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch2.mean_Theta(:,4));
title('Theta ERP-ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch2.mean_Theta(:,5));
title('Theta PLI > Ch2 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch2.mean_Theta(:,6));
title('Theta PLI > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

%% ALFA ROI definitions %%
        %BASAL%
f1(1) = 8;          %start of the freq. interval
f2(1) = 13;          %end of the freq. interval
        %P1/N1%
f1(2) = 8;          %start of the freq. interval
f2(2) = 13;          %end of the freq. interval
         %N2/P3%
f1(3) = 8;          %start of the freq. interval
f2(3) = 13;          %end of the freq. interval
        %P3%
f1(4) = 8;          %start of the freq. interval
f2(4) = 13;          %end of the freq. interval

for i = 1:n_roi         %over all ROIs - find indexes based on "f" variables
    [c f1_idx(i)] = min(abs(subject(1).f-f1(i)));
    [c f2_idx(i)] = min(abs(subject(1).f-f2(i)));
end

clear c f1 f2 i;

% ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,1};
        results.Ch1.ERO_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,1};
        results.Ch1.ERO_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,2};
        results.Ch2.ERO_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,2};
        results.Ch2.ERO_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% AVG_ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,1};
        results.Ch1.AVG_ERO_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,1};
        results.Ch1.AVG_ERO_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,2};
        results.Ch2.AVG_ERO_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,2};
        results.Ch2.AVG_ERO_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% PLI %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,1};
        results.Ch1.PLI_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,1};
        results.Ch1.PLI_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,2};
        results.Ch2.PLI_Alfa_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,2};
        results.Ch2.PLI_Alfa_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% Alfa ROI MEAN calculations %
for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Alfa_ERO_T(i,k) = mean2( results.Ch1.ERO_Alfa_T{i,k}); %for target
        results.Ch1.mean_Alfa_AVG_ERO_T(i,k) = mean2(results.Ch1.AVG_ERO_Alfa_T{i,k}); 
        results.Ch1.mean_Alfa_PLI_T(i,k) = mean2(results.Ch1.PLI_Alfa_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Alfa_ERO_T(i,k) = mean2( results.Ch2.ERO_Alfa_T{i,k}); %for target
        results.Ch2.mean_Alfa_AVG_ERO_T(i,k) = mean2(results.Ch2.AVG_ERO_Alfa_T{i,k}); 
        results.Ch2.mean_Alfa_PLI_T(i,k) = mean2(results.Ch2.PLI_Alfa_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Alfa_ERO_NT(i,k) = mean2( results.Ch1.ERO_Alfa_NT{i,k}); %for target
        results.Ch1.mean_Alfa_AVG_ERO_NT(i,k) = mean2(results.Ch1.AVG_ERO_Alfa_NT{i,k}); 
        results.Ch1.mean_Alfa_PLI_NT(i,k) = mean2(results.Ch1.PLI_Alfa_NT{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Alfa_ERO_NT(i,k) = mean2( results.Ch2.ERO_Alfa_NT{i,k}); %for target
        results.Ch2.mean_Alfa_AVG_ERO_NT(i,k) = mean2(results.Ch2.AVG_ERO_Alfa_NT{i,k}); 
        results.Ch2.mean_Alfa_PLI_NT(i,k) = mean2(results.Ch2.PLI_Alfa_NT{i,k}); 
    end
end
clear i k;

% AVERAGE OVER ALL MICE Target and Non-target %
for i=1:n_roi;
results.Ch1.mean_Alfa(i,1)=mean2(results.Ch1.mean_Alfa_ERO_T(i,:)); %target ERO Ch1
results.Ch2.mean_Alfa(i,1)=mean2(results.Ch2.mean_Alfa_ERO_T(i,:)); %target ERO Ch2

results.Ch1.mean_Alfa(i,2)=mean2(results.Ch1.mean_Alfa_ERO_NT(i,:)); % not-target ERO Ch1
results.Ch2.mean_Alfa(i,2)=mean2(results.Ch2.mean_Alfa_ERO_NT(i,:)); % not-target ERO Ch2

results.Ch1.mean_Alfa(i,3)=mean2(results.Ch1.mean_Alfa_AVG_ERO_T(i,:)); %target ERP_ERO Ch1
results.Ch2.mean_Alfa(i,3)=mean2(results.Ch2.mean_Alfa_AVG_ERO_T(i,:)); %target ERP_ERO Ch2

results.Ch1.mean_Alfa(i,4)=mean2(results.Ch1.mean_Alfa_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch1
results.Ch2.mean_Alfa(i,4)=mean2(results.Ch2.mean_Alfa_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch2

results.Ch1.mean_Alfa(i,5)=mean2(results.Ch1.mean_Alfa_PLI_T(i,:)); %target PLI Ch1
results.Ch2.mean_Alfa(i,5)=mean2(results.Ch2.mean_Alfa_PLI_T(i,:)); %target PLI Ch2

results.Ch1.mean_Alfa(i,6)=mean2(results.Ch1.mean_Alfa_PLI_NT(i,:)); % not-target PLI Ch1
results.Ch2.mean_Alfa(i,6)=mean2(results.Ch2.mean_Alfa_PLI_NT(i,:)); % not-target PLI Ch2
end
clear i f1_idx f2_idx;

% BARS PLOT %
figure; 
subplot (2,3,1); bar(results.Ch1.mean_Alfa(:,1));
title('Alfa ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch1.mean_Alfa(:,2));
title('Alfa ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch1.mean_Alfa(:,3));
title('Alfa ERP-ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch1.mean_Alfa(:,4));
title('Alfa ERP-ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch1.mean_Alfa(:,5));
title('Alfa PLI > Ch1 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch1.mean_Alfa(:,6));
title('Alfa PLI > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

figure; 
subplot (2,3,1); bar(results.Ch2.mean_Alfa(:,1));
title('Alfa ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch2.mean_Alfa(:,2));
title('Alfa ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch2.mean_Alfa(:,3));
title('Alfa ERP-ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch2.mean_Alfa(:,4));
title('Alfa ERP-ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch2.mean_Alfa(:,5));
title('Alfa PLI > Ch2 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch2.mean_Alfa(:,6));
title('Alfa PLI > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

%% BETA ROI definitions %%
        %BASAL%
f1(1) = 13;          %start of the freq. interval
f2(1) = 30;          %end of the freq. interval
        %P1/N1%
f1(2) = 13;          %start of the freq. interval
f2(2) = 30;          %end of the freq. interval
         %N2/P3%
f1(3) = 13;          %start of the freq. interval
f2(3) = 30;          %end of the freq. interval
        %P3%
f1(4) = 13;          %start of the freq. interval
f2(4) = 30;          %end of the freq. interval

for i = 1:n_roi         %over all ROIs - find indexes based on "f" variables    
    [c f1_idx(i)] = min(abs(subject(1).f-f1(i)));
    [c f2_idx(i)] = min(abs(subject(1).f-f2(i)));
end

clear c f1 f2 i;

% ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,1};
        results.Ch1.ERO_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,1};
        results.Ch1.ERO_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,2};
        results.Ch2.ERO_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,2};
        results.Ch2.ERO_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% AVG_ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,1};
        results.Ch1.AVG_ERO_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,1};
        results.Ch1.AVG_ERO_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,2};
        results.Ch2.AVG_ERO_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,2};
        results.Ch2.AVG_ERO_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% PLI %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,1};
        results.Ch1.PLI_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,1};
        results.Ch1.PLI_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,2};
        results.Ch2.PLI_Beta_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,2};
        results.Ch2.PLI_Beta_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% Beta ROI MEAN calculations %
for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Beta_ERO_T(i,k) = mean2( results.Ch1.ERO_Beta_T{i,k}); %for target
        results.Ch1.mean_Beta_AVG_ERO_T(i,k) = mean2(results.Ch1.AVG_ERO_Beta_T{i,k}); 
        results.Ch1.mean_Beta_PLI_T(i,k) = mean2(results.Ch1.PLI_Beta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Beta_ERO_T(i,k) = mean2( results.Ch2.ERO_Beta_T{i,k}); %for target
        results.Ch2.mean_Beta_AVG_ERO_T(i,k) = mean2(results.Ch2.AVG_ERO_Beta_T{i,k}); 
        results.Ch2.mean_Beta_PLI_T(i,k) = mean2(results.Ch2.PLI_Beta_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Beta_ERO_NT(i,k) = mean2( results.Ch1.ERO_Beta_NT{i,k}); %for target
        results.Ch1.mean_Beta_AVG_ERO_NT(i,k) = mean2(results.Ch1.AVG_ERO_Beta_NT{i,k}); 
        results.Ch1.mean_Beta_PLI_NT(i,k) = mean2(results.Ch1.PLI_Beta_NT{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Beta_ERO_NT(i,k) = mean2( results.Ch2.ERO_Beta_NT{i,k}); %for target
        results.Ch2.mean_Beta_AVG_ERO_NT(i,k) = mean2(results.Ch2.AVG_ERO_Beta_NT{i,k}); 
        results.Ch2.mean_Beta_PLI_NT(i,k) = mean2(results.Ch2.PLI_Beta_NT{i,k}); 
    end
end
clear i k;

% AVERAGE OVER ALL MICE Target and Non-target %
for i=1:n_roi;
results.Ch1.mean_Beta(i,1)=mean2(results.Ch1.mean_Beta_ERO_T(i,:)); %target ERO Ch1
results.Ch2.mean_Beta(i,1)=mean2(results.Ch2.mean_Beta_ERO_T(i,:)); %target ERO Ch2

results.Ch1.mean_Beta(i,2)=mean2(results.Ch1.mean_Beta_ERO_NT(i,:)); % not-target ERO Ch1
results.Ch2.mean_Beta(i,2)=mean2(results.Ch2.mean_Beta_ERO_NT(i,:)); % not-target ERO Ch2

results.Ch1.mean_Beta(i,3)=mean2(results.Ch1.mean_Beta_AVG_ERO_T(i,:)); %target ERP_ERO Ch1
results.Ch2.mean_Beta(i,3)=mean2(results.Ch2.mean_Beta_AVG_ERO_T(i,:)); %target ERP_ERO Ch2

results.Ch1.mean_Beta(i,4)=mean2(results.Ch1.mean_Beta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch1
results.Ch2.mean_Beta(i,4)=mean2(results.Ch2.mean_Beta_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch2

results.Ch1.mean_Beta(i,5)=mean2(results.Ch1.mean_Beta_PLI_T(i,:)); %target PLI Ch1
results.Ch2.mean_Beta(i,5)=mean2(results.Ch2.mean_Beta_PLI_T(i,:)); %target PLI Ch2

results.Ch1.mean_Beta(i,6)=mean2(results.Ch1.mean_Beta_PLI_NT(i,:)); % not-target PLI Ch1
results.Ch2.mean_Beta(i,6)=mean2(results.Ch2.mean_Beta_PLI_NT(i,:)); % not-target PLI Ch2
end
clear i f1_idx f2_idx;

% BARS PLOT %
figure; 
subplot (2,3,1); bar(results.Ch1.mean_Beta(:,1));
title('Beta ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch1.mean_Beta(:,2));
title('Beta ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch1.mean_Beta(:,3));
title('Beta ERP-ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch1.mean_Beta(:,4));
title('Beta ERP-ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch1.mean_Beta(:,5));
title('Beta PLI > Ch1 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch1.mean_Beta(:,6));
title('Beta PLI > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

figure; 
subplot (2,3,1); bar(results.Ch2.mean_Beta(:,1));
title('Beta ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch2.mean_Beta(:,2));
title('Beta ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch2.mean_Beta(:,3));
title('Beta ERP-ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch2.mean_Beta(:,4));
title('Beta ERP-ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch2.mean_Beta(:,5));
title('Beta PLI > Ch2 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch2.mean_Beta(:,6));
title('Beta PLI > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

%% GAMMA ROI definitions %%
        %BASAL%
f1(1) = 32;          %start of the freq. interval
f2(1) = 70;          %end of the freq. interval
        %P1/N1%
f1(2) = 32;          %start of the freq. interval
f2(2) = 70;          %end of the freq. interval
         %N2/P3%
f1(3) = 32;          %start of the freq. interval
f2(3) = 70;          %end of the freq. interval
        %P3%
f1(4) = 32;          %start of the freq. interval
f2(4) = 70;          %end of the freq. interval

for i = 1:n_roi         %over all ROIs - find indexes based on "f" variables  
    [c f1_idx(i)] = min(abs(subject(1).f-f1(i)));
    [c f2_idx(i)] = min(abs(subject(1).f-f2(i)));
end

clear c f1 f2 i;

% ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,1};
        results.Ch1.ERO_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,1};
        results.Ch1.ERO_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{1,2};
        results.Ch2.ERO_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).ERO{2,2};
        results.Ch2.ERO_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% AVG_ERO %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,1};
        results.Ch1.AVG_ERO_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,1};
        results.Ch1.AVG_ERO_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{1,2};
        results.Ch2.AVG_ERO_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).AVG_ERO{2,2};
        results.Ch2.AVG_ERO_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% PLI %
for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,1};
        results.Ch1.PLI_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,1};
        results.Ch1.PLI_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{1,2};
        results.Ch2.PLI_Gamma_NT{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

for i=1:n_roi;
    for k= 1:n_sbj;
        temp=subject(k).PLI{2,2};
        results.Ch2.PLI_Gamma_T{i,k} = temp(f1_idx(i):f2_idx(i),t1_idx(i):t2_idx(i));
    end
end
clear i k temp;

% Gamma ROI MEAN calculations %
for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Gamma_ERO_T(i,k) = mean2( results.Ch1.ERO_Gamma_T{i,k}); %for target
        results.Ch1.mean_Gamma_AVG_ERO_T(i,k) = mean2(results.Ch1.AVG_ERO_Gamma_T{i,k}); 
        results.Ch1.mean_Gamma_PLI_T(i,k) = mean2(results.Ch1.PLI_Gamma_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Gamma_ERO_T(i,k) = mean2( results.Ch2.ERO_Gamma_T{i,k}); %for target
        results.Ch2.mean_Gamma_AVG_ERO_T(i,k) = mean2(results.Ch2.AVG_ERO_Gamma_T{i,k}); 
        results.Ch2.mean_Gamma_PLI_T(i,k) = mean2(results.Ch2.PLI_Gamma_T{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch1.mean_Gamma_ERO_NT(i,k) = mean2( results.Ch1.ERO_Gamma_NT{i,k}); %for target
        results.Ch1.mean_Gamma_AVG_ERO_NT(i,k) = mean2(results.Ch1.AVG_ERO_Gamma_NT{i,k}); 
        results.Ch1.mean_Gamma_PLI_NT(i,k) = mean2(results.Ch1.PLI_Gamma_NT{i,k}); 
    end
end
clear i k;

for i=1:n_roi;
    for k= 1:n_sbj
        results.Ch2.mean_Gamma_ERO_NT(i,k) = mean2( results.Ch2.ERO_Gamma_NT{i,k}); %for target
        results.Ch2.mean_Gamma_AVG_ERO_NT(i,k) = mean2(results.Ch2.AVG_ERO_Gamma_NT{i,k}); 
        results.Ch2.mean_Gamma_PLI_NT(i,k) = mean2(results.Ch2.PLI_Gamma_NT{i,k}); 
    end
end
clear i k;

% AVERAGE OVER ALL MICE Target and Non-target %
for i=1:n_roi;
results.Ch1.mean_Gamma(i,1)=mean2(results.Ch1.mean_Gamma_ERO_T(i,:)); %target ERO Ch1
results.Ch2.mean_Gamma(i,1)=mean2(results.Ch2.mean_Gamma_ERO_T(i,:)); %target ERO Ch2

results.Ch1.mean_Gamma(i,2)=mean2(results.Ch1.mean_Gamma_ERO_NT(i,:)); % not-target ERO Ch1
results.Ch2.mean_Gamma(i,2)=mean2(results.Ch2.mean_Gamma_ERO_NT(i,:)); % not-target ERO Ch2

results.Ch1.mean_Gamma(i,3)=mean2(results.Ch1.mean_Gamma_AVG_ERO_T(i,:)); %target ERP_ERO Ch1
results.Ch2.mean_Gamma(i,3)=mean2(results.Ch2.mean_Gamma_AVG_ERO_T(i,:)); %target ERP_ERO Ch2

results.Ch1.mean_Gamma(i,4)=mean2(results.Ch1.mean_Gamma_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch1
results.Ch2.mean_Gamma(i,4)=mean2(results.Ch2.mean_Gamma_AVG_ERO_NT(i,:)); %not-target ERP_ERO Ch2

results.Ch1.mean_Gamma(i,5)=mean2(results.Ch1.mean_Gamma_PLI_T(i,:)); %target PLI Ch1
results.Ch2.mean_Gamma(i,5)=mean2(results.Ch2.mean_Gamma_PLI_T(i,:)); %target PLI Ch2

results.Ch1.mean_Gamma(i,6)=mean2(results.Ch1.mean_Gamma_PLI_NT(i,:)); % not-target PLI Ch1
results.Ch2.mean_Gamma(i,6)=mean2(results.Ch2.mean_Gamma_PLI_NT(i,:)); % not-target PLI Ch2
end
clear i n_roi f1_idx f2_idx t1_idx t2_idx

% BARS PLOT %
figure; 
subplot (2,3,1); bar(results.Ch1.mean_Gamma(:,1));
title('Gamma ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch1.mean_Gamma(:,2));
title('Gamma ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch1.mean_Gamma(:,3));
title('Gamma ERP-ERO > Ch1 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch1.mean_Gamma(:,4));
title('Gamma ERP-ERO > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch1.mean_Gamma(:,5));
title('Gamma PLI > Ch1 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch1.mean_Gamma(:,6));
title('Gamma PLI > Ch1 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

figure; 
subplot (2,3,1); bar(results.Ch2.mean_Gamma(:,1));
title('Gamma ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);
subplot (2,3,4); bar(results.Ch2.mean_Gamma(:,2));
title('Gamma ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 1.4]);

subplot (2,3,2); bar(results.Ch2.mean_Gamma(:,3));
title('Gamma ERP-ERO > Ch2 target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);
subplot (2,3,5); bar(results.Ch2.mean_Gamma(:,4));
title('Gamma ERP-ERO > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); % ylim([0 3.5]);

subplot (2,3,3); bar(results.Ch2.mean_Gamma(:,5));
title('Gamma PLI > Ch2 target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);
subplot (2,3,6); bar(results.Ch2.mean_Gamma(:,6));
title('Gamma PLI > Ch2 not-target');
xlabel('time point');
ylabel('ROI'); ylim([0 0.6]);

%% SAVE ALL DATA and CLEAR TEMPORAL FILES%%
% clear results.Ch1.ERO_T results.Ch1.ERO_NT results.Ch1.ERP_ERO_T results.Ch1.ERP_ERO_NT results.Ch1.PLI_T results.Ch1.PLI_NT
% clear results.Ch2.ERO_T results.Ch2.ERO_NT results.Ch2.ERP_ERO_T results.Ch2.ERP_ERO_NT results.Ch2.PLI_T results.Ch2.PLI_NT
%save ROI_NAME_results.mat results
%clear all