% clc; clear; close all;
clear
cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'

%readtable
evaldata = readtable('00042_evalDataset_notext.csv'); 
patchar = readtable('PatientCharacteristics_withDose.csv'); %Demographic Data

cd 'C:\Users\nadir\Desktop\Matching\Matlab'
%Changing relevant column names to "pid"
patchar.Properties.VariableNames([6]) = {'pid'};
% lesiond.Properties.VariableNames([1]) = {'pid'};

%Getting rid of bad features
evaldata(:,[1 3 4 7:119 120:127 129:207 209:245 248:259 261:276 278:364 369:375]) = []; % Removing unimportant/missing data 
patchar(:,[1:5 8 10:11 13 23:26 28 31 33])= [];
% lesiond(:,[2:7 8 15:27 28:142]) = []; % Removing unimportant/missing data

%Combining 2 tables by PID
datav = join(evaldata,patchar,'Keys','pid');
% datav1= Rjoin(datav1, lesiond,'Keys','pid'); %Doesn't work

datav1= extractVisitData(datav, 1);
datav2= extractVisitData(datav, 2);
% Visit 3 and 4 not used by Andrew and Alex

%Trimming visit-2 no-shows from visit-1 data
ev2_temp = datav2(:,2); datav1 = join(ev2_temp,datav1,'Keys','pid'); clear ev2_temp

%Rearranging columns so PID is dirst column
datav2= [datav2(:,2) datav2(:,1) datav2(:,3:29)];

dv1=table2numarray(datav1); %converting table to array of numbers
dv1= fixMissingValues(dv1); %fixing missing valuest
dv1= [dv1(:,1:15) dv1(:,17:29) dv1(:,16)];
% dv1(:,2)=[];

dv1= [normalize_data(dv1(:,1:28)) dv1(:,29)];  %Using custom function to normalize data
dv1(:,1)= [];

dv2=table2numarray(datav2); %converting table to array of numbers
dv2= fixMissingValues(dv2); %fixing missing valuest
dv2= [dv2(:,1:15) dv2(:,17:29) dv2(:,16)];
dv2= [normalize_data(dv2(:,1:28)) dv2(:,29)] ;  %Using custom function to normalize data
dv2(:,1)= [];

mean_std_means_e= []; mean_r_var_matched_e= [], mean_r_var_unmatched_e= [];
mean_std_means_p= []; mean_r_var_matched_p= [], mean_r_var_unmatched_p= [];
mean_std_means_m= []; mean_r_var_matched_m= [], mean_r_var_unmatched_m= [];

for kp= 1:100
dv1(:,28)=50*rand(300,1); %%Pseudo dose added
treat= find(dv1(:,28)<25);
control= find(dv1(:,28)>=25);

dv1t= dv1(treat,1:end-1);
dv1c= dv1(control,1:end-1);
%% EUCLIDIAN DISTANCE MATCHING-----------------------------------------------
for i= 1:size(dv1t,1)
    k=10000;
    for j = 1:size(dv1c,1)   %matching with replacement
        if pdist2(dv1t(i,:),dv1c(j,:), @nanhamdist)<k %nanhamdist accounts for NAN values
            matched_index(i,1)= j;
            k= pdist2(dv1t(i,:),dv1c(j,:), @nanhamdist);
            matched_index(i,2)= k;    %for checking error value
        end
    end
end
matched_control_e = control(matched_index(:,1),:);
%Evaluating matches
[mean_std_means_e(kp,:), mean_r_var_matched_e(kp), mean_r_var_unmatched_e(kp)]= ...
    find_std_means(treat, control, matched_control_e, dv1)
clear h p matched_index i j k;             % clearing varaiables

%% PROPENSITY SCORE MATCHING-----------------------------------------------
propdata= dv1;
%Assigning treatment 1 and control 0
propdata(treat, size(propdata,2))=1; 
propdata(control,size(propdata,2))=0;
X = propdata(:, 1:27); y = propdata(:, 28); % Creating X,y
%Performing logistic regression
[m, n] = size(X);
X = [ones(m, 1) X];     % Add intercept term to X
% Initialize fitting parameters and find initial cost and gradient
initial_theta = ones(n+1, 1);
[cost, grad] = costFunction(initial_theta, X, y);

%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);
%  Run fminunc to obtain the optimal theta
%  This function will return best theta and matching cost 
[theta, cost] = ...
	fminunc(@(t)(costFunction(t, X, y)), initial_theta, options);

%Calculating Propensity score from calculated Theta
prop_score= sigmoid(X*theta);

% Creating calipers
%Checking calipers
for i =1:300
logit(i) = log10(prop_score(i)/(1-prop_score(i)));
end
caliper= 0.2*std(logit);

matched_index_p=[]; incaliper=[];
%MATCHING
for i= 1:size(treat)
    k=1000;
    for j = 1:size(control)   %matching with replacement
        logiti= log10(prop_score(treat(i))/(1-prop_score(treat(i))));
        logitj= log10(prop_score(control(j))/(1-prop_score(control(j))));
        
        if pdist2(logiti,logitj)<k %nanhamdist accounts for NAN values
            matched_index_p(i,1)= j;
            k= pdist2(logiti,logitj);
            matched_index_p(i,2)= k;   %for checking error value
        end
    end
    %Finding out what values are within calipers
    if matched_index_p(i,2)<caliper 
        incaliper= [incaliper; i];
    end  
end
%Redifining new treatment and control groups with calipers
matched_index_p=matched_index_p(incaliper, :);
treat_new= treat(incaliper, :);
matched_control_p = control(matched_index_p(:,1),:);
%Evaluating matches
[mean_std_means_p(kp,:), mean_r_var_matched_p(kp), mean_r_var_unmatched_p(kp)]= ...
    find_std_means(treat, control, matched_control_p, dv1)
clear i j k p h propdata theta X y m n options cost grad logiti logitj

%% MAHALANOBIS DISTANCE MATCHING
matched_index_m=[];
for i= 1:size(treat)
    k=1000000;
    for j = 1:size(control)   %matching with replacement
        if sqrt((dv1t(i,:)-dv1c(j,:))*inv(cov(dv1(:,1:end-1)))*(dv1t(i,:)-dv1c(j,:))')<k %nanhamdist accounts for NAN values
            matched_index_m(i,1)= j;
            k= sqrt((dv1t(i,:)-dv1c(j,:))*inv(cov(dv1(:,1:end-1)))*(dv1t(i,:)-dv1c(j,:))');
            matched_index_m(i,2)= k;    %for checking error value
        end
    end
end
matched_control_m = control(matched_index_m(:,1),:);
%Evaluating matches
[mean_std_means_m(kp,:), mean_r_var_matched_m(kp), mean_r_var_unmatched_m(kp)]= ...
    find_std_means(treat, control, matched_control_m, dv1)
clear h p matched_index_m i j k;             % clearing varaiables
% unique_euclidean_controls= length(unique(matched_control_e))
% unique_mahalanobis_controls= length(unique(matched_control_m))
% unique_propensity_controls= length(unique(matched_control_p))

vtrt= var(prop_score(treat));
vcnt= var(prop_score(control));
vcnte= var(prop_score(matched_control_e));
vcntp= var(prop_score(matched_control_p));
vcntm= var(prop_score(matched_control_m));

var_prop_unmatched(kp)= max([vtrt/vcnt vcnt/vtrt]);
var_prop_matched_e(kp)= max([vtrt/vcnte vcnte/vtrt]);
var_prop_matched_p(kp)= max([vtrt/vcntp vcntp/vtrt]);
var_prop_matched_m(kp)= max([vtrt/vcntm vcntm/vtrt]);
end

mean_var_prop_unmatched= mean(var_prop_unmatched)
mean_var_prop_matched_euc= mean(var_prop_matched_e)
mean_var_prop_matched_mah= mean(var_prop_matched_m)
mean_var_prop_matched_prop= mean(var_prop_matched_p)

Std_means_euc_1000= mean(mean_std_means_e)
Std_means_mah_1000= mean(mean_std_means_m)
Std_means_prop_1000= mean(mean_std_means_p)
mean_ratio_variance_euc=mean(mean_r_var_matched_e)
mean_ratio_variance_mah=mean(mean_r_var_matched_m)
mean_ratio_variance_prop=mean(mean_r_var_matched_p)
mean_ratio_variance_unmatched=mean(mean_r_var_unmatched_e) %It's same for unmatched for all 3
