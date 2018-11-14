clc; clear; close all;
cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'

%readtable
evaldata = readtable('00042_evalDataset_notext.csv'); 
patchar = readtable('PatientCharacteristics.csv'); %Demographic Data

cd 'C:\Users\nadir\Desktop\Matching\Matlab'
%Changing relevant column names to "pid"
patchar.Properties.VariableNames([6]) = {'pid'};
% lesiond.Properties.VariableNames([1]) = {'pid'};

%Getting rid of bad features
evaldata(:,[1 3 4 7:119 120:127 129:207 209:245 248:259 261:276 278:364 369:375]) = []; % Removing unimportant/missing data 
patchar(:,[1:5 8 10:12 22:25 27 30])= [];
% lesiond(:,[2:7 8 15:27 28:142]) = []; % Removing unimportant/missing data

%Combining 2 tables by PIDS
datav = join(evaldata,patchar,'Keys','pid');


datav1= extractVisitData(datav, 1);
datav2= extractVisitData(datav, 2);
% Visit 3 and 4 not used by Andrew and Alex

%Trimming visit-2 no-shows from visit-1 data
ev2_temp = datav2(:,2); datav1 = join(ev2_temp,datav1,'Keys','pid'); clear ev2_temp

%Rearranging columns so PID is dirst column
datav2= [datav2(:,2) datav2(:,1) datav2(:,3:28)];

dv1=table2numarrayold(datav1); %converting table to array of numbers
dv1= fixMissingValues(dv1); %fixing missing valuest
dv1= normalize_data(dv1);   %Using custom function to normalize data
dv1(:,1)= [];

dv2=table2numarrayold(datav2); %converting table to array of numbers
dv2= fixMissingValues(dv2); %fixing missing valuest
dv2= normalize_data(dv2);   %Using custom function to normalize data
dv2(:,1)= [];
%% EUCLIDIAN DISTANCE MATCHING-----------------------------------------------
dv1t =dv1(1:100, :);  %splitting into groups- let first 100 subjects be treated
% matched_index=[];
for i= 1:100
    k=10000;
    for j = 101:300   %matching with replacement
        if pdist2(dv1t(i,:),dv1(j,:), @nanhamdist)<k %nanhamdist accounts for NAN values
            matched_index(i,1)= j;
            k= pdist2(dv1t(i,:),dv1(j,:), @nanhamdist);
            matched_index(i,2)= k;    %for checking error value
        end
    end
end
[h, p, treated_effect, control_effect]= ttest_wmft (matched_index, dv1, dv2);
effect_histogram(treated_effect, control_effect, "Eucladian")

plot_std_means(matched_index, "Euclidian", dv1) %See standardized difference of means; 
clear h p matched_index i j k;             % clearing varaiables
%% PROPENSITY SCORE MATCHING-----------------------------------------------
propdata= [dv1 zeros(size((dv1),1),1)];
propdata(1:100, 28)=1; %Assigning first 100 elements to group 1 and rest to group 0
X = propdata(:, 1:27); y = propdata(:, 28); % Creating X,y
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

%MATCHING
for i= 1:100
    k=10000;
    for j = 101:300   %matching with replacement
        if pdist2(prop_score(i),prop_score(j))<k %nanhamdist accounts for NAN values
            matched_index(i,1)= j;
            k= pdist2(prop_score(i),prop_score(j));
            matched_index(i,2)= k;    %for checking error value
        end
    end
end

[h, p, treated_effect, control_effect]= ttest_wmft (matched_index, dv1, dv2);
effect_histogram(treated_effect, control_effect, "Propensity")

plot_std_means(matched_index, "Propensity score", dv1) %See standardized difference of means; 
clear i j k p h propdata theta X y m n options cost grad matched_index

%Checking calipers
for i =1:300
logit(i,1) = log(prop_score(i)/(1-prop_score(i)));
end
caliper= 0.2*std(logit)
%% MAHALANOBIS DISTANCE MATCHING
for i= 1:100
    k=1000000;
    for j = 101:300   %matching with replacement
        if sqrt((dv1(i,:)-dv1(j,:))*inv(cov(dv1))*(dv1(i,:)-dv1(j,:))')<k %nanhamdist accounts for NAN values
            matched_index(i,1)= j;
            k= sqrt((dv1(i,:)-dv1(j,:))*inv(cov(dv1))*(dv1(i,:)-dv1(j,:))');
            matched_index(i,2)= k;    %for checking error value
        end
    end
end
[h, p, treated_effect, control_effect]= ttest_wmft (matched_index, dv1, dv2);
effect_histogram(treated_effect, control_effect, "Mahalanobis")
plot_std_means(matched_index, "Mahalanobis", dv1) %See standardized difference of means; 

clear p h i j k;