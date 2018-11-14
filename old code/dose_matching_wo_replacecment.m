clc; clear; close all;
cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'

%variables
cutoff= 27;
replace= 1;

%readtable
evaldata = readtable('00042_evalDataset_notext.csv'); 
patchar = readtable('PatientCharacteristics_withDose.csv'); %Demographic Data

cd 'C:\Users\nadir\Desktop\Matching\Matlab'
%Changing relevant column names to "pid"
patchar.Properties.VariableNames([6]) = {'pid'};
% lesiond.Properties.VariableNames([1]) = {'pid'};

%Getting rid of bad features
evaldata(:,[1:4 7:119 120:127 129:203 205:245 248:259 261:276 278:364 369:375]) = []; % Removing unimportant/missing data 
patchar(:,[1:5 7 8 10:11 13 23:26 28 31 33])= [];
% lesiond(:,[2:7 8 15:142]) = []; % Removing unimportant/missing data

%Combining 2 tables by PID
datav = join(evaldata,patchar,'Keys','pid');

datav1= extractVisitData(datav, 1);
datav2= extractVisitData(datav, 2);
% Visit 3 and 4 not used by Andrew and Alex

%Trimming visit-2 no-shows from visit-1 data
ev2_temp = datav2(:,1); 
datav1 = join(ev2_temp,datav1,'Keys','pid'); 
clear ev2_temp

writetable(datav2,"Data for visit 2.csv");

dv1=tabnum_reducedfeatures(datav1); %converting table to array of numbers
dv1= fixMissingValues(dv1); %fixing missing valuest
dv1= [dv1(:,1:13) dv1(:,15:end) dv1(:,14)]; %moving dose hours to end
dv1(:,1)= []; %Removing PID column
dv1_wo_normalize= dv1;
dv1= [normalize_data(dv1(:,1:size(dv1,2)-1)) dv1(:,size(dv1,2))];  %Using custom function to normalize data

dv2=tabnum_reducedfeatures(datav2); %converting table to array of numbers
dv2= fixMissingValues(dv2); %fixing missing valuest
dv2= [dv2(:,1:13) dv2(:,15:end) dv2(:,14)]; %moving dose hours to end
dv2(:,1)= []; %Removing PID column
dv2_wo_normalize= dv2;
dv2= [normalize_data(dv2(:,1:size(dv2,2)-1)) dv2(:,size(dv2,2))] ;  %Using custom function to normalize data

% dv1(:,28)=50*rand(300,1); %%Pseudo dose added
treat= find(dv1(:,size(dv1,2))>cutoff);
control= find(dv1(:,size(dv1,2))<=cutoff);

dv1t= dv1(treat,1:end-1);
dv1c= dv1(control,1:end-1);
%% EUCLIDIAN DISTANCE MATCHING-----------------------------------------------
for i= 1:size(dv1t,1)
    k=10000;
    for j = 1:size(dv1c,1)                              
        if pdist2(dv1t(i,:),dv1c(j,:))<k            
            matched_index(i,1)= j;
            k= pdist2(dv1t(i,:),dv1c(j,:));
            matched_index(i,2)= k;    %for checking error value
        end
    end
    dv1c( matched_index(i,1),:)= 5*ones(1,size(dv1c,2)); %Making picked sample large so doesn't get picked again
end
matched_control_e = control(matched_index(:,1));
%Evaluating matches
plot_std_means(treat, control, matched_control_e, "Euclidian", dv1) %See standardized difference of means; 
%Testing effect
[h, p, treated_effect, control_effect]= ttest_wmft (treat, matched_control_e, dv1, dv2);
median_treated_effect_e= median(treated_effect); %Finding the normalized improvement in treated group
median_control_effect_e= median(control_effect); %Finding normalized improvement in control group 
effect_histogram(treated_effect, control_effect, "Eucladian")
h
p
clear h p matched_index i j k;             % clearing varaiables

d1treat= dv1_wo_normalize(treat,:);
d1control= dv1_wo_normalize(control,:);
d1m_control_e= dv1_wo_normalize(matched_control_e,:);
%% PROPENSITY SCORE MATCHING-----------------------------------------------
propdata= dv1;
%Assigning treatment 1 and control 0
propdata(treat, size(propdata,2))=1; 
propdata(control,size(propdata,2))=0;
X = propdata(:, 1:size(propdata,2)-1); y = propdata(:, size(propdata,2)); % Creating X,y
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
logit= log10(prop_score./(1-prop_score));
caliper= 0.2*std(logit) %CHANGED

matched_index_p=[]; 
incaliper=[];
prop_score_original= prop_score
%MATCHING
for i= 1:size(treat)
    k=1000;
    for j = 1:size(control)   
        logiti= log10(prop_score(treat(i))/(1-prop_score(treat(i))));
        logitj= log10(prop_score(control(j))/(1-prop_score(control(j))));
        
        if pdist2(logiti,logitj)<k 
            matched_index_p(i,1)= j;
            k= pdist2(logiti,logitj);
            matched_index_p(i,2)= k;   %for checking error value
        end
    end
    %Finding out what values are within calipers
    if matched_index_p(i,2)<caliper 
        incaliper= [incaliper; i];
    end  
    prop_score(control(matched_index_p(i,1)))= 0; %setting prop score as big so not reused
end
%Redifining new treatment and control groups with calipers
matched_index_p=matched_index_p(incaliper, :);
treat_new= treat(incaliper, :);
matched_control_p = control(matched_index_p(:,1),:);
%Evaluating matches
plot_std_means(treat_new, control, matched_control_p, "Propensity score", dv1) %See standardized difference of means; 
%Testing effect
[h, p, treated_effect, control_effect]= ttest_wmft (treat_new, matched_control_p, dv1, dv2);
median_treated_effect_p= median(treated_effect); %Finding the normalized improvement in treated group
median_control_effect_p= median(control_effect); %Finding normalized improvement in control group 
effect_histogram(treated_effect, control_effect, "Propensity")
h
p
clear i j k p h propdata theta X y m n options cost grad logiti logitj
d1m_control_p= dv1_wo_normalize(matched_control_p,:);
%% 
unique_euclidean_controls= length(unique(matched_control_e))
unique_propensity_controls= length(unique(matched_control_p))

%Mean of variates
X1(:,1)= mean(d1treat);
X1(:,2)= mean(d1control);
X1(:,5)= mean(d1m_control_e);
X1(:,8)= mean(d1m_control_p);
for i=1:size(X1,1)
[h,p] = ttest2(d1treat(:,i),d1control(:,i));
X1(i,3)= h;
X1(i,4)= p;
[h,p] = ttest2(d1treat(:,i),d1m_control_e(:,i));
X1(i,6)= h;
X1(i,7)= p;
[h,p] = ttest2(d1treat(:,i),d1m_control_p(:,i));
X1(i,9)= h;
X1(i,10)= p;
end

multihist(d1m_control_e)
multihist(d1m_control_p)