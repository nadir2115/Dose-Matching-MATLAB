% Dose matching to evaluate impacts of high or low dose- Nadir Nibras
clc; clear; close all;

cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'
evalData = readtable('00042_evalDataset_notext.csv');               % unmatched + test data
patChar = readtable('PatientCharacteristics_withDose.csv');         % Demographic Data

% variables
thresh= 27;
replace= 0;
onehotencode= 0;

cd 'C:\Users\nadir\Desktop\Matching\Matlab'
patChar.Properties.VariableNames([6]) = {'pid'};                    % Change column name: "PID" to "pid"
evalData(:,[1:4 7:119 120:127 129:203 205:245 248:259 ...
    261:276 278:364 369:375]) = [];                                 % Remove sparse/useless features 
patChar(:,[1:5 7 8 10:11 13 23:26 28 31 33])= [];                   % Remove sparse/useless features
dataMerged = join(evalData,patChar,'Keys','pid');                   % Combine 2 tables by pid column
dataVisit1= extractVisitData(dataMerged, 1);                        % Extract visit 1 data
dataVisit2= extractVisitData(dataMerged, 2);                        % Extract visit 2 data

temp = dataVisit2(:,1);                                             % Trim visit-2 no-shows from visit-1
dataVisit1 = join(temp,dataVisit1,'Keys','pid'); 
clear temp

% Visit 1 data
visit1Mat= table2matrix(dataVisit1, onehotencode);                 	% convert table to numerical matrix
visit1Mat= fixMissingValues(visit1Mat);                             % impute missing values (mean)
visit1Mat(:,1)= [];                                                 % remove PID column
visit1MatNotNormalized= visit1Mat;                                  % non-normalized Dataset saved for later    
visit1Mat= [normalize(visit1Mat(:,1:end-1)),visit1Mat(:,end)];      % normalize all data except dose

% Visit 2 data
visit2Mat= table2matrix(dataVisit2, onehotencode);                  % convert table to numerical matrix
visit2Mat= fixMissingValues(visit2Mat);                             % impute missing values (mean)
visit2Mat(:,1)= [];                                                 % remove PID column
visit2MatNotNormalized= visit2Mat;                                  % dataset saved for plots    
visit2Mat= [normalize(visit2Mat(:,1:end-1)),visit2Mat(:,end)];      % normalize all data except dose                 

treatedIndex= find(visit1Mat(:,size(visit1Mat,2))<=thresh);         % Index for treated group
unmatchedControlIndex= find(visit1Mat(:,size(visit1Mat,2))>thresh); % Index for control group
vis1TreatedData= visit1Mat(treatedIndex,1:end-1);                   % Separate treated group data
vis1ControlData= visit1Mat(unmatchedControlIndex,1:end-1);          % Separate controls(to match to treated group)

%% EUCLIDIAN DISTANCE MATCHING-----------------------------------------------
controlDataTemp = vis1ControlData;
for i= 1:size(vis1TreatedData,1)
    dist=1000; 
    for j = 1:size(controlDataTemp,1)                              
       if pdist2(vis1TreatedData(i,:),controlDataTemp(j,:))<dist            
           dist= pdist2(vis1TreatedData(i,:),controlDataTemp(j,:));	% new lowest distance 
           matchedSubIndex_e(i,:)= [j dist];                       	% Matched control sub-index + match distance                           
       end
    end
    if replace==0                                                   % IF without replacement
       controlDataTemp(matchedSubIndex_e(i,1),:)= 50*...
           ones(1,size(controlDataTemp,2));                         % Make used control unmatchable
    end
end
matchedControlIndex_e = unmatchedControlIndex(matchedSubIndex_e(:,1));
assessBalance(treatedIndex, unmatchedControlIndex, matchedControlIndex_e,...
    "Euclidian", visit1Mat)                                         % Visually evaluate balance 
% Test effect
[h, p, effectTreated, effectControl]= ttestImprovement(treatedIndex,...
    matchedControlIndex_e, visit1Mat, visit2Mat);
medianTreatedEffect_e= median(effectTreated);                       % Median improvement in treated group
medianControlEffect_e= median(effectControl);                       % Median improvement in control group 
effectsHistogram(effectTreated, effectControl, "Eucladian")
h
p

clear h p i j dist controlDataTemp;                                                 
V1data_Treated= visit1MatNotNormalized(treatedIndex,:);
V1data_UnmatchedControl= visit1MatNotNormalized(unmatchedControlIndex,:);
V1data_MatchedControl_e= visit1MatNotNormalized(matchedControlIndex_e,:);

%% PROPENSITY SCORE MATCHING-----------------------------------------------
propData= visit1Mat;
propData(treatedIndex,end)=1;                                       % Set treated dose-value to 1
propData(unmatchedControlIndex,end)=0;                             	% Set control dose-value to 0
X = propData(:,1:end-1); y = propData(:,end);                       % Create log regression dataset 

%Logistic regression
[m, n] = size(X);
X = [ones(m, 1) X];                                                 % Add intercept term to X
initial_theta = ones(n+1, 1);                                       % Initialize fitting parameters
[cost, grad] = costFunction(initial_theta, X, y);                   % find initial cost & gradient
options = optimset('GradObj', 'on', 'MaxIter', 400);                % Set options for fminunc
[theta, cost] = ...
	fminunc(@(t)(costFunction(t, X, y)), initial_theta, options);   % Run fminunc: obtain optimal theta & matching cost 
propScore= sigmoid(X*theta);                                        % Calculate Prop. score from optimal Theta

caliper= 0.2*std(logit(propScore));                                 % Calculate calipers 
matchedSubIndex_incaliper=[];

tempPropScore=propScore;
for i= 1:size(treatedIndex)                                         % MATCHING
    dist=1000;
    logiti= logit(tempPropScore(treatedIndex(i)));
    for j = 1:size(unmatchedControlIndex)   
        logitj= logit(tempPropScore(unmatchedControlIndex(j)));
        if pdist2(logiti,logitj)<dist 
            dist= pdist2(logiti,logitj);
            matchedSubIndex_p(i,:)= [j dist];
        end
    end
    % Find out what values are within calipers
    if matchedSubIndex_p(i,2)<caliper 
        matchedSubIndex_incaliper= [matchedSubIndex_incaliper; i];	% store index of matches w/in calipers
    end  
    if replace==0                                                   % If without replacement
        tempPropScore(unmatchedControlIndex(...
            matchedSubIndex_p(i,1)))= 0;                            % make used prop score unmatchable
    end
end
% Redifine new treatment and control groups applying calipers
matchedSubIndex_p= matchedSubIndex_p(matchedSubIndex_incaliper, :); % subset control sub-index with calipers
treatedIndex_incaliper_p= treatedIndex(matchedSubIndex_incaliper,:);% subset treated index with calipers
matchedControlIndex_p = unmatchedControlIndex(...
    matchedSubIndex_p(:,1),:);                                      % main dataset index for matched control
assessBalance(treatedIndex_incaliper_p, unmatchedControlIndex, ...
    matchedControlIndex_p, "Propensity score", visit1Mat)           % Visually evaluate balance

% Test effect
[h, p, effectTreated, effectControl]= ttestImprovement (...
    treatedIndex_incaliper_p, matchedControlIndex_p, visit1Mat, visit2Mat);
medianTreatedEffect_p= median(effectTreated);                       % Median improvement in treated group
medianControlEffect_p= median(effectControl);                       % Median improvement in control group 
effectsHistogram(effectTreated, effectControl, "Propensity")
h
p
clear i j k p h propData theta initial_theta X y m n options cost grad logiti logitj dist tempPropScore
V1data_MatchedControl_p= visit1MatNotNormalized(matchedControlIndex_p,:);
V1data_Treated_p=  visit1MatNotNormalized(treatedIndex_incaliper_p,:);

%% MAHALANOBIS DISTANCE MATCHING
if onehotencode==0
    controlDataTemp = vis1ControlData;                            	
    for i= 1:size(treatedIndex)
        dist=100000;
        for j = 1:size(controlDataTemp)   
            if sqrt((vis1TreatedData(i,:)-controlDataTemp(j,:))*inv(cov(visit1Mat(:,1:end-1)))*...
                    (vis1TreatedData(i,:)-controlDataTemp(j,:))')<dist 
                dist= sqrt((vis1TreatedData(i,:)-controlDataTemp(j,:))*...
                    inv(cov(visit1Mat(:,1:end-1)))*(vis1TreatedData(i,:)-controlDataTemp(j,:))');
                matchedSubIndex_m(i,:)= [j,dist];   
            end
        end
       if replace==0                                            	% If without replacement
           controlDataTemp(matchedSubIndex_m(i,1),:)= 50*...
           ones(1,size(controlDataTemp,2));                       	% Make used control unmatchable
       end
        
    end
    matchedControlIndex_m = unmatchedControlIndex(matchedSubIndex_m(:,1),:);
    assessBalance(treatedIndex, unmatchedControlIndex,...
        matchedControlIndex_m,"Mahalanobis", visit1Mat)           	% Visually evaluate balance 
    % Test effect
    [h, p, effectTreated, effectControl]= ttestImprovement (...
        treatedIndex, matchedControlIndex_m, visit1Mat, visit2Mat);
    effectsHistogram(effectTreated, effectControl, "Mahalanobis")
    h
    p
    
    clear h p matchedSubIndex_m i j dist controlDataTemp;           	 
    V1data_MatchedControl_m= visit1MatNotNormalized(matchedControlIndex_m,:);
    unique_mahalanobis_controls= length(unique(matchedControlIndex_m)) 
end
%%
unique_euclidean_controls= length(unique(matchedControlIndex_e))
unique_propensity_controls= length(unique(matchedControlIndex_p))

% Mean of variates, and ttest results before and after matching
MatchEffect(:,1)= mean(V1data_Treated);
MatchEffect(:,2)= mean(V1data_UnmatchedControl);
[MatchEffect(:,3),MatchEffect(:,4)] = ...
    ttest2(V1data_Treated,V1data_UnmatchedControl);           
MatchEffect(:,5)= mean(V1data_MatchedControl_e);
[MatchEffect(:,6),MatchEffect(:,7)] = ...
    ttest2(V1data_Treated,V1data_MatchedControl_e); 
MatchEffect(:,8)= mean(V1data_Treated_p);
MatchEffect(:,9)= mean(V1data_MatchedControl_p);
[MatchEffect(:,10),MatchEffect(:,11)] = ...
    ttest2(V1data_Treated_p,V1data_MatchedControl_p); 
if onehotencode==0
    MatchEffect(:,12)= mean(V1data_Treated);
    MatchEffect(:,13)= mean(V1data_MatchedControl_m);
    [MatchEffect(:,14),MatchEffect(:,15)] = ...
        ttest2(V1data_Treated,V1data_MatchedControl_m); 
end

multihist(V1data_Treated)
multihist(V1data_MatchedControl_e)
multihist(V1data_MatchedControl_p)