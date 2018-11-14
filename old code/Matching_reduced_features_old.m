clc; clear; close all;
cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'

% %Reading relevant datasets using xlsread
% edataset = xlsread('00042_evalDataset_notext.csv');
% ptchar = xlsread('PatientCharacteristics.csv');

% [header, edata]= importCSVFile('C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis','00042_evalDataset_notext.csv',1)

%readtable
evaldata = readtable('00042_evalDataset_notext.csv'); 
patchar = readtable('PatientCharacteristics.csv'); %Demographic Data
% lesiond = readtable('Raw Data_Included Subjects Only_mod.csv'); %Lesion data-most fields missing
% updated = readtable('icare results update.csv');      % not important

%Changing relevant column names to "pid"
patchar.Properties.VariableNames([6]) = {'pid'};
% updated.Properties.VariableNames([1]) = {'pid'};
% lesiond.Properties.VariableNames([1]) = {'pid'};

%Getting rid of bad features
evaldata(:,[1 3 4 7:119 120:127 129:207 209:245 248:259 261:276 278:364 369:375]) = []; % Removing unimportant/missing data 
patchar(:,[1:5 8 10:12 22:25 27 30])= [];
% lesiond(:,[2:7 8 15:27 28:142]) = []; % Removing unimportant/missing data

% evaldata(evaldata.==-99)= "NA"    % Code doesn't work
ev1 = evaldata;
ev1(ev1.visitnum~=1,:) = [];
ev1.visitnum= []; %removing redundant column
ev2 = evaldata;
ev2(ev2.visitnum~=2,:) = [];
ev2.visitnum= []; %removing redundant column
% Visit 3 and 4 not used by Andrew and Alex

%Trimming visit-2 no-shows from visit-1 data
ev2_temp = ev2(:,2);
ev1 = join(ev2_temp,ev1,'Keys','pid');
clear ev2_temp

%Rearranging columns so ev1 and ev2 match
x=ev2.pid; ev2(:,2)=[]; ev2= [x ev2]; clear x; ev2.Properties.VariableNames{'Var1'} = 'pid';

%Combining 2 tables by PID
datav1 = join(ev1,patchar,'Keys','pid');
% datav1= Rjoin(datav1, lesiond,'Keys','pid'); %Doesn't work
datav2 = join(ev2,patchar,'Keys','pid'); %The key variable for B must contain all values in the key variable for A.

clear ev1; clear ev2; %data already in datav1/2

% Converting numerical values labeled as strings to double
for i =3:13     %Matlab can only do one variable at a time??? Review
datav1.(i)= str2double(datav1.(i));
datav2.(i)= str2double(datav2.(i));
end

%ORGANZING VISIT 1 DATA
% Splitting tables to 4 arrays with specific data types
dv1a= table2array(datav1(:,1));
dv1b= table2array(datav1(:,2:13));
dv1c= table2array(datav1(:,14:15));
dv1d= table2array(datav1(:,16:28));

for i=1:300
    x=char(dv1a(i)); d1a(i,1)=str2num(x(2:7)); clear x;
    if dv1c(i,1)=="D"   %Convert treat-group to int
        d1c(i,1)= 1;
    elseif dv1c(i,1)=="U"
        d1c(i,1)= 2;
    elseif dv1c(i,1)=="A"
        d1c(i,1)= 3;
    end
    if dv1c(i,2)=="H"   %Convert severity to in
        d1c(i,2)= 1;
    elseif dv1c(i,2)=="L"
        d1c(i,2)= 0;
    end
end

clear dv1a; 
clear dv1c; 
dv1 =[d1a dv1b d1c dv1d];
clear d1a; clear dv1b; clear d1c; clear dv1d;
dv1(dv1 == -99) = NaN; %getting rid of null values

for i=2:28
   dv1(:,i)= dv1(:,i)./max(dv1(:,i));
end
clear i %Clearing counter
% 
% writetable(datav1, "Visit1data_feature reduced.xls")
% xlswrite("Visit1data_processed.xls",dv1)




% Checking distribution of time of treatment since stroke
hist(datav1.wmft_mean_time_MA_PA);

