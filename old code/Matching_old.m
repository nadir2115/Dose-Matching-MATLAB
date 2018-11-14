clc; clear; close all;
cd 'C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis'

% %Reading relevant datasets
% %xlsread
% edataset = xlsread('00042_evalDataset_notext.csv');
% ptchar = xlsread('PatientCharacteristics.csv');

% [header, edata]= importCSVFile('C:\Users\nadir\Desktop\Matching\BKN 599 Advanced Data Analysis','00042_evalDataset_notext.csv',1)

%readtable
evaldata = readtable('00042_evalDataset_notext.csv'); 
patchar = readtable('PatientCharacteristics.csv'); %Demographic Data
lesiond = readtable('Raw Data_Included Subjects Only_mod.csv'); %Lesion data-most fields missing
updated = readtable('icare results update.csv');

%Changing relevant column names to "pid"
patchar.Properties.VariableNames([6]) = {'pid'};
updated.Properties.VariableNames([1]) = {'pid'};
lesiond.Properties.VariableNames([1]) = {'pid'};

%Getting rid of bad features
lesiond(:,[2:7 8 15:27 28:142]) = []; % Removing unimportant/missing data
evaldata(:,[1 3 4 7:94 119 127 190 191 193 195 197 199 201 203 205 207 243 245 248:259 261 263:265 275 276 278 281 283 285 287 289 291 295:298 301 303 305 307 309 311 316 318 320 322 324 326 328 330 332 336:339 342 344 346 348 350 352 357 359 361:364 369 370 375]) = []; % Removing unimportant/missing data 
patchar(:,[1 2 4 5 10:12 22:25])= [];

% evaldata(evaldata.==-99)= "NA"    % Code doesn't work

ev1 = evaldata;
ev1(ev1.visitnum~=1,:) = [];
ev2 = evaldata;
ev2(ev2.visitnum~=2,:) = [];
% Visit 3 and 4 not used by Andrew and Alex

%Trimming visit-2 no-shows from visit-1 data
ev2_temp = ev2(:,2)
ev1 = join(ev2_temp,ev1,'Keys','pid');
clear ev2_temp

%Combining 2 tables by PID
datav1 = join(ev1,patchar,'Keys','pid');
datav1= join(datav1, lesiond,'Keys','pid');

datav2 = join(ev2,patchar,'Keys','pid'); %The key variable for B must contain all values in the key variable for A.

% Checking distribution of time of treatment since stroke
% hist(merg2.onset_to_rand, 15);

