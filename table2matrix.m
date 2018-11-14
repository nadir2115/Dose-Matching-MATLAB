function array= table2matrix(datav1, onehotencode)
[m, n]=size(datav1);

 % Converting numerical values labeled as strings to double
for i = 2:12     %Matlab can only do one variable at a time??? Review
datav1.(i)= str2double(datav1.(i)); 
end

% Split tables to 4 arrays with specific data types
dv1a= table2array(datav1(:,1)); 
dv1b= table2array(datav1(:,2:12));
dv1c= table2array(datav1(:,13)); 
dv1d= table2array(datav1(:,14:n));

%Converting columns into double types
x=char(dv1a);       
d1a= str2num(x(:,2:7));                 %Convert PID to numbers removing first and last 2 letters
d1c(dv1c=="H")=1;   d1c=d1c';         	%Convert severity to integer, H=1, L=0
dv1d= str2double(dv1d);

d =[d1a dv1b d1c dv1d];                 %Combining into 1 table

if onehotencode==1 %One hot encoding
    % Site code
    d(d(:,15)==1,n+1)=1;
    d(d(:,15)==2,n+2)=1;
    d(d(:,15)==3,n+3)=1;

    % %Ethnicity- hispanic/non-his/unidentified (Feature 3(other) removed- 4 patients)
    d(d(:,17)==1,n+4)=1;

    %Race- 1=American Indian, Aleutian, or Eskimo~2=Asian~3=Pacific Islander~4=African-American~5=Caucasian~6=Other~7=Declined to report
    d(d(:,18)==1,n+5)=1; %Grouping 1,2,3 together as they make less than 10% of group total
    d(d(:,18)==2,n+5)=1;
    d(d(:,18)==3,n+5)=1;
    d(d(:,18)==4,n+6)=1;
    d(d(:,18)==5,n+7)=1;
    d(d(:,18)==6,n+8)=1; %Removing declined to report

    %Stroke type 1=ISCHEMIC (w/o hemorrhagic conversion)~2=ISCHEMIC (w/ hemorrhagic conversion)~3=INTRAPARENCHYMAL HEMORRHAGIC~4=OTHER 
    % Feature 1 and 4 made 5% of data
    % BINARY- 1=ISCHEMIC(w/o hemorrhagic conversion), 0=hemorrhagic/other(2-4)
    d(d(:,21)==1,n+9)=1;

    %Stroke location 1=R Hemisphere~2=L Hemisphere~3=Cerebellar 4=Brainstem~5=Other
    % (only 1 case of 3 so we ignore it)
    d(d(:,22)==1,n+10)=1;
    d(d(:,22)==2,n+11)=1;
    d(d(:,22)==4,n+12)=1;
    d(d(:,22)==5,n+13)=1;

    d(:,[15 17 18 21 22])=[];           %Removing redundant one-hot-encoded columns
    d= [d(:,1:13) d(:,15:end) d(:,14)]; %move dose hours to last column
elseif onehotencode==0
    d(:,[5 13 15:23 25 27])=[];         %Removing categorical columns
    d= [d(:,1:11) d(:,13:end) d(:,12)]; %move dose hours to last column
end
array= d;


    