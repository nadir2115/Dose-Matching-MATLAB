function array= table2numarray(datav1)
m=size(datav1,1); 
n =size(datav1,2);

 % Converting numerical values labeled as strings to double
for i = 2:12     %Matlab can only do one variable at a time??? Review
datav1.(i)= str2double(datav1.(i)); 
end

% Splitting tables to 4 arrays with specific data types
dv1a= table2array(datav1(:,1)); 
dv1b= table2array(datav1(:,2:12));
dv1c= table2array(datav1(:,13)); 
dv1d= table2array(datav1(:,14:n));

    for i=1:m
        %Convert PID to numers removing first and last 2 letters
        x=char(dv1a(i)); d1a(i,1)=str2num(x(2:7)); clear x; 
        %Convert severity to integer
        if dv1c(i)=="H"   
            d1c(i)= 1;
        elseif dv1c(i)=="L"
            d1c(i)= 0;
        end
    end
    
%Combining into 1 table
dv1d= str2double(dv1d);
% dv1d= str2double(dv1d);
d1c = d1c'; %Turning row vector to column vector
d =[d1a dv1b d1c dv1d]; 

%One hot encoding
% Site code
d(d(:,15)==1,n+1)=1;
d(d(:,15)==2,n+2)=1;
d(d(:,15)==3,n+3)=1;
%Ethnicity- hispanic/non-his/unidentified (Feature removed- only 4 patients)
d(d(:,17)==1,n+4)=1;
d(d(:,17)==2,n+5)=1;
d(d(:,17)==3,n+6)=1;
%Race- 1=American Indian, Aleutian, or Eskimo~2=Asian~3=Pacific Islander~4=African-American~5=Caucasian~6=Other~7=Declined to report
d(d(:,18)==1,n+7)=1;
d(d(:,18)==2,n+8)=1;
d(d(:,18)==3,n+9)=1;
d(d(:,18)==4,n+10)=1;
d(d(:,18)==5,n+11)=1;
d(d(:,18)==6,n+12)=1;
d(d(:,18)==7,n+13)=1;
%Stroke type 1=ISCHEMIC (w/o hemorrhagic conversion)~2=ISCHEMIC (w/ hemorrhagic conversion)~3=INTRAPARENCHYMAL HEMORRHAGIC~4=OTHER 
d(d(:,21)==1,n+14)=1;
d(d(:,21)==2,n+15)=1;
d(d(:,21)==3,n+16)=1;
d(d(:,21)==4,n+17)=1;
%Stroke location 1=R Hemisphere~2=L Hemisphere~3=Cerebellar (only 1 case so we ignore it)
% ~4=Brainstem~5=Other 
d(d(:,22)==1,n+18)=1;
d(d(:,22)==2,n+19)=1;
d(d(:,22)==4,n+20)=1;
d(d(:,22)==5,n+21)=1;

d(:,[15 17 18 21 22])=[];%Removing one hot encoded columns

array= d;


    