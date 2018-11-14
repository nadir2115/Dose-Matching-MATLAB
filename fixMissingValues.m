%function to impute missing values
function data= fixMissingValues(dv1)
dv1(dv1==-99)=1.0001;
dv1(isnan(dv1)) = 1.0001;
for i= 1:size(dv1,1)        
    for j= 1:size(dv1,2)
        if dv1(i,j) == 1.0001
            temp= dv1(:,j);
            temp(temp==1.001)= [];
            dv1(i,j)=median(temp); 
        end
    end
end
data= dv1;