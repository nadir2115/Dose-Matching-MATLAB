 %function to normalize values to between 0 and 1
function [datav] = normalize(datas)        
for i=1:size(datas,2)           
       datav(:,i)= (datas(:,i)-mean(datas(:,i)))./std(datas(:,i));
end