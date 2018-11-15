function assessBalance(treatedIndex, unmatchedControlIndex, matchedControlIndex, type, visit1Mat)

visit1Matt= visit1Mat(treatedIndex,1:end-1);
mean_t= mean(visit1Matt); 
std_treatedIndexed= std(visit1Matt);
unmatched_c= visit1Mat(unmatchedControlIndex,1:end-1);
mean_c_unmatched= mean(unmatched_c);

match_c= visit1Mat(matchedControlIndex(:,1),1:end-1);
mean_c_matched= mean(match_c);
std_means(:,1)= (mean_t-mean_c_unmatched)./std_treatedIndexed;
std_means(:,2)= (mean_t-mean_c_matched)./std_treatedIndexed;

% Remove cases where the std dev of the treatedIndexed is close to 0
std_means(std_means<-1000 |std_means>1000)=0;
std_means(std_means>1000)=0;

var_t= var(visit1Matt)';
var_uc= var(unmatched_c)';
var_mc= var(match_c)';
var1= [var_t var_uc];
var2= [var_t var_mc];

maxi= max([var_t; var_uc; var_mc]);
mini= min([var_t; var_uc; var_mc]);

figure
subplot(1,3,1)
parallelcoords(std_means); 
title("Change in Standardized means of difference after " +type+ " matching")

subplot(1,3,2)
parallelcoords(var1); 
ylim([mini maxi])
title("Var- treatedIndex vs control")

subplot(1,3,3)
parallelcoords(var2);
ylim([mini maxi])
title("Var- treatedIndex vs matched control")