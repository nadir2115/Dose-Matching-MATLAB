function assessBalance(treat, control, matched_control, type, dv1)
dv1t= dv1(treat,1:end-1);
mean_t= mean(dv1t); 
std_treated= std(dv1t);
unmatched_c= dv1(control,1:end-1);
mean_c_unmatched= mean(unmatched_c);

match_c= dv1(matched_control(:,1),1:end-1);
mean_c_matched= mean(match_c);
std_means(:,1)= (mean_t-mean_c_unmatched)./std_treated;
std_means(:,2)= (mean_t-mean_c_matched)./std_treated;

% Getting rid of cases where the standard deviation of the treated is close
% to 0
std_means(std_means<-1000 |std_means>1000)=0;
std_means(std_means>1000)=0;

var_t= var(dv1t)';
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
title("Var- treat vs control")

subplot(1,3,3)
parallelcoords(var2);
ylim([mini maxi])
title("Var- treat vs matched control")