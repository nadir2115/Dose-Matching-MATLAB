function [h,p, treated_effect, control_effect]= ttestImprovement (treat, matched_control, dv1, dv2)
%Assuming DV1/2 has 27 rows
dv1t =dv1(treat,:); 
dv2t =dv2(treat,:); 

%Picking out matched control values in Dataset 2
dv1_matched_control= dv1(matched_control(:),:);

%Picking out matched control values in Dataset 2
dv2_matched_control = dv2((matched_control(:)),:);

% %effect on log_mean_time_MA_PA
% treated_effect= dv1t(:,9)-dv2t(:,9);    
% control_effect= dv1_matched_control(:,9)-dv2_matched_control(:,9);

%effect on SIS hand
treated_effect= dv1t(:,2)-dv2t(:,2);    
control_effect= dv1_matched_control(:,2)-dv2_matched_control(:,2);

%T-test-----------------------------------------------------------------
% Test null hypothesis: pairwise difference between treat/control has mean equal to 0
[h,p] = ttest(treated_effect,control_effect);
% h=1 indicates that ttest rejects null hypothesis at the default 5% sig level
