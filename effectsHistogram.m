function effectsHistogram(treated_effect, control_effect, method)

range= [treated_effect; control_effect];
figure
subplot(2,1,1)
hist(treated_effect)
xlim([min(range) max(range)]);

title("Histogram of improvement for treated subjects")
subplot(2,1,2)
hist(control_effect)
xlim([min(range) max(range)]);
title("Histogram of improvement for " + method + " matched control subjects")