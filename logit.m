%function to find logit of number or vector elements
function y= logit(x)
y= log10(x./(1-x)); 