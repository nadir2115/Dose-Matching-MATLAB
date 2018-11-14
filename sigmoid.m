function g = sigmoid(z)
%SIGMOID Compute sigmoid function
%   g = SIGMOID(z) computes the sigmoid of z.

% You need to return the following variables correctly 
g = [];

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the sigmoid of each value of z (z can be a matrix,
%               vector or scalar).
% for i=1:size(z,1) %Clunnky code- DON'T USE
%     for j=1:size(z,2)
%         g(i,j) = 1/(1+exp(-z(i,j)));
%     end
% end

g = 1.0 ./ (1.0 + exp(-z));
% =============================================================

end
