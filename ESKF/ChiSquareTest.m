function [isOutlier,H_out] = ChiSquareTest(H_in,P_hat, R, delta_y)
%CHISQUARETEST Summary of this function goes here
%   Detailed explanation goes here
thresh = 0;
dim = length(delta_y); 

% Threshold corresponds to a confidence interval of 95%
if (dim == 1)
    thresh = 3.841;
elseif (dim == 2)
    thresh = 5.991;
elseif (dim == 3)
    thresh = 7.815;
end

S = H_in * P_hat * H_in' + R;
T = delta_y' * (S^(-1)) * delta_y;
if (T > thresh) 
    isOutlier = true;
    H_out = zeros(dim,18);
else
    isOutlier = false;
    H_out = H_in;
end

end

