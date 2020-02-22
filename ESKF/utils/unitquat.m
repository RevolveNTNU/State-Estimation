function q = unitquat(eps)
% UNITQUAT - computes the unit quaternion having the input epsilon vector as vector part 
% Utilizing the property that sqrt(eta^2 + eps1^2 + eps2^2 + eps3^2) = 1

eta = sqrt(1- eps'*eps);
eta = abs(eta);
q = [eta; eps];

end

