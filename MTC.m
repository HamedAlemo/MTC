function [RMSE, rho2] = MTC(X)
% Multiplicative Triple Collocation (MTC) calculates the root-mean-squared-
% error (RMSE) and correlation coefficient between three sets of data (from 
% different measurement systems) and the true. The difference between MTC
% and classical TC is that in MTC the error model is multiplicative which
% is more appropriate in case of e.g. precipitation data.
%
%
%%% Input:
%
%       X    is a matrix of N * 3 in which N is the number of spatially and 
%            temporally collocated data points. Each column of y represents one
%            data set. Zero values and missing values should have been removed.
%
%%% Outputs:    
%
%       RMSE is a 1 * 3 vector representing the RMSE between each of the data 
%            sets and the true.
%
%       rho2  is a 1 * 3 vector representing the squared correlation coefficient 
%             between each of the data sets and the true.
%
%
%%% Reference
%
% Alemohammad S.H., McColl K.A., Konings A.G., Entekhabi D., Stoffelen A. 
% (2015), Characterization of Precipitation Product Errors across the 
% United States using Multiplicative Triple Collocation, Hydrology and 
% Earth System Sciences, 19. 
%
%
%
%
% Version: 1.0, July 2015.
% Author:  S. Hamed Alemohammad, hamed_al@mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Check input data
if size(X, 2) ~= 3
    error('Error: Input data X must be an N x 3 matrix');
end

if any(isnan(X(:)))
    error('Error: Input data X must not contain NaNs');
end
    

% Log transformation of data
X = log(X);

% Calculating the covarinace matrix 
C = nancov(X);

% Bootstrap
  opt = statset('UseParallel', true);   % Use false if you don't want to use parallel processing toolbox.
  N_Boot = 1000;                        % Number of samples for Bootstrap
  C_bootstrap = bootstrp(N_Boot, 'nancov', X, 'Options', opt);
  sigma_t = zeros(N_Boot, 3);
  for n = 1 : N_Boot
	  CP = reshape(C_bootstrap(n, :), 3, 3);
	  sigma_t(n,1) = (CP(1,1) - ((CP(1,2) * CP(1,3)) / CP(2,3)));
	  sigma_t(n,2) = (CP(2,2) - ((CP(1,2) * CP(2,3)) / CP(1,3)));
	  sigma_t(n,3) = (CP(3,3) - ((CP(1,3) * CP(2,3)) / CP(1,2)));
  end
% End of Bootstrap

sigma_t = sqrt(sigma_t);

RMSE(1, 1) = mean(sigma_t(:, 1));
RMSE(1, 2) = mean(sigma_t(:, 2));
RMSE(1, 3) = mean(sigma_t(:, 3));

rho2(1, 1) = (C(1,2) * C(1,3)) / (C(1,1) * C(2,3));
rho2(1, 2) = (C(1,2) * C(2,3)) / (C(2,2) * C(1,3));
rho2(1, 3) = (C(1,3) * C(2,3)) / (C(3,3) * C(1,2));
