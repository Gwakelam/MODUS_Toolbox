function [reshape_v] = Tstat_reshape_v1(v, Fs, T)
%
% Reshape 2D arrays of continuous data into 3D arrays of individual
% Tstat (period of stationarity) ensembles. 
%
%  INPUTS
%   [v]: A 2D matrix containing velocity measurement data over a number of depth bins
%              (Size = nt x nz)
%   [Fs]: Sampling frequency (Hz)
%              (Size: 1 x 1)
%   [T]:    Length of Tstat period of stationarity (in seconds)   
%
%  OUTPUTS:
%   [reshape_v]: 3D matrix of T-stat ensembles
%               (Size = ns x ne x nz)
%
%  WHERE:   nt = total number of samples in data set
%           ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data set
%           nz = number of depth/distance cells from sensor
%
% v1    03/04/2018      First Version
% v1    28/05/2018      Added check for complete ensemble, set to discard incomplete set at end of data    
%################################################################################################

% Calc sizes/dimensions
ns = round(T*Fs);
nt = size(v, 1);
ne = (nt/ns);

% Check for 'fit' - if there is an incomplete Tstat at end of dataset, drop
% partial ensemble.
if rem(ne, 1) ~= 0
    v = v(1 : floor(ne)*ns, :) ;
end
    
% Reshape
reshape_v = reshape(v, ns, [], size(v, 2));


end

