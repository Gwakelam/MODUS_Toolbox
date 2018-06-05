function [ D, r ] = structureFunction_v1( fluc_v , z, options)
%
% Function to calculate the spatial structure function of a de-trended ADCP time series
% Based on function by J Thomson - see https://github.com/mguerrap/5Beam-Turbulence-Methods/blob/master/structureFunction.m
%
%  INPUTS:
%   [fluc_v]:   A 3D array containing detrended velocity measurements in 5-minute T-stat ensembles 
%              (Size = ns x ne x nz)
%
%   [z]:        A vector of z bin distances from sensor. 
%
%   [options]   Struct of optional parameters
%       options.surf_blank - distance below free surface to be 'ignored'
%
%
% OUTPUTS:
%   [D]:        Array of structure functions for each Tstat ensemble at
%               each depth bin (Size = nr x ne x nz)
%   [r]:        Array of radial separation values for each z bin
%               (Size = nr x nz)
%
%  WHERE:  
%           nr = number of radial separations ( range: 0 - nz, ie. nr = nz)
%           ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% J. Thomson, 7/2009, 
%       rev. 7/2010 (efficiency) 
%       rev. 8/2010 (explicitly remove mean [should be done already], allow nans in profiles) 
%       rev. 9/2011  (return signed r value, to allow upward or downward preference in D(z,r) fit)
% M. Guerra Dec 8 2015, change r deffinition to uplooking
% v1    18/05/2018  G Wakelam - Adapted for 3d matrix, aligned dimensions to other ReDAPT processing scripts
%       

[time, ensembles, bins] = size(fluc_v);

bins2 = length(find(z<= (z(end)- options.surf_blank)));


 for i = 1:bins
        
        r(:,i) = round(z(:) - z(i)); % For uplooking z1 is the smallest v(z+r)-v(z), then r=(z+r)-z;
        %For downlooking r(i,:) = z(i) - z(:);
 end

for j = 1:ensembles;

    temp_v = squeeze(fluc_v(:, j,:));
    
for i = [1:bins],
  
        for t = 1:time, 
            % Here the sign of the dzrt is not important since it is
            % squared: v(z+r)-v(z)
            
        dzrt(t,:) = ( temp_v(t,:) - temp_v(t,i) ).^2; % do the same for velocities for each time
        
        end
        
    D(:, j, i) = nanmean(dzrt,1);       % Ojo que da vuelta el vector nanmean(dzrt,2) para tener D(z,r)
    
end

end
