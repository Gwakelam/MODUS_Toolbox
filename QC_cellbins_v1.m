function [celldist, cellfail] = QC_cellbins_v1(cellwidth, cellblank, nz)

% Checks that the sensor configuration for distance bins does not change during a Tstat
% If the bin configuration changes, sets all measured values in the Tstat to NaN
%
%  INPUTS
%   [cellwidth]: Array containing data on cell/bin width for
%               each sample, arranged in Tstat ensembles
%              (Size = ns x ne)
%   [cellblank]: Array containing data on cell/bin blanking distance for
%               each sample, arranged in Tstat ensembles
%              (Size = ns x ne)
%   [nz]:       number of distance cells/bins in dataset
%               (Size = 1x1)
%
%  OUTPUTS:
%   [cells_out]: Array of width/blanking dis, ensembles with more than one 
%                   value are 'killed' (set to NaN)
%               (Size = ns x ne)
%
%  WHERE:   ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    07/06/2018      GW - First Version
%
[ns, ne] = size(cellwidth);
celldist = NaN(nz, ne);
cellfail = zeros(1, ne);
binvec = 0:1:nz-1;

for i = 1:ne
   if length(unique(cellwidth(:,i))) ~= 1  || length(unique(cellblank(:,i))) ~= 1
       cellfail(i) = true; 
   else
       celldist(:,i) = binvec.*cellwidth(1,i) + cellblank(1,i) ;
       
   end   
    
end