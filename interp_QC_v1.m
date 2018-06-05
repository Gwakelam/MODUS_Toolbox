function [Uout, Uf_out] = interp_QC_v1(U, badlim, maxgap)
% interpolates over individual gaps in sensor data, along 1st dimension of
% array. 'kills' Tstat ensembles with over (badlim)% of NaNs remaining after
% interpolation
%
%  INPUTS
%   [U]: Array containing sensor data arranged into Tstat ensembles
%              (Size = ns x ne x nz)
%   [badlim]: limit of 'bad' values (NaN, spikes??) in any one ensemble (as
%             a percent)
%               (Size = 1)
%   [maxgap]: maximum gap in data that will be interpolated (no. of samples)
%               (Size = 1)
%
%  OUTPUTS:
%   [Uout]: Array with single gaps in ensembles interpolated, ensembles
%           with more than badlim% 'bad values' are 'killed'
%               (Size = ns x ne x nz)
%   [Ubad]: number of 'bad' values in each ensemble (before interpolation)
%               (size = ne x nz)
%
%  WHERE:   ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    04/04/2018      First Version
% v1    28/05/2018      GW - Updated to chceck for ensembles with >20%
%                       NaN's. Will set all to NaN instead of tring to interpolate
%

% TO DO:
%
% Add filter for upper limit on interpolated values in a single Tstat?
% Allow selection of interpolation method?
%
%################################################################################################
%% MAIN

%testing
%badtest = sum(isnan(Data.ADCP.ADCP01_NW_Dep3.B1));
%badlim = 10;

% Variable Setup
L       = size(U);
badnum  = (badlim/100)*L(1);
Uout    = zeros(L);
maxgap  = maxgap + 1;            % Gap actually interpolated is n-1.... 

Ubad_in = squeeze(sum(isnan(U)));
Uf_in   = Ubad_in > badnum;

% Loop through all Ensembles... :'(
for i = 1:L(2)
    for j = 1:L(3)
        
        Utemp = U(:,i,j);
        
        % If missing val at start or end, replace with ensemble mean
        if isnan(Utemp(1))
            Utemp(1) = nanmean(Utemp);
        end
        if isnan(Utemp(end))
            Utemp(end) = nanmean(Utemp);
        end
        
        % Kill ensembles with >10% NaNs, esle interpolate single gaps
        if sum(isnan(Utemp)) > round(0.10 * length(Utemp))
            Uout(:,i,j) = NaN(size(Utemp));
        else
            Uout(:,i,j) = interp1gap(Utemp, maxgap);
        end
        
    end
end

% Sum NaNs in each Tstat ensemble, find those above limit
Ubad_out    = squeeze(sum(isnan(Uout)));
Uf_out      = Ubad_out > badnum;

% Get indices of failed ensembles
[ind2, ind3] = find(Uf_out);
 % 'kill' failed ensembles
Uout(:, ind2, ind3) = NaN;

%test plot
%imagesc(UF_out)

