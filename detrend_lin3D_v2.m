function [detr_v, mean_v, std_v] = detrend_lin3D_v2(v, dim)
%
% Linear detrend velocity data along any dimension ignoring NaN's, and find standard 
% deviation of T-stat ensembles. Based on detrendNaN3 function be Leander Moes, available at:
% https://uk.mathworks.com/matlabcentral/fileexchange/63749-detrendnan3
%
%  INPUTS:
%   [v]: A 3D matrix containing velocity measurements in 5-minute ensembles 
%              (Size = ns x ne x nz)
%   [dim]: Dimension of 3D matrix along which linear detrend will operate(Optional)
%              (Size: 1 x 1)
%         
%  OUTPUTS:
%   [detr_v]: 3D matrix of detrended T-stat ensembles
%               (Size = ns x ne x nz)
%   [mean_v]: array of mean velocity values for each T-stat ensemble
%               (Size = ne x nz)
%   [std_v]: array of standard deviation of values in each T-stat
%               (Size = ne x nz)
%
%  WHERE:   ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    21/03/2018      First Version
% v2    03/05/2018      Added mean and standard dev to outputs - GW
%################################################################################################
%

% Deal with dim being undefined
if nargin < 2
    dim = 1;
end

if isempty(dim)
    dim = 1;
end

% create default vector for sample spacing (ie. time - assume even)
t = 1:size(v,1);

% time vector to same format as v and same NaN positions
t = bsxfun(@times,permute(t,[2 1 3]),ones(size(v)));
t(isnan(v)) = NaN;

%% Calculation

% mean of time
xm = nanmean(t,dim);
% mean of A
ym = nanmean(v,dim);
% calculate slope using least squares
a = nansum(bsxfun(@times,bsxfun(@minus,t,xm),bsxfun(@minus,v,ym)),dim)./nansum(bsxfun(@minus,t,xm).^2,dim);
% calculate intercept
b = ym - a.*xm;
% calculate trend
trend = bsxfun(@plus,b,bsxfun(@times,a,t));
% remove trend
detr_v = v-trend;

mean_v = squeeze(nanmean(v, dim));

std_v = squeeze(std(v, 0, dim, 'omitnan'));
end


% figure
% plot(detrended(:, 200, 20))
