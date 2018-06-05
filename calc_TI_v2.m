function [TI] = calc_TI_v2(std_v, noise_v, mean_v, options)
%% [TI, metadata] = calc_TI(std_v,noise_v, mean_v, options)
%
% To calculate TI using variance method(?) 
%  - TI = (velocity variances - noise variances) / mean velocities
%  - TI = (std_v - noise_v) ./ (mean_v)
%
%  - NOTE: naming convention of velcoity directions must be consistent in
%    all input structs, but do not need to be in same order. 
%
%  INPUTS:
%    [std_v]: A structure containing the 5-minute ensemble standard deviations
%              of the detrended velocity (Size = ne x nz)        
% 
%    [noise_v]: A structure containing the signal standard deviation
%                   due to Doppler noise (Size = ne x nz, can also be 1 x nz 
%                   for depth varying or 1 x 1 for constant/fixed noise correction) 
%
%    [mean_v]: A structure containing the 5 minute ensemble mean velocity
%                   (Size = ne x nz)
% 
%    [options]: A structure containing optional parameters
%         options.slack         - set max velocity for slack water: TI values will not be 
%                                   calculated for Tstat's with mean velocity below this value
%         options.figure        - Flag to create figures (Default = 1)
%         options.hubheight     - set hub height of 'device' for figures
%
% OUTPUTS:
%    [TI]: A matrix of the TI values (Size = ne x nz)
%
% WHERE:
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% Versions
% Mar 12 2018   v1  Code created
% Apr 04 2018   v2  Removed struct format of input variables - only
%                   processes 1 dataset/direction at a time - GW
% May 24 2018   v2  Added minimum slack velocity cutoff to avoid large TI values - GW
%################################################################################################

%% Default Options
if ~isfield(options, 'slack')
    options.slack = 0.5;
    disp('Slack water not defined in options, setting to <= 0.5 m/s')
end

%% Dimension Check

    % Check dimensions match those of std_v data
    d_chk = size(std_v) ;
    
    if size(noise_v, 1) ~= d_chk(1) && size(noise_v, 1) ~=1
        error('[ERROR] - 1st dimension of noise_v is not 1 or equal to that of std_v')
    end
    
    if size(noise_v, 2) ~= d_chk(2) && size(noise_v, 2) ~=1
        error('[ERROR] - 2nd dimension of noise_v is is not 1 or equal to that of std_v')
    end
    
    if size(mean_v, 1) ~= d_chk(1)
        error('[ERROR] - 1st dimension of mean_v is not equal to that of std_v')
    end
    
    if size(mean_v, 2) ~= d_chk(2)
        error('[ERROR] - 2nd dimension of mean_v is not equal to that of std_v')
    end
   
    %% Calc TI
    
    % remove corrected variance values that are negative (occasionally
    % occur around free surface)
    corrected_v = (bsxfun(@minus, (std_v.^2), (noise_v.^2)));
    corrected_v(corrected_v < 0 ) = NaN;
    
    mean_v(abs(mean_v) <= options.slack) = NaN;
    
    % calculate TI - as standard dev has been squared to get variance, must
    % square root. Need to check that same numbers come out of std dev used
    % throughout
    TI =    (sqrt(corrected_v)./ abs(mean_v))*100;
    
end

% Add plotting options here
% if options.figure == 1
%   figure
%   plot( ... TI profile)
%   figure
%   plot( ... TI timeseries @ HH)
%   etc
% end




