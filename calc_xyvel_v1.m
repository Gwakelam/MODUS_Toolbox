function [sens, sens_mean] = calc_xyvel_v1(vb, tilt, options)
%
%  INPUTS:
%    [mean_vb]: A structure containing measured along-beam velocities (Size = ns x ne x nz)
%         vb.B1 - velocity along x-axis of ADCP 
%         vb.B2 - velocity along y-axis of ADCP
%         vb.B3 - velocity along x-axis of ADCP 
%         vb.B4 - velocity along y-axis of ADCP
%          **(Note: Need to verify that x and y are based on RDI convention)**
%
%    [tilt]: A structure containing values for pitch, roll & heading over
%    data period (size = ne)
%         tilt.p - pitch about x axis
%         tilt.r - roll about y axis
%         tilt.h - heading avout z axis
% 
%    [options]: A structure containing optional parameters
%         options.theta  - beam angle from vertical (Defult = 1)
%         options.figure - Flag to create figures (Default = 1)
%
% OUTPUTS:
%    [sens]: A matrix of the sensor-axes resolved velocities 
%         sens.x  - measured velocity along sensor x-axis (Size = ns x ne x nz)
%         sens.xm - mean velocity along to sensor x-axis  (Size = ne x nz)
%         sens.y  - measured velocity along sensor y-axis (Size = ns x ne x nz)
%         sens.ym - mean velocity along to sensor y-axis  (Size = ne x nz)
% 
%    [sens_mean]: signed magnitude of horizontal velocity
%
% WHERE:
%        ns = number of samples in a period of stationarity
%        ne = number of Tstat ensembles in data
%        nz = number of depth/distance cells from sensor
%
%CHANGES
% 24/05/18  GW - added sens_mean output struct for mean values
%

%% Set defaults

% beam angle
if ~isfield(options,'theta'); 
    options.theta = 20;
    disp(['beam angle set to ', num2str(options.theta)])
end


% figure
if ~isfield(options,'figure'); 
    options.figure = 0;
end

%% Calculate

theta = options.theta;

sens.x    = (vb.B2 - vb.B1)./(2*sind(theta)) - bsxfun(@times, tilt.p, ((vb.B1 + vb.B2)./(2*cosd(theta))));
sens.y    = (vb.B4 - vb.B3)./(2*sind(theta)) - bsxfun(@times, tilt.r, ((vb.B3 + vb.B4)./(2*cosd(theta))));
sens.z    = -1.*(vb.B2 + vb.B1 + vb.B3 + vb.B4)./(4*cosd(theta)) ...
            - bsxfun(@times, tilt.p, ((vb.B2 - vb.B1)./(2*sind(theta))))...
            + bsxfun(@times, tilt.r, ((vb.B4 - vb.B3)./(2*sind(theta))));
        
        
sens.mag       = sqrt(sens.x.^2 + sens.y.^2);      %JUST USE U FROM DATA - ALREADY CALCULATED

sens_mean.x = squeeze(nanmean(sens.x, 1));
sens_mean.y = squeeze(nanmean(sens.y, 1));
sens_mean.z = squeeze(nanmean(sens.z, 1));
sens_mean.mag = squeeze(nanmean(sens.mag, 1));

end