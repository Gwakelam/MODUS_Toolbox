function [R] = calc_Rstress_v1(std_vb, options)
% From Justine Mcmillan scripts
%
% Calculate TKE using a variance method approach
%  - TKE = C * sum(beam variances - noise variances)
%      
%    where C is 
%      
%    C = (alpha+1)/(1+2*alpha/tand(theta)^2)/(4*sind(theta)^2
%   
%    Note: C is a constant dependent on the 
%          anisotropy ratio (alpha) and beam angle (theta)
%
%   Calculate Reynolds Stresses Rxz & Ryz 
%
%  INPUTS:
%    [std_vb]: A structure containing the standard deviations
%              of the along-beam variances (Size = ne x nz)
%         std_vb.B1 - Beam 1      
%         std_vb.B2 - Beam 2 
%         std_vb.B3 - Beam 3      
%         std_vb.B4 - Beam 4   
% 
%    [options]: A structure containing optional parameters
%         options.theta     - Beam angle with respect to vertical (Default = 20)
%         options.heading   - array of sensor heading for each Tstat
%                             ensemble (or single value if constant) (Defalut = No correction)
%         options.figure    - Flag to create figures (Default = 0)
%
% OUTPUTS:  
%    [R]:  Structre of Reynolds stresses:
%           R.xz - xz Reynolds stress  (Size = ne x nz)
%           R.yz - yz Reynolds stress  (Size = ne x nz)
%
%   WHERE:
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    16/05/2018      First Version - created by separating
%                       calc_TKE_Rstress_v2.m
% v1    29/05/2018      GW - added transform to flow axes. See Guerra, 2017: 
%                       https://github.com/mguerrap/5Beam-Turbulence-Methods/blob/master/RS_VT.m
%
%################################################################################################

%% Set defaults

% beam angle
if ~isfield(options,'theta'); 
    options.theta = 20;
    disp(['beam angle set to ', num2str(options.theta)])
end

% beam angle
if ~isfield(options,'heading'); 
    flag_heading = 0;
    disp(['heading correction NOT applied.'])
elseif isfield(options,'heading'); 
    flag_heading = 1;
    disp(['heading correction applied.'])
end

% figure
if ~isfield(options,'figure'); 
    options.figure = 0;
end


%% Correlation calculations
vpxvpz = - (std_vb.B2.^2 - std_vb.B1.^2)/(4*cosd(options.theta)*sind(options.theta));
vpyvpz = - (std_vb.B4.^2 - std_vb.B3.^2)/(4*cosd(options.theta)*sind(options.theta));

R.xz = vpxvpz;
R.yz = vpyvpz;

% Transform from sensor axes to flow axes
if flag_heading
    % Heading correction as per Guerra: https://github.com/mguerrap/5Beam-Turbulence-Methods/blob/master/RS_VT.m 
    R.uw = bsxfun(@times, R.xz, cosd(heading)) + bsxfun(@times, R.yz, sind(heading));
    R.vw = bsxfun(@times, R.yz, cosd(heading)) - bsxfun(@times, R.xz, sind(heading));
end

end