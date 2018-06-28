function [TKE] = calc_TKE_Lu_v1(std_vb, std_noise, options)
% From Justine Mcmillan scripts
%
% Calculate TKE using a variance method approach as described by Lu & Lueck
% (1999)
%  - TKE = C * sum(beam variances - noise variances)
%      
%    where C is 
%      
%    C = (alpha+1)/(1+2*alpha/tand(theta)^2)/(4*sind(theta)^2
%   
%    Note: C is a constant dependent on the 
%          anisotropy ratio (alpha) and beam angle (phi)
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
%    [std_noise]: A structure containing the standard deviations
%                 of the Doppler noise (Size = 1 x nz)
%         std_noise.B1 - Beam 1      
%         std_noise.B2 - Beam 2 
%         std_noise.B3 - Beam 3      
%         std_noise.B4 - Beam 4 
%
%    [R]: A structure of Reynolds stresses for tilt correction - NEEDS TO BE ADDED
%               (Size = ne x nz)
%         R.xz - 
%         R.yz
% 
%    [options]: A structure containing optional parameters
%         options.phi - Beam angle with respect to vertical (Default = 20)
%         options.alpha0 - Anisotropic ratio (Default = 0.2)
%         options.figure - Flag to create figures (Default = 0)
%
% OUTPUTS:
%    [TKE]:     A matrix of the specific TKE density values [m^2/s^2] -
%               multiply by density of seawater to got [J/m^3]
%               (Size = ne x nz)
%
%   WHERE:
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    04/04/2018      First Version
% v2    03/05/2018      Changed input/output structure to match other
%                       functions (dim 2 is nz)
%################################################################################################

disp('Computing TKE...')
%% Set defaults

% beam angle
if ~isfield(options,'phi'); 
    options.phi = 20;
    disp(['beam angle set to ', num2str(options.phi)])
end

% anisotropy ratio
if ~isfield(options,'alpha0'); 
    options.alpha0 = 0.2;
    disp(['anisotropy ratio Alpha set to ', num2str(options.alpha0)])
end

% figure
if ~isfield(options,'figure'); 
    options.figure = 0;
end

%% 
[nt,nz] = size(std_vb.B1);
if nz ~= length(std_noise.B1);
    error('std_noise length does not match std_vb dim 2')
end


%% Sensitivity to alpha
alpha = [0:0.01:0.5];
const = (alpha+1)./(1+2*alpha/tand(options.phi)^2)/(4*sind(options.phi)^2);

figure(1),clf
plot(alpha,const)
hold all
plot(options.alpha0,interp1(alpha,const,options.alpha0),'ok')
xlabel('\alpha')
ylabel('C')

%% Compute from ADCP data
const0 = (options.alpha0+1)/(1+2*options.alpha0/tand(options.phi)^2)/(4*sind(options.phi)^2);

% Create matrix of noise bias values for easy removal
for bb = 1:4
    vb = ['B',num2str(bb)];
    std_noise_mat.(vb) = ones(nt,1) * std_noise.(vb);
end


% calc TKE (This part doesn't change between Oslusi & Lu/Lueck methods)
TKE = const0*(std_vb.B1.^2-std_noise_mat.B1.^2+...
    std_vb.B2.^2-std_noise_mat.B2.^2+...
    std_vb.B3.^2-std_noise_mat.B3.^2+...
    std_vb.B4.^2-std_noise_mat.B4.^2);
TKE0 = const0*(std_vb.B1.^2+ std_vb.B2.^2+...
    std_vb.B3.^2+std_vb.B4.^2);


if options.figure
    figure(3),clf
    subplot(211)
    pcolor(TKE)
    shading flat
    colorbar
    caxis([0 0.05])
    title('TKE - noise corrected')
    xlabel('year day'), ylabel('z [m]')

    subplot(212)
    pcolor(TKE0)
    shading flat
    colorbar
    caxis([0 0.05])
    title('TKE - No noise correction')
    xlabel('year day'), ylabel('z [m]')
end



%% Correlation calculations - alternative to passing R as input var for tilt correction?

% vpxvpz = - (std_vb.v2.^2 - std_vb.v1.^2)/(4*cosd(options.phi)*sind(options.phi));
% vpyvpz = - (std_vb.v4.^2 - std_vb.v3.^2)/(4*cosd(options.phi)*sind(options.phi));
% 
% R.xz = vpxvpz;
% R.yz = vpyvpz;


end