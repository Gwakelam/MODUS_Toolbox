function [TKE] = calc_TKE_Osalusi_v1(std_vb, std_noise, options)
% Based on function provided by Justine Mcmillan
%
% Calculate TKE using a variance method approach as described by Osalusi(2009)
%  - TKE = C * [sum(beam variances - noise variances)]^2
%      
%    where C is 
%      
%    C = 1/(16*sind(theta)^4*(1+xi*(2*cotd(theta)^2+1)))
%   
%    Note: C is a constant dependent on the 
%          empirical value xi (= 0.17 from Nezu & Nakagawa (1993)) and beam angle (phi)
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
%    [R]: A structure of Reynolds stresses for tilt correction - NOT YET IMPLEMENTED
%               (Size = ne x nz)
%         R.xz - 
%         R.yz - 
% 
%    [options]: A structure containing optional parameters
%         options.theta - Beam angle with respect to vertical (Default = 20)
%         options.xi0 - Anisotropic ratio (Default = 0.17)
%         options.tilt - Flag to enable tilt correction using reynolds stresses (Default = 0) (to be implemented)
%         options.figure - Flag to create figures (Default = 0)
%
% OUTPUTS:
%    [TKE]:     A matrix of the TKE density values (Size = ne x nz)
%
%   WHERE:
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    16/05/2018      First Version - adapted from TKE_Lu_v1
%
%################################################################################################

disp('Computing TKE...')
%% Set defaults

% beam angle
if ~isfield(options,'theta'); 
    options.theta = 20;
    disp(['beam angle set to ', num2str(options.theta)])
end

% Xi - semi empirical constant. Used by Osalusi, value from Nesun & nakagawa.
if ~isfield(options,'xi0'); 
    options.xi0 = 0.17;
    disp(['Empirical constant Xi set to ', num2str(options.xi0)])
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
xi = [0:0.01:0.5];
const = 1./(16*sind(options.theta)^4*(1+xi*(2*cotd(options.theta)^2-1)));

figure(1),clf
plot(xi,const)
hold all
plot(options.xi0,interp1(xi,const,options.xi0),'ok')
xlabel('\alpha')
ylabel('C')

%% Compute from ADCP data
const0 = 1/(16*sind(options.theta)^4*(1+options.xi0*(2*cotd(options.theta)^2-1)));

% Create matrix of noise bias values for easy removal
for bb = 1:4
    vb = ['B',num2str(bb)];
    std_noise_mat.(vb) = ones(nt,1) * std_noise.(vb);
end


% calc TKE (This part doesn't change between Oslusi & Lu/Lueck methods)
TKE = const0*((std_vb.B1.^2-std_noise_mat.B1.^2+...
    std_vb.B2.^2-std_noise_mat.B2.^2+...
    std_vb.B3.^2-std_noise_mat.B3.^2+...
    std_vb.B4.^2-std_noise_mat.B4.^2).^2);
% No noise correction
TKE0 = const0*((std_vb.B1.^2+ std_vb.B2.^2+...
    std_vb.B3.^2+std_vb.B4.^2).^2);


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



end