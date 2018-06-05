function [prod,shear] = calc_Prod_v1(mean_vb,R,z,options)
%
%  INPUTS:
%    [mean_vb]: A structure containing the mean values of the
%               relevant velocities (Size = ne x nz)
%         mean_vb.x - velocity along x-axis of ADCP 
%         mean_vb.y - velocity along y-axis of ADCP    
%         mean_vb.mag - A signed horizontal velocity (Can set to NaN)
%          **(Note: Need to verify that x and y are based on RDI convention)**
%
%    [R]:  Structre of Reynolds stresses:
%           R.xz - xz Reynolds stress  (Size = ne x nz)
%           R.yz - yz Reynolds stress  (Size = ne x nz)  
% 
%    [z]: A vector with the bin heights above bottom (Size = nz x 1)
%      Note: z values are assumed to be uniform
% 
%    [options]: A structure containing optional parameters
%         options.figure - Flag to create figures (Default = 1)
%
% OUTPUTS:
%    [prod]: A matrix of the production values (Size = ne x nz)
% 
%    [shear]: A structure with the shear in instrument coordinates (Size = ne x nz)
%         shear.dvxdz - Vertical shear of x-velocity
%         shear.dvydz - Vertical shear of y-velocity
%
% WHERE:
%        ne = number of Tstat ensembles in data
%        nz = number of depth/distance cells from sensor
%
% JMM
% Nov 25, 2016 -    Code created
% Sep 20, 2017 -    Code cleaned up for dissemination
% v1    16/05/2018  Adapted from original code, isolated production
%                   calculation - GW

disp('Computing production')

% figure
if ~isfield(options,'figure')
    options.figure = 1;
end


%% Define simple variables
U = mean_vb.mag;
vx = mean_vb.x;
vy = mean_vb.y;


%% Calculate shear

% initialize
[ne, nz] = size(U);
dUdz = NaN*ones(ne,nz);
dvxdz = NaN*ones(ne,nz);
dvydz = NaN*ones(ne,nz);

% calculate shear (central difference)
dz = mean(diff(z));
for zz = 2:nz-1 % neglect top and bottom bins
    dUdz(:,zz) = (U(:,zz+1)-U(:,zz-1))/(2*dz);
    dvxdz(:,zz) = (vx(:,zz+1)-vx(:,zz-1))/(2*dz);
    dvydz(:,zz) = (vy(:,zz+1)-vy(:,zz-1))/(2*dz);
end


%% compute production
Px = -R.xz.*dvxdz;
Py = -R.yz.*dvydz;

%% variables to export
prod = Px + Py;

shear.dvxdz = dvxdz;
shear.dvydz = dvydz;

disp('----------------------------')

%% plots
z = squeeze(z);
if options.figure
    %----------
    % plot shear
    %----------
    figure(1),clf
    ax(1)=subplot(411);
    pcolor(1:ne,z,U'), shading flat
    colorbar
    caxis(2.2*[-1 1])
    ylabel('U')
    title('veloctiy shear')
    
    ax(2)=subplot(412);
    pcolor(1:ne,z,dUdz'), shading flat
    colorbar
    caxis(0.12*[-1 1])
    ylabel('dU/dz')
    
    ax(3)=subplot(413);
    pcolor(1:ne,z,dvxdz'), shading flat
    colorbar
    caxis(0.12*[-1 1])
    ylabel('dv_x / dz')
    
    ax(4)=subplot(414);
    pcolor(1:ne,z,dvydz'), shading flat
    colorbar
    caxis(0.12*[-1 1])
    ylabel('dv_y / dz')
    
    linkaxes(ax,'x')
    
    %--------
    % production
    %--------
    figure(3),clf, clear ax
    
    ax(1)=subplot(211);
    pcolor(1:ne,z,Px'),shading flat
    colormap(gca,flipud(hot))
    caxis([-0.00005 0.3e-3])
    colorbar
    title('P_x')
    
    ax(2)=subplot(212);
    pcolor(1:ne,z,Py'),shading flat
    colormap(gca,flipud(hot))
    caxis([-0.00005 0.3e-3])
    colorbar
    title('P_y')
    
end
