function [LSt] = Lscale_transform_v1(LS, orient, theta)
%% Transform along-beam ADCP lengthscales to flow axes
% Use transform given by RDI instruments 
%
% Transform to sensor axes (x, y, z) THEN transform to flow axes using
% transforms in RDI manual.

%  INPUTS:
%   [LS]:   A struct containing along-beam ADCP lengthscales for each beam 
%           LS.B1 - Lengthscales along beam 1 (Size = ns x ne OR 1 x ne)
%           LS.B2 - Lengthscales along beam 2 (Size = ns x ne OR 1 x ne)
%           LS.B3 - Lengthscales along beam 3 (Size = ns x ne OR 1 x ne)
%           LS.B4 - Lengthscales along beam 4 (Size = ns x ne OR 1 x ne)
%
%
%   [orient]:   Struct containing orientations of the ADCP sensor
%           orient.p - pitch of sensor (Size = ne x1)
%           orient.r - roll of sensor (Size = ne x1)
%           orient.h - heading of sensor (Size = ne x1)
%
%   [theta]:    Beam angle from vertical (in degrees).
%              
%         
%  OUTPUTS:
%   [LSt]:   Struct of lengthscales transformed to flow axis 
%           Lscale.U - Lengthscale in U direction (Size = ne x nz)
%           Lscale.V - Lengthscale in V direction (Size = ne x nz)
%           Lscale.W - Lengthscale in W direction (Size = ne x nz)
%
%  WHERE:  
%           ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    29/05/2018      First Version

%% Check inputs

% pitch
if ~isfield(orient,'p'); 
    orient.p = zeros(size(LS.B1, 1));
    disp(['[INFO] - Pitch correction NOT applied.'])
end

% roll
if ~isfield(orient,'r'); 
    orient.r = zeros(size(LS.B1, 1));
    disp(['[INFO] - Roll correction NOT applied.'])
end

% heading
if ~isfield(orient,'h'); 
    orient.h = zeros(size(LS.B1, 1));
    disp(['[INFO] - heading correction NOT applied.'])
end


%% 1st Transform - to sensor axes
a = 1/(2*sind(theta));
b = 1/(4*cosd(theta));
d = a/(sqrt(2));

% % As perscribed
% LSt.X = (a.*LS.B1) - (a.*LS.B2);
% LSt.Y = -(a.*LS.B3) + (a.*LS.B4);
% LSt.Z = b.*(LS.B1 + LS.B2 + LS.B3 + LS.B4);
% LSt.Ze = d.*(LS.B1 + LS.B2) - d.*(LS.B3 + LS.B4);

% Magnitude only
LSt.X = (a.*abs(LS.B1)) + (a.*abs(LS.B2));
LSt.Y = (a.*abs(LS.B3)) + (a.*abs(LS.B4));
LSt.Z = b.*(abs(LS.B1) + abs(LS.B2) + abs(LS.B3) + abs(LS.B4));
LSt.Ze = d.*(abs(LS.B1) + abs(LS.B2)) + d.*(abs(LS.B3) + abs(LS.B4));

% Test
test = 0;
if test 
    LS.B1 = Data_out.B1.Lscale_ACF;
    LS.B2 = Data_out.B2.Lscale_ACF;
    LS.B3 = Data_out.B3.Lscale_ACF;
    LS.B4 = Data_out.B4.Lscale_ACF;
    
    orient.p = Data_out.ref.roll;
    orient.r = Data_out.ref.pitch;
    orient.h = Data_out.ref.heading;
        
    
    theta = 20;
    
end


%% 2nd Transform - to flow axes
% Rotation matrix [M] is of form:

%   (cosH*cosR + sinH*sinP*sinR)    (sinH*cosP)     (cosH*sinR - sinH*sinP*cosR)
%   (-sinH*cosR + cosH*sinP*sinR)   (cosH*cosP)     (-sinH*sinR - cosH*sinP*cosR)     = [M]
%          (-cosP*sinR)                (sinP)               (cosP*cosR)

% Note that for UPWARD looking sensor - add 180 degrees to measured roll
% (already done in CUBE data?)

% Loop through each ensemble, generate [M] for each step and find
% transformed lengthscales

if size(orient.h, 1) > 1
    orient.h = squeeze(nanmean(orient.h));
elseif size(orient.p, 1) > 1
    orient.p = squeeze(nanmean(orient.p));
elseif size(orient.r, 1) > 1 
    orient.r = squeeze(nanmean(orient.r));  
end

% Need to adjust heading to align with inflow? Or already done by Brian?

[ne nz] = size(LS.B1);

for j = 1:nz
    
    
    for i = 1:ne
        
        h = orient.h(i);
        p = orient.p(i);
        r = orient.r(i) ;
        
        M(1, 1:3) = [(cosd(h)*cosd(r)+sind(h)*sind(p)*sind(r)) (sind(h)*cosd(p)) (cosd(h)*sind(r) - sind(h)*sind(p)*cosd(r))];
        M(2, 1:3) = [(-sind(h)*cosd(r)+cosd(h)*sind(p)*sind(r)) (cosd(h)*cosd(p)) (-sind(h)*sind(r) - cosd(h)*sind(p)*cosd(r))];
        M(3, 1:3) = [(-cosd(p)*sind(r))  (sind(p))  (cosd(p)*cosd(r))];
        
        
        tempLS = [LSt.X(i,j); LSt.Y(i,j); LSt.Z(i,j)];
        
        transLS(:,i) = M*tempLS;
        
    end
    
    LSt.U(:,j) = transLS(1,:);
    LSt.V(:,j) = transLS(2,:);
    LSt.W(:,j) = transLS(3,:);
    
end

  figure
    pcolor(1:ne, 1:nz, abs(LSt.U'))
    colorbar; caxis([0 8]); shading flat;
    xlabel('Tstat ensemble'); ylabel('depth bin');
    title('transformed integral lengthscale [m]')
    
end


