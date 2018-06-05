function [diss, D, r] = calc_Diss_v1(fluc_v,z,options)
% Calcualte TKE dissipation using Osalusi (2009) structure function method.
% 
%  INPUTS:
%    [fluc_v]: A structure containing detrended along-beam velocity measurements (Size = ns x ne x nz)  
% 
%    [z]: bin values for distance from sensor, (Size = nz)
% 
%    [options]: A structure containing optional parameters
%         options.C      - constant for isotropic turbulence (Default 2.1)
%         options.disspts - number of points to use in sub-window (Default = 
%         options.figure - Flag to create figures (Default = 1)
%         options.surf_blank - distance to 'blank' below free surface [m]
%         options.bot_blank - distance to blank from sensor to first bin
%         [m]
%
% OUTPUTS:
%    [diss]: A matrix of TKE dissipation values through water column for each 
%               Tstat ensemble (Size = ne x nz)
%    [D]:    Matrix of structure function values for fluc_v ensembles
%               (Size = nr x ne x nz)
%    [r]:    Array of radial separation values for each z bin
%               (Size = nr x nz)
%
% WHERE:
%        nr = number of radial separations ( range: 0 - nz, ie. nr = nz)
%        ns = number of samples in a Tstat ensemble
%        ne = number of Tstat ensembles in data
%        nz = number of depth/distance cells from sensor
%
% v1    16/05/2018  Function crated - GW  
%                 

disp('Computing dissipation')

%% Check defaults

% C
if ~isfield(options,'C')
    options.C = 2.1;
end

% figure
if ~isfield(options,'figure')
    options.figure = 1;
end

% disspts
if ~isfield(options,'disspts')
    options.disspts = 4;
end

% surf_blank
if ~isfield(options,'surf_blank')
    options.surf_blank = 7;
end

% surf_blank
if ~isfield(options,'bot_blank')
    options.bot_blank = 7;
end

%% TESTING
test = 0;

if test
    fluc_v = Data_out.B1.detrended;
    z = Data_out.ref.zbin;
    options.surf_blank = 5;
    options.bot_blank = 5;
end

%% Centred Difference Structure Function
[D, r] = structureFunction_v1(fluc_v, z, options);



% cycle through Tstat ensembles
for k = 1:size(D, 2)
    
    D_temp = squeeze(D(:,k,:));

% cycle through depth bins
for j = 1:length(z)
    
%% Get 'Good' values from structure function

goodpts = ~isnan(D_temp(:,j));
surfpts = z' <= z(end) - options.surf_blank;
botpts = z' >= z(1) + options.bot_blank;
usepts = find(goodpts == 1 & surfpts == 1 & botpts == 1);

ruse = r(usepts, j);
Duse = D_temp(usepts, j);
Ur = unique(ruse);

if test
figure
plot(ruse, Duse, ':o' )
end
%% Least squares fit to structure function
% of the form A*r^(2/3) + n

[fitAux statsAux] = robustfit(ruse.^(2/3), Duse);
%Aaux( = 

if test
figure
plot(ruse, fitAux(2).*ruse.^(2/3)+fitAux(1), ruse, Duse, 'o')
%ylim([0 0.016])
end

% Check for flatness of structure func .^(-2/3) 
Dflat = Duse.*ruse.^(-2/3);
mDflat = nanmean(Dflat);
stdDflat = nanstd(Dflat);

% linear fit to D^(-2/3)
[fit stats] = robustfit(ruse, Dflat);
A(j) = fit(2);
N(j) = fit(1);
Dfit_flat(:, j) = N(j)+A(j)*ruse;
            RMSEfit(j)=sqrt(mean((Dflat-Dfit_flat(:,j)).^2));
            SlopeDiff(j)=(abs(A(j))); %Difference from zero slope
            
slopes(j)=A(j); 


residual(j) = stats.s;
Aerror(j) = stats.se(2);
Nerror(j) = stats.se(1);

if test
figure
plot(ruse, Dflat, ':o', ruse, Dfit_flat)
end

% end depth bin loop
end

%% Derive TKE dissipation from slope of fit, A
% as A = C*diss^(2/3)
posA = find(A >0);
diss(posA,k) = (A(posA) ./ options.C).^(3/2); % m^2/s^3 

% end ensemble loop
end




disp_linebreak('-')

%% plots

if options.figure
    %----------
    % plot diss
    %----------
    figure(1),clf
    ax(1)=subplot(411);
    pcolor(1:nt,z,U), shading flat
    colorbar
    caxis(2.2*[-1 1])
    ylabel('U')
    
    ax(2)=subplot(412);
    pcolor(1:nt,z,dUdz), shading flat
    colorbar
    caxis(0.12*[-1 1])
    ylabel('dU/dz')
    
    linkaxes(ax,'x')
    
    
end
