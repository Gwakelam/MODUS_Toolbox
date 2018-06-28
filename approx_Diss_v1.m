function [diss] = approx_Diss_v1(TKE,L,options)
% Approximate TKE dissipation. See George, link below, (p.65? )
% http://www.turbulence-online.com/Publications/Lecture_Notes/Turbulence_Lille/TB_16January2013.pdf
% 
%  INPUTS:
%    [TKE]: A structure containing TKE density values (Size = ne x nz)  
% 
%    [L]: Values of integral lengthscales (Size = ne x nz)
% 
%    [options]: A structure containing optional parameters
%       options.figure
%         [m]
%
% OUTPUTS:
%    [diss]: A matrix of TKE dissipation values through water column for each 
%               Tstat ensemble (Size = ne x nz)
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
if ~isfield(options,'figure')
    options.figure = 1;
end

%% TKE dissipation approximation

Diss = (2.*(TKE.^(3/2)))./L;

end
