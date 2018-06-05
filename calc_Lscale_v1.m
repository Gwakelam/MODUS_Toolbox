function [Lscale, ACF, Tscale] = calc_Lscale_v1(fluc_v, mean_v, Fs, lags, zc)
%% Calculate Integral Lengthscale
% Use autocorrelation method to calculate integral timescale (see
% Turbulence (Davidson, 2015)) 
%
%  INPUTS:
%   [fluc_v]:   A 3D array containing detrended velocity measurements in 5-minute T-stat ensembles 
%              (Size = ns x ne x nz)
%   [mean_v]:   A 2D array containing the mean velocity of 5-minute T-stat ensembles 
%              (Size = ne x nz)
%   [Fs]:       Sampling Frequency (Hz)
%              (Size = 1 x 1)
%   [lags]:     Number of lags up to which the autocorrelation function will process
%              (Size: 1 x 1)
%   [zc]:       Switch for zero-crossing interpolation (0 for off, 1 for
%               on)
%         
%  OUTPUTS:
%   [Lscale]:   2D matric of integral lengthscales 
%               (Size = ne x nz)
%   [ACF]:      3D matrix of autocorrelation functions for each 5 min Tstat
%               (Size = lags x ne x nz)
%   [Tscale]:   2D matrix of integral timescales 
%               (Size = ne x nz)
%
%  WHERE:  
%           ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%
% v1    03/04/2018      First Version
%################################################################################################

 %% Main Body
 
% test data setup
% fluc_v = Data_out.B1.detrended;
% mean_v = Data_out.B1.mean;
% Fs = Data_out.ref.Fs;
% zc = 0;
% lags = 125;
% end

ns = size(fluc_v, 1);
ne = size(fluc_v, 2);
nz = size(fluc_v, 3);
    
if ~exist('lags', 'var') || isempty(lags)
    lags = ns;
    disp('[INFO] - LAGS not specified, defaulting to N')
end
disp(strcat('lags = ', num2str(lags)));

if ~exist('zc', 'var') || isempty(zc)
    zc = 1;
    disp('[INFO] - ZC not specified, defaulting to ON')
end

% initialise results matrices
Lscale  = zeros( ne, nz);
ACF     = zeros( lags+1, ne, nz);
Tscale  = zeros( ne, nz);
nocross = 0;  


% Loop through Tstats & process
for i = 1:ne
    for j = 1:nz
        
        v_temp      = fluc_v(:, i, j);
        ACF_temp    = autocorr(v_temp, lags);
        
        % cumulative integral of ACF
        cumint = (1/Fs).*cumtrapz(ACF_temp);
        
        % Find zero down crossing (sign change + to -)
        snch  = diff(sign(ACF_temp));
        ind_down = find(snch < 0, 1);
        
        if isempty(ind_down)
            Tscale(i,j)   = NaN;
            Lscale(i,j)   = NaN;
            nocross = nocross+1;
            
        else
            if zc == 1
                %  Zero Crossing interpolation
                % ind_down is index of point BEFORE sign change
                yvals               = ACF_temp([ind_down ind_down+1]);
                yvals2(1, [1 3])   = ACF_temp([ind_down ind_down+1]);
                xvals               = [1 2];
                % interpolate timesteps over ACF values including 0 crossing
                % (gives fraction of a timestep to zero crossing)
                xvals2 = interp1(yvals, xvals, yvals2, 'linear');
                
                % Fraction of one sample period for zero crossing location
                zc_loc  = xvals2(2)-1;
                % trapezoidal method for area under zero crossing line - converted to seconds by
                % multiplying by sampling period (1/sample rate)
                zc_corr = 0.5*ACF_temp(ind_down)*zc_loc*(1/Fs);
                
                % corrected Tscale = cumul. sum up to last sample + zero crossing correction
                Tscale(i, j) = cumint(ind_down) + zc_corr;
            elseif zc == 0
                Tscale(i, j) = cumint(ind_down);
            end
            
            ACF(:,i,j) = ACF_temp;
        end
        
    end
end

% Add check for dimension match?
Lscale = Tscale .* mean_v;

if nocross >= 5
    fprintf('[WARNING - %d ACF''s do not cross zero! Consider increasing lags.\n', nocross)
end

%test = Tscale(:,:,1);

end
