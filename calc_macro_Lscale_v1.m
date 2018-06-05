function [macroLS, E] = calc_macro_Lscale_v1(vb, U, Fs, f_cut)
% Calculates integral length scale as per Roach (1986) and implemented by el-Gabry (2014) - 
%   - method used is as described by el-Gabry (2014) - "Procedure for
%   determining turbulence length scales...."

%  INPUTS:
% 
%    [vb]: Matrix of detrended velocity measurements fron a single beam, arranged in
%    ensembles of stationarity.
%           (Size = ns x ne x nz)
%
%    [Fs]: Sampling frequency of data, Hz (Size: 1)
%
%    [f_cut]: cutoff for low frequency averaging, Hz (Size = 1)
%
% OUTPUTS:
%    [macroLS]: Integral macro-length scale for each ensemble
%         (Size = ne x nz)
%
% WHERE:
%        ns = number of samples in a period of stationarity
%        ne = number of Tstat ensembles in data
%        nz = number of depth/distance cells from sensor
%        nF = number of frequencies the PSD/FFT function is evaluated at
%
% v1    30/05/2018  Dropped first value of Energy spectra (0 frequency - 0
%                   energy. Biases results low)

%% Check inputs

% Fs
if nargin < 3 
    disp('No cutoff for low freq supplied - defauting to 0.05 Hz')
    f_cut = 0.03;
end


%% Calculate

% Find energy as f >> 0, using fft
[ne, nz] = size(vb);
L   = 2^nextpow2(ne);     % More efficient when power of 2
f   = (Fs*(0:(L/2))/L)';          % Frequencies at which fft is evaluated
P   = abs(fft(vb, L, 1)).^2/L;    % Power
E   = P/Fs;                       % Energy

U = abs(U);
Uvar = squeeze(var(vb, 0, 1, 'omitnan'));

%figure
%loglog(f,(nanmean(E(1:length(f),:,20), 2)))

% Frquencies approaching 0
f0  = find( (f <= f_cut) & (f ~= 0));
% Ignore value for 0 freq.
Ef0_test = E(f0(2:end),:,:);
% Find mean energy
Ef0 = squeeze(nanmean(E(f0,:,:), 1));

% Calc macro length scale
macroLS = (Ef0.*U)./(4.*(Uvar));

%test
test = 0;
if test
    test_macroLS = squeeze(nanmean(macroLS, 1));
    z = 1:1: size(Ef0, 2);
    z = 2.1 + 1.*(z-1);
    figure
    plot(test_macroLS, z)
    xlabel('Macro Lengthscale [m]'); ylabel('mab [m]');
end

end
