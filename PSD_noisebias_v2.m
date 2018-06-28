function [noise_v, f, PSD_v, FFT_v, f_cut] = PSD_noisebias_v2(fluc_v, Fs, options)
%% Find Sensor Doppler Noise Bias
% 
% Based on clean_FFT & JBR_noise_code functions by Duncan Sutherland & Brian Sellar.
% Calcualtes acoustic sensor noise bias as per Richard (2013)
% 
%  INPUTS:
%   [fluc_v]:   A 3D array containing velocity fluctuations in 5-minute ensembles 
%               (Size = ns x ne x nz)
%   [Fs]:       Sampling frequency of acoustic sensor data - ***should be
%               frequency or period? - Confirm ***
%               (Size: 1 x 1)
%   [options]:  Struct of options for function operation (Optional)
%       options.range   - Range of PSD for fitting of f^-(5/3) slope. Leave
%                       empty for default.
%       options.c       - Weighting of each frequency component in fitting
%                       range
%       options.
%         
%  OUTPUTS:
%   [noise_v]: Noise bias values of detrended T-stat ensemble - VARIANCE  due to noise (m^2/s^2)
%               (Size = 1 x ne x nz)
%   [f]:     Frequency values at which FFT & PSS are evaluated
%               (Size = nf x 1)
%   [PSD_v]: Array of spectral density values
%               (Size = nf x ne x nz)
%   [FFT_v]: Array of FFT values - Single sided spectrum from FFT.
%               (Size = nf x ne x nz)
%   [f_cut]: Frequency at which noise begins to dominate PSS plot
%               (Size = 1 x ne x nz)
%
%  WHERE:   nf = number of frequencies evaluated in FFT
%           ns = number of samples in a T-stat period
%           ne = number of Tstat ensembles in data
%           nz = number of depth/distance cells from sensor
%           NOTE - single 5-min Tstat periods can be used if desired
%
%   v1  01/04/2018  First Version
%   v2  07/06/2018  GW - Changed method of calculating PSD to match clean_FFT -
%                   What method did first version use? Was different formula, gave different results...
%
%################################################################################################

%% Initialise 
L = size(fluc_v, 1);
n = 2^nextpow2(L);
% Frequency range
f = (Fs*(0:(n/2))/n)';

%% Use FFT to evaluate PSD

% Find FFT amplitude spectrum
P = fft(fluc_v, n, 1);
P2 = abs(P./n);                  % Compute double sided spectrum
P1 = P2(1 : n/2+1,:,:);             % Single-sided Amplitude spectrum
P1(2:end-1,:,:) = 2*P1(2:end-1,:,:);    % middle values must be doubled

% plot(f, P1)

% PSD spectra: 
PSD2 = (2.*abs(P).^2)/(n*Fs);               % double sided PSD spectra
PSD1 = PSD2(1:n/2+1, :, :);                      % single-sided PSD spectra
PSD1([1 end], :, :) = PSD1([1 end], :, :)/2;    % first & last values are not x2

% figure
% loglog(f, PSD1)

% Assign variables for output
FFT_v = P1;
PSD_v = PSD1;

% Averaged PSD for noise floor method (as investigated previously)
PSD_m = nanmean(PSD_v, 2);

% Range for fitting of f^-(5/3) slope
if ~isfield(options, 'range') || isempty(options.range) || isnan(options.range)
     [~,rangestart(1)]=min(abs(f-13^-1.0));
else 
    range = options.range;
end

range = rangestart(1)+1:1:length(f);

%% Find and treat large amplitudes at high freq end of PSD spectra.
% Replace with mean value of 'tail'

for i = 1:size(PSD_m, 3)
    
    assume_good = range(round(length(range)/2):end-15);
    tail_vals = PSD_m(assume_good,1,i);
    meanval = nanmean(tail_vals, 1);
    stdev = std(tail_vals); 
    
   % find 'out of bound' values
    oob = find(PSD_m(:, 1, i) > meanval+2*stdev);
    oob2 = oob(oob >= size(PSD_m, 1) - 15);
    
    PSD_m(oob2, 1, i) = meanval;
    
   % NaN's dont work - gives NaNs in results 
   % testnan = PSD_m;
   % testnan(oob2, 1,i) = NaN;
   
%     testmean = PSD_m;
%     testmean(oob2, 1, i) = meanval;
%     
    
%     %test plots
%     figure
%     loglog(f, PSD_m)
%     hold
%     loglog(f(assume_good), tail_vals)
%     loglog(f, repmat(meanval, length(f),1))
%     loglog(f, repmat(meanval + 2*stdev, length(f),1))
    
    
    
end

%% Use PSD to evaluate Sensor Noise Floor

% Set weigthing vector (Default is equal weighting of 1)
if ~isfield(options, 'c') || isempty(options.c) || sum(isnan(options.c))
    c=repmat(1,length(range),1);
elseif length(options.c) == length(range)
    c = options.c;
else
    c = ones(length(range));
end

% Set up LH matrix A
A11=sum(c);
A12=sum(c.*(f(range).^-(5/3)));
A21=sum(c.*(f(range).^-(5/3)));
A22=sum(c.*(f(range).^-(10/3)));
A=[A11,A12;A21,A22];

% Set up RH matrix B
f2=f(range).^-(5/3);

% CAN USE BSXFUN FOR 3D ARRAYS??? (seems to work - GW)
cPSD1 = bsxfun(@times, c, PSD_m(range,:,:)); 
cPSD1f2 = bsxfun(@times, cPSD1, f2);

B1=sum(cPSD1, 1);        
B2=sum(cPSD1f2, 1);

% to do 'all at once' with matrix functions...
B=[B1;B2];

for i = 1:size(B,2)
    for j = 1:size(B,3)
        
% Solve from A*X=B for X
X(:,i,j)=inv(A)*B(:,i,j);

N(i,j)=X(1, i, j);  % Why So many negative?

if N(i,j)<0 || X(2, i ,j)<0
N(i,j)=NaN;
end
K(i, j)=X(2, i, j);  % Integral subrange coeficient

int_range_fit=K(i,j)*f.^(-5/3);


% end loops
    end 
end

f_cut=(N./K).^(-3/5);

% outputs

noise_v = N.*f(end);

test_noise = sqrt(noise_v)*100;

% plot Results
plot_switch = 0;
if plot_switch==1
    dist = 64;
    figure
%     fig_h=figure;
    loglog(f(2:end),squeeze(PSD_m(2:end,1,dist))) % Original PSD
    hold all
    loglog(f(range),squeeze(PSD_m(range,1,dist))) % PSD of selected decade
    loglog(f,repmat(N(dist),1,length(f))) % Noise Floor
    loglog(f,int_range_fit) % 5/3 integral subrange fit
    loglog(f,repmat(N(dist),length(f),1)+int_range_fit) % Combined fit
    [~,loc]=min(abs(f-f_cut(dist))) ;
    loglog(f(assume_good),squeeze(PSD_m(assume_good,1,dist))) % PSD of selected decade
%     loglog([f_cut,f_cut],[N(1),Pss(loc)])
else
    fig_h=[];
end


end