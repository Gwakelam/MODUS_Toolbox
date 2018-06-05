function [Tref_2, Uref_2, hm0_2] = trim_refdata_v2(Tref, Uref, timestamp, Tstat)
% trim ref data to matchmeasured sensor data 
%
%  INPUTS
%   [Tref]: A 1D vector containing 5-min reference timestamps
%           (corresponding to each Tstat ensemmble)
%              (Size = nX) 
%   [Uref]: A 1D vector containing 5-min reference timestamps
%           (corresponding to each Tstat ensemmble)
%              (Size = nX) 
%   [timestamp]: A 1D vector of measurement timestamps
%              (size = nT)
%   [hm0]: hm0 data corresponding to each Tstat ensemble
%              (size = nX)
%   [Tstat]: Length of a period of stationarity, in minutes
%               (size = 1)
%
%  OUTPUTS:
%   [Tref_2]: 1D vector containing only values that correspond to a Tstat
%               ensemble for which measured data exists
%              (Size = ne) 
%   [Uref_2]: 1D vector containing only values that correspond to a Tstat
%               ensemble for which measured data exists
%              (Size = ne) 
%   [hm0_2]: hm0 data corresponding to each Tstat ensemble for which
%               measured data exists
%              (size = ne)
%
%  WHERE:   nX = number of 5-minute reference measurements before trimming
%           ne = number of Tstat ensembles in data set
%           nT = total number of measurements in dataset
%
% v1    04/04/2018      First Version
% v2    28/05/2018      Added check for correct number of Tref emnsemble GW
% v2    03/06/2018      Removed Hm0 trimming - now done in get_hm0 function GW

%################################################################################################

%% Ref Data Trimming

% USING INTERSECT WITH FLOATING POINT NUMBERS PROBLEMATIC - CONVERT DATENUMS TO DATEVECS USING DATEVEC(),
% THEN USE INTERSECT WITH 'ROWS' ARGUMENT TO MATCH TIMESTAMPS
Tref_dv       = datevec(Tref);
timestamp_dv  = datevec(timestamp);

% Find Tref values for which there is a matching data ensemble
[Match, indTref, indtst] = ...
    intersect(Tref_dv, unique(floor(timestamp_dv), 'rows'), 'rows');

% Sampling freq.
Fs =1/((timestamp(2) - timestamp(1))*(24*60*60));
% length of Tstat ensemble
Tlen = round((Tstat*60)*Fs);
% Number of Tstat ensembles
ne = numel(timestamp)/Tlen;

% If calculated number of Tstats doesnt match intersecting Tref values
if length(indTref) < ne 
    if length(Tref) >= indTref(end) +1
        indTref(end+1) = indTref(end) +1;
    else % For last data file produced by CUBE query, an extra Tref value generally needs to be added to end - workaround for now
        Tref(end+1) = Tref(end) + (Tstat/(60*24));
        Uref(end+1) = Uref(end);
        indTref(end+1) = indTref(end) +1;
    end
end
% If more Tstats, shorten. Should never be out by more than 1.
if length(indTref) > ne
    indTref = indTref(1:ne);
end

% Extract ref data matching sensor data - ref values are for stationary period preceding
% that timestamp? - in which case, these all need 'shifting' to index n-1
% in vector OR add extra value to BEGINNING of edges vector? - Needs to be
% clarified
Tref_2        = Tref(indTref);
Uref_2        = Uref(indTref);

end
    