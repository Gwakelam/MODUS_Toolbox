function [hm0] = get_hm0_v1(Tref)
% Finds Hm0 values from cut down CUBElite
%
%  INPUTS
%   [Tref]: A 1D vector containing 5-min reference timestamps
%           (corresponding to each Tstat ensemmble)
%              (Size = ne) 
%
%  OUTPUTS:
%   [hm0]: 1D vecrot of hm0 values corresponding to each T-stat ensemble
%               (Size = ns x ne x nz)
%
%  WHERE:   ne = number of Tstat ensembles in data set
%
% v1    04/04/2018      GW First Version
% v1    28/05/2018      GW Added intersect function to macth gaps internal to
%                       dataset (ie from missing data, other tide
%                       direction) 
%################################################################################################

% Load file with reference measurements 
load('C:\3. ReDAPT Matlab Work\Matlab_Query_Scripts\v1.1\input_data\cube_refs.mat');

% find start and end values, time range matching data
date_start = Tref(1);
date_end = Tref(end);

over = find(cube_refs.TS(:,1) >= datenum(date_start));
under = find(cube_refs.TS(:,1) <= datenum(date_end));

daterange = [over(1) under(end)];

% Get Hm0 values
hm0 = cube_refs.TS(daterange(1):daterange(2), 8);
% Corresponding timestamps
Thm0 = cube_refs.TS(daterange(1):daterange(2), 1);

% To eliminate additional Hm0 measurements in middle of dataset that dont match
% measured data (ie, on opposing tide, gaps in data etc.)
Tref_dv       = datevec(Tref);
Thm0_dv       = datevec(Thm0);

[Match, indTref, indHm0] = ...
    intersect(Tref_dv, unique(floor(Thm0_dv), 'rows'), 'rows');

hm0 = hm0(indHm0);

figure
plot(hm0)

